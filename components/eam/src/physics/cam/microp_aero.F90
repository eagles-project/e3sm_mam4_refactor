
module microp_aero

!---------------------------------------------------------------------------------
! Purpose:
!   CAM driver layer for aerosol activation processes.
!
! ***N.B.*** This module is currently hardcoded to recognize only the aerosols/modes that
!            affect the climate calculation.  This is implemented by using list
!            index 0 in all the calls to rad_constituent interfaces.
!
! Author: Andrew Gettelman
! Based on code from: Hugh Morrison, Xiaohong Liu and Steve Ghan
! May 2010
! Description in: Morrison and Gettelman, 2008. J. Climate (MG2008)
!                 Gettelman et al., 2010 J. Geophys. Res. - Atmospheres (G2010)         
! for questions contact Andrew Gettelman  (andrew@ucar.edu)
! Modifications: A. Gettelman Nov 2010  - changed to support separation of 
!                  microphysics and macrophysics and concentrate aerosol information here
!                B. Eaton, Sep 2014 - Refactored to move CAM interface code into the CAM
!                  interface modules and preserve just the driver layer functionality here.
!
!---------------------------------------------------------------------------------

use shr_kind_mod,     only: r8=>shr_kind_r8
use spmd_utils,       only: masterproc
use ppgrid,           only: pcols, pver, pverp, begchunk, endchunk
use ref_pres,         only: top_lev => trop_cloud_top_lev
use physconst,        only: rair, gravit, pi
use constituents,     only: cnst_get_ind
use physics_types,    only: physics_state, physics_ptend, physics_ptend_init
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field
use phys_control,     only: phys_getopts, cam_chempkg_is, use_hetfrz_classnuc
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_aer_mmr, rad_cnst_get_aer_props, &
                            rad_cnst_get_mode_num

use nucleate_ice_cam, only: use_preexisting_ice, nucleate_ice_cam_readnl, nucleate_ice_cam_register, &
                            nucleate_ice_cam_init, nucleate_ice_cam_calc
use ndrop,            only: ndrop_init, dropmixnuc
use ndrop_bam,        only: ndrop_bam_init, ndrop_bam_run, ndrop_bam_ccn

use hetfrz_classnuc_cam, only: hetfrz_classnuc_cam_readnl,  hetfrz_classnuc_cam_init, &
                               hetfrz_classnuc_cam_calc

use cam_history,      only: addfld, add_default, outfld
use cam_logfile,      only: iulog
use cam_abortutils,       only: endrun
use perf_mod,         only: t_startf, t_stopf

use error_messages, only: alloc_err

implicit none
private
save

public :: microp_aero_init, microp_aero_run, microp_aero_readnl, microp_aero_register

! Private module data

character(len=16)   :: eddy_scheme
logical             :: micro_do_icesupersat

!!   icenul_wsub_scheme = 1 : f(TKE) as default
!!                        2 : Mean updraft calculated from Gausssian PDF, with stddev=f(TKE)
integer             :: icenul_wsub_scheme = 1
 
! contact freezing due to dust
! dust number mean radius (m), Zender et al JGR 2003 assuming number mode radius of 0.6 micron, sigma=2
real(r8), parameter :: rn_dst1 = 0.258e-6_r8
real(r8), parameter :: rn_dst2 = 0.717e-6_r8
real(r8), parameter :: rn_dst3 = 1.576e-6_r8
real(r8), parameter :: rn_dst4 = 3.026e-6_r8

real(r8) :: bulk_scale    ! prescribed aerosol bulk sulfur scale factor

! smallest mixing ratio considered in microphysics
real(r8), parameter :: qsmall = 1.e-18_r8

! minimum allowed cloud fraction
real(r8), parameter :: mincld = 0.0001_r8

! indices in state%q and pbuf structures
integer :: cldliq_idx = -1
integer :: cldice_idx = -1
integer :: numliq_idx = -1
integer :: numice_idx = -1
integer :: kvh_idx = -1
integer :: tke_idx = -1
integer :: wp2_idx = -1
integer :: ast_idx = -1
integer :: alst_idx = -1
integer :: aist_idx = -1

integer :: cldo_idx = -1
integer :: dgnumwet_idx = -1
integer :: dgnum_idx = -1

integer :: naai_idx
integer :: naai_hom_idx

! pbuf indices for fields provided by heterogeneous freezing
integer :: frzimm_idx
integer :: frzcnt_idx
integer :: frzdep_idx

! Bulk aerosols
character(len=20), allocatable :: aername(:)
real(r8), allocatable :: num_to_mass_aer(:)

integer :: naer_all      ! number of aerosols affecting climate
integer :: idxsul   = -1 ! index in aerosol list for sulfate
integer :: idxdst2  = -1 ! index in aerosol list for dust2
integer :: idxdst3  = -1 ! index in aerosol list for dust3
integer :: idxdst4  = -1 ! index in aerosol list for dust4

! modal aerosols
logical :: clim_modal_aero

integer :: mode_accum_idx  = -1  ! index of accumulation mode
integer :: mode_aitken_idx = -1  ! index of aitken mode
integer :: mode_coarse_idx = -1  ! index of coarse mode
integer :: mode_coarse_dst_idx = -1  ! index of coarse dust mode
integer :: mode_coarse_slt_idx = -1  ! index of coarse sea salt mode
integer :: coarse_dust_idx = -1  ! index of dust in coarse mode
integer :: coarse_nacl_idx = -1  ! index of nacl in coarse mode
integer :: mode_fine_dst_idx = -1   ! index of dust in fine dust mode
integer :: mode_pcarbon_idx  = -1  ! index of dust in accum mode
integer :: accum_dust_idx    = -1  ! index of dust in accum mode
logical :: dem_in            = .false.

integer :: npccn_idx, rndst_idx, nacon_idx

logical  :: separate_dust = .false.
logical  :: liqcf_fix
real(r8), parameter :: unset_r8   = huge(1.0_r8)
real(r8) :: wsubmin = unset_r8 !PMA sets a much lower lower bound

integer, parameter :: ncnst = 20
integer :: hetfrz_aer_spec_idx(1:ncnst) = -1

character(len=8) :: hetfrz_aer_specname(ncnst)

! Copy of cloud borne aerosols before modification by droplet nucleation
! The basis is converted from mass to volume.
real(r8), allocatable :: aer_cb(:,:,:,:)

contains
!=========================================================================================

subroutine microp_aero_register
   !----------------------------------------------------------------------- 
   ! 
   ! Purpose: 
   ! Register pbuf fields for aerosols needed by microphysics
   ! 
   ! Author: Cheryl Craig October 2012
   ! 
   !-----------------------------------------------------------------------
   use ppgrid,         only: pcols
   use physics_buffer, only: pbuf_add_field, dtype_r8

   call pbuf_add_field('NPCCN',      'physpkg',dtype_r8,(/pcols,pver/), npccn_idx)

   call pbuf_add_field('RNDST',      'physpkg',dtype_r8,(/pcols,pver,4/), rndst_idx)
   call pbuf_add_field('NACON',      'physpkg',dtype_r8,(/pcols,pver,4/), nacon_idx)

   call pbuf_add_field('NAAI',     'physpkg', dtype_r8, (/pcols,pver/), naai_idx)
   call pbuf_add_field('NAAI_HOM', 'physpkg', dtype_r8, (/pcols,pver/), naai_hom_idx)   
 
   ! pbuf fields provided by hetfrz_classnuc
   call pbuf_add_field('FRZIMM', 'physpkg', dtype_r8, (/pcols,pver/), frzimm_idx)
   call pbuf_add_field('FRZCNT', 'physpkg', dtype_r8, (/pcols,pver/), frzcnt_idx)
   call pbuf_add_field('FRZDEP', 'physpkg', dtype_r8, (/pcols,pver/), frzdep_idx)

end subroutine microp_aero_register

!=========================================================================================

subroutine microp_aero_init

   !----------------------------------------------------------------------- 
   ! 
   ! Purpose: 
   ! Initialize constants for aerosols needed by microphysics
   ! 
   ! Author: Andrew Gettelman May 2010
   ! 
   !-----------------------------------------------------------------------

   ! local variables
   integer  :: iaer, ierr, ispec, istat
   integer  :: m, n, nmodes, nspec

   character(len=32) :: str32
   character(len=*), parameter :: routine = 'microp_aero_init'
   logical :: history_amwg
   !-----------------------------------------------------------------------

   ! Query the PBL eddy scheme
   call phys_getopts(eddy_scheme_out          = eddy_scheme, &
        history_amwg_out = history_amwg, &
        micro_do_icesupersat_out = micro_do_icesupersat, &
        liqcf_fix_out    = liqcf_fix,    & 
        demott_ice_nuc_out = dem_in      ) 
   
   if(masterproc)write(iulog,*)'DEMOTT is:', dem_in 

   hetfrz_aer_specname(1:ncnst) = (/'so4_c1  ', 'bc_c1   ', 'pom_c1  ', 'soa_c1  ', &
                                    'dst_c1  ', 'ncl_c1  ', 'mom_c1  ', 'num_c1  ', &
                                    'dst_c3  ', 'ncl_c3  ', 'so4_c3  ', 'bc_c3   ', &
                                    'pom_c3  ', 'soa_c3  ', 'mom_c3  ', 'num_c3  ', &
                                    'bc_c4   ', 'pom_c4  ', 'mom_c4  ', 'num_c4  '/)

   ! Access the physical properties of the aerosols that are affecting the climate
   ! by using routines from the rad_constituents module.

   ! get indices into state and pbuf structures
   call cnst_get_ind('CLDLIQ', cldliq_idx)
   call cnst_get_ind('CLDICE', cldice_idx)
   call cnst_get_ind('NUMLIQ', numliq_idx)
   call cnst_get_ind('NUMICE', numice_idx)

   select case(trim(eddy_scheme))
   case ('diag_TKE')
      tke_idx      = pbuf_get_index('tke')   
   case ('CLUBB_SGS')
      wp2_idx = pbuf_get_index('WP2_nadv')
   case default
      kvh_idx      = pbuf_get_index('kvh')
   end select

   ! clim_modal_aero determines whether modal aerosols are used in the climate calculation.
   ! The modal aerosols can be either prognostic or prescribed.
   call rad_cnst_get_info(0, nmodes=nmodes)
   clim_modal_aero = (nmodes > 0)

   ast_idx      = pbuf_get_index('AST')
   alst_idx      = pbuf_get_index('ALST')
   aist_idx      = pbuf_get_index('AIST')
   
   do ispec = 1, ncnst
      hetfrz_aer_spec_idx(ispec) = pbuf_get_index(hetfrz_aer_specname(ispec))
   enddo   

   ! Allocate space for copy of cloud borne aerosols before modification by
   ! droplet nucleation.
   !pw: zero-basing the lchunk index to work around PGI bug/feature.
   allocate(aer_cb(pcols,pver,ncnst,0:(endchunk-begchunk)), stat=istat)
   call alloc_err(istat, routine, 'aer_cb', pcols*pver*ncnst*(endchunk-begchunk+1))

   if (clim_modal_aero) then

      cldo_idx     = pbuf_get_index('CLDO')
      dgnumwet_idx = pbuf_get_index('DGNUMWET')
      dgnum_idx    = pbuf_get_index('DGNUM' )      

      call ndrop_init()

      ! Init indices for specific modes/species

      ! mode index for specified mode types
      do m = 1, nmodes
         call rad_cnst_get_info(0, m, mode_type=str32)
         select case (trim(str32))
         case ('accum')
            mode_accum_idx = m
         case ('aitken')
            mode_aitken_idx = m
         case ('coarse')
            mode_coarse_idx = m
         case ('coarse_dust')
            mode_coarse_dst_idx = m
         case ('coarse_seasalt')
            mode_coarse_slt_idx = m
         case ('fine_dust')
            mode_fine_dst_idx = m
         case ('primary_carbon')
            mode_pcarbon_idx  = m            
         end select
      end do

      ! check if coarse dust is in separate mode
      separate_dust = mode_coarse_dst_idx > 0

      ! for 3-mode 
      if ( mode_coarse_dst_idx<0 ) mode_coarse_dst_idx = mode_coarse_idx
      if ( mode_coarse_slt_idx<0 ) mode_coarse_slt_idx = mode_coarse_idx

      ! Check that required mode types were found
      if (mode_accum_idx == -1 .or. mode_aitken_idx == -1 .or. &
          mode_coarse_dst_idx == -1.or. mode_coarse_slt_idx == -1) then
         write(iulog,*) routine//': ERROR required mode type not found - mode idx:', &
            mode_accum_idx, mode_aitken_idx, mode_coarse_dst_idx, mode_coarse_slt_idx
         call endrun(routine//': ERROR required mode type not found')
      end if

      ! species indices for specified types
      ! find indices for the dust and seasalt species in the coarse mode
      call rad_cnst_get_info(0, mode_coarse_dst_idx, nspec=nspec)
      do n = 1, nspec
         call rad_cnst_get_info(0, mode_coarse_dst_idx, n, spec_type=str32)
         select case (trim(str32))
         case ('dust')
            coarse_dust_idx = n
         end select
      end do
      call rad_cnst_get_info(0, mode_coarse_slt_idx, nspec=nspec)
      do n = 1, nspec
         call rad_cnst_get_info(0, mode_coarse_slt_idx, n, spec_type=str32)
         select case (trim(str32))
         case ('seasalt')
            coarse_nacl_idx = n
         end select
      end do
      

      ! Check that required mode specie types were found
      if ( coarse_dust_idx == -1 .or. coarse_nacl_idx == -1) then
         write(iulog,*) routine//': ERROR required mode-species type not found - indicies:', &
            coarse_dust_idx, coarse_nacl_idx
         call endrun(routine//': ERROR required mode-species type not found')
      end if

   end if

   call addfld('LCLOUD', (/ 'lev' /), 'A', ' ', 'Liquid cloud fraction used in stratus activation')

   call addfld('WSUB',  (/ 'lev' /), 'A', 'm/s', 'Diagnostic sub-grid vertical velocity'                   )
   call addfld('WSUBI', (/ 'lev' /), 'A', 'm/s', 'Diagnostic sub-grid vertical velocity for ice'           )

   call addfld('WLARGE',(/ 'lev' /), 'A', 'm/s', 'Large-scale vertical velocity'                           )
   call addfld('WSIG',  (/ 'lev' /), 'A', 'm/s', 'Subgrid standard deviation of vertical velocity'         )
   call addfld('WSUBI2',(/ 'lev' /), 'A', 'm/s', 'Mean updraft, with stddev=f(TKE)'                        )
   call addfld('RHICE', (/ 'lev' /), 'A', '0-1', 'RHi for ice nucleation'                                  )

   if (history_amwg) then
      call add_default ('WSUB     ', 1, ' ')
   end if

   call nucleate_ice_cam_init(mincld, bulk_scale)
   call hetfrz_classnuc_cam_init(mincld)

end subroutine microp_aero_init

!=========================================================================================

subroutine microp_aero_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Namelist variables
   real(r8) :: microp_aero_bulk_scale = 2._r8  ! prescribed aerosol bulk sulfur scale factor
   real(r8) :: microp_aero_wsubmin  = 0.2_r8
   integer  :: microp_aero_wsub_scheme = 1     ! updraft velocity parameterization option for ice nucleation
 
   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'microp_aero_readnl'

   namelist /microp_aero_nl/ microp_aero_bulk_scale, microp_aero_wsub_scheme,microp_aero_wsubmin

   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'microp_aero_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, microp_aero_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variable
   call mpibcast(microp_aero_bulk_scale, 1, mpir8, 0, mpicom)
   call mpibcast(microp_aero_wsub_scheme, 1, mpiint, 0, mpicom)
   call mpibcast(microp_aero_wsubmin,     1, mpir8, 0, mpicom)
#endif

   ! set local variables
   bulk_scale = microp_aero_bulk_scale
   icenul_wsub_scheme = microp_aero_wsub_scheme
   wsubmin = microp_aero_wsubmin

   call nucleate_ice_cam_readnl(nlfile)
   call hetfrz_classnuc_cam_readnl(nlfile)

end subroutine microp_aero_readnl

!=========================================================================================

subroutine microp_aero_run ( &
   state, ptend, deltatin, pbuf, liqcldfo )
   
   use mam_support, only: min_max_bound
   use modal_aero_data,   only: qqcw_get_field, nspec_amode, numptrcw_amode, lmassptrcw_amode, maxd_aspectype, ntot_amode
   use ndrop,             only: ptr2d_t, mam_idx, ncnst_tot
  
   ! input arguments
   type(physics_state), target, intent(in)    :: state
   type(physics_ptend),         intent(out)   :: ptend
   real(r8),                    intent(in)    :: deltatin     ! time step (s)
   real(r8),                    intent(in)    :: liqcldfo(pcols,pver)  ! old liquid cloud fraction
   type(physics_buffer_desc),   pointer       :: pbuf(:)

   ! local workspace
   ! all units mks unless otherwise stated

   integer :: icol, kk, mm
   integer :: itim_old
   integer :: nmodes 
   integer :: lchnk_zb                  ! zero-based local chunk id
   integer :: ispec
   integer :: lspec   ! index for aerosol number / chem-mass / water-mass
   integer :: imode   ! aerosol mode index
   integer :: icnst   ! tracer index
   integer :: kvh_idx_dropmixnuc ! pbuf index of kvh needed for dropmixnuc input

   ! pbuf pointers 
   real(r8), pointer :: ast(:,:)        
   real(r8), pointer :: alst(:,:)        
   real(r8), pointer :: aist(:,:)        

   real(r8), pointer :: npccn(:,:)      ! number of CCN (liquid activated)

   real(r8), pointer :: kvh(:,:)        ! vertical eddy diff coef (m2 s-1)
   real(r8), pointer :: tke(:,:)        ! TKE from the UW PBL scheme (m2 s-2)
   real(r8), pointer :: wp2(:,:)        ! CLUBB vertical velocity variance

   real(r8), pointer :: cldn(:,:)       ! cloud fraction
   real(r8), pointer :: cldo(:,:)       ! old cloud fraction

   real(r8), pointer :: dgnumwet(:,:,:) ! aerosol mode diameter
   real(r8), pointer :: dgnum(:,:,:)

   ! naai and naai_hom are the outputs from nucleate_ice_cam_calc shared with the microphysics
   real(r8), pointer :: naai(:,:)       ! number of activated aerosol for ice nucleation [#/kg]
   real(r8), pointer :: naai_hom(:,:)   ! number of activated aerosol for ice nucleation (homogeneous freezing only) [#/kg]

   ! the following is used in droplet nucleation
   type(ptr2d_t), allocatable :: qqcw(:)     ! cloud-borne aerosol mass, number mixing ratios [#/kg or kg/kg]

   real(r8), pointer :: frzimm(:,:)
   real(r8), pointer :: frzcnt(:,:)
   real(r8), pointer :: frzdep(:,:)

   real(r8), pointer :: ptr2d(:,:)

   real(r8)          :: icecldf(pcols,pver)    ! ice cloud fraction   
   real(r8)          :: liqcldf(pcols,pver)    ! liquid cloud fraction

   real(r8) :: rho(pcols,pver)     ! air density (kg m-3)

   real(r8) :: lcldm(pcols,pver)   ! liq cloud fraction

   real(r8) :: lcldn(pcols,pver)   ! fractional coverage of new liquid cloud
   real(r8) :: lcldo(pcols,pver)   ! fractional coverage of old liquid cloud
   real(r8) :: qcld                ! total cloud water
   real(r8) :: nctend_mixnuc(pcols,pver)
   real(r8) :: dum, dum2           ! temporary dummy variable


   real(r8) :: qs(pcols)            ! liquid-ice weighted sat mixing rat (kg/kg)
   real(r8) :: es(pcols)            ! liquid-ice weighted sat vapor press (pa)
   real(r8) :: gammas(pcols)        ! parameter for cond/evap of cloud water

   real(r8) :: wsub(pcols,pver)    ! diagnosed sub-grid vertical velocity st. dev. (m/s)
   real(r8) :: wsubi(pcols,pver)   ! diagnosed sub-grid vertical velocity ice (m/s)
   real(r8) :: wsubice(pcols,pver) ! final updraft velocity for ice nucleation (m/s)
   real(r8) :: wsig(pcols,pver)    ! diagnosed standard deviation of vertical velocity ~ f(TKE)
   real(r8) :: nucboast

   real(r8) :: w0(pcols,pver)      ! large scale velocity (m/s) 
   real(r8) :: w2(pcols,pver)      ! subgrid mean updraft velocity, Gaussian PDF, stddev=f(tke)


   real(r8), allocatable :: factnum(:,:,:) ! activation fraction for aerosol number
   !-------------------------------------------------------------------------------

   associate( &
      lchnk => state%lchnk,             &
      ncol  => state%ncol,              &
      psetcols  => state%psetcols,      &
      temperature     => state%t,       &
      state_q         => state%q,       &
      qc    => state%q(:pcols,:pver,cldliq_idx), &
      qi    => state%q(:pcols,:pver,cldice_idx), &
      nc    => state%q(:pcols,:pver,numliq_idx), &
! BJG below might be equal to above but not sure about indices
      ncldwtr => state%q(:,:,numliq_idx),   &
      omega => state%omega,             &
      pmid     => state%pmid,           &
      pint     => state%pint,           &
      pdel     => state%pdel,           &
      rpdel    => state%rpdel,          &
      zm       => state%zm             )


   call t_startf('microp_aero_run_init')

   ! note for C++ porting, the following variables that are obtained using pbuf
   ! should be obtained from AD for the previous time step
   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, ast_idx,      ast, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, alst_idx,     alst, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, aist_idx,     aist, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, ast_idx,  cldn, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
   call pbuf_get_field(pbuf, cldo_idx, cldo, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) ) 
   call pbuf_get_field(pbuf, wp2_idx, wp2, start=(/1,1,itim_old/),kount=(/pcols,pverp,1/))  

   ! note for C++ porting, the following variables that are obtained using pbuf
   ! should be obtained from AD
   call pbuf_get_field(pbuf, npccn_idx, npccn) 
   call rad_cnst_get_info(0, nmodes=nmodes)
   call pbuf_get_field(pbuf, dgnumwet_idx, dgnumwet, start=(/1,1,1/), kount=(/pcols,pver,nmodes/) )
   call pbuf_get_field(pbuf, dgnum_idx, dgnum) 
 
   call pbuf_get_field(pbuf, naai_idx, naai)
   call pbuf_get_field(pbuf, naai_hom_idx, naai_hom)

   ! frzimm, frzcnt, frzdep are the outputs of hetfrz_classnuc_cam_calc used by the microphysics
   call pbuf_get_field(pbuf, frzimm_idx, frzimm)
   call pbuf_get_field(pbuf, frzcnt_idx, frzcnt)
   call pbuf_get_field(pbuf, frzdep_idx, frzdep)

   liqcldf(:ncol,:pver) = alst(:ncol,:pver) 
   icecldf(:ncol,:pver) = aist(:ncol,:pver)

   allocate(factnum(pcols,pver,nmodes))


   ! initialize output
   npccn(1:ncol,1:pver)    = 0._r8  

   lchnk_zb = lchnk - begchunk

   ! save copy of cloud borne aerosols for use in heterogeneous freezing
   !call hetfrz_classnuc_cam_save_cbaero(state, pbuf)
   do ispec = 1, ncnst
      call pbuf_get_field(pbuf, hetfrz_aer_spec_idx(ispec), ptr2d)
      
      aer_cb(:,:,ispec,lchnk_zb) = ptr2d
   enddo
  
   ! initialize time-varying parameters
   do kk = top_lev, pver
      do icol = 1, ncol
         rho(icol,kk) = pmid(icol,kk)/(rair*temperature(icol,kk))
      end do
   end do

   do ispec = 1, ncnst
      aer_cb(:ncol,:,ispec,lchnk_zb) = aer_cb(:ncol,:,ispec,lchnk_zb) * rho(:ncol,:)
   enddo

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! More refined computation of sub-grid vertical velocity 
   ! Set to be zero at the surface by initialization.

   allocate(tke(pcols,pverp))
   tke(:ncol,:) = (3._r8/2._r8)*wp2(:ncol,:)

   !PMA no longer needs the minimum value that is designed for CAM5-UW scheme which 
   !produces very low values

   wsub(:ncol,:top_lev-1)  = wsubmin
   wsubi(:ncol,:top_lev-1) = 0.001_r8
   wsig(:ncol,:top_lev-1)  = 0.001_r8

   do kk = top_lev, pver
      do icol = 1, ncol

         wsub(icol,kk)  = sqrt(0.5_r8*(tke(icol,kk) + tke(icol,kk+1))*(2._r8/3._r8))      
         wsig(icol,kk)  = min_max_bound(0.001_r8, 10._r8, wsub(icol,kk))        
         wsubi(icol,kk) = min_max_bound(0.2_r8, 10._r8, wsub(icol,kk))

         wsub(icol,kk)  = max(wsubmin, wsub(icol,kk))

      end do
   end do

   !!.......................................................... 
   !! Initialization
   !!.......................................................... 

   w0(1:ncol,1:pver) = 0._r8
   w2(1:ncol,1:pver) = 0._r8
   wsubice(1:ncol,1:pver) = 0._r8

   !!.......................................................... 
   !!  Convert from omega to w 
   !!  Negative omega means rising motion
   !!.......................................................... 

   do kk = top_lev, pver
      do icol = 1, ncol
         w0(icol,kk) = -1._r8*omega(icol,kk)/(rho(icol,kk)*gravit)
      enddo
   enddo

   call t_stopf('microp_aero_run_init')

   !!.......................................................... 
   !! icenul_wsub_scheme = 2 : Mean updraft calculated from Gausssian PDF, with
   !stddev=f(TKE)    
   !!.......................................................... 

   call t_startf('subgrid_mean_updraft')
   call subgrid_mean_updraft(ncol, w0, wsig, w2)
   call t_stopf('subgrid_mean_updraft')

   wsubice(1:ncol,1:pver) = wsubi(1:ncol,1:pver)

   call outfld('WSUB',   wsub, pcols, lchnk)
   call outfld('WSUBI',  wsubice, pcols, lchnk)
   call outfld('WSIG',   wsig, pcols, lchnk)
   call outfld('WLARGE', w0, pcols, lchnk)
   call outfld('WSUBI2', w2, pcols, lchnk)



   if (trim(eddy_scheme) == 'CLUBB_SGS') deallocate(tke)

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !ICE Nucleation

   call t_startf('nucleate_ice_cam_calc')
   call nucleate_ice_cam_calc(ncol, lchnk, temperature, state_q, pmid, &      ! input
                              rho, wsubice, ast, dgnum, &                     ! input
                              naai, naai_hom)                                 ! output
   call t_stopf('nucleate_ice_cam_calc')


   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! Droplet Activation

   ! for modal aerosol

   ! partition cloud fraction into liquid water part
   lcldn = 0._r8
   lcldo = 0._r8
   do kk = top_lev, pver
      do icol = 1, ncol
         qcld = qc(icol,kk) + qi(icol,kk)
         if (qcld > qsmall) then
            lcldn(icol,kk)=liqcldf(icol,kk)
            lcldo(icol,kk)=liqcldfo(icol,kk)
         end if
      end do
   end do
    
   call outfld('LCLOUD', lcldn, pcols, lchnk)

! BJG below could be moved to top, but safest to include here for now.

   ! Init qqcw (a flat 1d pointer array) to mode number and species mass mixing ratios for
   ! cloud borne phases.

   allocate(qqcw(ncnst_tot)) 
   do imode = 1, ntot_amode
      do lspec = 0, nspec_amode(imode)  ! loop through all species for mode 'imode'
          mm = mam_idx(imode,lspec)
          if (lspec == 0) then   ! number
             icnst = numptrcw_amode(imode)
          else ! aerosol mass
             icnst = lmassptrcw_amode(lspec,imode)
          endif
          qqcw(mm)%fld => qqcw_get_field(pbuf,icnst,lchnk,.true.) 
      enddo
   enddo

   ! NOTE FOR C++ porting, the following variable obtained using pbuf is used
   ! for vertical mixing in the activation subroutine dropmixnuc
   ! Do not use 'kvh_idx' since that index is not necessarily set in this module

   kvh_idx_dropmixnuc      = pbuf_get_index('kvh')
   call pbuf_get_field(pbuf, kvh_idx_dropmixnuc, kvh)

   call t_startf('dropmixnuc')
   call dropmixnuc( &
         lchnk,ncol,psetcols,deltatin,temperature,pmid,pint,pdel,rpdel,zm, &  ! in
         state_q,ncldwtr,kvh,wsub,lcldn, lcldo, &  ! in
         qqcw, &  ! inout
         ptend, nctend_mixnuc, factnum)  !out
   call t_stopf('dropmixnuc')

   npccn(:ncol,:) = nctend_mixnuc(:ncol,:)


   ! heterogeneous freezing
   call t_startf('hetfrz_classnuc_cam_calc')
   call hetfrz_classnuc_cam_calc(ncol, lchnk, temperature, pmid, rho, ast, &   ! in
                                 qc, nc, state_q, aer_cb(:,:,:,lchnk_zb), deltatin, factnum, & ! in
                                 frzimm, frzcnt, frzdep)                       ! out
   call t_stopf('hetfrz_classnuc_cam_calc')

   deallocate(factnum)

   end associate

end subroutine microp_aero_run

!=========================================================================================

subroutine subgrid_mean_updraft(ncol, w0, wsig, ww)

!---------------------------------------------------------------------------------
! Purpose: Calculate the mean updraft velocity inside a GCM grid assuming the 
!          vertical velocity distribution is Gaussian and peaks at the 
!          GCM resolved large-scale vertical velocity. 
!          When icenul_wsub_scheme = 2, the model uses the mean updraft velocity as the 
!          characteristic updraft velocity to calculate the ice nucleation rate. 
! Author:  Kai Zhang (kai.zhang@pnnl.gov) 
! Last Modified: Oct, 2015 
!---------------------------------------------------------------------------------

   !! interface 

   integer,  intent(in) :: ncol              ! number of cols 
   real(r8), intent(in) :: wsig(pcols,pver ) ! standard deviation (m/s)
   real(r8), intent(in) :: w0(pcols,pver ) ! large scale vertical velocity (m/s) 
   real(r8), intent(out):: ww(pcols,pver) ! mean updraft velocity(m/s) -> characteristic w*

   !! local 
   integer, parameter :: nbin = 50

   real(r8) :: wlarge,sigma
   real(r8) :: xx, yy 
   real(r8) :: zz(nbin) 
   real(r8) :: wa(nbin) 
   integer  :: kp(nbin) 
   integer  :: i, k
   integer  :: ibin

   !! program begins 

   do k = 1, pver
   do i = 1, ncol

      sigma  = max(0.001_r8, wsig(i,k))
      wlarge = w0(i,k)

      xx = 6._r8 * sigma / nbin

      do ibin = 1, nbin
         yy = wlarge - 3._r8*sigma + 0.5*xx
         yy = yy + (ibin-1)*xx
         !! wbar = integrator < w * f(w) * dw > 
         zz(ibin) = yy * exp(-1.*(yy-wlarge)**2/(2*sigma**2))/(sigma*sqrt(2*pi))*xx
      end do 

      kp(:) = 0 
      wa(:) = 0._r8 
 
      where(zz.gt.0._r8) 
         kp = 1 
         wa = zz
      elsewhere 
         kp = 0 
         wa = 0._r8 
      end where 

      if(sum(kp).gt.0) then 
         !! wbar = integrator < w * f(w) * dw > 
         ww(i,k) = sum(wa)
      else 
         ww(i,k) = 0.001_r8
      end if 

      !!write(6,*) 'i, k, w0, wsig, ww : ', i, k, w0(i,k), wsig(i,k), ww(i,k) 

  end do
  end do

end subroutine subgrid_mean_updraft
!================================================================================================

end module microp_aero
