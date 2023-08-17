
module modal_aer_opt

! parameterizes aerosol coefficients using chebychev polynomial
! parameterize aerosol radiative properties in terms of
! surface mode wet radius and wet refractive index

! Ghan and Zaveri, JGR 2007.

! uses Wiscombe's (1979) mie scattering code


use shr_kind_mod,      only: r8 => shr_kind_r8, shr_kind_cl
use ppgrid,            only: pcols, pver, pverp
use constituents,      only: pcnst
use spmd_utils,        only: masterproc
use phys_control,      only: cam_chempkg_is
use ref_pres,          only: top_lev => clim_modal_aero_top_lev
use physconst,         only: rhoh2o, rga, rair
use radconstants,      only: nswbands, nlwbands, idx_sw_diag, idx_uv_diag, idx_nir_diag
use rad_constituents,  only: n_diag, rad_cnst_get_call_list, rad_cnst_get_info, rad_cnst_get_aer_mmr, &
                             rad_cnst_get_aer_props, rad_cnst_get_mode_props
use physics_types,     only: physics_state

use physics_buffer, only : pbuf_get_index,physics_buffer_desc,pbuf_get_field,pbuf_old_tim_idx
use pio,               only: file_desc_t, var_desc_t, pio_inq_dimlen, pio_inq_dimid, pio_inq_varid, &
                             pio_get_var, pio_nowrite, pio_closefile
use cam_pio_utils,     only: cam_pio_openfile
use cam_history,       only:  addfld, horiz_only, add_default, outfld
use cam_history_support, only: fillvalue
use cam_logfile,       only: iulog
use perf_mod,          only: t_startf, t_stopf
use cam_abortutils,        only: endrun

use modal_aero_wateruptake, only: modal_aero_wateruptake_dr
use modal_aero_calcsize,    only: modal_aero_calcsize_diag,modal_aero_calcsize_sub
use shr_log_mod ,           only: errmsg => shr_log_errmsg

implicit none
private
save

public :: modal_aer_opt_readnl, modal_aer_opt_init, modal_aero_sw, modal_aero_lw


character(len=*), parameter :: unset_str = 'UNSET'

! Namelist variables:
character(shr_kind_cl)      :: modal_optics_file = unset_str   ! full pathname for modal optics dataset
character(shr_kind_cl)      :: water_refindex_file = unset_str ! full pathname for water refractive index dataset

! Dimension sizes in coefficient arrays used to parameterize aerosol radiative properties
! in terms of refractive index and wet radius
integer, parameter :: ncoef=5, prefr=7, prefi=10

real(r8) :: xrmin, xrmax

! refractive index for water read in read_water_refindex
complex(r8) :: crefwsw(nswbands) ! complex refractive index for water visible
complex(r8) :: crefwlw(nlwbands) ! complex refractive index for water infrared

! physics buffer indices
integer :: cld_idx      = 0
integer :: dgnumwet_idx = -1
integer :: qaerwat_idx  = -1

character(len=4) :: diag(0:n_diag) = (/'    ','_d1 ','_d2 ','_d3 ','_d4 ','_d5 ', &
                                       '_d6 ','_d7 ','_d8 ','_d9 ','_d10'/)

!Declare the following threadprivate variables to be used for calcsize and water uptake
!These are defined as module level variables to aviod allocation-deallocation in a loop
real(r8), allocatable, target :: dgnumdry_m(:,:,:) ! number mode dry diameter for all modes
real(r8), allocatable, target :: dgnumwet_m(:,:,:) ! number mode wet diameter for all modes
real(r8), allocatable, target :: qaerwat_m(:,:,:)  ! aerosol water (g/g) for all modes
!$OMP THREADPRIVATE(dgnumdry_m, dgnumwet_m, qaerwat_m)


!===============================================================================
CONTAINS
!===============================================================================

subroutine modal_aer_opt_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'modal_aer_opt_readnl'

   namelist /modal_aer_opt_nl/ water_refindex_file
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'modal_aer_opt_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, modal_aer_opt_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   call mpibcast(water_refindex_file, len(water_refindex_file), mpichar, 0, mpicom)
#endif


end subroutine modal_aer_opt_readnl

!===============================================================================

subroutine modal_aer_opt_init()

   use ioFileMod,        only: getfil
   use phys_control,     only: phys_getopts
   use shr_log_mod ,     only: errmsg => shr_log_errmsg

   ! Local variables

   integer  :: i, m
   real(r8) :: rmmin, rmmax       ! min, max aerosol surface mode radius treated (m)
   character(len=256) :: locfile

   logical           :: history_amwg            ! output the variables used by the AMWG diag package
   logical           :: history_verbose         ! produce verbose history output
   logical           :: history_aero_optics     ! output aerosol optics diagnostics

   logical :: call_list(0:n_diag)
   integer :: ilist, nmodes, m_ncoef, m_prefr, m_prefi
   integer :: errcode, istat

   character(len=*), parameter :: routine='modal_aer_opt_init'
   !----------------------------------------------------------------------------

   rmmin = 0.01e-6_r8
   rmmax = 25.e-6_r8
   xrmin = log(rmmin)
   xrmax = log(rmmax)

   ! set flags to check aerosol settings
   call rad_cnst_get_info(0, nmodes=nmodes)

   ! Check that dimension sizes in the coefficient arrays used to
   ! parameterize aerosol radiative properties are consistent between this
   ! module and the mode physprop files.
   call rad_cnst_get_call_list(call_list)
   do ilist = 0, n_diag
      if (call_list(ilist)) then
         call rad_cnst_get_info(ilist, nmodes=nmodes)
         do m = 1, nmodes
            call rad_cnst_get_mode_props(ilist, m, ncoef=m_ncoef, prefr=m_prefr, prefi=m_prefi)
            if (m_ncoef /= ncoef .or. m_prefr /= prefr .or. m_prefi /= prefi) then
               write(iulog,*) routine//': ERROR - file and module values do not match:'
               write(iulog,*) '   ncoef:', ncoef, m_ncoef
               write(iulog,*) '   prefr:', prefr, m_prefr
               write(iulog,*) '   prefi:', prefi, m_prefi
               call endrun(routine//': ERROR - file and module values do not match')
            end if
         end do
      end if
   end do

   cld_idx        = pbuf_get_index('CLD')


   ! Initialize physics buffer indices for dgnumwet and qaerwat.  Note the implicit assumption
   ! that the loops over modes in the optics calculations will use the values for dgnumwet and qaerwat
   ! that are set in the aerosol_wet_intr code.
   dgnumwet_idx = pbuf_get_index('DGNUMWET',errcode)
   if (errcode < 0) then
      call endrun(routine//' ERROR: cannot find physics buffer field DGNUMWET')
   end if
   qaerwat_idx  = pbuf_get_index('QAERWAT')
   if (errcode < 0) then
      call endrun(routine//' ERROR: cannot find physics buffer field QAERWAT')
   end if

   call getfil(water_refindex_file, locfile)
   call read_water_refindex(locfile)
   if (masterproc) write(iulog,*) "modal_aer_opt_init: read water refractive index file:", trim(locfile)

   call phys_getopts(history_amwg_out        = history_amwg, &
                     history_verbose_out = history_verbose, &
                     history_aero_optics_out = history_aero_optics )


   !obtain nmodes for the climate (ilist = 0) list
   call rad_cnst_get_info(0, nmodes=nmodes)

   !Allocate dry and wet size variables in an OMP PARRALLEL region as these
   !arrays are private for each thread and needs to be allocated for each OMP thread
   !$OMP PARALLEL
   allocate(dgnumdry_m(pcols,pver,nmodes),stat=istat)
   if (istat .ne. 0) call endrun("Unable to allocate dgnumdry_m: "//errmsg(__FILE__,__LINE__) )
   allocate(dgnumwet_m(pcols,pver,nmodes),stat=istat)
   if (istat .ne. 0) call endrun("Unable to allocate dgnumwet_m: "//errmsg(__FILE__,__LINE__) )
   allocate(qaerwat_m(pcols,pver,nmodes),stat=istat)
   if (istat .ne. 0) call endrun("Unable to allocate qaerwat_m: "//errmsg(__FILE__,__LINE__) )
   !$OMP END PARALLEL

   ! Add diagnostic fields to history output.

   call addfld ('EXTINCT',(/ 'lev' /),    'A','/m','Aerosol extinction', flag_xyfill=.true.)
   call addfld ('tropopause_m',horiz_only,    'A',' m  ','tropopause level in meters', flag_xyfill=.true.)
   call addfld ('ABSORB',(/ 'lev' /),    'A','/m','Aerosol absorption', flag_xyfill=.true.)
   call addfld ('AODVIS',horiz_only,    'A','  ','Aerosol optical depth 550 nm', flag_xyfill=.true., &
   standard_name='atmosphere_optical_thickness_due_to_ambient_aerosol_particles')
   call addfld ('AODALL',horiz_only,    'A','  ','AOD 550 nm for all time and species', flag_xyfill=.true.)
   call addfld ('AODUV',horiz_only,    'A','  ','Aerosol optical depth 350 nm', flag_xyfill=.true.)
   call addfld ('AODNIR',horiz_only,    'A','  ','Aerosol optical depth 850 nm', flag_xyfill=.true.)
   call addfld ('AODABS',horiz_only,    'A','  ','Aerosol absorption optical depth 550 nm', flag_xyfill=.true., &
   standard_name='atmosphere_absorption_optical_thickness_due_to_ambient_aerosol_particles')
   call addfld ('AODMODE1',horiz_only,    'A','  ','Aerosol optical depth 550 nm mode 1'           , flag_xyfill=.true.)
   call addfld ('AODMODE2',horiz_only,    'A','  ','Aerosol optical depth 550 nm mode 2'           , flag_xyfill=.true.)
   call addfld ('AODMODE3',horiz_only,    'A','  ','Aerosol optical depth 550 nm mode 3'           , flag_xyfill=.true.)
   call addfld ('AODDUST1',horiz_only,    'A','  ','Aerosol optical depth 550 nm model 1 from dust', flag_xyfill=.true.)
   call addfld ('AODDUST2',horiz_only,    'A','  ','Aerosol optical depth 550 nm model 2 from dust', flag_xyfill=.true.)
   call addfld ('AODDUST3',horiz_only,    'A','  ','Aerosol optical depth 550 nm model 3 from dust', flag_xyfill=.true.)
   call addfld ('AODDUST',horiz_only,    'A','  ','Aerosol optical depth 550 nm from dust', flag_xyfill=.true.)
   call addfld ('AODSO4',horiz_only,    'A','  ','Aerosol optical depth 550 nm from SO4', flag_xyfill=.true.)
   call addfld ('AODPOM',horiz_only,    'A','  ','Aerosol optical depth 550 nm from POM', flag_xyfill=.true.)
   call addfld ('AODSOA',horiz_only,    'A','  ','Aerosol optical depth 550 nm from SOA', flag_xyfill=.true.)
   call addfld ('AODBC',horiz_only,    'A','  ','Aerosol optical depth 550 nm from BC', flag_xyfill=.true.)
   call addfld ('AODSS',horiz_only,    'A','  ','Aerosol optical depth 550 nm from seasalt', flag_xyfill=.true.)
   call addfld ('AODABSBC',horiz_only, 'A','  ','Aerosol absorption optical depth 550 nm from BC', flag_xyfill=.true.)
   call addfld ('AODMOM',horiz_only,    'A','  ','Aerosol optical depth 550 nm from marine organic', flag_xyfill=.true.)
   call addfld ('BURDEN1',horiz_only,  'A','kg/m2'      ,'Aerosol burden mode 1'      , flag_xyfill=.true.)
   call addfld ('BURDEN2',horiz_only,  'A','kg/m2'      ,'Aerosol burden mode 2'      , flag_xyfill=.true.)
   call addfld ('BURDEN3',horiz_only,  'A','kg/m2'      ,'Aerosol burden mode 3'      , flag_xyfill=.true.)
   call addfld ('BURDENDUST',horiz_only,  'A','kg/m2'   ,'Dust aerosol burden'        , flag_xyfill=.true.)
   call addfld ('BURDENSO4',horiz_only,  'A','kg/m2'    ,'Sulfate aerosol burden'     , flag_xyfill=.true.)
   call addfld ('BURDENPOM',horiz_only,  'A','kg/m2'    ,'POM aerosol burden'         , flag_xyfill=.true.)
   call addfld ('BURDENSOA',horiz_only,  'A','kg/m2'    ,'SOA aerosol burden'         , flag_xyfill=.true.)
   call addfld ('BURDENBC',horiz_only,  'A','kg/m2'     ,'Black carbon aerosol burden', flag_xyfill=.true.)
   call addfld ('BURDENSEASALT',horiz_only,  'A','kg/m2','Seasalt aerosol burden'     , flag_xyfill=.true.)
   call addfld ('BURDENMOM',horiz_only,  'A','kg/m2'    ,'Marine organic aerosol burden', flag_xyfill=.true.)
   call addfld ('SSAVIS',horiz_only,    'A','  ','Aerosol singel-scatter albedo', flag_xyfill=.true.)


   if (history_amwg) then
      call add_default ('AODDUST1'     , 1, ' ')
      call add_default ('AODDUST3'     , 1, ' ')
      call add_default ('AODVIS'       , 1, ' ')
      call add_default ('AODALL'       , 1, ' ')
      call add_default ('BURDEN1'      , 1, ' ')
      call add_default ('BURDEN2'      , 1, ' ')
      call add_default ('BURDEN3'      , 1, ' ')
      if ( history_verbose ) then
      call add_default ('BURDENDUST'   , 1, ' ')
      call add_default ('BURDENSO4'    , 1, ' ')
      call add_default ('BURDENPOM'    , 1, ' ')
      call add_default ('BURDENSOA'    , 1, ' ')
      call add_default ('BURDENBC'     , 1, ' ')
      call add_default ('BURDENSEASALT', 1, ' ')
      call add_default ('BURDENMOM', 1, ' ')
      end if
   end if

   if (history_aero_optics) then
      call add_default ('AODDUST1'     , 1, ' ')
      call add_default ('AODDUST3'     , 1, ' ')
      call add_default ('AODMODE1'     , 1, ' ')
      call add_default ('AODMODE2'     , 1, ' ')
      call add_default ('AODMODE3'     , 1, ' ')
      call add_default ('AODVIS'       , 1, ' ')
      call add_default ('AODALL'       , 1, ' ')
      call add_default ('AODUV'        , 1, ' ')
      call add_default ('AODNIR'       , 1, ' ')
      call add_default ('AODABS'       , 1, ' ')
      call add_default ('AODABSBC'     , 1, ' ')
      call add_default ('AODDUST'      , 1, ' ')
      call add_default ('AODSO4'       , 1, ' ')
      call add_default ('AODPOM'       , 1, ' ')
      call add_default ('AODSOA'       , 1, ' ')
      call add_default ('AODBC'        , 1, ' ')
      call add_default ('AODSS'        , 1, ' ')
      call add_default ('BURDEN1'      , 1, ' ')
      call add_default ('BURDEN2'      , 1, ' ')
      call add_default ('BURDEN3'      , 1, ' ')
      if (history_verbose) then
      call add_default ('ABSORB'       , 1, ' ')
      call add_default ('BURDENDUST'   , 1, ' ')
      call add_default ('BURDENSO4'    , 1, ' ')
      call add_default ('BURDENPOM'    , 1, ' ')
      call add_default ('BURDENSOA'    , 1, ' ')
      call add_default ('BURDENBC'     , 1, ' ')
      call add_default ('BURDENSEASALT', 1, ' ')
      call add_default ('BURDENMOM'    , 1, ' ')
      end if
      call add_default ('SSAVIS'       , 1, ' ')
      call add_default ('EXTINCT'      , 1, ' ')
  end if
  if (cam_chempkg_is('trop_mam4').or.cam_chempkg_is('trop_mam4_resus').or. &
       cam_chempkg_is('trop_mam4_mom').or.cam_chempkg_is('trop_mam4_resus_mom').or. &
       cam_chempkg_is('trop_mam4_resus_soag').or.cam_chempkg_is('trop_mam7').or. &
       cam_chempkg_is('trop_mam9').or.cam_chempkg_is('trop_strat_mam7').or. &
       cam_chempkg_is('linoz_mam4_resus').or.cam_chempkg_is('linoz_mam4_resus_soag').or.&
       cam_chempkg_is('linoz_mam4_resus_mom').or. &
       cam_chempkg_is('linoz_mam4_resus_mom_soag').or.cam_chempkg_is('superfast_mam4_resus_mom_soag')) then
     call addfld ('AODDUST4',horiz_only,    'A','  ','Aerosol optical depth 550 nm model 4 from dust', flag_xyfill=.true.)
     call addfld ('AODMODE4',horiz_only,    'A','  ','Aerosol optical depth 550 nm mode 4', flag_xyfill=.true.)
     call addfld ('BURDEN4',horiz_only,    'A','kg/m2','Aerosol burden mode 4', flag_xyfill=.true.)


     if (history_aero_optics) then
        call add_default ('AODDUST4', 1, ' ')
        call add_default ('AODMODE4', 1, ' ')
        call add_default ('BURDEN4' , 1, ' ')
     end if
  end if
   if (cam_chempkg_is('trop_mam7').or.cam_chempkg_is('trop_mam9').or.cam_chempkg_is('trop_strat_mam7')) then
      call addfld ('AODDUST5',horiz_only,    'A','  ','Aerosol optical depth 550 nm model 5 from dust', flag_xyfill=.true.)
      call addfld ('AODDUST6',horiz_only,    'A','  ','Aerosol optical depth 550 nm model 6 from dust', flag_xyfill=.true.)
      call addfld ('AODDUST7',horiz_only,    'A','  ','Aerosol optical depth 550 nm model 7 from dust', flag_xyfill=.true.)
      call addfld ('AODMODE5',horiz_only,    'A','  ','Aerosol optical depth 550 nm mode 5', flag_xyfill=.true.)
      call addfld ('AODMODE6',horiz_only,    'A','  ','Aerosol optical depth 550 nm mode 6', flag_xyfill=.true.)
      call addfld ('AODMODE7',horiz_only,    'A','  ','Aerosol optical depth 550 nm mode 7', flag_xyfill=.true.)
      call addfld ('BURDEN5',horiz_only,    'A','kg/m2','Aerosol burden mode 5', flag_xyfill=.true.)
      call addfld ('BURDEN6',horiz_only,    'A','kg/m2','Aerosol burden mode 6', flag_xyfill=.true.)
      call addfld ('BURDEN7',horiz_only,    'A','kg/m2','Aerosol burden mode 7', flag_xyfill=.true.)

      if (history_aero_optics) then
         call add_default ('AODDUST5', 1, ' ')
         call add_default ('AODDUST6', 1, ' ')
         call add_default ('AODDUST7', 1, ' ')
         call add_default ('AODMODE5', 1, ' ')
         call add_default ('AODMODE6', 1, ' ')
         call add_default ('AODMODE7', 1, ' ')
         if (history_verbose) then
         call add_default ('BURDEN5', 1, ' ')
         call add_default ('BURDEN6', 1, ' ')
         call add_default ('BURDEN7', 1, ' ')
         end if
      end if
   end if
   if (cam_chempkg_is('trop_mam9')) then
      call addfld ('AODMODE8',horiz_only,    'A','  '  ,'Aerosol optical depth 550 nm mode 8', flag_xyfill=.true.)
      call addfld ('AODMODE9',horiz_only,    'A','  '  ,'Aerosol optical depth 550 nm mode 9', flag_xyfill=.true.)
      call addfld ('BURDEN8',horiz_only,    'A','kg/m2','Aerosol burden mode 8', flag_xyfill=.true.)
      call addfld ('BURDEN9',horiz_only,    'A','kg/m2','Aerosol burden mode 9', flag_xyfill=.true.)
      if (history_aero_optics) then
         call add_default ('AODMODE8', 1, ' ')
         call add_default ('AODMODE9', 1, ' ')
         call add_default ('BURDEN8', 1, ' ')
         call add_default ('BURDEN9', 1, ' ')
      end if
   end if

   do ilist = 1, n_diag
      if (call_list(ilist)) then

         call addfld ('EXTINCT'//diag(ilist), (/ 'lev' /), 'A','1/m', &
              'Aerosol extinction', flag_xyfill=.true.)
         call addfld ('ABSORB'//diag(ilist),  (/ 'lev' /), 'A','1/m', &
              'Aerosol absorption', flag_xyfill=.true.)
         call addfld ('AODVIS'//diag(ilist),       horiz_only, 'A','  ', &
              'Aerosol optical depth 550 nm', flag_xyfill=.true.)
         call addfld ('AODALL'//diag(ilist),       horiz_only, 'A','  ', &
              'AOD 550 nm all time', flag_xyfill=.true.)
         call addfld ('AODABS'//diag(ilist),       horiz_only, 'A','  ', &
              'Aerosol absorption optical depth 550 nm', flag_xyfill=.true.)

         call add_default ('EXTINCT'//diag(ilist), 1, ' ')
         call add_default ('ABSORB'//diag(ilist),  1, ' ')
         call add_default ('AODVIS'//diag(ilist),  1, ' ')
         call add_default ('AODALL'//diag(ilist),  1, ' ')
         call add_default ('AODABS'//diag(ilist),  1, ' ')

      end if
   end do

end subroutine modal_aer_opt_init

!===============================================================================

subroutine modal_aero_sw(dt, state, pbuf, nnite, idxnite, is_cmip6_volc, ext_cmip6_sw, trop_level,  &
                         tauxar, wa, ga, fa)
   ! calculates aerosol sw radiative properties

  use mam_support, only : min_max_bound
  use modal_aero_data,  only: ntot_amode, nspec_amode, specdens_amode, &
                              specname_amode, spechygro, & 
                              sigmag_amode, lspectype_amode, lmassptr_amode, &
                              specrefndxsw, refrtabsw, refitabsw, extpsw, abspsw, asmpsw

   real(r8),            intent(in) :: dt             !timestep (s)
   type(physics_state), intent(in), target :: state          ! state variables

   type(physics_buffer_desc), pointer :: pbuf(:)
   integer,             intent(in) :: nnite          ! number of night columns
   integer,             intent(in) :: idxnite(nnite) ! local column indices of night columns
   integer,             intent(in) :: trop_level(pcols)!tropopause level for each column
   real(r8),            intent(in) :: ext_cmip6_sw(pcols,pver)
   logical,             intent(in) :: is_cmip6_volc

   real(r8), intent(out) :: tauxar(pcols,0:pver,nswbands) ! layer extinction optical depth
   real(r8), intent(out) :: wa(pcols,0:pver,nswbands)     ! layer single-scatter albedo
   real(r8), intent(out) :: ga(pcols,0:pver,nswbands)     ! asymmetry factor
   real(r8), intent(out) :: fa(pcols,0:pver,nswbands)     ! forward scattered fraction

   ! Local variables
   integer :: i, ifld, isw, k, l, ll,m, nc, ns, ilev_tropp
   integer :: list_idx       ! index of the climate or a diagnostic list
   integer :: lchnk                    ! chunk id
   integer :: ncol                     ! number of active columns in the chunk
   integer :: nspec
   integer :: istat
   integer :: itim_old           ! index

   real(r8) :: mass(pcols,pver)        ! layer mass
   real(r8) :: air_density(pcols,pver) ! (kg/m3)

   real(r8),    pointer :: state_q(:,:,:)      ! state%q
   real(r8),    pointer :: temperature(:,:)    ! temperatures [K]
   real(r8),    pointer :: pmid(:,:)           ! layer pressure [Pa]
   real(r8),    pointer :: specmmr(:,:)        ! species mass mixing ratio
   character*32         :: spectype            ! species type
   real(r8)             :: hygro_aer           !

   real(r8), pointer :: cldn(:,:)        ! layer cloud fraction [fraction]
   real(r8), pointer :: dgnumwet(:,:)     ! number mode wet diameter
   real(r8), pointer :: qaerwat(:,:)      ! aerosol water (g/g)

   real(r8) :: sigma_logr_aer         ! geometric standard deviation of number distribution
   real(r8) :: radsurf(pcols,pver)    ! aerosol surface mode radius
   real(r8) :: logradsurf(pcols,pver) ! log(aerosol surface mode radius)
   real(r8) :: cheb(ncoef,pcols,pver)

   real(r8),allocatable :: volf(:,:)           ! volume fraction of insoluble aerosol
   real(r8),allocatable :: specdens(:)         ! species density for all species [kg/m3]
   complex(r8),allocatable :: specrefindex(:,:)     ! species refractive index

   real(r8)    :: refr(pcols)     ! real part of refractive index
   real(r8)    :: refi(pcols)     ! imaginary part of refractive index
   complex(r8) :: crefin(pcols)   ! complex refractive index

   real(r8) :: vol(pcols)      ! volume concentration of aerosol specie (m3/kg)
   real(r8) :: dryvol(pcols)   ! volume concentration of aerosol mode (m3/kg)
   real(r8) :: watervol(pcols) ! volume concentration of water in each mode (m3/kg)
   real(r8) :: wetvol(pcols)   ! volume concentration of wet mode (m3/kg)

   integer  :: itab(pcols), jtab(pcols)
   real(r8) :: ttab(pcols), utab(pcols)
   real(r8) :: cext(pcols,ncoef), cabs(pcols,ncoef), casm(pcols,ncoef)
   real(r8) :: pext(pcols)     ! parameterized specific extinction (m2/kg)
   real(r8) :: specpext(pcols) ! specific extinction (m2/kg)
   real(r8) :: dopaer(pcols)   ! aerosol optical depth in layer
   real(r8) :: pabs(pcols)     ! parameterized specific absorption (m2/kg)
   real(r8) :: pasm(pcols)     ! parameterized asymmetry factor
   real(r8) :: palb(pcols)     ! parameterized single scattering albedo

   ! Diagnostics
   real(r8) :: extinct(pcols,pver), tropopause_m(pcols)
   real(r8) :: absorb(pcols,pver)
   real(r8) :: aodvis(pcols)               ! extinction optical depth
   real(r8) :: aodall(pcols)               ! extinction optical depth
   real(r8) :: aodabs(pcols)               ! absorption optical depth

   real(r8) :: aodabsbc(pcols)             ! absorption optical depth of BC

   real(r8) :: ssavis(pcols)
   real(r8) :: dustvol(pcols)              ! volume concentration of dust in aerosol mode (m3/kg)

   real(r8) :: burden(pcols)
   real(r8) :: burdendust(pcols), burdenso4(pcols), burdenbc(pcols), &
               burdenpom(pcols), burdensoa(pcols), burdenseasalt(pcols)
   real(r8) :: burdenmom(pcols)

   real(r8) :: aodmode(pcols)
   real(r8) :: dustaodmode(pcols)          ! dust aod in aerosol mode

   real(r8) :: specrefr, specrefi
   real(r8) :: scatdust(pcols), scatso4(pcols), scatbc(pcols), &
               scatpom(pcols), scatsoa(pcols), scatseasalt(pcols)
   real(r8) :: scatmom(pcols)
   real(r8) :: absdust(pcols), absso4(pcols), absbc(pcols), &
               abspom(pcols), abssoa(pcols), absseasalt(pcols)
   real(r8) :: absmom(pcols)
   real(r8) :: hygrodust(pcols), hygroso4(pcols), hygrobc(pcols), &
               hygropom(pcols), hygrosoa(pcols), hygroseasalt(pcols)
   real(r8) :: hygromom(pcols)

   real(r8) :: scath2o, absh2o, sumscat, sumabs, sumhygro
   real(r8) :: aodc                        ! aod of component

   ! total species AOD
   real(r8) :: dustaod(pcols), so4aod(pcols), bcaod(pcols), &
               pomaod(pcols), soaaod(pcols), seasaltaod(pcols)
   real(r8) :: momaod(pcols)



   logical :: savaervis ! true if visible wavelength (0.55 micron)
   logical :: savaernir ! true if near ir wavelength (~0.88 micron)
   logical :: savaeruv  ! true if uv wavelength (~0.35 micron)

   real(r8) :: aoduv(pcols)               ! extinction optical depth in uv
   real(r8) :: aodnir(pcols)              ! extinction optical depth in nir


   character(len=32) :: outname

   ! debug output
   integer, parameter :: nerrmax_dopaer=1000
   integer  :: nerr_dopaer = 0
   character(len=*), parameter :: subname = 'modal_aero_sw'
   !----------------------------------------------------------------------------
   lchnk = state%lchnk
   ncol  = state%ncol
   state_q      => state%q
   temperature  => state%t
   pmid         => state%pmid

   mass(:ncol,:)        = state%pdeldry(:ncol,:)*rga
   air_density(:ncol,:) = state%pmid(:ncol,:)/(rair*state%t(:ncol,:))

   !FORTRAN refactoring: For prognostic aerosols only, other options are removed
   list_idx = 0   ! index of the climate or a diagnostic list
   itim_old    =  pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, cld_idx, cldn, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

   ! initialize output variables
   tauxar(:ncol,:,:) = 0._r8
   wa(:ncol,:,:)     = 0._r8
   ga(:ncol,:,:)     = 0._r8
   fa(:ncol,:,:)     = 0._r8
   ! zero'th layer does not contain aerosol
   tauxar(1:ncol,0,:)  = 0._r8
   wa(1:ncol,0,:)      = 0.925_r8
   ga(1:ncol,0,:)      = 0.850_r8
   fa(1:ncol,0,:)      = 0.7225_r8

   ! diagnostics for visible band summed over modes
   extinct(1:ncol,:)     = 0.0_r8
   absorb(1:ncol,:)      = 0.0_r8
   aodvis(1:ncol)        = 0.0_r8
   aodall(1:ncol)        = 0.0_r8
   aodabs(1:ncol)        = 0.0_r8
   burdendust(:ncol)     = 0.0_r8
   burdenso4(:ncol)      = 0.0_r8
   burdenpom(:ncol)      = 0.0_r8
   burdensoa(:ncol)      = 0.0_r8
   burdenbc(:ncol)       = 0.0_r8
   burdenseasalt(:ncol)  = 0.0_r8
   burdenmom(:ncol)      = 0.0_r8
   momaod(:ncol)         = 0.0_r8
   ssavis(1:ncol)        = 0.0_r8

   aodabsbc(:ncol)       = 0.0_r8
   dustaod(:ncol)        = 0.0_r8
   so4aod(:ncol)         = 0.0_r8
   pomaod(:ncol)         = 0.0_r8
   soaaod(:ncol)         = 0.0_r8
   bcaod(:ncol)          = 0.0_r8
   seasaltaod(:ncol)     = 0.0_r8


   ! diags for other bands
   aoduv(:ncol)          = 0.0_r8
   aodnir(:ncol)         = 0.0_r8

   ! Calculate aerosol size distribution parameters and aerosol water uptake
   !For prognostic aerosols
   call modal_aero_calcsize_sub(state, dt, pbuf, list_idx_in=list_idx, update_mmr_in = .false., &
           dgnumdry_m=dgnumdry_m)

   call modal_aero_wateruptake_dr(lchnk, ncol, state_q, temperature, pmid, & ! in 
                                  cldn, dgnumdry_m, & ! in
                                  dgnumwet_m, qaerwat_m, & ! inout
                                  list_idx_in=list_idx   ) ! optional in
   

   ! loop over all aerosol modes
   do m = 1, ntot_amode

      ! diagnostics for visible band for each mode
      burden(:ncol)       = 0._r8
      aodmode(1:ncol)     = 0.0_r8
      dustaodmode(1:ncol) = 0.0_r8

      dgnumwet => dgnumwet_m(:,:,m)
      qaerwat  => qaerwat_m(:,:,m)

      ! get mode info
      nspec = nspec_amode(m)
      sigma_logr_aer = sigmag_amode(m)

      ! calc size parameter for all columns
      call modal_size_parameters(ncol, sigma_logr_aer, dgnumwet, radsurf, logradsurf, cheb)


      allocate(volf(ncol,nspec),stat=istat)
      if (istat /= 0) call endrun("Unable to allocate volf: "//errmsg(__FILE__,__LINE__) )
      allocate(specdens(nspec),stat=istat)
      if (istat /= 0) call endrun("Unable to allocate specdens: "//errmsg(__FILE__,__LINE__) )
      allocate(specrefindex(nspec,nswbands),stat=istat)
      if (istat /= 0) call endrun("Unable to allocate specrefindex: "//errmsg(__FILE__,__LINE__) )

      do isw = 1, nswbands
         savaervis = (isw .eq. idx_sw_diag)
         savaeruv  = (isw .eq. idx_uv_diag)
         savaernir = (isw .eq. idx_nir_diag)

         do k = top_lev, pver

            ! form bulk refractive index
            crefin(:ncol) = (0._r8, 0._r8)
            dryvol(:ncol) = 0._r8
            dustvol(:ncol) = 0._r8

            scatdust(:ncol)     = 0._r8
            absdust(:ncol)      = 0._r8
            hygrodust(:ncol)    = 0._r8
            scatso4(:ncol)      = 0._r8
            absso4(:ncol)       = 0._r8
            hygroso4(:ncol)     = 0._r8
            scatbc(:ncol)       = 0._r8
            absbc(:ncol)        = 0._r8
            hygrobc(:ncol)      = 0._r8
            scatpom(:ncol)      = 0._r8
            abspom(:ncol)       = 0._r8
            hygropom(:ncol)     = 0._r8
            scatsoa(:ncol)      = 0._r8
            abssoa(:ncol)       = 0._r8
            hygrosoa(:ncol)     = 0._r8
            scatseasalt(:ncol)  = 0._r8
            absseasalt(:ncol)   = 0._r8
            hygroseasalt(:ncol) = 0._r8
            scatmom(:ncol)  = 0._r8
            absmom(:ncol)   = 0._r8
            hygromom(:ncol) = 0._r8

            ! aerosol species loop
            do l = 1, nspec

               ! get aerosol properties and save for each species
               specmmr => state_q(:,:,lmassptr_amode(l,m))
               spectype = specname_amode(lspectype_amode(l,m))
               hygro_aer = spechygro(lspectype_amode(l,m))
               specdens(l) = specdens_amode(lspectype_amode(l,m))
               specrefindex(l,:) = specrefndxsw(:,lspectype_amode(l,m))
               volf(:,l) = specmmr(:,k)/specdens(l)

               do i = 1, ncol
                  vol(i)      = specmmr(i,k)/specdens(l)
               end do

               ! compute some diagnostics for visible band only
               if (savaervis) then

                  specrefr = real(specrefindex(l,isw))
                  specrefi = aimag(specrefindex(l,isw))

                  do i = 1, ncol
                     burden(i) = burden(i) + specmmr(i,k)*mass(i,k)
                  enddo

                  if (trim(spectype) == 'dust') then
                        call calc_diag_spec ( ncol, specmmr(:,k), mass(:,k), & ! in
                                        vol, specrefr, specrefi, hygro_aer, & ! in
                                        burdendust, scatdust, absdust, hygrodust) ! out
                        dustvol(:)    = vol(:)
                  endif
                  if (trim(spectype) == 'sulfate') then
                        call calc_diag_spec ( ncol, specmmr(:,k), mass(:,k), & ! in
                                        vol, specrefr, specrefi, hygro_aer, & ! in
                                        burdenso4, scatso4, absso4, hygroso4) ! out
                  endif
                  if (trim(spectype) == 'black-c') then
                        call calc_diag_spec ( ncol, specmmr(:,k), mass(:,k), & ! in
                                        vol, specrefr, specrefi, hygro_aer, & ! in
                                        burdenbc, scatbc, absbc, hygrobc) ! out
                  endif
                  if (trim(spectype) == 'p-organic') then
                        call calc_diag_spec ( ncol, specmmr(:,k), mass(:,k), & ! in
                                        vol, specrefr, specrefi, hygro_aer, & ! in
                                        burdenpom, scatpom, abspom, hygropom) ! out
                  endif
                  if (trim(spectype) == 's-organic') then
                        call calc_diag_spec ( ncol, specmmr(:,k), mass(:,k), & ! in
                                        vol, specrefr, specrefi, hygro_aer, & ! in
                                        burdensoa, scatsoa, abssoa, hygrosoa) ! out
                  endif
                  if (trim(spectype) == 'seasalt') then
                        call calc_diag_spec ( ncol, specmmr(:,k), mass(:,k), & ! in
                                        vol, specrefr, specrefi, hygro_aer, & ! in
                                        burdenseasalt, scatseasalt, absseasalt, hygroseasalt) ! out
                  endif
                  if (trim(spectype) == 'm-organic') then
                        call calc_diag_spec ( ncol, specmmr(:,k), mass(:,k), & ! in
                                        vol, specrefr, specrefi, hygro_aer, & ! in
                                        burdenmom, scatmom, absmom, hygromom) ! out
                  endif

               endif ! if (savaervis)
            enddo ! species loop l

            call calc_refin_complex ('sw', ncol, isw,          & ! in
                   qaerwat(:,k), volf, specrefindex,           & ! in
                   dryvol, wetvol, watervol, crefin, refr, refi) ! out

            ! call t_startf('binterp')

            ! interpolate coefficients linear in refractive index
            ! first call calcs itab,jtab,ttab,utab
            itab(:ncol) = 0
            call binterp(extpsw(m,:,:,:,isw), ncol, ncoef, prefr, prefi, &
                         refr, refi, refrtabsw(m,:,isw), refitabsw(m,:,isw), &
                         itab, jtab, ttab, utab, cext)
            call binterp(abspsw(m,:,:,:,isw), ncol, ncoef, prefr, prefi, &
                         refr, refi, refrtabsw(m,:,isw), refitabsw(m,:,isw), &
                         itab, jtab, ttab, utab, cabs)
            call binterp(asmpsw(m,:,:,:,isw), ncol, ncoef, prefr, prefi, &
                         refr, refi, refrtabsw(m,:,isw), refitabsw(m,:,isw), &
                         itab, jtab, ttab, utab, casm)

            ! call t_stopf('binterp')

            ! parameterized optical properties
            call calc_parameterized (ncol, ncoef, cext, cheb(:,:,k), & ! in
                                     pext) ! out
            call calc_parameterized (ncol, ncoef, cabs, cheb(:,:,k), & ! in
                                     pabs) ! out
            call calc_parameterized (ncol, ncoef, casm, cheb(:,:,k), & ! in
                                     pasm) ! out
            do i=1,ncol

               if (logradsurf(i,k) <= xrmax) then
                  pext(i) = exp(pext(i))
               else
                  pext(i) = 1.5_r8/(radsurf(i,k)*rhoh2o) ! geometric optics
               endif
               ! convert from m2/kg water to m2/kg aerosol
               specpext(i) = pext(i)
               pext(i) = pext(i)*wetvol(i)*rhoh2o
               pabs(i) = pabs(i)*wetvol(i)*rhoh2o

               pabs(i) = min_max_bound(0._r8, pext(i), pabs(i))
               palb(i) = 1._r8-pabs(i)/max(pext(i),1.e-40_r8)
               dopaer(i) = pext(i)*mass(i,k)

            enddo

            if (savaeruv) then
               do i = 1, ncol
                  aoduv(i) = aoduv(i) + dopaer(i)
               enddo
            endif

            if (savaernir) then
               do i = 1, ncol
                  aodnir(i) = aodnir(i) + dopaer(i)
               enddo
            endif

            ! Save aerosol optical depth at longest visible wavelength
            ! sum over layers
            if (savaervis) then
               ! aerosol extinction (/m)
               do i = 1, ncol
                  extinct(i,k) = extinct(i,k) + dopaer(i)*air_density(i,k)/mass(i,k)
                  absorb(i,k)  = absorb(i,k) + pabs(i)*air_density(i,k)
                  aodvis(i)    = aodvis(i) + dopaer(i)
                  aodall(i)    = aodall(i) + dopaer(i)
                  aodabs(i)    = aodabs(i) + pabs(i)*mass(i,k)
                  aodmode(i)   = aodmode(i) + dopaer(i)
                  ssavis(i)    = ssavis(i) + dopaer(i)*palb(i)

                  if (wetvol(i) > 1.e-40_r8) then

                     dustaodmode(i) = dustaodmode(i) + dopaer(i)*dustvol(i)/wetvol(i)

                     ! partition optical depth into contributions from each constituent
                     ! assume contribution is proportional to refractive index X volume

                     scath2o        = watervol(i)*real(crefwsw(isw))
                     absh2o         = -watervol(i)*aimag(crefwsw(isw))
                     sumscat        = scatso4(i) + scatpom(i) + scatsoa(i) + scatbc(i) + &
                                      scatdust(i) + scatseasalt(i) + scath2o + &
                                      scatmom(i)
                     sumabs         = absso4(i) + abspom(i) + abssoa(i) + absbc(i) + &
                                      absdust(i) + absseasalt(i) + absh2o + &
                                      absmom(i)
                     sumhygro       = hygroso4(i) + hygropom(i) + hygrosoa(i) + hygrobc(i) + &
                                      hygrodust(i) + hygroseasalt(i) + &
                                      hygromom(i)

                     call update_aod_spec ( i, scath2o, absh2o, & ! in
                           sumhygro, sumscat, sumabs, & ! in
                           hygrodust(i), palb(i), dopaer(i), & ! in
                           scatdust, absdust, dustaod ) ! inout

                     call update_aod_spec ( i, scath2o, absh2o, & ! in
                           sumhygro, sumscat, sumabs, & ! in
                           hygroso4(i), palb(i), dopaer(i), & ! in
                           scatso4, absso4, so4aod ) ! inout

                     call update_aod_spec ( i, scath2o, absh2o, & ! in
                           sumhygro, sumscat, sumabs, & ! in
                           hygropom(i), palb(i), dopaer(i), & ! in
                           scatpom, abspom, pomaod ) ! inout

                     call update_aod_spec ( i, scath2o, absh2o, & ! in
                           sumhygro, sumscat, sumabs, & ! in
                           hygrosoa(i), palb(i), dopaer(i), & ! in
                           scatsoa, abssoa, soaaod ) ! inout

                     call update_aod_spec ( i, scath2o, absh2o, & ! in
                           sumhygro, sumscat, sumabs, & ! in
                           hygrobc(i), palb(i), dopaer(i), & ! in
                           scatbc, absbc, bcaod ) ! inout

                     call update_aod_spec ( i, scath2o, absh2o, & ! in
                           sumhygro, sumscat, sumabs, & ! in
                           hygroseasalt(i), palb(i), dopaer(i), & ! in
                           scatseasalt, absseasalt, seasaltaod ) ! inout

                     call update_aod_spec ( i, scath2o, absh2o, & ! in
                           sumhygro, sumscat, sumabs, & ! in
                           hygromom(i), palb(i), dopaer(i), & ! in
                           scatmom, absmom, momaod ) ! inout

                     aodabsbc(i)    = aodabsbc(i) + absbc(i)*dopaer(i)*(1.0_r8-palb(i))

                  endif

               enddo ! i
            endif

            do i = 1, ncol
                call check_error_warning('sw', i, k, m, isw, nspec,list_idx, & ! in
                        dopaer(i), pabs(i), dryvol, wetvol, watervol, crefin,cabs,& ! in
                        specdens, specrefindex, volf, & ! in
                        nerr_dopaer, & ! inout
                        pext(i), specpext(i) ) ! optional in

               tauxar(i,k,isw) = tauxar(i,k,isw) + dopaer(i)
               wa(i,k,isw)     = wa(i,k,isw)     + dopaer(i)*palb(i)
               ga(i,k,isw)     = ga(i,k,isw)     + dopaer(i)*palb(i)*pasm(i)
               fa(i,k,isw)     = fa(i,k,isw)     + dopaer(i)*palb(i)*pasm(i)*pasm(i)
            enddo

         enddo ! pver

      enddo ! sw bands

      ! mode diagnostics
      ! The diagnostics are currently only output for the climate list.  Code mods will
      ! be necessary to provide output for the rad_diag lists.
      if (list_idx == 0) then
         do i = 1, nnite
            aodmode(idxnite(i)) = fillvalue
            dustaodmode(idxnite(i)) = fillvalue
         enddo

         write(outname,'(a,i1)') 'BURDEN', m
         call outfld(trim(outname), burden, pcols, lchnk)

         write(outname,'(a,i1)') 'AODMODE', m
         call outfld(trim(outname), aodmode, pcols, lchnk)

         write(outname,'(a,i1)') 'AODDUST', m
         call outfld(trim(outname), dustaodmode, pcols, lchnk)
      endif

      deallocate(volf)
      deallocate(specdens)
      deallocate(specrefindex)

   enddo ! nmodes

   !Add contributions from volcanic aerosols directly read in extinction
   if(is_cmip6_volc) then
        call calc_volc_ext(ncol, trop_level, state%zm, ext_cmip6_sw, & ! in
                extinct, tropopause_m ) ! inout/out
   endif


   ! Output visible band diagnostics for quantities summed over the modes
   ! These fields are put out for diagnostic lists as well as the climate list.
   do i = 1, nnite
      extinct(idxnite(i),:) = fillvalue
      absorb(idxnite(i),:)  = fillvalue
      aodvis(idxnite(i))    = fillvalue
      aodabs(idxnite(i))    = fillvalue
   enddo

   call outfld('EXTINCT'//diag(list_idx),  extinct, pcols, lchnk)
   call outfld('tropopause_m', tropopause_m, pcols, lchnk)
   call outfld('ABSORB'//diag(list_idx),   absorb,  pcols, lchnk)
   call outfld('AODVIS'//diag(list_idx),   aodvis,  pcols, lchnk)
   call outfld('AODALL'//diag(list_idx),   aodall,  pcols, lchnk)
   call outfld('AODABS'//diag(list_idx),   aodabs,  pcols, lchnk)

   ! These diagnostics are output only for climate list
   if (list_idx == 0) then
      do i = 1, ncol
         if (aodvis(i) > 1.e-10_r8) then
            ssavis(i) = ssavis(i)/aodvis(i)
         else
            ssavis(i) = 0.925_r8
         endif
      enddo

      do i = 1, nnite
         ssavis(idxnite(i))     = fillvalue
         aoduv(idxnite(i))      = fillvalue
         aodnir(idxnite(i))     = fillvalue
         aodabsbc(idxnite(i))   = fillvalue
         dustaod(idxnite(i))    = fillvalue
         so4aod(idxnite(i))     = fillvalue
         pomaod(idxnite(i))     = fillvalue
         soaaod(idxnite(i))     = fillvalue
         bcaod(idxnite(i))      = fillvalue
         momaod(idxnite(i))     = fillvalue
       enddo

      call outfld('SSAVIS',        ssavis,        pcols, lchnk)
      call outfld('AODUV',         aoduv,         pcols, lchnk)
      call outfld('AODNIR',        aodnir,        pcols, lchnk)
      call outfld('BURDENDUST',    burdendust,    pcols, lchnk)
      call outfld('BURDENSO4' ,    burdenso4,     pcols, lchnk)
      call outfld('BURDENPOM' ,    burdenpom,     pcols, lchnk)
      call outfld('BURDENSOA' ,    burdensoa,     pcols, lchnk)
      call outfld('BURDENBC'  ,    burdenbc,      pcols, lchnk)
      call outfld('BURDENSEASALT', burdenseasalt, pcols, lchnk)
      call outfld('BURDENMOM',     burdenmom,     pcols, lchnk)
      call outfld('AODABSBC',      aodabsbc,      pcols, lchnk)
      call outfld('AODDUST',       dustaod,       pcols, lchnk)
      call outfld('AODSO4',        so4aod,        pcols, lchnk)
      call outfld('AODPOM',        pomaod,        pcols, lchnk)
      call outfld('AODSOA',        soaaod,        pcols, lchnk)
      call outfld('AODBC',         bcaod,         pcols, lchnk)
      call outfld('AODSS',         seasaltaod,    pcols, lchnk)
      call outfld('AODMOM',        momaod,        pcols, lchnk)
   endif  !  if (list_idx == 0)

end subroutine modal_aero_sw

!===============================================================================
subroutine modal_aero_lw(dt, state, pbuf, & ! in
                        tauxar            ) ! out

  use shr_log_mod ,     only: errmsg => shr_log_errmsg
  use modal_aero_data,  only: ntot_amode, nspec_amode, specdens_amode, &
                              sigmag_amode, lspectype_amode, lmassptr_amode, &
                              specrefndxlw, refrtablw, refitablw, absplw

   ! calculates aerosol lw radiative properties

   real(r8),            intent(in)  :: dt       ! time step [s]
   type(physics_state), intent(in), target :: state    ! state variables
   type(physics_buffer_desc), pointer :: pbuf(:)

   real(r8), intent(out) :: tauxar(pcols,pver,nlwbands) ! layer absorption optical depth

   ! Local variables
   integer :: icol, ilw, kk, ll, mm, nc
   integer :: lchnk
   integer :: list_idx                 ! index of the climate or a diagnostic list
   integer :: ncol                     ! number of active columns in the chunk
   integer :: nspec
   integer :: istat
   integer :: itim_old           ! index

   real(r8), pointer :: dgnumwet(:,:)  ! wet number mode diameter [m]
   real(r8), pointer :: qaerwat(:,:)   ! aerosol water [g/g]

   real(r8) :: sigma_logr_aer          ! geometric standard deviation of number distribution
   real(r8) :: alnsg_amode
   real(r8) :: cheby(ncoef,pcols,pver) ! chebychef polynomials
   real(r8) :: radsurf(pcols,pver)     ! aerosol surface mode radius
   real(r8) :: logradsurf(pcols,pver)  ! log(aerosol surface mode radius)

   real(r8) :: mass(pcols,pver) ! layer mass

   real(r8),allocatable :: volf(:,:)           ! volume fraction of insoluble aerosol
   real(r8),allocatable :: specdens(:)         ! species density for all species [kg/m3]
   complex(r8),allocatable :: specrefindex(:,:)     ! species refractive index

   real(r8) :: dryvol(pcols)    ! volume concentration of aerosol mode [m3/kg]
   real(r8) :: wetvol(pcols)    ! volume concentration of wet mode [m3/kg]
   real(r8) :: watervol(pcols)  ! volume concentration of water in each mode [m3/kg]
   real(r8) :: refr(pcols)      ! real part of refractive index
   real(r8) :: refi(pcols)      ! imaginary part of refractive index
   complex(r8) :: crefin(pcols) ! complex refractive index
   real(r8), pointer :: state_q(:,:,:)     ! state%q
   real(r8), pointer :: specmmr(:,:)       ! species mass mixing ratio [g/g]
   real(r8), pointer :: temperature(:,:)   ! temperatures [K]
   real(r8), pointer :: pmid(:,:)          ! layer pressure [Pa]
   real(r8), pointer :: cldn(:,:)        ! layer cloud fraction [fraction]

   integer  :: itab(pcols), jtab(pcols)
   real(r8) :: ttab(pcols), utab(pcols)
   real(r8) :: cabs(pcols,ncoef)
   real(r8) :: pabs(pcols)      ! parameterized specific absorption [m2/kg]
   real(r8) :: dopaer    ! aerosol optical depth in layer

   integer  :: nerr_dopaer = 0

   !----------------------------------------------------------------------------

   ncol  = state%ncol
   lchnk = state%lchnk
   state_q      => state%q
   temperature  => state%t
   pmid         => state%pmid
   ! dry mass in each cell
   mass(:ncol,:) = state%pdeldry(:ncol,:)*rga

   !FORTRAN refactoring: For prognostic aerosols only, other options are removed
   list_idx = 0   ! index of the climate or a diagnostic list
   itim_old    =  pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, cld_idx, cldn, start=(/1,1,itim_old/),   kount=(/pcols,pver,1/) )


   call modal_aero_calcsize_sub(state, dt, pbuf, list_idx_in=list_idx, update_mmr_in = .false., &
           dgnumdry_m=dgnumdry_m)

   call modal_aero_wateruptake_dr(lchnk, ncol, state_q, temperature, pmid, & ! in
                                  cldn, dgnumdry_m, & ! in
                                  dgnumwet_m, qaerwat_m, & ! inout
                                  list_idx_in=list_idx   ) ! optional in
   

   ! initialize output variables
   tauxar(:ncol,:,:) = 0._r8

   ! loop over all aerosol modes
   do mm = 1, ntot_amode

      dgnumwet => dgnumwet_m(:,:,mm)
      qaerwat  => qaerwat_m(:,:,mm)

      ! get mode info
      nspec = nspec_amode(mm)
      sigma_logr_aer = sigmag_amode(mm)
      
      ! calc size parameter for all columns
      ! FORTRAN refactoring: ismethod2 is tempararily used to ensure BFB test. 
      ! can be removed when porting to C++
      call modal_size_parameters(ncol, sigma_logr_aer, dgnumwet, & ! in
                                 radsurf, logradsurf, cheby, ismethod2=.true.) 

      allocate(volf(ncol,nspec),stat=istat)
      if (istat /= 0) call endrun("Unable to allocate volf: "//errmsg(__FILE__,__LINE__) )
      allocate(specdens(nspec),stat=istat)
      if (istat /= 0) call endrun("Unable to allocate specdens: "//errmsg(__FILE__,__LINE__) )
      allocate(specrefindex(nspec,nlwbands),stat=istat)
      if (istat /= 0) call endrun("Unable to allocate specrefindex: "//errmsg(__FILE__,__LINE__) )

      do ilw = 1, nlwbands

         do kk = top_lev, pver
            ! get aerosol properties and save for each species
            do ll = 1, nspec
               specmmr => state_q(:,:,lmassptr_amode(ll,mm))
               specdens(ll) = specdens_amode(lspectype_amode(ll,mm))
               specrefindex(ll,:) = specrefndxlw(:,lspectype_amode(ll,mm))
               volf(:,ll) = specmmr(:,kk)/specdens(ll)
            enddo

            ! calculate complex refractive index
            call calc_refin_complex('lw', ncol, ilw, & ! in
                   qaerwat(:,kk), volf, specrefindex,   & ! in
                   dryvol, wetvol, watervol, crefin, refr, refi) ! out

            ! interpolate coefficients linear in refractive index
            ! first call calcs itab,jtab,ttab,utab
            itab(:ncol) = 0
            call binterp(absplw(mm,:,:,:,ilw), ncol, ncoef, prefr, prefi, &
                         refr, refi, refrtablw(mm,:,ilw), refitablw(mm,:,ilw), &
                         itab, jtab, ttab, utab, cabs)

            ! parameterized optical properties
            call calc_parameterized (ncol, ncoef, cabs, cheby(:,:,kk), & ! in
                                     pabs) ! out
            do icol = 1, ncol
               pabs(icol)   = pabs(icol)*wetvol(icol)*rhoh2o
               pabs(icol)   = max(0._r8,pabs(icol))
               dopaer = pabs(icol)*mass(icol,kk)

               ! FORTRAN refactor: check and writeout error/warning message
               call check_error_warning('lw', icol, kk,mm, ilw, nspec, list_idx,& ! in
                        dopaer, pabs(icol), dryvol, wetvol, watervol, crefin,cabs,& ! in
                        specdens, specrefindex, volf, & ! in
                        nerr_dopaer) ! inout

               ! update absorption optical depth
               tauxar(icol,kk,ilw) = tauxar(icol,kk,ilw) + dopaer
            enddo

         enddo ! kk = top_lev, pver

      enddo  ! nlwbands

      deallocate(volf)
      deallocate(specdens)
      deallocate(specrefindex)

   enddo ! mm = 1, ntot_amode

end subroutine modal_aero_lw

!===============================================================================
! Private routines
!===============================================================================

subroutine calc_parameterized (ncol, ncoef, coef, cheb_k, & ! in
                               para) ! out
    ! calculate parameterized absorption, extinction or asymmetry factor
    ! further calculations are needed. see modal_aero_sw and modal_aero_lw

    implicit none
    integer,  intent(in) :: ncol,ncoef
    real(r8), intent(in) :: coef(pcols,ncoef)
    real(r8), intent(in) :: cheb_k(ncoef,pcols)
    real(r8), intent(out):: para(pcols)

    integer :: nc

    para(:ncol) = 0.5_r8*coef(:ncol,1)
    do nc = 2, ncoef
        para(:ncol) = para(:ncol) + cheb_k(nc,:ncol)*coef(:ncol,nc)
    enddo

end subroutine calc_parameterized
!===============================================================================
subroutine calc_diag_spec ( ncol, specmmr_k, mass_k, & ! in
                        vol, specrefr, specrefi, hygro_aer, & ! in
                        burden_s, scat_s, abs_s, hygro_s) ! out
   ! calculate some diagnostics for a species
   implicit none
   integer,  intent(in) :: ncol
   real(r8), intent(in),  pointer :: specmmr_k(:)
   real(r8), intent(in) :: mass_k(:)
   real(r8), intent(in) :: vol(:) ! volume concentration of aerosol species (m3/kg)
   real(r8), intent(in) :: specrefr, specrefi  ! real and image part of specrefindex
   real(r8), intent(in) :: hygro_aer  ! aerosol hygroscopicity
   real(r8), intent(out) :: burden_s(pcols) ! aerosol burden of species
   real(r8), intent(out) :: scat_s(pcols)   ! scattering of species
   real(r8), intent(out) :: abs_s(pcols)    ! absorption of species
   real(r8), intent(out) :: hygro_s(pcols)  ! hygroscopicity of species
   
   integer :: icol

   burden_s(:ncol) = 0._r8
   do icol = 1, ncol
       burden_s(icol) = burden_s(icol) + specmmr_k(icol)*mass_k(icol)
       scat_s(icol)   = vol(icol)*specrefr
       abs_s(icol)    = -vol(icol)*specrefi
       hygro_s(icol)  = vol(icol)*hygro_aer
   enddo

end subroutine calc_diag_spec

!===============================================================================
subroutine update_aod_spec ( icol, scath2o, absh2o, & ! in
                           sumhygro, sumscat, sumabs, & ! in
                           hygro_s, palb, dopaer, & ! in
                           scat_s, abs_s, aod_s ) ! inout
   ! update aerosol optical depth from scattering and absorption

   implicit none
   integer,  intent(in) :: icol
   real(r8), intent(in) :: scath2o, absh2o, sumscat, sumabs, sumhygro
   real(r8), intent(in) :: hygro_s, palb, dopaer
   real(r8), intent(inout) :: scat_s(:), abs_s(:), aod_s(:)  ! scatering, absorption and aod for a species
   ! local variables
   real(r8) :: aodc  ! aod component

   scat_s(icol)     = (scat_s(icol) + scath2o*hygro_s/sumhygro)/sumscat
   abs_s(icol)      = (abs_s(icol) + absh2o*hygro_s/sumhygro)/sumabs

   aodc           = (abs_s(icol)*(1.0_r8 - palb) + palb*scat_s(icol))*dopaer

   aod_s(icol)      = aod_s(icol) + aodc

end subroutine update_aod_spec

!===============================================================================
subroutine calc_volc_ext(ncol, trop_level, state_zm, ext_cmip6_sw, & ! in
                extinct, tropopause_m ) ! inout/out
   ! calculate contributions from volcanic aerosol extinction

   implicit none
   integer,  intent(in) :: ncol
   integer,  intent(in) :: trop_level(pcols)!tropopause level for each column
   real(r8), intent(in) :: state_zm(:,:) ! state%zm [m]
   real(r8), intent(in) :: ext_cmip6_sw(pcols,pver)
   real(r8), intent(inout) :: extinct(pcols,pver)
   real(r8), intent(out)   :: tropopause_m(pcols)

   ! local variables
   integer :: icol, kk_tropp

   do icol = 1, ncol
      kk_tropp = trop_level(icol)

      ! diagnose tropopause height
      tropopause_m(icol) = state_zm(icol,kk_tropp)!in meters

      !update tropopause layer first
      extinct(icol,kk_tropp) = 0.5_r8*( extinct(icol,kk_tropp) + ext_cmip6_sw(icol,kk_tropp) )
      !extinction is assigned read in values only for visible band above tropopause
      extinct(icol, 1:kk_tropp-1) = ext_cmip6_sw(icol, 1:kk_tropp-1)
   enddo

end subroutine calc_volc_ext

!===============================================================================
subroutine calc_refin_complex (lwsw, ncol, ilwsw,           & ! in
                qaerwat_kk,  volf, specrefindex,            & ! in
                dryvol, wetvol, watervol, crefin, refr, refi) ! out
    !-------------------------------------------------------------------
    ! calculate complex refractive index 
    ! also output wetvol and watervol
    !-------------------------------------------------------------------

    implicit none
    character(len=2), intent(in) :: lwsw   ! indicator if this is lw or sw
    integer,  intent(in) :: ncol, ilwsw       
    real(r8), intent(in) :: qaerwat_kk(:)   ! aerosol water at level kk [g/g]
    real(r8), intent(in) :: volf(:,:)       ! volume fraction of insoluble aerosol [fraction]
    complex(r8), intent(in) :: specrefindex(:,:)     ! species refractive index

    real(r8),intent(out) :: dryvol(pcols)    ! volume concentration of aerosol mode [m3/kg]
    real(r8),intent(out) :: wetvol(pcols)    ! volume concentration of wet mode [m3/kg]
    real(r8),intent(out) :: watervol(pcols)  ! volume concentration of water in each mode [m3/kg]
    real(r8),intent(out) :: refr(pcols)      ! real part of refractive index
    real(r8),intent(out) :: refi(pcols)      ! imaginary part of refractive index
    complex(r8),intent(out) :: crefin(pcols) ! complex refractive index

    integer :: icol

    if ((lwsw /= 'lw') .and. (lwsw /= 'sw')) call endrun()

    crefin(:ncol) = (0._r8, 0._r8)
    dryvol(:ncol) = 0._r8

    do icol = 1, ncol
       dryvol(icol) = sum(volf(icol,:))
       crefin(icol) = sum(volf(icol,:)*specrefindex(:,ilwsw))

       watervol(icol) = qaerwat_kk(icol)/rhoh2o
       wetvol(icol)   = watervol(icol) + dryvol(icol)

       if (watervol(icol) < 0.0_r8 .and. lwsw=='lw') then
          if (abs(watervol(icol)) > 1.e-1_r8*wetvol(icol)) then
              write(iulog,*) 'watervol,wetvol,dryvol=',watervol(icol), &
                         wetvol(icol),dryvol(icol)
          endif
          watervol(icol) = 0._r8
          wetvol(icol)   = dryvol(icol)
       endif

       ! some different treatments for lw and sw
       if (lwsw=='lw') then
          crefin(icol) = crefin(icol) + watervol(icol)*crefwlw(ilwsw)
          if (wetvol(icol) > 1.e-40_r8) crefin(icol) = crefin(icol)/wetvol(icol)
       elseif (lwsw=='sw') then
          crefin(icol) = crefin(icol) + watervol(icol)*crefwsw(ilwsw)
          crefin(icol) = crefin(icol)/max(wetvol(icol),1.e-60_r8)
       endif

       refr(icol) = real(crefin(icol))
       refi(icol) = aimag(crefin(icol))
    enddo

end subroutine calc_refin_complex

!===================================================================
subroutine check_error_warning(lwsw, icol, kk, mm, ilwsw, nspec,list_idx, & ! in
                        dopaer, pabs, dryvol, wetvol, watervol, crefin,cabs,& ! in
                        specdens, specrefindex, volf, & ! in
                        nerr_dopaer, & ! inout
                        pext, specpext ) ! optional in
    !------------------------------------------------------------
    ! check and writeout error and warning message
    !------------------------------------------------------------

   implicit none
   character(len=2), intent(in) :: lwsw   ! indicator if this is lw or sw
   integer,intent(in) :: icol, kk, mm, ilwsw, nspec, list_idx  ! indices
   real(r8),intent(in) :: dopaer    ! aerosol optical depth in layer
   real(r8),intent(in) :: pabs      ! parameterized specific absorption [m2/kg]
   real(r8),intent(in) :: dryvol(:)    ! volume concentration of aerosol mode [m3/kg]
   real(r8),intent(in) :: wetvol(:)    ! volume concentration of wet mode [m3/kg]
   real(r8),intent(in) :: watervol(:)  ! volume concentration of water in each mode [m3/kg]
   complex(r8),intent(in) :: crefin(:) ! complex refractive index
   real(r8),intent(in) :: cabs(:,:)
   real(r8),intent(in) :: volf(:,:)   ! volume fraction of insoluble aerosol [fraction]
   real(r8),intent(in) :: specdens(:) ! species density [kg/m3]
   complex(r8), intent(in) :: specrefindex(:, :)     ! species refractive index
   integer, intent(inout) :: nerr_dopaer    ! total number of error times
   real(r8),intent(in),optional :: pext         ! only write out for sw
   real(r8),intent(in),optional :: specpext     ! only write out for sw

   integer :: ll
   integer, parameter :: nerrmax_dopaer=1000

    if ((lwsw /= 'lw') .and. (lwsw /= 'sw')) call endrun()

    ! FORTRAN refactor: This if condition is never met in testing run ...
    if ((dopaer <= -1.e-10_r8) .or. (dopaer >= 20._r8)) then

        if (dopaer <= -1.e-10_r8) then
            write(iulog,*) "ERROR: Negative aerosol optical depth &
                          &in this layer."
        else
            write(iulog,*) "WARNING: Aerosol optical depth is &
                         &unreasonably high in this layer."
        endif

        write(iulog,*) 'dopaer(',icol,',',kk,',',mm,',)=', dopaer
        write(iulog,*) 'kk=',kk,' pabs=', pabs
        if ((lwsw == 'sw') .and. present(pext) .and. present(specpext)) then
            write(iulog,*) 'kk=', kk, ' pext=', pext, ' specext=', specpext
        endif
        write(iulog,*) 'wetvol=',wetvol(icol),' dryvol=',dryvol(icol),     &
                       ' watervol=',watervol(icol)
        write(iulog,*) 'cabs=', (cabs(icol,ll),ll=1,ncoef)
        write(iulog,*) 'crefin=', crefin(icol)
        write(iulog,*) 'nspec=', nspec
        do ll = 1,nspec
            write(iulog,*) 'll=',ll,'vol(l)=',volf(icol,ll)
            if (lwsw=='lw') then
               write(iulog,*) 'ilw=',ilwsw,' specrefindex(ilw)=',specrefindex(ll,ilwsw)
            elseif (lwsw=='sw') then
               write(iulog,*) 'isw=', ilwsw, 'specrefindex(isw)=', specrefindex(ll,ilwsw)
            endif
            write(iulog,*) 'specdens=',specdens(ll)
        enddo

        nerr_dopaer = nerr_dopaer + 1
        if (nerr_dopaer >= nerrmax_dopaer .or. dopaer < -1.e-10_r8) then
             write(iulog,*) '*** halting after nerr_dopaer =', nerr_dopaer
             call endrun()
         endif

     endif

end subroutine check_error_warning


!======================================================================
subroutine read_water_refindex(infilename)

   ! read water refractive index file and set module data

   character*(*), intent(in) :: infilename   ! modal optics filename

   ! Local variables

   integer            :: i, ierr
   type(file_desc_t)  :: ncid              ! pio file handle
   integer            :: did               ! dimension ids
   integer            :: dimlen            ! dimension lengths
   type(var_desc_t)   :: vid               ! variable ids
   real(r8) :: refrwsw(nswbands), refiwsw(nswbands) ! real, imaginary ref index for water visible
   real(r8) :: refrwlw(nlwbands), refiwlw(nlwbands) ! real, imaginary ref index for water infrared
   !----------------------------------------------------------------------------

   ! open file
   call cam_pio_openfile(ncid, infilename, PIO_NOWRITE)

   ! inquire dimensions.  Check that file values match parameter values.

   ierr = pio_inq_dimid(ncid, 'lw_band', did)
   ierr = pio_inq_dimlen(ncid, did, dimlen)
   if (dimlen .ne. nlwbands) then
      write(iulog,*) 'lw_band len=', dimlen, ' from ', infilename, ' ne nlwbands=', nlwbands
      call endrun('read_modal_optics: bad lw_band value')
   endif

   ierr = pio_inq_dimid(ncid, 'sw_band', did)
   ierr = pio_inq_dimlen(ncid, did, dimlen)
   if (dimlen .ne. nswbands) then
      write(iulog,*) 'sw_band len=', dimlen, ' from ', infilename, ' ne nswbands=', nswbands
      call endrun('read_modal_optics: bad sw_band value')
   endif

   ! read variables
   ierr = pio_inq_varid(ncid, 'refindex_real_water_sw', vid)
   ierr = pio_get_var(ncid, vid, refrwsw)

   ierr = pio_inq_varid(ncid, 'refindex_im_water_sw', vid)
   ierr = pio_get_var(ncid, vid, refiwsw)

   ierr = pio_inq_varid(ncid, 'refindex_real_water_lw', vid)
   ierr = pio_get_var(ncid, vid, refrwlw)

   ierr = pio_inq_varid(ncid, 'refindex_im_water_lw', vid)
   ierr = pio_get_var(ncid, vid, refiwlw)

   ! set complex representation of refractive indices as module data
   do i = 1, nswbands
      crefwsw(i)  = cmplx(refrwsw(i), abs(refiwsw(i)),kind=r8)
   end do
   do i = 1, nlwbands
      crefwlw(i)  = cmplx(refrwlw(i), abs(refiwlw(i)),kind=r8)
   end do

   call pio_closefile(ncid)

end subroutine read_water_refindex

!===============================================================================

subroutine modal_size_parameters(ncol, sigma_logr_aer, dgnumwet, & ! in
                                 radsurf, logradsurf, cheb,      & ! out
                                 ismethod2 ) ! optional in

   use mam_support, only : min_max_bound

   integer,  intent(in)  :: ncol
   real(r8), intent(in)  :: sigma_logr_aer  ! geometric standard deviation of number distribution
   real(r8), intent(in)  :: dgnumwet(:,:)   ! aerosol wet number mode diameter [m]
   real(r8), intent(out) :: radsurf(:,:)    ! aerosol surface mode radius [m]
   real(r8), intent(out) :: logradsurf(:,:) ! log(aerosol surface mode radius)
   real(r8), intent(out) :: cheb(:,:,:)     ! chebychev polynomial parameters

   ! FORTRAN refactoring: ismethod is tempararily used to ensure BFB test
   logical,intent(in),optional :: ismethod2

   integer  :: icol, kk, nc
   real(r8) :: alnsg_amode      ! log(sigma)
   real(r8) :: explnsigma
   real(r8) :: xrad ! normalized aerosol radius
   !-------------------------------------------------------------------------------

   alnsg_amode = log(sigma_logr_aer)
   explnsigma = exp(2.0_r8*alnsg_amode*alnsg_amode)

   do kk = top_lev, pver
      do icol = 1, ncol
         ! convert from number mode diameter to surface area
         radsurf(icol,kk) = 0.5_r8*dgnumwet(icol,kk)*explnsigma

         ! --------------- FORTRAN refactoring -------------------
         ! here two calculations are used to ensure passing BFB test
         ! can be simplified (there is only round-off difference
         if (present(ismethod2) .and. ismethod2) then
             logradsurf(icol,kk) = log(0.5_r8*dgnumwet(icol,kk)) + 2.0_r8*alnsg_amode*alnsg_amode
         else
             logradsurf(icol,kk) = log(radsurf(icol,kk))
         endif
         ! --------------- FORTRAN refactoring -------------------

         ! normalize size parameter
         xrad = min_max_bound(xrmin, xrmax, logradsurf(icol,kk))
         xrad = (2._r8*xrad-xrmax-xrmin)/(xrmax-xrmin)
         ! chebyshev polynomials
         cheb(1,icol,kk) = 1._r8
         cheb(2,icol,kk) = xrad
         do nc = 3, ncoef
            cheb(nc,icol,kk) = 2._r8*xrad*cheb(nc-1,icol,kk) - cheb(nc-2,icol,kk)
         enddo
      enddo
   enddo

end subroutine modal_size_parameters

!===============================================================================

      subroutine binterp(table,ncol,km,im,jm,x,y,xtab,ytab,ix,jy,t,u,out)

!     bilinear interpolation of table
!
      implicit none
      integer im,jm,km,ncol
      real(r8) table(km,im,jm),xtab(im),ytab(jm),out(pcols,km)
      integer i,ix(pcols),ip1,j,jy(pcols),jp1,k,ic
      real(r8) x(pcols),dx,t(pcols),y(pcols),dy,u(pcols), &
             tu(pcols),tuc(pcols),tcu(pcols),tcuc(pcols)

      if(ix(1).gt.0)go to 30
      if(im.gt.1)then
        do ic=1,ncol
          do i=1,im
            if(x(ic).lt.xtab(i))go to 10
          enddo
   10     ix(ic)=max0(i-1,1)
          ip1=min(ix(ic)+1,im)
          dx=(xtab(ip1)-xtab(ix(ic)))
          if(abs(dx).gt.1.e-20_r8)then
             t(ic)=(x(ic)-xtab(ix(ic)))/dx
          else
             t(ic)=0._r8
          endif
	end do
      else
        ix(:ncol)=1
        t(:ncol)=0._r8
      endif
      if(jm.gt.1)then
        do ic=1,ncol
          do j=1,jm
            if(y(ic).lt.ytab(j))go to 20
          enddo
   20     jy(ic)=max0(j-1,1)
          jp1=min(jy(ic)+1,jm)
          dy=(ytab(jp1)-ytab(jy(ic)))
          if(abs(dy).gt.1.e-20_r8)then
             u(ic)=(y(ic)-ytab(jy(ic)))/dy
             if(u(ic).lt.0._r8.or.u(ic).gt.1._r8)then
                write(iulog,*) 'u,y,jy,ytab,dy=',u(ic),y(ic),jy(ic),ytab(jy(ic)),dy
             endif
          else
            u(ic)=0._r8
          endif
	end do
      else
        jy(:ncol)=1
        u(:ncol)=0._r8
      endif
   30 continue
      do ic=1,ncol
         tu(ic)=t(ic)*u(ic)
         tuc(ic)=t(ic)-tu(ic)
         tcuc(ic)=1._r8-tuc(ic)-u(ic)
         tcu(ic)=u(ic)-tu(ic)
         jp1=min(jy(ic)+1,jm)
         ip1=min(ix(ic)+1,im)
         do k=1,km
            out(ic,k)=tcuc(ic)*table(k,ix(ic),jy(ic))+tuc(ic)*table(k,ip1,jy(ic))   &
               +tu(ic)*table(k,ip1,jp1)+tcu(ic)*table(k,ix(ic),jp1)
	 end do
      enddo
      return
      end subroutine binterp

end module modal_aer_opt
