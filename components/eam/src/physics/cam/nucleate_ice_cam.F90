module nucleate_ice_cam

!---------------------------------------------------------------------------------
!
!  CAM Interfaces for nucleate_ice module.
!
!  B. Eaton - Sept 2014
!---------------------------------------------------------------------------------

use shr_kind_mod,   only: r8=>shr_kind_r8
use spmd_utils,     only: masterproc
use ppgrid,         only: pcols, pver
use physconst,      only: pi, rair, tmelt
use constituents,   only: cnst_get_ind
use physics_types,  only: physics_state
use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field
use phys_control,   only: use_hetfrz_classnuc
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_aer_mmr, rad_cnst_get_aer_props, &
                            rad_cnst_get_mode_num, rad_cnst_get_mode_props

use physics_buffer, only: pbuf_add_field, dtype_r8, pbuf_old_tim_idx, &
                          pbuf_get_index, pbuf_get_field
use cam_history,    only: addfld, add_default, outfld

use ref_pres,       only: top_lev => trop_cloud_top_lev
use wv_saturation,  only: qsat_water
#ifndef HAVE_ERF_INTRINSICS
use shr_spfn_mod,   only: erf => shr_spfn_erf
#endif
use cam_logfile,    only: iulog
use cam_abortutils, only: endrun

use nucleate_ice,   only: nucleati_init, nucleati


implicit none
private
save

public :: &
   nucleate_ice_cam_readnl,   &
   nucleate_ice_cam_register, &
   nucleate_ice_cam_init,     &
   nucleate_ice_cam_calc
   

! Namelist variables
logical, public, protected :: use_preexisting_ice = .false.
logical                    :: hist_preexisting_ice = .false.
logical, public, protected :: use_nie_nucleate = .false.
logical, public, protected :: use_dem_nucleate = .false.
real(r8)                   :: nucleate_ice_subgrid
real(r8)                   :: so4_sz_thresh_icenuc = huge(1.0_r8) !ice nucleation SO2 size threshold for aitken mode

! Vars set via init method.
real(r8) :: mincld      ! minimum allowed cloud fraction
real(r8) :: bulk_scale  ! prescribed aerosol bulk sulfur scale factor

! constituent indices
integer :: &
   cldliq_idx = -1, &
   cldice_idx = -1, &
   numice_idx = -1

integer :: &
   naai_idx,     &
   naai_hom_idx

integer :: &
   ast_idx   = -1, &
   dgnum_idx = -1


! modal aerosols
logical :: clim_modal_aero

integer :: nmodes = -1
integer :: mode_accum_idx  = -1  ! index of accumulation mode
integer :: mode_aitken_idx = -1  ! index of aitken mode
integer :: mode_coarse_idx = -1  ! index of coarse mode
integer :: mode_coarse_dst_idx = -1  ! index of coarse dust mode
integer :: mode_coarse_slt_idx = -1  ! index of coarse sea salt mode
integer :: coarse_dust_idx = -1  ! index of dust in coarse mode
integer :: coarse_nacl_idx = -1  ! index of nacl in coarse mode

integer :: coarse_so4_idx = -1  ! index of so4 in coarse mode
integer :: coarse_mom_idx = -1  ! index of mom in coarse mode
integer :: coarse_bc_idx  = -1  ! index of bc in coarse mode
integer :: coarse_pom_idx = -1  ! index of pom in coarse mode
integer :: coarse_soa_idx = -1  ! index of soa in coarse mode

integer :: mode_fine_dst_idx = -1   ! index of dust in fine dust mode
integer :: mode_pcarbon_idx  = -1  ! index of dust in accum mode
integer :: accum_dust_idx    = -1  ! index of dust in accum mode
integer :: fine_dust_idx    = -1   ! index of dust in fine mode

logical  :: separate_dust = .false.
real(r8) :: sigmag_aitken
real(r8) :: sigmag_coarse


!===============================================================================
contains
!===============================================================================

subroutine nucleate_ice_cam_readnl(nlfile)

  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use mpishorthand

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

  ! Local variables
  integer :: unitn, ierr
  character(len=*), parameter :: subname = 'nucleate_ice_cam_readnl'

  namelist /nucleate_ice_nl/ use_preexisting_ice, hist_preexisting_ice, &
                             use_nie_nucleate, use_dem_nucleate,        &
                             nucleate_ice_subgrid, so4_sz_thresh_icenuc 

  !-----------------------------------------------------------------------------

  if (masterproc) then
     unitn = getunit()
     open( unitn, file=trim(nlfile), status='old' )
     call find_group_name(unitn, 'nucleate_ice_nl', status=ierr)
     if (ierr == 0) then
        read(unitn, nucleate_ice_nl, iostat=ierr)
        if (ierr /= 0) then
           call endrun(subname // ':: ERROR reading namelist')
        end if
     end if
     close(unitn)
     call freeunit(unitn)

  end if

#ifdef SPMD
  ! Broadcast namelist variables
  call mpibcast(use_preexisting_ice,  1, mpilog, 0, mpicom)
  call mpibcast(hist_preexisting_ice, 1, mpilog, 0, mpicom)
  call mpibcast(use_nie_nucleate, 1, mpilog, 0, mpicom)
  call mpibcast(use_dem_nucleate, 1, mpilog, 0, mpicom)
  call mpibcast(nucleate_ice_subgrid, 1, mpir8, 0, mpicom)
  call mpibcast(so4_sz_thresh_icenuc, 1, mpir8, 0, mpicom)
#endif

end subroutine nucleate_ice_cam_readnl

!================================================================================================

subroutine nucleate_ice_cam_register()

   call pbuf_add_field('NAAI',     'physpkg', dtype_r8, (/pcols,pver/), naai_idx)
   call pbuf_add_field('NAAI_HOM', 'physpkg', dtype_r8, (/pcols,pver/), naai_hom_idx)

end subroutine nucleate_ice_cam_register

!================================================================================================

subroutine nucleate_ice_cam_init(mincld_in, bulk_scale_in)

   real(r8), intent(in) :: mincld_in
   real(r8), intent(in) :: bulk_scale_in

   ! local variables
   integer  :: iaer
   integer  :: m, n, nspec

   character(len=32) :: str32
   character(len=*), parameter :: routine = 'nucleate_ice_cam_init'
   !--------------------------------------------------------------------------------------------

   mincld     = mincld_in
   bulk_scale = bulk_scale_in

   call cnst_get_ind('CLDLIQ', cldliq_idx)
   call cnst_get_ind('CLDICE', cldice_idx)
   call cnst_get_ind('NUMICE', numice_idx)

   call addfld('NIHF', (/ 'lev' /), 'A',  '1/m3', 'Activated Ice Number Concentation due to homogenous freezing')
   call addfld('NIDEP', (/ 'lev' /), 'A', '1/m3', 'Activated Ice Number Concentation due to deposition nucleation')
   call addfld('NIIMM', (/ 'lev' /), 'A', '1/m3', 'Activated Ice Number Concentation due to immersion freezing')
   call addfld('NIMEY', (/ 'lev' /), 'A', '1/m3', 'Activated Ice Number Concentation due to meyers deposition')


   ! clim_modal_aero determines whether modal aerosols are used in the climate calculation.
   ! The modal aerosols can be either prognostic or prescribed.
   call rad_cnst_get_info(0, nmodes=nmodes)
   clim_modal_aero = (nmodes > 0)

   if (clim_modal_aero) then

      dgnum_idx    = pbuf_get_index('DGNUM' )

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
         end select
      end do

      ! check if coarse dust is in separate mode
      separate_dust = mode_coarse_dst_idx > 0

      ! for 3-mode 
      if (mode_coarse_dst_idx < 0) mode_coarse_dst_idx = mode_coarse_idx
      if (mode_coarse_slt_idx < 0) mode_coarse_slt_idx = mode_coarse_idx
      if (mode_fine_dst_idx < 0 .and. (use_nie_nucleate .or. use_dem_nucleate)) mode_fine_dst_idx = mode_accum_idx


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

      if (mode_coarse_idx > 0) then
         call rad_cnst_get_info(0, mode_coarse_idx, nspec=nspec)
         do n = 1, nspec
            call rad_cnst_get_info(0, mode_coarse_idx, n, spec_type=str32)
            select case (trim(str32))
            case ('sulfate')
               coarse_so4_idx = n
            end select
         end do
      end if

      ! Check that required mode specie types were found
      if (mode_coarse_idx > 0) then
         if ( coarse_so4_idx == -1) then
            write(iulog,*) routine//': ERROR required mode-species type not found - indicies:', &
               coarse_so4_idx
            call endrun(routine//': ERROR required mode-species type not found')
         end if
      end if

      call rad_cnst_get_info(0, mode_coarse_idx, nspec=nspec)
      do n = 1, nspec
         call rad_cnst_get_info(0, mode_coarse_idx, n, spec_type=str32)
         select case (trim(str32))
         case ('m-organic')
            coarse_mom_idx = n
         end select
      end do

      ! Check that required mode specie types were found
      if ( coarse_mom_idx == -1) then
         write(iulog,*) routine//': ERROR required mode-species type not found - indicies:', &
            coarse_mom_idx
         call endrun(routine//': ERROR required mode-species type not found')
      end if

      call rad_cnst_get_info(0, mode_coarse_idx, nspec=nspec)
      do n = 1, nspec
         call rad_cnst_get_info(0, mode_coarse_idx, n, spec_type=str32)
         select case (trim(str32))
         case ('black-c')
            coarse_bc_idx = n
         end select
      end do

      call rad_cnst_get_info(0, mode_coarse_idx, nspec=nspec)
      do n = 1, nspec
         call rad_cnst_get_info(0, mode_coarse_idx, n, spec_type=str32)
         select case (trim(str32))
         case ('p-organic')
            coarse_pom_idx = n
         end select
      enddo

      call rad_cnst_get_info(0, mode_coarse_idx, nspec=nspec)
      do n = 1, nspec
         call rad_cnst_get_info(0, mode_coarse_idx, n, spec_type=str32)
         select case (trim(str32))
         case ('s-organic')
            coarse_soa_idx = n
         end select
      enddo

      ! Check that required mode specie types were found
      if ( coarse_bc_idx == -1 .or. coarse_pom_idx == -1 .or. coarse_soa_idx == -1 ) then
         write(iulog,*) routine//': ERROR required mode-species type not found - indicies:', &
            coarse_bc_idx, coarse_pom_idx, coarse_soa_idx 
         call endrun(routine//': ERROR required mode-species type not found')
      endif

      ! get specific mode properties
      call rad_cnst_get_mode_props(0, mode_aitken_idx, sigmag=sigmag_aitken)

   endif

   call nucleati_init(use_preexisting_ice, use_hetfrz_classnuc,  &
                      use_nie_nucleate, use_dem_nucleate,        &
                      iulog, pi, mincld, nucleate_ice_subgrid)

   ! get indices for fields in the physics buffer
   ast_idx      = pbuf_get_index('AST')

end subroutine nucleate_ice_cam_init

!================================================================================================

subroutine nucleate_ice_cam_calc( ncol, lchnk, temperature, q, pmid, rho, wsubi, pbuf)
   use modal_aero_data,   only: modeptr_accum, modeptr_aitken, modeptr_coarse
   use modal_aero_data,   only: numptr_amode, alnsg_amode
   use modal_aero_data,   only: lptr_dust_a_amode, lptr_nacl_a_amode, lptr_so4_a_amode, &
                                lptr_bc_a_amode, lptr_pom_a_amode, lptr_soa_a_amode, lptr_mom_a_amode

   ! arguments
   integer, intent(in) :: ncol
   integer, intent(in) :: lchnk
   real(r8), intent(in) :: temperature(:,:)     ! input temperature [K]
   real(r8), intent(in) :: q(:,:,:)             ! input mixing ratio for all water and chem species [#/kg,kg/kg]
   real(r8), intent(in) :: pmid(:,:)            ! pressure at layer midpoints [pa]
   real(r8), intent(in) :: rho(:,:)             ! air density [kg/m^3]
   real(r8),                    intent(in)    :: wsubi(:,:)   ! updraft velocity for ice nucleation [m/s]  
   type(physics_buffer_desc),   pointer       :: pbuf(:)
 
   ! local workspace

   ! naai and naai_hom are the outputs shared with the microphysics
   real(r8), pointer :: naai(:,:)       ! number of activated aerosol for ice nucleation [#/kg]
   real(r8), pointer :: naai_hom(:,:)   ! number of activated aerosol for ice nucleation (homogeneous freezing only) [#/kg]

   real(r8), pointer :: dgnum(:,:,:)    ! mode dry radius [m]
   real(r8), pointer :: ast(:,:)        ! cloud fraction [unitless]


   integer :: itim_old
   integer :: i, k, m

   real(r8) :: qn(pcols,pver)           ! water vapor mixing ratio [kg/kg]
   real(r8) :: num_aitken(pcols,pver)   ! number m.r. of aitken mode [#/kg]
   real(r8) :: num_coarse(pcols,pver)   ! number m.r. of coarse mode [#/kg]
   real(r8) :: coarse_dust(pcols,pver)  ! mass m.r. of coarse dust [kg/kg]
   real(r8) :: coarse_nacl(pcols,pver)  ! mass m.r. of coarse nacl [kg/kg]

   real(r8) :: coarse_so4(pcols,pver)   ! mass m.r. of coarse so4 [kg/kg]
   real(r8) :: coarse_mom(pcols,pver)   ! mass m.r. of coarse mom [kg/kg]
   real(r8) :: coarse_bc(pcols,pver)    ! mass m.r. of coarse bc [kg/kg]
   real(r8) :: coarse_pom(pcols,pver)   ! mass m.r. of coarse pom [kg/kg] 
   real(r8) :: coarse_soa(pcols,pver)   ! mass m.r. of coarse soa [kg/kg]
   real(r8) :: icecldf(pcols,pver)      ! ice cloud fraction [unitless]

   real(r8) :: qs(pcols)                ! liquid-ice weighted sat mixing rat [kg/kg]
   real(r8) :: es(pcols)                ! liquid-ice weighted sat vapor press [pa]
   real(r8) :: gammas(pcols)            ! parameter for cond/evap of cloud water

   real(r8) :: relhum(pcols,pver)       ! relative humidity [unitless]
   real(r8) :: icldm(pcols,pver)        ! ice cloud fraction [unitless]

   real(r8) :: so4_num                  ! so4 aerosol number [#/cm^3]
   real(r8) :: dst3_num                 ! dust aerosol number [#/cm^3]
   real(r8) :: dst_num                  ! total dust aerosol number [#/cm^3]
   real(r8) :: wght                     ! mass weight [unitless]
   real(r8) :: dmc                      ! dust mass concentration [kg/m^3]
   real(r8) :: ssmc                     ! sea salt mass concentration [kg/m^3]
   real(r8) :: so4mc                    ! sulfate mass concentration [kg/m^3]
   real(r8) :: mommc                    ! marine organic mass concentration [kg/m^3]
   real(r8) :: bcmc                     ! black carbon mass concentration [kg/m^3]
   real(r8) :: pommc                    ! primary organic carbon mass concentration [kg/m^3]
   real(r8) :: soamc                    ! secondary organic carbon mass concentration [kg/m^3]

   ! history output for ice nucleation
   real(r8) :: nihf(pcols,pver)  !output number conc of ice nuclei due to heterogenous freezing [1/m3]
   real(r8) :: niimm(pcols,pver) !output number conc of ice nuclei due to immersion freezing (hetero nuc) [1/m3]
   real(r8) :: nidep(pcols,pver) !output number conc of ice nuclei due to deoposion nucleation (hetero nuc) [1/m3]
   real(r8) :: nimey(pcols,pver) !output number conc of ice nuclei due to meyers deposition [1/m3]

   real(r8), parameter :: num_m3_to_cm3 = 1.0e-6_r8

   !-------------------------------------------------------------------------------

   ! note for C++ porting
   ! read ast, dgnum, naai, naai_hom from pbuf
   ! will need to change according to how pbuf variables are stored in C++ structure
   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, ast_idx, ast, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

   icecldf(:ncol,:pver) = ast(:ncol,:pver)

   call pbuf_get_field(pbuf, dgnum_idx, dgnum)

   ! naai and naai_hom are the outputs from this parameterization
   call pbuf_get_field(pbuf, naai_idx, naai)
   call pbuf_get_field(pbuf, naai_hom_idx, naai_hom)


   qn(:ncol,:pver) = q(:ncol,:pver,1)

   ! mode number mixing ratios
   num_aitken(:ncol,:pver) = q(:ncol,:pver,numptr_amode(modeptr_aitken))
   num_coarse(:ncol,:pver) = q(:ncol,:pver,numptr_amode(modeptr_coarse))

   ! mode specie mass m.r. 
   coarse_dust(:ncol,:pver) = q(:ncol,:pver,lptr_dust_a_amode(modeptr_coarse))
   coarse_nacl(:ncol,:pver) = q(:ncol,:pver,lptr_nacl_a_amode(modeptr_coarse))
   coarse_so4(:ncol,:pver)  = q(:ncol,:pver,lptr_so4_a_amode(modeptr_coarse))
   coarse_mom(:ncol,:pver)  = q(:ncol,:pver,lptr_mom_a_amode(modeptr_coarse))
   coarse_bc(:ncol,:pver)   = q(:ncol,:pver,lptr_bc_a_amode(modeptr_coarse))
   coarse_pom(:ncol,:pver)  = q(:ncol,:pver,lptr_pom_a_amode(modeptr_coarse))
   coarse_soa(:ncol,:pver)  = q(:ncol,:pver,lptr_soa_a_amode(modeptr_coarse))


   naai(1:ncol,1:pver)     = 0._r8  
   naai_hom(1:ncol,1:pver) = 0._r8  

   ! initialize history output fields for ice nucleation
   nihf(1:ncol,1:pver)  = 0._r8  
   niimm(1:ncol,1:pver) = 0._r8  
   nidep(1:ncol,1:pver) = 0._r8 
   nimey(1:ncol,1:pver) = 0._r8 

   do k = top_lev, pver

      ! Get humidity and saturation vapor pressures
      ! This subsoutine is also used in MG cloud microphysics
      ! Probably done by SCREAM team, leave it as it is.
      call qsat_water(temperature(:ncol,k), pmid(:ncol,k), &
           es(:ncol), qs(:ncol), gam=gammas(:ncol))

      do i = 1, ncol

         relhum(i,k) = qn(i,k)/qs(i)

         ! get cloud fraction, check for minimum
         icldm(i,k) = max(icecldf(i,k), mincld)

      enddo
   enddo


   do k = top_lev, pver
      do i = 1, ncol

         if (temperature(i,k) < tmelt - 5._r8) then

            ! compute aerosol number for so4, soot, and dust with units #/cm^3
            ! remove soot number, because it is set to zero
            so4_num  = 0._r8
            dst3_num = 0._r8
            dst_num  = 0._r8

            
            !For modal aerosols, assume for the upper troposphere:
            ! soot = accumulation mode
            ! sulfate = aiken mode
            ! dust = coarse mode
            ! since modal has internal mixtures.
            dmc  = coarse_dust(i,k)*rho(i,k)      
            ssmc = coarse_nacl(i,k)*rho(i,k)
            so4mc  = coarse_so4(i,k)*rho(i,k)

            mommc  = coarse_mom(i,k)*rho(i,k)
            bcmc   = coarse_bc(i,k)*rho(i,k)
            pommc  = coarse_pom(i,k)*rho(i,k)
            soamc  = coarse_soa(i,k)*rho(i,k)

            if (dmc > 0._r8) then
               wght = dmc/(ssmc + dmc + so4mc + bcmc + pommc + soamc + mommc)
               dst3_num = wght * num_coarse(i,k)*rho(i,k)*num_m3_to_cm3
            endif

            dst_num = dst3_num

            if (dgnum(i,k,mode_aitken_idx) > 0._r8) then
               ! only allow so4 with D>0.1 um in ice nucleation
               so4_num  = num_aitken(i,k)*rho(i,k)*num_m3_to_cm3 &
                          * (0.5_r8 - 0.5_r8*erf(log(so4_sz_thresh_icenuc/dgnum(i,k,mode_aitken_idx))/  &
                          (2._r8**0.5_r8*alnsg_amode(modeptr_aitken))))
            endif
            so4_num = max(0.0_r8, so4_num)

            call nucleati( &
               wsubi(i,k), temperature(i,k), pmid(i,k), relhum(i,k), icldm(i,k),   &
               rho(i,k), so4_num, dst_num,                               &
               naai(i,k), nihf(i,k), niimm(i,k), nidep(i,k), nimey(i,k), &
               dst3_num)


            naai_hom(i,k) = nihf(i,k)

            ! output activated ice (convert from #/kg -> #/m3)
            nihf(i,k)     = nihf(i,k) *rho(i,k)
            niimm(i,k)    = niimm(i,k)*rho(i,k)
            nidep(i,k)    = nidep(i,k)*rho(i,k)
            nimey(i,k)    = nimey(i,k)*rho(i,k)

         endif
      enddo
   enddo

   call outfld('NIHF',   nihf, pcols, lchnk)
   call outfld('NIIMM', niimm, pcols, lchnk)
   call outfld('NIDEP', nidep, pcols, lchnk)
   call outfld('NIMEY', nimey, pcols, lchnk)

end subroutine nucleate_ice_cam_calc

!================================================================================================

end module nucleate_ice_cam
