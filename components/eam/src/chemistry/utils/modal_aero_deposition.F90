module modal_aero_deposition

!------------------------------------------------------------------------------------------------
! Purpose:
!
! Partition the contributions from modal components of wet and dry 
! deposition at the surface into the fields passed to the coupler.
!
! *** N.B. *** Currently only a simple scheme for the 3-mode version
!              of MAM has been implemented.
!
! Revision history:
! Feb 2009  M. Flanner, B. Eaton   Original version for trop_mam3.
! Jul 2011  F Vitt -- made avaliable to be used in a prescribed modal aerosol mode (no prognostic MAM)
! Mar 2012  F Vitt -- made changes for to prevent abort when 7-mode aeroslol model is used
!                     some of the needed consituents do not exist in 7-mode so bin_fluxes will be false
! May 2014  F Vitt -- included contributions from MAM4 aerosols and added soa_a2 to the ocphiwet fluxes
!------------------------------------------------------------------------------------------------

use shr_kind_mod,     only: r8 => shr_kind_r8
use camsrfexch,       only: cam_out_t     
use constituents,     only: cnst_get_ind
use cam_abortutils,   only: endrun

implicit none
private
save

public :: &
   modal_aero_deposition_init, &
   set_srf_drydep,             &
   set_srf_wetdep

! Private module data
integer :: idx_bc1  = -1
integer :: idx_pom1 = -1
integer :: idx_soa1 = -1
integer :: idx_soa2 = -1
integer :: idx_dst1 = -1
integer :: idx_dst3 = -1
integer :: idx_ncl3 = -1
integer :: idx_so43 = -1
integer :: idx_bc4  = -1
integer :: idx_pom4 = -1
integer :: idx_bc3  = -1
integer :: idx_pom3 = -1

logical :: bin_fluxes = .false.

logical :: initialized = .false.

!==============================================================================
contains
!==============================================================================

subroutine modal_aero_deposition_init(bc1_ndx,pom1_ndx,soa1_ndx,soa2_ndx,dst1_ndx, &
                            dst3_ndx,ncl3_ndx,so43_ndx,num3_ndx,bc4_ndx,pom4_ndx)

! set aerosol indices for re-mapping surface deposition fluxes:
! *_a1 = accumulation mode
! *_a2 = aitken mode
! *_a3 = coarse mode

   ! can be initialized with user specified indices
   ! if called from aerodep_flx module (for prescribed modal aerosol fluxes) then these indices are specified

   integer, optional, intent(in) :: bc1_ndx,pom1_ndx,soa1_ndx,soa2_ndx,dst1_ndx,dst3_ndx,ncl3_ndx,so43_ndx,num3_ndx
   integer, optional, intent(in) :: bc4_ndx,pom4_ndx

   ! if already initialized abort the run
   if (initialized) then
     call endrun('modal_aero_deposition_init is already initialized')
   endif

   if (present(bc1_ndx)) then
      idx_bc1  = bc1_ndx
   else
      call cnst_get_ind('bc_a1',  idx_bc1)
   endif
   if (present(pom1_ndx)) then
      idx_pom1 = pom1_ndx
   else
      call cnst_get_ind('pom_a1', idx_pom1)
   endif
   if (present(soa1_ndx)) then
      idx_soa1 = soa1_ndx
   else
      call cnst_get_ind('soa_a1', idx_soa1)
   endif
   if (present(soa2_ndx)) then
      idx_soa2 = soa2_ndx
   else
      call cnst_get_ind('soa_a2', idx_soa2)
   endif
   if (present(dst1_ndx)) then
      idx_dst1 = dst1_ndx
   else
      call cnst_get_ind('dst_a1', idx_dst1,abrtf=.false.)
   endif
   if (present(dst3_ndx)) then
      idx_dst3 = dst3_ndx
   else
      call cnst_get_ind('dst_a3', idx_dst3,abrtf=.false.)
   endif
   if (present(ncl3_ndx)) then
      idx_ncl3 = ncl3_ndx
   else
      call cnst_get_ind('ncl_a3', idx_ncl3,abrtf=.false.)
   endif
   if (present(so43_ndx)) then
      idx_so43 = so43_ndx
   else
      call cnst_get_ind('so4_a3', idx_so43,abrtf=.false.)
   endif
   if (present(bc4_ndx)) then
      idx_bc4 = bc4_ndx
   else
      call cnst_get_ind('bc_a4', idx_bc4,abrtf=.false.)
   endif
   if (present(pom4_ndx)) then
      idx_pom4 = pom4_ndx
   else
      call cnst_get_ind('pom_a4', idx_pom4,abrtf=.false.)   
   endif

   bin_fluxes = .true.

   ! assign additional indices for resuspended BC and POM to coarse mode:
   call cnst_get_ind('bc_a3',  idx_bc3)
   call cnst_get_ind('pom_a3', idx_pom3)

   initialized = .true.

end subroutine modal_aero_deposition_init

!==============================================================================
subroutine set_srf_wetdep(aerdepwetis, aerdepwetcw, cam_out)

! Set surface wet deposition fluxes passed to coupler.

   ! Arguments:
   real(r8), intent(in) :: aerdepwetis(:,:)  ! aerosol wet deposition (interstitial)
   real(r8), intent(in) :: aerdepwetcw(:,:)  ! aerosol wet deposition (cloud water)
   type(cam_out_t), intent(inout) :: cam_out     ! cam export state

   ! Local variables:
   integer :: icol
   integer :: ncol                      ! number of columns
   !----------------------------------------------------------------------------

   if (.not.bin_fluxes) return

   ncol = cam_out%ncol

   do icol = 1, ncol

      ! MAM4

      ! in SNICAR+MAM, bcphiwet represents BC mixed internally within
      ! hydrometeors
      cam_out%bcphiwet(icol) = -(aerdepwetcw(icol,idx_bc1)+aerdepwetcw(icol,idx_bc4))

      ! bcphidry represents BC mixed externally to hydrometeors
      cam_out%bcphidry(icol) = -(aerdepwetis(icol,idx_bc1)+aerdepwetis(icol,idx_bc4))

      ! ocphiwet represents OC mixed internally within hydrometeors
      cam_out%ocphiwet(icol) = -(aerdepwetcw(icol,idx_pom1)+aerdepwetcw(icol,idx_pom4)+ &
                              aerdepwetcw(icol,idx_soa1)+aerdepwetcw(icol,idx_soa2))

      ! ocphidry represents OC mixed externally to hydrometeors
      cam_out%ocphidry(icol) = -(aerdepwetis(icol,idx_pom1)+aerdepwetis(icol,idx_pom4)+ &
                              aerdepwetis(icol,idx_soa1)+aerdepwetis(icol,idx_soa2))

       ! add resuspended coarse-mode BC and OC
         cam_out%bcphiwet(icol) = cam_out%bcphiwet(icol) -aerdepwetcw(icol,idx_bc3)
         cam_out%bcphidry(icol) = cam_out%bcphidry(icol) -aerdepwetis(icol,idx_bc3)
         cam_out%ocphiwet(icol) = cam_out%ocphiwet(icol) -aerdepwetcw(icol,idx_pom3)
         cam_out%ocphidry(icol) = cam_out%ocphidry(icol) -aerdepwetis(icol,idx_pom3)

      ! Four dust bins in SNICAR represent dust with dry diameters of
      ! 0.1-1.0um, 1.0-2.5um, 2.5-5.0um, 5.0-10um, respectively.  Dust
      ! mass is partitioned into these bins based on global-mean size
      ! distributions of MAM7 fine dust and coarse dust shown in Table
      ! 1 of Liu et al (2012, doi:10.5194/gmd-5-709-2012).  In MAM3,
      ! accumulation-mode dust is assumed to resemble fine dust
      cam_out%dstwet1(icol) = -(0.625_r8*(aerdepwetis(icol,idx_dst1)+aerdepwetcw(icol,idx_dst1))+ &
                             0.015_r8*(aerdepwetis(icol,idx_dst3)+aerdepwetcw(icol,idx_dst3)))

      cam_out%dstwet2(icol) = -(0.345_r8*(aerdepwetis(icol,idx_dst1)+aerdepwetcw(icol,idx_dst1))+ &
                             0.252_r8*(aerdepwetis(icol,idx_dst3)+aerdepwetcw(icol,idx_dst3)))

      cam_out%dstwet3(icol) = -(0.029_r8*(aerdepwetis(icol,idx_dst1)+aerdepwetcw(icol,idx_dst1))+ &
                             0.444_r8*(aerdepwetis(icol,idx_dst3)+aerdepwetcw(icol,idx_dst3)))

      cam_out%dstwet4(icol) = -(0.001_r8*(aerdepwetis(icol,idx_dst1)+aerdepwetcw(icol,idx_dst1))+ &
                             0.289_r8*(aerdepwetis(icol,idx_dst3)+aerdepwetcw(icol,idx_dst3)))



      ! in rare cases, integrated deposition tendency is upward
      if (cam_out%bcphiwet(icol) .lt. 0._r8) cam_out%bcphiwet(icol) = 0._r8
      if (cam_out%bcphidry(icol) .lt. 0._r8) cam_out%bcphidry(icol) = 0._r8
      if (cam_out%ocphiwet(icol) .lt. 0._r8) cam_out%ocphiwet(icol) = 0._r8
      if (cam_out%ocphidry(icol) .lt. 0._r8) cam_out%ocphidry(icol) = 0._r8
      if (cam_out%dstwet1(icol)  .lt. 0._r8) cam_out%dstwet1(icol)  = 0._r8
      if (cam_out%dstwet2(icol)  .lt. 0._r8) cam_out%dstwet2(icol)  = 0._r8
      if (cam_out%dstwet3(icol)  .lt. 0._r8) cam_out%dstwet3(icol)  = 0._r8
      if (cam_out%dstwet4(icol)  .lt. 0._r8) cam_out%dstwet4(icol)  = 0._r8
   enddo

end subroutine set_srf_wetdep

!==============================================================================
subroutine set_srf_drydep(aerdepdryis, aerdepdrycw, cam_out)

! Set surface dry deposition fluxes passed to coupler.
   
   ! Arguments:
   real(r8), intent(in) :: aerdepdryis(:,:)  ! aerosol dry deposition (interstitial)
   real(r8), intent(in) :: aerdepdrycw(:,:)  ! aerosol dry deposition (cloud water)
   type(cam_out_t), intent(inout) :: cam_out     ! cam export state

   ! Local variables:
   integer :: icol
   integer :: ncol                      ! number of columns
   !----------------------------------------------------------------------------
   if (.not.bin_fluxes) return

   ncol = cam_out%ncol

   do icol = 1, ncol


      ! MAM4

      ! in SNICAR+MAM, bcphodry represents BC mixed external to hydrometeors
      cam_out%bcphodry(icol) = aerdepdryis(icol,idx_bc1)+aerdepdryis(icol,idx_bc4)+ &
                            aerdepdrycw(icol,idx_bc1)+aerdepdrycw(icol,idx_bc4)

      ! ocphodry represents OC mixed external to hydrometeors
      cam_out%ocphodry(icol) = aerdepdryis(icol,idx_pom1)+aerdepdryis(icol,idx_pom4)+aerdepdryis(icol,idx_soa1)+aerdepdryis(icol,idx_soa2)+ &
                            aerdepdrycw(icol,idx_pom1)+aerdepdrycw(icol,idx_pom4)+aerdepdrycw(icol,idx_soa1)+aerdepdrycw(icol,idx_soa2)

       ! add resuspended coarse-mode BC and OC
         cam_out%bcphodry(icol) = cam_out%bcphodry(icol) + (aerdepdryis(icol,idx_bc3)+aerdepdrycw(icol,idx_bc3))
         cam_out%ocphodry(icol) = cam_out%ocphodry(icol) + (aerdepdryis(icol,idx_pom3)+aerdepdrycw(icol,idx_pom3))

      ! NOTE: drycw fluxes shown above ideally would be included as
      ! within-hydrometeor species, but this would require passing
      ! additional species through the coupler.  drycw fluxes are
      ! extremely small in the global-mean, so this will make little
      ! difference.


      ! Four dust bins in SNICAR represent dust with dry diameters of
      ! 0.1-1.0um, 1.0-2.5um, 2.5-5.0um, 5.0-10um, respectively.  Dust
      ! mass is partitioned into these bins based on global-mean size
      ! distributions of MAM7 fine dust and coarse dust shown in Table
      ! 1 of Liu et al (2012, doi:10.5194/gmd-5-709-2012).  In MAM3,
      ! accumulation-mode dust is assumed to resemble fine dust
      cam_out%dstdry1(icol) = (0.625_r8*(aerdepdryis(icol,idx_dst1)+aerdepdrycw(icol,idx_dst1))+ &
                            0.015_r8*(aerdepdryis(icol,idx_dst3)+aerdepdrycw(icol,idx_dst3)))

      cam_out%dstdry2(icol) = (0.345_r8*(aerdepdryis(icol,idx_dst1)+aerdepdrycw(icol,idx_dst1))+ &
                            0.252_r8*(aerdepdryis(icol,idx_dst3)+aerdepdrycw(icol,idx_dst3)))

      cam_out%dstdry3(icol) = (0.029_r8*(aerdepdryis(icol,idx_dst1)+aerdepdrycw(icol,idx_dst1))+ &
                            0.444_r8*(aerdepdryis(icol,idx_dst3)+aerdepdrycw(icol,idx_dst3)))

      cam_out%dstdry4(icol) = (0.001_r8*(aerdepdryis(icol,idx_dst1)+aerdepdrycw(icol,idx_dst1))+ &
                            0.289_r8*(aerdepdryis(icol,idx_dst3)+aerdepdrycw(icol,idx_dst3)))



      ! in rare cases, integrated deposition tendency is upward
      if (cam_out%bcphidry(icol) .lt. 0._r8) cam_out%bcphidry(icol) = 0._r8
      if (cam_out%bcphodry(icol) .lt. 0._r8) cam_out%bcphodry(icol) = 0._r8
      if (cam_out%ocphidry(icol) .lt. 0._r8) cam_out%ocphidry(icol) = 0._r8
      if (cam_out%ocphodry(icol) .lt. 0._r8) cam_out%ocphodry(icol) = 0._r8
      if (cam_out%dstdry1(icol)  .lt. 0._r8) cam_out%dstdry1(icol)  = 0._r8
      if (cam_out%dstdry2(icol)  .lt. 0._r8) cam_out%dstdry2(icol)  = 0._r8
      if (cam_out%dstdry3(icol)  .lt. 0._r8) cam_out%dstdry3(icol)  = 0._r8
      if (cam_out%dstdry4(icol)  .lt. 0._r8) cam_out%dstdry4(icol)  = 0._r8
   enddo
end subroutine set_srf_drydep
!==============================================================================

end module modal_aero_deposition
