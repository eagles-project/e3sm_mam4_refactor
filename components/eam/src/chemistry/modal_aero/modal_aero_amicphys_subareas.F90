module modal_aero_amicphys_subareas
!--------------------------------------------------------------------------------------
! This module contains subroutines that reconstruct cloudy and clear-air subareas 
! within a grid cell of a host model using information about grid cell mean values,
! cloud fraction, and additional assumptions about gases and aerosol distributions.
!
! History:
! All algorithms were developed by Richard (Dick) C. Easter and Steve Ghan.
! The code was separated from modal_aero_amicphys.F90, reorganized, and separated
! into smaller subroutines by Hui Wan (PNNL) in 2020-2022.
!--------------------------------------------------------------------------------------

  use shr_kind_mod, only: wp => shr_kind_r8
  use modal_aero_amicphys_control, only: ncnst=>gas_pcnst, maxsubarea

  implicit none

  public

contains

subroutine setup_subareas( cld,                                     &! in
                           nsubarea, ncldy_subarea, jclea, jcldy,   &! out
                           iscldy_subarea, afracsub, fclea, fcldy )! out
!--------------------------------------------------------------------------------------
! Purpose: Determine the number of sub-areas and their fractional areas.
!          Assign values to some bookkeeping variables.
!--------------------------------------------------------------------------------------

  implicit none

  real(wp), intent(in)  :: cld            ! cloud fraction in the grid cell [unitless]

  integer,  intent(out) :: nsubarea       ! total number of subareas to do calculations for
  integer,  intent(out) :: ncldy_subarea  ! total # of cloudy subareas
  integer,  intent(out) :: jclea, jcldy   ! indices of the clear and cloudy subareas
  logical,  intent(out) :: iscldy_subarea(maxsubarea) 
  real(wp), intent(out) :: afracsub(maxsubarea)  ! area fraction of each active subarea
  real(wp), intent(out) :: fclea, fcldy          ! area fraction of clear/cloudy subarea

  ! Cloud chemistry is only on when cld(i,k) >= 1.0e-5_wp
  ! It may be that the macrophysics has a higher threshold that this
  real(wp), parameter :: fcld_locutoff = 1.0e-5_wp

  real(wp), parameter :: fcld_hicutoff = 0.999_wp

  ! if cloud fraction ~= 0, the grid-cell has a single clear  sub-area      (nsubarea = 1)
  ! if cloud fraction ~= 1, the grid-cell has a single cloudy sub-area      (nsubarea = 1)
  ! otherwise,              the grid-cell has a clear and a cloudy sub-area (nsubarea = 2)

  if (cld < fcld_locutoff) then
     fcldy = 0.0_wp
     nsubarea = 1
     ncldy_subarea = 0
     jclea = 1
     jcldy = 0

  else if (cld > fcld_hicutoff) then
     fcldy = 1.0_wp
     nsubarea = 1
     ncldy_subarea = 1
     jclea = 0
     jcldy = 1

  else
     fcldy = cld
     nsubarea = 2
     ncldy_subarea = 1
     jclea = 1
     jcldy = 2
  end if

  fclea = 1.0_wp - fcldy

  ! Set up a logical array to indicate whether the subareas are clear or cloudy

  iscldy_subarea(:) = .false.
  if (jcldy>0) iscldy_subarea(jcldy) = .true.

  ! Save the area fractions to an array

  afracsub(:) = 0.0_wp
  if (jclea>0) afracsub(jclea) = fclea
  if (jcldy>0) afracsub(jcldy) = fcldy

end subroutine setup_subareas

subroutine set_subarea_relhum( ncldy_subarea,jclea,jcldy, &! in
                               afracsub,relhumgcm,        &! in
                               relhumsub                  )! out
!----------------------------------------------------------------------------
! Purpose: Set relative humidity in subareas.
!----------------------------------------------------------------------------

  integer,  intent(in) :: ncldy_subarea         ! # of cloudy subareas
  integer,  intent(in) :: jclea, jcldy          ! indices of clear and cloudy subareas
  real(wp), intent(in) :: afracsub(maxsubarea)  ! area fraction [unitless] of subareas
  real(wp), intent(in) :: relhumgcm             ! grid cell mean relative humidity [unitless]

  real(wp), intent(out) :: relhumsub(maxsubarea) ! relative humidity in subareas [unitless]

  real(wp) :: relhum_tmp

  if (ncldy_subarea <= 0) then
  ! Entire grid cell is cloud-free. RH in subarea = grid cell mean.

     relhumsub(:) = relhumgcm

#if ( defined( MAM_STANDALONE ) )
  else if (cldy_rh_sameas_clear > 0) then
     relhumsub(:) = relhumgcm
#endif

  else
     ! Grid cell has a cloudy subarea. Set RH in that part to 1.0.
     relhumsub(jcldy) = 1.0_wp

     ! If the grid cell also has a clear portion, back out the subarea RH from 
     ! the grid-cell mean RH and cloud fraction.

     if (jclea > 0) then
        relhum_tmp = (relhumgcm - afracsub(jcldy))/afracsub(jclea)
        relhumsub(jclea) = max( 0.0_wp, min( 1.0_wp, relhum_tmp ) )
     end if
  end if

end subroutine set_subarea_relhum

subroutine copy_cnst( q_in, q_copy, lcopy )
!------------------------------------------------------
! Purpose: copy values from one array to another
!------------------------------------------------------
  logical, intent(in)    :: lcopy(ncnst)  ! whether individual constituents should be copied
  real(wp),intent(in)    :: q_in(ncnst)   ! values to be copied
  real(wp),intent(inout) :: q_copy(ncnst) ! copy of input values 

  where( lcopy )
    q_copy = q_in
  end where

end subroutine copy_cnst

subroutine compute_qsub_from_gbm_and_qsub_of_other_subarea( lcompute, f_a, f_b, qgcm, &! in
                                                            qsub_a, qsub_b            )! inout
!-----------------------------------------------------------------------------------------
! Purpose: Calculate the value of qsub_b assuming qgcm is a weighted average defined as
!          qgcm = f_a*qsub_a + f_b*qsub_b.
!-----------------------------------------------------------------------------------------

  logical, intent(in)    :: lcompute(ncnst)
  real(wp),intent(in)    :: f_a, f_b          ! area fractions [unitless] of subareas
  real(wp),intent(in)    :: qgcm  (ncnst)     ! grid cell mean (known)
  real(wp),intent(inout) :: qsub_a(ncnst)     ! value in subarea A (known, but might get adjusted)
  real(wp),intent(inout) :: qsub_b(ncnst)     ! value in subarea B (to be calculated here)

  integer :: icnst    ! consitituent index

  do icnst = 1, ncnst
   if (lcompute(icnst)) then

     ! Calculate qsub_b

     qsub_b(icnst) = ( qgcm(icnst) - f_a*qsub_a(icnst) )/f_b

     ! Check that this does not produce a negative value.
     ! If so, set qsub_b to zero and adjust the value of qsub_a.

     if (qsub_b(icnst) < 0.0_wp) then
        qsub_b(icnst) = 0.0_wp
        qsub_a(icnst) = qgcm(icnst)/f_a
     end if

   end if
  end do

end subroutine compute_qsub_from_gbm_and_qsub_of_other_subarea

subroutine set_subarea_q_numb_for_cldbrn_aerosols( loffset, jclea, jcldy, fcldy, qqcwgcm, qqcwsub )
!-----------------------------------------------------------------------------------------
! Purpose: Set the number mixing ratios of cloud-borne aerosols in subareas:
!          - zero in clear air;
!          - grid-cell-mean divided by cloud-fraction in the cloudy subarea.
!          This is done for all lognormal modes.
!-----------------------------------------------------------------------------------------

   use modal_aero_data, only: ntot_amode, nspec_amode, numptrcw_amode

   integer, intent(in)    :: loffset 
   integer, intent(in)    :: jclea, jcldy               ! indices of subareas
   real(wp),intent(in)    :: fcldy                      ! area fraction [unitless] of the cloudy subarea
   real(wp),intent(in)    :: qqcwgcm(ncnst)             ! grid cell mean (unit does not matter for this subr.)
   real(wp),intent(inout) :: qqcwsub(ncnst,maxsubarea)  ! values in subareas (unit does not matter for this subr.)

   integer :: imode    ! mode index
   integer :: icnst    ! consitituent index 

   !----------------------------------------------------------------
   do imode = 1, ntot_amode

      icnst = numptrcw_amode(imode) - loffset

      qqcwsub(icnst,jclea) = 0.0_wp
      qqcwsub(icnst,jcldy) = qqcwgcm(icnst)/fcldy

   end do
   !----------------------------------------------------------------

end subroutine set_subarea_q_numb_for_cldbrn_aerosols

subroutine set_subarea_q_mass_for_cldbrn_aerosols( loffset, jclea, jcldy, fcldy, qqcwgcm, qqcwsub )
!-----------------------------------------------------------------------------------------
! Purpose: Set the mass mixing ratios of cloud-borne aerosols in subareas:
!          - zero in clear air;
!          - grid-cell-mean/cloud-fraction in the cloudy subarea.
!          This is done for all lognormal modes and all chemical species.
!-----------------------------------------------------------------------------------------

   use modal_aero_data, only: ntot_amode, nspec_amode, lmassptrcw_amode

   integer, intent(in)    :: loffset 
   integer, intent(in)    :: jclea, jcldy              ! subarea indices 
   real(wp),intent(in)    :: fcldy                     ! area fraction [unitless] of the cloudy subarea 
   real(wp),intent(in)    :: qqcwgcm(ncnst)            ! grid cell mean (unit does not matter for this subr.)
   real(wp),intent(inout) :: qqcwsub(ncnst,maxsubarea) ! values in subareas (unit does not matter for this subr.)

   integer :: imode    ! mode index
   integer :: ispec    ! aerosol species index
   integer :: icnst    ! consitituent index 

   !----------------------------------------------------------------
   do imode = 1, ntot_amode          ! loop thru all modes
    do ispec = 1, nspec_amode(imode) ! loop thru all species in a mode

         icnst = lmassptrcw_amode(ispec,imode) - loffset

         qqcwsub(icnst,jclea) = 0.0_wp
         qqcwsub(icnst,jcldy) = qqcwgcm(icnst)/fcldy

    end do ! ispec - species loop
   end do ! imode - mode loop
   !----------------------------------------------------------------

end subroutine set_subarea_q_mass_for_cldbrn_aerosols

subroutine set_subarea_q_numb_for_intrst_aerosols( loffset, jclea, jcldy, fclea, fcldy, &! in
                                                   qgcm, qqcwgcm,                       &! in
                                                   qgcmx,qsubx                          )! inout
!-----------------------------------------------------------------------------------------
! Purpose: Set the number mixing ratios of interstitial aerosols in subareas.
!          Interstitial aerosols can exist in both cloudy and clear subareas, so a
!          grid cell mean needs to be partitioned. Different lognormal modes are
!          partitioned differently based on the mode-specific number mixing ratios.
!-----------------------------------------------------------------------------------------

   use modal_aero_data, only: ntot_amode, numptr_amode, numptrcw_amode

   integer, intent(in)    :: loffset 
   integer, intent(in)    :: jclea, jcldy    ! subarea indices 
   real(wp),intent(in)    :: fclea, fcldy    ! area fraction [unitless] of the clear and cloudy subareas
   real(wp),intent(in)    :: qgcm   (ncnst) ! grid cell mean of interstitial aerosol mixing ratio [unit does not matter)
   real(wp),intent(in)    :: qqcwgcm(ncnst) ! grid cell mean of cloud-borne aerosol mixing ratio [unit does not matter)
   real(wp),intent(in)    :: qgcmx  (ncnst)

   real(wp),intent(inout) :: qsubx(ncnst,maxsubarea)

   integer :: imode    ! mode index
   integer :: icnst    ! consitituent index 

   real(wp) :: tmp_qa_gcav
   real(wp) :: tmp_qc_gcav

   real(wp) :: tmp_aa_clea
   real(wp) :: tmp_aa_cldy

   do imode = 1, ntot_amode

      ! calculate partitioning factors

      tmp_qa_gcav = qgcm   ( numptr_amode(imode)-loffset )
      tmp_qc_gcav = qqcwgcm( numptrcw_amode(imode)-loffset )

      call get_partition_factors( tmp_qa_gcav, tmp_qc_gcav, fcldy, fclea, &
                                  tmp_aa_clea, tmp_aa_cldy )

      ! apply partitioning factors

      icnst = numptr_amode(imode) - loffset

      qsubx(icnst,jclea) = qgcmx(icnst)*tmp_aa_clea
      qsubx(icnst,jcldy) = qgcmx(icnst)*tmp_aa_cldy

   end do ! imode

end subroutine set_subarea_q_numb_for_intrst_aerosols

subroutine set_subarea_q_mass_for_intrst_aerosols( loffset, jclea, jcldy, fclea, fcldy, &
                                                   qgcm,  qqcwgcm, &
                                                   qgcmx, qsubx    )! inout
!-----------------------------------------------------------------------------------------
! Purpose: Set the mass mixing ratios of interstitial aerosols in subareas.
!          Interstitial aerosols can exist in both cloudy and clear subareas, so a
!          grid cell mean needs to be partitioned. Different lognormal modes are
!          partitioned differently based on the mode-specific mixing ratios.
!          All species in the same mode are partitioned the same way, consistent
!          with the internal mixing assumption used in MAM.
!-----------------------------------------------------------------------------------------

  use modal_aero_data, only: ntot_amode, nspec_amode, lmassptr_amode, lmassptrcw_amode

  integer, intent(in)    :: loffset 
  integer, intent(in)    :: jclea, jcldy
  real(wp),intent(in)    :: fclea, fcldy
  real(wp),intent(in)    :: qgcm   (ncnst)
  real(wp),intent(in)    :: qqcwgcm(ncnst)

  real(wp),intent(in)    :: qgcmx(ncnst)
  real(wp),intent(inout) :: qsubx(ncnst,maxsubarea)

  integer :: imode    ! mode index
  integer :: ispec    ! aerosol species index
  integer :: icnst    ! consitituent index 

  real(wp) :: tmp_qa_gcav
  real(wp) :: tmp_qc_gcav

  real(wp) :: tmp_aa_clea
  real(wp) :: tmp_aa_cldy

  do imode = 1, ntot_amode

     ! calculcate partitioning factors

     tmp_qa_gcav = 0.0_wp
     tmp_qc_gcav = 0.0_wp

     do ispec = 1, nspec_amode(imode)
        tmp_qa_gcav = tmp_qa_gcav + qgcm( lmassptr_amode(ispec,imode) - loffset )
        tmp_qc_gcav = tmp_qc_gcav + qqcwgcm( lmassptrcw_amode(ispec,imode) - loffset )
     end do

     call get_partition_factors( tmp_qa_gcav, tmp_qc_gcav, fcldy, fclea, &
                                 tmp_aa_clea, tmp_aa_cldy )

     ! apply partitioning factors

     do ispec = 1, nspec_amode(imode)

        icnst = lmassptr_amode(ispec,imode) - loffset

        qsubx(icnst,jclea) = qgcmx(icnst)*tmp_aa_clea
        qsubx(icnst,jcldy) = qgcmx(icnst)*tmp_aa_cldy

     end do ! ispec 
  end do ! imode


end subroutine set_subarea_q_mass_for_intrst_aerosols

subroutine get_partition_factors(  q_intrst_gcm, q_cldbrn_gcm, fcldy, fclea, &! in
                                   part_fac_q_intrst_clea, part_fac_q_intrst_cldy  )! out
!------------------------------------------------------------------------------------
! Purpose: Calculate the partitioning factors for distributing interstitial aerosol
!          mixing ratios to cloudy and clear subareas in a grid box.
!          The partitioning factors depend on the grid cell mean mixing ratios of
!          both interstitial and cloud-borne aerosols.
!------------------------------------------------------------------------------------

  real(wp), intent(in)  ::  q_intrst_gcm  ! grid cell mean interstitial aerosol mixing ratio
  real(wp), intent(in)  ::  q_cldbrn_gcm  ! grid cell mean cloud-borne aerosol mixing ratio

  real(wp), intent(in)  ::  fcldy           ! cloudy fraction of the grid cell
  real(wp), intent(in)  ::  fclea           ! clear  fraction of the grid cell

  real(wp), intent(out) ::  part_fac_q_intrst_clea
  real(wp), intent(out) ::  part_fac_q_intrst_cldy

  real(wp) :: tmp_q_intrst_clea, tmp_q_intrst_cldy
  real(wp) :: tmp_q_cldbrn_cldy
  real(wp) :: tmp_aa

  real(wp),parameter :: eps = 1.e-35_wp

  ! Calculate mixing ratios of each subarea

  tmp_q_cldbrn_cldy = q_cldbrn_gcm/fcldy ! cloud-borne,  cloudy subarea
  tmp_q_intrst_cldy = max( 0.0_wp, ((q_intrst_gcm+q_cldbrn_gcm) - tmp_q_cldbrn_cldy) ) ! interstitial, cloudy subarea

  tmp_q_intrst_clea = (q_intrst_gcm - fcldy*tmp_q_intrst_cldy)/fclea ! interstitial, clear  subarea

  ! Calculate the corresponding paritioning factors for interstitial aerosols
  ! using the above-derived subarea mixing ratios plus the constraint that
  ! the cloud fraction weighted average of subarea mean need to match grid box mean.

  tmp_aa = max( eps, tmp_q_intrst_clea*fclea ) / max( eps, q_intrst_gcm )
  tmp_aa = max( 0.0_wp, min( 1.0_wp, tmp_aa ) )

  part_fac_q_intrst_clea = tmp_aa/fclea
  part_fac_q_intrst_cldy = (1.0_wp-tmp_aa)/fcldy

end subroutine get_partition_factors


end module modal_aero_amicphys_subareas
