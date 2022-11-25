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

  real(wp), intent(in)  :: cld                        ! cloud fraction in the grid cell [unitless]

  integer,  intent(out) :: nsubarea                   ! total # of subareas to do calculations for
  integer,  intent(out) :: ncldy_subarea              ! total # of cloudy subareas
  integer,  intent(out) :: jclea, jcldy               ! indices of the clear and cloudy subareas
  logical,  intent(out) :: iscldy_subarea(maxsubarea) ! whether a subarea is cloudy
  real(wp), intent(out) :: afracsub(maxsubarea)       ! area fraction of each active subarea [unitless]
  real(wp), intent(out) :: fclea, fcldy               ! area fraction of clear/cloudy subarea [unitless]

  ! Cloud chemistry is only active when cld(i,k) >= 1.0e-5_wp
  ! It may be that the macrophysics has a higher threshold that this
  real(wp), parameter :: fcld_locutoff = 1.0e-5_wp

  ! Grid cells with cloud fraction larger than this cutoff is considered to be overcast
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

subroutine set_subarea_rh( ncldy_subarea,jclea,jcldy,afracsub,relhumgcm, &! in
                           relhumsub                                     )! out
!----------------------------------------------------------------------------
! Purpose: Set relative humidity in subareas.
!----------------------------------------------------------------------------

  integer,  intent(in) :: ncldy_subarea         ! # of cloudy subareas
  integer,  intent(in) :: jclea, jcldy          ! indices of clear and cloudy subareas
  real(wp), intent(in) :: afracsub(maxsubarea)  ! area fraction in subareas [unitless]
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

end subroutine set_subarea_rh

subroutine set_subarea_gases_and_aerosols( loffset, nsubarea, jclea, jcldy, fclea, fcldy, &! in
                                           qgcm1, qgcm2, qqcwgcm2, qgcm3, qqcwgcm3,       &! in
                                           qsub1, qsub2, qqcwsub2, qsub3, qqcwsub3        )! out
!------------------------------------------------------------------------------------------------
! Purpose: Partition grid cell mean mixing ratios to clear/cloudy subareas.
!------------------------------------------------------------------------------------------------

  use cam_abortutils,              only: endrun
  use modal_aero_amicphys_control, only: lmapcc_all, lmapcc_val_gas

  integer, intent(in) :: loffset         ! # of tracers in the host model that are not part of MAM
  integer, intent(in) :: nsubarea        ! # of active subareas in the current grid cell
  integer, intent(in) :: jclea, jcldy    ! indices of the clear and cloudy subareas
  real(wp),intent(in) :: fclea, fcldy    ! area fractions of the clear and cloudy subareas [unitless]

  ! The next set of argument variables are tracer mixing ratios.
  !  - The units are different for gases, aerosol number, and aerosol mass. The exact units do not
  !    matter for this subroutine, as long as the grid cell mean values ("gcm") and the corresponding 
  !    subarea values ("sub") have the same units.
  !  - q* and qqcw* are correspond to the interstitial and cloud-borne species, respectively
  !  - The numbers 1-3 correspond to different locations in the host model's time integration loop.

  ! Grid cell mean mixing ratios

  real(wp),intent(in) :: qgcm1(ncnst)

  real(wp),intent(in) :: qgcm2(ncnst)
  real(wp),intent(in) :: qqcwgcm2(ncnst)

  real(wp),intent(in) :: qgcm3(ncnst)
  real(wp),intent(in) :: qqcwgcm3(ncnst)

  ! Subarea mixing ratios

  real(wp),intent(out) :: qsub1(ncnst,maxsubarea)
  real(wp),intent(out) :: qsub2(ncnst,maxsubarea)
  real(wp),intent(out) :: qsub3(ncnst,maxsubarea)
  real(wp),intent(out) :: qqcwsub2(ncnst,maxsubarea)
  real(wp),intent(out) :: qqcwsub3(ncnst,maxsubarea)
  !----

  logical :: grid_cell_has_only_clea_area
  logical :: grid_cell_has_only_cldy_area
  logical :: gird_cell_is_partly_cldy

  logical :: lcopy(ncnst)
  logical :: cnst_is_gas(ncnst)

  integer :: imode, ispec, icnst, jsub

  character(len=200) :: tmp_str

  !------------------------------------------------------------------------------------
  ! Initialize mixing ratios in subareas before the aerosol microphysics calculations
  !------------------------------------------------------------------------------------
  ! Gases and interstitial aerosols
  qsub1(:,:) = 0.0_wp
  qsub2(:,:) = 0.0_wp
  qsub3(:,:) = 0.0_wp

  ! Cloud-borne aerosols
  qqcwsub2(:,:) = 0.0_wp
  qqcwsub3(:,:) = 0.0_wp

  !---------------------------------------------------------------------------------------------------
  ! Determine which category the current grid cell belongs to: partly cloudy, all cloudy, or all clear
  !---------------------------------------------------------------------------------------------------
  grid_cell_has_only_clea_area = ((jclea == 1) .and. (jcldy == 0) .and. (nsubarea == 1))
  grid_cell_has_only_cldy_area = ((jclea == 0) .and. (jcldy == 1) .and. (nsubarea == 1))
  gird_cell_is_partly_cldy = (jclea > 0) .and. (jcldy > 0) .and.  (jclea+jcldy == 3) .and. (nsubarea == 2)

  ! Sanity check
  if ( (.not.grid_cell_has_only_clea_area) .and. &
       (.not.grid_cell_has_only_cldy_area) .and. &
       (.not.gird_cell_is_partly_cldy)           ) then

     write(tmp_str,'(a,3(1x,i10))') '*** modal_aero_amicphys - bad jclea, jcldy, nsubarea', jclea, jcldy, nsubarea
     call endrun( tmp_str )

  end if

  !*************************************************************************************************
  ! Category I: grid cell is either all clear or all cloudy. Copy the grid cell mean values.
  !*************************************************************************************************
  if (grid_cell_has_only_clea_area.or.grid_cell_has_only_cldy_area) then

     lcopy(1:ncnst) = lmapcc_all(1:ncnst) > 0     ! copy all gases and aerosols

     do jsub = 1,nsubarea
        call copy_cnst( qgcm1, qsub1(:,jsub), lcopy ) !from, to, flag
        call copy_cnst( qgcm2, qsub2(:,jsub), lcopy ) !from, to, flag
        call copy_cnst( qgcm3, qsub3(:,jsub), lcopy ) !from, to, flag

        call copy_cnst( qqcwgcm2, qqcwsub2(:,jsub), lcopy ) !from, to, flag
        call copy_cnst( qqcwgcm3, qqcwsub3(:,jsub), lcopy ) !from, to, flag
     end do

  !*************************************************************************************************
  ! Category II: partly cloudy grid cell. Tracer mixing ratios are generally assumed different
  ! in clear and cloudy subareas.  This is primarily because the
  ! interstitial aerosol mixing ratios are assumed to be lower in the cloudy sub-area than in
  ! the clear sub-area, as much of the aerosol is activated in the cloudy sub-area.
  !*************************************************************************************************
  else if ( gird_cell_is_partly_cldy ) then

     !===================================
     ! Set gas mixing ratios in subareas
     !===================================
     cnst_is_gas(:) = lmapcc_all(:).eq.lmapcc_val_gas 

     !------------------------------------------------------------------------------------------
     ! Before gas chemistry, gas mixing ratios are assumed to be the same in all subareas,
     ! i.e., they all equal the grid cell mean.
     !------------------------------------------------------------------------------------------
     do jsub = 1,nsubarea
        call copy_cnst( qgcm1, qsub1(:,jsub), cnst_is_gas ) !from, to, flag
     end do

     !------------------------------------------------------------------------------------------
     ! After gas chemistry, still assume gas mixing ratios are the same in all subareas.
     !------------------------------------------------------------------------------------------
     do jsub = 1,nsubarea
        call copy_cnst( qgcm2, qsub2(:,jsub), cnst_is_gas ) !from, to, flag
     end do

     !----------------------------------------------------------------------------------------
     ! After cloud chemistry, gas and aerosol mass mixing ratios in the clear subarea are 
     ! assumed to be the same as their values before cloud chemistry  (because by definition,
     ! cloud chemistry did not happen in clear sky), while the mixing ratios in the cloudy 
     ! subarea likely have changed.
     !----------------------------------------------------------------------------------------
     ! Gases in the clear subarea remain the same as their values before cloud chemistry.

     call copy_cnst( qsub2(:,jclea), qsub3(:,jclea), cnst_is_gas ) !from, to, flag

     ! Calculater the gas mixing ratios in the cloudy subarea using the grid-cell mean, 
     ! cloud fraction and the clear-sky values

     call compute_qsub_from_gcm_and_qsub_of_other_subarea( cnst_is_gas, fclea, fcldy, qgcm3, &! in
                                                           qsub3(:,jclea), qsub3(:,jcldy)    )! inout

     !=========================================================================
     ! Set AEROSOL mixing ratios in subareas.
     ! Only need to do this for points 2 and 3 in the time integraion loop,
     ! i.e., the before-cloud-chem and after-cloud-chem states.
     !=========================================================================
     ! Cloud-borne aerosols. (They are straightforward to partition, 
     ! as they only exist in the cloudy subarea.)
     !----------------------------------------------------------------------------------------
     ! Partition mass and number before cloud chemistry

     call set_subarea_qnumb_for_cldbrn_aerosols( loffset, jclea, jcldy, fcldy, qqcwgcm2, qqcwsub2 )
     call set_subarea_qmass_for_cldbrn_aerosols( loffset, jclea, jcldy, fcldy, qqcwgcm2, qqcwsub2 )

     ! Partition mass and number before cloud chemistry

     call set_subarea_qnumb_for_cldbrn_aerosols( loffset, jclea, jcldy, fcldy, qqcwgcm3, qqcwsub3 )
     call set_subarea_qmass_for_cldbrn_aerosols( loffset, jclea, jcldy, fcldy, qqcwgcm3, qqcwsub3 )

     !----------------------------------------------------------------------------------------
     ! Interstitial aerosols. (They can exist in both cloudy and clear subareas, and hence 
     ! need to be partitioned.)
     !----------------------------------------------------------------------------------------
     ! Partition mass and number before cloud chemistry

     call set_subarea_qnumb_for_intrst_aerosols( loffset, jclea, jcldy, fclea, fcldy, &! in
                                                 qgcm2, qqcwgcm2,                     &! in
                                                 qgcm2, qsub2                         )! in, inout

     call set_subarea_qmass_for_intrst_aerosols( loffset, jclea, jcldy, fclea, fcldy, &! in
                                                 qgcm2, qqcwgcm2,                     &! in
                                                 qgcm2, qsub2                         )! in, inout

     ! Partition mass and number before cloud chemistry

     call set_subarea_qnumb_for_intrst_aerosols( loffset, jclea, jcldy, fclea, fcldy, &! in
                                                 qgcm2, qqcwgcm2,                     &! in
                                                 qgcm3, qsub3                         )! in, inout

     call set_subarea_qmass_for_intrst_aerosols( loffset, jclea, jcldy, fclea, fcldy, &! in
                                                 qgcm2, qqcwgcm2,                     &! in
                                                 qgcm3, qsub3                         )! in, inout

  end if ! different categories

end subroutine set_subarea_gases_and_aerosols

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

subroutine compute_qsub_from_gcm_and_qsub_of_other_subarea( lcompute, f_a, f_b, qgcm, &! in
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

end subroutine compute_qsub_from_gcm_and_qsub_of_other_subarea

subroutine set_subarea_qnumb_for_cldbrn_aerosols( loffset, jclea, jcldy, fcldy, qqcwgcm, qqcwsub )
!-----------------------------------------------------------------------------------------
! Purpose: Set the number mixing ratios of cloud-borne aerosols in subareas:
!          - zero in clear air;
!          - grid-cell-mean divided by cloud-fraction in the cloudy subarea.
!          This is done for all lognormal modes.
!-----------------------------------------------------------------------------------------

   use modal_aero_data, only: ntot_amode, nspec_amode, numptrcw_amode

   integer, intent(in)    :: loffset                   ! # of tracers in the host model that are not part of MAM
   integer, intent(in)    :: jclea, jcldy              ! indices of subareas
   real(wp),intent(in)    :: fcldy                     ! area fraction [unitless] of the cloudy subarea
   real(wp),intent(in)    :: qqcwgcm(ncnst)            ! grid cell mean (unit does not matter for this subr.)
   real(wp),intent(inout) :: qqcwsub(ncnst,maxsubarea) ! values in subareas (unit does not matter for this subr.)

   integer :: imode    ! mode index
   integer :: icnst    ! consitituent index 

   !----------------------------------------------------------------
   do imode = 1, ntot_amode

      icnst = numptrcw_amode(imode) - loffset

      qqcwsub(icnst,jclea) = 0.0_wp
      qqcwsub(icnst,jcldy) = qqcwgcm(icnst)/fcldy

   end do
   !----------------------------------------------------------------

end subroutine set_subarea_qnumb_for_cldbrn_aerosols

subroutine set_subarea_qmass_for_cldbrn_aerosols( loffset, jclea, jcldy, fcldy, qqcwgcm, qqcwsub )
!-----------------------------------------------------------------------------------------
! Purpose: Set the mass mixing ratios of cloud-borne aerosols in subareas:
!          - zero in clear air;
!          - grid-cell-mean/cloud-fraction in the cloudy subarea.
!          This is done for all lognormal modes and all chemical species.
!-----------------------------------------------------------------------------------------

   use modal_aero_data, only: ntot_amode, nspec_amode, lmassptrcw_amode

   integer, intent(in)    :: loffset                   ! # of tracers in the host model that are not part of MAM
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

end subroutine set_subarea_qmass_for_cldbrn_aerosols

subroutine set_subarea_qnumb_for_intrst_aerosols( loffset, jclea, jcldy, fclea, fcldy, &! in
                                                  qgcm, qqcwgcm,                       &! in
                                                  qgcmx,qsubx                          )! inout
!-----------------------------------------------------------------------------------------
! Purpose: Set the number mixing ratios of interstitial aerosols in subareas.
!          Interstitial aerosols can exist in both cloudy and clear subareas, so a
!          grid cell mean needs to be partitioned. Different lognormal modes are
!          partitioned differently based on the mode-specific number mixing ratios.
!-----------------------------------------------------------------------------------------
  use modal_aero_data, only: ntot_amode, numptr_amode, numptrcw_amode

  integer, intent(in) :: loffset        ! # of tracers in the host model that are not part of MAM
  integer, intent(in) :: jclea, jcldy   ! subarea indices 
  real(wp),intent(in) :: fclea, fcldy   ! area fraction [unitless] of the clear and cloudy subareas
  real(wp),intent(in) :: qgcm   (ncnst) ! grid cell mean, interstitial constituents (unit does not matter)
  real(wp),intent(in) :: qqcwgcm(ncnst) ! grid cell mean, cloud-borne  constituents (unit does not matter)

  real(wp),intent(in) :: qgcmx  (ncnst) ! grid cell mean, interstitial constituents (unit does not matter)
  real(wp),intent(inout) :: qsubx(ncnst,maxsubarea)  ! subarea mixing ratios of interst. constituents 
                                                     ! (unit does not matter as long as they are consistent
                                                     ! with the grid cell mean values)

  ! Note: qgcm and qqcwgcm are used for calculating the patitioning factors.
  ! qgcmx is the actual grid cell mean that is partitioned into qsubx.

  integer :: imode    ! mode index
  integer :: icnst    ! consitituent index 

  real(wp) :: qgcm_intrst  ! grid cell mean of interstitial aerosol mixing ratio of a single mode
  real(wp) :: qgcm_cldbrn  ! grid cell mean of cloud-borne  aerosol mixing ratio of a single mode

  real(wp) :: factor_clea  ! partitioning factor for clear  subarea [unitless]
  real(wp) :: factor_cldy  ! partitioning factor for cloudy subarea [unitless]

  do imode = 1, ntot_amode

     ! calculate partitioning factors

     qgcm_intrst = qgcm   ( numptr_amode(imode)-loffset )
     qgcm_cldbrn = qqcwgcm( numptrcw_amode(imode)-loffset )

     call get_partition_factors( qgcm_intrst, qgcm_cldbrn, fcldy, fclea, &
                                 factor_clea, factor_cldy )

     ! apply partitioning factors

     icnst = numptr_amode(imode) - loffset

     qsubx(icnst,jclea) = qgcmx(icnst)*factor_clea
     qsubx(icnst,jcldy) = qgcmx(icnst)*factor_cldy

  end do ! imode

end subroutine set_subarea_qnumb_for_intrst_aerosols

subroutine set_subarea_qmass_for_intrst_aerosols( loffset, jclea, jcldy, fclea, fcldy, &
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

  integer, intent(in) :: loffset        ! # of tracers in the host model that are not part of MAM
  integer, intent(in) :: jclea, jcldy   ! subarea indices 
  real(wp),intent(in) :: fclea, fcldy   ! area fraction [unitless] of the clear and cloudy subareas
  real(wp),intent(in) :: qgcm   (ncnst) ! grid cell mean, interstitial constituents (unit does not matter)
  real(wp),intent(in) :: qqcwgcm(ncnst) ! grid cell mean, cloud-borne  constituents (unit does not matter)

  real(wp),intent(in) :: qgcmx  (ncnst) ! grid cell mean, interstitial constituents (unit does not matter)
  real(wp),intent(inout) :: qsubx(ncnst,maxsubarea)  ! subarea mixing ratios of interst. constituents 
                                                     ! (unit does not matter as long as they are consistent
                                                     ! with the grid cell mean values)

  ! Note: qgcm and qqcwgcm are used for calculating the patitioning factors.
  ! qgcmx is the actual grid cell mean that is partitioned into qsubx.

  integer :: imode    ! mode index
  integer :: ispec    ! species index
  integer :: icnst    ! consitituent index 

  real(wp) :: qgcm_intrst  ! grid cell mean of interstitial aerosol mixing ratio of a single mode
  real(wp) :: qgcm_cldbrn  ! grid cell mean of cloud-borne  aerosol mixing ratio of a single mode

  real(wp) :: factor_clea  ! partitioning factor for clear  subarea [unitless]
  real(wp) :: factor_cldy  ! partitioning factor for cloudy subarea [unitless]

  do imode = 1, ntot_amode

     ! calculcate partitioning factors

     qgcm_intrst = 0.0_wp
     qgcm_cldbrn = 0.0_wp

     do ispec = 1, nspec_amode(imode)
        qgcm_intrst = qgcm_intrst + qgcm( lmassptr_amode(ispec,imode) - loffset )
        qgcm_cldbrn = qgcm_cldbrn + qqcwgcm( lmassptrcw_amode(ispec,imode) - loffset )
     end do

     call get_partition_factors( qgcm_intrst, qgcm_cldbrn, fcldy, fclea, &
                                 factor_clea, factor_cldy )

     ! apply partitioning factors

     do ispec = 1, nspec_amode(imode)

        icnst = lmassptr_amode(ispec,imode) - loffset

        qsubx(icnst,jclea) = qgcmx(icnst)*factor_clea
        qsubx(icnst,jcldy) = qgcmx(icnst)*factor_cldy

     end do ! ispec 
  end do ! imode


end subroutine set_subarea_qmass_for_intrst_aerosols

subroutine get_partition_factors(  qgcm_intrst, qgcm_cldbrn, fcldy, fclea, &! in
                                   factor_clea, factor_cldy                )! out
!------------------------------------------------------------------------------------
! Purpose: Calculate the partitioning factors for distributing interstitial aerosol
!          mixing ratios to cloudy and clear subareas in a grid box.
!          The partitioning factors depend on the grid cell mean mixing ratios of
!          both interstitial and cloud-borne aerosols.
!------------------------------------------------------------------------------------

  real(wp), intent(in)  ::  qgcm_intrst  ! grid cell mean interstitial aerosol mixing ratio
  real(wp), intent(in)  ::  qgcm_cldbrn  ! grid cell mean cloud-borne aerosol mixing ratio

  real(wp), intent(in)  ::  fcldy        ! cloudy fraction of the grid cell [unitless]
  real(wp), intent(in)  ::  fclea        ! clear  fraction of the grid cell [unitless]

  real(wp), intent(out) ::  factor_clea  ! partitioning factor for clear  subarea
  real(wp), intent(out) ::  factor_cldy  ! partitioning factor for cloudy subarea

  real(wp) :: tmp_q_intrst_clea, tmp_q_intrst_cldy
  real(wp) :: tmp_q_cldbrn_cldy
  real(wp) :: clea2gcm_ratio

  real(wp),parameter :: eps = 1.e-35_wp

  ! Calculate subarea-mean mixing ratios

  tmp_q_cldbrn_cldy = qgcm_cldbrn/fcldy                                              ! cloud-borne,  cloudy subarea
  tmp_q_intrst_cldy = max( 0.0_wp, ((qgcm_intrst+qgcm_cldbrn) - tmp_q_cldbrn_cldy) ) ! interstitial, cloudy subarea

  tmp_q_intrst_clea = (qgcm_intrst - fcldy*tmp_q_intrst_cldy)/fclea                  ! interstitial, clear  subarea

  ! Calculate the corresponding paritioning factors for interstitial aerosols
  ! using the above-derived subarea-mean mixing ratios plus the constraint that
  ! the cloud fraction weighted average of subarea mean need to match grid box mean.

  clea2gcm_ratio = max( eps, tmp_q_intrst_clea*fclea ) / max( eps, qgcm_intrst )
  clea2gcm_ratio = max( 0.0_wp, min( 1.0_wp, clea2gcm_ratio ) )

  factor_clea = clea2gcm_ratio/fclea
  factor_cldy = (1.0_wp-clea2gcm_ratio)/fcldy

end subroutine get_partition_factors

!==================================================

subroutine form_gcm_of_gases_and_aerosols_from_subareas( nsubarea, ncldy_subarea, afracsub, &! in
                                                         qsub, qqcwsub, qqcwgcm_old,        &! in
                                                         qgcm, qqcwgcm )
!---------------------------------------------------- ------------------------------------------------------
! Purpose: Form grid cell mean values by calculating area-weighted averages of the subareas.
!           - For gases and interstitial aerosols, sum over all active subareas. 
!           - For cloud-borne aerosols, 
!----------------------------------------------------------------------------------------------------------

  integer,  intent(in) :: nsubarea                   ! # of active subareas
  integer,  intent(in) :: ncldy_subarea              ! # of cloudy subareas
  real(wp), intent(in) :: afracsub(maxsubarea)       ! area fraction of subareas [unitless]

  ! The following arguments are mixing ratios. Their units do not matter for this subroutine.

  real(wp), intent(in) :: qsub   (ncnst, maxsubarea)  ! gas and interst. aerosol mixing ratios in subareas
  real(wp), intent(in) :: qqcwsub(ncnst, maxsubarea)  ! cloud-borne aerosol mixing ratios in subareas
  real(wp), intent(in) :: qqcwgcm_old(ncnst)          ! grid cell mean cloud-borne aerosol mixing ratios
                                                      ! before aerosol microphysics calculations

  real(wp), intent(out):: qgcm   (ncnst)  ! grid cell mean gas and interst.  aerosol mixing ratios
  real(wp), intent(out):: qqcwgcm(ncnst)  ! cloud-borne aerosol mixing ratios in subareas
  !---

  integer :: jsub

  ! Gases and interstitial aerosols

  qgcm(:) = 0.0_wp
  do jsub = 1, nsubarea
     qgcm(:) = qgcm(:) + qsub(:,jsub)*afracsub(jsub)
  end do

  qgcm(:) = max( 0._wp, qgcm(:) )

  ! Cloud-borne aerosols

  if (ncldy_subarea <= 0) then
     qqcwgcm(:) = qqcwgcm_old(:)
  else
     qqcwgcm(:) = 0.0_wp
     do jsub = 1, nsubarea
        qqcwgcm(:) = qqcwgcm(:) + qqcwsub(:,jsub)*afracsub(jsub)
     end do
  end if

  qqcwgcm(:) = max( 0._wp, qqcwgcm(:) )

end subroutine form_gcm_of_gases_and_aerosols_from_subareas

end module modal_aero_amicphys_subareas
