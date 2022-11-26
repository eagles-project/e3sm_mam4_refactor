module modal_aero_amicphys_diags

  use modal_aero_amicphys_control, only: do_q_coltendaa, do_qqcw_coltendaa, &
                                          iqtend_cond, &
                                          iqtend_rnam, iqqcwtend_rnam, &
                                          iqtend_nnuc, iqtend_coag
!                                         suffix_q_coltendaa, suffix_qqcw_coltendaa
  implicit none

public

contains

subroutine amicphys_diags_init( do_cond, do_rename, do_newnuc, do_coag )
!-----------------------------------------------------------------------------------
! Purpose: set switches that turn on or off column-integrated tendency diagnostics.
!-----------------------------------------------------------------------------------
  logical,intent(in) :: do_cond, do_rename, do_newnuc, do_coag

  if ( (.not. do_cond) .and. (.not. do_rename) ) then
     do_q_coltendaa(:,iqtend_cond) = .false.
     do_q_coltendaa(:,iqtend_rnam) = .false.
     do_qqcw_coltendaa(:,iqqcwtend_rnam) = .false.
  end if
  if (.not.do_newnuc) do_q_coltendaa(:,iqtend_nnuc) = .false.
  if (.not.do_coag)   do_q_coltendaa(:,iqtend_coag) = .false.

end subroutine amicphys_diags_init

subroutine get_gcm_tend_diags_from_subareas( nsubarea, ncldy_subarea, afracsub, &! in
                                             qsub_tendaa, qqcwsub_tendaa,       &! in
                                             qgcm_tendaa, qqcwgcm_tendaa        )! out 
!--------------------------------------------------------------------------------
! Purpose: Get grid cell mean tendencies by calculating area-weighted averages
!          of the values in different subareas.
!--------------------------------------------------------------------------------
  use modal_aero_amicphys_control, only: wp=>r8, ncnst=>gas_pcnst, maxsubarea, nqtendaa, nqqcwtendaa

  integer,  intent(in) :: nsubarea                   ! # of active subareas
  integer,  intent(in) :: ncldy_subarea              ! # of cloudy subareas
  real(wp), intent(in) :: afracsub(maxsubarea)       ! area fraction of subareas [unitless]

  real(wp),intent(in) ::    qsub_tendaa(ncnst,   nqtendaa,maxsubarea)
  real(wp),intent(in) :: qqcwsub_tendaa(ncnst,nqqcwtendaa,maxsubarea)

  real(wp),intent(out) ::    qgcm_tendaa(ncnst,   nqtendaa)
  real(wp),intent(out) :: qqcwgcm_tendaa(ncnst,nqqcwtendaa)

  integer :: jsub

  ! Gases and interstitial aerosols

  qgcm_tendaa(:,:) = 0.0_wp
  do jsub = 1, nsubarea
     qgcm_tendaa(:,:) = qgcm_tendaa(:,:) &
                      + qsub_tendaa(:,:,jsub)*afracsub(jsub)
  end do

  ! Cloud-borne aerosols

  if (ncldy_subarea <= 0) then
     qqcwgcm_tendaa(:,:) = 0.0_wp
  else
     qqcwgcm_tendaa(:,:) = 0.0_wp
     do jsub = 1, nsubarea
        qqcwgcm_tendaa(:,:) = qqcwgcm_tendaa(:,:) &
                            + qqcwsub_tendaa(:,:,jsub)*afracsub(jsub)
     end do
  end if

end subroutine get_gcm_tend_diags_from_subareas

end module modal_aero_amicphys_diags
