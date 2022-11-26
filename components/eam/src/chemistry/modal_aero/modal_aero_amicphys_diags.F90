module modal_aero_amicphys_diags

  use modal_aero_amicphys_control, only:  nqtendaa, nqqcwtendaa, &
                                          do_q_coltendaa, do_qqcw_coltendaa, iqtend_cond, &
                                          iqtend_rnam, iqqcwtend_rnam, &
                                          iqtend_nnuc, iqtend_coag, &
                                          suffix_q_coltendaa, suffix_qqcw_coltendaa
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

end module modal_aero_amicphys_diags
