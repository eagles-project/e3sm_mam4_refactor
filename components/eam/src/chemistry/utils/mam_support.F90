module mam_support
  !---------------------------------------------------------------------
  !Purpose:
  !This module contains utlity variables/subroutines/functions which
  ! multiple MAM routinesuse.
  !---------------------------------------------------------------------

  implicit none

  private ! make everything private

  !explicitly declare public functions/variables
  public:: min_max_bound

contains

  pure function min_max_bound(minlim, maxlim, input) result(bounded)
    !Bound a quantity between a min limit and a max limit
    real(r8), intent(in) :: minlim, maxlim
    real(r8), intent(in) :: input

    !return value
    real(r8) :: bounded

    bounded = max(min(maxlim, input), minlim)

  end function min_max_bound

end module mam_support
