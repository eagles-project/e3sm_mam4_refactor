module kessler_autoconv

  use shr_kind_mod,   only: r8=>shr_kind_r8
  use cam_abortutils, only: endrun
  use cam_logfile,    only: iulog
  use ppgrid,         only: pver, pcols

  implicit none
  private

  public :: kessler_autoconv_tend 

contains

  !------------------------------------------------------------------------
  ! Initialization of the kessler simple microphysics model.
  ! Currently this contains just registration of output variables.
  !------------------------------------------------------------------------

  subroutine kessler_autoconv_init()

    use cam_history,    only: addfld

    implicit none

    call addfld ('KES_ql_excess',  (/'lev'/), 'I','kg/kg','Liquid water removed by kessler scheme')

  end subroutine kessler_autoconv_init


  !--------------------------------------------------------------------------------------
  ! This subroutine is an implementation of the Kessler autoconversion parameterization.
  ! For time stepping, autoconversion is treated in isolation from other processes,
  ! using an analytical solution to the autoconversion-only equation. 
  !--------------------------------------------------------------------------------------
  subroutine kessler_autoconv_tend(state, ptend, dtime, ixcldliq, ast, &
                                   kessler_autoconv_tau, kessler_autoconv_ql_crit)

  use constituents, only: pcnst
  use physics_types,only: physics_state, physics_ptend, physics_ptend_init
 !use cam_history,  only: outfld

  implicit none
  !
  ! arguments
  !
  type(physics_state), intent(in), target    :: state       ! State variables
  type(physics_ptend), intent(out)           :: ptend       ! Package tendencies

  real(r8), intent(in) :: dtime                      ! time step size (seconds) 
  integer,  intent(in) :: ixcldliq                   ! tracer index for cloud liquid 

  real(r8), intent(in) :: ast(:,:)                   ! cloud fraction
  real(r8), intent(in) :: kessler_autoconv_tau       ! autoconversion time scale
  real(r8), intent(in) :: kessler_autoconv_ql_crit   ! critical in-cloud ql for autoconversion to occur

  real(r8) :: zconst   ! ( exp(-dtime/kessler_autoconv_tau) - 1 )/dtime
  real(r8) :: ql_excess(pcols,pver)

  integer  :: ncol  
  logical  :: lq(pcnst)              ! logical array used when calling subroutine physics_ptend_init indicating 

  ncol = state%ncol
  zconst = ( exp(-dtime/kessler_autoconv_tau) - 1._r8 )/dtime

  !===================================================================
  ! Initialize ptend to be returned to the calling subroutine

  lq(:) = .FALSE.
  lq(ixcldliq) = .TRUE.
  call physics_ptend_init(ptend,state%psetcols, "Kessler autoconversion", ls=.false., lq=lq)

  !===================================================================
  ! Calculate grid-box mean liquid tendency due to autoconversion

  ptend%q(:ncol,:pver,ixcldliq) = 0._r8   ! initialize liquid tendency with zero

  ql_excess = state%q(:ncol,:pver,ixcldliq) - ast(:ncol,:pver)*kessler_autoconv_ql_crit
  where ( ql_excess(:ncol,:pver).gt.0._r8 ) 
    ptend%q(:ncol,:pver,ixcldliq) = ql_excess(:ncol,:pver)*zcost 
  end where

  call outfld('KES_ql_excess',    ql_excess,    pcols, lchnk)

  end subroutine kessler_autoconv_tend

end module kessler_autoconv
