!-------------------------------------------------------------------------
! $Id: sigma_sqd_w_module.F90 6849 2014-04-22 21:52:30Z charlass@uwm.edu $
!===============================================================================
module sigma_sqd_w_module

  implicit none

  public :: compute_sigma_sqd_w, compute_a1wp3_on_wp2, compute_a2wp3_on_wp2, compute_a3p3_wp2_sqd

  private ! Default scope

  contains
!---------------------------------------------------------------------------------------------------
  elemental function compute_sigma_sqd_w( gamma_Skw_fnc, wp2, thlp2, rtp2, wpthlp, wprtp ) &
    result( sigma_sqd_w )
! Description:
!   Compute the variable sigma_sqd_w (PDF width parameter)
!
! References:
!   Eqn 22 in ``Equations for CLUBB''
!---------------------------------------------------------------------------------------------------
    use constants_clubb, only: &
      w_tol,  & ! Constant(s)
      rt_tol, &
      thl_tol

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: min, max, sqrt

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      gamma_Skw_fnc, & ! Gamma as a function of skewness   [-]
      wp2,           & ! Variance of vertical velocity     [m^2/s^2]
      thlp2,         & ! Variance of liquid pot. temp.     [K^2]
      rtp2,          & ! Variance of total water           [kg^2/kg^2]
      wpthlp,        & ! Flux of liquid pot. temp.         [m/s K]
      wprtp            ! Flux of total water               [m/s kg/kg]

    ! Output Variable
    real( kind = core_rknd ) :: sigma_sqd_w ! PDF width parameter      [-]

    ! ---- Begin Code ----

    !----------------------------------------------------------------
    ! Compute sigma_sqd_w with new formula from Vince
    !----------------------------------------------------------------

    sigma_sqd_w = gamma_Skw_fnc * &
      ( 1.0_core_rknd - min( &
                  max( ( wpthlp / ( sqrt( wp2 * thlp2 )  &
                      + 0.01_core_rknd * w_tol * thl_tol ) )**2, &
                       ( wprtp / ( sqrt( wp2 * rtp2 )  &
                      + 0.01_core_rknd * w_tol * rt_tol ) )**2 &
                     ), & ! max
             1.0_core_rknd ) & ! min - Known magic number (eq. 22 from "Equations for CLUBB")
       )

    return
  end function compute_sigma_sqd_w

  elemental function compute_a1wp3_on_wp2( gamma_Skw_fnc, wp2, wp3, thlp2, rtp2, wpthlp, wprtp, sigma_sqd_wthl, sigma_sqd_wrt ) &
    result( a1wp3_on_wp2 )
! Description:
!   Compute the variable sigma_sqd_w (PDF width parameter)
!
! References:
!   Eqn 22 in ``Equations for CLUBB''
!---------------------------------------------------------------------------------------------------
    use constants_clubb, only: &
      w_tol,  & ! Constant(s)
      rt_tol, &
      thl_tol

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: min, max, sqrt

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      gamma_Skw_fnc, & ! Gamma as a function of skewness   [-]
      wp2,           & ! Variance of vertical velocity     [m^2/s^2]
      wp3,           & ! Third order moment of vertical velocity [m^3/s^3]
      thlp2,         & ! Variance of liquid pot. temp.     [K^2]
      rtp2,          & ! Variance of total water           [kg^2/kg^2]
      wpthlp,        & ! Flux of liquid pot. temp.         [m/s K]
      wprtp,         & ! Flux of total water               [m/s kg/kg]
      sigma_sqd_wthl,& ! parameters for the contional statement [m^2/s^2 K^2 kg^2/kg^2] 
      sigma_sqd_wrt, & ! parameters for the contional statement [m^2/s^2 K^2 kg^2/kg^2]
            
    ! Output Variable
    real( kind = core_rknd ) :: a1wp3_on_wp2   ! PDF width parameter      [-]

    ! ---- Begin Code ----

    !----------------------------------------------------------------
    ! Compute a1wp3_on_wp2 with new formula from SZhang and HWan
    ! a1wp3_on_wp2 = a1*wp3/wp2
    !----------------------------------------------------------------
    !sigma_sqd_wthl > sigma_sqd_wrt: a1*wp3/wp2 = wp3*thlp2/( (1-gamma)*wp2*thlp2 + gamma*wpthlp*wpthlp )
    !sigma_sqd_wthl < sigma_sqd_wrt: a1*wp3/wp2 = wp3*rtp2/( (1-gamma)*wp2*rtp2 + gamma*wprtp*wprtp )
    !sigma_sqd_wthl = sigma_sqd_wrt: a1*wp3/wp2 = wp3*thlp2/( (1-gamma)*wp2*thlp2 + gamma*wpthlp*wpthlp )
    !                                           = wp3*rtp2/( (1-gamma)*wp2*rtp2 + gamma*wprtp*wprtp )

    a1wp3_on_wp2 = 0.0_core_rknd

    where ( ( (wp2 .ne. 0.0_core_rknd) .and. (rtp2 .ne. 0.0_core_rknd) ) .or. (wprtp .ne. 0.0_core_rknd) )
       a1wp3_on_wp2 = ( wp3 * rtp2 ) / &
                      ( (1.0_core_rknd - gamma_Skw_fnc) * wp2 * rtp2 + &
                         gamma_Skw_fnc * wprtp * wprtp &
                      )
    endwhere

    where( ( (wp2 .eq. 0.0_core_rknd) .or. (rtp2 .eq. 0.0_core_rknd) ) .and. (wprtp .eq. 0.0_core_rknd) )
       a1wp3_on_wp2 = ( wp3 * rtp2 ) / &
                      ( (1.0_core_rknd - gamma_Skw_fnc) * (wp2 + 0.01_core_rknd * w_tol * w_tol) * (rtp2 + 0.01_core_rknd * rt_tol * rt_tol)+ &
                         gamma_Skw_fnc * (wprtp + 0.01_core_rknd * w_tol * rt_tol) * (wprtp + 0.01_core_rknd * w_tol * rt_tol) &
                      )
    endwhere

    where ( (sigma_sqd_wthl > sigma_sqd_wrt) .and. ( ( (wp2 .ne. 0.0_core_rknd) .and. (thlp2 .ne. 0.0_core_rknd) ) .or. (wpthlp .ne. 0.0_core_rknd) ) )
       a1wp3_on_wp2 = ( wp3 * thlp2 ) / &
                      ( (1.0_core_rknd - gamma_Skw_fnc) * wp2 * thlp2 + &
                         gamma_Skw_fnc * wpthlp * wpthlp &
                      )
    endwhere

    where ( (sigma_sqd_wthl > sigma_sqd_wrt) .and. ( ( (wp2 .eq. 0.0_core_rknd) .or. (thlp2 .eq. 0.0_core_rknd) ) .and. (wpthlp .eq. 0.0_core_rknd) ) )
       a1wp3_on_wp2 = ( wp3 * thlp2 ) / &
                      ( (1.0_core_rknd - gamma_Skw_fnc) * (wp2 + 0.01_core_rknd * w_tol * w_tol) * (thlp2 + 0.01_core_rknd * thl_tol * thl_tol) + &
                         gamma_Skw_fnc * (wpthlp + 0.01_core_rknd * w_tol * thl_tol) * (wpthlp + 0.01_core_rknd * w_tol * thl_tol) &
                      )
    endwhere

    return
  end function compute_a1wp3_on_wp2

  elemental function compute_a2wp3_on_wp2( gamma_Skw_fnc, wp2, wp3, thlp2, rtp2, wpthlp, wprtp, sigma_sqd_wthl, sigma_sqd_wrt ) &
    result( a2wp3_on_wp2 )
! Description:
!   Compute the variable sigma_sqd_w (PDF width parameter)
!
! References:
!   Eqn 22 in ``Equations for CLUBB''
!---------------------------------------------------------------------------------------------------
    use constants_clubb, only: &
      w_tol,  & ! Constant(s)
      rt_tol, &
      thl_tol

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: min, max, sqrt

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      gamma_Skw_fnc, & ! Gamma as a function of skewness   [-]
      wp2,           & ! Variance of vertical velocity     [m^2/s^2]
      wp3,           & ! Third order moment of vertical velocity [m^3/s^3]
      thlp2,         & ! Variance of liquid pot. temp.     [K^2]
      rtp2,          & ! Variance of total water           [kg^2/kg^2]
      wpthlp,        & ! Flux of liquid pot. temp.         [m/s K]
      wprtp,         & ! Flux of total water               [m/s kg/kg]
      sigma_sqd_wthl,& ! parameters for the contional statement [m^2/s^2 K^2 kg^2/kg^2] 
      sigma_sqd_wrt, & ! parameters for the contional statement [m^2/s^2 K^2 kg^2/kg^2]

    ! Output Variable
    real( kind = core_rknd ) :: a2wp3_on_wp2   ! PDF width parameter      [-]

    ! ---- Begin Code ----

    !----------------------------------------------------------------
    ! Compute a1wp3_on_wp2 with new formula from SZhang and HWan
    ! a1wp3_on_wp2 = a1*wp3/wp2
    !----------------------------------------------------------------
    !sigma_sqd_wthl > sigma_sqd_wrt: a1*wp3/wp2 = wp3*thlp2*thlp2/( (1-gamma)*wp2*thlp2 + gamma*wpthlp*wpthlp )^2
    !sigma_sqd_wthl < sigma_sqd_wrt: a1*wp3/wp2 = wp3*rtp2*rtp2/( (1-gamma)*wp2*rtp2 + gamma*wprtp*wprtp )^2
    !sigma_sqd_wthl = sigma_sqd_wrt: a1*wp3/wp2 = wp3*thlp2*thlp2/( (1-gamma)*wp2*thlp2 + gamma*wpthlp*wpthlp )^2
    !                                           = wp3*rtp2*rtp2/( (1-gamma)*wp2*rtp2 + gamma*wprtp*wprtp )^2        


    a2wp3_on_wp2 = 0.0_core_rknd

    where ( ((wp2 .ne. 0.0_core_rknd) .and. (rtp2 .ne. 0.0_core_rknd)) .or. (wprtp .ne. 0.0_core_rknd) )
       a2wp3_on_wp2 = ( wp3 * rtp2 * rtp2 ) / &
                      ( (1.0_core_rknd - gamma_Skw_fnc) * wp2 * rtp2 + &
                         gamma_Skw_fnc * wprtp * wprtp &
                      ) / &
                      ( (1.0_core_rknd - gamma_Skw_fnc) * wp2 * rtp2 + &
                         gamma_Skw_fnc * wprtp * wprtp &
 
    endwhere

    where ( ((wp2 .eq. 0.0_core_rknd) .or. (rtp2 .eq. 0.0_core_rknd)) .and. (wprtp .eq. 0.0_core_rknd) )
       a2wp3_on_wp2 = ( wp3 * rtp2 * rtp2 ) / &
                      ( (1.0_core_rknd - gamma_Skw_fnc) * (wp2 + 0.01_core_rknd * w_tol * w_tol) * (rtp2 + 0.01_core_rknd * rt_tol * rt_tol)+ &
                         gamma_Skw_fnc * (wprtp + 0.01_core_rknd * w_tol * rt_tol) * (wprtp + 0.01_core_rknd * w_tol * rt_tol) &
                      ) / &
                      ( (1.0_core_rknd - gamma_Skw_fnc) * (wp2 + 0.01_core_rknd * w_tol * w_tol) * (rtp2 + 0.01_core_rknd * rt_tol * rt_tol)+ &
                         gamma_Skw_fnc * (wprtp + 0.01_core_rknd * w_tol * rt_tol) * (wprtp + 0.01_core_rknd * w_tol * rt_tol) &
                      )
    endwhere

    where ( (sigma_sqd_wthl > sigma_sqd_wrt) .and. ( ( (wp2 .ne. 0.0_core_rknd) .and. (thlp2 .ne. 0.0_core_rknd) ) .or. (wpthlp .ne. 0.0_core_rknd) ) )
       a2wp3_on_wp2 = ( wp3 * thlp2 * thlp2 ) / &
                      ( (1.0_core_rknd - gamma_Skw_fnc) * wp2 * thlp2 + &
                         gamma_Skw_fnc * wpthlp * wpthlp &
                      ) / &
                      ( (1.0_core_rknd - gamma_Skw_fnc) * wp2 * thlp2 + &
                         gamma_Skw_fnc * wpthlp * wpthlp &
                      )
    endwhere

    where ( (sigma_sqd_wthl > sigma_sqd_wrt) .and. ( ( (wp2 .eq. 0.0_core_rknd) .or. (thlp2 .eq. 0.0_core_rknd) ) .and. (wpthlp .eq. 0.0_core_rknd) ) )
       a2wp3_on_wp2 = ( wp3 * thlp2  * thlp2 ) / &
                      ( (1.0_core_rknd - gamma_Skw_fnc) * (wp2 + 0.01_core_rknd * w_tol * w_tol) * (thlp2 + 0.01_core_rknd * thl_tol * thl_tol) + &
                         gamma_Skw_fnc * (wpthlp + 0.01_core_rknd * w_tol * thl_tol) * (wpthlp + 0.01_core_rknd * w_tol * thl_tol) &
                      ) / &
                      ( (1.0_core_rknd - gamma_Skw_fnc) * (wp2 + 0.01_core_rknd * w_tol * w_tol) * (thlp2 + 0.01_core_rknd * thl_tol * thl_tol) + &
                         gamma_Skw_fnc * (wpthlp + 0.01_core_rknd * w_tol * thl_tol) * (wpthlp + 0.01_core_rknd * w_tol * thl_tol) &
                      ) 
    endwhere

    return
  end function compute_a2wp3_on_wp2

  elemental function compute_a3p3_wp2_sqd( gamma_Skw_fnc, wp2, wp3, thlp2, rtp2, wpthlp, wprtp, sigma_sqd_wthl, sigma_sqd_wrt ) &
    result( a3p3_wp2_sqd )
! Description:
!   Compute the variable sigma_sqd_w (PDF width parameter)
!
! References:
!   Eqn 22 in ``Equations for CLUBB''
!---------------------------------------------------------------------------------------------------
    use constants_clubb, only: &
      w_tol,  & ! Constant(s)
      rt_tol, &
      thl_tol

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: min, max, sqrt

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      gamma_Skw_fnc, & ! Gamma as a function of skewness   [-]
      wp2,           & ! Variance of vertical velocity     [m^2/s^2]
      wp3,           & ! Third order moment of vertical velocity [m^3/s^3]
      thlp2,         & ! Variance of liquid pot. temp.     [K^2]
      rtp2,          & ! Variance of total water           [kg^2/kg^2]
      wpthlp,        & ! Flux of liquid pot. temp.         [m/s K]
      wprtp,         & ! Flux of total water               [m/s kg/kg]
      sigma_sqd_wthl,& ! parameters for the contional statement [m^2/s^2 K^2 kg^2/kg^2] 
      sigma_sqd_wrt, & ! parameters for the contional statement [m^2/s^2 K^2 kg^2/kg^2]

    ! Output Variable
    real( kind = core_rknd ) :: a3p3_wp2_sqd   ! PDF width parameter      [-]

    ! ---- Begin Code ----

    !----------------------------------------------------------------
    ! Compute a3p3_wp2_sqd with new formula from SZhang and HWan
    ! a3p3_wp2_sqd = wp2*wp2*(a3+3)
    !----------------------------------------------------------------
    !sigma_sqd_wthl > sigma_sqd_wrt: (a3+3)*wp2*wp2 = (-2*gamma^2*(wp2-wpthlp*wpthlp/thlp2)^2 +4*gamma*(wp2*wp2-wp2*wpthlp*wpthlp/thlp2) - wp2*wp2 )
    !sigma_sqd_wthl < sigma_sqd_wrt: (a3+3)*wp2*wp2 = (-2*gamma^2*(wp2-wprtp*wprp/rtp2)^2 +4*gamma*(wp2*wp2-wp2*wprtp*wprtp/rtp2) - wp2*wp2 ) 
    !sigma_sqd_wthl = sigma_sqd_wrt: (a3+3)*wp2*wp2 = (-2*gamma^2*(wp2-wpthlp*wpthlp/thlp2)^2 +4*gamma*(wp2*wp2-wp2*wpthlp*wpthlp/thlp2) - wp2*wp2 )
    !                                               = (-2*gamma^2*(wp2-wprtp*wprp/rtp2)^2 +4*gamma*(wp2*wp2-wp2*wprtp*wprtp/rtp2) - wp2*wp2 ) 
    
    a3p3_wp2_sqd = 0.0_core_rknd
    
    where ( rtp2 .ne. 0.0_core_rknd )
       a3p3_wp2_sqd = ( - 2.0_core_rknd * gamma_Skw_fnc * gamma_Skw_fnc &
                          * ( wp2 - wprtp * wprtp / rtp2 ) * ( wp2 - wprtp * wprtp / rtp2 )
                        + 4.0_core_rknd * gamma_Skw_fnc * ( wp2 * wp2 - wp2 * wprtp * wprtp / rtp2 )
                        - wp2 * wp2 ) 
    endwhere
   
    where ( rtp2 .eq. 0.0_core_rknd)
       a3p3_wp2_sqd = ( - 2.0_core_rknd * gamma_Skw_fnc * gamma_Skw_fnc &
                          * ( wp2 - wprtp * wprtp / ( rtp2 + 0.01_core_rknd * rt_tol * rt_tol ) ) &
                          * ( wp2 - wprtp * wprtp / ( rtp2 + 0.01_core_rknd * rt_tol * rt_tol) ) 
                        + 4.0_core_rknd * gamma_Skw_fnc * ( wp2 * wp2 - wp2 * wprtp * wprtp / (rtp2 + 0.01_core_rknd * rt_tol * rt_tol) )
                        - wp2 * wp2 ) 
    endwhere

    where ( (sigma_sqd_wthl > sigma_sqd_wrt) .and. ( thlp2 .ne. 0.0_core_rknd) )
       a3p3_wp2_sqd = ( - 2.0_core_rknd * gamma_Skw_fnc * gamma_Skw_fnc &
                          * ( wp2 - wpthlp * wpthlp / thlp2 ) * ( wp2 - wpthlp * wpthlp / thlp2 )
                        + 4.0_core_rknd * gamma_Skw_fnc * ( wp2 * wp2 - wp2 * wpthlp * wpthlp / thlp2 )
                        - wp2 * wp2 )
    endwhere

    where ( (sigma_sqd_wthl > sigma_sqd_wrt) .and. (thlp2 .eq. 0.0_core_rknd) )
       a3p3_wp2_sqd = ( - 2.0_core_rknd * gamma_Skw_fnc * gamma_Skw_fnc &
                          * ( wp2 - wpthlp * wpthlp / (thlp2 + 0.01_core_rknd * thl_tol * thl_tol) ) &
                          * ( wp2 - wpthlp * wpthlp / (thlp2 + 0.01_core_rknd * thl_tol * thl_tol) )
                        + 4.0_core_rknd * gamma_Skw_fnc * ( wp2 * wp2 - wp2 * wpthlp * wpthlp / (thlp2 + 0.01_core_rknd * thl_tol * thl_tol) )
                        - wp2 * wp2 ) 
    endwhere

    return
  end function compute_a3p3_wp2_sqd

end module sigma_sqd_w_module
