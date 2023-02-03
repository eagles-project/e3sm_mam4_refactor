module mo_spitfire_transport

!---------------------------------------------------------------------------------
! Purpose:
!
! Contains routines to compute tendencies from sedimentation
!
! Author: Phil Rasch
!
!---------------------------------------------------------------------------------

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use ppgrid,        only: pcols, pver, pverp
  use physconst,     only: gravit, rair
  use cam_logfile,   only: iulog
  use cam_abortutils,only: endrun

  private

  public :: getflx

contains

!===============================================================================
! Calculate tracer fluxes across cell boundaries using the 1D SPITFIRE
! (SPlit Implementation of Transport using Flux Integral REpresentation) 
! algorithm of Rasch and Lawrence (1998):
!   Rasch, P. J., and M. Lawrence, Recent development in transport methods
!   at NCAR, MPI-Rep. 265, pp. 65 – 75, Max-Planck-Inst. fuer Meteorol.,
!   Hamburg, Germany, 1998. 
!===============================================================================
  subroutine getflx(ncol, xw, phi, vel, deltat, flux)

!-----------------------------------------------------------------
! Assumed grid staggering:
!
! Input:
!
!     xw1.......xw2.......xw3.......xw4.......xw5.......xw6
!    vel1......vel2......vel3......vel4......vel5......vel6
!     ....phi1......phi2.......phi3.....phi4.......phi5....
!
! Work arrays:
!
!    psi1......psi2......psi3......psi4......psi5......psi6
!
! Output:
!
!   flux1.....flux2.....flux3.....flux4.....flux5.....flux6
!-----------------------------------------------------------------

    implicit none

    integer, intent(in) :: ncol     ! number of colums to process

    real(r8),intent(in) ::  xw(pcols,pverp)     ! coordinate variable, values at layer interfaces. In EAM this is pint.
    real(r8),intent(in) :: phi(pcols,pverp-1)   ! grid cell mean tracer mixing ratio
    real(r8),intent(in) :: vel(pcols,pverp)     ! velocity in the xw coordinate. In EAM this is grav * rho * v
                                                ! where v is velocity in the height (z) coordinate 
    real(r8),intent(in) :: deltat

    real(r8),intent(out) :: flux(pcols,pverp)

    real(r8) psi(pcols,pverp)    ! integral of phi along the xw coordinate
    real(r8) xxk(pcols,pver)     ! departure point in the xw coordinate
    real(r8) psistar(pcols)      ! integral of phi at departure point
    real(r8) fdot(pcols,pverp)   ! derivative of interpolating polynomial

    integer :: kk  ! vertical layer index

    ! Set fluxes at boundaries to zero

    flux(:ncol,1) = 0
    flux(:ncol,pverp) = 0._r8

    ! Get the vertical integral of phi. 
    ! See Rasch and Lawrence (1998), Eq (3) but note we are using a pressure coordinate here.

    psi(:ncol,1) = 0._r8
    do kk = 2,pverp
       psi(:ncol,kk) = phi(:ncol,kk-1)*( xw(:ncol,kk)-xw(:ncol,kk-1) ) + psi(:ncol,kk-1)
    end do

    ! Calculate the derivatives for the interpolating polynomial

    call cfdotmc_pro (ncol, xw, psi, &! in
                      fdot           )! out

    ! Calculate fluxes at interior interfaces

    do kk = 2,pver

       ! Find departure point. Rasch and Lawrence (1998), Eq (4)

       xxk(:ncol,kk) = xw(:ncol,kk)-vel(:ncol,kk)*deltat   

       ! Calculate the integral, psistar, at the departure point xxk.

       call cfint2(ncol, xw, psi, fdot, xxk(:,kk), &! in
                   psistar                         )! out

       ! Calculate the flux at interface kk. Rasch and Lawrence (1998), Eq (5)

       flux(:ncol,kk) = psi(:ncol,kk)-psistar(:ncol)
    end do

  end subroutine getflx



!##############################################################################

  subroutine cfint2 (ncol, x, f, fdot, xin, psistar)


    implicit none

! input
    integer ncol                      ! number of colums to process

    real (r8) x(pcols, pverp)
    real (r8) f(pcols, pverp)
    real (r8) fdot(pcols, pverp)
    real (r8) xin(pcols)

! output
    real (r8) fxdot(pcols)
    real (r8) fxdd(pcols)
    real (r8) psistar(pcols)

    integer i
    integer k
    integer intz(pcols)
    real (r8) dx
    real (r8) s
    real (r8) c2
    real (r8) c3
    real (r8) xx
    real (r8) xinf
    real (r8) psi1, psi2, psi3, psim
    real (r8) cfint
    real (r8) cfnew
    real (r8) xins(pcols)

!     the minmod function 
    real (r8) a, b, c
    real (r8) minmod
    real (r8) medan
    minmod(a,b) = 0.5_r8*(sign(1._r8,a) + sign(1._r8,b))*min(abs(a),abs(b))
    medan(a,b,c) = a + minmod(b-a,c-a)

    do i = 1,ncol
       xins(i) = medan(x(i,1), xin(i), x(i,pverp))
       intz(i) = 0
    end do

! first find the interval 
    do k =  1,pverp-1
       do i = 1,ncol
          if ((xins(i)-x(i,k))*(x(i,k+1)-xins(i)).ge.0._r8) then
             intz(i) = k
          endif
       end do
    end do

    do i = 1,ncol
       if (intz(i).eq.0) then
          write(iulog,*) ' interval was not found for col i ', i
          call endrun('DUST_SEDIMENT_MOD:cfint2 -- interval was not found ')
       endif
    end do

! now interpolate
    do i = 1,ncol
       k = intz(i)
       dx = (x(i,k+1)-x(i,k))
       s = (f(i,k+1)-f(i,k))/dx
       c2 = (3*s-2*fdot(i,k)-fdot(i,k+1))/dx
       c3 = (fdot(i,k)+fdot(i,k+1)-2*s)/dx**2
       xx = (xins(i)-x(i,k))
       fxdot(i) =  (3*c3*xx + 2*c2)*xx + fdot(i,k)
       fxdd(i) = 6*c3*xx + 2*c2
       cfint = ((c3*xx + c2)*xx + fdot(i,k))*xx + f(i,k)

!        limit the interpolant
       psi1 = f(i,k)+(f(i,k+1)-f(i,k))*xx/dx
       if (k.eq.1) then
          psi2 = f(i,1)
       else
          psi2 = f(i,k) + (f(i,k)-f(i,k-1))*xx/(x(i,k)-x(i,k-1))
       endif
       if (k+1.eq.pverp) then
          psi3 = f(i,pverp)
       else
          psi3 = f(i,k+1) - (f(i,k+2)-f(i,k+1))*(dx-xx)/(x(i,k+2)-x(i,k+1))
       endif
       psim = medan(psi1, psi2, psi3)
       cfnew = medan(cfint, psi1, psim)
       if (abs(cfnew-cfint)/(abs(cfnew)+abs(cfint)+1.e-36_r8)  .gt..03_r8) then
!     CHANGE THIS BACK LATER!!!
!     $        .gt..1) then


!     UNCOMMENT THIS LATER!!!
!            write(iulog,*) ' cfint2 limiting important ', cfint, cfnew


       endif
       psistar(i) = cfnew
    end do

    return
  end subroutine cfint2



!##############################################################################

  subroutine cfdotmc_pro (ncol, x, f, fdot)

!     prototype version; eventually replace with final SPITFIRE scheme

!     calculate the derivative for the interpolating polynomial
!     multi column version


    implicit none

! input
    integer ncol                      ! number of colums to process

    real (r8) x(pcols, pverp)
    real (r8) f(pcols, pverp)
! output
    real (r8) fdot(pcols, pverp)          ! derivative at nodes

! assumed variable distribution
!     x1.......x2.......x3.......x4.......x5.......x6     1,pverp points
!     f1.......f2.......f3.......f4.......f5.......f6     1,pverp points
!     ...sh1.......sh2......sh3......sh4......sh5....     1,pver points
!     .........d2.......d3.......d4.......d5.........     2,pver points
!     .........s2.......s3.......s4.......s5.........     2,pver points
!     .............dh2......dh3......dh4.............     2,pver-1 points
!     .............eh2......eh3......eh4.............     2,pver-1 points
!     ..................e3.......e4..................     3,pver-1 points
!     .................ppl3......ppl4................     3,pver-1 points
!     .................ppr3......ppr4................     3,pver-1 points
!     .................t3........t4..................     3,pver-1 points
!     ................fdot3.....fdot4................     3,pver-1 points


! work variables


    integer i
    integer k

    real (r8) a                    ! work var
    real (r8) b                    ! work var
    real (r8) c                    ! work var
    real (r8) s(pcols,pverp)             ! first divided differences at nodes
    real (r8) sh(pcols,pverp)            ! first divided differences between nodes
    real (r8) d(pcols,pverp)             ! second divided differences at nodes
    real (r8) dh(pcols,pverp)            ! second divided differences between nodes
    real (r8) e(pcols,pverp)             ! third divided differences at nodes
    real (r8) eh(pcols,pverp)            ! third divided differences between nodes
    real (r8) pp                   ! p prime
    real (r8) ppl(pcols,pverp)           ! p prime on left
    real (r8) ppr(pcols,pverp)           ! p prime on right
    real (r8) qpl
    real (r8) qpr
    real (r8) ttt
    real (r8) t
    real (r8) tmin
    real (r8) tmax
    real (r8) delxh(pcols,pverp)


!     the minmod function 
    real (r8) minmod
    real (r8) medan
    minmod(a,b) = 0.5_r8*(sign(1._r8,a) + sign(1._r8,b))*min(abs(a),abs(b))
    medan(a,b,c) = a + minmod(b-a,c-a)

    do k = 1,pver


!        first divided differences between nodes
       do i = 1, ncol
          delxh(i,k) = (x(i,k+1)-x(i,k))
          sh(i,k) = (f(i,k+1)-f(i,k))/delxh(i,k)
       end do

!        first and second divided differences at nodes
       if (k.ge.2) then
          do i = 1,ncol
             d(i,k) = (sh(i,k)-sh(i,k-1))/(x(i,k+1)-x(i,k-1))
             s(i,k) = minmod(sh(i,k),sh(i,k-1))
          end do
       endif
    end do

!     second and third divided diffs between nodes
    do k = 2,pver-1
       do i = 1, ncol
          eh(i,k) = (d(i,k+1)-d(i,k))/(x(i,k+2)-x(i,k-1))
          dh(i,k) = minmod(d(i,k),d(i,k+1))
       end do
    end do

!     treat the boundaries
    do i = 1,ncol
       e(i,2) = eh(i,2)
       e(i,pver) = eh(i,pver-1)
!        outside level
       fdot(i,1) = sh(i,1) - d(i,2)*delxh(i,1)  &
            - eh(i,2)*delxh(i,1)*(x(i,1)-x(i,3))
       fdot(i,1) = minmod(fdot(i,1),3*sh(i,1))
       fdot(i,pverp) = sh(i,pver) + d(i,pver)*delxh(i,pver)  &
            + eh(i,pver-1)*delxh(i,pver)*(x(i,pverp)-x(i,pver-1))
       fdot(i,pverp) = minmod(fdot(i,pverp),3*sh(i,pver))
!        one in from boundary
       fdot(i,2) = sh(i,1) + d(i,2)*delxh(i,1) - eh(i,2)*delxh(i,1)*delxh(i,2)
       fdot(i,2) = minmod(fdot(i,2),3*s(i,2))
       fdot(i,pver) = sh(i,pver) - d(i,pver)*delxh(i,pver)   &
            - eh(i,pver-1)*delxh(i,pver)*delxh(i,pver-1)
       fdot(i,pver) = minmod(fdot(i,pver),3*s(i,pver))
    end do


    do k = 3,pver-1
       do i = 1,ncol
          e(i,k) = minmod(eh(i,k),eh(i,k-1))
       end do
    end do



    do k = 3,pver-1

       do i = 1,ncol

!           p prime at k-0.5
          ppl(i,k)=sh(i,k-1) + dh(i,k-1)*delxh(i,k-1)  
!           p prime at k+0.5
          ppr(i,k)=sh(i,k)   - dh(i,k)  *delxh(i,k)

          t = minmod(ppl(i,k),ppr(i,k))

!           derivate from parabola thru f(i,k-1), f(i,k), and f(i,k+1)
          pp = sh(i,k-1) + d(i,k)*delxh(i,k-1) 

!           quartic estimate of fdot
          fdot(i,k) = pp                            &
               - delxh(i,k-1)*delxh(i,k)            &
               *(  eh(i,k-1)*(x(i,k+2)-x(i,k  ))    &
               + eh(i,k  )*(x(i,k  )-x(i,k-2))      &
               )/(x(i,k+2)-x(i,k-2))

!           now limit it
          qpl = sh(i,k-1)       &
               + delxh(i,k-1)*minmod(d(i,k-1)+e(i,k-1)*(x(i,k)-x(i,k-2)), &
               d(i,k)  -e(i,k)*delxh(i,k))
          qpr = sh(i,k)         &
               + delxh(i,k  )*minmod(d(i,k)  +e(i,k)*delxh(i,k-1),        &
               d(i,k+1)+e(i,k+1)*(x(i,k)-x(i,k+2)))

          fdot(i,k) = medan(fdot(i,k), qpl, qpr)

          ttt = minmod(qpl, qpr)
          tmin = min(0._r8,3*s(i,k),1.5_r8*t,ttt)
          tmax = max(0._r8,3*s(i,k),1.5_r8*t,ttt)

          fdot(i,k) = fdot(i,k) + minmod(tmin-fdot(i,k), tmax-fdot(i,k))

       end do

    end do

    return
  end subroutine cfdotmc_pro

end module mo_spitfire_transport
