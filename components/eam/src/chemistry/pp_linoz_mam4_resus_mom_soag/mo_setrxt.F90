
      module mo_setrxt

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: setrxt
      public :: setrxt_hrates

      contains

      subroutine setrxt( rate, temp, ncol )

      use ppgrid,       only : pver, pcols
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : rxntot

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol    ! number of columns
      real(r8), intent(in)    :: temp(pcols,pver)  ! midpoint temperature [K]
      real(r8), intent(inout) :: rate(ncol,pver,rxntot) ! reaction rate 

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      real(r8)  ::  itemp(ncol,pver)  ! inverse midpoint temperature [K^-1]

      itemp(:ncol,:) = 1._r8 / temp(:ncol,:)
      rate(:,:,3) = 2.9e-12_r8 * exp( -160._r8 * itemp(:,:) )
      rate(:,:,5) = 9.6e-12_r8 * exp( -234._r8 * itemp(:,:) )
      rate(:,:,7) = 1.9e-13_r8 * exp( 520._r8 * itemp(:,:) )

      end subroutine setrxt


      subroutine setrxt_hrates( rate, temp, m, ncol, kbot )

      use ppgrid,       only : pver, pcols
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      integer, intent(in) :: kbot
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer  ::  n
      real(r8)  ::  itemp(ncol,kbot)
      real(r8)  ::  exp_fac(ncol,kbot)


      end subroutine setrxt_hrates

      end module mo_setrxt
