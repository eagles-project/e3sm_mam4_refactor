
      module mo_setrxt

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: setrxt

      contains

      subroutine setrxt( rate, & ! inout
                         temp, ncol )   ! in

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

      end module mo_setrxt
