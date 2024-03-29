module mo_indprd
#include "../yaml/common_files/common_uses.ymlf90"

   use shr_kind_mod, only : r8 => shr_kind_r8

   private
   public :: indprd

   contains

!==========================================================================
   subroutine indprd( class, & ! in
                       prod, & ! inout
                nprod, extfrc, rxt, ncol ) ! in

      use chem_mods, only : gas_pcnst, extcnt, rxntot
      use ppgrid, only : pver

      implicit none

!--------------------------------------------------------------------
! ... dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: class  ! 1 for explicit, 4 for implicit
      integer, intent(in) :: ncol
      integer, intent(in) :: nprod
      real(r8), intent(in) :: rxt(ncol,pver,rxntot)
      real(r8), intent(in) :: extfrc(ncol,pver,extcnt)
      real(r8), intent(inout) :: prod(ncol,pver,nprod)
#include "../yaml/mo_indprd/f90_yaml/indprd_beg_yml.f90"

!--------------------------------------------------------------------
! ... "independent" production for Explicit species
!--------------------------------------------------------------------
      if( class == 1 ) then     ! FORTRAN refactor: this if condition is never reached as the explicit solver is removed
         prod(:,:,1) = 0._r8
!--------------------------------------------------------------------
! ... "independent" production for Implicit species
!--------------------------------------------------------------------
      elseif( class == 4 ) then
         prod(:,:,1) =rxt(:,:,2)
         prod(:,:,2) = 0._r8
         prod(:,:,3) = + extfrc(:,:,1)
         prod(:,:,4) = 0._r8
         prod(:,:,5) = + extfrc(:,:,9)
         prod(:,:,6) = + extfrc(:,:,2)
         prod(:,:,7) = 0._r8
         prod(:,:,8) = 0._r8
         prod(:,:,9) = 0._r8
         prod(:,:,10) = 0._r8
         prod(:,:,11) = 0._r8
         prod(:,:,12) = 0._r8
         prod(:,:,13) = + extfrc(:,:,6)
         prod(:,:,14) = + extfrc(:,:,3)
         prod(:,:,15) = 0._r8
         prod(:,:,16) = 0._r8
         prod(:,:,17) = 0._r8
         prod(:,:,18) = + extfrc(:,:,7)
         prod(:,:,19) = 0._r8
         prod(:,:,20) = 0._r8
         prod(:,:,21) = 0._r8
         prod(:,:,22) = 0._r8
         prod(:,:,23) = 0._r8
         prod(:,:,24) = 0._r8
         prod(:,:,25) = 0._r8
         prod(:,:,26) = 0._r8
         prod(:,:,27) = + extfrc(:,:,4)
         prod(:,:,28) = + extfrc(:,:,5)
         prod(:,:,29) = 0._r8
         prod(:,:,30) = + extfrc(:,:,8)

      endif

#include "../yaml/mo_indprd/f90_yaml/indprd_end_yml.f90"
   end subroutine indprd
!==========================================================================

end module mo_indprd
