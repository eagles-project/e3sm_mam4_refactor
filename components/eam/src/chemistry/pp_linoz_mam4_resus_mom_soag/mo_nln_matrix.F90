
module mo_nln_matrix
#include "../yaml/common_files/common_uses.ymlf90"

   use shr_kind_mod, only : r8 => shr_kind_r8

   private
   public :: nlnmat

   contains

!=========================================================
   subroutine nlnmat( mat,      & ! out
                      lmat, dti ) ! in

      use chem_mods, only : nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: dti
      real(r8), intent(in) :: lmat(nzcnt)
      real(r8), intent(out) :: mat(nzcnt)
#include "../yaml/mo_nln_matrix/f90_yaml/nlnmat_beg_yml.f90"
!----------------------------------------------
! ... local variables
!----------------------------------------------
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
         mat( 1) = lmat( 1)
         mat( 2) = lmat( 2)
         mat( 3) = lmat( 3)
         mat( 4) = lmat( 4)
         mat( 5) = lmat( 5)
         mat( 6) = lmat( 6)
         mat( 7) = lmat( 7)
         mat( 8) = lmat( 8)
         mat( 9) = lmat( 9)
         mat( 10) = lmat( 10)
         mat( 11) = lmat( 11)
         mat( 12) = lmat( 12)
         mat( 13) = lmat( 13)
         mat( 14) = lmat( 14)
         mat( 15) = lmat( 15)
         mat( 16) = lmat( 16)
         mat( 17) = lmat( 17)
         mat( 18) = lmat( 18)
         mat( 19) = lmat( 19)
         mat( 20) = lmat( 20)
         mat( 21) = lmat( 21)
         mat( 22) = lmat( 22)
         mat( 23) = lmat( 23)
         mat( 24) = lmat( 24)
         mat( 25) = lmat( 25)
         mat( 26) = lmat( 26)
         mat( 27) = lmat( 27)
         mat( 28) = lmat( 28)
         mat( 29) = lmat( 29)
         mat( 30) = lmat( 30)
         mat( 31) = lmat( 31)
         mat( 32) = lmat( 32)
         mat( 1) = mat( 1) - dti
         mat( 2) = mat( 2) - dti
         mat( 4) = mat( 4) - dti
         mat( 6) = mat( 6) - dti
         mat( 7) = mat( 7) - dti
         mat( 8) = mat( 8) - dti
         mat( 9) = mat( 9) - dti
         mat( 10) = mat( 10) - dti
         mat( 11) = mat( 11) - dti
         mat( 12) = mat( 12) - dti
         mat( 13) = mat( 13) - dti
         mat( 14) = mat( 14) - dti
         mat( 15) = mat( 15) - dti
         mat( 16) = mat( 16) - dti
         mat( 17) = mat( 17) - dti
         mat( 18) = mat( 18) - dti
         mat( 19) = mat( 19) - dti
         mat( 20) = mat( 20) - dti
         mat( 21) = mat( 21) - dti
         mat( 22) = mat( 22) - dti
         mat( 23) = mat( 23) - dti
         mat( 24) = mat( 24) - dti
         mat( 25) = mat( 25) - dti
         mat( 26) = mat( 26) - dti
         mat( 27) = mat( 27) - dti
         mat( 28) = mat( 28) - dti
         mat( 29) = mat( 29) - dti
         mat( 30) = mat( 30) - dti
         mat( 31) = mat( 31) - dti
         mat( 32) = mat( 32) - dti
#include "../yaml/mo_nln_matrix/f90_yaml/nlnmat_end_yml.f90"
   end subroutine nlnmat
!=========================================================
end module mo_nln_matrix
