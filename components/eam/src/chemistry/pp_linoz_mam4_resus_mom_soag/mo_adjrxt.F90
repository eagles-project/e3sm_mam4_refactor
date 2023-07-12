






      module mo_adjrxt
#include "../yaml/common_files/common_uses.ymlf90"
      private
      public :: adjrxt

      contains

      subroutine adjrxt( rate, &  ! inout
                          inv, mtot, ncol )  ! in

      use ppgrid, only : pver
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : nfs, rxntot

      implicit none

!--------------------------------------------------------------------
! ... dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: ncol  ! number of columns
      real(r8), intent(in) :: inv(ncol,pver,nfs) ! invariant molec. density [molecules/cm^3]
      real(r8), intent(in) :: mtot(ncol,pver) ! total atm molec. density [molecules/cm^3]
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)  ! reaction rate

!--------------------------------------------------------------------
! ... local variables
!--------------------------------------------------------------------
      real(r8) :: imtot(ncol,pver)  ! inverse of total atm molec. density [cm^3/molecules]
#include "../yaml/mo_adjrxt/f90_yaml/adjrxt_beg_yml.f90"

         rate(:,:, 3) = rate(:,:, 3) * inv(:,:, 5)
         rate(:,:, 4) = rate(:,:, 4) * inv(:,:, 5)
         rate(:,:, 5) = rate(:,:, 5) * inv(:,:, 5)
         rate(:,:, 6) = rate(:,:, 6) * inv(:,:, 5)
         rate(:,:, 7) = rate(:,:, 7) * inv(:,:, 6)
         imtot(:,:) = 1._r8 / mtot(:,:)
         rate(:,:, 2) = rate(:,:, 2) * inv(:,:, 7) * inv(:,:, 7) * imtot(:,:)
#include "../yaml/mo_adjrxt/f90_yaml/adjrxt_end_yml.f90"
      end subroutine adjrxt

      end module mo_adjrxt
