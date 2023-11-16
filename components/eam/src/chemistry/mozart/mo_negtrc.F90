
      module mo_negtrc
#include "../yaml/common_files/common_uses.ymlf90"
      private
      public :: negtrc

      contains

      subroutine negtrc( fld, ncol )
!-----------------------------------------------------------------------
!  	... Check for negative constituent values and
!	    replace with zero value
!-----------------------------------------------------------------------

      use shr_kind_mod, only: r8 => shr_kind_r8
      use chem_mods,   only : gas_pcnst
      use ppgrid,      only : pver

      implicit none

!-----------------------------------------------------------------------
!  	... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in)          :: ncol
      real(r8), intent(inout)      :: fld(ncol,pver,gas_pcnst) ! field to check [vmr]

!-----------------------------------------------------------------------
!  	... Local variables
!-----------------------------------------------------------------------
      integer :: icnst
      integer :: nneg                       ! flag counter
#include "../yaml/mo_negtrc/f90_yaml/negtrc_beg_yml.f90"
      do icnst  = 1, gas_pcnst
         nneg = count( fld(:,:,icnst) < 0._r8 )
	 if( nneg > 0 ) then

            where( fld(:,:,icnst) < 0._r8 )
	       fld(:,:,icnst) = 0._r8
	    endwhere
	 endif
      enddo
#include "../yaml/mo_negtrc/f90_yaml/negtrc_end_yml.f90"
      end subroutine negtrc

      end module mo_negtrc
