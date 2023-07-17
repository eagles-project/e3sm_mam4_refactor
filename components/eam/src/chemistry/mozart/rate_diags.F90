!--------------------------------------------------------------------------------
! Manages writing reaction rates to history
!--------------------------------------------------------------------------------
module rate_diags

  use shr_kind_mod, only : r8 => shr_kind_r8
  use cam_history,  only : fieldname_len
  use cam_history,  only : addfld
  use cam_history,  only : outfld
  use chem_mods,    only : rxt_tag_cnt, rxt_tag_lst, rxt_tag_map

  implicit none
  private 
  public :: rate_diags_init
  public :: rate_diags_calc

  character(len=fieldname_len) :: rate_names(rxt_tag_cnt)

contains

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
  subroutine rate_diags_init

    integer :: ii, len, pos

    character(len=64) :: name

    do ii = 1,rxt_tag_cnt
       pos = 0
       pos = index(rxt_tag_lst(ii),'tag_')
       if (pos <= 0) pos = index(rxt_tag_lst(ii),'usr_')
       if (pos <= 0) pos = index(rxt_tag_lst(ii),'cph_')
       if (pos <= 0) pos = index(rxt_tag_lst(ii),'ion_')
       if (pos>0) then
          name = 'r_'//trim(rxt_tag_lst(ii)(5:))
       else
          name = 'r_'//trim(rxt_tag_lst(ii)(1:))
       endif
       len = min(fieldname_len,len_trim(name))
       rate_names(ii) = trim(name(1:len))
       call addfld(rate_names(ii), (/ 'lev' /),'A', 'molecules/cm3/sec','reaction rate')
    enddo

  end subroutine rate_diags_init

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
  subroutine rate_diags_calc( rxt_rates, & !inout
      vmr, air_density, ncol, lchnk ) !in

    use mo_rxt_rates_conv, only: set_rates

    real(r8), intent(inout) :: rxt_rates(:,:,:) ! reaction rates[molec/cm3/sec]
    real(r8), intent(in)    :: vmr(:,:,:)       ! volule mixing ratio [mol/mol]
    real(r8), intent(in)    :: air_density(:,:) ! air density [molecules/cm3]
    integer,  intent(in)    :: ncol, lchnk

    integer :: ii

    call set_rates( rxt_rates, vmr, ncol )
    
    do ii = 1, rxt_tag_cnt
       ! convert from vmr/sec to molecules/cm3/sec
       rxt_rates(:ncol,:,rxt_tag_map(ii)) = rxt_rates(:ncol,:,rxt_tag_map(ii)) *  air_density(:,:)
       call outfld( rate_names(ii), rxt_rates(:ncol,:,rxt_tag_map(ii)), ncol, lchnk )

    enddo
  end subroutine rate_diags_calc

end module rate_diags
