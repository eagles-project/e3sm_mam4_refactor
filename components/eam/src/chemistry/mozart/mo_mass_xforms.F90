module mo_mass_xforms
#include "../yaml/common_files/common_uses.ymlf90"
  use ppgrid,       only : pcols, pver
  use shr_kind_mod, only : r8 => shr_kind_r8


  private
  public :: mmr2vmr, vmr2mmr, h2o_to_vmr, init_mass_xforms
  save

  real(r8) :: adv_mass_h2o = 18._r8

contains

  subroutine init_mass_xforms

    implicit none

    adv_mass_h2o = 18._r8

  endsubroutine init_mass_xforms

  subroutine mmr2vmr( mmr, vmr, mbar, ncol )
    !-----------------------------------------------------------------
    !	... Xfrom from mass to volume mixing ratio
    !-----------------------------------------------------------------
    use chem_mods, only : adv_mass, gas_pcnst


    implicit none

    !-----------------------------------------------------------------
    !	... Dummy args
    !-----------------------------------------------------------------
    integer, intent(in)     :: ncol ! number of columns
    real(r8), intent(in)    :: mbar(ncol,pver)  ! atmos mean atomic mass [g/mol or amu]
    real(r8), intent(in)    :: mmr(pcols,pver,gas_pcnst)  ! mass mixing ratio [kg/kg air]
    real(r8), intent(inout) :: vmr(ncol,pver,gas_pcnst) ! volume mixing ratio [mol/mol air]

    !-----------------------------------------------------------------
    !	... Local variables
    !-----------------------------------------------------------------
    integer :: kk, mm   ! indices for vertical level and  gas constituent
#include "../yaml/mo_mass_xforms/f90_yaml/mmr2vmr_beg_yml.f90"
    do mm = 1,gas_pcnst
       if( adv_mass(mm) /= 0._r8 ) then
          do kk = 1,pver
             vmr(:ncol,kk,mm) = mbar(:ncol,kk) * mmr(:ncol,kk,mm) / adv_mass(mm)
          enddo
       endif
    enddo
#include "../yaml/mo_mass_xforms/f90_yaml/mmr2vmr_end_yml.f90"
  end subroutine mmr2vmr

  subroutine vmr2mmr( vmr, mmr, mbar, ncol )
    !-----------------------------------------------------------------
    !	... Xfrom from volume to mass mixing ratio
    !-----------------------------------------------------------------

    use chem_mods, only : adv_mass, gas_pcnst

    implicit none

    !-----------------------------------------------------------------
    !	... Dummy args
    !-----------------------------------------------------------------
    integer, intent(in)     :: ncol ! number of columns
    real(r8), intent(in)    :: mbar(ncol,pver)  ! atmos mean atomic mass [g/mol or amu]
    real(r8), intent(in)    :: vmr(ncol,pver,gas_pcnst) ! volume mixing ratio [mol/mol air]
    real(r8), intent(inout) :: mmr(pcols,pver,gas_pcnst) ! mass mixing ratio [kg/kg air]

    !-----------------------------------------------------------------
    !	... Local variables
    !-----------------------------------------------------------------
    integer :: kk, mm ! indices for vertical level and  gas constituent

    !-----------------------------------------------------------------
    !	... The non-group species
    !-----------------------------------------------------------------
    do mm = 1,gas_pcnst
       if( adv_mass(mm) /= 0._r8 ) then
          do kk = 1,pver
             mmr(:ncol,kk,mm) = adv_mass(mm) * vmr(:ncol,kk,mm) / mbar(:ncol,kk)
          end do
       end if
    end do

  end subroutine vmr2mmr

  subroutine h2o_to_vmr( h2o_mmr, h2o_vmr, mbar, ncol )
    !-----------------------------------------------------------------------
    !     ... Transform water vapor from mass to volumetric mixing ratio
    !-----------------------------------------------------------------------

    implicit none

    !-----------------------------------------------------------------------
    !	... Dummy arguments
    !-----------------------------------------------------------------------
    integer, intent(in) ::    ncol  !  number of columns
    real(r8), dimension(pcols,pver), intent(in) :: &
         h2o_mmr                ! specific humidity [kg h2o / kg air]
    real(r8), dimension(ncol,pver), intent(in)  :: &
         mbar                   ! atmos mean atomic mass [g/mol or amu]
    real(r8), dimension(ncol,pver), intent(out) :: &
         h2o_vmr                ! water vapor vmr [mol h2o/ mol air]

    !-----------------------------------------------------------------------
    !	... Local variables
    !-----------------------------------------------------------------------
    integer ::   kk  ! vertical level index
#include "../yaml/mo_mass_xforms/f90_yaml/h2o_to_vmr_beg_yml.f90"

    do kk = 1,pver
       h2o_vmr(:ncol,kk) = mbar(:ncol,kk) * h2o_mmr(:ncol,kk) / adv_mass_h2o
    enddo
#include "../yaml/mo_mass_xforms/f90_yaml/h2o_to_vmr_end_yml.f90"
  end subroutine h2o_to_vmr

end module mo_mass_xforms
