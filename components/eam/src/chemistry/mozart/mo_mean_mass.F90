
module mo_mean_mass

  implicit none

  private
  public :: set_mean_mass, init_mean_mass

  integer :: id_o2, id_o, id_h, id_n

contains

  subroutine init_mean_mass
    use mo_chem_utls, only : get_spc_ndx

    implicit none

    id_o2 = get_spc_ndx('O2')
    id_o  = get_spc_ndx('O')
    id_h  = get_spc_ndx('H')
    id_n  = get_spc_ndx('N')

  endsubroutine init_mean_mass

  subroutine set_mean_mass( ncol, mbar )
    !-----------------------------------------------------------------
    !        ... Set the invariant densities (molecules/cm**3)
    !-----------------------------------------------------------------

    use shr_kind_mod, only : r8 => shr_kind_r8
    use ppgrid,       only : pver
    use physconst,    only : mwdry                   ! molecular weight of dry air

    implicit none

    !-----------------------------------------------------------------
    !        ... Dummy arguments
    !-----------------------------------------------------------------
    integer, intent(in)   ::      ncol                 ! number of columns
    real(r8), intent(out) ::      mbar(:,:)            ! atmos mean atomic mass [g/mol or amu]

    !-----------------------------------------------------------------
    !	... use CAM meam molecular weight 
    !-----------------------------------------------------------------
    mbar(:ncol,:pver) = mwdry  

  end subroutine set_mean_mass

end module mo_mean_mass
