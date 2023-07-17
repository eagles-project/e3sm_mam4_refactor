
module mo_mean_mass

  implicit none

  private
  public :: set_mean_mass

contains
  subroutine set_mean_mass( ncol, &! in
       mbar ) !out
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
