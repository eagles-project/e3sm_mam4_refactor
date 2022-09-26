module yaml_input_file_io

  !Following two module are not used here but they are added for convenience so that
  !we just need to include yaml_input_file_io to get access to icolprnt, getunit
  !and freeunit it may be a bad practice but we are doing it here for temporary
  !refactoring work

  use module_perturb, only: icolprnt
  use units,          only: getunit, freeunit

  use shr_kind_mod, only: r8 => shr_kind_r8 !real kind
  use shr_log_mod , only: errMsg => shr_log_errMsg
  use cam_abortutils, only: endrun

  implicit none

  private ! make everything private in this module

  !explicit public subroutines and variables
  public :: write_var_with_levs, icolprnt, getunit, freeunit

contains
  subroutine write_var_with_levs(unit,fld_name,dim,field)
    !------------------------------------------------------------------
    !Purpose: Writes a 1D field in a YAML file format
    !for a given column
    !------------------------------------------------------------------
    implicit none

    integer, intent(in)   :: unit            ! stream unit number
    character(len=*), intent(in) :: fld_name ! name of the field
    integer, intent(in)   :: dim             ! dimensions of the field
    real(r8), intent(in)  :: field(:)        ! field values in r8

    integer :: k

    !check if file is open to write or not
    call is_file_open(unit)

    !format statement to write in double precision
10  format(d17.10)

    write(unit,'(3A)',advance="no")'    ',trim(adjustl(fld_name)),':['
    do k = 1, dim
       write(unit,10,advance="no")field(k)
    enddo

    write(unit,'(A)')']'

    end subroutine write_var_with_levs

    !================================================================================
    !================================================================================

    subroutine is_file_open(unit)
    !------------------------------------------------------------------
    !Purpose: Writes a 1D field in a YAML file format
    !for a given column
    !------------------------------------------------------------------
      integer, intent(in) :: unit

      !local
      logical ::  is_open
      character(len=2000) :: err_msg

      inquire(unit=unit, opened=is_open) !inquire if the unit is open
      if ( .not. is_open ) then ! error out if not open
         write(err_msg,*)"Unit ", unit," is not open to write"//errmsg(__FILE__,__LINE__)
         call endrun(err_msg)
      end if

    end subroutine is_file_open

end module yaml_input_file_io
