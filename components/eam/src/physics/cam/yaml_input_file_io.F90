module yaml_input_file_io

  !Following three module are not used here but they are added for convenience so that
  !we just need to include yaml_input_file_io to get access to icolprnt, getunit
  !and freeunit. It may be a bad practice but we are doing it here for temporary
  !refactoring work

  use module_perturb, only: icolprnt
  use units,          only: getunit, freeunit
  use time_manager,   only: get_nstep

  use shr_kind_mod, only: r8 => shr_kind_r8 !real kind
  use shr_log_mod , only: errMsg => shr_log_errMsg
  use cam_abortutils, only: endrun

  implicit none

  private ! make everything private in this module

  !explicit public subroutines and variables
  public :: open_files
  public :: write_var_with_levs
  public :: write_output_var_with_levs
  public :: write_input_header
  public :: write_output_header
  public :: icolprnt, getunit, freeunit, get_nstep

contains

  subroutine open_files(name, & ! input
       unit_input, unit_output) ! output

    !intent-ins
    character(len=*), intent(in) :: name

    !intent outs
    integer, intent(out) ::unit_input, unit_output

    !local
    character(len=2000) :: finp, fout

    write(finp,'(2A,I1,A)')trim(adjustl(name)),'_input_ts_',get_nstep(),'.yaml'
    write(fout,'(2A,I1,A)')trim(adjustl(name)),'_output_ts_',get_nstep(),'.py'

    !get a unused unit numbers to write
    unit_input = getunit()
    unit_output = getunit()

    !open I/O files to write in the append mode
    open( unit_input, file=trim(adjustl(finp)), action='write', position='append')
    open( unit_output, file=trim(adjustl(fout)), action='write', position='append')

  end subroutine open_files


  subroutine write_var_with_levs(unit_input, unit_output, fld_name,dim,field)
    !------------------------------------------------------------------
    !Purpose: Writes a 1D input and output field in a YAML file format
    !for a given column
    !------------------------------------------------------------------
    implicit none

    integer, intent(in)   :: unit_input      ! input stream unit number
    integer, intent(in)   :: unit_output     ! output stream unit number
    character(len=*), intent(in) :: fld_name ! name of the field
    integer, intent(in)   :: dim             ! dimensions of the field
    real(r8), intent(in)  :: field(:)        ! field values in r8

    integer :: k

    !check if file is open to write or not
    call is_file_open(unit_input)

    !format statement to write in double precision
10  format(E17.10)
11  format(A,E17.10)

    write(unit_input,'(3A)',advance="no")'    ',trim(adjustl(fld_name)),':['

    write(unit_input,10,advance="no")field(1)

    do k = 2, dim
       write(unit_input,11,advance="no")',',field(k)
    enddo

    write(unit_input,'(A)')']'

    call write_output_var_with_levs(unit_output, fld_name, dim, field, "input")

  end subroutine write_var_with_levs

  !================================================================================
  !================================================================================

  subroutine write_output_var_with_levs(unit_output, fld_name, dim, field, inp_out_str)
    !------------------------------------------------------------------
    !Purpose: Writes a 1D output field in a YAML file format
    !for a given column
    !------------------------------------------------------------------
    implicit none

    integer, intent(in)   :: unit_output     ! output stream unit number
    character(len=*), intent(in) :: fld_name ! name of the field
    integer, intent(in)   :: dim             ! dimensions of the field
    real(r8), intent(in)  :: field(:)        ! field values in r8

    !optional input
    character(len=*), intent(in), optional :: inp_out_str ! input or output

    !local
    integer :: k
    character(len=20) :: object

    !check if file is open to write or not
    call is_file_open(unit_output)

    !format statement to write in double precision
12  format(E17.10,A)

    object = "output"
    if (present(inp_out_str)) then
       object = trim(adjustl(inp_out_str))
    endif

    write(unit_output,'(4A)',advance="no")trim(adjustl(object)),'.',trim(adjustl(fld_name)),'=['

    do k = 1, dim
       write(unit_output,12,advance="no"),field(k),','
    enddo

    write(unit_output,'(A)')']'

  end subroutine write_output_var_with_levs

  !================================================================================
  !================================================================================

  subroutine write_input_header(unit_input, unit_output)
    !------------------------------------------------------------------
    !Purpose: Write input header for the YAML input and the output files
    !------------------------------------------------------------------
    integer, intent(in) :: unit_input, unit_output

    !check if file is open to write or not
    call is_file_open(unit_input)
    call is_file_open(unit_output)

    write(unit_input,'(A)')'input:'
    write(unit_input,'(A)')'  enumerated:'


    write(unit_output,'(A)')'# Object is just a dynamic container that stores input/output data.'
    write(unit_output,'(A)')'class Object(object):'
    write(unit_output,'(A)')'    pass'

    write(unit_output,'(A)')'# Settings are stored here.'
    write(unit_output,'(A)')'settings = Object()'
    write(unit_output,'(A)')'# Input is stored here.'
    write(unit_output,'(A)')'input = Object()'

  end subroutine write_input_header

  !================================================================================
  !================================================================================

  subroutine write_output_header(unit_output)
    !------------------------------------------------------------------
    !Purpose: Write output header for the YAML input
    !------------------------------------------------------------------
    integer, intent(in) :: unit_output

    !check if file is open to write or not
    call is_file_open(unit_output)

    !write output header
    write(unit_output,'(A)')'# Output data is stored here.'
    write(unit_output,'(A)')'output = Object()'

  end subroutine write_output_header

  !================================================================================
  !================================================================================

  subroutine is_file_open(unit)
    !------------------------------------------------------------------
    !Purpose: Check if the file is open or not
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
