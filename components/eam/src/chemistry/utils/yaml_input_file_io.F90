module yaml_input_file_io

  !Following three module are not used here but they are added for convenience so that
  !we just need to include yaml_input_file_io to get access to icolprnt, getunit
  !and freeunit. It may be a bad practice but we are doing it here for temporary
  !refactoring work

  use module_perturb, only: icolprnt
  use units,          only: getunit, freeunit
  use time_manager,   only: get_nstep

  use shr_kind_mod,   only: r8 => shr_kind_r8 !real kind
  use shr_log_mod ,   only: errMsg => shr_log_errMsg
  use cam_abortutils, only: endrun
  use constituents,   only: pcnst, cnst_name
  use phys_grid,      only: get_rlat_p, get_rlon_p


  implicit none

  private ! make everything private in this module

  !explicit public subroutines and variables
  public :: open_files   !opens input and output files to write
  public :: write_input_output_header !writes input related  header for both input and output files
  public :: write_output_header       !writes output related  header for output file
  public :: one_print_ts !enables only one sample write at each time step
  public :: yaml_vars ! data structure to carry yaml vars

  public :: write_var    !write a single variable
  public :: write_output_var  !write a single output variable
  public :: write_1d_var ! writes 1D variables of any dimension
  public :: write_1d_output_var !writes 1D variables of any dimension in the output python module
  public :: write_2d_var ! writes 2D variables of any dimension
  public :: write_2d_output_var !writes 2D variables of any dimension in the output python module
  public :: write_aerosol_mmr_from_stateq   !writes aerosol mmr from state%q or q vector (cloud borne and interstitial)
  public :: write_output_aerosol_mmr_from_stateq !writes aerosol mmr from state%q or q vector(cloud borne and interstitial) in the output python module

  !Share data from other modules
  public :: icolprnt, getunit, freeunit, get_nstep

  !A structure to carry all the yaml input/output variables
  type yaml_vars
     integer :: col_print, lev_print, nstep_print
     integer :: lchnk_print = huge(1)
     logical :: flag_print = .false.
  end type yaml_vars


  !interfaces for datatypes (integers and reals)
  interface write_var
     module procedure write_var_int
     module procedure write_var_real
     module procedure write_var_logical
  end interface write_var

  interface write_output_var
     module procedure write_output_var_int
     module procedure write_output_var_real
     module procedure write_output_var_logical
  end interface write_output_var

  interface write_aerosol_mmr_from_stateq
     module procedure write_aerosol_mmr_from_stateq_real
     module procedure write_aerosol_mmr_from_stateq_logical
  end interface write_aerosol_mmr_from_stateq

  interface write_output_aerosol_mmr_from_stateq
     module procedure write_output_aerosol_mmr_from_stateq_real
     module procedure write_output_aerosol_mmr_from_stateq_logical
  end interface write_output_aerosol_mmr_from_stateq

contains

  !================================================================================
  !------------------- * Utilities functions/subroutines * -----------------------
  !================================================================================

  subroutine open_files(name, & ! input
       unit_input, unit_output) ! output

    !intent-ins
    character(len=*), intent(in) :: name

    !intent outs
    integer, intent(out) ::unit_input, unit_output

    !local
    character(len=2000) :: finp, fout

    write(finp,'(2A,I0,A)')trim(adjustl(name)),'_input_ts_',get_nstep(),'.yaml'
    write(fout,'(2A,I0,A)')trim(adjustl(name)),'_output_ts_',get_nstep(),'.py'

    !get a unused unit numbers to write
    unit_input = getunit()
    unit_output = getunit()

    !open I/O files to write in the append mode
    open( unit_input, file=trim(adjustl(finp)), action='write', position='append')
    open( unit_output, file=trim(adjustl(fout)), action='write', position='append')

  end subroutine open_files

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

  !================================================================================

  function one_print_ts(print_time_step, do_update_print_time)
    !output
    integer, intent(inout) :: print_time_step

    !input
    !flag to decide if print time should be updated or not
    !(print_time_step is updated by 1 if do_update_print_time is .true.
    !so that the next printout happens on the next time step)
    logical, intent(in), optional :: do_update_print_time

    !local (return value)
    logical :: one_print_ts

    if(print_time_step == get_nstep()) then
       if(present(do_update_print_time) .and. do_update_print_time) print_time_step = print_time_step + 1
       one_print_ts = .true.
    else
       one_print_ts = .false.
    endif
    return

  end function one_print_ts

  !================================================================================
  !-------------------------- * Write headers * -----------------------------------
  !================================================================================

  subroutine write_input_output_header(unit_input, unit_output,lchunk, icol, subr_name, nstep, lev, dt)
    !------------------------------------------------------------------
    !Purpose: Write input header for the YAML input and the output files
    !------------------------------------------------------------------
    integer,  intent(in) :: unit_input, unit_output, lchunk, icol, nstep, lev
    real(r8), intent(in), optional :: dt !time step

    character(len = *), intent(in) :: subr_name

    !local
    character(len = 2000) :: meta_msg
    integer, parameter :: RAD_TO_DEG = 57.296_r8 !180/pi
    real(r8) :: time_step

    time_step = 0.0_r8 !defaut time step
    if (present(dt)) time_step = dt


    !check if file is open to write or not
    call is_file_open(unit_input)
    call is_file_open(unit_output)

    !write some meta data
    write(meta_msg,*)'# Data at [lat=', get_rlat_p(lchunk, icol)*RAD_TO_DEG, &
         ', lon=',get_rlon_p(lchunk, icol)*RAD_TO_DEG,', k=',lev,', time step=', nstep,']'

    !write header for the input file
    write(unit_input,'(A)')adjustl(trim(meta_msg))
    write(unit_input,'(A)')'mam4xx:'
    write(unit_input,'(A,A)')'  function: ',trim(adjustl(subr_name))
    write(unit_input,'(A)')'input:'
    write(unit_input,'(A)' )'  fixed:'
    write(unit_input,'(A,F8.2)')'    dt:', time_step

    !write header for the output file
    write(unit_output,'(A)')adjustl(trim(meta_msg))
    write(unit_output,'(A)')'# This file was generated by E3SM.'
    write(unit_output,'(A)')
    write(unit_output,'(A)')'from math import nan as nan, inf as inf'
    write(unit_output,'(A)')

    write(unit_output,'(A)')'# Object is just a dynamic container that stores input/output data.'
    write(unit_output,'(A)')'class Object(object):'
    write(unit_output,'(A)')'    pass'

    write(unit_output,'(A)')'# Settings are stored here.'
    write(unit_output,'(A)')'settings = Object()'
    write(unit_output,'(A)')'# Input is stored here.'
    write(unit_output,'(A)')'input = Object()'
    write(unit_output,'(A,F8.2,A)')'input.dt = [',time_step,', ]'

  end subroutine write_input_output_header

  !================================================================================

  subroutine write_output_header(unit_output)
    !------------------------------------------------------------------
    !Purpose: Write output header in the output python module
    !------------------------------------------------------------------
    integer, intent(in) :: unit_output

    !check if file is open to write or not
    call is_file_open(unit_output)

    !write output header
    write(unit_output,'(A)')'# Output data is stored here.'
    write(unit_output,'(A)')'output = Object()'

  end subroutine write_output_header

  !================================================================================
  !-------------------------- * write data * --------------------------------------
  !================================================================================

  subroutine write_aerosol_mmr_from_stateq_real(unit_input, unit_output, fld_name,field,aer_num_only)
    !------------------------------------------------------------------
    !Purpose: write interstitial aerosols mmr from state%q vector
    !------------------------------------------------------------------
    implicit none

    integer,  intent(in)   :: unit_input      ! input stream unit number
    integer,  intent(in)   :: unit_output     ! output stream unit number
    character(len=*), intent(in) :: fld_name ! name of the field
    real(r8), intent(in)  :: field(:)        ! field values in r8

    !optional input
    logical,  intent(in), optional :: aer_num_only!to print only aerosol numbers

    integer :: ispec, aer_spec
    logical :: mass_or_num


    !format statement to write in double precision
10  format(E17.10)
11  format(A,E17.10)

    write(unit_input,'(3A)',advance="no")'    ',trim(adjustl(fld_name)),': ['
    !In MAM4,first 15 species are non-aerosols, so we start with 16th species
    !We store first specie mmr seprately to adjust comma(",") we need
    !in the output array
    aer_spec = 16
    if(present(aer_num_only) .and. aer_num_only)aer_spec = 23 !23rd index is the first num aerosol
    write(unit_input,10,advance="no")field(aer_spec)
    aer_spec = aer_spec + 1

    do ispec = aer_spec, pcnst
       !ignore species based on "aer_num_only" optional input
       mass_or_num = (index(trim(adjustl(cnst_name(ispec))),'num')==0)
       if(present(aer_num_only) .and. aer_num_only) mass_or_num = (index(trim(adjustl(cnst_name(ispec))),'num').ne.0) !if we need to print only aerosol numbers

       if(mass_or_num) then !ignore the aerosol number
          write(unit_input,11,advance="no")',',field(ispec)
       endif

    enddo
    write(unit_input,'(A)')']'

    call write_output_aerosol_mmr_from_stateq_real(unit_output, fld_name, field, aer_num_only, "input")

  end subroutine write_aerosol_mmr_from_stateq_real

  !================================================================================

  subroutine write_output_aerosol_mmr_from_stateq_real(unit_output, fld_name, field, aer_num_only, inp_out_str)
    !------------------------------------------------------------------
    !Purpose: Writes aerosol mmr from state%q or q vector for the input and the output yaml file
    !------------------------------------------------------------------
    implicit none

    integer, intent(in)   :: unit_output     ! output stream unit number
    character(len=*), intent(in) :: fld_name ! name of the field
    real(r8), intent(in)  :: field(:)        ! field values in r8

    !optional input
    logical,  intent(in), optional :: aer_num_only!to print only aerosol numbers
    character(len=*), intent(in), optional :: inp_out_str ! input or output

    !local
    integer :: ispec, aer_spec
    logical :: mass_or_num
    character(len=20) :: object

    !check if file is open to write or not
    call is_file_open(unit_output)

    !format statement to write in double precision
12  format(E17.10,A)

    object = "output"
    if (present(inp_out_str)) then
       object = trim(adjustl(inp_out_str))
    endif

    write(unit_output,'(4A)',advance="no")trim(adjustl(object)),'.',trim(adjustl(fld_name)),'=[['

    !In MAM4,first 15 species are non-aerosols, so we start with 16th species
    aer_spec = 16

    !if we need to print only aerosol number mixing ratios, we should start from 23rd species
    if(present(aer_num_only) .and. aer_num_only)aer_spec = 23 !23rd index is the first num aerosol

    do ispec = aer_spec, pcnst
       !ignore species based on "aer_num_only" optional input
       mass_or_num = (index(trim(adjustl(cnst_name(ispec))),'num')==0)
       if(present(aer_num_only) .and. aer_num_only) &
            mass_or_num = (index(trim(adjustl(cnst_name(ispec))),'num').ne.0) !if we need to print only aerosol numbers

       if(mass_or_num) then !ignore the aerosol number
          write(unit_output,12,advance="no")field(ispec),','
       endif

    enddo

    write(unit_output,'(A)')'],]'

  end subroutine write_output_aerosol_mmr_from_stateq_real

  !================================================================================

  subroutine write_aerosol_mmr_from_stateq_logical(unit_input, unit_output, fld_name,field,aer_num_only)
    !------------------------------------------------------------------
    !Purpose: write interstitial aerosols mmr from state%q or q vector
    !------------------------------------------------------------------
    implicit none

    integer,  intent(in)   :: unit_input      ! input stream unit number
    integer,  intent(in)   :: unit_output     ! output stream unit number
    character(len=*), intent(in) :: fld_name ! name of the field
    logical, intent(in)  :: field(:)        ! field values in boolean

    !optional input
    logical,  intent(in), optional :: aer_num_only!to print only aerosol numbers

    integer :: ispec, aer_spec
    logical :: mass_or_num


    !format statement to write in double precision
10  format(E17.10)
11  format(A,L10)

    write(unit_input,'(3A)',advance="no")'    ',trim(adjustl(fld_name)),': ['
    !In MAM4,first 15 species are non-aerosols, so we start with 16th species
    !We store first specie mmr seprately to adjust comma(",") we need
    !in the output array
    aer_spec = 16
    if(present(aer_num_only) .and. aer_num_only)aer_spec = 23 !23rd index is the first num aerosol
    write(unit_input,10,advance="no")field(aer_spec)
    aer_spec = aer_spec + 1

    do ispec = aer_spec, pcnst
       !ignore species based on "aer_num_only" optional input
       mass_or_num = (index(trim(adjustl(cnst_name(ispec))),'num')==0)
       if(present(aer_num_only) .and. aer_num_only) mass_or_num = (index(trim(adjustl(cnst_name(ispec))),'num').ne.0) !if we need to print only aerosol numbers

       if(mass_or_num) then !ignore the aerosol number
          write(unit_input,11,advance="no")',',field(ispec)
       endif

    enddo
    write(unit_input,'(A)')']'

    call write_output_aerosol_mmr_from_stateq_logical(unit_output, fld_name, field, aer_num_only, "input")

  end subroutine write_aerosol_mmr_from_stateq_logical

  !================================================================================

  subroutine write_output_aerosol_mmr_from_stateq_logical(unit_output, fld_name, field, aer_num_only, inp_out_str)
    !------------------------------------------------------------------
    !Purpose: Writes aerosol mmr for the input and the output yaml file
    !------------------------------------------------------------------
    implicit none

    integer, intent(in)   :: unit_output     ! output stream unit number
    character(len=*), intent(in) :: fld_name ! name of the field
    logical, intent(in)  :: field(:)        ! field values in boolean

    !optional input
    logical,  intent(in), optional :: aer_num_only!to print only aerosol numbers
    character(len=*), intent(in), optional :: inp_out_str ! input or output

    !local
    integer :: ispec, aer_spec
    logical :: mass_or_num
    character(len=20) :: object

    !check if file is open to write or not
    call is_file_open(unit_output)

    !format statement to write in double precision
12  format(L10,A)

    object = "output"
    if (present(inp_out_str)) then
       object = trim(adjustl(inp_out_str))
    endif

    write(unit_output,'(4A)',advance="no")trim(adjustl(object)),'.',trim(adjustl(fld_name)),'=[['

    !In MAM4,first 15 species are non-aerosols, so we start with 16th species
    aer_spec = 16

    !if we need to print only aerosol number mixing ratios, we should start from 23rd species
    if(present(aer_num_only) .and. aer_num_only)aer_spec = 23 !23rd index is the first num aerosol

    do ispec = aer_spec, pcnst
       !ignore species based on "aer_num_only" optional input
       mass_or_num = (index(trim(adjustl(cnst_name(ispec))),'num')==0)
       if(present(aer_num_only) .and. aer_num_only) &
            mass_or_num = (index(trim(adjustl(cnst_name(ispec))),'num').ne.0) !if we need to print only aerosol numbers

       if(mass_or_num) then !ignore the aerosol number
          write(unit_output,12,advance="no")field(ispec),','
       endif

    enddo

    write(unit_output,'(A)')'],]'

  end subroutine write_output_aerosol_mmr_from_stateq_logical

  !================================================================================

  subroutine write_var_int(unit_input, unit_output, fld_name,field)
    !------------------------------------------------------------------
    !Purpose: Writes a 1D input and output field in a YAML file format
    !for a given column
    !------------------------------------------------------------------
    implicit none

    integer, intent(in)   :: unit_input      ! input stream unit number
    integer, intent(in)   :: unit_output     ! output stream unit number
    character(len=*), intent(in) :: fld_name ! name of the field
    integer, intent(in)  :: field        ! field values in int

    !check if file is open to write or not
    call is_file_open(unit_input)

    write(unit_input,'(3A,I10,A)')'    ',trim(adjustl(fld_name)),': [', field,']'

    call write_output_var_int(unit_output, fld_name, field, "input")

  end subroutine write_var_int

  !================================================================================

  subroutine write_output_var_int(unit_output, fld_name, field, inp_out_str)
    !------------------------------------------------------------------
    !Purpose: Writes a 1D output field in a YAML file format
    !for a given column
    !------------------------------------------------------------------
    implicit none

    integer, intent(in)   :: unit_output     ! output stream unit number
    character(len=*), intent(in) :: fld_name ! name of the field
    integer, intent(in)  :: field           ! field values in int

    !optional input
    character(len=*), intent(in), optional :: inp_out_str ! input or output

    !local
    character(len=20) :: object

    !check if file is open to write or not
    call is_file_open(unit_output)

    object = "output"
    if (present(inp_out_str)) then
       object = trim(adjustl(inp_out_str))
    endif

    write(unit_output,'(4A,I10,A)')trim(adjustl(object)),'.',trim(adjustl(fld_name)),'=[[', field,'],]'

  end subroutine write_output_var_int

  !================================================================================

  subroutine write_var_logical(unit_input, unit_output, fld_name,field)
    !------------------------------------------------------------------
    !Purpose: Writes a 1D input and output field in a YAML file format
    !for a given column
    !------------------------------------------------------------------
    implicit none

    integer, intent(in)   :: unit_input      ! input stream unit number
    integer, intent(in)   :: unit_output     ! output stream unit number
    character(len=*), intent(in) :: fld_name ! name of the field
    logical, intent(in)  :: field        ! field values in logical

    !check if file is open to write or not
    call is_file_open(unit_input)

    write(unit_input,'(3A,L10,A)')'    ',trim(adjustl(fld_name)),': [', field,']'

    call write_output_var_logical(unit_output, fld_name, field, "input")

  end subroutine write_var_logical

  !================================================================================

  subroutine write_output_var_logical(unit_output, fld_name, field, inp_out_str)
    !------------------------------------------------------------------
    !Purpose: Writes a 1D output field in a YAML file format
    !for a given column
    !------------------------------------------------------------------
    implicit none

    integer, intent(in)   :: unit_output     ! output stream unit number
    character(len=*), intent(in) :: fld_name ! name of the field
    logical, intent(in)  :: field           ! field values in int

    !optional input
    character(len=*), intent(in), optional :: inp_out_str ! input or output

    !local
    character(len=20) :: object

    !check if file is open to write or not
    call is_file_open(unit_output)

    object = "output"
    if (present(inp_out_str)) then
       object = trim(adjustl(inp_out_str))
    endif

    write(unit_output,'(4A,L10,A)')trim(adjustl(object)),'.',trim(adjustl(fld_name)),'=[[', field,'],]'

  end subroutine write_output_var_logical

  !================================================================================

  subroutine write_var_real(unit_input, unit_output, fld_name,field)
    !------------------------------------------------------------------
    !Purpose: Writes a 1D input and output field in a YAML file format
    !for a given column
    !------------------------------------------------------------------
    implicit none

    integer, intent(in)   :: unit_input      ! input stream unit number
    integer, intent(in)   :: unit_output     ! output stream unit number
    character(len=*), intent(in) :: fld_name ! name of the field
    real(r8), intent(in)  :: field        ! field values in r8

    !check if file is open to write or not
    call is_file_open(unit_input)

    write(unit_input,'(3A,E17.10,A)')'    ',trim(adjustl(fld_name)),': [', field,']'

    call write_output_var_real(unit_output, fld_name, field, "input")

  end subroutine write_var_real

  !================================================================================

  subroutine write_output_var_real(unit_output, fld_name, field, inp_out_str)
    !------------------------------------------------------------------
    !Purpose: Writes a 1D output field in a YAML file format
    !for a given column
    !------------------------------------------------------------------
    implicit none

    integer, intent(in)   :: unit_output     ! output stream unit number
    character(len=*), intent(in) :: fld_name ! name of the field
    real(r8), intent(in)  :: field           ! field values in r8

    !optional input
    character(len=*), intent(in), optional :: inp_out_str ! input or output

    !local
    character(len=20) :: object

    !check if file is open to write or not
    call is_file_open(unit_output)

    object = "output"
    if (present(inp_out_str)) then
       object = trim(adjustl(inp_out_str))
    endif

    write(unit_output,'(4A,E17.10,A)')trim(adjustl(object)),'.',trim(adjustl(fld_name)),'=[[', field,'],]'

  end subroutine write_output_var_real

  !================================================================================

  subroutine write_1d_var(unit_input, unit_output, fld_name,dim,field)
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

    write(unit_input,'(3A)',advance="no")'    ',trim(adjustl(fld_name)),': ['

    write(unit_input,10,advance="no")field(1)

    do k = 2, dim
       write(unit_input,11,advance="no")',',field(k)
    enddo

    write(unit_input,'(A)')']'

    call write_1d_output_var(unit_output, fld_name, dim, field, "input")

  end subroutine write_1d_var

  !================================================================================

  subroutine write_1d_output_var(unit_output, fld_name, dim, field, inp_out_str)
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
    !Since we capture for input and output of a subroutine in python format, this
    !subroutine is called for writing both the input and the output
    !Do not provide this optional argument if it is called from
    !an external subroutine external to this file
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

    write(unit_output,'(4A)',advance="no")trim(adjustl(object)),'.',trim(adjustl(fld_name)),'=[['

    do k = 1, dim
       write(unit_output,12,advance="no"),field(k),','
    enddo

    write(unit_output,'(A)')'],]'

  end subroutine write_1d_output_var

  !================================================================================

  subroutine write_2d_var(unit_input, unit_output, fld_name, dim1, dim2, field)
    !------------------------------------------------------------------
    !Purpose: Writes a 2D input and output field in a YAML file format
    !for a given column
    !------------------------------------------------------------------
    implicit none

    integer, intent(in)   :: unit_input      ! input stream unit number
    integer, intent(in)   :: unit_output     ! output stream unit number
    integer, intent(in)   :: dim1, dim2      ! dimensions of the field
    character(len=*), intent(in) :: fld_name ! name of the field
    real(r8), intent(in)  :: field(:,:)        ! field values in r8

    !local
    integer :: d1, d2

    !check if file is open to write or not
    call is_file_open(unit_input)

    !format statement to write in double precision
10  format(E17.10)
11  format(A,E17.10)

    write(unit_input,'(3A)',advance="no")'    ',trim(adjustl(fld_name)),': ['

    ! For maintaining format in the YAML inout file we have to  print first element of the array first
    write(unit_input,10,advance="no")field(1,1)

    !print rest of the column for d2=1
    d2 = 1
    do d1 = 2, dim1 !first element is already printed, start from the 2nd
       write(unit_input,11,advance="no")',',field(d1,d2)
    enddo

    do d2 = 2, dim2 !First column is already printed, start from the 2nd
       do d1 = 1, dim1
          write(unit_input,11,advance="no")',',field(d1,d2)
       enddo
    enddo

    write(unit_input,'(A)')']'

    call write_2d_output_var(unit_output, fld_name, dim1, dim2, field, "input")

  end subroutine write_2d_var

  !================================================================================

  subroutine write_2d_output_var(unit_output, fld_name, dim1, dim2, field, inp_out_str)
    !------------------------------------------------------------------
    !Purpose: Writes a 2D output field in a YAML file format
    !for a given column
    !------------------------------------------------------------------
    implicit none

    integer, intent(in)   :: unit_output     ! output stream unit number
    character(len=*), intent(in) :: fld_name ! name of the field
    integer, intent(in)   :: dim1, dim2      ! dimensions of the field
    real(r8), intent(in)  :: field(:,:)        ! field values in r8

    !optional input
    !Since we capture for input and output of a subroutine in python format, this
    !subroutine is called for writing both the input and the output
    !Do not provide this optional argument if it is called from
    !an external subroutine external to this file
    character(len=*), intent(in), optional :: inp_out_str ! input or output

    !local
    integer :: d1, d2
    character(len=20) :: object

    !check if file is open to write or not
    call is_file_open(unit_output)

    !format statement to write in double precision
12  format(E17.10,A)

    object = "output"
    if (present(inp_out_str)) then
       object = trim(adjustl(inp_out_str))
    endif

    write(unit_output,'(4A)',advance="no")trim(adjustl(object)),'.',trim(adjustl(fld_name)),'=[['

    do d2 = 1, dim2
       do d1 = 1, dim1
          write(unit_output,12,advance="no"),field(d1,d2),','
       enddo
    enddo

    write(unit_output,'(A)')'],]'

  end subroutine write_2d_output_var

end module yaml_input_file_io
