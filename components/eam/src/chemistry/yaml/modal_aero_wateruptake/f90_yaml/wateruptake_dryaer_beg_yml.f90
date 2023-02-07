#ifdef YAML_WUPTAKE
  !<"lchnk" is needed for the following code to work,
  ! temporarily pass it along from upper level subroutines
  ! as y_lchnk and uncomment the following code:
  integer, intent(in) :: y_lchnk

  ! **OR** if this subroutine is called in a nested loop of columns and levels,
  ! we also might need column index (y_i or icol) and level
  ! index (y_k or klev) to be passed to this routine and
  ! uncomment the following code:
  ! integer, intent(in) :: y_i, y_k, y_lchnk

  !>

  ! character(len=200) :: ext_str  ! this is used when multiple sets of yaml output is needed to cover different options (e.g., true and false)
  type(yaml_vars) :: yaml
  integer  :: unit_input, unit_output, y_nstep
  integer,save :: n_calls=0   ! some subroutines are called multiple times in one timestep, record the number of calls

  !populate YAML structure
  yaml%lev_print = 51 !level (**remove these if generating data for a dependent subroutines**)
  yaml%nstep_print = 355 !time step(**remove these if generating data for a dependent subroutines**)

  ! set ext_str if there are multiple sets of yaml output to be write out
  ! here gives an example that "flag" in the code can be 0, 1, or 2:
  ! write(ext_str,'(I2)') flag
  ! ext_str = 'flag_'//adjustl(ext_str)


  !YAML file input generation code- DO NOT PORT to C++
  !print all inputs one-by-one at column "yaml%col_print"
  yaml%col_print = icolprnt(y_lchnk) !column to write data(**remove these if generating data for a dependent subroutines**)
  y_nstep = get_nstep() !time step (**remove these if generating data for a dependent subroutines**)

  yaml%flag_print = .false. ! to write or not to write data (**remove these if generating data for a dependent subroutines**)
!  if(yaml%col_print == y_i .and. y_nstep==yaml%nstep_print .and. y_k == yaml%lev_print) then ! if this column exists in y_lchnk
     !<
     !In the case of y_i or y_k are not passed as arguments:
     if(yaml%col_print >0 .and. y_nstep==yaml%nstep_print) then

     !For generating data for a dependent subroutines where "yaml" derived type is already initialized:
     !if(yaml%flag_print) then
     !>

     ! record number of calls that can output yaml file. you only need to write one set of input/output
     n_calls = n_calls+1
     if (n_calls==1) then ! output the first call. change this if writing out other calls
     ! if ((n_calls==1 .and. flag==0) .or. (n_calls==3 .and. flag==1) .or. (n_calls==5 .and. flag==2)) then ! this is an example of writing multiple sets of YAML files

        yaml%lchnk_print = y_lchnk !(**remove these if generating data for a dependent subroutines**)
        yaml%flag_print  = .true.!(**remove these if generating data for a dependent subroutines**)

        !open I/O yaml files (it can have an extra optional argument to pass a unique string to differentiate file names)
        call open_files('wateruptake_dryaer', &  !intent-in
             unit_input, unit_output) !intent-out
!             unit_input, unit_output, trim(ext_str)) !intent-out, with the use of ext_str


        !start by adding an input string
        call write_input_output_header(unit_input, unit_output,yaml%lchnk_print,yaml%col_print, &
             'compute_tendencies',yaml%nstep_print, yaml%lev_print)

        !< add code for writing data here>
        call write_var(unit_input, unit_output, 'ncol', ncol)
        call write_1d_var(unit_input, unit_output, 'dgncur_a', nmodes,dgncur_a(yaml%col_print, yaml%lev_print,:)) 
        call write_1d_var(unit_input, unit_output, 'state_q', nmodes, state_q(yaml%col_print, yaml%lev_print,:))

        !call write_var(unit_input, unit_output, fld_name,field)!write a single variable
        !call write_1d_var(unit_input, unit_output, fld_name,dim,field) ! writes 1D variables of any dimension
        !call write_2d_var(unit_input, unit_output, fld_name, dim1, dim2, field) ! writes 2D variables of any dimension: field(dim1,dim2)

        !writes aerosol mmr from state%q or q vector (cloud borne and interstitial)
        !"aer_num_only" is .ture. if printing aerosol num only
        !call write_aerosol_mmr_from_stateq(unit_input, unit_output, fld_name,field,aer_num_only)

        !close only the input file, not the output file
        close(unit_input)
        call freeunit(unit_input)
     endif
  endif
#endif
