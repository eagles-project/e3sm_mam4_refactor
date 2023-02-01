#ifdef YAML_RENAME
  !<"lchnk" is needed for the following code to work,
  ! temporarily pass it along from upper level subroutines
  ! as y_lchnk and uncomment the following code:
  ! integer, intent(in) :: y_lchnk

  ! **OR** if this subroutine is called in a nested loop of columns and levels,
  ! we also might need column index (y_i or icol) and level
  ! index (y_k or klev) to be passed to this routine and
  ! uncomment the following code:
  ! integer, intent(in) :: y_i, y_k, y_lchnk

  !>

  type(yaml_vars), intent(in) :: yaml
  integer  :: unit_input, unit_output, y_nstep
  integer,save :: n_calls=0   ! some subroutines are called multiple times in one timestep, record the number of calls

     !For generating data for a dependent subroutines where "yaml" derived type is already initialized:
     if(yaml%flag_print) then
     n_calls = n_calls+1
     if (n_calls==1) then
     !>
     !open I/O yaml files (it can have an extra optional argument to pass a unique string to differentiate file names)
     call open_files('compute_dryvol_change_in_src_mode', &  !intent-in
          unit_input, unit_output) !intent-out

     !start by adding an input string
     call write_input_output_header(unit_input, unit_output,yaml%lchnk_print,yaml%col_print, &
          'compute_tendencies',yaml%nstep_print, yaml%lev_print)

     !< add code for writing data here>
call write_var(unit_input, unit_output, 'nmode', nmode)
call write_var(unit_input, unit_output, 'nspec', nspec)
call write_1d_var(unit_input, unit_output,'dest_mode_of_mode',max_mode,dest_mode_of_mode)
call write_2d_var(unit_input, unit_output,'q_mmr',max_aer,max_mode,q_mmr)
call write_2d_var(unit_input, unit_output,'q_del_growth',max_aer,max_mode,q_del_growth)

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
