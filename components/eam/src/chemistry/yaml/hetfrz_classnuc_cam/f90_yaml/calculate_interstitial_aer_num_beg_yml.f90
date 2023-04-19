#ifdef YAML_HETFRZ_CLASSNUC
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

  ! character(len=200) :: ext_str  ! this is used when multiple sets of yaml output is needed to cover different options (e.g., true and false)
  type(yaml_vars) :: yaml
  integer  :: unit_input, unit_output, y_nstep

  if(yaml%flag_print) then

        !open I/O yaml files (it can have an extra optional argument to pass a unique string to differentiate file names)
        call open_files('calculate_interstitial_aer_num', &  !intent-in
             unit_input, unit_output) !intent-out


        !start by adding an input string
        call write_input_output_header(unit_input, unit_output,yaml%lchnk_print,yaml%col_print, &
             'calculate_interstitial_aer_num',yaml%nstep_print, yaml%lev_print)

        !< add code for writing data here>
        call write_var(unit_input,unit_output,'ncnst',ncnst)
        call write_var(unit_input,unit_output,'aer',aer)


        !close only the input file, not the output file
        close(unit_input)
        call freeunit(unit_input)
     endif
#endif
