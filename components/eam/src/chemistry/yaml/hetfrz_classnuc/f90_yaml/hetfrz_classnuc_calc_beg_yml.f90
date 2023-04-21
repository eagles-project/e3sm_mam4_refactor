#ifdef YAML_HETFRZ_CLASSNUC
  !-----------------------------------------------------------------------------------------
  !"lchnk" is needed for the following code to work,
  ! temporarily pass it along from upper level subroutines
  ! as y_lchnk and uncomment the following code:
  ! integer, intent(in) :: y_lchnk
  !-----------------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------------------
  ! **OR** if this subroutine is called in a nested loop of columns and levels,
  ! we also might need column index (y_i or icol) and level
  ! index (y_k or klev) to be passed to this routine and
  ! uncomment the following code:
  ! integer, intent(in) :: y_i, y_k, y_lchnk
  !-----------------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------------------
  ! This is used when multiple sets of yaml output is needed
  !to cover different options (e.g., true and false)
  character(len=200) :: ext_str
  !-----------------------------------------------------------------------------------------

  type(yaml_vars) :: yaml
  integer  :: unit_input, unit_output, y_nstep, y_k

  ! some subroutines are called multiple times in one timestep, record the number of calls
  integer,save :: n_calls=0


  !-----------------------------------------------------------------------------------------
  !For generating data for a dependent subroutines where "yaml" derived type is already
  !initialized, use the following if condition
  if(yaml%flag_print .and. mod(y_k,10) == 0 ) then
     !-----------------------------------------------------------------------------------------

     !-----------------------------------------------------------------------------------------
     ! Set "ext_str" if there are multiple sets of yaml output to be written out
     ! Example:"flag" in the code can be 0, 1, or 2, we can update "ext_str" as:
     write(ext_str,'(I2)') y_k
     ext_str = 'lev_'//adjustl(ext_str)
     yaml%lev_print = y_k
     !-----------------------------------------------------------------------------------------


        !open I/O yaml files
        !(with an optional argument to pass a unique string to differentiate file names)
        call open_files('hetfrz_classnuc_calc', &  !intent-in
             unit_input, unit_output) !intent-out
        !    unit_input, unit_output, trim(ext_str)) !intent-out, with the use of ext_str


        !start by adding an input string
        call write_input_output_header(unit_input, unit_output,yaml%lchnk_print,yaml%col_print, &
             'hetfrz_classnuc_calc',yaml%nstep_print, yaml%lev_print)

        ! add code for writing data here

        call write_var(unit_input,unit_output,'deltat',deltat)
        call write_var(unit_input,unit_output,'temperature',temperature)
        call write_var(unit_input,unit_output,'pressure',pressure)
        call write_var(unit_input,unit_output,'supersatice',supersatice)
        call write_var(unit_input,unit_output,'fn',fn)
        call write_var(unit_input,unit_output,'r3lx',r3lx)
        call write_var(unit_input,unit_output,'icnlx',icnlx)
        call write_var(unit_input,unit_output,'hetraer',hetraer)
        call write_var(unit_input,unit_output,'awcam',awcam)
        call write_var(unit_input,unit_output,'awfacm',awfacm)
        call write_var(unit_input,unit_output,'dstcoat',dstcoat)
        call write_var(unit_input,unit_output,'total_aer_num',total_aer_num)
        call write_var(unit_input,unit_output,'coated_aer_num',coated_aer_num)
        call write_var(unit_input,unit_output,'uncoated_aer_num',uncoated_aer_num)
        call write_var(unit_input,unit_output,'total_interstitial_aer_num',total_interstitial_aer_num)
        call write_var(unit_input,unit_output,'total_cloudborne_aer_num',total_cloudborne_aer_num)


        !writes aerosol mmr from state%q or q vector (cloud borne and interstitial)
        !"aer_num_only" is .ture. if printing aerosol num only
        !call write_aerosol_mmr_from_stateq(unit_input, unit_output, fld_name,field,aer_num_only)

        !close only the input file, not the output file
        close(unit_input)
        call freeunit(unit_input)
  endif
#endif
