#ifdef YAML_JLONG
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
  ! character(len=200) :: ext_str
  !-----------------------------------------------------------------------------------------

  type(yaml_vars), intent(in) :: yaml
  integer  :: unit_input, unit_output, y_nstep

  ! some subroutines are called multiple times in one timestep, record the number of calls
  integer,save :: n_calls=0

  real(r8),allocatable :: rsf_tab_1d(:)  ! change rsf_tab from r4 to r8 and dimension to 1d
  !-----------------------------------------------------------------------------------------
  !In the case of y_i or y_k are not passed as arguments, use the following if condition:
  !if(yaml%col_print >0 .and. y_nstep==yaml%nstep_print) then
  !-----------------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------------------
  !For generating data for a dependent subroutines where "yaml" derived type is already
  !initialized, use the following if condition
  if(yaml%flag_print) then
  !-----------------------------------------------------------------------------------------

     !-----------------------------------------------------------------------------------------
     ! Set "ext_str" if there are multiple sets of yaml output to be written out
     ! Example:"flag" in the code can be 0, 1, or 2, we can update "ext_str" as:
     ! write(ext_str,'(I2)') flag
     ! ext_str = 'flag_'//adjustl(ext_str)
     !-----------------------------------------------------------------------------------------


     n_calls = n_calls+1

     if (n_calls==1) then ! output at the first call only, modify this if condition (see below) if writing out other calls
     !if ((n_calls==1 .and. flag==0) .or. (n_calls==3 .and. flag==1) .or. (n_calls==5 .and. flag==2)) then


        !open I/O yaml files
        !(with an optional argument to pass a unique string to differentiate file names)
        call open_files('calc_sum_wght', &  !intent-in
             unit_input, unit_output) !intent-out
        !    unit_input, unit_output, trim(ext_str)) !intent-out, with the use of ext_str


        !start by adding an input string
        call write_input_output_header(unit_input, unit_output,yaml%lchnk_print,yaml%col_print, &
             'calc_sum_wght',yaml%nstep_print, yaml%lev_print)

        ! add code for writing data here
        
        call write_var(unit_input,unit_output,'dels',dels)
        call write_var(unit_input,unit_output,'wrk0',wrk0)
        call write_var(unit_input,unit_output,'iz',iz)
        call write_var(unit_input,unit_output,'is',is)
        call write_var(unit_input,unit_output,'iv',iv)
        call write_var(unit_input,unit_output,'ial',ial)
        call write_var(unit_input,unit_output,'shape_rsf_tab',[67,1,2,2,2])


allocate(rsf_tab_1d(67*2*2*2))  ! only write out part of the array to reduce size
       rsf_tab_1d(:) = reshape(rsf_tab(:,iz, is:isp1, iv:ivp1, ial:ialp1),&
                        (/size(rsf_tab_1d)/)) 
        call write_var(unit_input,unit_output,'rsf_tab',rsf_tab_1d)
deallocate(rsf_tab_1d)

        !writes aerosol mmr from state%q or q vector (cloud borne and interstitial)
        !"aer_num_only" is .ture. if printing aerosol num only
        !call write_aerosol_mmr_from_stateq(unit_input, unit_output, fld_name,field,aer_num_only)

        !close only the input file, not the output file
        close(unit_input)
        call freeunit(unit_input)
     endif
  endif
#endif
