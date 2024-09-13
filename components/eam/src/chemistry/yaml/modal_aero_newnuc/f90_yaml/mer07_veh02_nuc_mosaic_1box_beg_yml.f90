#ifdef YAML_VEHPBLNUC
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
   integer, intent(in) :: y_i, y_k, y_lchnk
  !-----------------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------------------
  ! This is used when multiple sets of yaml output is needed
  !to cover different options (e.g., true and false)
  ! character(len=200) :: ext_str
  !-----------------------------------------------------------------------------------------

  type(yaml_vars) :: yaml
  integer  :: unit_input, unit_output, y_nstep

  ! some subroutines are called multiple times in one timestep, record the number of calls
  integer,save :: n_calls=0


  !populate YAML structure
  !(**remove yaml%lev_print, nstep_print, col_print if generating data for a dependent subroutines**)
  yaml%lev_print = 68       !BJG level
  yaml%nstep_print = 379       !time step

  yaml%col_print = icolprnt(y_lchnk)                !column to write data

  !current_time step
  y_nstep = get_nstep()

  !Flag to decide to write or not to write data
  yaml%flag_print = .false. !(**remove these if generating data for a dependent subroutines**)

  if(yaml%col_print == y_i .and. y_nstep==yaml%nstep_print .and. y_k == yaml%lev_print) then ! if this column exists in y_lchnk

  !-----------------------------------------------------------------------------------------
  !In the case of y_i or y_k are not passed as arguments, use the following if condition:
  !if(yaml%col_print >0 .and. y_nstep==yaml%nstep_print) then
  !-----------------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------------------
  !For generating data for a dependent subroutines where "yaml" derived type is already
  !initialized, use the following if condition
  !if(yaml%flag_print) then
  !-----------------------------------------------------------------------------------------

     !-----------------------------------------------------------------------------------------
     ! Set "ext_str" if there are multiple sets of yaml output to be written out
     ! Example:"flag" in the code can be 0, 1, or 2, we can update "ext_str" as:
     ! write(ext_str,'(I2)') flag
     ! ext_str = 'flag_'//adjustl(ext_str)
     !-----------------------------------------------------------------------------------------


     !Record number of calls that can output yaml file if you only need to write one set of input/output
     n_calls = n_calls+1

     if (n_calls==1) then ! output at the first call only, modify this if condition (see below) if writing out other calls

     !if ((n_calls==1 .and. flag==0) .or. (n_calls==3 .and. flag==1) .or. (n_calls==5 .and. flag==2)) then

        !(**remove these yaml% variables if generating data for a dependent subroutines**)
        yaml%lchnk_print = y_lchnk
        yaml%flag_print  = .true.


        !open I/O yaml files
        !(with an optional argument to pass a unique string to differentiate file names)
        call open_files('mer07_veh02_nuc_mosaic_1box', &  !intent-in
             unit_input, unit_output) !intent-out
        !    unit_input, unit_output, trim(ext_str)) !intent-out, with the use of ext_str


        !start by adding an input string
        call write_input_output_header(unit_input, unit_output,yaml%lchnk_print,yaml%col_print, &
             'mer07_veh02_nuc_mosaic_1box',yaml%nstep_print, yaml%lev_print)

        ! add code for writing data here
        
        call write_var(unit_input,unit_output,'dtnuc',dtnuc)
        call write_var(unit_input,unit_output,'temp_in',temp_in)
        call write_var(unit_input,unit_output,'rh_in',rh_in)
        call write_var(unit_input,unit_output,'press_in',press_in)
        call write_var(unit_input,unit_output,'zm_in',zm_in)
        call write_var(unit_input,unit_output,'pblh_in',pblh_in)
        call write_var(unit_input,unit_output,'qh2so4_cur',qh2so4_cur)
        call write_var(unit_input,unit_output,'qh2so4_avg',qh2so4_avg)
        call write_var(unit_input,unit_output,'qnh3_cur',qnh3_cur)
        call write_var(unit_input,unit_output,'h2so4_uptkrate',h2so4_uptkrate)
        call write_var(unit_input,unit_output,'mw_so4a_host',mw_so4a_host)
        call write_var(unit_input,unit_output,'newnuc_method_flagaa',newnuc_method_flagaa)
        call write_var(unit_input,unit_output,'nsize',nsize)
        call write_var(unit_input,unit_output,'maxd_asize',maxd_asize)
        call write_var(unit_input,unit_output,'dplom_sect',dplom_sect(:))
        call write_var(unit_input,unit_output,'dphim_sect',dphim_sect(:))
        call write_var(unit_input,unit_output,'ldiagaa',ldiagaa)

        ! below are external module inputs to mer07_veh02_nuc_mosaic_1box

        call write_var(unit_input,unit_output,'rgas',rgas)
        call write_var(unit_input,unit_output,'avogad',avogad)
        call write_var(unit_input,unit_output,'mw_so4a',mw_so4a)
        call write_var(unit_input,unit_output,'mw_nh4a',mw_nh4a)

        !writes aerosol mmr from state%q or q vector (cloud borne and interstitial)
        !"aer_num_only" is .ture. if printing aerosol num only
        !call write_aerosol_mmr_from_stateq(unit_input, unit_output, fld_name,field,aer_num_only)

        !close only the input file, not the output file
        close(unit_input)
        call freeunit(unit_input)
     endif
  endif
#endif
