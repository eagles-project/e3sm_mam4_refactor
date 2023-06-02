#ifdef YAML_NDROP
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

 ! This within nested column loop but not level loop, hence below
  integer, intent(in) :: y_i, y_lchnk
  !-----------------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------------------
  ! This is used when multiple sets of yaml output is needed
  !to cover different options (e.g., true and false)
  ! character(len=200) :: ext_str
  !-----------------------------------------------------------------------------------------

  type(yaml_vars) :: yaml
  integer  :: unit_input, unit_output, y_nstep
   !  print output from timestep nstep_print_lo to yaml%nstep_print
  integer  :: nstep_print_lo
  ! some subroutines are called multiple times in one timestep, record the number of calls
!  integer,save :: n_calls=0
  integer,save :: y_nstep_old = 0  

  !populate YAML structure
  !(**remove yaml%lev_print, nstep_print, col_print if generating data for a dependent subroutines**)
 ! used with phys_debug_lat = 39.0553303768561, phys_debug_lon = 262.904774388088

  yaml%lev_print = 54       !level
  yaml%nstep_print = 1417 !time step
  nstep_print_lo = 1400
  yaml%col_print = icolprnt(y_lchnk)                !column to write data

  !current_time step
  y_nstep = get_nstep()

  !Flag to decide to write or not to write data
  yaml%flag_print = .false. !(**remove these if generating data for a dependent subroutines**)

!!  if(yaml%col_print == y_i .and. y_nstep==yaml%nstep_print .and. y_k == yaml%lev_print) then ! if this column exists in y_lchnk
  !  print output for timesteps between nstep_print_lo and yaml%nstep_print
  if(yaml%col_print == y_i .and. y_nstep<=yaml%nstep_print .and. y_nstep>=nstep_print_lo) then

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
!     n_calls = n_calls+1

!     if (n_calls==1) then ! output at the first call only, modify this if condition (see below) if writing out other calls

     !if ((n_calls==1 .and. flag==0) .or. (n_calls==3 .and. flag==1) .or. (n_calls==5 .and. flag==2)) then

     !We only want to print the first data occurrence within a timestep.
     !If current timestep (y_nstep) is the same as the timestep during the last pass of this code segment
     !(y_nstep_old), then we do not print
     !If y_nstep does not equal y_nstep_old, then we can reset all flags, and print next occurrence of everything.

     if(y_nstep /= y_nstep_old) then
        y_nstep_old = y_nstep

        !(**remove these yaml% variables if generating data for a dependent subroutines**)
        yaml%lchnk_print = y_lchnk
        yaml%flag_print  = .true.


        !open I/O yaml files
        !(with an optional argument to pass a unique string to differentiate file names)
        call open_files('update_from_explmix', &  !intent-in
             unit_input, unit_output) !intent-out
        !    unit_input, unit_output, trim(ext_str)) !intent-out, with the use of ext_str


        !start by adding an input string
        !Use current timestep, given by y_nstep, to form name of input / output files
        call write_input_output_header(unit_input, unit_output,yaml%lchnk_print,yaml%col_print, &
             'update_from_explmix',y_nstep, yaml%lev_print)

        ! add code for writing data here
        ! below are explicit input arguments to update_from_explmix
        ! subroutine operates on multiple vertical levels, so choose to print all vertical levels for given column
       
        call write_var(unit_input,unit_output,'dtmicro',dtmicro)
        call write_var(unit_input,unit_output,'csbot',csbot(:))
        call write_var(unit_input,unit_output,'cldn_col',cldn_col(:))
        call write_var(unit_input,unit_output,'zn',zn(:))
        call write_var(unit_input,unit_output,'zs',zs(:))
        call write_var(unit_input,unit_output,'ekd',ekd(:))

        call write_var(unit_input,unit_output,'nact',nact(:,:))
        call write_var(unit_input,unit_output,'mact',mact(:,:))
        call write_var(unit_input,unit_output,'qcld',qcld(:))
        call write_var(unit_input,unit_output,'nnew',nnew)
        call write_var(unit_input,unit_output,'nsav',nsav)

        ! below, do separate prints for two time levels of raercol and raercol_cw, otherwise 3d variables
        ! nnew is index of time level being updated, nsav is index of time level with previous update

        call write_var(unit_input,unit_output,'raercol_1',raercol(:,:,1))
        call write_var(unit_input,unit_output,'raercol_cw_1',raercol_cw(:,:,1))
        call write_var(unit_input,unit_output,'raercol_2',raercol(:,:,2))
        call write_var(unit_input,unit_output,'raercol_cw_2',raercol_cw(:,:,2))

        ! below are external module variable inputs to update_from_explmix

        call write_var(unit_input,unit_output,'pver',pver)
        call write_var(unit_input,unit_output,'top_lev',top_lev)

        ! below are ndrop module variable inputs to update_from_explmix

        call write_var(unit_input,unit_output,'ntot_amode',ntot_amode)
        call write_var(unit_input,unit_output,'mam_idx',mam_idx)
        call write_var(unit_input,unit_output,'nspec_amode',nspec_amode)

        !writes aerosol mmr from state%q or q vector (cloud borne and interstitial)
        !"aer_num_only" is .ture. if printing aerosol num only
        !call write_aerosol_mmr_from_stateq(unit_input, unit_output, fld_name,field,aer_num_only)

        !close only the input file, not the output file
        close(unit_input)
        call freeunit(unit_input)
     endif
  endif
#endif
