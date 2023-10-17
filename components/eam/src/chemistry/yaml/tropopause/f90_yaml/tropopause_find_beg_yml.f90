#ifdef YAML_TROPFIND
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
  integer  :: unit_input, unit_output, y_nstep
 !  print output from timestep nstep_print_lo to yaml%nstep_print
  integer  :: nstep_print_lo

  ! some subroutines are called multiple times in one timestep, record the number of calls
  integer,save :: n_calls=0
  integer,save :: y_nstep_old=0
  integer      :: itype
  logical,save :: hadtype(0:6) = .false.


  !populate YAML structure
  !(**remove yaml%lev_print, nstep_print, col_print if generating data for a dependent subroutines**)
  ! used with phys_debug_lat = 39.0553303768561, phys_debug_lon = 262.904774388088

  yaml%lev_print = 30       !print level above tropopause
  yaml%nstep_print = 1417 !time step
  nstep_print_lo = 1400
  yaml%col_print = icolprnt(lchnk)                !column to write data

  ! below are levels / timeranges needed to test tropopause_climate option for the above column

  !yaml%lev_print = 54       !print level above tropopause
  !yaml%nstep_print = 574 !time step
  !nstep_print_lo = 552
  !yaml%col_print = icolprnt(lchnk)                !column to write data

  !current_time step
  y_nstep = get_nstep()

  !Flag to decide to write or not to write data
  yaml%flag_print = .false. !(**remove these if generating data for a dependent subroutines**)

 !! if(yaml%col_print == y_i .and. y_nstep==yaml%nstep_print .and. y_k == yaml%lev_print) then ! if this column exists in y_lchnk
  !  print output for timesteps between nstep_print_lo and yaml%nstep_print
  !  output depends on multiple vertical levels with internal kk loop, so use below condition
  if(yaml%col_print >0 .and. y_nstep<=yaml%nstep_print .and. y_nstep>=nstep_print_lo) then
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
     !
     !  itype = 0, not called from any of the following (should not occur)  
     !        = 1, called from mozart/chemistry.F90
     !        = 2, called from utils/prescribed_volcaero.F90
     !        = 3, third call from cam/tropopause_f90 (tropopause_output subroutine)
     !        = 4, second call from cam/tropopause_f90 (tropopause_output subroutine)
     !        = 5, first call from cam/tropopause_f90 (tropopause_output subroutine)
     !        = 6, called from cam/aer_rad_props.F90

     itype = 0
     if(present(primary)) then
        if(primary == TROP_ALG_HYBSTOB) then
          itype = 1
        elseif(primary == TROP_ALG_TWMO) then
          itype = 2
        elseif(primary == TROP_ALG_CPP) then
          itype = 3
        endif
     elseif(present(backup)) then
        if(backup == TROP_ALG_NONE) then
           itype = 4
        endif
     elseif(present(tropP)) then
        itype = 5
     else
        itype = 6
     endif
     write(ext_str,'(A5,I1)') 'itype',itype

  
     !-----------------------------------------------------------------------------------------

     !Record number of calls that can output yaml file if you only need to write one set of input/output
     !Increment n_calls if current timestep, y_nstep, is same as for previous pass through this code.
     !Otherwise, we are at a new timestep, can reset n_calls to 1.

     if(y_nstep /= y_nstep_old) then
         hadtype(:) = .false.
     endif
     y_nstep_old = y_nstep

     if(.not.(hadtype(itype))) then
        hadtype(itype) = .true.

     !if (n_calls==1) then ! output at the first call only, modify this if condition (see below) if writing out other calls

     !if ((n_calls==1 .and. flag==0) .or. (n_calls==3 .and. flag==1) .or. (n_calls==5 .and. flag==2)) then

        !(**remove these yaml% variables if generating data for a dependent subroutines**)
        yaml%lchnk_print = lchnk
        yaml%flag_print  = .true.


        !open I/O yaml files
        !(with an optional argument to pass a unique string to differentiate file names)
        call open_files('tropopause_find', &  !intent-in
        !     unit_input, unit_output) !intent-out
            unit_input, unit_output, trim(ext_str)) !intent-out, with the use of ext_str


        !start by adding an input string
        !Use current timestep, given by y_nstep, to form name of input / output files
        call write_input_output_header(unit_input, unit_output,yaml%lchnk_print,yaml%col_print, &
             'tropopause_find',y_nstep, yaml%lev_print)

        ! add code for writing data here
        ! below are explicit argument inputs to tropopause_find
        
        call write_var(unit_input,unit_output,'lchnk',lchnk)
        call write_var(unit_input,unit_output,'ncol',ncol)
        call write_var(unit_input,unit_output,'pmid',pmid(yaml%col_print,:))
        call write_var(unit_input,unit_output,'pint',pint(yaml%col_print,:))
        call write_var(unit_input,unit_output,'temp',temp(yaml%col_print,:))
        call write_var(unit_input,unit_output,'zm',zm(yaml%col_print,:))
        call write_var(unit_input,unit_output,'zi',zi(yaml%col_print,:))
        
        if(present(primary)) then
           call write_var(unit_input,unit_output,'primary',primary)
        endif
        if(present(backup)) then
           call write_var(unit_input,unit_output,'backup',backup)
        endif

        ! below are external module inputs to tropopause_find
 
        call write_var(unit_input,unit_output,'pcols',pcols)
        call write_var(unit_input,unit_output,'fillvalue',fillvalue)

        ! below are internal module inputs to tropopause_find
        
        call write_var(unit_input,unit_output,'NOTFOUND',NOTFOUND)
        call write_var(unit_input,unit_output,'TROP_ALG_NONE',TROP_ALG_NONE)
        call write_var(unit_input,unit_output,'default_primary',default_primary)
        call write_var(unit_input,unit_output,'default_backup',default_backup)

        ! below are internal module fields used in this data writing code

        call write_var(unit_input,unit_output,'TROP_ALG_TWMO',TROP_ALG_TWMO)
        call write_var(unit_input,unit_output,'TROP_ALG_HYBSTOB',TROP_ALG_HYBSTOB)
        call write_var(unit_input,unit_output,'TROP_ALG_CPP',TROP_ALG_CPP)



        !writes aerosol mmr from state%q or q vector (cloud borne and interstitial)
        !"aer_num_only" is .ture. if printing aerosol num only
        !call write_aerosol_mmr_from_stateq(unit_input, unit_output, fld_name,field,aer_num_only)

        !close only the input file, not the output file
        close(unit_input)
        call freeunit(unit_input)
     endif
  endif
#endif
