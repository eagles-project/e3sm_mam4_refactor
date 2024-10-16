#ifdef YAML_GASAEREXCH
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

  type(yaml_vars) :: yaml
  integer  :: unit_input, unit_output, y_nstep

  ! some subroutines are called multiple times in one timestep, record the number of calls
  integer,save :: n_calls=0


  !populate YAML structure
  !(**remove yaml%lev_print, nstep_print, col_print if generating data for a dependent subroutines**)
  yaml%lev_print = k       !level
  yaml%nstep_print = 379 !time step

  yaml%col_print = icolprnt(lchnk)                !column to write data

  !current_time step
  y_nstep = get_nstep()

  !Flag to decide to write or not to write data
  yaml%flag_print = .false. !(**remove these if generating data for a dependent subroutines**)

  if(yaml%col_print == i .and. y_nstep==yaml%nstep_print .and. k == yaml%lev_print) then ! if this column exists in y_lchnk

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
        yaml%lchnk_print = lchnk
        yaml%flag_print  = .true.


        !open I/O yaml files
        !(with an optional argument to pass a unique string to differentiate file names)
        call open_files('mam_gasaerexch_1subarea', &  !intent-in
             unit_input, unit_output) !intent-out
        !    unit_input, unit_output, trim(ext_str)) !intent-out, with the use of ext_str


        !start by adding an input string
        call write_input_output_header(unit_input, unit_output,yaml%lchnk_print,yaml%col_print, &
             'mam_gasaerexch_1subarea',yaml%nstep_print, yaml%lev_print)

        ! add code for writing data here
        
        call write_var(unit_input,unit_output,'nstep',nstep)
        call write_var(unit_input,unit_output,'lchnk',lchnk)
        call write_var(unit_input,unit_output,'i',i)
        call write_var(unit_input,unit_output,'k',k)
        call write_var(unit_input,unit_output,'jsub',jsub)
        call write_var(unit_input,unit_output,'jtsubstep',jtsubstep)
        call write_var(unit_input,unit_output,'ntsubstep',ntsubstep)
        call write_var(unit_input,unit_output,'latndx',latndx)
        call write_var(unit_input,unit_output,'lonndx',lonndx)
        call write_var(unit_input,unit_output,'lund',lund)
        call write_var(unit_input,unit_output,'n_mode',n_mode)
        call write_var(unit_input,unit_output,'dtsubstep',dtsubstep)
        call write_var(unit_input,unit_output,'temp',temp)
        call write_var(unit_input,unit_output,'pmid',pmid)
        call write_var(unit_input,unit_output,'aircon',aircon)

        call write_var(unit_input,unit_output,'qgas_cur',qgas_cur(:))
        call write_var(unit_input,unit_output,'qgas_avg',qgas_avg(:))
        call write_var(unit_input,unit_output,'qgas_netprod_otrproc',qgas_netprod_otrproc(:))
        call write_var(unit_input,unit_output,'qaer_cur',qaer_cur(:,:))
        call write_var(unit_input,unit_output,'qnum_cur',qnum_cur(:))
        call write_var(unit_input,unit_output,'qwtr_cur',qwtr_cur(:))
        call write_var(unit_input,unit_output,'dgn_a',dgn_a(:))
        call write_var(unit_input,unit_output,'dgn_awet',dgn_awet(:))
        call write_var(unit_input,unit_output,'wetdens',wetdens(:))
        call write_var(unit_input,unit_output,'uptkaer',uptkaer(:,:))
        call write_var(unit_input,unit_output,'uptkrate_h2so4',uptkrate_h2so4)

!  below are external module inputs to mam_gasaerexch_1subarea
        call write_var(unit_input,unit_output,'max_gas',max_gas)
        call write_var(unit_input,unit_output,'max_aer',max_aer)
        call write_var(unit_input,unit_output,'max_mode',max_mode)
        call write_var(unit_input,unit_output,'ntot_amode',ntot_amode)
        call write_var(unit_input,unit_output,'pver',pver)
        call write_var(unit_input,unit_output,'igas_nh3',igas_nh3)
        call write_var(unit_input,unit_output,'igas_h2so4',igas_h2so4)
        call write_var(unit_input,unit_output,'lmap_aer',lmap_aer(:,:))
        call write_var(unit_input,unit_output,'nsoa',nsoa)
        call write_var(unit_input,unit_output,'ngas',ngas)
        call write_var(unit_input,unit_output,'igas_hno3',igas_hno3)
        call write_var(unit_input,unit_output,'igas_hcl',igas_hcl)


        !writes aerosol mmr from state%q or q vector (cloud borne and interstitial)
        !"aer_num_only" is .ture. if printing aerosol num only
        !call write_aerosol_mmr_from_stateq(unit_input, unit_output, fld_name,field,aer_num_only)

        !close only the input file, not the output file
        close(unit_input)
        call freeunit(unit_input)
     endif
  endif
#endif
