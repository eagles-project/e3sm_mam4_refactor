#ifdef YAML_AMIC2
  !-----------------------------------------------------------------------------------------
  !"lchnk" is needed for the following code to work,
  ! temporarily pass it along from upper level subroutines
  ! as y_lchnk and uncomment the following code:
  ! integer, intent(in) :: lchnk
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
  yaml%lev_print = 64       !level
  yaml%nstep_print = 379 !time step

  yaml%col_print = icolprnt(lchnk)                !column to write data

  !current_time step
  y_nstep = get_nstep()

  !Flag to decide to write or not to write data
  yaml%flag_print = .false. !(**remove these if generating data for a dependent subroutines**)

  !!!if(yaml%col_print == y_i .and. y_nstep==yaml%nstep_print .and. y_k == yaml%lev_print) then ! if this column exists in y_lchnk

  !-----------------------------------------------------------------------------------------
  !In the case of y_i or y_k are not passed as arguments, use the following if condition:
  if(yaml%col_print >0 .and. y_nstep==yaml%nstep_print) then
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
        call open_files('modal_aero_amicphys_intr', &  !intent-in
             unit_input, unit_output) !intent-out
        !    unit_input, unit_output, trim(ext_str)) !intent-out, with the use of ext_str


        !start by adding an input string
        call write_input_output_header(unit_input, unit_output,yaml%lchnk_print,yaml%col_print, &
             'modal_aero_amicphys_intr',yaml%nstep_print, yaml%lev_print)

        ! add code for writing data here
        
        call write_var(unit_input,unit_output,'mdo_gasaerexch',mdo_gasaerexch)
        call write_var(unit_input,unit_output,'mdo_rename',mdo_rename)
        call write_var(unit_input,unit_output,'mdo_newnuc',mdo_newnuc)
        call write_var(unit_input,unit_output,'mdo_coag',mdo_coag)
        call write_var(unit_input,unit_output,'lchnk',lchnk)
        call write_var(unit_input,unit_output,'ncol',ncol)
        call write_var(unit_input,unit_output,'nstep',nstep)
        call write_var(unit_input,unit_output,'loffset',loffset)
        call write_var(unit_input,unit_output,'deltat',deltat)
        call write_var(unit_input,unit_output,'latndx',latndx(yaml%col_print))
        call write_var(unit_input,unit_output,'lonndx',lonndx(yaml%col_print))
        call write_var(unit_input,unit_output,'temp',temp(yaml%col_print,yaml%lev_print))
        call write_var(unit_input,unit_output,'pmid',pmid(yaml%col_print,yaml%lev_print))
        call write_var(unit_input,unit_output,'pdel',pdel(yaml%col_print,yaml%lev_print))
        call write_var(unit_input,unit_output,'zm',zm(yaml%col_print,yaml%lev_print))
        call write_var(unit_input,unit_output,'pblh',pblh(yaml%col_print))
        call write_var(unit_input,unit_output,'qv',qv(yaml%col_print,yaml%lev_print))
        call write_var(unit_input,unit_output,'cld',cld(yaml%col_print,yaml%lev_print))
        call write_var(unit_input,unit_output,'qq',qq(yaml%col_print,yaml%lev_print,:))
        call write_var(unit_input,unit_output,'qqcw',qqcw(yaml%col_print,yaml%lev_print,:))
        call write_var(unit_input,unit_output,'q_pregaschem',q_pregaschem(yaml%col_print,yaml%lev_print,:))
        call write_var(unit_input,unit_output,'q_precldchem',q_precldchem(yaml%col_print,yaml%lev_print,:))
        call write_var(unit_input,unit_output,'qqcw_precldchem',qqcw_precldchem(yaml%col_print,yaml%lev_print,:))
        call write_var(unit_input,unit_output,'dgncur_a',dgncur_a(yaml%col_print,yaml%lev_print,:))
        call write_var(unit_input,unit_output,'dgncur_awet',dgncur_awet(yaml%col_print,yaml%lev_print,:))
        call write_var(unit_input,unit_output,'wetdens_host',wetdens_host(yaml%col_print,yaml%lev_print,:))

        ! below are optional input arguments to modal_aero_amicphys_intr
        if(present(qaerwat)) then
           call write_var(unit_input,unit_output,'qaerwat',qaerwat(yaml%col_print,yaml%lev_print,:))
        endif

        ! below are external module inputs to modal_aero_amicphys_intr
        call write_var(unit_input,unit_output,'pcols',pcols)
        call write_var(unit_input,unit_output,'pver',pver)
        call write_var(unit_input,unit_output,'gas_pcnst',gas_pcnst)
        call write_var(unit_input,unit_output,'ntot_amode',ntot_amode)
        call write_var(unit_input,unit_output,'max_mode',max_mode)
        call write_var(unit_input,unit_output,'ntot_amode_extd',ntot_amode_extd)
        call write_var(unit_input,unit_output,'maxsubarea',maxsubarea)
        call write_var(unit_input,unit_output,'nqtendaa',nqtendaa)
        call write_var(unit_input,unit_output,'nqqcwtendaa',nqqcwtendaa)
        call write_var(unit_input,unit_output,'top_lev',top_lev)     
        call write_var(unit_input,unit_output,'gravit',gravit)     


        !writes aerosol mmr from state%q or q vector (cloud borne and interstitial)
        !"aer_num_only" is .ture. if printing aerosol num only
        !call write_aerosol_mmr_from_stateq(unit_input, unit_output, fld_name,field,aer_num_only)

        !close only the input file, not the output file
        close(unit_input)
        call freeunit(unit_input)
     endif
  endif
#endif
