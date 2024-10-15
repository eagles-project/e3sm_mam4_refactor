#ifdef YAML_AMICPHYS_1SUBAREA
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
  yaml%lev_print = 7       !level
  yaml%nstep_print = 379 !time step

  yaml%col_print = icolprnt(lchnk)                !column to write data

  !current_time step
  y_nstep = get_nstep()

  !Flag to decide to write or not to write data
  yaml%flag_print = .false. !(**remove these if generating data for a dependent subroutines**)

  if(yaml%col_print == ii .and. y_nstep==yaml%nstep_print .and. kk == yaml%lev_print) then ! if this column exists in y_lchnk

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
        call open_files('mam_amicphys_1subarea', &  !intent-in
             unit_input, unit_output) !intent-out
        !    unit_input, unit_output, trim(ext_str)) !intent-out, with the use of ext_str


        !start by adding an input string
        call write_input_output_header(unit_input, unit_output,yaml%lchnk_print,yaml%col_print, &
             'mam_amicphys_1subarea',yaml%nstep_print, yaml%lev_print)

        ! add code for writing data here
        
        call write_var(unit_input,unit_output,'do_cond_sub',do_cond_sub)
        call write_var(unit_input,unit_output,'do_rename_sub',do_rename_sub)
        call write_var(unit_input,unit_output,'do_newnuc_sub',do_newnuc_sub)
        call write_var(unit_input,unit_output,'do_coag_sub',do_coag_sub)
        call write_var(unit_input,unit_output,'nstep',nstep)
        call write_var(unit_input,unit_output,'lchnk',lchnk)
        call write_var(unit_input,unit_output,'ii',ii)
        call write_var(unit_input,unit_output,'kk',kk)
        call write_var(unit_input,unit_output,'latndx',latndx)
        call write_var(unit_input,unit_output,'lonndx',lonndx)
        call write_var(unit_input,unit_output,'lund',lund)
        call write_var(unit_input,unit_output,'loffset',loffset)
        call write_var(unit_input,unit_output,'deltat',deltat)
        call write_var(unit_input,unit_output,'jsubarea',jsubarea)
        call write_var(unit_input,unit_output,'nsubarea',nsubarea)
        call write_var(unit_input,unit_output,'iscldy_subarea',iscldy_subarea)
        call write_var(unit_input,unit_output,'afracsub',afracsub)
        call write_var(unit_input,unit_output,'temp',temp)
        call write_var(unit_input,unit_output,'pmid',pmid)
        call write_var(unit_input,unit_output,'pdel',pdel)
        call write_var(unit_input,unit_output,'zmid',zmid)
        call write_var(unit_input,unit_output,'pblh',pblh)
        call write_var(unit_input,unit_output,'relhum',relhum)
        call write_var(unit_input,unit_output,'dgn_a',dgn_a)
        call write_var(unit_input,unit_output,'dgn_awet',dgn_awet)
        call write_var(unit_input,unit_output,'wetdens',wetdens)
        call write_var(unit_input,unit_output,'qgas1',qgas1)
        call write_var(unit_input,unit_output,'qgas3',qgas3)
        call write_var(unit_input,unit_output,'qgas_cur',qgas_cur)
        call write_var(unit_input,unit_output,'qgas_delaa',qgas_delaa)
        call write_var(unit_input,unit_output,'qnum3',qnum3)
        call write_var(unit_input,unit_output,'qnum_cur',qnum_cur)
        call write_var(unit_input,unit_output,'qnum_delaa',qnum_delaa)
        call write_var(unit_input,unit_output,'qaer2',qaer2)
        call write_var(unit_input,unit_output,'qaer3',qaer3)
        call write_var(unit_input,unit_output,'qaer_cur',qaer_cur)
        call write_var(unit_input,unit_output,'qaer_delaa',qaer_delaa)
        call write_var(unit_input,unit_output,'qwtr3',qwtr3)
        call write_var(unit_input,unit_output,'qwtr_cur',qwtr_cur)
        call write_var(unit_input,unit_output,'qnumcw3',qnumcw3)
        call write_var(unit_input,unit_output,'qnumcw_cur',qnumcw_cur)
        call write_var(unit_input,unit_output,'qnumcw_delaa',qnumcw_delaa)
        call write_var(unit_input,unit_output,'qaercw2',qaercw2)
        call write_var(unit_input,unit_output,'qaercw3',qaercw3)
        call write_var(unit_input,unit_output,'qaercw_cur',qaercw_cur)
        call write_var(unit_input,unit_output,'qaercw_delaa',qaercw_delaa)
        call write_var(unit_input,unit_output,'misc_vars_aa_sub',misc_vars_aa_sub%ncluster_tend_nnuc_1grid)

        call write_var(unit_input,unit_output,'n_mode',n_mode)
        call write_var(unit_input,unit_output,'r_universal',r_universal)
        call write_var(unit_input,unit_output,'ntot_poaspec',ntot_poaspec)
        call write_var(unit_input,unit_output,'ntot_soaspec',ntot_soaspec)
        call write_var(unit_input,unit_output,'max_mode',max_mode)
        call write_var(unit_input,unit_output,'max_gas',max_gas)
        call write_var(unit_input,unit_output,'max_aer',max_aer)
        call write_var(unit_input,unit_output,'gaexch_h2so4_uptake_optaa',gaexch_h2so4_uptake_optaa)
        call write_var(unit_input,unit_output,'ngas',ngas)
        call write_var(unit_input,unit_output,'igas_h2so4',igas_h2so4)
        call write_var(unit_input,unit_output,'igas_nh3',igas_nh3)
        call write_var(unit_input,unit_output,'nait',nait)
        call write_var(unit_input,unit_output,'nacc',nacc)
        call write_var(unit_input,unit_output,'max_agepair',max_agepair)
        call write_var(unit_input,unit_output,'n_agepair',n_agepair)
        call write_var(unit_input,unit_output,'nqtendaa',nqtendaa)
        call write_var(unit_input,unit_output,'nqqcwtendaa',nqqcwtendaa)
        call write_var(unit_input,unit_output,'iqtend_cond',iqtend_cond)
        call write_var(unit_input,unit_output,'iqtend_rnam',iqtend_rnam)
        call write_var(unit_input,unit_output,'iqtend_nnuc',iqtend_nnuc)
        call write_var(unit_input,unit_output,'iqtend_coag',iqtend_coag)
        call write_var(unit_input,unit_output,'iqqcwtend_rnam',iqqcwtend_rnam)
        call write_var(unit_input,unit_output,'qwtr3',qwtr3)
        call write_var(unit_input,unit_output,'qwtr3',qwtr3)
  

        !writes aerosol mmr from state%q or q vector (cloud borne and interstitial)
        !"aer_num_only" is .ture. if printing aerosol num only
        !call write_aerosol_mmr_from_stateq(unit_input, unit_output, fld_name,field,aer_num_only)

        !close only the input file, not the output file
        close(unit_input)
        call freeunit(unit_input)
     endif
  endif
#endif
