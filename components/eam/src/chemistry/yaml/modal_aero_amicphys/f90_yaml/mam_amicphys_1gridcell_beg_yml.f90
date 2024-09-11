#ifdef YAML_AMICPHYS_1GRIDCELL
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
        call open_files('mam_amicphys_1gridcell', &  !intent-in
             unit_input, unit_output) !intent-out
        !    unit_input, unit_output, trim(ext_str)) !intent-out, with the use of ext_str


        !start by adding an input string
        call write_input_output_header(unit_input, unit_output,yaml%lchnk_print,yaml%col_print, &
             'mam_amicphys_1gridcell',yaml%nstep_print, yaml%lev_print)

        ! add code for writing data here
        
        call write_var(unit_input,unit_output,'do_cond',do_cond)
        call write_var(unit_input,unit_output,'do_rename',do_rename)
        call write_var(unit_input,unit_output,'do_newnuc',do_newnuc)
        call write_var(unit_input,unit_output,'do_coag',do_coag)
        call write_var(unit_input,unit_output,'nstep',nstep)
        call write_var(unit_input,unit_output,'lchnk',lchnk)
        call write_var(unit_input,unit_output,'i',i)
        call write_var(unit_input,unit_output,'k',k)
        call write_var(unit_input,unit_output,'latndx',latndx)
        call write_var(unit_input,unit_output,'lonndx',lonndx)
        call write_var(unit_input,unit_output,'lund',lund)
        call write_var(unit_input,unit_output,'loffset',loffset)
        call write_var(unit_input,unit_output,'deltat',deltat)
        call write_var(unit_input,unit_output,'nsubarea',nsubarea)
        call write_var(unit_input,unit_output,'ncldy_subarea',ncldy_subarea)
        call write_var(unit_input,unit_output,'iscldy_subarea',iscldy_subarea)
        call write_var(unit_input,unit_output,'afracsub',afracsub)
        call write_var(unit_input,unit_output,'temp',temp)
        call write_var(unit_input,unit_output,'pmid',pmid)
        call write_var(unit_input,unit_output,'pdel',pdel)
        call write_var(unit_input,unit_output,'zmid',zmid)
        call write_var(unit_input,unit_output,'pblh',pblh)
        call write_var(unit_input,unit_output,'relhumsub',relhumsub)
        call write_var(unit_input,unit_output,'dgn_a',dgn_a)
        call write_var(unit_input,unit_output,'dgn_awet',dgn_awet)
        call write_var(unit_input,unit_output,'wetdens',wetdens)
        call write_var(unit_input,unit_output,'qsub1',qsub1)
        call write_var(unit_input,unit_output,'qsub2',qsub2)
        call write_var(unit_input,unit_output,'qqcwsub2',qqcwsub2)
        call write_var(unit_input,unit_output,'qsub3',qsub3)
        call write_var(unit_input,unit_output,'qqcwsub3',qqcwsub3)
        call write_var(unit_input,unit_output,'qaerwatsub3',qaerwatsub3)
        call write_var(unit_input,unit_output,'qsub4',qsub4)
        call write_var(unit_input,unit_output,'qqcwsub4',qqcwsub4)
        call write_var(unit_input,unit_output,'qaerwatsub4',qaerwatsub4)
        call write_var(unit_input,unit_output,'qsub_tendaa',qsub_tendaa)
        call write_var(unit_input,unit_output,'qqcwsub_tendaa',qqcwsub_tendaa)
        call write_var(unit_input,unit_output,'misc_vars_aa',misc_vars_aa%ncluster_tend_nnuc_1grid)


        !writes aerosol mmr from state%q or q vector (cloud borne and interstitial)
        !"aer_num_only" is .ture. if printing aerosol num only
        !call write_aerosol_mmr_from_stateq(unit_input, unit_output, fld_name,field,aer_num_only)

        !close only the input file, not the output file
        close(unit_input)
        call freeunit(unit_input)
     endif
  endif
#endif
