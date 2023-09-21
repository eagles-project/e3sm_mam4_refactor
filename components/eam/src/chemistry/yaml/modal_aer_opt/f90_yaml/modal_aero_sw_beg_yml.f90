#ifdef YAML_OPT
  !-----------------------------------------------------------------------------------------
  !"lchnk" is needed for the following code to work,
  ! temporarily pass it along from upper level subroutines
  ! as y_lchnk and uncomment the following code:
  integer :: y_lchnk
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
  integer  :: imm
  real(r8), allocatable :: qqcw_in(:,:)
  real(r8), allocatable :: qqcw_out(:,:)

  ! some subroutines are called multiple times in one timestep, record the number of calls
  integer,save :: n_calls=0

  y_lchnk = lchnk

  !populate YAML structure
  !(**remove yaml%lev_print, nstep_print, col_print if generating data for a dependent subroutines**)
  yaml%lev_print = 0       !level
  yaml%nstep_print = 355 !time step

  yaml%col_print = icolprnt(y_lchnk)                !column to write data

  !current_time step
  y_nstep = get_nstep()

  !Flag to decide to write or not to write data
  yaml%flag_print = .false. !(**remove these if generating data for a dependent subroutines**)

!  if(yaml%col_print == y_i .and. y_nstep==yaml%nstep_print .and. y_k == yaml%lev_print) then ! if this column exists in y_lchnk

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
        yaml%lchnk_print = y_lchnk
        yaml%flag_print  = .true.


        !open I/O yaml files
        !(with an optional argument to pass a unique string to differentiate file names)
        call open_files('modal_aero_sw', &  !intent-in
             unit_input, unit_output) !intent-out
        !    unit_input, unit_output, trim(ext_str)) !intent-out, with the use of ext_str

        !start by adding an input string
        call write_input_output_header(unit_input, unit_output,yaml%lchnk_print,yaml%col_print, &
             'modal_aero_sw',yaml%nstep_print, yaml%lev_print)

        ! add code for writing data here
        
        call write_var(unit_input,unit_output,'dt',dt)
        call write_var(unit_input,unit_output,'ncol',ncol)
        call write_var(unit_input,unit_output,'state_q',state_q(yaml%col_print,:,:))
        call write_var(unit_input,unit_output,'state_zm',state_zm(yaml%col_print,:))
        call write_var(unit_input,unit_output,'temperature',temperature(yaml%col_print,:))
        call write_var(unit_input,unit_output,'pmid',pmid(yaml%col_print,:))
        call write_var(unit_input,unit_output,'pdel',pdel(yaml%col_print,:))
        call write_var(unit_input,unit_output,'pdeldry',pdeldry(yaml%col_print,:))
        call write_var(unit_input,unit_output,'cldn',cldn(yaml%col_print,:))
        call write_var(unit_input,unit_output,'nnite',nnite)
        call write_var(unit_input,unit_output,'idxnite',idxnite)
        call write_var(unit_input,unit_output,'is_cmip6_volc',is_cmip6_volc)
        call write_var(unit_input,unit_output,'ext_cmip6_sw',ext_cmip6_sw(yaml%col_print,:))
        call write_var(unit_input,unit_output,'trop_level',trop_level(yaml%col_print))
!do imm=1,size(qqcw)
!write(102,*) shape(qqcw(imm)%fld)
!enddo
        if (allocated(qqcw_in)) deallocate(qqcw_in)
        allocate(qqcw_in(size(qqcw),pver))
        do imm=1,size(qqcw)
           if (imm<=15) then
                qqcw_in(imm,:) = -9999.9
           else
                qqcw_in(imm,:) = qqcw(imm)%fld(yaml%col_print,:)
           endif
        enddo
        call write_var(unit_input,unit_output,'qqcw',qqcw_in)
        deallocate(qqcw_in)

        ! module variables
        call write_var(unit_input,unit_output,'top_lev',top_lev)
        call write_var(unit_input,unit_output,'nswbands',nswbands)
        call write_var(unit_input,unit_output,'extpsw',reshape(extpsw,(/size(extpsw)/)))
        call write_var(unit_input,unit_output,'shape_of_extpsw', shape(extpsw))
        call write_var(unit_input,unit_output,'abspsw',reshape(abspsw,(/size(abspsw)/)))
        call write_var(unit_input,unit_output,'shape_of_abspsw', shape(abspsw))
        call write_var(unit_input,unit_output,'asmpsw',reshape(asmpsw,(/size(asmpsw)/)))
        call write_var(unit_input,unit_output,'shape_of_asmpsw', shape(asmpsw))
        call write_var(unit_input,unit_output,'specrefndxsw_real',real(specrefndxsw))
        call write_var(unit_input,unit_output,'specrefndxsw_imag',aimag(specrefndxsw))
        call write_var(unit_input,unit_output,'refrtabsw',reshape(refrtabsw,(/size(refrtabsw)/)))
        call write_var(unit_input,unit_output,'refitabsw',reshape(refitabsw,(/size(refitabsw)/)))
        call write_var(unit_input,unit_output,'crefwsw_real',real(crefwsw))
        call write_var(unit_input,unit_output,'crefwsw_imag',aimag(crefwsw))


        !writes aerosol mmr from state%q or q vector (cloud borne and interstitial)
        !"aer_num_only" is .ture. if printing aerosol num only
        !call write_aerosol_mmr_from_stateq(unit_input, unit_output, fld_name,field,aer_num_only)

        !close only the input file, not the output file
        close(unit_input)
        call freeunit(unit_input)
     endif
  endif
#endif
