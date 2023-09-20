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
        call open_files('modal_aero_lw', &  !intent-in
             unit_input, unit_output) !intent-out
        !    unit_input, unit_output, trim(ext_str)) !intent-out, with the use of ext_str


        !start by adding an input string
        call write_input_output_header(unit_input, unit_output,yaml%lchnk_print,yaml%col_print, &
             'modal_aero_lw',yaml%nstep_print, yaml%lev_print)

        ! add code for writing data here
        
        call write_var(unit_input,unit_output,'dt',dt)
        call write_var(unit_input,unit_output,'ncol',ncol)
        call write_var(unit_input,unit_output,'state_q',state_q(yaml%col_print,:,:))
        call write_var(unit_input,unit_output,'temperature',temperature(yaml%col_print,:))
        call write_var(unit_input,unit_output,'pmid',pmid(yaml%col_print,:))
        call write_var(unit_input,unit_output,'pdel',pdel(yaml%col_print,:))
        call write_var(unit_input,unit_output,'pdeldry',pdeldry(yaml%col_print,:))
        call write_var(unit_input,unit_output,'cldn',cldn(yaml%col_print,:))

        ! qqcw is an array of pointers that hold references
        ! to 2D arrays (column,vertical level) of number / species mass.
        ! Convert qqcw for a given column to a 2D array (number/mass index,
        ! vertical level)
        ! and then write that 2D array to the yaml / python files.
        if (allocated(qqcw_in)) deallocate(qqcw_in)
        allocate(qqcw_in(size(qqcw),pver))
        do imm=1,size(qqcw)
           if (imm<=15) then
                qqcw_in(imm,:) = -9999.9
           else
                qqcw_in(imm,:) = qqcw(imm)%fld(yaml%col_print,:)
           endif
        enddo
        call write_var(unit_input,unit_output,'qqcw',qqcw_in(:,:))
        deallocate(qqcw_in)

        ! other module variables
        call write_var(unit_input,unit_output,'top_lev',top_lev)
        call write_var(unit_input,unit_output,'nlwbands',nlwbands)
        call write_var(unit_input,unit_output,'absplw',reshape(absplw,(/size(absplw)/)))
        call write_var(unit_input,unit_output,'shape_of_absplw', shape(absplw))
        call write_var(unit_input,unit_output,'specrefndxlw_real',real(specrefndxlw))
        call write_var(unit_input,unit_output,'specrefndxlw_imag',aimag(specrefndxlw))
        call write_var(unit_input,unit_output,'refrtablw',reshape(refrtablw,(/size(refrtablw)/)))
        call write_var(unit_input,unit_output,'refitablw',reshape(refitablw,(/size(refitablw)/)))


        !writes aerosol mmr from state%q or q vector (cloud borne and interstitial)
        !"aer_num_only" is .ture. if printing aerosol num only
        !call write_aerosol_mmr_from_stateq(unit_input, unit_output, fld_name,field,aer_num_only)

        !close only the input file, not the output file
        close(unit_input)
        call freeunit(unit_input)
     endif
  endif
#endif
