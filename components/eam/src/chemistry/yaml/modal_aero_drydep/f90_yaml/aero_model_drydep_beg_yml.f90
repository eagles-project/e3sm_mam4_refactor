#ifdef YAML_AERO_DRYDEP
  
  integer :: y_lchnk

  type(yaml_vars) :: yaml
  integer  :: unit_input, unit_output, y_nstep

  character(len=200)::vname !name for qqcw constituents
  y_lchnk = lchnk

  !populate YAML structure
  yaml%lev_print = 52      !level

  yaml%col_print = icolprnt(y_lchnk)                !column to write data

  !current_time step
  y_nstep = get_nstep()

  !Flag to decide to write or not to write data
  yaml%flag_print = .false.

  !-----------------------------------------------------------------------------------------
  !In the case of y_i or y_k are not passed as arguments, use the following if condition:
  if(yaml%col_print >0 .and. y_nstep>1380) then
     !-----------------------------------------------------------------------------------------
     yaml%nstep_print = y_nstep
     yaml%lchnk_print = y_lchnk
     yaml%flag_print  = .true.


     !open I/O yaml files
     !(with an optional argument to pass a unique string to differentiate file names)
     call open_files('aero_model_drydep', &  !intent-in
          unit_input, unit_output) !intent-out

     !start by adding an input string
     call write_input_output_header(unit_input, unit_output,yaml%lchnk_print,yaml%col_print, &
          'aero_model_drydep',yaml%nstep_print, yaml%lev_print)

     ! add code for writing data here
     call write_var(unit_input,unit_output,'ntot_amode',ntot_amode)
     call write_var(unit_input,unit_output,'lchnk',lchnk)
     call write_var(unit_input,unit_output,'ncol',ncol)
     call write_var(unit_input,unit_output,'psetcols',psetcols)
     call write_var(unit_input,unit_output,'tair',tair(yaml%col_print,:))
     call write_var(unit_input,unit_output,'pmid',pmid(yaml%col_print,:))
     call write_var(unit_input,unit_output,'pint',pint(yaml%col_print,:))
     call write_var(unit_input,unit_output,'pdel',pdel(yaml%col_print,:))
     do icnst = 16,pcnst
        write(vname,'(A,I2)')'statq_',icnst
        call write_var(unit_input,unit_output,vname,state_q(yaml%col_print,:,icnst))
     enddo
     call write_var(unit_input,unit_output,'dgncur_awet',dgncur_awet(yaml%col_print,:,:))
     call write_var(unit_input,unit_output,'wetdens',wetdens(yaml%col_print,:,:))
     do icnst = 16,pcnst
        write(vname,'(A,I2)')'qqcw_',icnst
        call write_var(unit_input,unit_output,trim(vname),qqcw(icnst)%fld(yaml%col_print,:))
     enddo
     call write_var(unit_input,unit_output,'obklen',obklen(yaml%col_print))
     call write_var(unit_input,unit_output,'ustar',ustar(yaml%col_print))
     call write_var(unit_input,unit_output,'landfrac',landfrac(yaml%col_print))
     call write_var(unit_input,unit_output,'icefrac',icefrac(yaml%col_print))
     call write_var(unit_input,unit_output,'ocnfrac',ocnfrac(yaml%col_print))
     call write_var(unit_input,unit_output,'fricvelin',fricvelin(yaml%col_print))
     call write_var(unit_input,unit_output,'ram1in',ram1in(yaml%col_print))
     call write_var(unit_input,unit_output,'dt',dt)

     !close only the input file, not the output file
     close(unit_input)
     call freeunit(unit_input)
  endif
#endif
