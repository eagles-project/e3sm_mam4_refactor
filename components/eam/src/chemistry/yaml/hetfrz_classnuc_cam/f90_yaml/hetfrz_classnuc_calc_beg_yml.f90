#ifdef YAML_HETFRZ_CLASSNUC

  !Following code was used to get lat-lon and time step where
  !there are max non-zero values in factnum array.
  !do icol = 1, ncol
  !   if (any(factnum(icol,:,:)>0)) then
  !      write(102,*)''
  !      write(102,*)'BALLI:',get_lat(lchnk,icol),get_lon(lchnk,icol),get_nstep()
  !      write(102,*)count(factnum(icol,:,:)>0.0)
  !   endif
  !enddo


  character(len=200) :: ext_str  ! this is used when multiple sets of yaml output is needed to cover different options (e.g., true and false)
  type(yaml_vars) :: yaml
  integer  :: unit_input, unit_output, y_nstep, icolp, y_lchnk
  integer,save :: n_calls=0   ! some subroutines are called multiple times in one timestep, record the number of calls

  y_lchnk = lchnk

  !populate YAML structure
  yaml%lev_print = -1 !level (**remove these if generating data for a dependent subroutines**)

  ! set ext_str if there are multiple sets of yaml output to be write out
  ! here gives an example that "flag" in the code can be 0, 1, or 2:
  ! write(ext_str,'(I2)') flag
  ! ext_str = 'flag_'//adjustl(ext_str)


  !YAML file input generation code- DO NOT PORT to C++
  !print all inputs one-by-one at column "yaml%col_print"
  yaml%col_print = icolprnt(y_lchnk) !column to write data(**remove these if generating data for a dependent subroutines**)
  y_nstep = get_nstep() !time step (**remove these if generating data for a dependent subroutines**)

  yaml%flag_print = .false. ! to write or not to write data (**remove these if generating data for a dependent subroutines**)
  if(yaml%col_print >0 .and. y_nstep>1400) then ! if this column exists in y_lchnk
     icolp = yaml%col_print
     yaml%nstep_print = y_nstep !time step(**remove these if generating data for a dependent subroutines**)
     !<
     !In the case of y_i or y_k are not passed as arguments:
     !if(yaml%col_print >0 .and. y_nstep==yaml%nstep_print) then

     !For generating data for a dependent subroutines where "yaml" derived type is already initialized:
     !if(yaml%flag_print) then
     !>

     ! record number of calls that can output yaml file. you only need to write one set of input/output
     !n_calls = n_calls+1
     !if (n_calls==1) then ! output the first call. change this if writing out other calls
     ! if ((n_calls==1 .and. flag==0) .or. (n_calls==3 .and. flag==1) .or. (n_calls==5 .and. flag==2)) then ! this is an example of writing multiple sets of YAML files

        yaml%lchnk_print = y_lchnk !(**remove these if generating data for a dependent subroutines**)
        yaml%flag_print  = .true.!(**remove these if generating data for a dependent subroutines**)

        !open I/O yaml files (it can have an extra optional argument to pass a unique string to differentiate file names)
        call open_files('hetfrz_classnuc_calc', &  !intent-in
             unit_input, unit_output) !intent-out
!             unit_input, unit_output, trim(ext_str)) !intent-out, with the use of ext_str


        !start by adding an input string
        call write_input_output_header(unit_input, unit_output,yaml%lchnk_print,yaml%col_print, &
             'compute_tendencies',yaml%nstep_print, yaml%lev_print)

        !< add code for writing data here>
        call write_var(unit_input, unit_output, 'ncnst', ncnst)
        call write_var(unit_input, unit_output, 'pi', pi)
        call write_var(unit_input, unit_output, 'rhoh2o',rhoh2o)
        call write_var(unit_input, unit_output, 'deltatin',deltatin)
        call write_var(unit_input, unit_output, 'rair',rair)
        call write_var(unit_input, unit_output, 'mincld',mincld)
        call write_var(unit_input, unit_output, 'temperature',pver,temperature(icolp,:))
        call write_var(unit_input, unit_output, 'pmid',       pver,pmid(icolp,:))
        call write_var(unit_input, unit_output, 'rho',        pver,rho(icolp,:))
        call write_var(unit_input, unit_output, 'ast',        pver,ast(icolp,:))
        call write_var(unit_input, unit_output, 'qc',         pver,qc(icolp,:))
        call write_var(unit_input, unit_output, 'nc',         pver,nc(icolp,:))

        call write_var(unit_input, unit_output, 'state_q_bc_accum', pver, state_q(icolp, :pver, lptr_bc_a_amode(modeptr_accum)))
        call write_var(unit_input, unit_output, 'state_q_pom_accum', pver, state_q(icolp, :pver, lptr_pom_a_amode(modeptr_accum)))
        call write_var(unit_input, unit_output, 'state_q_soa_accum', pver, state_q(icolp, :pver, lptr_soa_a_amode(modeptr_accum)))
        call write_var(unit_input, unit_output, 'state_q_dust_accum', pver, state_q(icolp, :pver, lptr_dust_a_amode(modeptr_accum)))
        call write_var(unit_input, unit_output, 'state_q_nacl_accum', pver, state_q(icolp, :pver, lptr_nacl_a_amode(modeptr_accum)))
        call write_var(unit_input, unit_output, 'state_q_mom_accum', pver, state_q(icolp, :pver, lptr_mom_a_amode(modeptr_accum)))
        call write_var(unit_input, unit_output, 'state_q_num_accum', pver, state_q(icolp, :pver, numptr_amode(modeptr_accum)))

        call write_var(unit_input, unit_output, 'state_q_dust_coarse', pver, state_q(icolp, :pver, lptr_dust_a_amode(modeptr_coarse)))
        call write_var(unit_input, unit_output, 'state_q_nacl_coarse', pver, state_q(icolp, :pver, lptr_nacl_a_amode(modeptr_coarse)))
        call write_var(unit_input, unit_output, 'state_q_so4_coarse', pver, state_q(icolp, :pver, lptr_so4_a_amode(modeptr_coarse)))
        call write_var(unit_input, unit_output, 'state_q_bc_coarse', pver, state_q(icolp, :pver, lptr_bc_a_amode(modeptr_coarse)))
        call write_var(unit_input, unit_output, 'state_q_pom_coarse', pver, state_q(icolp, :pver, lptr_pom_a_amode(modeptr_coarse)))
        call write_var(unit_input, unit_output, 'state_q_soa_coarse', pver, state_q(icolp, :pver, lptr_soa_a_amode(modeptr_coarse)))
        call write_var(unit_input, unit_output, 'state_q_mom_coarse', pver, state_q(icolp, :pver, lptr_mom_a_amode(modeptr_coarse)))
        call write_var(unit_input, unit_output, 'state_q_num_coarse', pver, state_q(icolp, :pver, numptr_amode(modeptr_coarse)))

        call write_var(unit_input, unit_output, 'state_q_bc_pcarbon', pver, state_q(icolp, :pver, lptr_bc_a_amode(modeptr_pcarbon)))
        call write_var(unit_input, unit_output, 'state_q_pom_pcarbon', pver, state_q(icolp, :pver, lptr_pom_a_amode(modeptr_pcarbon)))
        call write_var(unit_input, unit_output, 'state_q_mom_pcarbon', pver, state_q(icolp, :pver, lptr_mom_a_amode(modeptr_pcarbon)))
        call write_var(unit_input, unit_output, 'state_q_num_pcarbon', pver, state_q(icolp, :pver, numptr_amode(modeptr_pcarbon)))

        call write_var(unit_input, unit_output, 'aer_cb', pver, ncnst, aer_cb(icolp,:,:))

        call write_var(unit_input, unit_output, 'factnum'    ,pver, nmodes, factnum(icolp,:,:))

        close(unit_input)
        call freeunit(unit_input)
     !endif
     endif
#endif
