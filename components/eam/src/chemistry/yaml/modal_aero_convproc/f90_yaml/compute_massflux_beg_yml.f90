#ifdef YAML_CONV

  !>

  ! character(len=200) :: ext_str  ! this is used when multiple sets of yaml output is needed to cover different options (e.g., true and false)
  type(yaml_vars), intent(in) :: yaml
  integer  :: unit_input, unit_output, y_nstep
  integer,save :: n_calls=0   ! some subroutines are called multiple times in one timestep, record the number of calls

  !YAML file input generation code- DO NOT PORT to C++
     !<
     !In the case of y_i or y_k are not passed as arguments:
     !if(yaml%col_print >0 .and. y_nstep==yaml%nstep_print) then

     !For generating data for a dependent subroutines where "yaml" derived type is already initialized:
     if(yaml%flag_print) then
     !>

     ! record number of calls that can output yaml file. you only need to write one set of input/output
     n_calls = n_calls+1
     if (n_calls==1) then ! output the first call. change this if writing out other calls
     ! if ((n_calls==1 .and. flag==0) .or. (n_calls==3 .and. flag==1) .or. (n_calls==5 .and. flag==2)) then ! this is an example of writing multiple sets of YAML files


        !open I/O yaml files (it can have an extra optional argument to pass a unique string to differentiate file names)
        call open_files('compute_massflux', &  !intent-in
             unit_input, unit_output) !intent-out
!             unit_input, unit_output, trim(ext_str)) !intent-out, with the use of ext_str


        !start by adding an input string
        call write_input_output_header(unit_input, unit_output,yaml%lchnk_print,yaml%col_print, &
             'compute_tendencies',yaml%nstep_print, yaml%lev_print)

        !< add code for writing data here>
        call write_var(unit_input, unit_output, 'ii', ii)
        call write_var(unit_input, unit_output, 'icol', icol)
        call write_var(unit_input, unit_output, 'ktop', ktop)
        call write_var(unit_input, unit_output, 'kbot', kbot)
        call write_1d_var(unit_input,unit_output,'dpdry_i', pver,dpdry_i(:))
        call write_1d_var(unit_input,unit_output,'du', pver,du(yaml%col_print,:))
        call write_1d_var(unit_input,unit_output,'eu', pver,eu(yaml%col_print,:))
        call write_1d_var(unit_input,unit_output,'ed', pver,ed(yaml%col_print,:))
        call write_var(unit_input, unit_output, 'xx_mfup_max', xx_mfup_max(yaml%col_print))

        !call write_var(unit_input, unit_output, fld_name,field)!write a single variable
        !call write_1d_var(unit_input, unit_output, fld_name,dim,field) ! writes 1D variables of any dimension
        !call write_2d_var(unit_input, unit_output, fld_name, dim1, dim2, field) ! writes 2D variables of any dimension: field(dim1,dim2)

        !writes aerosol mmr from state%q or q vector (cloud borne and interstitial)
        !"aer_num_only" is .ture. if printing aerosol num only
        !call write_aerosol_mmr_from_stateq(unit_input, unit_output, fld_name,field,aer_num_only)

        !close only the input file, not the output file
        close(unit_input)
        call freeunit(unit_input)
     endif
  endif
#endif
