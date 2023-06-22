#ifdef YAML_CONV

  !>

  ! character(len=200) :: ext_str  ! this is used when multiple sets of yaml output is needed to cover different options (e.g., true and false)
  type(yaml_vars), intent(in) :: yaml
  integer  :: unit_input, unit_output, y_nstep
  integer,save :: n_calls=0   ! some subroutines are called multiple times in one timestep, record the number of calls


!  if(yaml%col_print == y_i .and. y_nstep==yaml%nstep_print .and. y_k == yaml%lev_print) then ! if this column exists in y_lchnk
     !<
     !In the case of y_i or y_k are not passed as arguments:
     !if(yaml%col_print >0 .and. y_nstep==yaml%nstep_print) then

     !For generating data for a dependent subroutines where "yaml" derived type is already initialized:
     if(yaml%flag_print .and. convtype=='deep') then
     !>

     ! record number of calls that can output yaml file. you only need to write one set of input/output
     n_calls = n_calls+1
     if (n_calls==1) then ! output the first call. change this if writing out other calls
     ! if ((n_calls==1 .and. flag==0) .or. (n_calls==3 .and. flag==1) .or. (n_calls==5 .and. flag==2)) then ! this is an example of writing multiple sets of YAML files

        !open I/O yaml files (it can have an extra optional argument to pass a unique string to differentiate file names)
        call open_files('ma_convproc_tend', &  !intent-in
             unit_input, unit_output) !intent-out
!             unit_input, unit_output, trim(ext_str)) !intent-out, with the use of ext_str


        !start by adding an input string
        call write_input_output_header(unit_input, unit_output,yaml%lchnk_print,yaml%col_print, &
             'compute_tendencies',yaml%nstep_print, yaml%lev_print)

        !< add code for writing data here>
        call write_var(unit_input,unit_output,'ncnst',ncnst)
        call write_var(unit_input,unit_output,'dt',dt)
        call write_1d_var(unit_input,unit_output,'temperature',pver,temperature(yaml%col_print,:))
        call write_1d_var(unit_input,unit_output,'pmid', pver,pmid(yaml%col_print,:))
        call write_2d_var(unit_input, unit_output,'qnew', pver,ncnst, qnew(yaml%col_print, :,:))
        call write_1d_var(unit_input,unit_output,'mu', pver,mu(1,:))
        call write_1d_var(unit_input,unit_output,'md', pver,md(1,:))
        call write_1d_var(unit_input,unit_output,'du', pver,du(1,:))
        call write_1d_var(unit_input,unit_output,'eu', pver,eu(1,:))
        call write_1d_var(unit_input,unit_output,'ed', pver,ed(1,:))
        call write_1d_var(unit_input,unit_output,'dp', pver,dp(1,:))
        call write_1d_var(unit_input,unit_output,'dpdry', pver,dpdry(1,:))
        call write_var(unit_input,unit_output,'jt', jt(1))
        call write_var(unit_input,unit_output,'mx', mx(1))
        call write_var(unit_input,unit_output,'ideep', ideep(1))
        call write_var(unit_input,unit_output,'il1g',il1g)
        call write_var(unit_input,unit_output,'il2g',il2g)
        call write_1d_var(unit_input,unit_output,'cldfrac', pver,cldfrac(yaml%col_print,:))
        call write_1d_var(unit_input,unit_output,'icwmr', pver,icwmr(yaml%col_print,:))
        call write_1d_var(unit_input,unit_output,'rprd', pver,rprd(yaml%col_print,:))
        call write_1d_var(unit_input,unit_output,'evapc',pver,evapc(yaml%col_print,:))
        call write_1d_var(unit_input,unit_output,'doconvproc', ncnst,doconvproc(:))
        call write_1d_var(unit_input, unit_output,'species_class', pcnst, species_class(:))
        call write_var(unit_input,unit_output,'nsrflx',nsrflx)
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
