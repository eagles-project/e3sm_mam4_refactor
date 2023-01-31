#ifdef YAML_WETDEP
  ! YAML file output generation code- DO NOT PORT to C++
  if(yaml%flag_print) then ! if this column exists in lchnk

     !write output header
     call write_output_header(unit_output)

     !start writing data
     !<add code for writing data here>
     call write_output_var(unit_output, 'precabx_base_new', precabx_base_new)
     call write_output_var(unit_output, 'precabx_new', precabx_new)
     call write_output_var(unit_output, 'scavabx_new', scavabx_new)
     call write_output_var(unit_output, 'resusp_x', resusp_x)

     !call write_output_var(unit_output, fld_name, field, inp_out_str)  !write a single output variable
     !call write_1d_output_var(unit_output, fld_name, dim, field, inp_out_str) !writes 1D variables of any dimension in the output python module
     !call write_2d_output_var(unit_output, fld_name, dim1, dim2, field, inp_out_str) !writes 2D variables of any dimension in the output python module

     !writes aerosol mmr from state%q or q vector(cloud borne and interstitial) in the output python module
     !"aer_num_only" is .ture. if printing aerosol num only
     !call write_output_aerosol_mmr_from_stateq(unit_output, fld_name, field, aer_num_only, inp_out_str)

     !close the output file
     close(unit_output)
     call freeunit(unit_output)
  endif

!! find lag/lon/y_nstep/y_k
!if(masterproc .and. y_nstep>353 .and. y_nstep<360 .and. (resusp_x/=0.0 .or. scavabx_new/=0.0)) then
!       write(102,*)'phys_debug_lat = ',get_lat(y_lchnk, y_i), &
!      ' phys_debug_lon = ', get_lon(y_lchnk, y_i), 'kk=',y_k, y_nstep, &
!        mam_prevap_resusp_optcc,scavabx_new,resusp_x
!endif
!! find n_calls. the first condition should be consistent with *_beg_yml.f90 file
!if((yaml%col_print == y_i .and. y_nstep==yaml%nstep_print .and. y_k == yaml%lev_print) .and. (resusp_x/=0.0 .or. scavabx_new/=0.0)) then
!       write(102,*) 'n_calls = ', n_calls,mam_prevap_resusp_optcc,scavabx_new,resusp_x
!endif

#endif
