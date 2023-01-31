#ifdef YAML_WETDEP
  ! YAML file output generation code- DO NOT PORT to C++
  if(yaml%flag_print) then ! if this column exists in lchnk

     !write output header
     call write_output_header(unit_output)

     !start writing data
     !<add code for writing data here>
     call write_1d_output_var(unit_output, 'fracis', pver,fracis(yaml%col_print,:))
     call write_1d_output_var(unit_output, 'scavt',pver,scavt(yaml%col_print,:))
     call write_1d_output_var(unit_output, 'iscavt',pver,iscavt(yaml%col_print,:))
     call write_1d_output_var(unit_output, 'icscavt', pver,icscavt(yaml%col_print,:))
     call write_1d_output_var(unit_output, 'isscavt', pver,isscavt(yaml%col_print,:))
     call write_1d_output_var(unit_output, 'bcscavt', pver,bcscavt(yaml%col_print,:))
     call write_1d_output_var(unit_output, 'bsscavt', pver,bsscavt(yaml%col_print,:))
     call write_1d_output_var(unit_output, 'rcscavt', pver,rcscavt(yaml%col_print,:))
     call write_1d_output_var(unit_output, 'rsscavt', pver,rsscavt(yaml%col_print,:))

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
!if(masterproc .and. y_nstep==355 .and. any(scavt(yaml%col_print,:)/=0.0) .and. any(bcscavt(yaml%col_print,:)/=0.0) .and. any(bsscavt(yaml%col_print,:)/=0.0) .and. any(rcscavt(yaml%col_print,:)/=0.0) .and. any(rsscavt(yaml%col_print,:)/=0.0)) then
!       write(102,*) 'phys_debug_lat = ',get_lat(y_lchnk, yaml%col_print), &
!      ' phys_debug_lon = ', get_lon(y_lchnk, yaml%col_print),  mam_prevap_resusp_optcc
!endif
!! find n_calls. the first condition should be consistent with *_beg_yml.f90 file
!if((yaml%col_print >0 .and. y_nstep==yaml%nstep_print)) then
!       write(102,*) 'n_calls = ',n_calls,mam_prevap_resusp_optcc,maxval(abs(scavt(yaml%col_print,:))),maxval(abs(iscavt(yaml%col_print,:))),maxval(abs(isscavt(yaml%col_print,:))),maxval(abs(icscavt(yaml%col_print,:))),maxval(abs(bcscavt(yaml%col_print,:))),maxval(abs(bsscavt(yaml%col_print,:))),maxval(abs(rcscavt(yaml%col_print,:))),maxval(abs(rsscavt(yaml%col_print,:)))
!endif
#endif
