#ifdef YAML_WETDEP
  ! YAML file output generation code- DO NOT PORT to C++
  if(yaml%flag_print) then ! if this column exists in lchnk
  
     !write output header
     call write_output_header(unit_output)

     !start writing data
     !<add code for writing data here>
     call write_1d_output_var(unit_output, 'cldv',pver, cldv(yaml%col_print,:))
     call write_1d_output_var(unit_output, 'sumppr_all',pver, sumppr_all(yaml%col_print,:))

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
#endif
