#ifdef YAML_AERO_DRYDEP
  if(yaml%flag_print) then ! if this column exists in lchnk

     !write output header
     call write_output_header(unit_output)

     !start writing data
     
        call write_var(unit_output,'vlc_dry',vlc_dry(yaml%col_print,:))
        call write_var(unit_output,'vlc_trb',vlc_trb(yaml%col_print))
        call write_var(unit_output,'vlc_grv',vlc_grv(yaml%col_print,:))

     !writes aerosol mmr from state%q or q vector(cloud borne and interstitial) in the output python module
     !"aer_num_only" is .ture. if printing aerosol num only
     !call write_output_aerosol_mmr_from_stateq(unit_output, fld_name, field, aer_num_only, inp_out_str)

     !close the output file
     close(unit_output)
     call freeunit(unit_output)
  endif


#endif
