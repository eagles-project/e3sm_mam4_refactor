#ifdef YAML_AERO_DRYDEP
  if(yaml%flag_print) then ! if this column exists in lchnk

     !write output header
     call write_output_header(unit_output)

     !start writing data
     
        call write_var(unit_output,'vlc_grv',vlc_grv(yaml%col_print,:))
        
     !close the output file
     close(unit_output)
     call freeunit(unit_output)
  endif
#endif
