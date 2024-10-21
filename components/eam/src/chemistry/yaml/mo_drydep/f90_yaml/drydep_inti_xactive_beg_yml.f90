#ifdef YAML_INTI_DRYDEP

        !open I/O yaml files
        !(with an optional argument to pass a unique string to differentiate file names)
        call open_files('drydep_inti_xactive', &  !intent-in
!             unit_input, unit_output, trim(ext_str)) !intent-out, with the use of ext_str
             unit_input, unit_output) !intent-out

        !start by adding an input string
!!!        call write_input_output_header(unit_input, unit_output,yaml%lchnk_print,yaml%col_print, &
!!!             'drydep_inti_xactive',yaml%nstep_print, yaml%lev_print)

        ! add code for writing data here
        
        call write_var(unit_input,unit_output,'clat',clat(:))
        call write_var(unit_input,unit_output,'wk_lai',wk_lai(:,:,:))
        call write_var(unit_input,unit_output,'lat_lai',lat_lai(:))

        !writes aerosol mmr from state%q or q vector (cloud borne and interstitial)
        !"aer_num_only" is .ture. if printing aerosol num only
        !call write_aerosol_mmr_from_stateq(unit_input, unit_output, fld_name,field,aer_num_only)

        !close only the input file, not the output file
        close(unit_input)
        call freeunit(unit_input)

#endif
