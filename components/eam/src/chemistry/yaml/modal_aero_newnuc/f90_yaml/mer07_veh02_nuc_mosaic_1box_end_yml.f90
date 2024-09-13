#ifdef YAML_VEHPBLNUC
  if(yaml%flag_print) then ! if this column exists in lchnk

     !write output header
     call write_output_header(unit_output)

     !start writing data
     
        call write_var(unit_output,'isize_nuc',isize_nuc)
        call write_var(unit_output,'qnuma_del',qnuma_del)
        call write_var(unit_output,'qso4a_del',qso4a_del)
        call write_var(unit_output,'qnh4a_del',qnh4a_del)
        call write_var(unit_output,'qh2so4_del',qh2so4_del)
        call write_var(unit_output,'dens_nh4so4a',dens_nh4so4a)

!  below are optional output arguments to mer07_veh02_nuc_mosaic_1box
        if(present(dnclusterdt)) then
           call write_var(unit_output,'dnclusterdt',dnclusterdt)
        endif

     !writes aerosol mmr from state%q or q vector(cloud borne and interstitial) in the output python module
     !"aer_num_only" is .ture. if printing aerosol num only
     !call write_output_aerosol_mmr_from_stateq(unit_output, fld_name, field, aer_num_only, inp_out_str)

     !close the output file
     close(unit_output)
     call freeunit(unit_output)
  endif

!------------------------------------------------------------------------------------------------------------------------------
! This part of the code is for finding the best lat/lon and n_calls.

! find lag/lon/y_nstep/y_k
!if(masterproc .and. some condition) then
!  write(102,*)'phys_debug_lat = ',get_lat(y_lchnk, yaml%col_print), &
!  ' phys_debug_lon = ', get_lon(y_lchnk, yaml%col_print), get_nstep(), yaml%lev_print
!endif

! find n_calls. the if statement below should be consistent with the if condition in the *_beg_yml.f90 file
! so that the corresponding n_calls are printed out here.
! Note that n_calls may change when lat/lon/y_nstep/y_k change
!if(yaml%col_print == y_i .and. y_nstep==yaml%nstep_print .and. y_k == yaml%lev_print) then
!      write(103,*) 'n_calls = ',n_calls, <any of the input/output variables to choose n_calls>
!endif
!------------------------------------------------------------------------------------------------------------------------------

#endif
