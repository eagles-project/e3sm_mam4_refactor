#ifdef YAML_SOXMOD
  if(yaml%flag_print) then ! if this column exists in lchnk

     !write output header
     call write_output_header(unit_output)

     !start writing data
     
        call write_var(unit_output,'conc_obj%xlwc',conc_obj%xlwc(yaml%col_print,yaml%lev_print))
        call write_var(unit_output,'conc_obj%so4c',conc_obj%so4c(yaml%col_print,yaml%lev_print))

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
!if(yaml%flag_print) then
!    do kk = 1,73
!      write(103,*) 'kk = ',kk, cldfrc(yaml%col_print,kk), lwc(yaml%col_print,kk), conc_obj%so4c(yaml%col_print,kk)
!    enddo
!endif
!------------------------------------------------------------------------------------------------------------------------------

#endif
