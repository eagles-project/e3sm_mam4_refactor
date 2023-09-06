#ifdef YAML_LINSTRAT
  if(yaml%flag_print) then ! if this column exists in lchnk

     !write output header
     call write_output_header(unit_output)

     !start writing data

     if(multicol) then
     
        call write_var(unit_output,'o3_vmr',o3_vmr(:,yaml%lev_print))

     !write data sent to outfld

        call write_var(unit_output,'do3_linoz',do3_linoz(:,yaml%lev_print))
        call write_var(unit_output,'do3_linoz_psc',do3_linoz_psc(:,yaml%lev_print))
        call write_var(unit_output,'ss_o3',ss_o3(:,yaml%lev_print))
        call write_var(unit_output,'o3col_du_diag',o3col_du_diag(:,yaml%lev_print))
        call write_var(unit_output,'o3clim_linoz_diag',o3clim_linoz_diag(:,yaml%lev_print))

     else

        call write_var(unit_output,'o3_vmr',o3_vmr(yaml%col_print,yaml%lev_print))

     !write data sent to outfld

        call write_var(unit_output,'do3_linoz',do3_linoz(yaml%col_print,yaml%lev_print))
        call write_var(unit_output,'do3_linoz_psc',do3_linoz_psc(yaml%col_print,yaml%lev_print))
        call write_var(unit_output,'ss_o3',ss_o3(yaml%col_print,yaml%lev_print))
        call write_var(unit_output,'o3col_du_diag',o3col_du_diag(yaml%col_print,yaml%lev_print))
        call write_var(unit_output,'o3clim_linoz_diag',o3clim_linoz_diag(yaml%col_print,yaml%lev_print))

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
