#ifdef YAML_AMICPHYS_1GRIDCELL
  if(yaml%flag_print) then ! if this column exists in lchnk

     !write output header
     call write_output_header(unit_output)

     !start writing data
     
        call write_var(unit_output,'dgn_a',dgn_a)
        call write_var(unit_output,'dgn_awet',dgn_awet)
        call write_var(unit_output,'wetdens',wetdens)
        call write_var(unit_output,'qsub4',qsub4)
        call write_var(unit_output,'qqcwsub4',qqcwsub4)
        call write_var(unit_output,'qaerwatsub3',qaerwatsub3)
        call write_var(unit_output,'qaerwatsub4',qaerwatsub4)
        call write_var(unit_output,'qsub_tendaa',qsub_tendaa)
        call write_var(unit_output,'qqcwsub_tendaa',qqcwsub_tendaa)
        call write_var(unit_output,'misc_vars_aa',misc_vars_aa%ncluster_tend_nnuc_1grid)

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
