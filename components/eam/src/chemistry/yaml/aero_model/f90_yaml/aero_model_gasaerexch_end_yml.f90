#ifdef YAML_AERO_MODEL_GASAEREXCH
  if(yaml%flag_print) then ! if this column exists in lchnk

     !write output header
     call write_output_header(unit_output)

     !start writing data
     
        if (allocated(qqcw_out)) deallocate(qqcw_out)
        allocate(qqcw_out(size(qqcw),pver))
        do imm=1,size(qqcw)
           if (imm<=15) then
                qqcw_out(imm,:) = -9999.9
           else
                qqcw_out(imm,:) = qqcw(imm)%fld(yaml%col_print,:)
           endif
        enddo
        call write_var(unit_output,'qqcw',qqcw_out)
        deallocate(qqcw_out)

        call write_var(unit_output,'vmr',vmr(yaml%col_print,:,:))

        ! diagnostics fields (outfld)
        call write_var(unit_output,'AQ_'//trim(solsym(gas_pcnst))//'_dvmrdt_col',dvmrdt_col(yaml%col_print))

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
