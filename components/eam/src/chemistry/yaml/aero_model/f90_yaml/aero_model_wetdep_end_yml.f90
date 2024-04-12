#ifdef YAML_AERO_MODEL_WETDEP
  if(yaml%flag_print) then ! if this column exists in lchnk

     !write output header
     call write_output_header(unit_output)

     !start writing data

        call write_var(unit_output,'dgncur_a',dgncur_a(yaml%col_print,:,:))
        call write_var(unit_output,'wetdens',wetdens(yaml%col_print,:,:))
        call write_var(unit_output,'qaerwat',qaerwat(yaml%col_print,:,:))
        call write_var(unit_output,'dgnumwet',dgnumwet(yaml%col_print,:,:))
        call write_var(unit_output,'fracis',fracis(yaml%col_print,:,:))  

        call write_var(unit_output,'aerdepwetis',aerdepwetis(yaml%col_print,:))
        call write_var(unit_output,'aerdepwetcw',aerdepwetcw(yaml%col_print,:))


        ! write out fields that are set within intent-out structure 'cam_out'

        call write_var(unit_output,'cam_out_bcphiwet',cam_out%bcphiwet(yaml%col_print))
        call write_var(unit_output,'cam_out_bcphidry',cam_out%bcphidry(yaml%col_print))
        call write_var(unit_output,'cam_out_ocphiwet',cam_out%ocphiwet(yaml%col_print))
        call write_var(unit_output,'cam_out_ocphidry',cam_out%ocphidry(yaml%col_print))
        call write_var(unit_output,'cam_out_dstwet1',cam_out%dstwet1(yaml%col_print))
        call write_var(unit_output,'cam_out_dstwet2',cam_out%dstwet2(yaml%col_print))
        call write_var(unit_output,'cam_out_dstwet3',cam_out%dstwet3(yaml%col_print))
        call write_var(unit_output,'cam_out_dstwet4',cam_out%dstwet4(yaml%col_print))

        ! ptend is passed as intent-out to aero_model_wetdep, so only include as output data.
        ! The following extract the fields from ptend used as output and write them to the data files.
        ! Use -9999.9 to fill in unused portions of the constituent arrays

        if (allocated(ptend_lq_1darr)) deallocate(ptend_lq_1darr)
        allocate(ptend_lq_1darr(size(ptend%lq)))

        if (allocated(ptend_q_3darr)) deallocate(ptend_q_3darr)
        allocate(ptend_q_3darr(size(ptend%q,1),size(ptend%q,2),size(ptend%q,3)))

        do imm=1,size(ptend%lq)
           if (imm<=15) then
                ptend_lq_1darr(imm) = -9999.9
                ptend_q_3darr(:,:,imm) = -9999.9
           else
                ptend_lq_1darr(imm) = ptend%lq(imm)
                ptend_q_3darr(:,:,imm) = ptend%q(:,:,imm)
           endif
        enddo
        call write_var(unit_output,'ptend_lq',ptend_lq_1darr(:))
        call write_var(unit_output,'ptend_lq',ptend_q_3darr(yaml%col_print,:,:))
        deallocate(ptend_lq_1darr)
        deallocate(ptend_q_3darr)

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
