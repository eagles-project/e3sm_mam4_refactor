#ifdef YAML_OPT
  if(yaml%flag_print) then ! if this column exists in lchnk

     !write output header
     call write_output_header(unit_output)

     !start writing data
        if (allocated(qqcw_out)) deallocate(qqcw_out)
        allocate(qqcw_out(size(qqcw),pver))
        do imm=1,size(qqcw)
           if (size(qqcw(imm)%fld,1) > 10000) then
                qqcw_out(imm,:) = -9999.9
           else
                qqcw_out(imm,:) = qqcw(imm)%fld(yaml%col_print,:)
           endif
        enddo     
        call write_var(unit_output,'qqcw',qqcw_out)
        deallocate(qqcw_out)

        call write_var(unit_output,'tauxar',tauxar(yaml%col_print,0:pver,:))
        call write_var(unit_output,'wa',wa(yaml%col_print,0:pver,:))
        call write_var(unit_output,'ga',ga(yaml%col_print,0:pver,:))
        call write_var(unit_output,'fa',fa(yaml%col_print,0:pver,:))

        call write_var(unit_output,'extinct', extinct(yaml%col_print,:))
        call write_var(unit_output,'absorb', absorb(yaml%col_print,:))
        call write_var(unit_output,'aodabs', aodabs(yaml%col_print))
        call write_var(unit_output,'aodvis', aodvis(yaml%col_print))
        call write_var(unit_output,'aodall', aodall(yaml%col_print))
        call write_var(unit_output,'burdendust', burdendust(yaml%col_print))
        call write_var(unit_output,'burdenso4', burdenso4(yaml%col_print))
        call write_var(unit_output,'burdenpom', burdenpom(yaml%col_print))
        call write_var(unit_output,'burdensoa', burdensoa(yaml%col_print))
        call write_var(unit_output,'burdenbc', burdenbc(yaml%col_print))
        call write_var(unit_output,'burdenseasalt', burdenseasalt(yaml%col_print))
        call write_var(unit_output,'burdenmom', burdenmom(yaml%col_print))
        call write_var(unit_output,'dustaod', dustaod(yaml%col_print))
        call write_var(unit_output,'so4aod', so4aod(yaml%col_print))
        call write_var(unit_output,'pomaod', pomaod(yaml%col_print))
        call write_var(unit_output,'soaaod', soaaod(yaml%col_print))
        call write_var(unit_output,'bcaod', bcaod(yaml%col_print))
        call write_var(unit_output,'seasaltaod', seasaltaod(yaml%col_print))
        call write_var(unit_output,'momaod', momaod(yaml%col_print))

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
!if(yaml%flag_print) write(101,*) shape(extpsw), shape(abspsw), shape(asmpsw)
!------------------------------------------------------------------------------------------------------------------------------

#endif
