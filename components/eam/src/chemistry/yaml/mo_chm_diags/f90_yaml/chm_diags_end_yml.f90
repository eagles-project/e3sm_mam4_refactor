#ifdef YAML_CHMDIAGS
  if(yaml%flag_print) then ! if this column exists in lchnk

     !write output header
     call write_output_header(unit_output)

     !start writing data
     
        call write_var(unit_output,'area',area(yaml%col_print))
        call write_var(unit_output,'mass',mass(yaml%col_print,:))
        call write_var(unit_output,'drymass',drymass(yaml%col_print,:))
        call write_var(unit_output,'ozone_col',ozone_col(yaml%col_print))
        call write_var(unit_output,'ozone_strat',ozone_strat(yaml%col_print))
        call write_var(unit_output,'ozone_trop',ozone_trop(yaml%col_print))
        call write_var(unit_output,'net_chem',net_chem(yaml%col_print,:))
        call write_var(unit_output,'mass_bc',mass_bc(yaml%col_print,:))
        call write_var(unit_output,'mass_dst',mass_dst(yaml%col_print,:))
        call write_var(unit_output,'mass_mom',mass_mom(yaml%col_print,:))
        call write_var(unit_output,'mass_ncl',mass_ncl(yaml%col_print,:))
        call write_var(unit_output,'mass_pom',mass_pom(yaml%col_print,:))
        call write_var(unit_output,'mass_so4',mass_so4(yaml%col_print,:))
        call write_var(unit_output,'mass_soa',mass_soa(yaml%col_print,:))
        call write_var(unit_output,'vmr_nox',vmr_nox(yaml%col_print,:))
        call write_var(unit_output,'vmr_noy',vmr_noy(yaml%col_print,:))
        call write_var(unit_output,'vmr_clox',vmr_clox(yaml%col_print,:))
        call write_var(unit_output,'vmr_cloy',vmr_cloy(yaml%col_print,:))
        call write_var(unit_output,'vmr_brox',vmr_brox(yaml%col_print,:))
        call write_var(unit_output,'vmr_broy',vmr_broy(yaml%col_print,:))
        call write_var(unit_output,'vmr_tcly',vmr_tcly(yaml%col_print,:))
        call write_var(unit_output,'mmr_noy',mmr_noy(yaml%col_print,:))
        call write_var(unit_output,'mmr_sox',mmr_sox(yaml%col_print,:))
        call write_var(unit_output,'mmr_nhx',mmr_nhx(yaml%col_print,:))
        call write_var(unit_output,'df_noy',df_noy(yaml%col_print))
        call write_var(unit_output,'df_sox',df_sox(yaml%col_print))
        call write_var(unit_output,'df_nhx',df_nhx(yaml%col_print))

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
