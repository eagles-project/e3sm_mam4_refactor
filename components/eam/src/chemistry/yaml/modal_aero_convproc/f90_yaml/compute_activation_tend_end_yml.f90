#ifdef YAML_CONV
  ! YAML file output generation code- DO NOT PORT to C++
  if(yaml%flag_print .and. (yaml%lev_print==kk .or. kk==70) .and. (n_calls<=2) ) then ! if this column exists in lchnk

     !write output header
     call write_output_header(unit_output)

     !start writing data
     !<add code for writing data here>
     call write_output_var(unit_output,'kactcnt',kactcnt)
     call write_output_var(unit_output,'kactfirst',kactfirst)
     call write_2d_output_var(unit_output, 'conu', pcnst_extd,pverp, conu)
     call write_2d_output_var(unit_output, 'dconudt_activa', pcnst_extd,pverp, dconudt_activa)
     call write_output_var(unit_output,'xx_wcldbase',xx_wcldbase(yaml%col_print))
     call write_output_var(unit_output,'xx_kcldbase',xx_kcldbase(yaml%col_print))
     !call write_output_var(unit_output, fld_name, field, inp_out_str)  !write a single output variable
     !call write_1d_output_var(unit_output, fld_name, dim, field, inp_out_str) !writes 1D variables of any dimension in the output python module
     !call write_2d_output_var(unit_output, fld_name, dim1, dim2, field, inp_out_str) !writes 2D variables of any dimension in the output python module

     !writes aerosol mmr from state%q or q vector(cloud borne and interstitial) in the output python module
     !"aer_num_only" is .ture. if printing aerosol num only
     !call write_output_aerosol_mmr_from_stateq(unit_output, fld_name, field, aer_num_only, inp_out_str)

     !close the output file
     close(unit_output)
     call freeunit(unit_output)
  endif

! this part of the code is for finding the best lat/lon and n_calls.

! find lag/lon/y_nstep/y_k
!if(yaml%flag_print .and. any(dconudt_activa/=0.0)) then
!  write(102,*) n_calls,kk, do_act_this_lev
!endif

! find n_calls. the if statement below should be consistent with the if condition in the *_beg_yml.f90 file
! so that the corresponding n_calls are printed out here. 
! Note that n_calls may change when lat/lon/y_nstep/y_k change 
!if(yaml%flag_print ) then
!      write(103,*) 'n_calls = ',n_calls, kk, do_act_this_lev
!endif
     
#endif
