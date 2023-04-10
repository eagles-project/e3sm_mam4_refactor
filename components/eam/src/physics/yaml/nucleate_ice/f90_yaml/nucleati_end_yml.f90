#ifdef YAML_NUCICE
  ! YAML file output generation code- DO NOT PORT to C++
  if(yaml%flag_print) then ! if this column exists in lchnk

     !write output header
     call write_output_header(unit_output)

     !start writing data
     call write_output_var(unit_output, 'nuci', nuci)
     call write_output_var(unit_output, 'onihf', onihf)
     call write_output_var(unit_output, 'oniimm', oniimm)
     call write_output_var(unit_output, 'onidep', onidep)
     call write_output_var(unit_output, 'onimey', onimey)

     !<add code for writing data here>
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
!  if(so4_num >= num_threshold .and. dst3_num >= num_threshold .and. cldn > 0._r8) then
!     if((tc <= -35.0_r8) .and. ((relhum*svp_water(tair)/svp_ice(tair)*subgrid) >= 1.2_r8)) then ! use higher RHi threshold

!        if (tc > regm) then
!           if(tc < -40._r8 .and. wbar > 1._r8) then
!              write(102,*)'sub hf phys_debug_lat = ',get_lat(y_lchnk, y_i), &
!                          ' phys_debug_lon = ', get_lon(y_lchnk, y_i), get_nstep(), y_k, nihf
!           else
!              write(103,*)'sub hetero phys_debug_lat = ',get_lat(y_lchnk, y_i), &
!                          ' phys_debug_lon = ', get_lon(y_lchnk, y_i), get_nstep(), y_k, niimm
!           endif

!        else if (tc < regm-5._r8) then
!           write(104,*)'sub hf phys_debug_lat = ',get_lat(y_lchnk, y_i), &
!                       ' phys_debug_lon = ', get_lon(y_lchnk, y_i), get_nstep(), y_k, nihf
!        else
!           if (tc < -40._r8 .and. wbar > 1._r8) then
!              write(105,*)'sub hf phys_debug_lat = ',get_lat(y_lchnk, y_i), &
!                          ' phys_debug_lon = ', get_lon(y_lchnk, y_i), get_nstep(), y_k, nihf
!           else
!              write(106,*)'sub hf and hetero phys_debug_lat = ',get_lat(y_lchnk, y_i), &
!                          ' phys_debug_lon = ', get_lon(y_lchnk, y_i), get_nstep(), y_k, nihf, niimm
!           endif
!        endif
         
!     endif
!  endif


! find n_calls. the if statement below should be consistent with the if condition in the *_beg_yml.f90 file
! so that the corresponding n_calls are printed out here. 
! Note that n_calls may change when lat/lon/y_nstep/y_k change 
!if(yaml%col_print == y_i .and. y_nstep==yaml%nstep_print .and. y_k == yaml%lev_print) then
!      write(103,*) 'n_calls = ',n_calls, <any of the input/output variables to choose n_calls>
!endif
     
#endif
