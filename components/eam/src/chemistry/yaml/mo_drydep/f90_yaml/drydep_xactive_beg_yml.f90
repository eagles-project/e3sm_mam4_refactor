#ifdef YAML_GAS_DRYDEP
  !-----------------------------------------------------------------------------------------
  !"lchnk" is needed for the following code to work,
  ! temporarily pass it along from upper level subroutines
  ! as y_lchnk and uncomment the following code:
  ! integer, intent(in) :: y_lchnk
  !-----------------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------------------
  ! **OR** if this subroutine is called in a nested loop of columns and levels,
  ! we also might need column index (y_i or icol) and level
  ! index (y_k or klev) to be passed to this routine and
  ! uncomment the following code:
  ! integer, intent(in) :: y_i, y_k, y_lchnk
  !-----------------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------------------
  ! This is used when multiple sets of yaml output is needed
  !to cover different options (e.g., true and false)
  character(len=200) :: ext_str
  !-----------------------------------------------------------------------------------------

  type(yaml_vars) :: yaml
  integer  :: unit_input, unit_output, y_nstep

  ! some subroutines are called multiple times in one timestep, record the number of calls
  integer,save :: n_calls=0


  call seq_drydep_setHCoeff( ncol, sfc_temp, & ! in
                               heff )            ! out

  !populate YAML structure
  !(**remove yaml%lev_print, nstep_print, col_print if generating data for a dependent subroutines**)
  yaml%lev_print = 72       !level
  yaml%nstep_print = 1000 !time step
!  yaml%nstep_print = 100 !time step
!  yaml%nstep_print = 300 !time step

  yaml%col_print = icolprnt(lchnk)                !column to write data

  !current_time step
  y_nstep = get_nstep()

  !Flag to decide to write or not to write data
  yaml%flag_print = .false. !(**remove these if generating data for a dependent subroutines**)

  !if(yaml%col_print == y_i .and. y_nstep==yaml%nstep_print .and. y_k == yaml%lev_print) then ! if this column exists in y_lchnk

  !do icol = 1, ncol
  !   if (snow(icol) > .02_r8 .and. y_nstep == 300) then
  !      write(102,*)'phys_debug_lat = ',get_lat(y_lchnk, icol), &
  !                  ' phys_debug_lon = ', get_lon(y_lchnk, icol), &
  !                  ' snow = ', snow(icol)
  !   endif
  !enddo

  !-----------------------------------------------------------------------------------------
  !In the case of y_i or y_k are not passed as arguments, use the following if condition:
  if(yaml%col_print >0 .and. y_nstep==yaml%nstep_print) then
  !-----------------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------------------
  !For generating data for a dependent subroutines where "yaml" derived type is already
  !initialized, use the following if condition
  !if(yaml%flag_print) then
  !-----------------------------------------------------------------------------------------

     !-----------------------------------------------------------------------------------------
     ! Set "ext_str" if there are multiple sets of yaml output to be written out
     ! Example:"flag" in the code can be 0, 1, or 2, we can update "ext_str" as:
     ! write(ext_str,'(I2)') flag
     ! ext_str = 'flag_'//adjustl(ext_str)
      ext_str = "loc3_over_water"
     ! ext_str = "loc2_has_snow"
     !-----------------------------------------------------------------------------------------


     !Record number of calls that can output yaml file if you only need to write one set of input/output
     n_calls = n_calls+1

     if (n_calls==1) then ! output at the first call only, modify this if condition (see below) if writing out other calls

     !if ((n_calls==1 .and. flag==0) .or. (n_calls==3 .and. flag==1) .or. (n_calls==5 .and. flag==2)) then

        !(**remove these yaml% variables if generating data for a dependent subroutines**)
        yaml%lchnk_print = lchnk
        yaml%flag_print  = .true.


        !open I/O yaml files
        !(with an optional argument to pass a unique string to differentiate file names)
        call open_files('drydep_xactive', &  !intent-in
             unit_input, unit_output, trim(ext_str)) !intent-out, with the use of ext_str
!             unit_input, unit_output) !intent-out


        !start by adding an input string
        call write_input_output_header(unit_input, unit_output,yaml%lchnk_print,yaml%col_print, &
             'drydep_xactive',yaml%nstep_print, yaml%lev_print)

        ! add code for writing data here
        
        call write_var(unit_input,unit_output,'lchnk',lchnk)
        call write_var(unit_input,unit_output,'ncol',ncol)
        call write_var(unit_input,unit_output,'latndx',latndx(yaml%col_print))
        call write_var(unit_input,unit_output,'ncdate',ncdate)
        call write_var(unit_input,unit_output,'sfc_temp',sfc_temp(yaml%col_print))
        call write_var(unit_input,unit_output,'air_temp',air_temp(yaml%col_print))
        call write_var(unit_input,unit_output,'tv',tv(yaml%col_print))
        call write_var(unit_input,unit_output,'pressure_sfc',pressure_sfc(yaml%col_print))
        call write_var(unit_input,unit_output,'pressure_10m',pressure_10m(yaml%col_print))
        call write_var(unit_input,unit_output,'spec_hum',spec_hum(yaml%col_print))
        call write_var(unit_input,unit_output,'wind_speed',wind_speed(yaml%col_print))
        call write_var(unit_input,unit_output,'rain',rain(yaml%col_print))
        call write_var(unit_input,unit_output,'snow',snow(yaml%col_print))
        call write_var(unit_input,unit_output,'solar_flux',solar_flux(yaml%col_print))
        call write_var(unit_input,unit_output,'mmr',mmr(yaml%col_print,yaml%lev_print,:))
        call write_var(unit_input,unit_output,'rain_threshold',rain_threshold)
        call write_var(unit_input,unit_output,'temp_highbound',temp_highbound)
        call write_var(unit_input,unit_output,'ric',ric)
        call write_var(unit_input,unit_output,'heff',heff(yaml%col_print,:))

        
        call write_var(unit_input,unit_output,'fraction_landuse',fraction_landuse(yaml%col_print,:,lchnk))
        call write_var(unit_input,unit_output,'n_land_type',n_land_type)
        call write_var(unit_input,unit_output,'index_season_lai',index_season_lai(:,:))

        !writes aerosol mmr from state%q or q vector (cloud borne and interstitial)
        !"aer_num_only" is .ture. if printing aerosol num only
        !call write_aerosol_mmr_from_stateq(unit_input, unit_output, fld_name,field,aer_num_only)

        !close only the input file, not the output file
        close(unit_input)
        call freeunit(unit_input)
     endif
  endif
#endif
