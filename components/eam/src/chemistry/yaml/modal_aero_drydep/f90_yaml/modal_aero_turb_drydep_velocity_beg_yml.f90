#ifdef YAML_AERO_DRYDEP
  !call
  !call modal_aero_turb_drydep_velocity( imode,yaml, jvlc, moment, ncol, pcols, lchnk, radius_max,    &! in
  !subroutine:
  !subroutine modal_aero_turb_drydep_velocity( imode,yaml, jvlc,moment, ncol, pcols, lchnk, radius_max,  &! in

 integer, intent(in) :: imode, jvlc
  !-----------------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------------------
  ! This is used when multiple sets of yaml output is needed
  !to cover different options (e.g., true and false)
 character(len=200) :: ext_str, mode_str
  !-----------------------------------------------------------------------------------------

  type(yaml_vars) :: yaml
  integer  :: unit_input, unit_output
  !-----------------------------------------------------------------------------------------
  !For generating data for a dependent subroutines where "yaml" derived type is already
  !initialized, use the following if condition
  if(yaml%flag_print) then
  !-----------------------------------------------------------------------------------------

     !-----------------------------------------------------------------------------------------
     ! Set "ext_str" if there are multiple sets of yaml output to be written out
     ! Example:"flag" in the code can be 0, 1, or 2, we can update "ext_str" as:
   if (imode == -999) then
      write(mode_str,'(A)')'_imode_No_mode'
   else
      write(mode_str,'(A,I1)')'_imode_',imode
   endif
   write(ext_str,'(A,I1,A,I1,A,I1)')'jvlc_',jvlc,'_imnt_',moment,trim(mode_str)
     !-----------------------------------------------------------------------------------------

    !open I/O yaml files
        !(with an optional argument to pass a unique string to differentiate file names)
        call open_files('modal_aero_turb_drydep_velocity', &  !intent-in
         !    unit_input, unit_output) !intent-out
            unit_input, unit_output, trim(ext_str)) !intent-out, with the use of ext_str


        !start by adding an input string
        call write_input_output_header(unit_input, unit_output,yaml%lchnk_print,yaml%col_print, &
             'modal_aero_turb_drydep_velocity',yaml%nstep_print, yaml%lev_print)

        ! add code for writing data here
        call write_var(unit_input,unit_output,'fraction_landuse',fraction_landuse(yaml%col_print,:,lchnk))
        call write_var(unit_input,unit_output,'n_land_type',n_land_type)
        call write_var(unit_input,unit_output,'moment',moment)
        call write_var(unit_input,unit_output,'ncol',ncol)
        call write_var(unit_input,unit_output,'pcols',pcols)
        call write_var(unit_input,unit_output,'lchnk',lchnk)
        call write_var(unit_input,unit_output,'radius_max',radius_max)
        call write_var(unit_input,unit_output,'tair',tair(yaml%col_print))
        call write_var(unit_input,unit_output,'pmid',pmid(yaml%col_print))
        call write_var(unit_input,unit_output,'radius_part',radius_part(yaml%col_print))
        call write_var(unit_input,unit_output,'density_part',density_part(yaml%col_print))
        call write_var(unit_input,unit_output,'sig_part',sig_part(yaml%col_print))
        call write_var(unit_input,unit_output,'fricvel',fricvel(yaml%col_print))
        call write_var(unit_input,unit_output,'ram1',ram1(yaml%col_print))
        call write_var(unit_input,unit_output,'vlc_grv',vlc_grv(yaml%col_print))

        !close only the input file, not the output file
        close(unit_input)
        call freeunit(unit_input)
     endif
#endif