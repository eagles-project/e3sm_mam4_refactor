#ifdef YAML_CONV
  !<"lchnk" is needed for the following code to work,
  ! temporarily pass it along from upper level subroutines
  ! as y_lchnk and uncomment the following code:
  ! integer, intent(in) :: y_lchnk

  ! **OR** if this subroutine is called in a nested loop of columns and levels,
  ! we also might need column index (y_i or icol) and level
  ! index (y_k or klev) to be passed to this routine and
  ! uncomment the following code:
  ! integer, intent(in) :: y_i, y_k, y_lchnk

  !>

  character(len=200) :: ext_str  ! this is used when multiple sets of yaml output is needed to cover different options (e.g., true and false)
  type(yaml_vars),intent(in) :: yaml
  integer  :: unit_input, unit_output, y_nstep
  integer,save :: n_calls=0   ! some subroutines are called multiple times in one timestep, record the number of calls

  !populate YAML structure

  ! set ext_str if there are multiple sets of yaml output to be write out
  ! here gives an example that "flag" in the code can be 0, 1, or 2:
  ! write(ext_str,'(I2)') flag
  ! ext_str = 'flag_'//adjustl(ext_str)
  if (do_act_this_lev) then
        ext_str = 'do_act_true'
  else
        ext_str = 'do_act_false'
  endif

  !YAML file input generation code- DO NOT PORT to C++
  !print all inputs one-by-one at column "yaml%col_print"

     !<
     !In the case of y_i or y_k are not passed as arguments:
     !if(yaml%col_print >0 .and. y_nstep==yaml%nstep_print) then

     !For generating data for a dependent subroutines where "yaml" derived type is already initialized:
     if(yaml%flag_print .and. (yaml%lev_print==kk .or. kk==70) ) then
     !>

     ! record number of calls that can output yaml file. you only need to write one set of input/output
     n_calls = n_calls+1
     if (n_calls<=2) then ! output the first call. change this if writing out other calls
     ! if ((n_calls==1 .and. flag==0) .or. (n_calls==3 .and. flag==1) .or. (n_calls==5 .and. flag==2)) then ! this is an example of writing multiple sets of YAML files


        !open I/O yaml files (it can have an extra optional argument to pass a unique string to differentiate file names)
        call open_files('compute_activation_tend', &  !intent-in
!             unit_input, unit_output) !intent-out
             unit_input, unit_output, trim(ext_str)) !intent-out, with the use of ext_str


        !start by adding an input string
        call write_input_output_header(unit_input, unit_output,yaml%lchnk_print,yaml%col_print, &
             'compute_tendencies',yaml%nstep_print, kk)

        !< add code for writing data here>
        call write_var(unit_input, unit_output, 'icol', icol)
        call write_var(unit_input, unit_output, 'kk', kk)
        call write_var(unit_input, unit_output, 'f_ent', f_ent)
        call write_1d_var(unit_input,unit_output,'cldfrac_i', pver,cldfrac_i)
        call write_1d_var(unit_input,unit_output,'rhoair_i', pver,rhoair_i)
        call write_1d_var(unit_input,unit_output,'mu_i', pverp,mu_i)
        call write_1d_var(unit_input,unit_output,'dt_u', pver,dt_u)
        call write_1d_var(unit_input,unit_output,'wup', pver,wup)
        call write_1d_var(unit_input,unit_output,'temperature',pver,temperature(yaml%col_print,:))
        call write_1d_var(unit_input,unit_output,'icwmr',pver,icwmr(yaml%col_print,:))
        call write_var(unit_input, unit_output,'kactcnt',kactcnt)
        call write_var(unit_input, unit_output,'kactfirst',kactfirst)
        call write_2d_var(unit_input, unit_output,'conu', pcnst_extd, pverp,conu)
        call write_2d_var(unit_input, unit_output,'dconudt_activa', pcnst_extd, pverp,dconudt_activa)
        call write_var(unit_input,unit_output,'xx_wcldbase',xx_wcldbase(yaml%col_print))
        call write_var(unit_input,unit_output,'xx_kcldbase',xx_kcldbase(yaml%col_print))
        !call write_var(unit_input, unit_output, fld_name,field)!write a single variable
        !call write_1d_var(unit_input, unit_output, fld_name,dim,field) ! writes 1D variables of any dimension
        !call write_2d_var(unit_input, unit_output, fld_name, dim1, dim2, field) ! writes 2D variables of any dimension: field(dim1,dim2)

        !writes aerosol mmr from state%q or q vector (cloud borne and interstitial)
        !"aer_num_only" is .ture. if printing aerosol num only
        !call write_aerosol_mmr_from_stateq(unit_input, unit_output, fld_name,field,aer_num_only)

        !close only the input file, not the output file
        close(unit_input)
        call freeunit(unit_input)
     endif
  endif
#endif
