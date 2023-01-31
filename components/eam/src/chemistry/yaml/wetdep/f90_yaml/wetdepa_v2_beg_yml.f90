#ifdef YAML_WETDEP
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
  character(len=200) :: ext_str
  type(yaml_vars) :: yaml
  integer  :: unit_input, unit_output, y_nstep
  integer,save :: n_calls=0   ! some subroutines are called multiple times in one timestep, record the number of calls

  !populate YAML structure
  yaml%lev_print = 0 !level (**remove these if generating data for a dependent !  subroutines**)
  yaml%nstep_print = 355 !time step(**remove these if generating data for a dependent subroutines**)


  write(ext_str,'(I3)') mam_prevap_resusp_optcc
  ext_str = 'mam_prevap_resusp_optcc_'//adjustl(ext_str)
!  write(ext_str, '(A,I3)') 'mam_prevap_resusp_optcc_', mam_prevap_resusp_optcc
!  if (mam_prevap_resusp_optcc == 0) then
!       ext_str = 'mam_prevap_resusp_optcc_0'
!  elseif (mam_prevap_resusp_optcc == 130) then
!       ext_str = 'mam_prevap_resusp_optcc_130'
!  else
!       ext_str = 'mam_prevap_resusp_optcc_230'
!  endif

  !YAML file input generation code- DO NOT PORT to C++
  !print all inputs one-by-one at column "yaml%col_print"
  yaml%col_print = icolprnt(y_lchnk) !column to write data(**remove these if generating data for a dependent subroutines**)
  y_nstep = get_nstep() !time step (**remove these if generating data for a dependent subroutines**)

  yaml%flag_print = .false. ! to write or not to write data (**remove these if generating data for a dependent subroutines**)
!  if(yaml%col_print == y_i .and. y_nstep==yaml%nstep_print .and. y_k == yaml%lev_print) then ! if this column exists in y_lchnk
     !<
     !In the case of y_i or y_k are not passed as arguments:
  if(yaml%col_print >0 .and. y_nstep==yaml%nstep_print) then
     !For generating data for a dependent subroutines where "yaml" derived type is already initialized:
     !if(yaml%flag_print) then
     !>
     n_calls = n_calls+1
     if ((n_calls==8 .and. mam_prevap_resusp_optcc==0) .or. (n_calls==1 .and. mam_prevap_resusp_optcc==130) .or. (n_calls==50 .and. mam_prevap_resusp_optcc==230)) then
!     if ((n_calls==1) then  ! output the first call
!     if (n_calls>=0) then  ! output all the calls

     yaml%lchnk_print = y_lchnk !(**remove these if generating data for a dependent subroutines**)
     yaml%flag_print  = .true.!(**remove these if generating data for a dependent subroutines**)

     !open I/O yaml files (it can have an extra optional argument to pass a unique string to differentiate file names)
     call open_files('wetdepa_v2', &  !intent-in
          unit_input, unit_output, trim(ext_str)) !intent-out

     !start by adding an input string
     call write_input_output_header(unit_input, unit_output,yaml%lchnk_print,yaml%col_print, &
          'compute_tendencies',yaml%nstep_print, yaml%lev_print)

     !< add code for writing data here>
call write_var(unit_input, unit_output, 'ncol', ncol)
call write_var(unit_input, unit_output, 'deltat', deltat)
call write_1d_var(unit_input, unit_output, 'pdel',pver, pdel(yaml%col_print,:))
call write_1d_var(unit_input, unit_output, 'cmfdqr',pver, cmfdqr(yaml%col_print,:))
call write_1d_var(unit_input, unit_output, 'evapc',pver, evapc(yaml%col_print,:))
call write_1d_var(unit_input, unit_output, 'dlf',pver, dlf(yaml%col_print,:))
call write_1d_var(unit_input, unit_output, 'conicw',pver, conicw(yaml%col_print,:))
call write_1d_var(unit_input, unit_output, 'precs',pver, precs(yaml%col_print,:))
call write_1d_var(unit_input, unit_output, 'evaps',pver, evaps(yaml%col_print,:))
call write_1d_var(unit_input, unit_output, 'cwat',pver, cwat(yaml%col_print,:))
call write_1d_var(unit_input, unit_output, 'cldt',pver, cldt(yaml%col_print,:))
call write_1d_var(unit_input, unit_output, 'cldc',pver, cldc(yaml%col_print,:))
call write_1d_var(unit_input, unit_output, 'cldvcu',pver, cldvcu(yaml%col_print,:))
call write_1d_var(unit_input, unit_output, 'cldvst',pver, cldvst(yaml%col_print,:))
call write_1d_var(unit_input, unit_output, 'tracer',pver, tracer(yaml%col_print,:))
call write_var(unit_input, unit_output, 'mam_prevap_resusp_optcc',mam_prevap_resusp_optcc)
call write_var(unit_input, unit_output, 'is_strat_cloudborne', is_strat_cloudborne)
call write_1d_var(unit_input, unit_output, 'scavcoef',pver, scavcoef(yaml%col_print,:))
call write_1d_var(unit_input, unit_output, 'f_act_conv',pver, f_act_conv(yaml%col_print,:))
call write_1d_var(unit_input, unit_output, 'qqcw',pver, qqcw(yaml%col_print,:))
call write_var(unit_input, unit_output, 'sol_factb', sol_factb)
call write_var(unit_input, unit_output, 'sol_facti', sol_facti)
call write_1d_var(unit_input, unit_output, 'sol_factic',pver, sol_factic(yaml%col_print,:))


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
