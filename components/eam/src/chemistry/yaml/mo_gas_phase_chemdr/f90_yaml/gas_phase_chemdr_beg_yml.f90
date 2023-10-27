#ifdef YAML_CHEMDR
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
  ! character(len=200) :: ext_str
  !-----------------------------------------------------------------------------------------

  type(yaml_vars) :: yaml
  integer  :: unit_input, unit_output, y_nstep
 !  print output from timestep nstep_print_lo to yaml%nstep_print
  integer  :: nstep_print_lo

  ! some subroutines are called multiple times in one timestep, record the number of calls
  integer,save :: n_calls=0
  integer,save :: y_nstep_old=0

  ! variables to store data referenced by inout pointer qqcw
  ! qqcw_out is declared here but assigned in the _end_yml.f90 include file 
  ! these are valid for a particular comment, so dimensions are ( ,pver)
  ! for first dimension, elements beyond 15 are set to -9999

  real(r8), allocatable :: qqcw_in(:,:)
  real(r8), allocatable :: qqcw_out(:,:)
  integer imm

  !populate YAML structure
  !(**remove yaml%lev_print, nstep_print, col_print if generating data for a dependent subroutines**)
  ! used with phys_debug_lat = 39.0553303768561, phys_debug_lon = 262.904774388088

  yaml%lev_print = 30       !print level above tropopause
  yaml%nstep_print = 1417 !time step
  nstep_print_lo = 1400
  yaml%col_print = icolprnt(lchnk)                !column to write data

  !current_time step
  y_nstep = get_nstep()

  !Flag to decide to write or not to write data
  yaml%flag_print = .false. !(**remove these if generating data for a dependent subroutines**)

!!  if(yaml%col_print == y_i .and. y_nstep==yaml%nstep_print .and. y_k == yaml%lev_print) then ! if this column exists in y_lchnk

  !-----------------------------------------------------------------------------------------
  !In the case of y_i or y_k are not passed as arguments, use the following if condition:
  !if(yaml%col_print >0 .and. y_nstep==yaml%nstep_print) then
  !  print output for timesteps between nstep_print_lo and yaml%nstep_print
  !  output depends on multiple vertical levels with internal kk loop, so use below condition
  if(yaml%col_print >0 .and. y_nstep<=yaml%nstep_print .and. y_nstep>=nstep_print_lo) then
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
     !-----------------------------------------------------------------------------------------


     !Record number of calls that can output yaml file if you only need to write one set of input/output
     !Increment n_calls if current timestep, y_nstep, is same as for previous pass through this code.
     !Otherwise, we are at a new timestep, can reset n_calls to 1.
     if(y_nstep == y_nstep_old) then
        n_calls = n_calls+1
     else
        n_calls = 1
     endif
     y_nstep_old = y_nstep

     if (n_calls==1) then ! output at the first call only, modify this if condition (see below) if writing out other calls

     !if ((n_calls==1 .and. flag==0) .or. (n_calls==3 .and. flag==1) .or. (n_calls==5 .and. flag==2)) then

        !(**remove these yaml% variables if generating data for a dependent subroutines**)
        yaml%lchnk_print = lchnk
        yaml%flag_print  = .true.


        !open I/O yaml files
        !(with an optional argument to pass a unique string to differentiate file names)
        call open_files('gas_phase_chemdr', &  !intent-in
             unit_input, unit_output) !intent-out
        !    unit_input, unit_output, trim(ext_str)) !intent-out, with the use of ext_str


        !start by adding an input string
        !Use current timestep, given by y_nstep, to form name of input / output files
        call write_input_output_header(unit_input, unit_output,yaml%lchnk_print,yaml%col_print, &
             'gas_phase_chemdr',y_nstep, yaml%lev_print)

        ! add code for writing data here
        ! below are explicit input arguments to gas_phase_chemdr
        
        call write_var(unit_input,unit_output,'lchnk',lchnk)
        call write_var(unit_input,unit_output,'ncol',ncol)
        call write_var(unit_input,unit_output,'imozart',imozart)
        call write_var(unit_input,unit_output,'state_q',state_q(yaml%col_print,:,:))
        call write_var(unit_input,unit_output,'phis',phis(yaml%col_print))
        call write_var(unit_input,unit_output,'zm',zm(yaml%col_print,:))
        call write_var(unit_input,unit_output,'zi',zi(yaml%col_print,:))
        call write_var(unit_input,unit_output,'calday',calday)
        call write_var(unit_input,unit_output,'tfld',tfld(yaml%col_print,:))
        call write_var(unit_input,unit_output,'pmid',pmid(yaml%col_print,:))
        call write_var(unit_input,unit_output,'pdel',pdel(yaml%col_print,:))
        call write_var(unit_input,unit_output,'pdeldry',pdeldry(yaml%col_print,:))
        call write_var(unit_input,unit_output,'pint',pint(yaml%col_print,:))
        call write_var(unit_input,unit_output,'cldw',cldw(yaml%col_print,:))
        call write_var(unit_input,unit_output,'troplev',troplev(yaml%col_print))
        call write_var(unit_input,unit_output,'ncldwtr',ncldwtr(yaml%col_print,:))
        call write_var(unit_input,unit_output,'ufld',ufld(yaml%col_print,:))
        call write_var(unit_input,unit_output,'vfld',vfld(yaml%col_print,:))
        call write_var(unit_input,unit_output,'prain',prain(yaml%col_print,:))
        call write_var(unit_input,unit_output,'cldfr',cldfr(yaml%col_print,:))
        call write_var(unit_input,unit_output,'cmfdqr',cmfdqr(yaml%col_print,:))
        call write_var(unit_input,unit_output,'nevapr',nevapr(yaml%col_print,:))
        call write_var(unit_input,unit_output,'delt',delt)
        call write_var(unit_input,unit_output,'ps',ps(yaml%col_print))
        call write_var(unit_input,unit_output,'linoz_o3_clim',linoz_o3_clim(yaml%col_print,:))
        call write_var(unit_input,unit_output,'linoz_t_clim',linoz_t_clim(yaml%col_print,:))
        call write_var(unit_input,unit_output,'linoz_o3col_clim',linoz_o3col_clim(yaml%col_print,:))
        call write_var(unit_input,unit_output,'linoz_PmL_clim',linoz_PmL_clim(yaml%col_print,:))
        call write_var(unit_input,unit_output,'linoz_dPmL_dO3',linoz_dPmL_dO3(yaml%col_print,:))
        call write_var(unit_input,unit_output,'linoz_dPmL_dT',linoz_dPmL_dT(yaml%col_print,:))
        call write_var(unit_input,unit_output,'linoz_dPmL_dO3col',linoz_dPmL_dO3col(yaml%col_print,:))
        call write_var(unit_input,unit_output,'linoz_cariolle_psc',linoz_cariolle_psc(yaml%col_print,:))
        call write_var(unit_input,unit_output,'fsds',fsds(yaml%col_print))
        call write_var(unit_input,unit_output,'ts',ts(yaml%col_print))
        call write_var(unit_input,unit_output,'asdir',asdir(yaml%col_print))
        call write_var(unit_input,unit_output,'precc',precc(yaml%col_print))
        call write_var(unit_input,unit_output,'precl',precl(yaml%col_print))
        call write_var(unit_input,unit_output,'snowhland',snowhland(yaml%col_print))
        call write_var(unit_input,unit_output,'pblh',pblh(yaml%col_print))
        call write_var(unit_input,unit_output,'drydepflx',drydepflx(yaml%col_print,:))
        call write_var(unit_input,unit_output,'cflx',cflx(yaml%col_print,:))
        call write_var(unit_input,unit_output,'qtend',qtend(yaml%col_print,:,:))
        call write_var(unit_input,unit_output,'dgnum',dgnum(yaml%col_print,:,:))
        call write_var(unit_input,unit_output,'dgnumwet',dgnumwet(yaml%col_print,:,:))
        call write_var(unit_input,unit_output,'wetdens',wetdens(yaml%col_print,:,:))

        if (allocated(qqcw_in)) deallocate(qqcw_in)
        allocate(qqcw_in(size(qqcw),pver))
        do imm=1,size(qqcw)
           if (imm<=15) then
                qqcw_in(imm,:) = -9999.9
           else
                qqcw_in(imm,:) = qqcw(imm)%fld(yaml%col_print,:)
           endif
        enddo
        call write_var(unit_input,unit_output,'qqcw',qqcw_in)
        deallocate(qqcw_in)

        ! below are external module variable inputs to gas_phase_chemdr
        call write_var(unit_input,unit_output,'pi',pi)
        call write_var(unit_input,unit_output,'fieldname_len',fieldname_len)
        call write_var(unit_input,unit_output,'pcols',pcols)
        call write_var(unit_input,unit_output,'pver',pver)
        call write_var(unit_input,unit_output,'pcnst',pcnst)
        call write_var(unit_input,unit_output,'nabscol',nabscol)
        call write_var(unit_input,unit_output,'nfs',nfs)
        call write_var(unit_input,unit_output,'indexm',indexm)
        call write_var(unit_input,unit_output,'rga',rga)
        call write_var(unit_input,unit_output,'lambm0',lambm0)
        call write_var(unit_input,unit_output,'eccen',eccen)
        call write_var(unit_input,unit_output,'mvelpp',mvelpp)
        call write_var(unit_input,unit_output,'obliqr',obliqr)
        call write_var(unit_input,unit_output,'phtcnt',phtcnt)
        call write_var(unit_input,unit_output,'rxntot',rxntot)
        call write_var(unit_input,unit_output,'gas_pcnst',gas_pcnst)
        call write_var(unit_input,unit_output,'rxt_tag_cnt',rxt_tag_cnt)
        call write_var(unit_input,unit_output,'rxt_tag_map',rxt_tag_map)
        call write_var(unit_input,unit_output,'extcnt',extcnt)

        ! below are internal module variable inputs to gas_phase_chemdr
        call write_var(unit_input,unit_output,'map2chm',map2chm(:))
        call write_var(unit_input,unit_output,'synoz_ndx',synoz_ndx)    
        call write_var(unit_input,unit_output,'o3_ndx',o3_ndx)    
        call write_var(unit_input,unit_output,'ndx_h2so4',ndx_h2so4)    
        call write_var(unit_input,unit_output,'inv_ndx_cnst_o3',inv_ndx_cnst_o3)    
        call write_var(unit_input,unit_output,'rxn_names',rxn_names(:))    
        call write_var(unit_input,unit_output,'pht_names',pht_names(:))    
        call write_var(unit_input,unit_output,'tag_names',tag_names(:))    
        call write_var(unit_input,unit_output,'extfrc_name',extfrc_name(:))    


        !writes aerosol mmr from state%q or q vector (cloud borne and interstitial)
        !"aer_num_only" is .ture. if printing aerosol num only
        !call write_aerosol_mmr_from_stateq(unit_input, unit_output, fld_name,field,aer_num_only)

        !close only the input file, not the output file
        close(unit_input)
        call freeunit(unit_input)
     endif
  endif
#endif
