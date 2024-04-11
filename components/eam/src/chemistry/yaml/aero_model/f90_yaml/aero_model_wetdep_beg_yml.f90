#ifdef YAML_AERO_MODEL_WETDEP

  
  !-----------------------------------------------------------------------------------------
  !"lchnk" is needed for the following code to work,
  ! temporarily pass it along from upper level subroutines
  ! as y_lchnk and uncomment the following code:
  ! integer, intent(in) :: y_lchnk
    integer :: y_lchnk
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

  ! some subroutines are called multiple times in one timestep, record the number of calls
  integer,save :: n_calls=0

  logical, allocatable :: ptend_lq_1darr(:)   ! logical flag for tracer tendency (ptend%lq)
  real(r8), allocatable :: ptend_q_3darr(:,:,:)        ! tracer tendency (ptend%q) [kg/kg/s]

  real(r8), pointer :: dgncur_a_yaml(:,:,:)
  real(r8), pointer :: wetdens_yaml(:,:,:)
  real(r8), pointer :: qaerwat_yaml(:,:,:)
  real(r8), pointer :: dgnumwet_yaml(:,:,:)
  real(r8), pointer :: fracis_yaml(:,:,:)
  real(r8), pointer :: rprddp_yaml(:,:)
  real(r8), pointer :: rprdsh_yaml(:,:)
  real(r8), pointer :: evapcsh_yaml(:,:)
  real(r8), pointer :: evapcdp_yaml(:,:)

  type(ptr2d_t) :: qqcw_yaml(pcnst)                 !cloud-borne aerosols mass and number mixing rations
  real(r8), allocatable :: qqcw_yaml_2darr(:,:)
  integer :: imm

  !populate YAML structure
  !(**remove yaml%lev_print, nstep_print, col_print if generating data for a dependent subroutines**)
  yaml%lev_print = 68       !level
  yaml%nstep_print = 379 !time step

  y_lchnk = state%lchnk
  yaml%col_print = icolprnt(y_lchnk)                !column to write data

  !current_time step
  y_nstep = get_nstep()

  !Flag to decide to write or not to write data
  yaml%flag_print = .false. !(**remove these if generating data for a dependent subroutines**)

  !if(yaml%col_print == y_i .and. y_nstep==yaml%nstep_print .and. y_k == yaml%lev_print) then ! if this column exists in y_lchnk
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
     !-----------------------------------------------------------------------------------------


     !Record number of calls that can output yaml file if you only need to write one set of input/output
     n_calls = n_calls+1

     if (n_calls==1) then ! output at the first call only, modify this if condition (see below) if writing out other calls

     !if ((n_calls==1 .and. flag==0) .or. (n_calls==3 .and. flag==1) .or. (n_calls==5 .and. flag==2)) then

        !(**remove these yaml% variables if generating data for a dependent subroutines**)
        yaml%lchnk_print = y_lchnk
        yaml%flag_print  = .true.


        !open I/O yaml files
        !(with an optional argument to pass a unique string to differentiate file names)
        call open_files('aero_model_wetdep', &  !intent-in
             unit_input, unit_output) !intent-out
        !    unit_input, unit_output, trim(ext_str)) !intent-out, with the use of ext_str


        !start by adding an input string
        call write_input_output_header(unit_input, unit_output,yaml%lchnk_print,yaml%col_print, &
             'aero_model_wetdep',yaml%nstep_print, yaml%lev_print)

        ! add code for writing data here
        
        ! below are explicit input arguments to aero_model_wetdep

        call write_var(unit_input,unit_output,'dt',dt)
        call write_var(unit_input,unit_output,'dlf',dlf(yaml%col_print,yaml%lev_print)) 
        call write_var(unit_input,unit_output,'dlf2',dlf2(yaml%col_print,yaml%lev_print))
        call write_var(unit_input,unit_output,'cmfmc2',cmfmc2(yaml%col_print,yaml%lev_print))
        call write_var(unit_input,unit_output,'sh_e_ed_ratio',sh_e_ed_ratio(yaml%col_print,yaml%lev_print))
        call write_var(unit_input,unit_output,'mu',mu(yaml%col_print,yaml%lev_print))
        call write_var(unit_input,unit_output,'md',md(yaml%col_print,yaml%lev_print))
        call write_var(unit_input,unit_output,'du',du(yaml%col_print,yaml%lev_print))
        call write_var(unit_input,unit_output,'eu',eu(yaml%col_print,yaml%lev_print))
        call write_var(unit_input,unit_output,'ed',ed(yaml%col_print,yaml%lev_print))
        call write_var(unit_input,unit_output,'dp',dp(yaml%col_print,yaml%lev_print))
        call write_var(unit_input,unit_output,'jt',jt(yaml%col_print))
        call write_var(unit_input,unit_output,'maxg',maxg(yaml%col_print))
        call write_var(unit_input,unit_output,'ideep',ideep(yaml%col_print))
        call write_var(unit_input,unit_output,'lengath',lengath)
        call write_var(unit_input,unit_output,'species_class',species_class(:))


        !  the following are input arguments that are fields within the state vector.
        !  they are written out to datafiles here before being de-referenced in the subroutine body, using the same nomenclature

        call write_var(unit_input,unit_output,'lchnk',y_lchnk) 
        call write_var(unit_input,unit_output,'ncol',state%ncol)
        call write_var(unit_input,unit_output,'state_q',state%q(yaml%col_print,yaml%lev_print,:)) ! ??? distinct from state_q used later on?
        call write_var(unit_input,unit_output,'temperature',state%t(yaml%col_print,yaml%lev_print))
        call write_var(unit_input,unit_output,'pmid',state%pmid(yaml%col_print,yaml%lev_print))
        call write_var(unit_input,unit_output,'pdel',state%pdel(yaml%col_print,yaml%lev_print)) 

        ! the following are input arguments that are fields within pbuf.  We thus need to duplicate the pbuf extraction subroutine calls in order to write out the fields.

        call pbuf_get_field(pbuf, dgnum_idx,      dgncur_a_yaml )  ! note:  this is called again later on before passed into dgnumdry_m in modal_aero_calcsize_sub
        call pbuf_get_field(pbuf, wetdens_ap_idx, wetdens_yaml)
        call pbuf_get_field(pbuf, qaerwat_idx,    qaerwat_yaml)
        call pbuf_get_field(pbuf, dgnumwet_idx, dgnumwet_yaml, start=(/1,1,1/), kount=(/pcols,pver,nmodes/) )
        call pbuf_get_field(pbuf, fracis_idx,   fracis_yaml,   start=(/1,1,1/), kount=(/pcols,pver, pcnst/) )
        call pbuf_get_field(pbuf, rprddp_idx,      rprddp_yaml  )
        call pbuf_get_field(pbuf, rprdsh_idx,      rprdsh_yaml  )
        call pbuf_get_field(pbuf, nevapr_shcu_idx, evapcsh_yaml )
        call pbuf_get_field(pbuf, nevapr_dpcu_idx, evapcdp_yaml )

        call write_var(unit_input,unit_output,'dgncur_a',dgncur_a_yaml(yaml%col_print,yaml%lev_print,:)) 
        call write_var(unit_input,unit_output,'wetdens',wetdens_yaml(yaml%col_print,yaml%lev_print,:))  
        call write_var(unit_input,unit_output,'qaerwat',qaerwat_yaml(yaml%col_print,yaml%lev_print,:)) 
        call write_var(unit_input,unit_output,'dgnumwet',dgnumwet_yaml(yaml%col_print,yaml%lev_print,:))  
        call write_var(unit_input,unit_output,'fracis',fracis_yaml(yaml%col_print,yaml%lev_print,:)) 
        call write_var(unit_input,unit_output,'rprddp',rprddp_yaml(:,:)) 
        call write_var(unit_input,unit_output,'rprdsh',rprdsh_yaml(:,:)) 
        call write_var(unit_input,unit_output,'evapcsh',evapcsh_yaml(:,:)) 
        call write_var(unit_input,unit_output,'evapcdp',evapcdp_yaml(:,:)) 

        call get_cldbrn_mmr(lchnk, pbuf, &! in
         qqcw_yaml) !out      
        
        ! write out qqcw fields only over constituents that are set, for others set value to -9999.9
  
        if (allocated(qqcw_yaml_2darr)) deallocate(qqcw_yaml_2darr)
        allocate(qqcw_yaml_2darr(size(qqcw_yaml),pver))
        do imm=1,size(qqcw_yaml)
           if (imm<=15) then
                qqcw_yaml_2darr(imm,:) = -9999.9
           else
                qqcw_yaml_2darr(imm,:) = qqcw_yaml(imm)%fld(yaml%col_print,:)
           endif
        enddo
        call write_var(unit_input,unit_output,'qqcw',qqcw_yaml_2darr)
        deallocate(qqcw_yaml_2darr)

        !writes aerosol mmr from state%q or q vector (cloud borne and interstitial)
        !"aer_num_only" is .ture. if printing aerosol num only
        !call write_aerosol_mmr_from_stateq(unit_input, unit_output, fld_name,field,aer_num_only)

        !close only the input file, not the output file
        close(unit_input)
        call freeunit(unit_input)



     endif
  endif
#endif
