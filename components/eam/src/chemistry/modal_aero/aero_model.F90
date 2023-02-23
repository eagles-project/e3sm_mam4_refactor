!===============================================================================
! Modal Aerosol Model
! JUly 2015 B.Singh Added unified convection code
!===============================================================================
module aero_model
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use constituents,   only: pcnst, cnst_name, cnst_get_ind
  use ppgrid,         only: pcols, pver, pverp
  use cam_abortutils, only: endrun
  use cam_logfile,    only: iulog
  use perf_mod,       only: t_startf, t_stopf
  use camsrfexch,     only: cam_in_t, cam_out_t
  use aerodep_flx,    only: aerodep_flx_prescribed
  use physics_types,  only: physics_state, physics_ptend, physics_ptend_init
  use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_get_index, pbuf_set_field
  use phys_control,   only: phys_getopts
  use physconst,      only: gravit, rair, rhoh2o, spec_class_gas
  use mo_constants,   only: pi
  use spmd_utils,     only: masterproc

  use cam_history,    only: outfld, fieldname_len
  use chem_mods,      only: gas_pcnst, adv_mass
  use mo_tracname,    only: solsym

  use modal_aero_data,only: cnst_name_cw, ntot_amode

  implicit none
  private

  public :: aero_model_readnl
  public :: aero_model_register
  public :: aero_model_init
  public :: aero_model_gasaerexch ! create, grow, change, and shrink aerosols.
  public :: aero_model_wetdep     ! aerosol wet removal
  public :: aero_model_emissions  ! aerosol emissions
  public :: calc_1_impact_rate

  public :: drydep_lq, dgnumwet_idx, nmodes, wetdens_ap_idx

 ! Misc private data 
  integer :: nmodes  ! number of modes
  integer :: pblh_idx            = 0
  integer :: dgnum_idx           = 0
  integer :: dgnumwet_idx        = 0
  integer :: rate1_cw2pr_st_idx  = 0  

  integer :: wetdens_ap_idx      = 0
  integer :: qaerwat_idx         = 0

  integer :: fracis_idx          = 0
  integer :: prain_idx           = 0
  integer :: nevapr_idx          = 0
  integer :: rprddp_idx          = 0 
  integer :: rprdsh_idx          = 0 
  integer :: nevapr_shcu_idx     = 0
  integer :: nevapr_dpcu_idx     = 0

  integer :: icwmrdp_idx        = 0
  integer :: icwmrsh_idx        = 0
  integer :: sh_frac_idx        = 0
  integer :: dp_frac_idx        = 0

  integer :: imozart             = -1 

  ! variables for table lookup of aerosol impaction/interception scavenging rates
  ! These are made public to be used by MMF w/ ECPP
  integer, public, parameter :: nimptblgrow_mind=-7, nimptblgrow_maxd=12
  real(r8), public :: dlndg_nimptblgrow  = log( 1.25_r8 )
  real(r8), public :: scavimptblnum(nimptblgrow_mind:nimptblgrow_maxd, ntot_amode)
  real(r8), public :: scavimptblvol(nimptblgrow_mind:nimptblgrow_maxd, ntot_amode)

  character(len=fieldname_len) :: dgnum_name(ntot_amode)

  ! Namelist variables
  logical :: sscav_tuning, convproc_do_aer, resus_fix  
  character(len=16) :: wetdep_list(pcnst) = ' '
  character(len=16) :: drydep_list(pcnst) = ' '
  real(r8)          :: sol_facti_cloud_borne = 1._r8
  real(r8)          :: sol_factb_interstitial  = 0.1_r8
  real(r8)          :: sol_factic_interstitial = 0.4_r8
  real(r8)          :: seasalt_emis_scale

  integer :: ndrydep = 0
  integer,allocatable :: drydep_indices(:)
  integer :: nwetdep = 0
  integer,allocatable :: wetdep_indices(:)
  logical :: drydep_lq(pcnst)
  logical :: wetdep_lq(pcnst)


contains
  
  !=============================================================================
  ! reads aerosol namelist options
  !=============================================================================
  subroutine aero_model_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'aero_model_readnl'

    ! Namelist variables
    character(len=16) :: aer_wetdep_list(pcnst) = ' '
    character(len=16) :: aer_drydep_list(pcnst) = ' '

    namelist /aerosol_nl/ aer_wetdep_list, aer_drydep_list, sol_facti_cloud_borne, seasalt_emis_scale,&
       sscav_tuning,sol_factb_interstitial, sol_factic_interstitial
    ! FIXME: many of these namelists are not used, can be cleaned up.
    !-----------------------------------------------------------------------------

    ! Read namelist
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'aerosol_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, aerosol_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          endif
       endif
       close(unitn)
       call freeunit(unitn)

    endif

#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast(aer_wetdep_list,   len(aer_wetdep_list(1))*pcnst, mpichar, 0, mpicom)
    call mpibcast(aer_drydep_list,   len(aer_drydep_list(1))*pcnst, mpichar, 0, mpicom)
    call mpibcast(sol_facti_cloud_borne, 1,                         mpir8,   0, mpicom)
    call mpibcast(sol_factb_interstitial, 1,                        mpir8,   0, mpicom)
    call mpibcast(sol_factic_interstitial, 1,                       mpir8,   0, mpicom)
    call mpibcast(sscav_tuning,          1,                         mpilog,  0, mpicom)
    call mpibcast(seasalt_emis_scale, 1, mpir8,   0, mpicom)
#endif

    if (sscav_tuning == .false.) then
       call endrun('ERROR: sscav_tuning==.false. option is removed in MAM4xx')
    endif

    wetdep_list = aer_wetdep_list
    drydep_list = aer_drydep_list

  end subroutine aero_model_readnl

  !=============================================================================
  !=============================================================================
  subroutine aero_model_register(imozart_in, species_class)
    use modal_aero_initialize_data, only : modal_aero_register
    integer, intent(in) :: imozart_in
    integer, intent(inout) :: species_class(:) 

    imozart = imozart_in
    call modal_aero_register(species_class)

  end subroutine aero_model_register

  !=============================================================================
  !=============================================================================
  subroutine aero_model_init( pbuf2d, species_class, iflagaa )

    use mo_chem_utls,    only: get_inv_ndx
    use cam_history,     only: addfld, horiz_only, add_default
    use mo_chem_utls,    only: get_rxt_ndx, get_spc_ndx
    use modal_aero_initialize_data, only: modal_aero_initialize
    use rad_constituents,           only: rad_cnst_get_info
    use dust_model,      only: dust_init, dust_names, dust_active, dust_nbin, dust_nnum
    use seasalt_model,   only: seasalt_init, seasalt_names, seasalt_active,seasalt_nbin
    use drydep_mod,      only: inidrydep
    use wetdep,          only: wetdep_init
    use mo_chem_utls,    only: get_het_ndx
    use gas_wetdep_opts, only: gas_wetdep_cnt, gas_wetdep_list, gas_wetdep_method ! REASTER 08/04/2015

    ! args
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    integer, intent(inout) :: species_class(:)  
    integer, intent(in) :: iflagaa

    ! local vars
    integer :: id, l, m, n, nspc

    logical  :: history_aerosol ! Output MAM or SECT aerosol tendencies
    logical  :: history_verbose ! produce verbose history output

    character(len=*), parameter :: subrname = 'aero_model_init'
    character(len=20) :: dummy
    character(len=fieldname_len) :: wetdep_name, depflx_name
    character(len=6) :: test_name
    character(len=100) :: errmes
    character(len=2)  :: unit_basename  ! Units 'kg' or '1' 

    if ( masterproc ) write(iulog,'(a,i5)') 'aero_model_init iflagaa=', iflagaa ! REASTER 08/04/2015

    call phys_getopts( history_aerosol_out=history_aerosol, &
         history_verbose_out=history_verbose, &
         convproc_do_aer_out = convproc_do_aer, & 
         resus_fix_out       = resus_fix    ) 

    ! This section cannot execute until chemini, ..., chm_diags_inti have been called
    if ( iflagaa == 2 ) then
       if ( masterproc ) then
          write(iulog,'(a,i5,2x,a)') 'gas_wetdep_cnt,meth', gas_wetdep_cnt, gas_wetdep_method
          do m = 1, gas_wetdep_cnt
          write(iulog,'(a,i5,2x,a)') 'gas_wetdep_list    ', m, trim(gas_wetdep_list(m))
          enddo
       endif

       ! These WD_ and DF_ fields should always been in a MAM history file, 
       do m = 1,gas_pcnst
          call cnst_get_ind( solsym(m), l, abrtf=.false. )
          if ( ( history_aerosol ) .and. (l > 0) ) then
             if  ( species_class(l) == spec_class_gas ) then !RCE - only output WD_xxx and DF_xxx for gases
                wetdep_name = 'WD_'//trim(solsym(m))
                depflx_name = 'DF_'//trim(solsym(m)) 
                nspc = get_het_ndx(solsym(m)) 
                if (nspc > 0) call add_default( wetdep_name, 1, ' ' )
                call add_default( depflx_name, 1, ' ' )
             endif
          endif
       enddo ! m = 1,gas_pcnst
       return
    endif ! ( iflagaa == 2 )


    ! The unified convective transport/removal for aerosols does not 
    ! do gases yet, and convproc_do_gas is just a place holder.  For that reason, 
    !    (1) All of the "if ( convproc_do_aer .or. convproc_do_gas ) then" statements 
    !        in aero_model.F90 have been changed to "if ( convproc_do_aer ) then"
    !    (2) convproc_do_aer=.false. and convproc_do_gas=.true. is no longer allowed.
    ! for C++ porting: All convproc_do_aer=.false. and resus_fix=.false. conditions
    ! are removed as not been used/tested. These conditions are then not allowed.
    ! convproc_do_gas is not used
    if ( ( .not. convproc_do_aer ) .or. ( .not. resus_fix ) ) then
       errmes = 'aero_model_init - convproc_do_aer and resus_fix MUST BE .true.' 
       call endrun( errmes )
    endif

    dgnum_idx      = pbuf_get_index('DGNUM')
    dgnumwet_idx   = pbuf_get_index('DGNUMWET')
    
    call rad_cnst_get_info(0, nmodes=nmodes)

    call modal_aero_initialize(pbuf2d, imozart, species_class) 
    call modal_aero_bcscavcoef_init()
    call mam_prevap_resusp_init( ) ! REASTER 08/04/2015

    call dust_init()
    call seasalt_init()
    call wetdep_init()

    fracis_idx      = pbuf_get_index('FRACIS') 
    prain_idx       = pbuf_get_index('PRAIN')  
    nevapr_idx      = pbuf_get_index('NEVAPR') 
    rprddp_idx      = pbuf_get_index('RPRDDP')  
    rprdsh_idx      = pbuf_get_index('RPRDSH')  
    
    nevapr_shcu_idx = pbuf_get_index('NEVAPR_SHCU')
    nevapr_dpcu_idx = pbuf_get_index('NEVAPR_DPCU')

    icwmrdp_idx      = pbuf_get_index('ICWMRDP')
    icwmrsh_idx      = pbuf_get_index('ICWMRSH')
    sh_frac_idx      = pbuf_get_index('SH_FRAC')
    dp_frac_idx      = pbuf_get_index('DP_FRAC')


    nwetdep = 0
    ndrydep = 0

    count_species: do m = 1,pcnst
       if ( len_trim(wetdep_list(m)) /= 0 ) then
          nwetdep = nwetdep+1
       endif
       if ( len_trim(drydep_list(m)) /= 0 ) then
          ndrydep = ndrydep+1
       endif
    enddo count_species
    
    if (nwetdep>0) &
         allocate(wetdep_indices(nwetdep))
    if (ndrydep>0) &
         allocate(drydep_indices(ndrydep))

    do m = 1,ndrydep
       call cnst_get_ind ( drydep_list(m), id, abrtf=.false. )
       if (id>0) then
          drydep_indices(m) = id
       else
          call endrun(subrname//': invalid drydep species: '//trim(drydep_list(m)) )
       endif

       if (masterproc) then
          write(iulog,*) subrname//': '//drydep_list(m)//' will have drydep applied'
       endif
    enddo
    do m = 1,nwetdep
       call cnst_get_ind ( wetdep_list(m), id, abrtf=.false. )
       if (id>0) then
          wetdep_indices(m) = id
       else
          call endrun(subrname//': invalid wetdep species: '//trim(wetdep_list(m)) )
       endif
       
       if (masterproc) then
          write(iulog,*) subrname//': '//wetdep_list(m)//' will have wet removal'
       endif
    enddo

    if (ndrydep>0) then

       call inidrydep(rair, gravit)

       dummy = 'RAM1'
       call addfld (dummy,horiz_only, 'A','frac','RAM1')
       if ( history_aerosol ) then  
          call add_default (dummy, 1, ' ')
       endif
       dummy = 'airFV'
       call addfld (dummy,horiz_only, 'A','frac','FV')
       if ( history_aerosol ) then  
          call add_default (dummy, 1, ' ')
       endif

    endif

    if (dust_active) then
       ! emissions diagnostics ....

       do m = 1, dust_nbin+dust_nnum
          dummy = trim(dust_names(m)) // 'SF'
          call addfld (dummy,horiz_only, 'A','kg/m2/s',trim(dust_names(m))//' dust surface emission')
          if (history_aerosol) then
             call add_default (dummy, 1, ' ')
          endif
       enddo

       dummy = 'DSTSFMBL'
       call addfld (dummy,horiz_only, 'A','kg/m2/s','Mobilization flux at surface')
       if (history_aerosol) then
          call add_default (dummy, 1, ' ')
       endif

       dummy = 'LND_MBL'
       call addfld (dummy,horiz_only, 'A','1','Soil erodibility factor')
       if (history_aerosol) then
          call add_default (dummy, 1, ' ')
       endif

    endif

    if (seasalt_active) then
       
       dummy = 'SSTSFMBL'
       call addfld (dummy,horiz_only, 'A','kg/m2/s','Mobilization flux at surface')
       if (history_aerosol) then
          call add_default (dummy, 1, ' ')
       endif

       do m = 1, seasalt_nbin
          dummy = trim(seasalt_names(m)) // 'SF'
          call addfld (dummy,horiz_only, 'A','kg/m2/s',trim(seasalt_names(m))//' seasalt surface emission')
          if (history_aerosol) then
             call add_default (dummy, 1, ' ')
          endif
       enddo

#if (defined MODAL_AERO_9MODE || MODAL_AERO_4MODE_MOM)
       dummy = 'SSTSFMBL_OM'
       call addfld (dummy,horiz_only, 'A','kg/m2/s','Mobilization flux of marine organic matter at surface')
       if (history_aerosol) then
          call add_default (dummy, 1, ' ')
       endif

       dummy = 'F_eff'
       call addfld (dummy,horiz_only, 'A','1','Effective enrichment factor of marine organic matter')
       if (history_aerosol) then
          call add_default (dummy, 1, ' ')
       endif
#endif

    endif

    
    ! set flags for drydep tendencies
    drydep_lq(:) = .false.
    do m=1,ndrydep 
       id = drydep_indices(m)
       drydep_lq(id) =  .true.
    enddo

    ! set flags for wetdep tendencies
    wetdep_lq(:) = .false.
    do m=1,nwetdep
       id = wetdep_indices(m)
       wetdep_lq(id) = .true.
    enddo

    wetdens_ap_idx = pbuf_get_index('WETDENS_AP')
    qaerwat_idx    = pbuf_get_index('QAERWAT')
    pblh_idx       = pbuf_get_index('pblh')

    rate1_cw2pr_st_idx  = pbuf_get_index('RATE1_CW2PR_ST') 
    call pbuf_set_field(pbuf2d, rate1_cw2pr_st_idx, 0.0_r8)

    do m = 1,ndrydep
       
       ! units 
       if (drydep_list(m)(1:3) == 'num') then
          unit_basename = ' 1'
       else
          unit_basename = 'kg'  
       endif

       call addfld (trim(drydep_list(m))//'DDF',   horiz_only, 'A',unit_basename//'/m2/s ', &
            trim(drydep_list(m))//' dry deposition flux at bottom (grav + turb)')
       call addfld (trim(drydep_list(m))//'TBF',   horiz_only, 'A',unit_basename//'/m2/s', &
            trim(drydep_list(m))//' turbulent dry deposition flux')
       call addfld (trim(drydep_list(m))//'GVF',   horiz_only, 'A',unit_basename//'/m2/s ', &
            trim(drydep_list(m))//' gravitational dry deposition flux')
       call addfld (trim(drydep_list(m))//'DTQ',(/ 'lev' /), 'A',unit_basename//'/kg/s ', &
            trim(drydep_list(m))//' dry deposition')
       call addfld (trim(drydep_list(m))//'DDV',(/ 'lev' /), 'A','m/s', &
            trim(drydep_list(m))//' deposition velocity')

       if ( history_aerosol ) then 
          call add_default (trim(drydep_list(m))//'DDF', 1, ' ')
          if ( history_verbose ) then
             call add_default (trim(drydep_list(m))//'TBF', 1, ' ')
             call add_default (trim(drydep_list(m))//'GVF', 1, ' ')
          endif
       endif

    enddo

    do m = 1,nwetdep
       if ( masterproc ) write(iulog,'(a,i3,2x,a)') 'm, wetdep_list', m, trim(wetdep_list(m)) ! REASTER 08/04/2015
       
       ! units 
       if (wetdep_list(m)(1:3) == 'num') then
          unit_basename = ' 1'
       else
          unit_basename = 'kg'  
       endif

       call addfld (trim(wetdep_list(m))//'SFWET', &
            horiz_only,  'A',unit_basename//'/m2/s ','Wet deposition flux at surface')
       call addfld (trim(wetdep_list(m))//'SFSIC', &
            horiz_only,  'A',unit_basename//'/m2/s ','Wet deposition flux (incloud, convective) at surface')
       call addfld (trim(wetdep_list(m))//'SFSIS', &
            horiz_only,  'A',unit_basename//'/m2/s ','Wet deposition flux (incloud, stratiform) at surface')
       call addfld (trim(wetdep_list(m))//'SFSBC', &
            horiz_only,  'A',unit_basename//'/m2/s ','Wet deposition flux (belowcloud, convective) at surface')
       call addfld (trim(wetdep_list(m))//'SFSBS', &
            horiz_only,  'A',unit_basename//'/m2/s ','Wet deposition flux (belowcloud, stratiform) at surface')
       call addfld (trim(wetdep_list(m))//'SFSEC', &
            horiz_only,  'A','kg/m2/s','Wet deposition flux (precip evap, convective) at surface')
       call addfld (trim(wetdep_list(m))//'SFSES', &
            horiz_only,  'A','kg/m2/s','Wet deposition flux (precip evap, stratiform) at surface')
       call addfld (trim(wetdep_list(m))//'SFSED', &
            horiz_only,  'A','kg/m2/s','Wet deposition flux (precip evap, deep convective) at surface')
       call addfld (trim(wetdep_list(m))//'SFSID', &
            horiz_only,  'A','kg/m2/s','Wet deposition flux (incloud, deep convective) at surface')
       call addfld (trim(wetdep_list(m))//'SFSBD', &
            horiz_only,  'A','kg/m2/s','Wet deposition flux (belowcloud, deep convective) at surface')

       call addfld (trim(wetdep_list(m))//'WET',(/ 'lev' /), 'A',unit_basename//'/kg/s ','wet deposition tendency')
       call addfld (trim(wetdep_list(m))//'SIC',(/ 'lev' /), 'A',unit_basename//'/kg/s ', &
            trim(wetdep_list(m))//' ic wet deposition')
       call addfld (trim(wetdep_list(m))//'SIS',(/ 'lev' /), 'A',unit_basename//'/kg/s ', &
            trim(wetdep_list(m))//' is wet deposition')
       call addfld (trim(wetdep_list(m))//'SBC',(/ 'lev' /), 'A',unit_basename//'/kg/s ', &
            trim(wetdep_list(m))//' bc wet deposition')
       call addfld (trim(wetdep_list(m))//'SBS',(/ 'lev' /), 'A',unit_basename//'/kg/s ', &
            trim(wetdep_list(m))//' bs wet deposition')
       
       if ( history_aerosol ) then          
          call add_default (trim(wetdep_list(m))//'SFWET', 1, ' ')
          if ( history_verbose ) then
             call add_default (trim(wetdep_list(m))//'SFSIC', 1, ' ')
             call add_default (trim(wetdep_list(m))//'SFSIS', 1, ' ')
             call add_default (trim(wetdep_list(m))//'SFSBC', 1, ' ')
             call add_default (trim(wetdep_list(m))//'SFSBS', 1, ' ')
             call add_default (trim(wetdep_list(m))//'SFSEC', 1, ' ')
             call add_default (trim(wetdep_list(m))//'SFSES', 1, ' ')
          endif
       endif

    enddo ! m = 1,nwetdep

    do m = 1,gas_pcnst

       if  ( solsym(m)(1:3) == 'num') then
          unit_basename = ' 1'  ! Units 'kg' or '1' 
       else
          unit_basename = 'kg'  ! Units 'kg' or '1' 
       end if

       call addfld( 'GS_'//trim(solsym(m)),horiz_only,  'A', unit_basename//'/m2/s ', &
                    trim(solsym(m))//' gas chemistry/wet removal (for gas species)')
       call addfld( 'AQ_'//trim(solsym(m)),horiz_only,  'A', unit_basename//'/m2/s ', &
                    trim(solsym(m))//' aqueous chemistry (for gas species)')
       if ( history_aerosol ) then 
          if ( history_verbose ) then
             call add_default( 'GS_'//trim(solsym(m)), 1, ' ')
             call add_default( 'AQ_'//trim(solsym(m)), 1, ' ')
          else
             select case (trim(solsym(m)))
             case ('O3','H2O2','H2SO4','SO2','DMS','SOAG')
                  call add_default( 'AQ_'//trim(solsym(m)), 1, ' ')
             end select
          end if
       endif
       
       call cnst_get_ind(trim(solsym(m)), nspc, abrtf=.false. ) ! REASTER 08/04/2015
       if( nspc > 0 ) then                                      ! REASTER 08/04/2015
        if ( .not. cnst_name_cw(nspc) == ' ') then              ! REASTER 08/04/2015
          call addfld (trim(cnst_name_cw(nspc))//'SFSEC',horiz_only,  'A','kg/m2/s', &
               trim(cnst_name_cw(nspc))//' wet deposition flux (precip evap, convective) at surface')  !RCE
          call addfld (trim(cnst_name_cw(nspc))//'SFSES',horiz_only,  'A','kg/m2/s', &
               trim(cnst_name_cw(nspc))//' wet deposition flux (precip evap, stratiform) at surface')  !RCE             
          if(history_aerosol .and. history_verbose) then
             call add_default (trim(cnst_name_cw(nspc))//'SFSEC', 1, ' ')  !RCE
             call add_default (trim(cnst_name_cw(nspc))//'SFSES', 1, ' ')  !RCE
          endif
        endif
       endif

    enddo

    do n = 1,pcnst
       if( .not. (cnst_name_cw(n) == ' ') ) then

          if (cnst_name_cw(n)(1:3) == 'num') then
             unit_basename = ' 1'
          else
             unit_basename = 'kg'  
          endif

          call addfld( cnst_name_cw(n), (/ 'lev' /), 'A',                unit_basename//'/kg ', &
               trim(cnst_name_cw(n))//' in cloud water')
          call addfld (trim(cnst_name_cw(n))//'SFWET',horiz_only,  'A', unit_basename//'/m2/s ', &
               trim(cnst_name_cw(n))//' wet deposition flux at surface')
          call addfld (trim(cnst_name_cw(n))//'SFSIC',horiz_only,  'A', unit_basename//'/m2/s ', &
               trim(cnst_name_cw(n))//' wet deposition flux (incloud, convective) at surface')
          call addfld (trim(cnst_name_cw(n))//'SFSIS',horiz_only,  'A', unit_basename//'/m2/s ', &
               trim(cnst_name_cw(n))//' wet deposition flux (incloud, stratiform) at surface')
          call addfld (trim(cnst_name_cw(n))//'SFSBC',horiz_only,  'A', unit_basename//'/m2/s ', &
               trim(cnst_name_cw(n))//' wet deposition flux (belowcloud, convective) at surface')
          call addfld (trim(cnst_name_cw(n))//'SFSBS',horiz_only,  'A', unit_basename//'/m2/s ', &
               trim(cnst_name_cw(n))//' wet deposition flux (belowcloud, stratiform) at surface')
          call addfld (trim(cnst_name_cw(n))//'DDF',   horiz_only, 'A',   unit_basename//'/m2/s ', &
               trim(cnst_name_cw(n))//' dry deposition flux at bottom (grav + turb)')
          call addfld (trim(cnst_name_cw(n))//'TBF',   horiz_only, 'A',   unit_basename//'/m2/s ', &
               trim(cnst_name_cw(n))//' turbulent dry deposition flux')
          call addfld (trim(cnst_name_cw(n))//'GVF',   horiz_only, 'A',   unit_basename//'/m2/s ', &
               trim(cnst_name_cw(n))//' gravitational dry deposition flux')     

          if ( history_aerosol ) then 
             if (history_verbose) then
                call add_default( cnst_name_cw(n), 1, ' ' )
                call add_default (trim(cnst_name_cw(n))//'GVF', 1, ' ')
                call add_default (trim(cnst_name_cw(n))//'TBF', 1, ' ')
                call add_default (trim(cnst_name_cw(n))//'SFSBS', 1, ' ')      
                call add_default (trim(cnst_name_cw(n))//'SFSIC', 1, ' ')
                call add_default (trim(cnst_name_cw(n))//'SFSBC', 1, ' ')
                call add_default (trim(cnst_name_cw(n))//'SFSIS', 1, ' ')
             endif
             call add_default (trim(cnst_name_cw(n))//'SFWET', 1, ' ') 
             call add_default (trim(cnst_name_cw(n))//'DDF', 1, ' ')
          endif
       endif
    enddo

    do n=1,ntot_amode
       dgnum_name(n) = ' '
       write(dgnum_name(n),fmt='(a,i1)') 'dgnumwet',n
       call addfld( dgnum_name(n), (/ 'lev' /), 'I', 'm', 'Aerosol mode wet diameter' )
       if ( history_aerosol .and. history_verbose ) then 
          call add_default( dgnum_name(n), 1, ' ' )
       endif
    enddo

  end subroutine aero_model_init

  !=============================================================================
  !=============================================================================
    subroutine mam_prevap_resusp_init( )

    use modal_aero_data, only: &
       lmassptr_amode, lspectype_amode, modeptr_coarse, &
       nspec_amode, ntot_amode, numptr_amode, &
       mmtoo_prevap_resusp, ntoo_prevap_resusp

    integer :: lspec, lspec2
    integer :: mm, mmtoo, mm2
    integer :: n, ntoo, nch
    character(len=100) :: msg


! calculate pointers for resuspension
! mmtoo_prevap_resusp values are
!    >0 for aerosol mass species with    coarse mode counterpart
!    -1 for aerosol mass species WITHOUT coarse mode counterpart
!    -2 for aerosol number species
!     0 for other species

    mmtoo_prevap_resusp(:) = 0
    ntoo_prevap_resusp(:) = 0


       ntoo = modeptr_coarse
       do n = 1, ntot_amode   ! loop over aerosol modes that was wet-removed

          do lspec = 1, nspec_amode(n)   ! loop over chem constituents that was wet-removed
             mm = lmassptr_amode(lspec,n)  ! q-array index of the species that was wet-removed
             nch = len( trim( cnst_name(mm) ) ) - 1

             mmtoo = -1   ! q-array index of the coarse mode species that gets the resuspension
             do lspec2 = 1, nspec_amode(ntoo)
!               match based on the cnst_name (except for the last 1-2 characters)
                mm2 = lmassptr_amode(lspec2,ntoo)
                if ( cnst_name(mm)(1:nch) == cnst_name(mm2)(1:nch) ) then
                   mmtoo = mm2
                   exit
                endif
             enddo

             if ( masterproc ) then
                write(iulog,'(a,3(2x,a))') 'modal_aero_wetscav_init mmfrm/too:  ', &
                   cnst_name(mm), cnst_name(mmtoo), cnst_name(numptr_amode(ntoo))
             endif

             mmtoo_prevap_resusp(mm) = mmtoo
             ntoo_prevap_resusp(mm)  = ntoo
          enddo ! lspec

          mm = numptr_amode(n)
          mmtoo_prevap_resusp(mm) = -2
          ntoo_prevap_resusp(mm)  = ntoo
       enddo ! n


    if ( masterproc ) then
       do mm = 1, pcnst
          mmtoo = mmtoo_prevap_resusp(mm)
          ntoo = ntoo_prevap_resusp(mm)
          msg = ' '
          if (mmtoo > 0) msg = cnst_name(mmtoo)
          write(iulog,'(2a,3(1x,i9),2x,a)') 'name, mm, mmtoo, ntoo =  ', &
             cnst_name(mm), mm, mmtoo, ntoo, trim(msg)
       enddo
    endif

    end subroutine mam_prevap_resusp_init


  !=============================================================================
  !=============================================================================
  subroutine aero_model_wetdep(dt, dlf, dlf2, cmfmc2, state,                    & ! in
       sh_e_ed_ratio, mu, md, du, eu, ed, dp, jt, maxg, ideep, lengath,         & ! in
       species_class,                                                           & ! in
       cam_out,                                                                 & ! inout
       pbuf,                                                                    & ! Pointer
       ptend                                                                    ) ! out

    use modal_aero_deposition, only: set_srf_wetdep
    use wetdep,                only: wetdepa_v2, wetdep_inputs_set, &
                                     wetdep_inputs_unset, wetdep_inputs_t
    use modal_aero_data
    use modal_aero_calcsize,   only: modal_aero_calcsize_sub
    use modal_aero_wateruptake,only: modal_aero_wateruptake_dr
    use modal_aero_convproc,   only: ma_convproc_intr

    ! args
    type(physics_state), intent(in)    :: state       ! Physics state variables
    real(r8),            intent(in)    :: dt          ! time step
    real(r8),            intent(in)    :: dlf(:,:)    ! shallow+deep convective detrainment [kg/kg/s]
    real(r8),            intent(in)    :: dlf2(:,:)   ! Shal conv cldwtr detrainment (kg/kg/s - grid avg)
    real(r8),            intent(in)    :: cmfmc2(pcols,pverp) ! Shal conv mass flux (kg/m2/s)
    real(r8),            intent(in)    :: sh_e_ed_ratio(pcols,pver)  ! shallow conv [ent/(ent+det)] ratio
    ! mu, md, ..., ideep, lengath are all deep conv variables
    ! *** AND ARE GATHERED ***
    ! eu, ed, du are "d(massflux)/dp" and are all positive
    real(r8),            intent(in)    :: mu(pcols,pver)   ! Updraft mass flux (positive)
    real(r8),            intent(in)    :: md(pcols,pver)   ! Downdraft mass flux (negative)
    real(r8),            intent(in)    :: du(pcols,pver)   ! Mass detrain rate from updraft
    real(r8),            intent(in)    :: eu(pcols,pver)   ! Mass entrain rate into updraft
    real(r8),            intent(in)    :: ed(pcols,pver)   ! Mass entrain rate into downdraft
    real(r8),            intent(in)    :: dp(pcols,pver)   ! Delta pressure between interfaces
    
    integer,             intent(in)    :: jt(pcols)         ! Index of cloud top for each column
    integer,             intent(in)    :: maxg(pcols)       ! Index of cloud top for each column
    integer,             intent(in)    :: ideep(pcols)      ! Gathering array
    integer,             intent(in)    :: lengath           ! Gathered min lon indices over which to operate
    integer,             intent(in)    :: species_class(:)
    
    type(cam_out_t),     intent(inout) :: cam_out     ! export state
    type(physics_buffer_desc), pointer :: pbuf(:)

    type(physics_ptend), intent(out)   :: ptend       ! indivdual parameterization tendencies



    ! local vars

    integer :: jnv ! index for scavcoefnv 3rd dimension
    integer :: jnummaswtr  ! indicates current aerosol species type (0 = number, 1 = dry mass, 2 = water)
    integer :: imode
    integer :: lchnk ! chunk identifier
    integer :: lphase ! index for interstitial / cloudborne aerosol
    integer :: lspec ! index for aerosol number / chem-mass / water-mass
    integer :: mtmp ! mode index
    integer :: mm ! tracer (q-array) index
    integer :: ncol ! number of atmospheric columns
    integer :: mam_prevap_resusp_optcc

    real(r8) :: iscavt(pcols, pver)
    real(r8) :: icscavt(pcols, pver)
    real(r8) :: isscavt(pcols, pver)
    real(r8) :: bcscavt(pcols, pver)
    real(r8) :: bsscavt(pcols, pver)
    real(r8) :: sol_factb, sol_facti
    real(r8) :: sol_factic(pcols,pver)

    real(r8) :: sflx(pcols) ! deposition flux

    real(r8) :: dqdt_tmp(pcols,pver)      ! temporary array to hold tendency for the "current" aerosol species
    real(r8) :: f_act_conv(pcols,pver) ! prescribed aerosol activation fraction for convective cloud 
    real(r8) :: f_act_conv_coarse(pcols,pver) ! similar but for coarse mode 
    real(r8) :: f_act_conv_coarse_dust, f_act_conv_coarse_nacl
    real(r8) :: fracis_cw(pcols,pver)
    real(r8) :: prec(pcols) ! precipitation rate
    real(r8) :: q_tmp(pcols,pver) ! temporary array to hold "most current" mixing ratio for 1 species
    real(r8) :: qqcw_tmp(pcols,pver), qqcw_in(pcols,pver) ! temporary array to hold qqcw 
    real(r8) :: qqcw_sav(pcols,pver,0:maxd_aspectype) ! temporary array to hold qqcw for the current mode  !RCE
    real(r8) :: scavcoefnv(pcols,pver,0:2) ! Dana and Hales coefficient (/mm) for
                                           ! cloud-borne num & vol (0),
                                           ! interstitial num (1), interstitial vol (2)

    logical  :: isprx(pcols,pver) ! true if precipation

    real(r8) :: sflxec(pcols), sflxecdp(pcols)  ! deposition flux  !RCE
    real(r8) :: sflxic(pcols), sflxicdp(pcols)  ! deposition flux  !RCE
    real(r8) :: sflxbc(pcols), sflxbcdp(pcols)  ! deposition flux  !RCE
    real(r8) :: rcscavt(pcols, pver)  !RCE
    real(r8) :: rsscavt(pcols, pver)  !RCE
    real(r8) :: rtscavt_sv(pcols, pver, pcnst) ! REASTER 08/12/2015
    
    real(r8), pointer :: fldcw(:,:)
    real(r8), pointer :: dgnumwet(:,:,:)

    real(r8), pointer :: fracis(:,:,:)   ! fraction of transported species that are insoluble

    integer, parameter:: nsrflx_mzaer2cnvpr = 2  !RCE 2012/01/12 bgn
    real(r8)          :: aerdepwetis(pcols,pcnst) ! aerosol wet deposition (interstitial) 
    real(r8)          :: aerdepwetcw(pcols,pcnst) ! aerosol wet deposition (cloud water)  
    real(r8)          :: qsrflx_mzaer2cnvpr(pcols,pcnst,nsrflx_mzaer2cnvpr)
    real(r8)          :: rprddpsum(pcols),  rprdshsum(pcols)  
    real(r8)          :: evapcdpsum(pcols), evapcshsum(pcols) 
    real(r8), pointer :: rprddp(:,:)     ! rain production, deep convection
    real(r8), pointer :: rprdsh(:,:)     ! rain production, deep convection
    real(r8), pointer :: evapcsh(:,:)    ! Evaporation rate of shallow convective precipitation >=0.
    real(r8), pointer :: evapcdp(:,:)    ! Evaporation rate of deep    convective precipitation >=0.

    real(r8), pointer :: icwmrdp(:,:)    ! in cloud water mixing ratio, deep convection
    real(r8), pointer :: icwmrsh(:,:)    ! in cloud water mixing ratio, deep convection
    real(r8), pointer :: sh_frac(:,:)    ! Shallow convective cloud fraction
    real(r8), pointer :: dp_frac(:,:)    ! Deep convective cloud fraction

    type(wetdep_inputs_t) :: dep_inputs


    lchnk = state%lchnk
    ncol  = state%ncol

    call physics_ptend_init(ptend, state%psetcols, 'aero_model_wetdep_ma', lq=wetdep_lq)

    ! Do calculations of mode radius and water uptake if:
    ! 1) modal aerosols are affecting the climate, or
    ! 2) prognostic modal aerosols are enabled
    ! If not using prognostic aerosol call the diagnostic version

    ! Calculate aerosol size distribution parameters
    ! for prognostic modal aerosols the transfer of mass between aitken and 
    ! accumulation modes is done in conjunction with the dry radius calculation
    call t_startf('calcsize')
    call modal_aero_calcsize_sub(state, dt, pbuf, ptend)
    call t_stopf('calcsize')
    
    ! Aerosol water uptake
    call t_startf('wateruptake')
    call modal_aero_wateruptake_dr(state, pbuf)
    call t_stopf('wateruptake')

    if (nwetdep<1) return

    call wetdep_inputs_set( state, pbuf, dep_inputs )

    call pbuf_get_field(pbuf, dgnumwet_idx, dgnumwet, start=(/1,1,1/), kount=(/pcols,pver,nmodes/) )
    call pbuf_get_field(pbuf, fracis_idx,   fracis,   start=(/1,1,1/), kount=(/pcols,pver, pcnst/) )

    !Compute variables needed for convproc unified convective transport
    call pbuf_get_field(pbuf, rprddp_idx,      rprddp  )
    call pbuf_get_field(pbuf, rprdsh_idx,      rprdsh  )
    call pbuf_get_field(pbuf, nevapr_shcu_idx, evapcsh )
    call pbuf_get_field(pbuf, nevapr_dpcu_idx, evapcdp )

    call calc_sfc_flux(rprdsh,  state%pdel, rprdshsum)
    call calc_sfc_flux(rprddp,  state%pdel, rprddpsum)
    call calc_sfc_flux(evapcsh, state%pdel, evapcshsum)
    call calc_sfc_flux(evapcdp, state%pdel, evapcdpsum)

    ! initiate variables
    qsrflx_mzaer2cnvpr(:,:,:) = 0.0_r8  
    aerdepwetis(:,:)          = 0.0_r8  
    aerdepwetcw(:,:)          = 0.0_r8  
    qqcw_tmp(:,:)             = 0.0_r8  
    ! below-cloud scavcoef = 0.0 for cloud-borne species
    scavcoefnv(:,:,0)         = 0.0_r8
    ! resuspension goes to a different phase or mode
    rtscavt_sv(:,:,:)         = 0.0_r8

    ! examine if there is precipitation falling from above in each grid
    call examine_prec_exist ( ncol,  state%pdel,      & ! in
                 dep_inputs%prain,  dep_inputs%cmfdqr,& ! in
                 dep_inputs%evapr,                    & ! in
                 isprx                                ) ! out

    ! calculate the mass-weighted sol_factic for coarse mode species
    call set_f_act_coarse(      ncol,                           & ! in
                state%q,        ptend%q,        dt,             & ! in
                f_act_conv_coarse, f_act_conv_coarse_dust,      & ! out
                f_act_conv_coarse_nacl                          ) ! out

mmode_loop_aa: &
    do mtmp = 1, ntot_amode ! main loop over aerosol modes
       imode = mtmp
       ! for mam4, do accum, aitken, pcarbon, then coarse 
       ! so change the order of 3 and 4 here
       if (mtmp == modeptr_coarse) then
             imode = ntot_amode
       elseif (mtmp > modeptr_coarse) then
             imode = mtmp - 1
       endif

! loop over interstitial (1) and cloud-borne (2) forms         
!BSINGH (09/12/2014):Do cloudborne first for unified convection scheme so
!that the resuspension of cloudborne can be saved then applied to interstitial (RCE) 
lphase_loop_aa: &
       do lphase = 2,1,-1  ! do cloudborne (2) first then interstitial (1)
          if (lphase == 1) then ! interstial aerosol
             call modal_aero_bcscavcoef_get( imode, ncol, isprx, dgnumwet, &
                  scavcoefnv(:,:,1), scavcoefnv(:,:,2) )
          endif

          call define_act_frac ( lphase,     imode,         & ! in
                sol_facti, sol_factic, sol_factb, f_act_conv) ! out

! REASTER 08/12/2015 - changed ordering (mass then number) for prevap resuspend to coarse
lspec_loop_aa: &
          do lspec = 1, nspec_amode(imode)+2 ! loop over number + chem constituents + water

             call index_ordering (                 &
                        lspec, imode,  lphase,     & ! in
                        mm,    jnv, jnummaswtr     ) ! out

             if (mm <= 0 .or. jnummaswtr == 2) cycle  ! by pass wet aerosols

! mam_prevap_resusp_optcc values control the prevap_resusp calculations in wetdepa_v2:
!     0 = no resuspension
!   130 = non-linear resuspension of aerosol mass   based on scavenged aerosol mass
!   230 = non-linear resuspension of aerosol number based on raindrop number
!   the 130 thru 230 all use the new prevap_resusp code block in subr wetdepa_v2
!
             mam_prevap_resusp_optcc = 0

             if ( jnummaswtr == 1 ) then  ! dry mass
                   mam_prevap_resusp_optcc = 130
             elseif ( jnummaswtr == 0 .and. lphase == 1 .and. imode == modeptr_coarse ) then ! number
                   mam_prevap_resusp_optcc = 230
             endif

             ! set f_act_conv for interstitial (lphase=1) coarse mode species
             ! for the convective in-cloud, we conceptually treat the coarse dust and seasalt
             ! as being externally mixed, and apply 
             ! f_act_conv = f_act_conv_coarse_dust/nacl to dust/seasalt
             ! number and sulfate are conceptually partitioned to the dust and seasalt
             ! on a mass basis, so the f_act_conv for number and sulfate are
             ! mass-weighted averages of the values used for dust/seasalt
             if ((lphase == 1) .and. (imode == modeptr_coarse)) then
                f_act_conv = f_act_conv_coarse 
                if (jnummaswtr == 1) then
                   if (lmassptr_amode(lspec,imode) == lptr_dust_a_amode(imode)) then
                      f_act_conv = f_act_conv_coarse_dust 
                   elseif (lmassptr_amode(lspec,imode) == lptr_nacl_a_amode(imode)) then
                      f_act_conv = f_act_conv_coarse_nacl 
                   endif
                endif
             endif

lphase_jnmw_conditional: &
             if (lphase == 1) then
                ptend%lq(mm) = .true.
                ! q_tmp reflects changes from modal_aero_calcsize and is the "most current" q
                q_tmp(1:ncol,:) = state%q(1:ncol,:,mm) + ptend%q(1:ncol,:,mm)*dt
                !Feed in the saved cloudborne mixing ratios from phase 2
                qqcw_in(:,:) = qqcw_sav(:,:,lspec)

                call wetdepa_v2( &
                     ncol, dt, state%pdel,                                      & ! in
                     dep_inputs%cmfdqr, dep_inputs%evapc, dlf, dep_inputs%conicw, & ! in
                     dep_inputs%prain, dep_inputs%evapr, dep_inputs%totcond,    & ! in
                     dep_inputs%cldt, dep_inputs%cldcu,                         & ! in
                     dep_inputs%cldvcu, dep_inputs%cldvst,                      & ! in
                     sol_factb, sol_facti, sol_factic,                          & ! in
                     mam_prevap_resusp_optcc, .false., scavcoefnv(:,:,jnv), f_act_conv, & ! in
                     q_tmp, qqcw_in(:,:),                                       & ! in
                     fracis(:,:,mm), dqdt_tmp, iscavt,                          & ! out
                     icscavt, isscavt, bcscavt, bsscavt, rcscavt, rsscavt       ) ! out

                ! resuspension goes to coarse mode
                call calc_resusp_to_coarse(     ncol,   mm,     & ! in
                        mmtoo_prevap_resusp,    .true.,         & ! in
                        rcscavt,        rsscavt,                & ! in
                        dqdt_tmp,       rtscavt_sv              ) ! inout

                ptend%q(1:ncol,:,mm) = ptend%q(1:ncol,:,mm) + dqdt_tmp(1:ncol,:)

                call outfld( trim(cnst_name(mm))//'WET', dqdt_tmp(:,:), pcols, lchnk)
                call outfld( trim(cnst_name(mm))//'SIC', icscavt, pcols, lchnk)
                call outfld( trim(cnst_name(mm))//'SIS', isscavt, pcols, lchnk)
                call outfld( trim(cnst_name(mm))//'SBC', bcscavt, pcols, lchnk)
                call outfld( trim(cnst_name(mm))//'SBS', bsscavt, pcols, lchnk)

                call calc_sfc_flux(dqdt_tmp, state%pdel, sflx)
                aerdepwetis(:ncol,mm) = sflx(:ncol)

                call calc_sfc_flux(icscavt, state%pdel, sflx)
                sflxic = sflx

                call calc_sfc_flux(isscavt, state%pdel, sflx)
                call outfld( trim(cnst_name(mm))//'SFSIS', sflx, pcols, lchnk)

                call calc_sfc_flux(bcscavt, state%pdel, sflx)
                call outfld( trim(cnst_name(mm))//'SFSBC', sflx, pcols, lchnk)
                sflxbc = sflx

                call calc_sfc_flux(bsscavt, state%pdel, sflx)
                call outfld( trim(cnst_name(mm))//'SFSBS', sflx, pcols, lchnk)
               
                ! here the prevap resuspension is in rcscavt & rsscavt and column integral is written to history
                !BSINGH(09/15/2014):Following two nested do-loops are new additions for unified convection 
                !BSINGH(09/15/2014):After these do-loops, code was added by RCE, the comments by RCE are kept as it is
                call calc_sfc_flux(rcscavt, state%pdel, sflx)
                sflxec = sflx
                   
                call calc_sfc_flux(rsscavt, state%pdel, sflx)
                call outfld( trim(cnst_name(mm))//'SFSES', sflx, pcols, lchnk)                   
                   
                ! apportion convective surface fluxes to deep and shallow conv
                ! this could be done more accurately in subr wetdepa
                ! since deep and shallow rarely occur simultaneously, and these
                !    fields are just diagnostics, this approximate method is adequate
                ! only do this for interstitial aerosol, because conv clouds to not
                !    affect the stratiform-cloudborne aerosol
                call apportion_sfc_flux_deep ( ncol,              & ! in
                        rprddpsum,rprdshsum,evapcdpsum,evapcshsum,& ! in
                        sflxbc,             sflxec,               & ! in
                        sflxbcdp,           sflxecdp              ) ! out

                call outfld( trim(cnst_name(mm))//'SFSBD', sflxbcdp, pcols, lchnk)
                ! when ma_convproc_intr is used, convective in-cloud wet removal is done there
                ! the convective (total and deep) precip-evap-resuspension includes in- and below-cloud
                ! contributions, so pass the below-cloud contribution to ma_convproc_intr
                qsrflx_mzaer2cnvpr(1:ncol,mm,1) = sflxec(  1:ncol)
                qsrflx_mzaer2cnvpr(1:ncol,mm,2) = sflxecdp(1:ncol)

             elseif (lphase == 2) then lphase_jnmw_conditional
! There is no cloud-borne aerosol water in the model, so this code block
! should NEVER execute for lspec = nspec_amode(m)+1 (i.e., jnummaswtr = 2).
! The code only worked because the "do lspec" loop cycles when lspec = nspec_amode(m)+1,
! but that does not make the code correct.
                fldcw => qqcw_get_field(pbuf,mm,lchnk)
                qqcw_sav(1:ncol,:,lspec) = fldcw(1:ncol,:)  !RCE 2012/01/12
            
                ! FIXME: Not sure if this is a bug or not as qqcw_tmp seem different
                ! from the previous call and qqcw_tmp is always zero. May need
                ! further check.  - Shuaiqi Tang in refactoring for MAM4xx
                call wetdepa_v2( &
                   ncol, dt, state%pdel,                                        & ! in
                   dep_inputs%cmfdqr, dep_inputs%evapc, dlf, dep_inputs%conicw, & ! in
                   dep_inputs%prain, dep_inputs%evapr, dep_inputs%totcond,      & ! in
                   dep_inputs%cldt, dep_inputs%cldcu,                           & ! in
                   dep_inputs%cldvcu, dep_inputs%cldvst,                        & ! in
                   sol_factb, sol_facti, sol_factic,                            & ! in
                   mam_prevap_resusp_optcc, .true., scavcoefnv(:,:,jnv), f_act_conv, & ! in
                   fldcw, qqcw_tmp,                                             & ! in 
                   fracis_cw, dqdt_tmp, iscavt,                                 & ! out
                   icscavt, isscavt, bcscavt, bsscavt, rcscavt, rsscavt         ) ! out 

                ! resuspension goes to coarse mode
                call calc_resusp_to_coarse(    ncol,   mm,      & ! in
                        mmtoo_prevap_resusp,   .false.,         & ! in
                        rcscavt,        rsscavt,                & ! in
                        dqdt_tmp,       rtscavt_sv              ) ! inout
                   
                fldcw(1:ncol,:) = fldcw(1:ncol,:) + dqdt_tmp(1:ncol,:) * dt

                call calc_sfc_flux(dqdt_tmp, state%pdel, sflx)
                call outfld( trim(cnst_name_cw(mm))//'SFWET', sflx, pcols, lchnk)
                aerdepwetcw(:ncol,mm) = sflx(:ncol)
                   
                call calc_sfc_flux(icscavt, state%pdel, sflx)
                call outfld( trim(cnst_name_cw(mm))//'SFSIC', sflx, pcols, lchnk)

                call calc_sfc_flux(isscavt, state%pdel, sflx)
                call outfld( trim(cnst_name_cw(mm))//'SFSIS', sflx, pcols, lchnk)

                call calc_sfc_flux(bcscavt, state%pdel, sflx)
                call outfld( trim(cnst_name_cw(mm))//'SFSBC', sflx, pcols, lchnk)

                call calc_sfc_flux(bsscavt, state%pdel, sflx)
                call outfld( trim(cnst_name_cw(mm))//'SFSBS', sflx, pcols, lchnk)

                call calc_sfc_flux(rcscavt, state%pdel, sflx)
                call outfld( trim(cnst_name_cw(mm))//'SFSEC', sflx, pcols, lchnk)
                      
                call calc_sfc_flux(rsscavt, state%pdel, sflx)
                call outfld( trim(cnst_name_cw(mm))//'SFSES', sflx, pcols, lchnk)

             endif lphase_jnmw_conditional

          enddo lspec_loop_aa ! lspec = 1, nspec_amode(m)+2
       enddo lphase_loop_aa ! lphase = 1, 2
    enddo mmode_loop_aa ! m = 1, ntot_amode

    ! if the user has specified prescribed aerosol dep fluxes then
    ! do not set cam_out dep fluxes according to the prognostic aerosols
    if (.not.aerodep_flx_prescribed()) then
       call set_srf_wetdep(aerdepwetis, aerdepwetcw, cam_out)
    endif

       
    call pbuf_get_field(pbuf, icwmrdp_idx,     icwmrdp )
    call pbuf_get_field(pbuf, icwmrsh_idx,     icwmrsh )
    call pbuf_get_field(pbuf, sh_frac_idx,     sh_frac )
    call pbuf_get_field(pbuf, dp_frac_idx,     dp_frac )

    call t_startf('ma_convproc')
    call ma_convproc_intr( state, dt,                      & ! in
            dp_frac, icwmrdp, rprddp, evapcdp,                & ! in
            sh_frac, icwmrsh, rprdsh, evapcsh,                & ! in
            dlf, dlf2, cmfmc2, sh_e_ed_ratio,                 & ! in
            nsrflx_mzaer2cnvpr, qsrflx_mzaer2cnvpr,           & ! in
            mu, md, du, eu, ed, dp, jt, maxg,                 & ! in
            ideep, lengath,  species_class,                   & ! in
            ptend, aerdepwetis                                ) ! inout
    call t_stopf('ma_convproc')       

    call wetdep_inputs_unset(dep_inputs)

  end subroutine aero_model_wetdep

!=============================================================================
   subroutine examine_prec_exist ( ncol,  pdel,      & ! in
                        prain,  cmfdqr,   evapr,     & ! in
                        isprx                        ) ! out
     !-----------------------------------------------------------------------
     ! examine if each grid has precipitation exist from all layers above
     !-----------------------------------------------------------------------
     integer,  intent(in) :: ncol
     real(r8), intent(in) :: pdel(:,:)      ! pressure difference between two layers [Pa]
     real(r8), intent(in) :: prain(pcols,pver)    ! rain production rate from stratiform clouds [kg/kg/s]
     real(r8), intent(in) :: cmfdqr(pcols,pver)   ! dq/dt due to convective rainout [kg/kg/s]
     real(r8), intent(in) :: evapr(pcols,pver)    ! rain evaporation rate [kg/kg/s]
     logical, intent(out) :: isprx(pcols,pver)    ! if there is precipitation falling into this grid

     ! local variables
     integer  :: kk
     real(r8) :: prec(pcols)    ! precipitation falling from layers above [kg/m2/s]
     real(r8), parameter :: small_value_7=1.0e-7_r8

     ! initiate precipitation at the top level
     prec(:ncol)=0._r8
     do kk=1,pver  ! check from the top level downward
       where (prec(:ncol) >= small_value_7)
          isprx(:ncol,kk) = .true.
       elsewhere
          isprx(:ncol,kk) = .false.
       endwhere
       ! update precipitation to the level below kk
       prec(:ncol) = prec(:ncol) + &
           (prain(:ncol,kk)+cmfdqr(:ncol,kk)-evapr(:ncol,kk)) * pdel(:ncol,kk)/gravit
     enddo

   end subroutine examine_prec_exist

!=============================================================================
   subroutine set_f_act_coarse( ncol,                           & ! in
                state_q,        ptend_q,        dt,             & ! in
                f_act_conv_coarse, f_act_conv_coarse_dust,      & ! out
                f_act_conv_coarse_nacl                          ) ! out
     !-----------------------------------------------------------------------
     ! set the mass-weighted sol_factic for coarse mode species
     !-----------------------------------------------------------------------
     use modal_aero_data, only: modeptr_coarse, &
                   lptr_dust_a_amode, lptr_nacl_a_amode

     integer, intent(in)  :: ncol
     ! state_q and ptend_q only use dust and seasalt in this subroutine
     real(r8),intent(in)  :: state_q(:,:,:)     ! tracer of state%q [kg/kg]
     real(r8),intent(in)  :: ptend_q(:,:,:)     ! tracer tendency (ptend%q) [kg/kg/s]
     real(r8),intent(in)  :: dt                 ! time step [s]
     real(r8),intent(out) :: f_act_conv_coarse(pcols,pver) ! prescribed coarse mode aerosol activation fraction for convective cloud [fraction]
     real(r8),intent(out) :: f_act_conv_coarse_dust ! prescribed dust aerosol activation fraction for convective cloud [fraction]
     real(r8),intent(out) :: f_act_conv_coarse_nacl ! prescribed seasalt aerosol activation fraction for convective cloud [fraction]
    
     ! local variables
     integer  :: lcoardust, lcoarnacl ! indices for coarse mode dust and seasalt masses
     integer  :: kk,icol
     real(r8) :: tmpdust, tmpnacl    ! dust and seasalt mass concentration [kg/kg] 
     real(r8), parameter :: small_value_30 = 1.0e-30_r8

     ! initial value
     f_act_conv_coarse(:,:) = 0.60_r8
     f_act_conv_coarse_dust = 0.40_r8
     f_act_conv_coarse_nacl = 0.80_r8

     if (modeptr_coarse > 0) then
       lcoardust = lptr_dust_a_amode(modeptr_coarse)
       lcoarnacl = lptr_nacl_a_amode(modeptr_coarse)
       if ((lcoardust > 0) .and. (lcoarnacl > 0)) then
          do kk = 1, pver
             do icol = 1, ncol
                tmpdust = max( 0.0_r8, state_q(icol,kk,lcoardust) + ptend_q(icol,kk,lcoardust)*dt )
                tmpnacl = max( 0.0_r8, state_q(icol,kk,lcoarnacl) + ptend_q(icol,kk,lcoarnacl)*dt )
                if ((tmpdust+tmpnacl) > small_value_30) then
                   f_act_conv_coarse(icol,kk) = (f_act_conv_coarse_dust*tmpdust &
                        + f_act_conv_coarse_nacl*tmpnacl)/(tmpdust+tmpnacl)
                endif
             enddo
          enddo
       endif
     endif

   end subroutine set_f_act_coarse
!=============================================================================
   subroutine define_act_frac ( lphase, imode,          & ! in
                sol_facti, sol_factic, sol_factb, f_act_conv) ! out
     !-----------------------------------------------------------------------
     ! define sol_factb and sol_facti values, and f_act_conv
     ! sol_factb - currently this is basically a tuning factor
     ! sol_facti & sol_factic - currently has a physical basis, and
     ! reflects activation fraction
     ! f_act_conv is the activation fraction
     !
     ! 2008-mar-07 rce - sol_factb (interstitial) changed from 0.3 to 0.1
     ! - sol_factic (interstitial, dust modes) changed from 1.0 to 0.5
     ! - sol_factic (cloud-borne, pcarb modes) no need to set it to 0.0
     ! because the cloud-borne pcarbon == 0 (no activation)
     !
     ! rce 2010/05/02
     ! prior to this date, sol_factic was used for convective in-cloud wet removal,
     ! and its value reflected a combination of an activation fraction
     ! (which varied between modes) and a tuning factor
     ! from this date forward, two parameters are used for convective
     ! in-cloud wet removal
     !
     ! note that "non-activation" of aerosol in air entrained into updrafts should
     ! be included here
     ! eventually we might use the activate routine (with w ~= 1 m/s) to calculate
     ! this, but there is still the entrainment issue
     !
     ! sol_factic is strictly a tuning factor
     !-----------------------------------------------------------------------
     use modal_aero_data, only: modeptr_pcarbon

     integer, intent(in) :: lphase ! index for interstitial / cloudborne aerosol
     integer, intent(in) :: imode  ! index for aerosol mode

     real(r8), intent(out) :: sol_facti  ! in-cloud scavenging fraction
     real(r8), intent(out) :: sol_factb  ! below-cloud scavenging fraction
     real(r8), intent(out) :: sol_factic(pcols,pver) ! in-cloud convective scavenging fraction
     real(r8), intent(out) :: f_act_conv(pcols,pver) ! convection activation fraction

     if (lphase == 1) then ! interstial aerosol
        sol_facti = 0.0_r8 ! strat in-cloud scav totally OFF for institial
        ! if modal aero convproc is turned on for aerosols, then
        ! turn off the convective in-cloud removal for interstitial aerosols
        ! (but leave the below-cloud on, as convproc only does in-cloud)
        ! and turn off the outfld SFWET, SFSIC, SFSID, SFSEC, and SFSED calls
        ! for (stratiform)-cloudborne aerosols, convective wet removal
        ! (all forms) is zero, so no action is needed
        sol_factic(:,:) = 0.0_r8
        sol_factb  = 0.03_r8   ! all below-cloud scav ON (0.1 "tuning factor") 
        if (imode == modeptr_pcarbon) then
           f_act_conv(:,:) = 0.0_r8 
        else
           f_act_conv(:,:) = 0.4_r8 
        endif
     else ! cloud-borne aerosol (borne by stratiform cloud drops)
        sol_factb  = 0.0_r8   ! all below-cloud scav OFF (anything cloud-borne is located "in-cloud")
        sol_facti  = min(0.6_r8, sol_facti_cloud_borne)  ! strat  in-cloud scav totally ON for cloud-borne
        sol_factic(:,:) = 0.0_r8   ! conv   in-cloud scav OFF (having this on would mean
                              ! that conv precip collects strat droplets)
        f_act_conv(:,:) = 0.0_r8   ! conv   in-cloud scav OFF (having this on would mean
     endif

   end subroutine define_act_frac

!=============================================================================
   subroutine index_ordering (                    &
                lspec,     imode,     lphase,     & ! in
                mm,        jnv,       jnummaswtr  ) ! out
     !-----------------------------------------------------------------------
     ! changed ordering (mass then number) for prevap resuspend to coarse
     !-----------------------------------------------------------------------
     use modal_aero_data, only: &
       lmassptr_amode, lmassptrcw_amode, &
       nspec_amode,  numptr_amode, numptrcw_amode

     integer, intent(in) :: lspec  ! index for aerosol number / chem-mass / water-mass
     integer, intent(in) :: imode  ! index for aerosol mode
     integer, intent(in) :: lphase ! index for interstitial / cloudborne aerosol 
     integer, intent(out) :: mm         ! index of the tracers
     integer, intent(out) :: jnv        ! index for scavcoefnv 3rd dimension
     integer, intent(out) :: jnummaswtr ! indicates current aerosol species type (0 = number, 1 = dry mass, 2 = water) 

     if (lspec <= nspec_amode(imode)) then ! non-water mass
          jnummaswtr = 1
          if (lphase == 1) then
             mm = lmassptr_amode(lspec,imode)
             jnv = 2
          else
             mm = lmassptrcw_amode(lspec,imode)
             jnv = 0
          endif
     elseif (lspec == nspec_amode(imode)+1) then ! number
          jnummaswtr = 0
          if (lphase == 1) then
             mm = numptr_amode(imode)
             jnv = 1
          else
             mm = numptrcw_amode(imode)
             jnv = 0
          endif
     else ! water mass
          jnummaswtr = 2
     endif

   end subroutine index_ordering

!=============================================================================
   subroutine apportion_sfc_flux_deep ( ncol,                     & ! in
                        rprddpsum,rprdshsum,evapcdpsum,evapcshsum,& ! in
                        sflxbc,   sflxec,                         & ! in
                        sflxbcdp, sflxecdp                        ) ! out
     !-----------------------------------------------------------------------
     ! apportion convective surface fluxes to deep and shallow conv
     ! this could be done more accurately in subr wetdepa
     ! since deep and shallow rarely occur simultaneously, and these
     ! fields are just diagnostics, this approximate method is adequate
     ! only do this for interstitial aerosol, because conv clouds to not
     ! affect the stratiform-cloudborne aerosol
     !-----------------------------------------------------------------------
     use mam_support,     only: min_max_bound

     integer, intent(in) :: ncol
     real(r8),intent(in) :: rprddpsum(:)  ! vertical integration of deep rain production
     real(r8),intent(in) :: rprdshsum(:)  ! vertical integration of shallow rain production
     real(r8),intent(in) :: evapcdpsum(:) ! vertical integration of deep rain evaporation
     real(r8),intent(in) :: evapcshsum(:) ! vertical integration of shallow rain evaporation
     real(r8),intent(in) :: sflxbc(:)     ! surface flux of resuspension from bcscavt [kg/m2/s]
     real(r8),intent(in) :: sflxec(:)     ! surface flux of resuspension from rcscavt [kg/m2/s]
     real(r8),intent(out):: sflxbcdp(:)   ! surface flux of resuspension from bcscavt in deep conv. [kg/m2/s]
     real(r8),intent(out):: sflxecdp(:)   ! surface flux of resuspension from rcscavt in deep conv. [kg/m2/s]


     integer :: ii
     real(r8) :: tmp_precdp, tmp_precsh, tmp_evapdp, tmp_evapsh ! working variables for precipitation and evaporation from deep and shallow convection
     real(r8) :: tmp_resudp, tmp_resush         ! working variables for resuspension from deep and shallow convection
     real(r8) :: tmpa, tmpb   ! working variables of deep fraction
     real(r8), parameter :: small_value_35 = 1.0e-35_r8
     real(r8), parameter :: small_value_36 = 1.0e-36_r8

     do ii = 1, ncol
          tmp_precdp = max( rprddpsum(ii),  small_value_35 )
          tmp_precsh = max( rprdshsum(ii),  small_value_35 )
          tmp_evapdp = max( evapcdpsum(ii), small_value_36 )
          tmp_evapsh = max( evapcshsum(ii), small_value_36 )

          ! assume that in- and below-cloud removal are proportional to
          ! column precip production
          tmpa = tmp_precdp / (tmp_precdp + tmp_precsh)
          tmpa = min_max_bound(0.0_r8, 1.0_r8, tmpa)
          sflxbcdp(ii) = sflxbc(ii)*tmpa

          ! assume that resuspension is proportional to
          ! (wet removal)*[(precip evap)/(precip production)]
          tmp_resudp =           tmpa  * min( (tmp_evapdp/tmp_precdp),1.0_r8 )
          tmp_resush = (1.0_r8 - tmpa) * min( (tmp_evapsh/tmp_precsh),1.0_r8 )
          tmpb = max( tmp_resudp, small_value_35 ) / &
                 max((tmp_resudp+tmp_resush), small_value_35)
          tmpb = min_max_bound(0.0_r8, 1.0_r8, tmpb)

          sflxecdp(ii) = sflxec(ii)*tmpb
      enddo

   end subroutine apportion_sfc_flux_deep

!=============================================================================
   subroutine calc_resusp_to_coarse(    ncol,   mm,     & ! in
                mmtoo_prevap_resusp,    update_dqdt,    & ! in
                rcscavt,        rsscavt,                & ! in
                dqdt_tmp,       rtscavt_sv              ) ! inout
     !-----------------------------------------------------------------------
     ! resuspension goes to coarse mode
     !-----------------------------------------------------------------------
     integer, intent(in) :: ncol, mm
     integer, intent(in) :: mmtoo_prevap_resusp(:)  ! pointers for resuspension
     logical, intent(in) :: update_dqdt  ! if update dqdt_tmp with rtscavt_sv
     real(r8),intent(in) :: rcscavt(:,:) ! resuspention from convective [kg/kg/s]
     real(r8),intent(in) :: rsscavt(:,:) ! resuspention from stratiform [kg/kg/s]

     real(r8), intent(inout) :: dqdt_tmp(:,:) ! temporary array to hold tendency for the "current" aerosol species [kg/kg/s]
     real(r8), intent(inout) :: rtscavt_sv(:,:,:) ! resuspension that goes to coarse mode [kg/kg/s]

     integer :: mmtoo

     mmtoo = mmtoo_prevap_resusp(mm)

     ! first deduct the current resuspension from the dqdt_tmp of the current species
     dqdt_tmp(1:ncol,:) = dqdt_tmp(1:ncol,:) - ( rcscavt(1:ncol,:) + rsscavt(1:ncol,:) )
     ! then add the current resuspension to the rtscavt_sv of the appropriate coarse mode species
     if (mmtoo > 0) then
        rtscavt_sv(1:ncol,:,mmtoo) = rtscavt_sv(1:ncol,:,mmtoo) &
                                   + ( rcscavt(1:ncol,:) + rsscavt(1:ncol,:) )
     endif
     ! then add the rtscavt_sv of the current species to the dqdt_tmp
     ! of the current species. This is not called when lphase==2
     ! note that for so4_a3 and mam3, the rtscavt_sv at this point will have
     !  resuspension contributions from so4_a1/2/3 and so4c1/2/3
     if (update_dqdt) then
        dqdt_tmp(1:ncol,:) = dqdt_tmp(1:ncol,:) + rtscavt_sv(1:ncol,:,mm)
     endif

  end subroutine calc_resusp_to_coarse

  !=============================================================================
  subroutine calc_sfc_flux(layer_tend, pdel, sflx)
    !-----------------------------------------------------------------------
    ! calculate surface fluxes of wet deposition from vertical integration of tendencies 
    !-----------------------------------------------------------------------
    real(r8), intent(in) :: pdel(:,:)      ! pressure difference between two layers [Pa]
    real(r8), intent(in) :: layer_tend(:,:)! physical tendencies in each layer [kg/kg/s]
    real(r8), intent(out):: sflx(:)        ! integrated surface fluxes [kg/m2/s]

    integer :: kk

     sflx(:)=0.0_r8
     do kk=1,pver
        sflx(:) = sflx(:) + layer_tend(:,kk)*pdel(:,kk)/gravit
     enddo

  end subroutine calc_sfc_flux
  !=============================================================================
  !=============================================================================
  subroutine aero_model_gasaerexch( loffset, ncol, lchnk, delt, &
                                    latndx, lonndx, reaction_rates, &
                                    tfld, pmid, pdel, mbar, relhum, &
                                    zm,  qh2o, cwat, cldfr, cldnum, &
                                    airdens, invariants, del_h2so4_gasprod,  &
                                    vmr0, vmr, pbuf )

    use time_manager,          only : get_nstep
    use modal_aero_amicphys,   only : modal_aero_amicphys_intr
    use mo_setsox,             only : setsox, has_sox
    use modal_aero_data,       only : qqcw_get_field

    !-----------------------------------------------------------------------
    !      ... dummy arguments
    !-----------------------------------------------------------------------
    integer,  intent(in) :: loffset                ! offset applied to modal aero "pointers"
    integer,  intent(in) :: ncol                   ! number columns in chunk
    integer,  intent(in) :: lchnk                  ! chunk index
    integer,  intent(in) :: latndx(pcols)          ! latitude indices
    integer,  intent(in) :: lonndx(pcols)          ! longitude indices
    real(r8), intent(in) :: delt                   ! time step size (sec)
    real(r8), intent(in) :: reaction_rates(:,:,:)  ! reaction rates
    real(r8), intent(in) :: tfld(:,:)              ! temperature (K)
    real(r8), intent(in) :: pmid(:,:)              ! pressure at model levels (Pa)
    real(r8), intent(in) :: pdel(:,:)              ! pressure thickness of levels (Pa)
    real(r8), intent(in) :: mbar(:,:)              ! mean wet atmospheric mass ( amu )
    real(r8), intent(in) :: relhum(:,:)            ! relative humidity
    real(r8), intent(in) :: airdens(:,:)           ! total atms density (molec/cm**3)
    real(r8), intent(in) :: invariants(:,:,:)
    real(r8), intent(in) :: del_h2so4_gasprod(:,:) 
    real(r8), intent(in) :: zm(:,:) 
    real(r8), intent(in) :: qh2o(:,:) 
    real(r8), intent(in) :: cwat(:,:)          ! cloud liquid water content (kg/kg)
    real(r8), intent(in) :: cldfr(:,:) 
    real(r8), intent(in) :: cldnum(:,:)       ! droplet number concentration (#/kg)
    real(r8), intent(in) :: vmr0(:,:,:)       ! initial mixing ratios (before gas-phase chem changes)
    real(r8), intent(inout) :: vmr(:,:,:)         ! mixing ratios ( vmr )
    type(physics_buffer_desc), pointer :: pbuf(:)
    
    ! local vars 
    
    integer :: n, m
    integer :: i,k
    integer :: nstep

    real(r8) :: del_h2so4_aeruptk(ncol,pver)

    real(r8), pointer :: dgnum(:,:,:), dgnumwet(:,:,:), wetdens(:,:,:)
    real(r8), pointer :: pblh(:)                    ! pbl height (m)

    real(r8), dimension(ncol) :: wrk
    character(len=32)         :: name
    real(r8) :: dvmrcwdt(ncol,pver,gas_pcnst)
    real(r8) :: dvmrdt(ncol,pver,gas_pcnst)
    real(r8) :: vmrcw(ncol,pver,gas_pcnst)            ! cloud-borne aerosol (vmr)

    real(r8), pointer :: fldcw(:,:)

    call pbuf_get_field(pbuf, dgnum_idx,      dgnum,  start=(/1,1,1/), kount=(/pcols,pver,ntot_amode/) )
    call pbuf_get_field(pbuf, dgnumwet_idx,   dgnumwet )
    call pbuf_get_field(pbuf, wetdens_ap_idx, wetdens )
    call pbuf_get_field(pbuf, pblh_idx,       pblh)

    do n=1,ntot_amode
       call outfld(dgnum_name(n),dgnumwet(1:ncol,1:pver,n), ncol, lchnk )
    enddo

! do gas-aerosol exchange (h2so4, msa, nh3 condensation)

    nstep = get_nstep()

    ! calculate tendency due to gas phase chemistry and processes
    dvmrdt(:ncol,:,:) = (vmr(:ncol,:,:) - vmr0(:ncol,:,:)) / delt
    do m = 1, gas_pcnst
      wrk(:) = 0._r8
      do k = 1,pver
        wrk(:ncol) = wrk(:ncol) + dvmrdt(:ncol,k,m)*adv_mass(m)/mbar(:ncol,k)*pdel(:ncol,k)/gravit
      enddo
      name = 'GS_'//trim(solsym(m))
      call outfld( name, wrk(:ncol), ncol, lchnk )
    enddo

!
! Aerosol processes ...
!
    call qqcw2vmr( lchnk, vmrcw, mbar, ncol, loffset, pbuf )

    !------------------------------------------------------

      dvmrdt(:ncol,:,:) = vmr(:ncol,:,:)
      dvmrcwdt(:ncol,:,:) = vmrcw(:ncol,:,:)

      ! aqueous chemistry ...

      if( has_sox ) then
         call setsox(   &
              ncol,     &
              lchnk,    &
              loffset,  &
              delt,     &
              pmid,     &
              pdel,     &
              tfld,     &
              mbar,     &
              cwat,     &
              cldfr,    &
              cldnum,   &
              airdens,  &
              invariants, &
              vmrcw,    &
              vmr       &
              )
      endif

      ! Tendency due to aqueous chemistry 
      ! before aqueous chemistry, and cannot be used to hold aq. chem. tendencies
      ! ***Note - should calc & output tendencies for cloud-borne aerosol species 
      !           rather than interstitial here
      do m = 1, gas_pcnst
        wrk(:) = 0._r8
        do k = 1,pver
            ! here dvmrdt is vmr before aqueous chemistry, so need to calculate (delta vmr)/(delt)
            wrk(:ncol) = wrk(:ncol) + ((vmr(:ncol,k,m)-dvmrdt(:ncol,k,m))/delt) &
                                                        * adv_mass(m)/mbar(:ncol,k)*pdel(:ncol,k)/gravit
        enddo
        name = 'AQ_'//trim(solsym(m))
        call outfld( name, wrk(:ncol), ncol, lchnk )
      enddo

    ! do gas-aerosol exchange, nucleation, and coagulation using new routines

       call t_startf('modal_aero_amicphys')

       ! note that:
       !     vmr0 holds vmr before gas-phase chemistry
       !     dvmrdt and dvmrcwdt hold vmr and vmrcw before aqueous chemistry
       call modal_aero_amicphys_intr(                &
            1,                  1,                   &
            1,                  1,                   &
            lchnk,     ncol,    nstep,               &
            loffset,   delt,                         &
            latndx,    lonndx,                       &
            tfld,      pmid,    pdel,                &
            zm,        pblh,                         &
            qh2o,      cldfr,                        &
            vmr,                vmrcw,               &
            vmr0,                                    &
            dvmrdt,             dvmrcwdt,            &
            dgnum,              dgnumwet,            &
            wetdens                                  )

       call t_stopf('modal_aero_amicphys')


    call vmr2qqcw( lchnk, vmrcw, mbar, ncol, loffset, pbuf )

    ! diagnostics for cloud-borne aerosols... 
    do n = 1,pcnst
       fldcw => qqcw_get_field(pbuf,n,lchnk,errorhandle=.true.)
       if(associated(fldcw)) then
          call outfld( cnst_name_cw(n), fldcw(:,:), pcols, lchnk )
       endif
    end do

  end subroutine aero_model_gasaerexch

  !=============================================================================
  !=============================================================================
  subroutine aero_model_emissions( state, cam_in )
    use seasalt_model, only: seasalt_emis, seasalt_names, seasalt_indices, seasalt_active,seasalt_nbin, &
         has_mam_mom, nslt_om
    use dust_model,    only: dust_emis, dust_names, dust_indices, dust_active,dust_nbin, dust_nnum

    ! Arguments:

    type(physics_state),    intent(in)    :: state   ! Physics state variables
    type(cam_in_t),         intent(inout) :: cam_in  ! import state

    ! local vars

    integer :: lchnk, ncol
    integer :: m, mm
    real(r8) :: soil_erod_tmp(pcols)
    real(r8) :: sflx(pcols)   ! accumulate over all bins for output
    real(r8) :: u10cubed(pcols)
    real(r8) :: u10(pcols)               ! Needed in Gantt et al. calculation of organic mass fraction
    real(r8) :: F_eff(pcols) ! optional diagnostic output -- organic enrichment ratio

    real (r8), parameter :: z0=0.0001_r8  ! m roughness length over oceans--from ocean model

    lchnk = state%lchnk
    ncol = state%ncol

    if (dust_active) then

       call dust_emis( ncol, lchnk, cam_in%dstflx, cam_in%cflx, soil_erod_tmp )

       ! some dust emis diagnostics ...
       sflx(:)=0._r8
       do m=1,dust_nbin+dust_nnum
          mm = dust_indices(m)
          if (m<=dust_nbin) sflx(:ncol)=sflx(:ncol)+cam_in%cflx(:ncol,mm)
          call outfld(trim(dust_names(m))//'SF',cam_in%cflx(:,mm),pcols, lchnk)
       enddo
       call outfld('DSTSFMBL',sflx(:),pcols,lchnk)
       call outfld('LND_MBL',soil_erod_tmp(:),pcols, lchnk )
    endif

    if (seasalt_active) then
       u10(:ncol)=sqrt(state%u(:ncol,pver)**2+state%v(:ncol,pver)**2)
       ! move the winds to 10m high from the midpoint of the gridbox:
       ! follows Tie and Seinfeld and Pandis, p.859 with math.

       u10cubed(:ncol)=u10(:ncol)*log(10._r8/z0)/log(state%zm(:ncol,pver)/z0)

       ! we need them to the 3.41 power, according to Gong et al., 1997:
       u10cubed(:ncol)=u10cubed(:ncol)**3.41_r8

       sflx(:)=0._r8
       F_eff(:)=0._r8

       call seasalt_emis(u10, u10cubed, lchnk, cam_in%sst, cam_in%ocnfrac, ncol, cam_in%cflx, seasalt_emis_scale, F_eff)

       ! Write out salt mass fluxes to history files
       do m=1,seasalt_nbin-nslt_om
          mm = seasalt_indices(m)
          sflx(:ncol)=sflx(:ncol)+cam_in%cflx(:ncol,mm)
          call outfld(trim(seasalt_names(m))//'SF',cam_in%cflx(:,mm),pcols,lchnk)
       enddo
       ! accumulated flux
       call outfld('SSTSFMBL',sflx(:),pcols,lchnk)

       ! Write out marine organic mass fluxes to history files
       if ( has_mam_mom ) then
          sflx(:)=0._r8
          do m=seasalt_nbin-nslt_om+1,seasalt_nbin
             mm = seasalt_indices(m)
             sflx(:ncol)=sflx(:ncol)+cam_in%cflx(:ncol,mm)
             call outfld(trim(seasalt_names(m))//'SF',cam_in%cflx(:,mm),pcols,lchnk)
          enddo
          ! accumulated flux
          call outfld('SSTSFMBL_OM',sflx(:),pcols,lchnk)

       endif

    endif

  end subroutine aero_model_emissions

  !===============================================================================
  ! private methods


  !===============================================================================
  !===============================================================================
  subroutine modal_aero_bcscavcoef_init
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Computes lookup table for aerosol impaction/interception scavenging rates
    !
    ! Authors: R. Easter
    !
    !-----------------------------------------------------------------------
    
    use modal_aero_data, only: dgnum_amode, sigmag_amode, lspectype_amode, specdens_amode

    implicit none

    !   local variables
    integer :: jgrow, ll, imode ! indexes
    real(r8) :: dg0             ! aerosol diameter [m]
    real(r8) :: dg0_cgs         ! aerosol diameter in CGS unit [cm]
    real(r8) :: sigmag          ! standard deviation of aerosol size distribution 
    real(r8) :: rhodryaero, rhowetaero   ! dry and wet aerosol density [kg/m3] 
    real(r8) :: rhowetaero_cgs  ! wet aerosol density in CGS unit [g/cm3]
    real(r8) :: scavratenum     ! scavenging rate of aerosol number [1/s]
    real(r8) :: scavratevol     ! scavenging rate of aerosol volume [1/s]
    real(r8) :: wetdiaratio, wetvolratio ! ratio of diameter and volume for wet/dry aerosols [fraction]

    ! set up temperature-pressure pair to compute impaction scavenging rates
    real(r8), parameter :: temp_0C = 273.16_r8        ! K
    real(r8), parameter :: press_750hPa = 0.75e6_r8   ! dynes/cm2
    
    modeloop: do imode = 1, ntot_amode

       sigmag = sigmag_amode(imode)
       ll = lspectype_amode(1,imode)
       rhodryaero = specdens_amode(ll)

       growloop: do jgrow = nimptblgrow_mind, nimptblgrow_maxd

          wetdiaratio = exp( jgrow*dlndg_nimptblgrow )
          dg0 = dgnum_amode(imode)*wetdiaratio
          wetvolratio = exp( jgrow*dlndg_nimptblgrow*3._r8 )
          rhowetaero = 1.0_r8 + (rhodryaero-1.0_r8)/wetvolratio
          rhowetaero = min( rhowetaero, rhodryaero )

! FIXME: not sure why wet aerosol density is set as dry aerosol density here
! but the above calculation of rhowetaero is incorrect. 
! I think the number 1.0_r8 should be 1000._r8 as the unit is kg/m3
! the above calculation gives wet aerosol density very small number (a few kg/m3)
! this may cause some problem. I guess this is the reason of using dry density.
! should be better if fix the wet density bug and use it. Keep it for now for BFB testing
! -- (commented by Shuaiqi Tang when refactoring for MAM4xx)
          rhowetaero = rhodryaero

          ! compute impaction scavenging rates at 1 temp-press pair and save
          ! note that the subroutine calc_1_impact_rate uses CGS units
          dg0_cgs = dg0*1.0e2_r8   ! m to cm
          rhowetaero_cgs = rhowetaero*1.0e-3_r8   ! kg/m3 to g/cm3
          call calc_1_impact_rate( &
               dg0_cgs, sigmag, rhowetaero_cgs, temp_0C, press_750hPa, &
               scavratenum, scavratevol )

          scavimptblnum(jgrow,imode) = log( scavratenum )
          scavimptblvol(jgrow,imode) = log( scavratevol )

       enddo growloop
    enddo modeloop
    return

  end subroutine modal_aero_bcscavcoef_init

  !===============================================================================
  subroutine modal_aero_bcscavcoef_get( imode, ncol, isprx, dgn_awet, & ! in
                                        scavcoefnum, scavcoefvol      ) ! out
    !-----------------------------------------------------------------------
    ! compute impaction scavenging removal amount for aerosol volume and number
    !-----------------------------------------------------------------------

    use modal_aero_data, only: dgnum_amode

    implicit none

    integer,  intent(in) :: imode, ncol
    logical,  intent(in) :: isprx(pcols,pver)           ! if there is precip
    real(r8), intent(in) :: dgn_awet(pcols,pver,ntot_amode)  ! wet aerosol diameter [m]
    real(r8), intent(out):: scavcoefnum(pcols,pver)     ! scavenging removal for aerosol number [1/h]
    real(r8), intent(out):: scavcoefvol(pcols,pver)     ! scavenging removal for aerosol volume [1/h]

    ! local variables
    integer  :: icol, kk, jgrow         ! index
    real(r8) :: wetdiaratio             ! ratio of wet and dry aerosol diameter [fraction]
    real(r8) :: xgrow, dumfhi, dumflo   ! working variables
    real(r8) :: scavimpvol, scavimpnum  ! log of scavenging rates for volume and number


    do kk = 1, pver
       do icol = 1, ncol

          ! do only if there is precip
          if ( isprx(icol,kk) ) then
             
             ! interpolate table values using log of (actual-wet-size)/(base-dry-size)
             wetdiaratio = dgn_awet(icol,kk,imode)/dgnum_amode(imode)

             if ((wetdiaratio >= 0.99_r8) .and. (wetdiaratio <= 1.01_r8)) then
                scavimpvol = scavimptblvol(0,imode)
                scavimpnum = scavimptblnum(0,imode)
             else
                xgrow = log( wetdiaratio ) / dlndg_nimptblgrow
                jgrow = int( xgrow )
                if (xgrow < 0._r8) jgrow = jgrow - 1

                if (jgrow < nimptblgrow_mind) then
                   jgrow = nimptblgrow_mind
                   xgrow = jgrow
                else
                   jgrow = min( jgrow, nimptblgrow_maxd-1 )
                endif

                dumfhi = xgrow - jgrow
                dumflo = 1._r8 - dumfhi
                scavimpvol = dumflo*scavimptblvol(jgrow,imode) + &
                     dumfhi*scavimptblvol(jgrow+1,imode)
                scavimpnum = dumflo*scavimptblnum(jgrow,imode) + &
                     dumfhi*scavimptblnum(jgrow+1,imode)

             endif

             ! impaction scavenging removal amount for volume
             scavcoefvol(icol,kk) = exp( scavimpvol )
             ! impaction scavenging removal amount to number
             scavcoefnum(icol,kk) = exp( scavimpnum )

          else ! if no precip
             scavcoefvol(icol,kk) = 0._r8
             scavcoefnum(icol,kk) = 0._r8
          endif

       enddo
    enddo

    return
  end subroutine modal_aero_bcscavcoef_get

  !===============================================================================
  subroutine calc_1_impact_rate(                        &
                     dg0, sigmag, rhoaero, temp, press, & ! in
                     scavratenum, scavratevol           ) ! out
   !
   !   this subroutine computes a single impaction scavenging rate
   !   for precipitation rate of 1 mm/h
   !
   !-----------------------------------------------------------------
   
   use mo_constants, only: rgas => rgas_cgs

   implicit none

   !   subroutine parameters
   real(r8), intent(in)  :: dg0         ! geometric mean diameter of aerosol [cm]
   real(r8), intent(in)  :: sigmag      ! geometric standard deviation of size distribution
   real(r8), intent(in)  :: rhoaero     ! aerosol density [g/cm^3]
   real(r8), intent(in)  :: temp        ! temperature [K]
   real(r8), intent(in)  :: press       ! pressure [dyne/cm^2]
   real(r8), intent(out) :: scavratenum, scavratevol  ! scavenging rate for aerosol number and volume [1/hour]

   !   local variables
   integer, parameter :: lunerr = 6     ! logical unit for error message

   integer, parameter :: nrainsvmax=50  ! maximum bin number for rain
   real(r8) :: rrainsv(nrainsvmax)      ! rain radius for each bin [cm]
   real(r8) :: xnumrainsv(nrainsvmax)   ! rain number for each bin [#/cm3]
   real(r8) :: vfallrainsv(nrainsvmax)  ! rain falling velocity for each bin [cm/s]

   integer, parameter :: naerosvmax=51  ! maximum bin number for aerosol
   real(r8) :: raerosv(naerosvmax)      ! aerosol particle radius in each bin [cm] 
   real(r8) :: fnumaerosv(naerosvmax)   ! fraction of total number in the bin [fraction]
   real(r8) :: fvolaerosv(naerosvmax)   ! fraction of total volume in the bin [fraction]

   integer  :: ja, jr                   ! index for aerosol and rain
   integer  :: na, nr                   ! number of aerosol and rain bins
   real(r8) :: airkinvisc               ! air kinematic viscosity [cm^2/s]
   real(r8) :: freepath                 ! molecular freepath [cm]
   real(r8) :: etotal                   ! efficiency of total scavenging effects [fraction]
   real(r8) :: precip                   ! precipitation rate, fix as 1 mm/hr in this subroutine [cm/s]
   real(r8) :: r_rain, r_aer            ! rain droplet and aerosol particle radius [cm]
   real(r8) :: cair                     ! air molar density [dyne/cm^2/erg*mol = mol/cm^3]
   real(r8) :: rhoair                   ! air mass density [g/cm^3]
   real(r8) :: rlo, rhi, dr             ! rain droplet bin information [cm]
   real(r8) :: rainsweepout             ! rain droplet sweep out volume [cm3/cm3/s]
   real(r8) :: scavsumnum, scavsumnumbb ! scavenging rate of aerosol number, "*bb" is for each rain droplet radius bin [1/s] 
   real(r8) :: scavsumvol, scavsumvolbb ! scavenging rate of aerosol volume, "*bb" is for each rain droplet radius bin [1/s]
   real(r8) :: sx                       ! standard deviation (log-normal distribution)
   real(r8) :: vfall                    ! rain droplet fall speed [cm/s]
   real(r8) :: ag0, xg0, xg3, xhi, xlo,dx    ! aerosol bin information
   !-----------------------------------------------------------------

   ! this subroutine is calculated for a fix rainrate of 1 mm/hr
   precip = 1.0_r8/36000._r8        ! 1 mm/hr in cm/s

   ! set the iteration radius for rain droplet
   rlo = .005_r8
   rhi = .250_r8
   dr = 0.005_r8
   nr = 1 + nint( (rhi-rlo)/dr )
   if (nr .gt. nrainsvmax) then
      write(lunerr,*) '*** subr. calc_1_impact_rate -- nr > nrainsvmax'
      call endrun()
   endif

   ! aerosol modal information
   ag0 = dg0/2._r8      ! mean radius of aerosol
   sx = log( sigmag )   ! standard deviation (log-normal distribution)
   xg0 = log( ag0 )     ! log(mean radius) (log-normal distribution)
   xg3 = xg0 + 3._r8*sx*sx      ! mean + 3*std^2

   ! set the iteration radius for aerosol particles
   dx = max( 0.2_r8*sx, 0.01_r8 )
   xlo = xg3 - max( 4._r8*sx, 2._r8*dx )
   xhi = xg3 + max( 4._r8*sx, 2._r8*dx )
   na = 1 + nint( (xhi-xlo)/dx )
   if (na .gt. naerosvmax) then
      write(lunerr,*) '*** subr. calc_1_impact_rate -- na > naerosvmax'
      call endrun()
   endif

   !   air molar density [dyne/cm^2/erg*mol = mol/cm^3]
   cair = press/(rgas*temp)   
   !   air mass density [g/cm^3]
   rhoair = 28.966_r8*cair
   !   molecular freepath [cm]
   freepath = 2.8052e-10_r8/cair
   !   air kinematic viscosity
   airkinvisc = air_kinematic_viscosity( temp, rhoair )

   !
   !   compute rain drop number concentrations
   !
   call calc_rain_drop_conc( nr, rlo, dr, rhoair, precip, & ! in
                         rrainsv, xnumrainsv, vfallrainsv ) ! out

   !
   !   compute aerosol concentrations
   !
   call calc_aer_conc_frac( na, xlo, dx, xg0, sx,      & ! in
                        raerosv, fnumaerosv, fvolaerosv) ! out

   !
   !   compute scavenging
   !
   scavsumnum = 0._r8
   scavsumvol = 0._r8
   !   outer loop for rain drop radius
   jr_loop: do jr = 1, nr

      ! rain droplet radius
      r_rain = rrainsv(jr)
      vfall = vfallrainsv(jr)

      !   inner loop for aerosol particle radius
      scavsumnumbb = 0._r8
      scavsumvolbb = 0._r8
      ja_loop: do ja = 1, na
         ! aerosol particle radius
         r_aer = raerosv(ja)
         call calc_impact_efficiency( r_aer, r_rain, temp,  & ! in
                                 freepath, rhoaero, rhoair, & ! in
                                 vfall, airkinvisc,         & ! in
                                 etotal                     ) ! out

         rainsweepout = xnumrainsv(jr)*4._r8*pi*r_rain*r_rain*vfall
         scavsumnumbb = scavsumnumbb + rainsweepout*etotal*fnumaerosv(ja)
         scavsumvolbb = scavsumvolbb + rainsweepout*etotal*fvolaerosv(ja)
      enddo ja_loop

      scavsumnum = scavsumnum + scavsumnumbb
      scavsumvol = scavsumvol + scavsumvolbb

   enddo jr_loop

   scavratenum = scavsumnum*3600._r8
   scavratevol = scavsumvol*3600._r8

   return
  end subroutine calc_1_impact_rate
 
  !=============================================================================
  subroutine calc_rain_drop_conc( nr, rlo, dr, rhoair, precip, & ! in
                              rrainsv, xnumrainsv, vfallrainsv ) ! out
   !-----------------------------------------------------------------
   !   compute rain drop number concentrations, radius and falling velocity
   !-----------------------------------------------------------------
   
   integer,  intent(in) :: nr           ! number of rain bins
   real(r8), intent(in) :: rlo          ! lower limit of rain radius [cm]
   real(r8), intent(in) :: dr           ! rain radius bin width [cm]
   real(r8), intent(in) :: rhoair       ! air mass density [g/cm^3]
   real(r8), intent(in) :: precip       ! precipitation [cm/s] 

   real(r8), intent(out) :: rrainsv(:)  ! rain radius in each bin [cm]
   real(r8), intent(out) :: xnumrainsv(:)  ! rain number concentration in each bin [#/cm3]
   real(r8), intent(out) :: vfallrainsv(:) ! rain droplet falling velocity [cm/s]

   ! local variables
   integer  :: ii                       ! index of cloud bins
   real(r8) :: precipsum                ! sum of precipitation in all bins
   real(r8) :: rr                       ! rain radius in the bin [cm]
   real(r8) :: dd                       ! rain diameter in the bin [cm]
   real(r8) :: vfall, vfallstp          ! rain droplet falling speed [cm/s]

   precipsum = 0._r8
   do ii = 1, nr
      rr = rlo + (ii-1)*dr
      rrainsv(ii) = rr
      xnumrainsv(ii) = exp( -rr/2.7e-2_r8 )

      dd = 2._r8*rr
      if (dd .le. 0.007_r8) then
         vfallstp = 2.88e5_r8 * dd**2._r8
      elseif (dd .le. 0.025_r8) then
         vfallstp = 2.8008e4_r8 * dd**1.528_r8
      elseif (dd .le. 0.1_r8) then
         vfallstp = 4104.9_r8 * dd**1.008_r8
      elseif (dd .le. 0.25_r8) then
         vfallstp = 1812.1_r8 * dd**0.638_r8
      else
         vfallstp = 1069.8_r8 * dd**0.235_r8
      endif

      vfall = vfallstp * sqrt(1.204e-3_r8/rhoair)
      vfallrainsv(ii) = vfall
      precipsum = precipsum + vfall*(rr**3)*xnumrainsv(ii)
   enddo
  ! 1.333333 is simplified 4/3 for sphere volume calculation
   precipsum = precipsum*pi*1.333333_r8

   do ii = 1, nr
      xnumrainsv(ii) = xnumrainsv(ii)*(precip/precipsum)
   enddo

  end subroutine calc_rain_drop_conc

  !=============================================================================
  subroutine calc_aer_conc_frac( na, xlo, dx, xg0, sx,      & ! in
                             raerosv, fnumaerosv, fvolaerosv) ! out
  !-----------------------------------------------------------------
  !   compute aerosol concentration, radius and volume in each bin
  !-----------------------------------------------------------------
  
   integer,  intent(in) :: na           ! number of aerosol bins
   real(r8), intent(in) :: xlo          ! lower limit of aerosol radius (log)
   real(r8), intent(in) :: dx           ! aerosol radius bin width (log)
   real(r8), intent(in) :: xg0          ! log(mean radius)
   real(r8), intent(in) :: sx           ! standard deviation (log)

   real(r8), intent(out):: raerosv(:)   ! aerosol radius [cm]
   real(r8), intent(out):: fnumaerosv(:)! fraction of total number in the bin [fraction]
   real(r8), intent(out):: fvolaerosv(:)! fraction of total volume in the bin [fraction]

   ! local variables
   integer  :: ii                       ! index of aerosol bins
   real(r8) :: xx                       ! aerosol radius in the bin (log)
   real(r8) :: aa                       ! aerosol radius in the bin [cm]
   real(r8) :: dum                      ! working variable
   real(r8) :: anumsum, avolsum         ! total aerosol number and volume
 
   ! calculate total aerosol number and volume
   anumsum = 0._r8
   avolsum = 0._r8
   do ii = 1, na
      xx = xlo + (ii-1)*dx
      aa = exp( xx )
      raerosv(ii) = aa
      dum = (xx - xg0)/sx
      fnumaerosv(ii) = exp( -0.5_r8*dum*dum )
      ! 1.3333 is simplified 4/3 for sphere volume calculation
      fvolaerosv(ii) = fnumaerosv(ii)*1.3333_r8*pi*aa*aa*aa 
      anumsum = anumsum + fnumaerosv(ii)
      avolsum = avolsum + fvolaerosv(ii)
   enddo

   ! calculate fraction in each aerosol bin
   do ii = 1, na
      fnumaerosv(ii) = fnumaerosv(ii)/anumsum
      fvolaerosv(ii) = fvolaerosv(ii)/avolsum
   enddo

  end subroutine calc_aer_conc_frac
 
  !=============================================================================
  subroutine calc_impact_efficiency( r_aer,  r_rain,  temp,     & ! in
                                     freepath, rhoaero, rhoair, & ! in
                                     vfall, airkinvisc,         & ! in
                                     etotal                     ) ! out
  !-----------------------------------------------------------------
  ! calculate aerosol-collection efficiency for a given radius of rain and aerosol particles
  !-----------------------------------------------------------------
 
      real(r8), intent(in)  :: r_aer         ! aerosol radius [cm]
      real(r8), intent(in)  :: r_rain        ! rain radius [cm]
      real(r8), intent(in)  :: temp          ! temperature [K]
      real(r8), intent(in)  :: freepath      ! molecular freepath [cm]
      real(r8), intent(in)  :: rhoaero       ! density of aerosol particles [g/cm^3]
      real(r8), intent(in)  :: rhoair        ! air mass density [g/cm^3]
      real(r8), intent(in)  :: airkinvisc    ! air kinematic viscosity [cm^2/s]
      real(r8), intent(in)  :: vfall         ! rain droplet falling speed [cm/s]
      real(r8), intent(out) :: etotal        ! efficiency of total effects [fraction]

      ! local variables
      real(r8)  :: chi  ! ratio of aerosol and rain radius [fraction]
      real(r8)  :: dum, sstar   ! working variables [unitless]
      real(r8)  :: taurelax     ! Stokes number relaxation time [s]
      real(r8)  :: schmidt      ! Schmidt number [unitless]
      real(r8)  :: stokes       ! Stokes number [unitless]
      real(r8)  :: reynolds     ! Raynold number [unitless]
      real(r8)  :: sqrtreynolds ! sqrt of Raynold number [unitless]
      real(r8)  :: ebrown, eintercept, eimpact ! efficiency of aerosol-collection in different processes

      ! ratio of water viscosity to air viscosity (from Slinn)
      real(r8), parameter ::  xmuwaterair = 60.0_r8 ! [fraction]


      chi = r_aer/r_rain

      ! ---------- calcualte Brown effect ------------

      ! calculate unitless numbers
      call calc_schmidt_number( freepath, r_aer, temp,        & ! in
                              rhoaero, rhoair, airkinvisc,    & ! in
                              schmidt,   taurelax             ) ! out
      stokes = vfall*taurelax/r_rain
      reynolds = r_rain * vfall / airkinvisc
      sqrtreynolds = sqrt( reynolds )

      ebrown = 4._r8*(1._r8 + 0.4_r8*sqrtreynolds*(schmidt**0.3333333_r8)) / (reynolds*schmidt)

      ! ------------ calculate intercept effect ------------
      dum = (1._r8 + 2._r8*xmuwaterair*chi) / (1._r8 + xmuwaterair/sqrtreynolds)
      eintercept = 4._r8*chi*(chi + dum)

      ! ------------ calculate impact effect ------------
      dum = log( 1._r8 + reynolds )
      sstar = (1.2_r8 + dum/12._r8) / (1._r8 + dum)
      eimpact = 0._r8
      if (stokes > sstar) then
            dum = stokes - sstar
            eimpact = (dum/(dum+0.6666667_r8)) ** 1.5_r8
      endif

      ! ------------ calculate total effects ------------
      etotal = ebrown + eintercept + eimpact
      etotal = min( etotal, 1.0_r8 )

  end subroutine calc_impact_efficiency

  !=============================================================================
  subroutine calc_schmidt_number( freepath, r_aer, temp,        & ! in
                                rhoaero, rhoair, airkinvisc,    & ! in
                                schmidt,   taurelax             ) ! out
    !-----------------------------------------------------------------
    ! calculate Schmidt number
    ! also output relaxation time for Stokes number
    !
    ! note that there is a similar calculation of Schmidt number in dry
    ! depositon (in modal_aero_drydep.F90) but the calculation of dumfuchs (or
    ! slip_correction_factor) looks differently
    !-----------------------------------------------------------------

    use mo_constants, only: boltz_cgs

      real(r8), intent(in)  :: freepath      ! molecular freepath [cm]
      real(r8), intent(in)  :: r_aer         ! aerosol radius [cm]
      real(r8), intent(in)  :: temp          ! temperature [K]
      real(r8), intent(in)  :: rhoaero       ! density of aerosol particles [g/cm^3]
      real(r8), intent(in)  :: rhoair        ! air mass density [g/cm^3]
      real(r8), intent(in)  :: airkinvisc    ! air kinematic viscosity [cm2/s]

      real(r8), intent(out) :: schmidt       ! Schmidt number [unitless]
      real(r8), intent(out) :: taurelax      ! relaxation time for Stokes number [s]

      ! local variables
      real(r8)  :: dum          ! working variables [unitless]
      real(r8)  :: dumfuchs     ! slip correction factor [unitless]
      real(r8)  :: aeromass     ! single-particle aerosol mass [g]
      real(r8)  :: aerodiffus   ! aerosol diffusivity [cm^2/s]

      dum = freepath/r_aer
      dumfuchs = 1._r8 + 1.246_r8*dum + 0.42_r8*dum*exp(-0.87_r8/dum)
      taurelax = 2._r8*rhoaero*r_aer*r_aer*dumfuchs/(9._r8*rhoair*airkinvisc)

      aeromass = 4._r8*pi*r_aer*r_aer*r_aer*rhoaero/3._r8 ![g]
      aerodiffus = boltz_cgs*temp*taurelax/aeromass  ! [cm^2/s]

      schmidt = airkinvisc/aerodiffus

  end subroutine calc_schmidt_number

  !=============================================================================
  subroutine qqcw2vmr(lchnk, vmr, mbar, ncol, im, pbuf)
    use modal_aero_data, only : qqcw_get_field
    !-----------------------------------------------------------------
    !	... Xfrom from mass to volume mixing ratio
    !-----------------------------------------------------------------

    implicit none

    !-----------------------------------------------------------------
    !	... Dummy args
    !-----------------------------------------------------------------
    integer, intent(in)     :: lchnk, ncol, im
    real(r8), intent(in)    :: mbar(ncol,pver)
    real(r8), intent(inout) :: vmr(ncol,pver,gas_pcnst)
    type(physics_buffer_desc), pointer :: pbuf(:)

    !-----------------------------------------------------------------
    !	... Local variables
    !-----------------------------------------------------------------
    integer :: k, m
    real(r8), pointer :: fldcw(:,:)

    do m=1,gas_pcnst
       if( adv_mass(m) /= 0._r8 ) then
          fldcw => qqcw_get_field(pbuf, m+im,lchnk,errorhandle=.true.)
          if(associated(fldcw)) then
             do k=1,pver
                vmr(:ncol,k,m) = mbar(:ncol,k) * fldcw(:ncol,k) / adv_mass(m)
             enddo
          else
             vmr(:,:,m) = 0.0_r8
          endif
       endif
    enddo
  end subroutine qqcw2vmr


  !=============================================================================
  !=============================================================================
  subroutine vmr2qqcw( lchnk, vmr, mbar, ncol, im, pbuf )
    !-----------------------------------------------------------------
    !	... Xfrom from volume to mass mixing ratio
    !-----------------------------------------------------------------

    use m_spc_id
    use modal_aero_data, only : qqcw_get_field

    implicit none

    !-----------------------------------------------------------------
    !	... Dummy args
    !-----------------------------------------------------------------
    integer, intent(in)     :: lchnk, ncol, im
    real(r8), intent(in)    :: mbar(ncol,pver)
    real(r8), intent(in)    :: vmr(ncol,pver,gas_pcnst)
    type(physics_buffer_desc), pointer :: pbuf(:)

    !-----------------------------------------------------------------
    !	... Local variables
    !-----------------------------------------------------------------
    integer :: k, m
    real(r8), pointer :: fldcw(:,:)
    !-----------------------------------------------------------------
    !	... The non-group species
    !-----------------------------------------------------------------
    do m = 1,gas_pcnst
       fldcw => qqcw_get_field(pbuf, m+im,lchnk,errorhandle=.true.)
       if( adv_mass(m) /= 0._r8 .and. associated(fldcw)) then
          do k = 1,pver
             fldcw(:ncol,k) = adv_mass(m) * vmr(:ncol,k,m) / mbar(:ncol,k)
          enddo
       endif
    enddo

  end subroutine vmr2qqcw

  !=============================================================================
  real(r8) function air_dynamic_viscosity( temp )
    !-----------------------------------------------------------------
    ! Calculate dynamic viscosity of air, unit [g/cm/s]
    !
    ! note that this calculation is different with that used in dry deposition
    ! see the same-name function in modal_aero_drydep.F90 
    !-----------------------------------------------------------------

    real(r8),intent(in) :: temp   ! air temperature [K]

    air_dynamic_viscosity = 1.8325e-4_r8 * (416.16_r8/(temp+120._r8)) *    &
                            ((temp/296.16_r8)**1.5_r8)

  end function air_dynamic_viscosity

  !=============================================================================
  real(r8) function air_kinematic_viscosity( temp, rhoair )
    !-----------------------------------------------------------------
    ! Calculate kinematic viscosity of air, unit [cm^2/s]
    !-----------------------------------------------------------------

    real(r8),intent(in) :: temp     ! air temperature [K]
    real(r8),intent(in) :: rhoair   ! air density [g/cm3]

    real(r8) :: vsc_dyn_atm  ! dynamic viscosity of air [g/cm/s]

    vsc_dyn_atm = air_dynamic_viscosity( temp)
    air_kinematic_viscosity = vsc_dyn_atm/rhoair

  end function air_kinematic_viscosity

  !=============================================================================

end module aero_model
