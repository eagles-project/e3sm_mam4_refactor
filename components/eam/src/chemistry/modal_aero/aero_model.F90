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
  subroutine aero_model_wetdep(dt, dlf, dlf2, cmfmc2, state,                    &!Intent-ins
       sh_e_ed_ratio, mu, md, du, eu, ed, dp, dsubcld, jt, maxg, ideep, lengath,&
       species_class,                                                           &
       cam_out,                                                                 & !Intent-inout
       pbuf,                                                                    & !Pointer
       ptend                                                                    ) !Intent-out

    use modal_aero_deposition, only: set_srf_wetdep
    use wetdep,                only: wetdepa_v2, wetdep_inputs_set, &
                                     wetdep_inputs_unset, wetdep_inputs_t
    use modal_aero_data
    use modal_aero_calcsize,   only: modal_aero_calcsize_sub
    use modal_aero_wateruptake,only: modal_aero_wateruptake_dr
    use modal_aero_convproc,   only: ma_convproc_intr
    use mo_constants,          only: pi

    ! args
    type(physics_state), intent(in)    :: state       ! Physics state variables
    real(r8),            intent(in)    :: dt          ! time step
    real(r8),            intent(in)    :: dlf(:,:)    ! shallow+deep convective detrainment [kg/kg/s]
    real(r8),            intent(in)    :: dlf2(:,:)   ! Shal conv cldwtr detrainment (kg/kg/s - grid avg)
    real(r8),            intent(in)    :: cmfmc2(pcols,pverp) ! Shal conv mass flux (kg/m2/s)
    real(r8),            intent(in)    :: sh_e_ed_ratio(pcols,pver)  ! shallow conv [ent/(ent+det)] ratio
                                                ! mu, md, ..., ideep, lengath are all deep conv variables
                                                ! *** AND ARE GATHERED ***
    real(r8),            intent(in)    :: mu(pcols,pver)   ! Updraft mass flux (positive)
    real(r8),            intent(in)    :: md(pcols,pver)   ! Downdraft mass flux (negative)
    real(r8),            intent(in)    :: du(pcols,pver)   ! Mass detrain rate from updraft
    real(r8),            intent(in)    :: eu(pcols,pver)   ! Mass entrain rate into updraft
    real(r8),            intent(in)    :: ed(pcols,pver)   ! Mass entrain rate into downdraft
    ! eu, ed, du are "d(massflux)/dp" and are all positive
    real(r8),            intent(in)    :: dp(pcols,pver)   ! Delta pressure between interfaces
    real(r8),            intent(in)    :: dsubcld(pcols)   ! Delta pressure from cloud base to sfc
    
    integer,             intent(in)    :: jt(pcols)         ! Index of cloud top for each column
    integer,             intent(in)    :: maxg(pcols)       ! Index of cloud top for each column
    integer,             intent(in)    :: ideep(pcols)      ! Gathering array
    integer,             intent(in)    :: lengath           ! Gathered min lon indices over which to operate
    integer,             intent(in)    :: species_class(:)
    
    type(cam_out_t),     intent(inout) :: cam_out     ! export state
    type(physics_buffer_desc), pointer :: pbuf(:)

    type(physics_ptend), intent(out)   :: ptend       ! indivdual parameterization tendencies



    ! local vars

    integer :: i
    integer :: jnv ! index for scavcoefnv 3rd dimension
    integer :: jnummaswtr  ! indicates current aerosol species type (0 = number, 1 = dry mass, 2 = water)
    integer, parameter :: jaeronumb=0, jaeromass=1, jaerowater=2
    integer :: k
    integer :: lchnk ! chunk identifier
    integer :: lphase ! index for interstitial / cloudborne aerosol
    integer :: lspec ! index for aerosol number / chem-mass / water-mass
    integer :: lspectype
    integer :: lcoardust, lcoarnacl ! indices for coarse mode dust and seasalt masses
    integer :: m, mtmp ! mode index
    integer :: mm, mmai, mmtoo ! tracer (q-array) index
    integer :: ncol ! number of atmospheric columns
    integer :: mam_prevap_resusp_optcc
    integer :: strt_loop, end_loop, stride_loop !loop indices for the lphase loop

    real(r8) :: iscavt(pcols, pver)
    real(r8) :: icscavt(pcols, pver)
    real(r8) :: isscavt(pcols, pver)
    real(r8) :: bcscavt(pcols, pver)
    real(r8) :: bsscavt(pcols, pver)
    real(r8) :: sol_factb, sol_facti
    real(r8) :: sol_factic(pcols,pver)
    real(r8) :: sol_factbi, sol_factii, sol_factiic

    real(r8) :: sflx(pcols) ! deposition flux

    real(r8) :: d1p_prevap_resusp, v1p_prevap_resusp
    real(r8) :: dqdt_tmp(pcols,pver)      ! temporary array to hold tendency for the "current" aerosol species
    real(r8) :: dqdt_sv(pcols,pver,pcnst) ! temporary array to hold tendency for all interstitial aerosol species
    real(r8) :: f_act_conv(pcols,pver) ! prescribed aerosol activation fraction for convective cloud ! rce 2010/05/01
    real(r8) :: f_act_conv_coarse(pcols,pver) ! similar but for coarse mode ! rce 2010/05/02
    real(r8) :: f_act_conv_coarse_dust, f_act_conv_coarse_nacl ! rce 2010/05/02
    real(r8) :: fracis_cw(pcols,pver)
    real(r8) :: hygro_sum_old(pcols,pver) ! before removal [sum of (mass*hydro/dens)]
    real(r8) :: hygro_sum_del(pcols,pver) ! removal change to [sum of (mass*hydro/dens)]
    real(r8) :: hygro_sum_old_ik, hygro_sum_new_ik
    real(r8) :: prec(pcols) ! precipitation rate
    real(r8) :: q_tmp(pcols,pver) ! temporary array to hold "most current" mixing ratio for 1 species
    real(r8) :: qqcw_tmp(pcols,pver) ! temporary array to hold qqcw ! rce 2010/05/01
    real(r8) :: scavcoefnv(pcols,pver,0:2) ! Dana and Hales coefficient (/mm) for
                                           ! cloud-borne num & vol (0),
                                           ! interstitial num (1), interstitial vol (2)
    real(r8) :: tmpa, tmpb
    real(r8) :: tmpdust, tmpnacl
    real(r8) :: water_old, water_new ! temporary old/new aerosol water mix-rat

    logical  :: isprx(pcols,pver) ! true if precipation
    logical :: do_lphase2

    real(r8) :: tmp_evapdp, tmp_evapsh  !RCE
    real(r8) :: tmp_precdp, tmp_precsh  !RCE
    real(r8) :: tmp_resudp, tmp_resush  !RCE
    real(r8) :: sflxec(pcols), sflxecdp(pcols)  ! deposition flux  !RCE
    real(r8) :: sflxic(pcols), sflxicdp(pcols)  ! deposition flux  !RCE
    real(r8) :: sflxbc(pcols), sflxbcdp(pcols)  ! deposition flux  !RCE
    real(r8) :: rcscavt(pcols, pver)  !RCE
    real(r8) :: rsscavt(pcols, pver)  !RCE
    real(r8) :: qqcw_in(pcols,pver), qqcw_sav(pcols,pver,0:maxd_aspectype)       ! temporary array to hold qqcw for the current mode  !RCE
    real(r8) :: rtscavt_sv(pcols, pver, pcnst) ! REASTER 08/12/2015
    real(r8) :: rcscavt_cn_sv(pcols, pver)     ! REASTER 08/12/2015
    real(r8) :: rsscavt_cn_sv(pcols, pver)     ! REASTER 08/12/2015
    
    real(r8), pointer :: fldcw(:,:)

    real(r8), pointer :: dgnumwet(:,:,:)
    real(r8), pointer :: qaerwat(:,:,:)  ! aerosol water
    real(r8), pointer :: rate1ord_cw2pr_st(:,:)

    real(r8), pointer :: fracis(:,:,:)   ! fraction of transported species that are insoluble

    integer, parameter:: nsrflx_mzaer2cnvpr = 2  !RCE 2012/01/12 bgn
    real(r8)          :: aerdepwetis(pcols,pcnst) ! aerosol wet deposition (interstitial) 
    real(r8)          :: aerdepwetcw(pcols,pcnst) ! aerosol wet deposition (cloud water)  
    real(r8)          :: qsrflx_mzaer2cnvpr(pcols,pcnst,nsrflx_mzaer2cnvpr)
    real(r8)          :: rprddpsum(pcols),  rprdshsum(pcols)   ! RCE 2012/01/12
    real(r8)          :: evapcdpsum(pcols), evapcshsum(pcols)  ! RCE 2012/01/12
    real(r8), pointer :: rprddp(:,:)     ! rain production, deep convection
    real(r8), pointer :: rprdsh(:,:)     ! rain production, deep convection
    real(r8), pointer :: evapcsh(:,:)    ! Evaporation rate of shallow convective precipitation >=0.
    real(r8), pointer :: evapcdp(:,:)    ! Evaporation rate of deep    convective precipitation >=0.

    real(r8), pointer :: icwmrdp(:,:)    ! in cloud water mixing ratio, deep convection
    real(r8), pointer :: icwmrsh(:,:)    ! in cloud water mixing ratio, deep convection
    real(r8), pointer :: sh_frac(:,:)    ! Shallow convective cloud fraction
    real(r8), pointer :: dp_frac(:,:)    ! Deep convective cloud fraction

    character(len=100) :: msg

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

    call pbuf_get_field(pbuf, dgnumwet_idx,       dgnumwet, start=(/1,1,1/), kount=(/pcols,pver,nmodes/) )
    call pbuf_get_field(pbuf, qaerwat_idx,        qaerwat,  start=(/1,1,1/), kount=(/pcols,pver,nmodes/) )
    call pbuf_get_field(pbuf, rate1_cw2pr_st_idx, rate1ord_cw2pr_st)
    call pbuf_get_field(pbuf, fracis_idx,         fracis, start=(/1,1,1/), kount=(/pcols, pver, pcnst/) )

    !Compute variables needed for convproc unified convective transport
    call pbuf_get_field(pbuf, rprddp_idx,      rprddp  )
    call pbuf_get_field(pbuf, rprdsh_idx,      rprdsh  )
    call pbuf_get_field(pbuf, nevapr_shcu_idx, evapcsh )
    call pbuf_get_field(pbuf, nevapr_dpcu_idx, evapcdp )
    evapcdpsum(:) = 0.0_r8
    rprddpsum(:)  = 0.0_r8  
    evapcshsum(:) = 0.0_r8
    rprdshsum(:)  = 0.0_r8
    do k = 1, pver
       rprddpsum(:ncol)  = rprddpsum(:ncol)  +  rprddp(:ncol,k)*state%pdel(:ncol,k)/gravit
       rprdshsum(:ncol)  = rprdshsum(:ncol)  +  rprdsh(:ncol,k)*state%pdel(:ncol,k)/gravit
       evapcdpsum(:ncol) = evapcdpsum(:ncol) + evapcdp(:ncol,k)*state%pdel(:ncol,k)/gravit
       evapcshsum(:ncol) = evapcshsum(:ncol) + evapcsh(:ncol,k)*state%pdel(:ncol,k)/gravit
    enddo  !RCE 2012/01/12 end


    prec(:ncol)=0._r8
    do k=1,pver
       where (prec(:ncol) >= 1.e-7_r8)
          isprx(:ncol,k) = .true.
       elsewhere
          isprx(:ncol,k) = .false.
       endwhere
       prec(:ncol) = prec(:ncol) + (dep_inputs%prain(:ncol,k) + dep_inputs%cmfdqr(:ncol,k) - dep_inputs%evapr(:ncol,k)) &
            *state%pdel(:ncol,k)/gravit
    enddo
    
    qsrflx_mzaer2cnvpr(:,:,:) = 0.0_r8  !RCE
    aerdepwetis(:,:)          = 0.0_r8  !RCE
    aerdepwetcw(:,:)          = 0.0_r8  !RCE
    qqcw_tmp(:,:)             = 0.0_r8  !RCE

    ! calculate the mass-weighted sol_factic for coarse mode species
    f_act_conv_coarse(:,:) = 0.60_r8 ! rce 2010/05/02
    f_act_conv_coarse_dust = 0.40_r8 ! rce 2010/05/02
    f_act_conv_coarse_nacl = 0.80_r8 ! rce 2010/05/02
    if (modeptr_coarse > 0) then
       lcoardust = lptr_dust_a_amode(modeptr_coarse)
       lcoarnacl = lptr_nacl_a_amode(modeptr_coarse)
       if ((lcoardust > 0) .and. (lcoarnacl > 0)) then
          do k = 1, pver
             do i = 1, ncol
                tmpdust = max( 0.0_r8, state%q(i,k,lcoardust) + ptend%q(i,k,lcoardust)*dt )
                tmpnacl = max( 0.0_r8, state%q(i,k,lcoarnacl) + ptend%q(i,k,lcoarnacl)*dt )
                if ((tmpdust+tmpnacl) > 1.0e-30_r8) then
                   f_act_conv_coarse(i,k) = (f_act_conv_coarse_dust*tmpdust &
                        + f_act_conv_coarse_nacl*tmpnacl)/(tmpdust+tmpnacl) ! rce 2010/05/02
                endif
             enddo
          enddo
       endif
    endif

    scavcoefnv(:,:,0) = 0.0_r8 ! below-cloud scavcoef = 0.0 for cloud-borne species

       ! resuspension goes to a different phase or mode
    rtscavt_sv(:,:,:) = 0.0_r8
    rcscavt_cn_sv(:,:) = 0.0_r8
    rsscavt_cn_sv(:,:) = 0.0_r8

    !BSINGH: Decide the loop counters for the lphase loop 
    !for cases with and without the unified convective transport
    !Counters for "without" unified convective treatment (i.e. default case)
    !BSINGH (09/12/2014):Do cloudborne first for unified convection scheme so
    !that the resuspension of cloudborne
    !can be saved then applied to interstitial (RCE)
    strt_loop   =  2
    end_loop    =  1
    stride_loop = -1

mmode_loop_aa: &
    do mtmp = 1, ntot_amode ! main loop over aerosol modes
       m = mtmp
       if (ntot_amode == 4) then
          ! for mam4, do accum, aitken, pcarbon, then coarse 
          if (mtmp == modeptr_coarse) then
             m = ntot_amode
          else if (mtmp > modeptr_coarse) then
             m = mtmp - 1
          endif
       endif
          
lphase_loop_aa: &
       do lphase = strt_loop,end_loop, stride_loop ! loop over interstitial (1) and cloud-borne (2) forms

          ! sol_factb and sol_facti values
          ! sol_factb - currently this is basically a tuning factor
          ! sol_facti & sol_factic - currently has a physical basis, and reflects activation fraction
          !
          ! 2008-mar-07 rce - sol_factb (interstitial) changed from 0.3 to 0.1
          ! - sol_factic (interstitial, dust modes) changed from 1.0 to 0.5
          ! - sol_factic (cloud-borne, pcarb modes) no need to set it to 0.0
          ! because the cloud-borne pcarbon == 0 (no activation)
          !
          ! rce 2010/05/02
          ! prior to this date, sol_factic was used for convective in-cloud wet removal,
          ! and its value reflected a combination of an activation fraction (which varied between modes)
          ! and a tuning factor
          ! from this date forward, two parameters are used for convective in-cloud wet removal
          ! f_act_conv is the activation fraction
          ! note that "non-activation" of aerosol in air entrained into updrafts should
          ! be included here
          ! eventually we might use the activate routine (with w ~= 1 m/s) to calculate
          ! this, but there is still the entrainment issue
          ! sol_factic is strictly a tuning factor
          !
          if (lphase == 1) then ! interstial aerosol
             hygro_sum_old(:,:) = 0.0_r8
             hygro_sum_del(:,:) = 0.0_r8
             call modal_aero_bcscavcoef_get( m, ncol, isprx, dgnumwet, &
                  scavcoefnv(:,:,1), scavcoefnv(:,:,2) )

             sol_facti = 0.0_r8 ! strat in-cloud scav totally OFF for institial
             sol_factic = 0.4_r8 ! xl 2010/05/20
             if (sscav_tuning) then
                sol_factb  = 0.03_r8   ! all below-cloud scav ON (0.1 "tuning factor")  ! tuned 1/6
             else
                sol_factb  = 0.1_r8    ! all below-cloud scav ON (0.1 "tuning factor")
             endif
             if (m == modeptr_pcarbon) then
                f_act_conv = 0.0_r8 ! rce 2010/05/02
             elseif ((m == modeptr_finedust) .or. (m == modeptr_coardust)) then
                f_act_conv = 0.4_r8 ! rce 2010/05/02
             else
                if (sscav_tuning) then
                   f_act_conv = 0.4_r8   ! rce 2010/05/02
                else
                   f_act_conv = 0.8_r8   ! rce 2010/05/02
                endif
             endif

          else ! cloud-borne aerosol (borne by stratiform cloud drops)

             sol_factb  = 0.0_r8   ! all below-cloud scav OFF (anything cloud-borne is located "in-cloud")
             if (sscav_tuning) then 
                sol_facti  = min(0.6_r8, sol_facti_cloud_borne)  ! strat  in-cloud scav totally ON for cloud-borne  ! tuned 1/6
             else
                sol_facti  = sol_facti_cloud_borne   ! strat  in-cloud scav cloud-borne tuning factor
             endif
             sol_factic = 0.0_r8   ! conv   in-cloud scav OFF (having this on would mean
                                   !        that conv precip collects strat droplets)
             f_act_conv = 0.0_r8   ! conv   in-cloud scav OFF (having this on would mean

          endif

          if( lphase == 1 ) then  ! FORTRAN refactor: can be merged above
             ! if modal aero convproc is turned on for aerosols, then
             !    turn off the convective in-cloud removal for interstitial aerosols
             !    (but leave the below-cloud on, as convproc only does in-cloud)
             !    and turn off the outfld SFWET, SFSIC, SFSID, SFSEC, and SFSED calls 
             ! for (stratiform)-cloudborne aerosols, convective wet removal
             !    (all forms) is zero, so no action is needed
             sol_factic = 0.0_r8
          endif
          !
          ! rce 2010/05/03
          ! wetdepa has 6 "sol_fact" parameters:
          ! sol_facti, sol_factic, sol_factb for liquid cloud
          ! sol_factii, sol_factiic, sol_factbi for ice cloud
          ! the ice cloud parameters are optional, and if not provided, they default to
          ! one of the other sol_fact parameters (see subr. wetdepa about this)
          ! for now, we set the ice cloud parameters equal
          ! to their liquid cloud counterparts
          ! currently the ice parameters are not used in wetdepa as
          ! wetdepa sets "weight" (the ice cloud fraction) to 0.0
          ! if this changes, we will have to give more thought to
          ! the ice cloud parameter values
          !
          sol_factbi = sol_factb
          sol_factii = sol_facti
          sol_factiic = sol_factic(1,1)


! REASTER 08/12/2015 - changed ordering (mass then number) for prevap resuspend to coarse
lspec_loop_aa: &
          do lspec = 1, nspec_amode(m)+2 ! loop over number + chem constituents + water

             mmai = 0
             if (lspec <= nspec_amode(m)) then ! non-water mass
                jnummaswtr = jaeromass
                if (lphase == 1) then
                   mm = lmassptr_amode(lspec,m)
                   jnv = 2
                else
                   mm = lmassptrcw_amode(lspec,m)
                   mmai = lmassptr_amode(lspec,m)
                   jnv = 0
                endif
             elseif (lspec == nspec_amode(m)+1) then ! number
                jnummaswtr = jaeronumb
                if (lphase == 1) then
                   mm = numptr_amode(m)
                   jnv = 1
                else
                   mm = numptrcw_amode(m)
                   mmai = numptr_amode(m)
                   jnv = 0
                endif
             else ! water mass
                ! bypass wet removal of aerosol water
                jnummaswtr = jaerowater
                cycle 
             endif

             if (mm <= 0) cycle



! mam_prevap_resusp_optcc values control the prevap_resusp calculations in wetdepa_v2:
!     0 = no resuspension
!   130 = non-linear resuspension of aerosol mass   based on scavenged aerosol mass
!   230 = non-linear resuspension of aerosol number based on raindrop number
!   the 130 thru 230 all use the new prevap_resusp code block in subr wetdepa_v2
!
             mam_prevap_resusp_optcc = 0

             if ( jnummaswtr == jaeromass ) then
                   mam_prevap_resusp_optcc = 130
             elseif ( jnummaswtr == jaeronumb .and. lphase == 1 .and. m == modeptr_coarse ) then
                   mam_prevap_resusp_optcc = 230
             endif


             ! set f_act_conv for interstitial (lphase=1) coarse mode species
             ! for the convective in-cloud, we conceptually treat the coarse dust and seasalt
             ! as being externally mixed, and apply f_act_conv = f_act_conv_coarse_dust/nacl to dust/seasalt
             ! number and sulfate are conceptually partitioned to the dust and seasalt
             ! on a mass basis, so the f_act_conv for number and sulfate are
             ! mass-weighted averages of the values used for dust/seasalt
             if ((lphase == 1) .and. (m == modeptr_coarse)) then
                f_act_conv = f_act_conv_coarse ! rce 2010/05/02
                if (jnummaswtr == jaeromass) then
                   if (lmassptr_amode(lspec,m) == lptr_dust_a_amode(m)) then
                      f_act_conv = f_act_conv_coarse_dust ! rce 2010/05/02
                   elseif (lmassptr_amode(lspec,m) == lptr_nacl_a_amode(m)) then
                      f_act_conv = f_act_conv_coarse_nacl ! rce 2010/05/02
                   endif
                endif
             endif


lphase_jnmw_conditional: &
             if ((lphase == 1) .and. (jnummaswtr /= jaerowater)) then
                ptend%lq(mm) = .true.
                dqdt_tmp(:,:) = 0.0_r8
                ! q_tmp reflects changes from modal_aero_calcsize and is the "most current" q
                q_tmp(1:ncol,:) = state%q(1:ncol,:,mm) + ptend%q(1:ncol,:,mm)*dt
                !Feed in the saved cloudborne mixing ratios from phase 2
                qqcw_in(:,:) = qqcw_sav(:,:,lspec)

                call wetdepa_v2( &
                     ncol, dt, state%pdel, &
                     dep_inputs%cmfdqr, dep_inputs%evapc, dlf, dep_inputs%conicw, &
                     dep_inputs%prain, dep_inputs%evapr, dep_inputs%totcond, &
                     dep_inputs%cldt, dep_inputs%cldcu, &
                     dep_inputs%cldvcu, dep_inputs%cldvst, &
                     sol_factb, sol_facti, sol_factic, &
                     mam_prevap_resusp_optcc, .false., scavcoefnv(:,:,jnv), f_act_conv, &
                     q_tmp, qqcw_in(:,:), &
                     fracis(:,:,mm), dqdt_tmp, iscavt, &
                     icscavt, isscavt, bcscavt, bsscavt, rcscavt, rsscavt )

                ! resuspension goes to coarse mode
                ! first deduct the current resuspension from the dqdt_tmp of the current species
                dqdt_tmp(1:ncol,:) = dqdt_tmp(1:ncol,:) - ( rcscavt(1:ncol,:) + rsscavt(1:ncol,:) )
                ! then add the current resuspension to the rtscavt_sv of the appropriate coarse mode species
                mmtoo = mmtoo_prevap_resusp(mm)
                if (mmtoo > 0) rtscavt_sv(1:ncol,:,mmtoo) = rtscavt_sv(1:ncol,:,mmtoo) & 
                                  + ( rcscavt(1:ncol,:) + rsscavt(1:ncol,:) )
                ! then add the rtscavt_sv of the current species to the dqdt_tmp of the current species
                ! note that for so4_a3 and mam3, the rtscavt_sv at this point will have resuspension contributions
                !    from so4_a1/2/3 and so4c1/2/3
                dqdt_tmp(1:ncol,:) = dqdt_tmp(1:ncol,:) + rtscavt_sv(1:ncol,:,mm)

                ptend%q(1:ncol,:,mm) = ptend%q(1:ncol,:,mm) + dqdt_tmp(1:ncol,:)

                call outfld( trim(cnst_name(mm))//'WET', dqdt_tmp(:,:), pcols, lchnk)
                call outfld( trim(cnst_name(mm))//'SIC', icscavt, pcols, lchnk)
                call outfld( trim(cnst_name(mm))//'SIS', isscavt, pcols, lchnk)
                call outfld( trim(cnst_name(mm))//'SBC', bcscavt, pcols, lchnk)
                call outfld( trim(cnst_name(mm))//'SBS', bsscavt, pcols, lchnk)

                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+dqdt_tmp(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                aerdepwetis(:ncol,mm) = sflx(:ncol)

                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+icscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                sflxic = sflx

                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+isscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                call outfld( trim(cnst_name(mm))//'SFSIS', sflx, pcols, lchnk)

                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+bcscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                call outfld( trim(cnst_name(mm))//'SFSBC', sflx, pcols, lchnk)
                sflxbc = sflx

                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+bsscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                call outfld( trim(cnst_name(mm))//'SFSBS', sflx, pcols, lchnk)
               
                ! here the prevap resuspension is in rcscavt & rsscavt and column integral is written to history
                !BSINGH(09/15/2014):Following two nested do-loops are new additions for unified convection 
                !BSINGH(09/15/2014):After these do-loops, code was added by RCE, the comments by RCE are kept as it is
                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+rcscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                sflxec = sflx
                   
                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+rsscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                call outfld( trim(cnst_name(mm))//'SFSES', sflx, pcols, lchnk)                   
                   
                ! apportion convective surface fluxes to deep and shallow conv
                ! this could be done more accurately in subr wetdepa
                ! since deep and shallow rarely occur simultaneously, and these
                !    fields are just diagnostics, this approximate method is adequate
                ! only do this for interstitial aerosol, because conv clouds to not
                !    affect the stratiform-cloudborne aerosol
                do i = 1, ncol
                  tmp_precdp = max( rprddpsum(i),  1.0e-35_r8 )
                  tmp_precsh = max( rprdshsum(i),  1.0e-35_r8 )
                  tmp_evapdp = max( evapcdpsum(i), 0.1e-35_r8 )
                  tmp_evapsh = max( evapcshsum(i), 0.1e-35_r8 )
                      
                  ! assume that in- and below-cloud removal are proportional to column precip production
                  tmpa = tmp_precdp / (tmp_precdp + tmp_precsh)
                  tmpa = max( 0.0_r8, min( 1.0_r8, tmpa ) )
                  sflxicdp(i) = sflxic(i)*tmpa
                  sflxbcdp(i) = sflxbc(i)*tmpa
                      
                  ! assume that resuspension is proportional to (wet removal)*[(precip evap)/(precip production)]
                  tmp_resudp =           tmpa  * min( (tmp_evapdp/tmp_precdp), 1.0_r8 )
                  tmp_resush = (1.0_r8 - tmpa) * min( (tmp_evapsh/tmp_precsh), 1.0_r8 )
                  tmpb = max( tmp_resudp, 1.0e-35_r8 ) / max( (tmp_resudp+tmp_resush), 1.0e-35_r8 )
                  tmpb = max( 0.0_r8, min( 1.0_r8, tmpb ) )
                  sflxecdp(i) = sflxec(i)*tmpb
                enddo
                call outfld( trim(cnst_name(mm))//'SFSBD', sflxbcdp, pcols, lchnk)
                ! when ma_convproc_intr is used, convective in-cloud wet removal is done there
                ! the convective (total and deep) precip-evap-resuspension includes in- and below-cloud
                ! contributions, so pass the below-cloud contribution to ma_convproc_intr
                qsrflx_mzaer2cnvpr(1:ncol,mm,1) = sflxec(  1:ncol)
                qsrflx_mzaer2cnvpr(1:ncol,mm,2) = sflxecdp(1:ncol)

             elseif (lphase == 2) then lphase_jnmw_conditional ! lphase == 2
                do_lphase2 = .true.
! There is no cloud-borne aerosol water in the model, so the do_lphase2 code block
! should NEVER execute for lspec = nspec_amode(m)+1 (i.e., jnummaswtr = jaerowater).
! The code only worked because the "do lspec" loop cycles when lspec = nspec_amode(m)+1,
! but that does not make the code correct.
                if (jnummaswtr == jaerowater) do_lphase2 = .false.
do_lphase2_conditional: &
                if ( do_lphase2 ) then
                   dqdt_tmp(:,:) = 0.0_r8
                   fldcw => qqcw_get_field(pbuf,mm,lchnk)
                   qqcw_sav(1:ncol,:,lspec) = fldcw(1:ncol,:)  !RCE 2012/01/12
                   
                call wetdepa_v2( &
                     ncol, dt, state%pdel, &
                     dep_inputs%cmfdqr, dep_inputs%evapc, dlf, dep_inputs%conicw, &
                     dep_inputs%prain, dep_inputs%evapr, dep_inputs%totcond, &
                     dep_inputs%cldt, dep_inputs%cldcu, &
                     dep_inputs%cldvcu, dep_inputs%cldvst, &
                     sol_factb, sol_facti, sol_factic, &
                     mam_prevap_resusp_optcc, .true., scavcoefnv(:,:,jnv), f_act_conv, &
                     fldcw, qqcw_tmp,  & ! FIXIT: is it a bug using qqcw_tmp here? Different with the previous call
                     fracis_cw, dqdt_tmp, iscavt, &
                     icscavt, isscavt, bcscavt, bsscavt, rcscavt, rsscavt ) 

                   ! resuspension goes to coarse mode
                   ! first deduct the current resuspension from the dqdt_tmp of the current species
                   dqdt_tmp(1:ncol,:) = dqdt_tmp(1:ncol,:) - ( rcscavt(1:ncol,:) + rsscavt(1:ncol,:) )
                   ! then add the current resuspension to the rtscavt_sv of the appropriate coarse mode species
                   mmtoo = mmtoo_prevap_resusp(mm)
                   if (mmtoo > 0) rtscavt_sv(1:ncol,:,mmtoo) = rtscavt_sv(1:ncol,:,mmtoo) & 
                                  + ( rcscavt(1:ncol,:) + rsscavt(1:ncol,:) )
                   
                   fldcw(1:ncol,:) = fldcw(1:ncol,:) + dqdt_tmp(1:ncol,:) * dt

                   sflx(:)=0._r8
                   do k=1,pver
                      do i=1,ncol
                         sflx(i)=sflx(i)+dqdt_tmp(i,k)*state%pdel(i,k)/gravit
                      enddo
                   enddo
                   call outfld( trim(cnst_name_cw(mm))//'SFWET', sflx, pcols, lchnk)
                   aerdepwetcw(:ncol,mm) = sflx(:ncol)
                   
                   sflx(:)=0._r8
                   do k=1,pver
                      do i=1,ncol
                         sflx(i)=sflx(i)+icscavt(i,k)*state%pdel(i,k)/gravit
                      enddo
                   enddo
                   call outfld( trim(cnst_name_cw(mm))//'SFSIC', sflx, pcols, lchnk)

                   sflx(:)=0._r8
                   do k=1,pver
                      do i=1,ncol
                         sflx(i)=sflx(i)+isscavt(i,k)*state%pdel(i,k)/gravit
                      enddo
                   enddo
                   call outfld( trim(cnst_name_cw(mm))//'SFSIS', sflx, pcols, lchnk)

                   sflx(:)=0._r8
                   do k=1,pver
                      do i=1,ncol
                         sflx(i)=sflx(i)+bcscavt(i,k)*state%pdel(i,k)/gravit
                      enddo
                   enddo
                   call outfld( trim(cnst_name_cw(mm))//'SFSBC', sflx, pcols, lchnk)

                   sflx(:)=0._r8
                   do k=1,pver
                      do i=1,ncol
                         sflx(i)=sflx(i)+bsscavt(i,k)*state%pdel(i,k)/gravit
                      enddo
                   enddo
                   call outfld( trim(cnst_name_cw(mm))//'SFSBS', sflx, pcols, lchnk)

                   sflx(:)=0.0_r8
                   do k=1,pver
                      sflx(1:ncol)=sflx(1:ncol)+rcscavt(1:ncol,k)*state%pdel(1:ncol,k)/gravit
                   enddo
                   call outfld( trim(cnst_name_cw(mm))//'SFSEC', sflx, pcols, lchnk)
                      
                   sflx(:)=0.0_r8
                   do k=1,pver
                      sflx(1:ncol)=sflx(1:ncol)+rsscavt(1:ncol,k)*state%pdel(1:ncol,k)/gravit
                   enddo
                   call outfld( trim(cnst_name_cw(mm))//'SFSES', sflx, pcols, lchnk)

                endif do_lphase2_conditional

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
    integer :: lunerr           ! logical unit for error message
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
    
    lunerr = 6

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
               scavratenum, scavratevol, lunerr )

          scavimptblnum(jgrow,imode) = log( scavratenum )
          scavimptblvol(jgrow,imode) = log( scavratevol )

       enddo growloop
    enddo modeloop
    return

  end subroutine modal_aero_bcscavcoef_init

  !===============================================================================
  subroutine modal_aero_bcscavcoef_get( m, ncol, isprx, dgn_awet, scavcoefnum, scavcoefvol )

    use modal_aero_data, only: dgnum_amode
    !-----------------------------------------------------------------------
    implicit none

    integer,intent(in) :: m, ncol
    logical,intent(in):: isprx(pcols,pver)
    real(r8), intent(in) :: dgn_awet(pcols,pver,ntot_amode)
    real(r8), intent(out) :: scavcoefnum(pcols,pver), scavcoefvol(pcols,pver)

    integer i, k, jgrow
    real(r8) dumdgratio, xgrow, dumfhi, dumflo, scavimpvol, scavimpnum


    do k = 1, pver
       do i = 1, ncol

          ! do only if no precip
          if ( isprx(i,k) ) then
             !
             ! interpolate table values using log of (actual-wet-size)/(base-dry-size)

             dumdgratio = dgn_awet(i,k,m)/dgnum_amode(m)

             if ((dumdgratio .ge. 0.99_r8) .and. (dumdgratio .le. 1.01_r8)) then
                scavimpvol = scavimptblvol(0,m)
                scavimpnum = scavimptblnum(0,m)
             else
                xgrow = log( dumdgratio ) / dlndg_nimptblgrow
                jgrow = int( xgrow )
                if (xgrow .lt. 0._r8) jgrow = jgrow - 1
                if (jgrow .lt. nimptblgrow_mind) then
                   jgrow = nimptblgrow_mind
                   xgrow = jgrow
                else
                   jgrow = min( jgrow, nimptblgrow_maxd-1 )
                endif

                dumfhi = xgrow - jgrow
                dumflo = 1._r8 - dumfhi

                scavimpvol = dumflo*scavimptblvol(jgrow,m) + &
                     dumfhi*scavimptblvol(jgrow+1,m)
                scavimpnum = dumflo*scavimptblnum(jgrow,m) + &
                     dumfhi*scavimptblnum(jgrow+1,m)

             endif

             ! impaction scavenging removal amount for volume
             scavcoefvol(i,k) = exp( scavimpvol )
             ! impaction scavenging removal amount to number
             scavcoefnum(i,k) = exp( scavimpnum )

             ! scavcoef = impaction scav rate (1/h) for precip = 1 mm/h
             ! scavcoef = impaction scav rate (1/s) for precip = pfx_inrain
             ! (scavcoef/3600) = impaction scav rate (1/s) for precip = 1 mm/h
             ! (pfx_inrain*3600) = in-rain-area precip rate (mm/h)
             ! impactrate = (scavcoef/3600) * (pfx_inrain*3600)
          else
             scavcoefvol(i,k) = 0._r8
             scavcoefnum(i,k) = 0._r8
          endif

       enddo
    enddo

    return
  end subroutine modal_aero_bcscavcoef_get

  !===============================================================================
	subroutine calc_1_impact_rate(             &
     		dg0, sigmag, rhoaero, temp, press, &
     		scavratenum, scavratevol, lunerr )
   !
   !   routine computes a single impaction scavenging rate
   !	for precipitation rate of 1 mm/h
   !
   !   dg0 = geometric mean diameter of aerosol number size distrib. (cm)
   !   sigmag = geometric standard deviation of size distrib.
   !   rhoaero = density of aerosol particles (g/cm^3)
   !   temp = temperature (K)
   !   press = pressure (dyne/cm^2)
   !   scavratenum = number scavenging rate (1/h)
   !   scavratevol = volume or mass scavenging rate (1/h)
   !   lunerr = logical unit for error message
   !
   use mo_constants, only: boltz_cgs, pi, rhowater => rhoh2o_cgs, &
                           gravity => gravity_cgs, rgas => rgas_cgs

   implicit none

   !   subr. parameters
   integer lunerr
   real(r8) dg0, sigmag, rhoaero, temp, press, scavratenum, scavratevol

   !   local variables
   integer nrainsvmax
   parameter (nrainsvmax=50)
   real(r8) rrainsv(nrainsvmax), xnumrainsv(nrainsvmax),&
        vfallrainsv(nrainsvmax)

   integer naerosvmax
   parameter (naerosvmax=51)
   real(r8) aaerosv(naerosvmax), &
     	ynumaerosv(naerosvmax), yvolaerosv(naerosvmax)

   integer i, ja, jr, na, nr
   real(r8) a, aerodiffus, aeromass, ag0, airdynvisc, airkinvisc
   real(r8) anumsum, avolsum, cair, chi
   real(r8) d, dr, dum, dumfuchs, dx
   real(r8) ebrown, eimpact, eintercept, etotal, freepath
   real(r8) precip, precipmmhr, precipsum
   real(r8) r, rainsweepout, reynolds, rhi, rhoair, rlo, rnumsum
   real(r8) scavsumnum, scavsumnumbb
   real(r8) scavsumvol, scavsumvolbb
   real(r8) schmidt, sqrtreynolds, sstar, stokes, sx              
   real(r8) taurelax, vfall, vfallstp
   real(r8) x, xg0, xg3, xhi, xlo, xmuwaterair                     

   
   rlo = .005_r8
   rhi = .250_r8
   dr = 0.005_r8
   nr = 1 + nint( (rhi-rlo)/dr )
   if (nr .gt. nrainsvmax) then
      write(lunerr,9110)
      call endrun()
   endif
9110 format( '*** subr. calc_1_impact_rate -- nr > nrainsvmax' )

   precipmmhr = 1.0_r8
   precip = precipmmhr/36000._r8

   ag0 = dg0/2._r8
   sx = log( sigmag )
   xg0 = log( ag0 )
   xg3 = xg0 + 3._r8*sx*sx

   xlo = xg3 - 4._r8*sx
   xhi = xg3 + 4._r8*sx
   dx = 0.2_r8*sx

   dx = max( 0.2_r8*sx, 0.01_r8 )
   xlo = xg3 - max( 4._r8*sx, 2._r8*dx )
   xhi = xg3 + max( 4._r8*sx, 2._r8*dx )

   na = 1 + nint( (xhi-xlo)/dx )
   if (na .gt. naerosvmax) then
      write(lunerr,9120)
      call endrun()
   endif
9120 format( '*** subr. calc_1_impact_rate -- na > naerosvmax' )

   !   air molar density
   cair = press/(rgas*temp)
   !   air mass density
   rhoair = 28.966_r8*cair
   !   molecular freepath
   freepath = 2.8052e-10_r8/cair
   !   air dynamic viscosity
   airdynvisc = 1.8325e-4_r8 * (416.16_r8/(temp+120._r8)) *    &
        ((temp/296.16_r8)**1.5_r8)
   !   air kinemaic viscosity
   airkinvisc = airdynvisc/rhoair
   !   ratio of water viscosity to air viscosity (from Slinn)
   xmuwaterair = 60.0_r8

   !
   !   compute rain drop number concentrations
   !	rrainsv = raindrop radius (cm)
   !	xnumrainsv = raindrop number concentration (#/cm^3)
   !		(number in the bin, not number density)
   !	vfallrainsv = fall velocity (cm/s)
   !
   precipsum = 0._r8
   do i = 1, nr
      r = rlo + (i-1)*dr
      rrainsv(i) = r
      xnumrainsv(i) = exp( -r/2.7e-2_r8 )

      d = 2._r8*r
      if (d .le. 0.007_r8) then
         vfallstp = 2.88e5_r8 * d**2._r8
      elseif (d .le. 0.025_r8) then
         vfallstp = 2.8008e4_r8 * d**1.528_r8
      elseif (d .le. 0.1_r8) then
         vfallstp = 4104.9_r8 * d**1.008_r8
      elseif (d .le. 0.25_r8) then
         vfallstp = 1812.1_r8 * d**0.638_r8
      else
         vfallstp = 1069.8_r8 * d**0.235_r8
      endif

      vfall = vfallstp * sqrt(1.204e-3_r8/rhoair)
      vfallrainsv(i) = vfall
      precipsum = precipsum + vfall*(r**3)*xnumrainsv(i)
   enddo
   precipsum = precipsum*pi*1.333333_r8

   rnumsum = 0._r8
   do i = 1, nr
      xnumrainsv(i) = xnumrainsv(i)*(precip/precipsum)
      rnumsum = rnumsum + xnumrainsv(i)
   enddo

   !
   !   compute aerosol concentrations
   !	aaerosv = particle radius (cm)
   !	fnumaerosv = fraction of total number in the bin (--)
   !	fvolaerosv = fraction of total volume in the bin (--)
   !
   anumsum = 0._r8
   avolsum = 0._r8
   do i = 1, na
      x = xlo + (i-1)*dx
      a = exp( x )
      aaerosv(i) = a
      dum = (x - xg0)/sx
      ynumaerosv(i) = exp( -0.5_r8*dum*dum )
      yvolaerosv(i) = ynumaerosv(i)*1.3333_r8*pi*a*a*a
      anumsum = anumsum + ynumaerosv(i)
      avolsum = avolsum + yvolaerosv(i)
   enddo

   do i = 1, na
      ynumaerosv(i) = ynumaerosv(i)/anumsum
      yvolaerosv(i) = yvolaerosv(i)/avolsum
   enddo


   !
   !   compute scavenging
   !
   scavsumnum = 0._r8
   scavsumvol = 0._r8
   !
   !   outer loop for rain drop radius
   !
   jr_loop: do jr = 1, nr

      r = rrainsv(jr)
      vfall = vfallrainsv(jr)

      reynolds = r * vfall / airkinvisc
      sqrtreynolds = sqrt( reynolds )

      !
      !   inner loop for aerosol particle radius
      !
      scavsumnumbb = 0._r8
      scavsumvolbb = 0._r8

      ja_loop: do ja = 1, na

         a = aaerosv(ja)

         chi = a/r

         dum = freepath/a
         dumfuchs = 1._r8 + 1.246_r8*dum + 0.42_r8*dum*exp(-0.87_r8/dum)
         taurelax = 2._r8*rhoaero*a*a*dumfuchs/(9._r8*rhoair*airkinvisc)

         aeromass = 4._r8*pi*a*a*a*rhoaero/3._r8
         aerodiffus = boltz_cgs*temp*taurelax/aeromass

         schmidt = airkinvisc/aerodiffus
         stokes = vfall*taurelax/r

         ebrown = 4._r8*(1._r8 + 0.4_r8*sqrtreynolds*(schmidt**0.3333333_r8)) /  &
              (reynolds*schmidt)

         dum = (1._r8 + 2._r8*xmuwaterair*chi) /         &
              (1._r8 + xmuwaterair/sqrtreynolds)
         eintercept = 4._r8*chi*(chi + dum)

         dum = log( 1._r8 + reynolds )
         sstar = (1.2_r8 + dum/12._r8) / (1._r8 + dum)
         eimpact = 0._r8
         if (stokes .gt. sstar) then
	    dum = stokes - sstar
	    eimpact = (dum/(dum+0.6666667_r8)) ** 1.5_r8
         endif

         etotal = ebrown + eintercept + eimpact
         etotal = min( etotal, 1.0_r8 )

         rainsweepout = xnumrainsv(jr)*4._r8*pi*r*r*vfall

         scavsumnumbb = scavsumnumbb + rainsweepout*etotal*ynumaerosv(ja)
         scavsumvolbb = scavsumvolbb + rainsweepout*etotal*yvolaerosv(ja)

      enddo ja_loop

      scavsumnum = scavsumnum + scavsumnumbb
      scavsumvol = scavsumvol + scavsumvolbb

   enddo jr_loop

   scavratenum = scavsumnum*3600._r8
   scavratevol = scavsumvol*3600._r8

   return
 end subroutine calc_1_impact_rate
  
  !=============================================================================
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

end module aero_model
