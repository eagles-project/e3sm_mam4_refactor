
module modal_aer_opt

  ! parameterizes aerosol coefficients using chebychev polynomial
  ! parameterize aerosol radiative properties in terms of
  ! surface mode wet radius and wet refractive index

  ! Ghan and Zaveri, JGR 2007.

  ! uses Wiscombe's (1979) mie scattering code


  use shr_kind_mod,      only: r8 => shr_kind_r8, shr_kind_cl
  use ppgrid,            only: pcols, pver
  use spmd_utils,        only: masterproc
  use phys_control,      only: cam_chempkg_is
  use ref_pres,          only: top_lev => clim_modal_aero_top_lev
  use physconst,         only: rhoh2o, rga, rair
  use radconstants,      only: nswbands, nlwbands, idx_sw_diag, idx_uv_diag, idx_nir_diag
  use rad_constituents,  only: n_diag, rad_cnst_get_call_list, rad_cnst_get_info, &
       rad_cnst_get_mode_props

  use pio,               only: file_desc_t, var_desc_t, pio_inq_dimlen, pio_inq_dimid, pio_inq_varid, &
       pio_get_var, pio_nowrite, pio_closefile
  use cam_pio_utils,     only: cam_pio_openfile
  use cam_history,       only:  addfld, horiz_only, add_default, outfld
  use cam_history_support, only: fillvalue
  use cam_logfile,       only: iulog
  use cam_abortutils,        only: endrun

  use modal_aero_data,  only: ntot_amode, nspec_amode, specdens_amode, &
       specname_amode, spechygro, &
       sigmag_amode, lspectype_amode, lmassptr_amode, &
       specrefndxlw, refrtablw, refitablw, absplw, &
       specrefndxsw, refrtabsw, refitabsw, extpsw, abspsw, asmpsw

  use modal_aero_wateruptake, only: modal_aero_wateruptake_dr
  use modal_aero_calcsize,    only: modal_aero_calcsize_sub
  use shr_log_mod ,           only: errmsg => shr_log_errmsg
  use mam_support,            only: ptr2d_t

  implicit none
  private
  save

  public :: modal_aer_opt_readnl, modal_aer_opt_init, modal_aero_sw, modal_aero_lw


  character(len=*), parameter :: unset_str = 'UNSET'

  ! Namelist variables:
  character(shr_kind_cl)      :: modal_optics_file = unset_str   ! full pathname for modal optics dataset
  character(shr_kind_cl)      :: water_refindex_file = unset_str ! full pathname for water refractive index dataset

  ! Dimension sizes in coefficient arrays used to parameterize aerosol radiative properties
  ! in terms of refractive index and wet radius
  integer, parameter :: ncoef=5, prefr=7, prefi=10

  ! refractive index for water read in read_water_refindex
  complex(r8) :: crefwsw(nswbands) ! complex refractive index for water visible
  complex(r8) :: crefwlw(nlwbands) ! complex refractive index for water infrared

  character(len=4) :: diag(0:n_diag) = (/'    ','_d1 ','_d2 ','_d3 ','_d4 ','_d5 ', &
       '_d6 ','_d7 ','_d8 ','_d9 ','_d10'/)

  !Declare the following threadprivate variables to be used for calcsize and water uptake
  !These are defined as module level variables to aviod allocation-deallocation in a loop
  real(r8), allocatable, target :: dgnumdry_m(:,:,:) ! number mode dry diameter for all modes
  real(r8), allocatable, target :: dgnumwet_m(:,:,:) ! number mode wet diameter for all modes
  real(r8), allocatable, target :: qaerwat_m(:,:,:)  ! aerosol water (g/g) for all modes
  !$OMP THREADPRIVATE(dgnumdry_m, dgnumwet_m, qaerwat_m)

  ! small values treated as zero
  real(r8), parameter :: small_value_40 = 1.e-40_r8

  ! min, max aerosol surface mode radius treated [m]
  real(r8), parameter :: rmmin = 0.01e-6_r8
  real(r8), parameter :: rmmax = 25.e-6_r8
  real(r8), parameter :: xrmin = log(rmmin)
  real(r8), parameter :: xrmax = log(rmmax)

  !===============================================================================
CONTAINS
  !===============================================================================

  subroutine modal_aer_opt_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'modal_aer_opt_readnl'

    namelist /modal_aer_opt_nl/ water_refindex_file
    !-----------------------------------------------------------------------------

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'modal_aer_opt_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, modal_aer_opt_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    call mpibcast(water_refindex_file, len(water_refindex_file), mpichar, 0, mpicom)
#endif


  end subroutine modal_aer_opt_readnl

  !===============================================================================

  subroutine modal_aer_opt_init()

    use ioFileMod,        only: getfil
    use phys_control,     only: phys_getopts

    ! Local variables

    integer  :: i, m
    character(len=256) :: locfile

    logical           :: history_amwg            ! output the variables used by the AMWG diag package
    logical           :: history_verbose         ! produce verbose history output
    logical           :: history_aero_optics     ! output aerosol optics diagnostics

    logical :: call_list(0:n_diag)
    integer :: ilist, nmodes, m_ncoef, m_prefr, m_prefi
    integer :: errcode, istat

    character(len=*), parameter :: routine='modal_aer_opt_init'
    !----------------------------------------------------------------------------

    ! set flags to check aerosol settings
    call rad_cnst_get_info(0, nmodes=nmodes)

    ! Check that dimension sizes in the coefficient arrays used to
    ! parameterize aerosol radiative properties are consistent between this
    ! module and the mode physprop files.
    call rad_cnst_get_call_list(call_list)
    do ilist = 0, n_diag
       if (call_list(ilist)) then
          call rad_cnst_get_info(ilist, nmodes=nmodes)
          do m = 1, nmodes
             call rad_cnst_get_mode_props(ilist, m, ncoef=m_ncoef, prefr=m_prefr, prefi=m_prefi)
             if (m_ncoef /= ncoef .or. m_prefr /= prefr .or. m_prefi /= prefi) then
                write(iulog,*) routine//': ERROR - file and module values do not match:'
                write(iulog,*) '   ncoef:', ncoef, m_ncoef
                write(iulog,*) '   prefr:', prefr, m_prefr
                write(iulog,*) '   prefi:', prefi, m_prefi
                call endrun(routine//': ERROR - file and module values do not match')
             end if
          end do
       end if
    end do

    call getfil(water_refindex_file, locfile)
    call read_water_refindex(locfile)
    if (masterproc) write(iulog,*) "modal_aer_opt_init: read water refractive index file:", trim(locfile)

    call phys_getopts(history_amwg_out        = history_amwg, &
         history_verbose_out = history_verbose, &
         history_aero_optics_out = history_aero_optics )


    !obtain nmodes for the climate (ilist = 0) list
    call rad_cnst_get_info(0, nmodes=nmodes)

    !Allocate dry and wet size variables in an OMP PARRALLEL region as these
    !arrays are private for each thread and needs to be allocated for each OMP thread
    !$OMP PARALLEL
    allocate(dgnumdry_m(pcols,pver,nmodes),stat=istat)
    if (istat .ne. 0) call endrun("Unable to allocate dgnumdry_m: "//errmsg(__FILE__,__LINE__) )
    allocate(dgnumwet_m(pcols,pver,nmodes),stat=istat)
    if (istat .ne. 0) call endrun("Unable to allocate dgnumwet_m: "//errmsg(__FILE__,__LINE__) )
    allocate(qaerwat_m(pcols,pver,nmodes),stat=istat)
    if (istat .ne. 0) call endrun("Unable to allocate qaerwat_m: "//errmsg(__FILE__,__LINE__) )
    !$OMP END PARALLEL

    ! Add diagnostic fields to history output.

    call addfld ('EXTINCT',(/ 'lev' /),    'A','/m','Aerosol extinction', flag_xyfill=.true.)
    call addfld ('tropopause_m',horiz_only,    'A',' m  ','tropopause level in meters', flag_xyfill=.true.)
    call addfld ('ABSORB',(/ 'lev' /),    'A','/m','Aerosol absorption', flag_xyfill=.true.)
    call addfld ('AODVIS',horiz_only,    'A','  ','Aerosol optical depth 550 nm', flag_xyfill=.true., &
         standard_name='atmosphere_optical_thickness_due_to_ambient_aerosol_particles')
    call addfld ('AODALL',horiz_only,    'A','  ','AOD 550 nm for all time and species', flag_xyfill=.true.)
    call addfld ('AODUV',horiz_only,    'A','  ','Aerosol optical depth 350 nm', flag_xyfill=.true.)
    call addfld ('AODNIR',horiz_only,    'A','  ','Aerosol optical depth 850 nm', flag_xyfill=.true.)
    call addfld ('AODABS',horiz_only,    'A','  ','Aerosol absorption optical depth 550 nm', flag_xyfill=.true., &
         standard_name='atmosphere_absorption_optical_thickness_due_to_ambient_aerosol_particles')
    call addfld ('AODMODE1',horiz_only,    'A','  ','Aerosol optical depth 550 nm mode 1'           , flag_xyfill=.true.)
    call addfld ('AODMODE2',horiz_only,    'A','  ','Aerosol optical depth 550 nm mode 2'           , flag_xyfill=.true.)
    call addfld ('AODMODE3',horiz_only,    'A','  ','Aerosol optical depth 550 nm mode 3'           , flag_xyfill=.true.)
    call addfld ('AODDUST1',horiz_only,    'A','  ','Aerosol optical depth 550 nm model 1 from dust', flag_xyfill=.true.)
    call addfld ('AODDUST2',horiz_only,    'A','  ','Aerosol optical depth 550 nm model 2 from dust', flag_xyfill=.true.)
    call addfld ('AODDUST3',horiz_only,    'A','  ','Aerosol optical depth 550 nm model 3 from dust', flag_xyfill=.true.)
    call addfld ('AODDUST',horiz_only,    'A','  ','Aerosol optical depth 550 nm from dust', flag_xyfill=.true.)
    call addfld ('AODSO4',horiz_only,    'A','  ','Aerosol optical depth 550 nm from SO4', flag_xyfill=.true.)
    call addfld ('AODPOM',horiz_only,    'A','  ','Aerosol optical depth 550 nm from POM', flag_xyfill=.true.)
    call addfld ('AODSOA',horiz_only,    'A','  ','Aerosol optical depth 550 nm from SOA', flag_xyfill=.true.)
    call addfld ('AODBC',horiz_only,    'A','  ','Aerosol optical depth 550 nm from BC', flag_xyfill=.true.)
    call addfld ('AODSS',horiz_only,    'A','  ','Aerosol optical depth 550 nm from seasalt', flag_xyfill=.true.)
    call addfld ('AODABSBC',horiz_only, 'A','  ','Aerosol absorption optical depth 550 nm from BC', flag_xyfill=.true.)
    call addfld ('AODMOM',horiz_only,    'A','  ','Aerosol optical depth 550 nm from marine organic', flag_xyfill=.true.)
    call addfld ('BURDEN1',horiz_only,  'A','kg/m2'      ,'Aerosol burden mode 1'      , flag_xyfill=.true.)
    call addfld ('BURDEN2',horiz_only,  'A','kg/m2'      ,'Aerosol burden mode 2'      , flag_xyfill=.true.)
    call addfld ('BURDEN3',horiz_only,  'A','kg/m2'      ,'Aerosol burden mode 3'      , flag_xyfill=.true.)
    call addfld ('BURDENDUST',horiz_only,  'A','kg/m2'   ,'Dust aerosol burden'        , flag_xyfill=.true.)
    call addfld ('BURDENSO4',horiz_only,  'A','kg/m2'    ,'Sulfate aerosol burden'     , flag_xyfill=.true.)
    call addfld ('BURDENPOM',horiz_only,  'A','kg/m2'    ,'POM aerosol burden'         , flag_xyfill=.true.)
    call addfld ('BURDENSOA',horiz_only,  'A','kg/m2'    ,'SOA aerosol burden'         , flag_xyfill=.true.)
    call addfld ('BURDENBC',horiz_only,  'A','kg/m2'     ,'Black carbon aerosol burden', flag_xyfill=.true.)
    call addfld ('BURDENSEASALT',horiz_only,  'A','kg/m2','Seasalt aerosol burden'     , flag_xyfill=.true.)
    call addfld ('BURDENMOM',horiz_only,  'A','kg/m2'    ,'Marine organic aerosol burden', flag_xyfill=.true.)
    call addfld ('SSAVIS',horiz_only,    'A','  ','Aerosol singel-scatter albedo', flag_xyfill=.true.)


    if (history_amwg) then
       call add_default ('AODDUST1'     , 1, ' ')
       call add_default ('AODDUST3'     , 1, ' ')
       call add_default ('AODVIS'       , 1, ' ')
       call add_default ('AODALL'       , 1, ' ')
       call add_default ('BURDEN1'      , 1, ' ')
       call add_default ('BURDEN2'      , 1, ' ')
       call add_default ('BURDEN3'      , 1, ' ')
       if ( history_verbose ) then
          call add_default ('BURDENDUST'   , 1, ' ')
          call add_default ('BURDENSO4'    , 1, ' ')
          call add_default ('BURDENPOM'    , 1, ' ')
          call add_default ('BURDENSOA'    , 1, ' ')
          call add_default ('BURDENBC'     , 1, ' ')
          call add_default ('BURDENSEASALT', 1, ' ')
          call add_default ('BURDENMOM', 1, ' ')
       end if
    end if

    if (history_aero_optics) then
       call add_default ('AODDUST1'     , 1, ' ')
       call add_default ('AODDUST3'     , 1, ' ')
       call add_default ('AODMODE1'     , 1, ' ')
       call add_default ('AODMODE2'     , 1, ' ')
       call add_default ('AODMODE3'     , 1, ' ')
       call add_default ('AODVIS'       , 1, ' ')
       call add_default ('AODALL'       , 1, ' ')
       call add_default ('AODUV'        , 1, ' ')
       call add_default ('AODNIR'       , 1, ' ')
       call add_default ('AODABS'       , 1, ' ')
       call add_default ('AODABSBC'     , 1, ' ')
       call add_default ('AODDUST'      , 1, ' ')
       call add_default ('AODSO4'       , 1, ' ')
       call add_default ('AODPOM'       , 1, ' ')
       call add_default ('AODSOA'       , 1, ' ')
       call add_default ('AODBC'        , 1, ' ')
       call add_default ('AODSS'        , 1, ' ')
       call add_default ('BURDEN1'      , 1, ' ')
       call add_default ('BURDEN2'      , 1, ' ')
       call add_default ('BURDEN3'      , 1, ' ')
       if (history_verbose) then
          call add_default ('ABSORB'       , 1, ' ')
          call add_default ('BURDENDUST'   , 1, ' ')
          call add_default ('BURDENSO4'    , 1, ' ')
          call add_default ('BURDENPOM'    , 1, ' ')
          call add_default ('BURDENSOA'    , 1, ' ')
          call add_default ('BURDENBC'     , 1, ' ')
          call add_default ('BURDENSEASALT', 1, ' ')
          call add_default ('BURDENMOM'    , 1, ' ')
       end if
       call add_default ('SSAVIS'       , 1, ' ')
       call add_default ('EXTINCT'      , 1, ' ')
    end if
    if (cam_chempkg_is('trop_mam4').or.cam_chempkg_is('trop_mam4_resus').or. &
         cam_chempkg_is('trop_mam4_mom').or.cam_chempkg_is('trop_mam4_resus_mom').or. &
         cam_chempkg_is('trop_mam4_resus_soag').or.cam_chempkg_is('trop_mam7').or. &
         cam_chempkg_is('trop_mam9').or.cam_chempkg_is('trop_strat_mam7').or. &
         cam_chempkg_is('linoz_mam4_resus').or.cam_chempkg_is('linoz_mam4_resus_soag').or.&
         cam_chempkg_is('linoz_mam4_resus_mom').or. &
         cam_chempkg_is('linoz_mam4_resus_mom_soag').or.cam_chempkg_is('superfast_mam4_resus_mom_soag')) then
       call addfld ('AODDUST4',horiz_only,    'A','  ','Aerosol optical depth 550 nm model 4 from dust', flag_xyfill=.true.)
       call addfld ('AODMODE4',horiz_only,    'A','  ','Aerosol optical depth 550 nm mode 4', flag_xyfill=.true.)
       call addfld ('BURDEN4',horiz_only,    'A','kg/m2','Aerosol burden mode 4', flag_xyfill=.true.)


       if (history_aero_optics) then
          call add_default ('AODDUST4', 1, ' ')
          call add_default ('AODMODE4', 1, ' ')
          call add_default ('BURDEN4' , 1, ' ')
       end if
    end if
    if (cam_chempkg_is('trop_mam7').or.cam_chempkg_is('trop_mam9').or.cam_chempkg_is('trop_strat_mam7')) then
       call addfld ('AODDUST5',horiz_only,    'A','  ','Aerosol optical depth 550 nm model 5 from dust', flag_xyfill=.true.)
       call addfld ('AODDUST6',horiz_only,    'A','  ','Aerosol optical depth 550 nm model 6 from dust', flag_xyfill=.true.)
       call addfld ('AODDUST7',horiz_only,    'A','  ','Aerosol optical depth 550 nm model 7 from dust', flag_xyfill=.true.)
       call addfld ('AODMODE5',horiz_only,    'A','  ','Aerosol optical depth 550 nm mode 5', flag_xyfill=.true.)
       call addfld ('AODMODE6',horiz_only,    'A','  ','Aerosol optical depth 550 nm mode 6', flag_xyfill=.true.)
       call addfld ('AODMODE7',horiz_only,    'A','  ','Aerosol optical depth 550 nm mode 7', flag_xyfill=.true.)
       call addfld ('BURDEN5',horiz_only,    'A','kg/m2','Aerosol burden mode 5', flag_xyfill=.true.)
       call addfld ('BURDEN6',horiz_only,    'A','kg/m2','Aerosol burden mode 6', flag_xyfill=.true.)
       call addfld ('BURDEN7',horiz_only,    'A','kg/m2','Aerosol burden mode 7', flag_xyfill=.true.)

       if (history_aero_optics) then
          call add_default ('AODDUST5', 1, ' ')
          call add_default ('AODDUST6', 1, ' ')
          call add_default ('AODDUST7', 1, ' ')
          call add_default ('AODMODE5', 1, ' ')
          call add_default ('AODMODE6', 1, ' ')
          call add_default ('AODMODE7', 1, ' ')
          if (history_verbose) then
             call add_default ('BURDEN5', 1, ' ')
             call add_default ('BURDEN6', 1, ' ')
             call add_default ('BURDEN7', 1, ' ')
          end if
       end if
    end if
    if (cam_chempkg_is('trop_mam9')) then
       call addfld ('AODMODE8',horiz_only,    'A','  '  ,'Aerosol optical depth 550 nm mode 8', flag_xyfill=.true.)
       call addfld ('AODMODE9',horiz_only,    'A','  '  ,'Aerosol optical depth 550 nm mode 9', flag_xyfill=.true.)
       call addfld ('BURDEN8',horiz_only,    'A','kg/m2','Aerosol burden mode 8', flag_xyfill=.true.)
       call addfld ('BURDEN9',horiz_only,    'A','kg/m2','Aerosol burden mode 9', flag_xyfill=.true.)
       if (history_aero_optics) then
          call add_default ('AODMODE8', 1, ' ')
          call add_default ('AODMODE9', 1, ' ')
          call add_default ('BURDEN8', 1, ' ')
          call add_default ('BURDEN9', 1, ' ')
       end if
    end if

    do ilist = 1, n_diag
       if (call_list(ilist)) then

          call addfld ('EXTINCT'//diag(ilist), (/ 'lev' /), 'A','1/m', &
               'Aerosol extinction', flag_xyfill=.true.)
          call addfld ('ABSORB'//diag(ilist),  (/ 'lev' /), 'A','1/m', &
               'Aerosol absorption', flag_xyfill=.true.)
          call addfld ('AODVIS'//diag(ilist),       horiz_only, 'A','  ', &
               'Aerosol optical depth 550 nm', flag_xyfill=.true.)
          call addfld ('AODALL'//diag(ilist),       horiz_only, 'A','  ', &
               'AOD 550 nm all time', flag_xyfill=.true.)
          call addfld ('AODABS'//diag(ilist),       horiz_only, 'A','  ', &
               'Aerosol absorption optical depth 550 nm', flag_xyfill=.true.)

          call add_default ('EXTINCT'//diag(ilist), 1, ' ')
          call add_default ('ABSORB'//diag(ilist),  1, ' ')
          call add_default ('AODVIS'//diag(ilist),  1, ' ')
          call add_default ('AODALL'//diag(ilist),  1, ' ')
          call add_default ('AODABS'//diag(ilist),  1, ' ')

       end if
    end do

  end subroutine modal_aer_opt_init

  !===============================================================================
  subroutine modal_aero_sw(dt, lchnk, ncol, state_q, state_zm, temperature, pmid, pdel, &  ! in
       pdeldry, cldn, nnite, idxnite, is_cmip6_volc, &  ! in
       ext_cmip6_sw, trop_level,  & ! in
       qqcw, tauxar, wa, ga, fa) ! out
    ! calculates aerosol sw radiative properties

    use mam_support, only : min_max_bound

    real(r8),         intent(in) :: dt               !timestep [s]
    integer,          intent(in) :: lchnk            ! chunk id
    integer,          intent(in) :: ncol             ! number of active columns in the chunk
    real(r8), target, intent(in) :: state_q(:,:,:)
    real(r8),         intent(in) :: state_zm(:,:)
    real(r8),         intent(in) :: temperature(:,:)
    real(r8),         intent(in) :: pmid(:,:)
    real(r8),         intent(in) :: pdel(:,:)
    real(r8),         intent(in) :: pdeldry(:,:)

    real(r8), target, intent(in) :: cldn(:,:)         ! layer cloud fraction [fraction]

    integer,          intent(in) :: nnite          ! number of night columns
    integer,          intent(in) :: idxnite(nnite) ! local column indices of night columns
    integer,          intent(in) :: trop_level(pcols)!tropopause level for each column
    real(r8),         intent(in) :: ext_cmip6_sw(pcols,pver) ! aerosol shortwave extinction [1/m]
    logical,          intent(in) :: is_cmip6_volc

    type(ptr2d_t),    intent(inout) :: qqcw(:)               ! Cloud borne aerosols mixing ratios [kg/kg or 1/kg]
    real(r8),         intent(out) :: tauxar(pcols,0:pver,nswbands) ! layer extinction optical depth [1]
    real(r8),         intent(out) :: wa(pcols,0:pver,nswbands)     ! layer single-scatter albedo [1]
    real(r8),         intent(out) :: ga(pcols,0:pver,nswbands)     ! asymmetry factor [1]
    real(r8),         intent(out) :: fa(pcols,0:pver,nswbands)     ! forward scattered fraction [1]

    ! Local variables
    integer :: icol, jj, isw, kk, ll, mm    ! indices                   
    integer :: nspec
    integer :: istat
    integer :: itim_old           ! index

    real(r8) :: mass(pcols,pver)        ! layer mass [kg]
    real(r8) :: air_density(pcols,pver) ! [kg/m3]

    real(r8),    pointer :: specmmr(:,:)        ! species mass mixing ratio [kg/kg]
    character*32         :: spectype            ! species type
    real(r8)             :: hygro_aer           ! hygroscopicity [1]

    real(r8) :: sigma_logr_aer         ! geometric standard deviation of number distribution
    real(r8) :: radsurf(pcols,pver)    ! aerosol surface mode radius
    real(r8) :: logradsurf(pcols,pver) ! log(aerosol surface mode radius)
    real(r8) :: cheb(ncoef,pcols,pver) ! chebychev polynomial parameters

    real(r8),allocatable :: specvol(:,:)        ! volume concentration of aerosol specie [m3/kg]
    real(r8),allocatable :: specdens(:)         ! species density for all species [kg/m3]
    complex(r8),allocatable :: specrefindex(:,:)     ! species refractive index

    real(r8)    :: refr(pcols)     ! real part of refractive index
    real(r8)    :: refi(pcols)     ! imaginary part of refractive index
    complex(r8) :: crefin(pcols)   ! complex refractive index

    real(r8) :: dryvol(pcols)   ! volume concentration of aerosol mode [m3/kg]
    real(r8) :: watervol(pcols) ! volume concentration of water in each mode [m3/kg]
    real(r8) :: wetvol(pcols)   ! volume concentration of wet mode [m3/kg]

    integer  :: itab(pcols), jtab(pcols)
    real(r8) :: ttab(pcols), utab(pcols)
    real(r8) :: cext(pcols,ncoef), cabs(pcols,ncoef), casm(pcols,ncoef)
    real(r8) :: pext(pcols)     ! parameterized specific extinction [m2/kg]
    real(r8) :: specpext(pcols) ! specific extinction [m2/kg]
    real(r8) :: dopaer(pcols)   ! aerosol optical depth in layer
    real(r8) :: pabs(pcols)     ! parameterized specific absorption [m2/kg]
    real(r8) :: pasm(pcols)     ! parameterized asymmetry factor
    real(r8) :: palb(pcols)     ! parameterized single scattering albedo

    ! Diagnostics
    real(r8) :: tropopause_m(pcols)
    real(r8) :: extinct(pcols,pver)         ! aerosol extinction [1/m]
    real(r8) :: absorb(pcols,pver)          ! aerosol absorption [1/m]
    real(r8) :: aodvis(pcols)               ! extinction optical depth
    real(r8) :: aodall(pcols)               ! extinction optical depth
    real(r8) :: aodabs(pcols)               ! absorption optical depth

    real(r8) :: aodabsbc(pcols)             ! absorption optical depth of BC

    real(r8) :: burdenmode(pcols)           ! aerosol burden for each mode
    real(r8) :: aodmode(pcols)

    real(r8) :: dustaodmode(pcols)          ! dust aod in aerosol mode
    real(r8) :: dustvol(pcols)              ! volume concentration of dust in aerosol mode (m3/kg)

    real(r8) :: burdendust(pcols), burdenso4(pcols), burdenbc(pcols), &
         burdenpom(pcols), burdensoa(pcols), burdenseasalt(pcols), burdenmom(pcols)
    real(r8) :: scatdust(pcols), scatso4(pcols), scatbc(pcols), &
         scatpom(pcols), scatsoa(pcols), scatseasalt(pcols), scatmom(pcols)
    real(r8) :: absdust(pcols), absso4(pcols), absbc(pcols), &
         abspom(pcols), abssoa(pcols), absseasalt(pcols), absmom(pcols)
    real(r8) :: hygrodust(pcols), hygroso4(pcols), hygrobc(pcols), &
         hygropom(pcols), hygrosoa(pcols), hygroseasalt(pcols), hygromom(pcols)

    real(r8) :: ssavis(pcols)
    real(r8) :: specrefr, specrefi
    real(r8) :: scath2o, absh2o, sumscat, sumabs, sumhygro

    ! total species AOD
    real(r8) :: dustaod(pcols), so4aod(pcols), bcaod(pcols), &
         pomaod(pcols), soaaod(pcols), seasaltaod(pcols), momaod(pcols)

    logical :: savaervis ! true if visible wavelength (0.55 micron)
    logical :: savaernir ! true if near ir wavelength (~0.88 micron)
    logical :: savaeruv  ! true if uv wavelength (~0.35 micron)

    real(r8) :: aoduv(pcols)               ! extinction optical depth in uv
    real(r8) :: aodnir(pcols)              ! extinction optical depth in nir

    character(len=32) :: outname

    ! debug output
    integer  :: nerr_dopaer = 0
    !----------------------------------------------------------------------------

    mass(:ncol,:)        = pdeldry(:ncol,:)*rga
    air_density(:ncol,:) = pmid(:ncol,:)/(rair*temperature(:ncol,:))

    ! initialize output variables
    tauxar(:ncol,:,:) = 0._r8
    wa(:ncol,:,:)     = 0._r8
    ga(:ncol,:,:)     = 0._r8
    fa(:ncol,:,:)     = 0._r8
    ! zero'th layer does not contain aerosol
    tauxar(1:ncol,0,:)  = 0._r8
    wa(1:ncol,0,:)      = 0.925_r8
    ga(1:ncol,0,:)      = 0.850_r8
    fa(1:ncol,0,:)      = 0.7225_r8

    ! diagnostics for visible band summed over modes
    extinct(:ncol,:)      = 0.0_r8
    absorb(:ncol,:)       = 0.0_r8
    aodvis(:ncol)         = 0.0_r8
    aodall(:ncol)         = 0.0_r8
    aodabs(:ncol)         = 0.0_r8
    burdendust(:ncol)     = 0.0_r8
    burdenso4(:ncol)      = 0.0_r8
    burdenpom(:ncol)      = 0.0_r8
    burdensoa(:ncol)      = 0.0_r8
    burdenbc(:ncol)       = 0.0_r8
    burdenseasalt(:ncol)  = 0.0_r8
    burdenmom(:ncol)      = 0.0_r8
    momaod(:ncol)         = 0.0_r8
    ssavis(:ncol)         = 0.0_r8

    aodabsbc(:ncol)       = 0.0_r8
    dustaod(:ncol)        = 0.0_r8
    so4aod(:ncol)         = 0.0_r8
    pomaod(:ncol)         = 0.0_r8
    soaaod(:ncol)         = 0.0_r8
    bcaod(:ncol)          = 0.0_r8
    seasaltaod(:ncol)     = 0.0_r8

    ! diags for other bands
    aoduv(:ncol)          = 0.0_r8
    aodnir(:ncol)         = 0.0_r8

    ! Calculate aerosol size distribution parameters and aerosol water uptake
    !For prognostic aerosols
    call modal_aero_calcsize_sub(ncol, lchnk, state_q, pdel, dt, qqcw, list_idx_in=0, update_mmr_in = .false., & ! in
         dgnumdry_m=dgnumdry_m) ! out

    call modal_aero_wateruptake_dr(lchnk, ncol, state_q, temperature, pmid, & ! in 
         cldn, dgnumdry_m, & ! in
         dgnumwet_m, qaerwat_m, & ! inout
         list_idx_in=0   ) ! optional in


    ! loop over all aerosol modes
    do mm = 1, ntot_amode

       ! diagnostics for visible band for each mode
       burdenmode(1:ncol)  = 0.0_r8
       aodmode(1:ncol)     = 0.0_r8
       dustaodmode(1:ncol) = 0.0_r8

       ! get mode info
       nspec = nspec_amode(mm)
       sigma_logr_aer = sigmag_amode(mm)

       ! calc size parameter for all columns
       call modal_size_parameters(ncol, sigma_logr_aer, dgnumwet_m(:,:,mm), & ! in
            radsurf,  logradsurf, cheb      ) ! out


       allocate(specvol(pcols,nspec),stat=istat)
       if (istat /= 0) call endrun("Unable to allocate specvol: "//errmsg(__FILE__,__LINE__) )
       allocate(specdens(nspec),stat=istat)
       if (istat /= 0) call endrun("Unable to allocate specdens: "//errmsg(__FILE__,__LINE__) )
       allocate(specrefindex(nspec,nswbands),stat=istat)
       if (istat /= 0) call endrun("Unable to allocate specrefindex: "//errmsg(__FILE__,__LINE__) )

       do isw = 1, nswbands
          savaervis = (isw == idx_sw_diag)
          savaeruv  = (isw == idx_uv_diag)
          savaernir = (isw == idx_nir_diag)

          do kk = top_lev, pver

             dustvol(:ncol) = 0._r8

             scatdust(:ncol)     = 0._r8
             absdust(:ncol)      = 0._r8
             hygrodust(:ncol)    = 0._r8
             scatso4(:ncol)      = 0._r8
             absso4(:ncol)       = 0._r8
             hygroso4(:ncol)     = 0._r8
             scatbc(:ncol)       = 0._r8
             absbc(:ncol)        = 0._r8
             hygrobc(:ncol)      = 0._r8
             scatpom(:ncol)      = 0._r8
             abspom(:ncol)       = 0._r8
             hygropom(:ncol)     = 0._r8
             scatsoa(:ncol)      = 0._r8
             abssoa(:ncol)       = 0._r8
             hygrosoa(:ncol)     = 0._r8
             scatseasalt(:ncol)  = 0._r8
             absseasalt(:ncol)   = 0._r8
             hygroseasalt(:ncol) = 0._r8
             scatmom(:ncol)      = 0._r8
             absmom(:ncol)       = 0._r8
             hygromom(:ncol)     = 0._r8

             ! aerosol species loop
             do ll = 1, nspec

                ! get aerosol properties and save for each species
                specmmr => state_q(:,:,lmassptr_amode(ll,mm))
                spectype = specname_amode(lspectype_amode(ll,mm))
                hygro_aer = spechygro(lspectype_amode(ll,mm))
                specdens(ll) = specdens_amode(lspectype_amode(ll,mm))
                specrefindex(ll,:) = specrefndxsw(:,lspectype_amode(ll,mm))
                specvol(:,ll) = specmmr(:,kk)/specdens(ll)

                ! compute some diagnostics for visible band only
                if (savaervis) then

                   specrefr = real(specrefindex(ll,isw))
                   specrefi = aimag(specrefindex(ll,isw))

                   do icol = 1, ncol
                      burdenmode(icol) = burdenmode(icol) + specmmr(icol,kk)*mass(icol,kk)
                   enddo

                   if (trim(spectype) == 'dust') then
                      call calc_diag_spec ( ncol, specmmr(:,kk), mass(:,kk), & ! in
                           specvol(:,ll), specrefr, specrefi, hygro_aer, & ! in
                           burdendust, scatdust, absdust, hygrodust) ! out
                      dustvol(:)    = specvol(:,ll)
                   endif
                   if (trim(spectype) == 'sulfate') then
                      call calc_diag_spec ( ncol, specmmr(:,kk), mass(:,kk), & ! in
                           specvol(:,ll), specrefr, specrefi, hygro_aer, & ! in
                           burdenso4, scatso4, absso4, hygroso4) ! out
                   endif
                   if (trim(spectype) == 'black-c') then
                      call calc_diag_spec ( ncol, specmmr(:,kk), mass(:,kk), & ! in
                           specvol(:,ll), specrefr, specrefi, hygro_aer, & ! in
                           burdenbc, scatbc, absbc, hygrobc) ! out
                   endif
                   if (trim(spectype) == 'p-organic') then
                      call calc_diag_spec ( ncol, specmmr(:,kk), mass(:,kk), & ! in
                           specvol(:,ll), specrefr, specrefi, hygro_aer, & ! in
                           burdenpom, scatpom, abspom, hygropom) ! out
                   endif
                   if (trim(spectype) == 's-organic') then
                      call calc_diag_spec ( ncol, specmmr(:,kk), mass(:,kk), & ! in
                           specvol(:,ll), specrefr, specrefi, hygro_aer, & ! in
                           burdensoa, scatsoa, abssoa, hygrosoa) ! out
                   endif
                   if (trim(spectype) == 'seasalt') then
                      call calc_diag_spec ( ncol, specmmr(:,kk), mass(:,kk), & ! in
                           specvol(:,ll), specrefr, specrefi, hygro_aer, & ! in
                           burdenseasalt, scatseasalt, absseasalt, hygroseasalt) ! out
                   endif
                   if (trim(spectype) == 'm-organic') then
                      call calc_diag_spec ( ncol, specmmr(:,kk), mass(:,kk), & ! in
                           specvol(:,ll), specrefr, specrefi, hygro_aer, & ! in
                           burdenmom, scatmom, absmom, hygromom) ! out
                   endif

                endif ! if (savaervis)
             enddo ! species loop ll

             call calc_refin_complex ('sw', ncol, isw,          & ! in
                  qaerwat_m(:,kk,mm), specvol, specrefindex,  & ! in
                  dryvol, wetvol, watervol, crefin, refr, refi) ! out

             ! interpolate coefficients linear in refractive index
             ! first call calcs itab,jtab,ttab,utab
             itab(:ncol) = 0
             call binterp(extpsw(mm,:,:,:,isw), ncol,                   & !in
                  refr, refi, refrtabsw(mm,:,isw), refitabsw(mm,:,isw), & !in
                  itab, jtab, ttab, utab, cext)                           !inout/out

             call binterp(abspsw(mm,:,:,:,isw), ncol,                   & !in
                  refr, refi, refrtabsw(mm,:,isw), refitabsw(mm,:,isw), & !in
                  itab, jtab, ttab, utab, cabs)                           !inout/out

             call binterp(asmpsw(mm,:,:,:,isw), ncol,                   & !in
                  refr, refi, refrtabsw(mm,:,isw), refitabsw(mm,:,isw), & !in
                  itab, jtab, ttab, utab, casm)                           !inout/out

             ! parameterized optical properties
             call calc_parameterized (ncol, ncoef, cext, cheb(:,:,kk), & ! in
                  pext) ! out
             call calc_parameterized (ncol, ncoef, cabs, cheb(:,:,kk), & ! in
                  pabs) ! out
             call calc_parameterized (ncol, ncoef, casm, cheb(:,:,kk), & ! in
                  pasm) ! out

             do icol=1,ncol

                if (logradsurf(icol,kk) <= xrmax) then
                   pext(icol) = exp(pext(icol))
                else
                   pext(icol) = 1.5_r8/(radsurf(icol,kk)*rhoh2o) ! geometric optics
                endif
                ! convert from m2/kg water to m2/kg aerosol
                specpext(icol) = pext(icol)
                pext(icol) = pext(icol)*wetvol(icol)*rhoh2o
                pabs(icol) = pabs(icol)*wetvol(icol)*rhoh2o

                pabs(icol) = min_max_bound(0._r8, pext(icol), pabs(icol))
                palb(icol) = 1._r8-pabs(icol)/max(pext(icol), small_value_40)
                dopaer(icol) = pext(icol)*mass(icol,kk)

                if (savaeruv) then
                   aoduv(icol) = aoduv(icol) + dopaer(icol)
                endif
                if (savaernir) then
                   aodnir(icol) = aodnir(icol) + dopaer(icol)
                endif

             enddo

             ! Save aerosol optical depth at longest visible wavelength
             ! sum over layers
             if (savaervis) then
                ! aerosol extinction (/m)
                do icol = 1, ncol
                   extinct(icol,kk) = extinct(icol,kk) + dopaer(icol)*air_density(icol,kk)/mass(icol,kk)
                   absorb(icol,kk)  = absorb(icol,kk) + pabs(icol)*air_density(icol,kk)
                   aodvis(icol )    = aodvis(icol) + dopaer(icol)
                   aodall(icol )    = aodall(icol) + dopaer(icol)
                   aodabs(icol )    = aodabs(icol) + pabs(icol)*mass(icol,kk)
                   aodmode(icol)    = aodmode(icol) + dopaer(icol)
                   ssavis(icol )    = ssavis(icol) + dopaer(icol)*palb(icol)

                   if (wetvol(icol) > small_value_40) then

                      dustaodmode(icol) = dustaodmode(icol) + dopaer(icol)*dustvol(icol)/wetvol(icol)

                      ! partition optical depth into contributions from each constituent
                      ! assume contribution is proportional to refractive index X volume

                      scath2o  = watervol(icol)*real(crefwsw(isw))
                      absh2o   = -watervol(icol)*aimag(crefwsw(isw))
                      sumscat  = scatso4(icol) + scatpom(icol) + scatsoa(icol) + scatbc(icol) + &
                           scatdust(icol) + scatseasalt(icol) + scath2o + &
                           scatmom(icol)
                      sumabs   = absso4(icol) + abspom(icol) + abssoa(icol) + absbc(icol) + &
                           absdust(icol) + absseasalt(icol) + absh2o + &
                           absmom(icol)
                      sumhygro = hygroso4(icol) + hygropom(icol) + hygrosoa(icol) + hygrobc(icol) + &
                           hygrodust(icol) + hygroseasalt(icol) + &
                           hygromom(icol)

                      call update_aod_spec ( &
                           scath2o,absh2o, sumhygro, sumscat, sumabs, & ! in
                           hygrodust(icol), palb(icol), dopaer(icol), & ! in
                           scatdust(icol), absdust(icol), dustaod(icol) ) ! inout

                      call update_aod_spec ( &
                           scath2o,absh2o, sumhygro,sumscat, sumabs, & ! in
                           hygroso4(icol), palb(icol), dopaer(icol), & ! in
                           scatso4(icol), absso4(icol), so4aod(icol) ) ! inout

                      call update_aod_spec ( &
                           scath2o,absh2o, sumhygro,sumscat, sumabs, & ! in
                           hygropom(icol), palb(icol), dopaer(icol), & ! in
                           scatpom(icol), abspom(icol), pomaod(icol) ) ! inout

                      call update_aod_spec ( &
                           scath2o,absh2o, sumhygro,sumscat, sumabs, & ! in
                           hygrosoa(icol), palb(icol), dopaer(icol), & ! in
                           scatsoa(icol), abssoa(icol), soaaod(icol) ) ! inout

                      call update_aod_spec ( &
                           scath2o,absh2o,sumhygro, sumscat,sumabs, & ! in
                           hygrobc(icol), palb(icol), dopaer(icol), & ! in
                           scatbc(icol), absbc(icol), bcaod(icol) ) ! inout

                      call update_aod_spec ( &
                           scath2o, absh2o,   sumhygro, sumscat, sumabs, & ! in
                           hygroseasalt(icol), palb(icol), dopaer(icol), & ! in
                           scatseasalt(icol), absseasalt(icol), seasaltaod(icol) ) ! inout

                      call update_aod_spec ( & 
                           scath2o,absh2o, sumhygro,sumscat, sumabs, & ! in
                           hygromom(icol), palb(icol), dopaer(icol), & ! in
                           scatmom(icol), absmom(icol), momaod(icol) ) ! inout

                      aodabsbc(icol) = aodabsbc(icol) + absbc(icol)*dopaer(icol)*(1.0_r8-palb(icol))

                   endif

                enddo ! icol
             endif

             do icol = 1, ncol
                call check_error_warning('sw', icol, kk, mm, isw, nspec, & ! in
                     dopaer(icol), pabs(icol), dryvol, wetvol, watervol, crefin,cabs,& ! in
                     specdens, specrefindex, specvol, & ! in
                     nerr_dopaer, & ! inout
                     pext(icol), specpext(icol) ) ! optional in

                tauxar(icol,kk,isw) = tauxar(icol,kk,isw) + dopaer(icol)
                wa(icol,kk,isw)     = wa(icol,kk,isw)     + dopaer(icol)*palb(icol)
                ga(icol,kk,isw)     = ga(icol,kk,isw)     + dopaer(icol)*palb(icol)*pasm(icol)
                fa(icol,kk,isw)     = fa(icol,kk,isw)     + dopaer(icol)*palb(icol)*pasm(icol)*pasm(icol)
             enddo

          enddo ! pver

       enddo ! sw bands

       ! mode diagnostics
       ! The diagnostics are currently only output for the climate list.  Code mods will
       ! be necessary to provide output for the rad_diag lists.

       do jj = 1, nnite
          aodmode(idxnite(jj)) = fillvalue
          dustaodmode(idxnite(jj)) = fillvalue
       enddo

       write(outname,'(a,i1)') 'BURDEN', mm
       call outfld(trim(outname), burdenmode, pcols, lchnk)

       write(outname,'(a,i1)') 'AODMODE', mm
       call outfld(trim(outname), aodmode, pcols, lchnk)

       write(outname,'(a,i1)') 'AODDUST', mm
       call outfld(trim(outname), dustaodmode, pcols, lchnk)


       deallocate(specvol)
       deallocate(specdens)
       deallocate(specrefindex)

    enddo ! nmodes

    !Add contributions from volcanic aerosols directly read in extinction
    if(is_cmip6_volc) then
       call calc_volc_ext(ncol, trop_level, state_zm, ext_cmip6_sw, & ! in
            extinct, tropopause_m ) ! inout/out
    endif


    ! Output visible band diagnostics for quantities summed over the modes
    ! These fields are put out for diagnostic lists as well as the climate list.
    do jj = 1, nnite
       extinct(idxnite(jj),:) = fillvalue
       absorb(idxnite(jj),:)  = fillvalue
       aodvis(idxnite(jj))    = fillvalue
       aodabs(idxnite(jj))    = fillvalue
    enddo

    call outfld('EXTINCT'//diag(0),  extinct, pcols, lchnk)
    call outfld('tropopause_m', tropopause_m, pcols, lchnk)
    call outfld('ABSORB'//diag(0),   absorb,  pcols, lchnk)
    call outfld('AODVIS'//diag(0),   aodvis,  pcols, lchnk)
    call outfld('AODALL'//diag(0),   aodall,  pcols, lchnk)
    call outfld('AODABS'//diag(0),   aodabs,  pcols, lchnk)

    ! These diagnostics are output only for climate list
    do icol = 1, ncol
       if (aodvis(icol) > 1.e-10_r8) then
          ssavis(icol) = ssavis(icol)/aodvis(icol)
       else
          ssavis(icol) = 0.925_r8
       endif
    enddo

    do jj = 1, nnite
       ssavis(idxnite(jj))     = fillvalue
       aoduv(idxnite(jj))      = fillvalue
       aodnir(idxnite(jj))     = fillvalue
       aodabsbc(idxnite(jj))   = fillvalue
       dustaod(idxnite(jj))    = fillvalue
       so4aod(idxnite(jj))     = fillvalue
       pomaod(idxnite(jj))     = fillvalue
       soaaod(idxnite(jj))     = fillvalue
       bcaod(idxnite(jj))      = fillvalue
       momaod(idxnite(jj))     = fillvalue
    enddo

    call outfld('SSAVIS',        ssavis,        pcols, lchnk)
    call outfld('AODUV',         aoduv,         pcols, lchnk)
    call outfld('AODNIR',        aodnir,        pcols, lchnk)
    call outfld('BURDENDUST',    burdendust,    pcols, lchnk)
    call outfld('BURDENSO4' ,    burdenso4,     pcols, lchnk)
    call outfld('BURDENPOM' ,    burdenpom,     pcols, lchnk)
    call outfld('BURDENSOA' ,    burdensoa,     pcols, lchnk)
    call outfld('BURDENBC'  ,    burdenbc,      pcols, lchnk)
    call outfld('BURDENSEASALT', burdenseasalt, pcols, lchnk)
    call outfld('BURDENMOM',     burdenmom,     pcols, lchnk)
    call outfld('AODABSBC',      aodabsbc,      pcols, lchnk)
    call outfld('AODDUST',       dustaod,       pcols, lchnk)
    call outfld('AODSO4',        so4aod,        pcols, lchnk)
    call outfld('AODPOM',        pomaod,        pcols, lchnk)
    call outfld('AODSOA',        soaaod,        pcols, lchnk)
    call outfld('AODBC',         bcaod,         pcols, lchnk)
    call outfld('AODSS',         seasaltaod,    pcols, lchnk)
    call outfld('AODMOM',        momaod,        pcols, lchnk)


  end subroutine modal_aero_sw

  !===============================================================================
  subroutine modal_aero_lw(dt, lchnk, ncol, state_q, temperature, pmid, pdel, pdeldry, cldn, & ! in
       qqcw, tauxar            ) ! out

    ! calculates aerosol lw radiative properties

    real(r8),          intent(in)  :: dt       ! time step [s]
    integer,           intent(in) :: lchnk            ! chunk id
    integer,           intent(in) :: ncol             ! number of active columns in the chunk
    real(r8), target,  intent(in) :: state_q(:,:,:)
    real(r8),          intent(in) :: temperature(:,:)
    real(r8),          intent(in) :: pmid(:,:)
    real(r8),          intent(in) :: pdel(:,:)
    real(r8),          intent(in) :: pdeldry(:,:)
    real(r8),  target, intent(in) :: cldn(:,:)        ! layer cloud fraction [fraction]

    type(ptr2d_t), intent(inout) :: qqcw(:)               ! Cloud borne aerosols mixing ratios [kg/kg or 1/kg]
    real(r8), intent(out) :: tauxar(pcols,pver,nlwbands) ! layer absorption optical depth

    ! Local variables
    integer :: icol, ilw, kk, ll, mm
    integer :: nspec
    integer :: istat
    integer :: itim_old           ! index

    real(r8) :: sigma_logr_aer          ! geometric standard deviation of number distribution
    real(r8) :: alnsg_amode
    real(r8) :: cheby(ncoef,pcols,pver) ! chebychef polynomials
    real(r8) :: radsurf(pcols,pver)     ! aerosol surface mode radius
    real(r8) :: logradsurf(pcols,pver)  ! log(aerosol surface mode radius)

    real(r8) :: mass(pcols,pver) ! layer mass

    real(r8),allocatable :: specvol(:,:)        ! volume concentration of aerosol specie [m3/kg]
    real(r8),allocatable :: specdens(:)         ! species density for all species [kg/m3]
    complex(r8),allocatable :: specrefindex(:,:)     ! species refractive index

    real(r8) :: dryvol(pcols)    ! volume concentration of aerosol mode [m3/kg]
    real(r8) :: wetvol(pcols)    ! volume concentration of wet mode [m3/kg]
    real(r8) :: watervol(pcols)  ! volume concentration of water in each mode [m3/kg]
    real(r8) :: refr(pcols)      ! real part of refractive index
    real(r8) :: refi(pcols)      ! imaginary part of refractive index
    complex(r8) :: crefin(pcols) ! complex refractive index
    real(r8), pointer :: specmmr(:,:)       ! species mass mixing ratio [g/g]

    integer  :: itab(pcols), jtab(pcols)
    real(r8) :: ttab(pcols), utab(pcols)
    real(r8) :: cabs(pcols,ncoef)
    real(r8) :: pabs(pcols)      ! parameterized specific absorption [m2/kg]
    real(r8) :: dopaer    ! aerosol optical depth in layer

    integer  :: nerr_dopaer = 0

    !----------------------------------------------------------------------------

    ! dry mass in each cell
    mass(:ncol,:) = pdeldry(:ncol,:)*rga

    call modal_aero_calcsize_sub(ncol, lchnk, state_q, pdel, dt, qqcw, list_idx_in=0, update_mmr_in = .false., &
         dgnumdry_m=dgnumdry_m)

    call modal_aero_wateruptake_dr(lchnk, ncol, state_q, temperature, pmid, & ! in
         cldn, dgnumdry_m, & ! in
         dgnumwet_m, qaerwat_m, & ! inout
         list_idx_in=0   ) ! optional in


    ! initialize output variables
    tauxar(:ncol,:,:) = 0._r8

    ! loop over all aerosol modes
    do mm = 1, ntot_amode

       ! get mode info
       nspec = nspec_amode(mm)
       sigma_logr_aer = sigmag_amode(mm)

       ! calc size parameter for all columns
       ! FORTRAN refactoring: ismethod2 is tempararily used to ensure BFB test. 
       ! can be removed when porting to C++
       call modal_size_parameters(ncol, sigma_logr_aer, dgnumwet_m(:,:,mm), & ! in
            radsurf, logradsurf, cheby, ismethod2=.true.) 

       allocate(specvol(ncol,nspec),stat=istat)
       if (istat /= 0) call endrun("Unable to allocate specvol: "//errmsg(__FILE__,__LINE__) )
       allocate(specdens(nspec),stat=istat)
       if (istat /= 0) call endrun("Unable to allocate specdens: "//errmsg(__FILE__,__LINE__) )
       allocate(specrefindex(nspec,nlwbands),stat=istat)
       if (istat /= 0) call endrun("Unable to allocate specrefindex: "//errmsg(__FILE__,__LINE__) )

       do ilw = 1, nlwbands

          do kk = top_lev, pver
             ! get aerosol properties and save for each species
             do ll = 1, nspec
                specmmr => state_q(:,:,lmassptr_amode(ll,mm))
                specdens(ll) = specdens_amode(lspectype_amode(ll,mm))
                specrefindex(ll,:) = specrefndxlw(:,lspectype_amode(ll,mm))
                specvol(:,ll) = specmmr(:,kk)/specdens(ll)
             enddo

             ! calculate complex refractive index
             call calc_refin_complex('lw', ncol, ilw, & ! in
                  qaerwat_m(:,kk,mm), specvol, specrefindex,   & ! in
                  dryvol, wetvol, watervol, crefin, refr, refi) ! out

             ! interpolate coefficients linear in refractive index
             ! first call calcs itab,jtab,ttab,utab
             itab(:ncol) = 0
             call binterp(absplw(mm,:,:,:,ilw), ncol,                   & !in
                  refr, refi, refrtablw(mm,:,ilw), refitablw(mm,:,ilw), & !in
                  itab, jtab, ttab, utab, cabs)                           !inout/out

             ! parameterized optical properties
             call calc_parameterized (ncol, ncoef, cabs, cheby(:,:,kk), & ! in
                  pabs) ! out
             do icol = 1, ncol
                pabs(icol)   = pabs(icol)*wetvol(icol)*rhoh2o
                pabs(icol)   = max(0._r8,pabs(icol))
                dopaer = pabs(icol)*mass(icol,kk)

                ! FORTRAN refactor: check and writeout error/warning message
                call check_error_warning('lw', icol, kk,mm, ilw, nspec,& ! in
                     dopaer, pabs(icol), dryvol, wetvol, watervol, crefin,cabs,& ! in
                     specdens, specrefindex, specvol, & ! in
                     nerr_dopaer) ! inout

                ! update absorption optical depth
                tauxar(icol,kk,ilw) = tauxar(icol,kk,ilw) + dopaer
             enddo

          enddo ! kk = top_lev, pver

       enddo  ! nlwbands

       deallocate(specvol)
       deallocate(specdens)
       deallocate(specrefindex)

    enddo ! mm = 1, ntot_amode

  end subroutine modal_aero_lw

  !===============================================================================
  ! Private routines
  !===============================================================================

  subroutine calc_parameterized (ncol, ncoef, coef, cheb_k, & ! in
       para) ! out
    ! calculate parameterized absorption, extinction or asymmetry factor
    ! further calculations are needed. see modal_aero_sw and modal_aero_lw

    implicit none
    integer,  intent(in) :: ncol,ncoef
    real(r8), intent(in) :: coef(pcols,ncoef)
    real(r8), intent(in) :: cheb_k(ncoef,pcols)
    real(r8), intent(out):: para(pcols)

    integer :: nc

    para(:ncol) = 0.5_r8*coef(:ncol,1)
    do nc = 2, ncoef
       para(:ncol) = para(:ncol) + cheb_k(nc,:ncol)*coef(:ncol,nc)
    enddo

  end subroutine calc_parameterized
  !===============================================================================
  subroutine calc_diag_spec ( ncol, specmmr_k, mass_k,            & ! in
       vol, specrefr, specrefi, hygro_aer, & ! in
       burden_s, scat_s, abs_s, hygro_s    ) ! out
    ! calculate some diagnostics for a species
    implicit none
    integer,  intent(in) :: ncol
    real(r8), intent(in),  pointer :: specmmr_k(:)
    real(r8), intent(in) :: mass_k(:)
    real(r8), intent(in) :: vol(:) ! volume concentration of aerosol species [m3/kg]
    real(r8), intent(in) :: specrefr, specrefi  ! real and image part of specrefindex
    real(r8), intent(in) :: hygro_aer        ! aerosol hygroscopicity
    real(r8), intent(out) :: burden_s(pcols) ! aerosol burden of species
    real(r8), intent(out) :: scat_s(pcols)   ! scattering of species
    real(r8), intent(out) :: abs_s(pcols)    ! absorption of species
    real(r8), intent(out) :: hygro_s(pcols)  ! hygroscopicity of species

    integer :: icol

    burden_s(:ncol) = 0._r8
    do icol = 1, ncol
       burden_s(icol) = burden_s(icol) + specmmr_k(icol)*mass_k(icol)
       scat_s(icol)   = vol(icol)*specrefr
       abs_s(icol)    = -vol(icol)*specrefi
       hygro_s(icol)  = vol(icol)*hygro_aer
    enddo

  end subroutine calc_diag_spec

  !===============================================================================
  subroutine update_aod_spec ( scath2o, absh2o,           & ! in
       sumhygro, sumscat, sumabs, & ! in
       hygro_s, palb, dopaer,     & ! in
       scat_s, abs_s, aod_s       ) ! inout
    ! update aerosol optical depth from scattering and absorption

    implicit none
    real(r8), intent(in) :: scath2o, absh2o, sumscat, sumabs, sumhygro
    real(r8), intent(in) :: hygro_s, palb, dopaer
    real(r8), intent(inout) :: scat_s, abs_s, aod_s  ! scatering, absorption and aod for a species
    ! local variables
    real(r8) :: aodc  ! aod component

    scat_s     = (scat_s + scath2o*hygro_s/sumhygro)/sumscat
    abs_s      = (abs_s + absh2o*hygro_s/sumhygro)/sumabs

    aodc           = (abs_s*(1.0_r8 - palb) + palb*scat_s)*dopaer

    aod_s      = aod_s + aodc

  end subroutine update_aod_spec

  !===============================================================================
  subroutine calc_volc_ext(ncol, trop_level, state_zm, ext_cmip6_sw, & ! in
       extinct, tropopause_m ) ! inout/out
    ! calculate contributions from volcanic aerosol extinction

    implicit none
    integer,  intent(in) :: ncol
    integer,  intent(in) :: trop_level(pcols)!tropopause level for each column
    real(r8), intent(in) :: state_zm(:,:) ! state%zm [m]
    real(r8), intent(in) :: ext_cmip6_sw(pcols,pver) ! aerosol shortwave extinction [1/m]
    real(r8), intent(inout) :: extinct(pcols,pver) ! aerosol extinction [1/m]
    real(r8), intent(out)   :: tropopause_m(pcols)

    ! local variables
    integer :: icol, kk_tropp

    do icol = 1, ncol
       kk_tropp = trop_level(icol)

       ! diagnose tropopause height
       tropopause_m(icol) = state_zm(icol,kk_tropp)!in meters

       !update tropopause layer first
       extinct(icol,kk_tropp) = 0.5_r8*( extinct(icol,kk_tropp) + ext_cmip6_sw(icol,kk_tropp) )
       !extinction is assigned read in values only for visible band above tropopause
       extinct(icol, 1:kk_tropp-1) = ext_cmip6_sw(icol, 1:kk_tropp-1)
    enddo

  end subroutine calc_volc_ext

  !===============================================================================
  subroutine calc_refin_complex (lwsw, ncol, ilwsw,           & ! in
       qaerwat_kk, specvol, specrefindex,          & ! in
       dryvol, wetvol, watervol, crefin, refr, refi) ! out
    !-------------------------------------------------------------------
    ! calculate complex refractive index 
    ! also output wetvol and watervol
    !-------------------------------------------------------------------

    implicit none
    character(len=2), intent(in) :: lwsw   ! indicator if this is lw or sw
    integer,  intent(in) :: ncol, ilwsw       
    real(r8), intent(in) :: qaerwat_kk(:)   ! aerosol water at level kk [g/g]
    real(r8), intent(in) :: specvol(:,:)    ! volume concentration of aerosol specie [m3/kg]
    complex(r8), intent(in) :: specrefindex(:,:)     ! species refractive index

    real(r8),intent(out) :: dryvol(pcols)    ! volume concentration of aerosol mode [m3/kg]
    real(r8),intent(out) :: wetvol(pcols)    ! volume concentration of wet mode [m3/kg]
    real(r8),intent(out) :: watervol(pcols)  ! volume concentration of water in each mode [m3/kg]
    real(r8),intent(out) :: refr(pcols)      ! real part of refractive index
    real(r8),intent(out) :: refi(pcols)      ! imaginary part of refractive index
    complex(r8),intent(out) :: crefin(pcols) ! complex refractive index

    integer :: icol
    real(r8), parameter :: small_value_60 = 1.e-60_r8

    if ((lwsw /= 'lw') .and. (lwsw /= 'sw')) then
       call endrun('calc_refin_complex is called with '// lwsw// ', it should be called with either lw or sw')
    endif

    crefin(:ncol) = (0._r8, 0._r8)
    dryvol(:ncol) = 0._r8

    do icol = 1, ncol
       dryvol(icol) = sum(specvol(icol,:))
       crefin(icol) = sum(specvol(icol,:)*specrefindex(:,ilwsw))

       watervol(icol) = qaerwat_kk(icol)/rhoh2o
       wetvol(icol)   = watervol(icol) + dryvol(icol)

       if (watervol(icol) < 0.0_r8 .and. lwsw=='lw') then
          if (abs(watervol(icol)) > 1.e-1_r8*wetvol(icol)) then
             write(iulog,*) 'watervol,wetvol,dryvol=',watervol(icol), &
                  wetvol(icol),dryvol(icol)
          endif
          watervol(icol) = 0._r8
          wetvol(icol)   = dryvol(icol)
       endif

       ! some different treatments for lw and sw
       if (lwsw=='lw') then
          crefin(icol) = crefin(icol) + watervol(icol)*crefwlw(ilwsw)
          if (wetvol(icol) > small_value_40) crefin(icol) = crefin(icol)/wetvol(icol)
       elseif (lwsw=='sw') then
          crefin(icol) = crefin(icol) + watervol(icol)*crefwsw(ilwsw)
          crefin(icol) = crefin(icol)/max(wetvol(icol), small_value_60)
       endif

       refr(icol) = real(crefin(icol))
       refi(icol) = aimag(crefin(icol))
    enddo

  end subroutine calc_refin_complex

  !===================================================================
  subroutine check_error_warning(lwsw, icol, kk, mm, ilwsw, nspec,   & ! in
       dopaer, pabs, dryvol, wetvol, watervol, crefin,cabs,& ! in
       specdens, specrefindex, specvol, & ! in
       nerr_dopaer,   & ! inout
       pext, specpext ) ! optional in
    !------------------------------------------------------------
    ! check and writeout error and warning message
    !------------------------------------------------------------

    implicit none
    character(len=2), intent(in) :: lwsw   ! indicator if this is lw or sw
    integer,intent(in) :: icol, kk, mm, ilwsw, nspec  ! indices
    real(r8),intent(in) :: dopaer    ! aerosol optical depth in layer
    real(r8),intent(in) :: pabs      ! parameterized specific absorption [m2/kg]
    real(r8),intent(in) :: dryvol(:)    ! volume concentration of aerosol mode [m3/kg]
    real(r8),intent(in) :: wetvol(:)    ! volume concentration of wet mode [m3/kg]
    real(r8),intent(in) :: watervol(:)  ! volume concentration of water in each mode [m3/kg]
    complex(r8),intent(in) :: crefin(:) ! complex refractive index
    real(r8),intent(in) :: cabs(:,:)
    real(r8),intent(in) :: specvol(:,:) ! volume concentration of aerosol specie [m3/kg]
    real(r8),intent(in) :: specdens(:)  ! species density [kg/m3]
    complex(r8), intent(in) :: specrefindex(:, :)     ! species refractive index
    integer, intent(inout) :: nerr_dopaer    ! total number of error times
    real(r8),intent(in),optional :: pext         ! only write out for sw
    real(r8),intent(in),optional :: specpext     ! only write out for sw

    integer :: ll
    integer, parameter :: nerrmax_dopaer=1000
    integer, parameter :: small_value_neg = -1.e-10_r8

    if ((lwsw /= 'lw') .and. (lwsw /= 'sw')) then
       call endrun('check_error_warning is called with '// lwsw// ', it should be called with either lw or sw')
    endif

    ! FORTRAN refactor: This if condition is never met in testing run ...
    if ((dopaer <= small_value_neg) .or. (dopaer >= 20._r8)) then

       if (dopaer <= small_value_neg) then
          write(iulog,*) "ERROR: Negative aerosol optical depth &
               &in this layer."
       else
          write(iulog,*) "WARNING: Aerosol optical depth is &
               &unreasonably high in this layer."
       endif

       write(iulog,*) 'dopaer(',icol,',',kk,',',mm,',)=', dopaer
       write(iulog,*) 'kk=',kk,' pabs=', pabs
       if ((lwsw == 'sw') .and. present(pext) .and. present(specpext)) then
          write(iulog,*) 'kk=', kk, ' pext=', pext, ' specext=', specpext
       endif
       write(iulog,*) 'wetvol=',wetvol(icol),' dryvol=',dryvol(icol),     &
            ' watervol=',watervol(icol)
       write(iulog,*) 'cabs=', (cabs(icol,ll),ll=1,ncoef)
       write(iulog,*) 'crefin=', crefin(icol)
       write(iulog,*) 'nspec=', nspec
       do ll = 1,nspec
          write(iulog,*) 'll=',ll,'vol(l)=',specvol(icol,ll)
          if (lwsw=='lw') then
             write(iulog,*) 'ilw=',ilwsw,' specrefindex(ilw)=',specrefindex(ll,ilwsw)
          elseif (lwsw=='sw') then
             write(iulog,*) 'isw=', ilwsw, 'specrefindex(isw)=', specrefindex(ll,ilwsw)
          endif
          write(iulog,*) 'specdens=',specdens(ll)
       enddo

       nerr_dopaer = nerr_dopaer + 1
       if (nerr_dopaer >= nerrmax_dopaer .or. dopaer < small_value_neg) then
          write(iulog,*) '*** halting after nerr_dopaer =', nerr_dopaer
          call endrun()
       endif

    endif

  end subroutine check_error_warning


  !======================================================================
  subroutine read_water_refindex(infilename)

    ! read water refractive index file and set module data

    character*(*), intent(in) :: infilename   ! modal optics filename

    ! Local variables

    integer            :: i, ierr
    type(file_desc_t)  :: ncid              ! pio file handle
    integer            :: did               ! dimension ids
    integer            :: dimlen            ! dimension lengths
    type(var_desc_t)   :: vid               ! variable ids
    real(r8) :: refrwsw(nswbands), refiwsw(nswbands) ! real, imaginary ref index for water visible
    real(r8) :: refrwlw(nlwbands), refiwlw(nlwbands) ! real, imaginary ref index for water infrared
    !----------------------------------------------------------------------------

    ! open file
    call cam_pio_openfile(ncid, infilename, PIO_NOWRITE)

    ! inquire dimensions.  Check that file values match parameter values.

    ierr = pio_inq_dimid(ncid, 'lw_band', did)
    ierr = pio_inq_dimlen(ncid, did, dimlen)
    if (dimlen .ne. nlwbands) then
       write(iulog,*) 'lw_band len=', dimlen, ' from ', infilename, ' ne nlwbands=', nlwbands
       call endrun('read_modal_optics: bad lw_band value')
    endif

    ierr = pio_inq_dimid(ncid, 'sw_band', did)
    ierr = pio_inq_dimlen(ncid, did, dimlen)
    if (dimlen .ne. nswbands) then
       write(iulog,*) 'sw_band len=', dimlen, ' from ', infilename, ' ne nswbands=', nswbands
       call endrun('read_modal_optics: bad sw_band value')
    endif

    ! read variables
    ierr = pio_inq_varid(ncid, 'refindex_real_water_sw', vid)
    ierr = pio_get_var(ncid, vid, refrwsw)

    ierr = pio_inq_varid(ncid, 'refindex_im_water_sw', vid)
    ierr = pio_get_var(ncid, vid, refiwsw)

    ierr = pio_inq_varid(ncid, 'refindex_real_water_lw', vid)
    ierr = pio_get_var(ncid, vid, refrwlw)

    ierr = pio_inq_varid(ncid, 'refindex_im_water_lw', vid)
    ierr = pio_get_var(ncid, vid, refiwlw)

    ! set complex representation of refractive indices as module data
    do i = 1, nswbands
       crefwsw(i)  = cmplx(refrwsw(i), abs(refiwsw(i)),kind=r8)
    end do
    do i = 1, nlwbands
       crefwlw(i)  = cmplx(refrwlw(i), abs(refiwlw(i)),kind=r8)
    end do

    call pio_closefile(ncid)

  end subroutine read_water_refindex

  !===============================================================================

  subroutine modal_size_parameters(ncol, sigma_logr_aer, dgnumwet, & ! in
       radsurf, logradsurf, cheb,      & ! out
       ismethod2 ) ! optional in

    use mam_support, only : min_max_bound

    integer,  intent(in)  :: ncol
    real(r8), intent(in)  :: sigma_logr_aer  ! geometric standard deviation of number distribution
    real(r8), intent(in)  :: dgnumwet(:,:)   ! aerosol wet number mode diameter [m]
    real(r8), intent(out) :: radsurf(:,:)    ! aerosol surface mode radius [m]
    real(r8), intent(out) :: logradsurf(:,:) ! log(aerosol surface mode radius)
    real(r8), intent(out) :: cheb(:,:,:)     ! chebychev polynomial parameters

    ! FORTRAN refactoring: ismethod is tempararily used to ensure BFB test
    logical,intent(in),optional :: ismethod2

    integer  :: icol, kk, nc
    real(r8) :: alnsg_amode      ! log(sigma)
    real(r8) :: explnsigma
    real(r8) :: xrad ! normalized aerosol radius
    !-------------------------------------------------------------------------------

    alnsg_amode = log(sigma_logr_aer)
    explnsigma = exp(2.0_r8*alnsg_amode*alnsg_amode)

    do kk = top_lev, pver
       do icol = 1, ncol
          ! convert from number mode diameter to surface area
          radsurf(icol,kk) = 0.5_r8*dgnumwet(icol,kk)*explnsigma

          ! --------------- FORTRAN refactoring -------------------
          ! here two calculations are used to ensure passing BFB test
          ! can be simplified (there is only round-off difference
          if (present(ismethod2) .and. ismethod2) then
             logradsurf(icol,kk) = log(0.5_r8*dgnumwet(icol,kk)) + 2.0_r8*alnsg_amode*alnsg_amode
          else
             logradsurf(icol,kk) = log(radsurf(icol,kk))
          endif
          ! --------------- FORTRAN refactoring -------------------

          ! normalize size parameter
          xrad = min_max_bound(xrmin, xrmax, logradsurf(icol,kk))
          xrad = (2._r8*xrad-xrmax-xrmin)/(xrmax-xrmin)
          ! chebyshev polynomials
          cheb(1,icol,kk) = 1._r8
          cheb(2,icol,kk) = xrad
          do nc = 3, ncoef
             cheb(nc,icol,kk) = 2._r8*xrad*cheb(nc-1,icol,kk) - cheb(nc-2,icol,kk)
          enddo
       enddo
    enddo

  end subroutine modal_size_parameters

  !===============================================================================


  subroutine binterp(table, ncol, ref_real, ref_img, ref_real_tab, ref_img_tab, &! in
       itab, jtab, ttab, utab, coef) !inout/out

    !------------------------------------------------------------------------------
    ! Bilinear interpolation along the refractive index dimensions
    ! of the table to estimate Chebyshev coefficients at an
    ! intermediate refractive index.

    ! In short wave, the first call computes itab, jtab, ttab, utab and coef.
    ! The subsequent calls use itab, jtab, ttab and utab as inputs and compute coef

    ! In long wave, we have just one call to compute itab,jtab,ttab, utab and coef
    !------------------------------------------------------------------------------
    implicit none

    !intent-ins
    integer, intent(in)  :: ncol
    real(r8), intent(in) :: table(ncoef,prefr,prefi)
    real(r8), intent(in) :: ref_real(pcols), ref_img(pcols) !real and imganinary parts of refractive indices [unitless]
    real(r8), intent(in) :: ref_real_tab(prefr), ref_img_tab(prefi) !real and imganinary table refractive indices [unitless]

    !intent-inouts/outs
    integer,  intent(inout) :: itab(pcols), jtab(pcols)
    real(r8), intent(inout) :: ttab(pcols), utab(pcols)
    real(r8), intent(out)   :: coef(pcols,ncoef) !coefficient interpolated bilinearly

    !local
    integer  :: ip1, jp1, icoef, ic
    real(r8) :: tu(pcols), tuc(pcols),tcu(pcols), tcuc(pcols)

    if(itab(1) <= 0) then

       !compute factors for the real part
       call compute_factors(prefr, ncol, ref_real, ref_real_tab, & !in
            itab, ttab) !out

       !compute factors for the imaginary part
       call compute_factors(prefi, ncol, ref_img, ref_img_tab, & !in
            jtab, utab, .true.) !out
    endif

    !Interpolate coef
    do ic = 1, ncol
       tu(ic)   = ttab(ic)*utab(ic)
       tuc(ic)  = ttab(ic)-tu(ic)
       tcuc(ic) = 1._r8-tuc(ic)-utab(ic)
       tcu(ic)  = utab(ic)-tu(ic)
       jp1      = min(jtab(ic)+1,prefi)
       ip1      = min(itab(ic)+1,prefr)
       do icoef = 1, ncoef
          coef(ic,icoef) = tcuc(ic)* table(icoef,itab(ic),jtab(ic)) + &
               tuc(ic) * table(icoef,ip1,   jtab(ic)) + &
               tu(ic)  * table(icoef,ip1,   jp1)    + &
               tcu(ic) * table(icoef,itab(ic),jp1)
       enddo
    enddo
  end subroutine binterp


  subroutine compute_factors(prefri, ncol, ref_ind, ref_table, & !in
       ix, tt, & !out
       do_print) !in-optional

    ! Compute factors for the real or imaginary parts
    implicit none

    integer, intent(in)  :: prefri, ncol
    real(r8), intent(in) :: ref_table(:) !refractive index table [unitless]
    real(r8), intent(in) :: ref_ind(:)   !refractive index       [unitless]

    !FORTRAN refactor note: "do_print" is kept to maintain code consistenty with the
    !original code. THis can be removed from the C++ ported code
    logical, intent(in), optional :: do_print ! to print log msg or not

    !intent-inouts/outs
    integer,  intent(out) :: ix(:)
    real(r8), intent(out) :: tt(:)

    !local
    integer  :: ii, ip1, ic
    real(r8) :: dx

    real(r8), parameter :: threshold = 1.e-20_r8

    ix(:ncol) = 1
    tt(:ncol) = 0._r8

    if(prefri > 1) then
       do ic = 1, ncol
          do ii = 1, prefri
             if(ref_ind(ic) < ref_table(ii)) exit
          enddo
          ix(ic) = max(ii-1,1)
          ip1    = min(ix(ic)+1,prefri)
          dx     = (ref_table(ip1)-ref_table(ix(ic)))
          if(abs(dx) > threshold) then
             tt(ic) = (ref_ind(ic)-ref_table(ix(ic)))/dx
             if (present(do_print) .and. do_print) then
                if(tt(ic) < 0._r8 .or. tt(ic) > 1._r8) then
                   write(iulog,*) 'tt,ref_ind,ix,ref_table,dx=',tt(ic),ref_ind(ic),ix(ic),ref_table(ix(ic)),dx
                endif
             endif
          endif
       enddo
    endif
  end subroutine compute_factors

end module modal_aer_opt
