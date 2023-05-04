
module ndrop
#include "../../chemistry/yaml/common_files/common_uses.ymlf90"

  !---------------------------------------------------------------------------------
  ! Purpose:
  !   CAM Interface for droplet activation by modal aerosols
  !
  ! ***N.B.*** This module is currently hardcoded to recognize only the modes that
  !            affect the climate calculation.  This is implemented by using list
  !            index 0 in all the calls to rad_constituent interfaces.
  !---------------------------------------------------------------------------------

  use shr_kind_mod,     only: r8 => shr_kind_r8
  use spmd_utils,       only: masterproc
  use ppgrid,           only: pcols, pver, pverp
  use physconst,        only: pi, rhoh2o, mwh2o, r_universal, rh2o, &
       gravit, latvap, cpair, rair
  use constituents,     only: pcnst, cnst_get_ind
  use physics_types,    only: physics_state, physics_ptend, physics_ptend_init
  use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_get_field

  use wv_saturation,    only: qsat
  use phys_control,     only: phys_getopts
  use ref_pres,         only: top_lev => trop_cloud_top_lev
#ifndef HAVE_ERF_INTRINSICS
  use shr_spfn_mod,     only: erf => shr_spfn_erf
#endif
  use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_mode_num, rad_cnst_get_aer_mmr, &
       rad_cnst_get_aer_props, rad_cnst_get_mode_props,                &
       rad_cnst_get_mam_mmr_idx, rad_cnst_get_mode_num_idx
  use cam_history,      only: addfld, horiz_only, add_default, fieldname_len, outfld
  use cam_abortutils,       only: endrun
  use cam_logfile,      only: iulog
  use modal_aero_data,   only: lmassptrcw_amode, numptrcw_amode, lmassptr_amode, numptr_amode, &
       lspectype_amode, specdens_amode, spechygro

  implicit none
  private
  save

  public ndrop_init, dropmixnuc, activate_modal
  public ptr2d_t, mam_idx, ncnst_tot  !needed by microp_aero to extract dropmixnuc input fields from pbuf

  real(r8), allocatable :: alogsig(:)     ! natl log of geometric standard dev of aerosol
  real(r8), allocatable :: exp45logsig(:)
  real(r8), allocatable :: f1(:)          ! abdul-razzak functions of width
  real(r8), allocatable :: f2(:)          ! abdul-razzak functions of width

  real(r8) :: t0            ! reference temperature
  real(r8) :: aten
  real(r8) :: surften       ! surface tension of water w/respect to air (N/m)
  real(r8) :: alog2, alog3, alogaten
  real(r8) :: third, twothird, sixth, zero
  real(r8) :: sq2, sqpi

  ! CCN diagnostic fields
  integer,  parameter :: psat=6    ! number of supersaturations to calc ccn concentration
  real(r8), parameter :: supersat(psat)= & ! supersaturation (%) to determine ccn concentration
       (/ 0.02_r8, 0.05_r8, 0.1_r8, 0.2_r8, 0.5_r8, 1.0_r8 /)
  character(len=8) :: ccn_name(psat)= &
       (/'CCN1','CCN2','CCN3','CCN4','CCN5','CCN6'/)

  ! indices in state and pbuf structures
  integer :: numliq_idx = -1
  integer :: kvh_idx    = -1

  ! description of modal aerosols
  integer               :: ntot_amode     ! number of aerosol modes
  integer,  allocatable :: nspec_amode(:) ! number of chemical species in each aerosol mode
  real(r8), allocatable :: sigmag_amode(:)! geometric standard deviation for each aerosol mode
  real(r8), allocatable :: dgnumlo_amode(:)
  real(r8), allocatable :: dgnumhi_amode(:)
  real(r8), allocatable :: voltonumblo_amode(:)
  real(r8), allocatable :: voltonumbhi_amode(:)

  logical :: history_aerosol      ! Output the MAM aerosol tendencies
  character(len=fieldname_len), allocatable :: fieldname(:)    ! names for drop nuc tendency output fields
  character(len=fieldname_len), allocatable :: fieldname_cw(:) ! names for drop nuc tendency output fields

  ! local indexing for MAM
  integer, allocatable :: mam_idx(:,:) ! table for local indexing of modal aero number and mmr
  integer :: ncnst_tot                  ! total number of mode number conc + mode species

  ! Indices for MAM species in the ptend%q array.  Needed for prognostic aerosol case.
  integer, allocatable :: mam_cnst_idx(:,:)


  ! ptr2d_t is used to create arrays of pointers to 2D fields
  type ptr2d_t
     real(r8), pointer :: fld(:,:)
  end type ptr2d_t

  ! modal aerosols
  logical :: prog_modal_aero     ! true when modal aerosols are prognostic
  logical :: lq(pcnst) = .false. ! set flags true for constituents with non-zero tendencies
  ! in the ptend object

  logical :: regen_fix

  !===============================================================================
contains
  !===============================================================================

  subroutine ndrop_init

    integer  :: ii, l, lptr, m, mm
    integer  :: nspec_max            ! max number of species in a mode
    character(len=32)   :: tmpname
    character(len=32)   :: tmpname_cw
    character(len=128)  :: long_name
    character(len=8)    :: unit
    logical :: history_amwg         ! output the variables used by the AMWG diag package
    logical :: history_verbose      ! produce verbose history output

    !-------------------------------------------------------------------------------

    ! get indices into state%q and pbuf structures
    call cnst_get_ind('NUMLIQ', numliq_idx)

    kvh_idx      = pbuf_get_index('kvh')

    zero     = 0._r8
    third    = 1._r8/3._r8
    twothird = 2._r8*third
    sixth    = 1._r8/6._r8
    sq2      = sqrt(2._r8)
    sqpi     = sqrt(pi)

    t0       = 273._r8
    surften  = 0.076_r8
    aten     = 2._r8*mwh2o*surften/(r_universal*t0*rhoh2o)
    alogaten = log(aten)
    alog2    = log(2._r8)
    alog3    = log(3._r8)

    ! get info about the modal aerosols
    ! get ntot_amode
    call rad_cnst_get_info(0, nmodes=ntot_amode)

    allocate( &
         nspec_amode(ntot_amode),  &
         sigmag_amode(ntot_amode), &
         dgnumlo_amode(ntot_amode), &
         dgnumhi_amode(ntot_amode), &
         alogsig(ntot_amode),      &
         exp45logsig(ntot_amode),  &
         f1(ntot_amode),           &
         f2(ntot_amode),           &
         voltonumblo_amode(ntot_amode), &
         voltonumbhi_amode(ntot_amode)  )

    do m = 1, ntot_amode
       ! use only if width of size distribution is prescribed

       ! get mode info
       call rad_cnst_get_info(0, m, nspec=nspec_amode(m))

       ! get mode properties
       call rad_cnst_get_mode_props(0, m, sigmag=sigmag_amode(m),  &
            dgnumhi=dgnumhi_amode(m), dgnumlo=dgnumlo_amode(m))

       alogsig(m)     = log(sigmag_amode(m))
       exp45logsig(m) = exp(4.5_r8*alogsig(m)*alogsig(m))
       f1(m)          = 0.5_r8*exp(2.5_r8*alogsig(m)*alogsig(m))
       f2(m)          = 1._r8 + 0.25_r8*alogsig(m)

       voltonumblo_amode(m) = 1._r8 / ( (pi/6._r8)*                          &
            (dgnumlo_amode(m)**3._r8)*exp(4.5_r8*alogsig(m)**2._r8) )
       voltonumbhi_amode(m) = 1._r8 / ( (pi/6._r8)*                          &
            (dgnumhi_amode(m)**3._r8)*exp(4.5_r8*alogsig(m)**2._r8) )
    end do

    ! Init the table for local indexing of mam number conc and mmr.
    ! This table uses species index 0 for the number conc.

    ! Find max number of species in all the modes, and the total
    ! number of mode number concentrations + mode species
    nspec_max = nspec_amode(1)
    ncnst_tot = nspec_amode(1) + 1
    do m = 2, ntot_amode
       nspec_max = max(nspec_max, nspec_amode(m))
       ncnst_tot = ncnst_tot + nspec_amode(m) + 1
    end do

    allocate( &
         mam_idx(ntot_amode,0:nspec_max),      &
         mam_cnst_idx(ntot_amode,0:nspec_max), &
         fieldname(ncnst_tot),                 &
         fieldname_cw(ncnst_tot)               )

    ! Local indexing compresses the mode and number/mass indicies into one index.
    ! This indexing is used by the pointer arrays used to reference state and pbuf
    ! fields.
    ii = 0
    do m = 1, ntot_amode
       do l = 0, nspec_amode(m)
          ii = ii + 1
          mam_idx(m,l) = ii
       end do
    end do

    ! Add dropmixnuc tendencies for all modal aerosol species

    call phys_getopts(history_amwg_out = history_amwg, &
         history_verbose_out = history_verbose, &
         history_aerosol_out = history_aerosol, &
         prog_modal_aero_out=prog_modal_aero, &
         regen_fix_out=regen_fix                )


    do m = 1, ntot_amode
       do l = 0, nspec_amode(m)   ! loop over number + chem constituents

          mm = mam_idx(m,l)

          unit = 'kg/m2/s'
          if (l == 0) then   ! number
             unit = '#/m2/s'
          end if

          if (l == 0) then   ! number
             call rad_cnst_get_info(0, m, num_name=tmpname, num_name_cw=tmpname_cw)
          else
             call rad_cnst_get_info(0, m, l, spec_name=tmpname, spec_name_cw=tmpname_cw)
          end if

          fieldname(mm)    = trim(tmpname) // '_mixnuc1'
          fieldname_cw(mm) = trim(tmpname_cw) // '_mixnuc1'

          if (prog_modal_aero) then

             ! To set tendencies in the ptend object need to get the constituent indices
             ! for the prognostic species
             if (l == 0) then   ! number
                call rad_cnst_get_mode_num_idx(m, lptr)
             else
                call rad_cnst_get_mam_mmr_idx(m, l, lptr)
             end if
             mam_cnst_idx(m,l) = lptr
             lq(lptr)          = .true.

             ! Add tendency fields to the history only when prognostic MAM is enabled.
             long_name = trim(tmpname) // ' dropmixnuc mixnuc column tendency'
             call addfld(fieldname(mm), horiz_only, 'A', unit, long_name)

             long_name = trim(tmpname_cw) // ' dropmixnuc mixnuc column tendency'
             call addfld(fieldname_cw(mm), horiz_only, 'A', unit, long_name)

             if (history_aerosol .and. history_verbose) then
                call add_default(fieldname(mm), 1, ' ')
                call add_default(fieldname_cw(mm), 1, ' ')
             end if



          end if

       end do
    end do

    call addfld('CCN1',(/ 'lev' /), 'A','1/cm3','CCN concentration at S=0.02%')
    call addfld('CCN2',(/ 'lev' /), 'A','1/cm3','CCN concentration at S=0.05%')
    call addfld('CCN3',(/ 'lev' /), 'A','1/cm3','CCN concentration at S=0.1%')
    call addfld('CCN4',(/ 'lev' /), 'A','1/cm3','CCN concentration at S=0.2%')
    call addfld('CCN5',(/ 'lev' /), 'A','1/cm3','CCN concentration at S=0.5%')
    call addfld('CCN6',(/ 'lev' /), 'A','1/cm3','CCN concentration at S=1.0%')


    call addfld('WTKE', (/ 'lev' /), 'A', 'm/s', 'Standard deviation of updraft velocity')
    call addfld('NDROPMIX', (/ 'lev' /), 'A', '1/kg/s', 'Droplet number mixing')
    call addfld('NDROPSRC', (/ 'lev' /), 'A', '1/kg/s', 'Droplet number source')
    call addfld('NDROPSNK', (/ 'lev' /), 'A', '1/kg/s', 'Droplet number loss by microphysics')
    call addfld('NDROPCOL', horiz_only,    'A', '1/m2', 'Column droplet number')

    ! set the add_default fields
    if (history_amwg) then
       call add_default('CCN3', 1, ' ')
    endif

    if (history_aerosol .and. prog_modal_aero) then
       do m = 1, ntot_amode
          do l = 0, nspec_amode(m)   ! loop over number + chem constituents
             mm = mam_idx(m,l)
             if (l == 0) then   ! number
                call rad_cnst_get_info(0, m, num_name=tmpname, num_name_cw=tmpname_cw)
             else
                call rad_cnst_get_info(0, m, l, spec_name=tmpname, spec_name_cw=tmpname_cw)
             end if
             fieldname(mm)    = trim(tmpname) // '_mixnuc1'
             fieldname_cw(mm) = trim(tmpname_cw) // '_mixnuc1'
          end do
       end do
    endif
  end subroutine ndrop_init

  !===============================================================================

  subroutine dropmixnuc( &
       lchnk,ncol,psetcols,dtmicro,temp,pmid,pint,pdel,rpdel,zm,   &  ! in
       state_q,ncldwtr,kvh,wsub,cldn,cldo, &  ! in
       qqcw, & ! inout
       ptend, tendnd, factnum)  !out

    ! vertical diffusion and nucleation of cloud droplets
    ! assume cloud presence controlled by cloud fraction
    ! doesn't distinguish between warm, cold clouds

    use modal_aero_data,   only: qqcw_get_field, maxd_aspectype
    use mam_support, only: min_max_bound

    ! input arguments
    integer, intent(in)  :: lchnk               ! chunk identifier
    integer, intent(in)  :: ncol                ! number of columns
    integer, intent(in)  :: psetcols            ! maximum number of columns
    real(r8), intent(in) :: dtmicro     ! time step for microphysics [s]
    real(r8), intent(in) :: temp(:,:)    ! temperature [K]
    real(r8), intent(in) :: pmid(:,:)    ! mid-level pressure [Pa]
    real(r8), intent(in) :: pint(:,:)    ! pressure at layer interfaces [Pa]
    real(r8), intent(in) :: pdel(:,:)    ! pressure thickess of layer [Pa]
    real(r8), intent(in) :: rpdel(:,:)   ! inverse of pressure thickess of layer [/Pa]
    real(r8), intent(in) :: zm(:,:)      ! geopotential height of level [m]
    real(r8), intent(in) :: state_q(:,:,:) ! aerosol mmrs [kg/kg]
    real(r8), intent(in) :: ncldwtr(:,:) ! initial droplet number mixing ratio [#/kg]
    real(r8), intent(in) :: kvh(:,:)     ! vertical diffusivity [m^2/s]
    real(r8), intent(in) :: wsub(pcols,pver)    ! subgrid vertical velocity [m/s]
    real(r8), intent(in) :: cldn(pcols,pver)    ! cloud fraction [fraction]
    real(r8), intent(in) :: cldo(pcols,pver)    ! cloud fraction on previous time step [fraction]

    ! inout arguments
    type(ptr2d_t), intent(inout) :: qqcw(:)     ! cloud-borne aerosol mass, number mixing ratios [#/kg or kg/kg]

    ! output arguments
    type(physics_ptend), intent(out)   :: ptend
    real(r8), intent(out) :: tendnd(pcols,pver) ! tendency in droplet number mixing ratio [#/kg/s]
    real(r8), intent(out) :: factnum(:,:,:)     ! activation fraction for aerosol number [fraction]

    !--------------------Local storage-------------------------------------

    integer  :: mm                  ! local array index for MAM number, species
    integer  :: nnew, nsav          ! indices for old, new time levels in substepping
    integer  :: lptr
    integer  :: ccn3d_idx   ! index of ccn3d in pbuf
    integer  :: icol        ! column index
    integer  :: imode       ! mode index
    integer  :: kk          ! level index
    integer  :: lspec      ! species index for given mode
    integer  :: lsat       !  level of supersaturation
    integer  :: spc_idx, num_idx  ! species, number indices

    real(r8), parameter :: zkmin = 0.01_r8, zkmax = 100._r8  ! min, max vertical diffusivity [m^2/s]
    real(r8), parameter :: wmixmin = 0.1_r8        ! minimum turbulence vertical velocity [m/s]

    real(r8) :: dtinv      ! inverse time step for microphysics [s^-1]
    real(r8) :: raertend(pver)  ! tendency of interstitial aerosol mass, number mixing ratios [#/kg/s or kg/kg/s]
    real(r8) :: qqcwtend(pver)  ! tendency of cloudborne aerosol mass, number mixing ratios [#/kg/s or kg/kg/s]
    real(r8) :: zs(pver) ! inverse of distance between levels [m^-1]
    real(r8) :: qcld(pver) ! cloud droplet number mixing ratio [#/kg]
    real(r8) :: csbot(pver)       ! air density at bottom (interface) of layer [kg/m^3]
    real(r8) :: csbot_cscen(pver) ! inverse normalized air density csbot(i)/cs(i,k) [dimensionless]
    real(r8) :: zn(pver)   ! g/pdel for layer [m^2/kg]
    real(r8) :: ekd(pver)       ! diffusivity for droplets [m^2/s]
    real(r8) :: ndropcol(pcols)               ! column-integrated droplet number [#/m2]
    real(r8) :: cs(pcols,pver)      ! air density [kg/m^3]
    real(r8) :: dz(pver)      ! geometric thickness of layers [m]
    real(r8) :: wtke(pcols,pver)     ! turbulent vertical velocity at base of layer k [m/s]
    real(r8) :: nsource(pcols,pver)            ! droplet number mixing ratio source tendency [#/kg/s]
    real(r8) :: ndropmix(pcols,pver)           ! droplet number mixing ratio tendency due to mixing [#/kg/s]
    real(r8) :: ccn(pcols,pver,psat)    ! number conc of aerosols activated at supersat [#/m^3]
    real(r8) :: qcldbrn(pcols,pver,maxd_aspectype,ntot_amode) ! ! cloud-borne aerosol mass mixing ratios [kg/kg]
    real(r8) :: qcldbrn_num(pcols,pver,ntot_amode) ! ! cloud-borne aerosol number mixing ratios [#/kg]


    real(r8), pointer :: ccn3d(:, :)  !  CCN at 0.2% supersat [#/m^3]
    real(r8), allocatable :: nact(:,:)  ! fractional aero. number  activation rate [/s]
    real(r8), allocatable :: mact(:,:)  ! fractional aero. mass    activation rate [/s]
    real(r8), allocatable :: raercol(:,:,:)    ! single column of aerosol mass, number mixing ratios [#/kg or kg/kg]
    real(r8), allocatable :: raercol_cw(:,:,:) ! same as raercol but for cloud-borne phase [#/kg or kg/kg]

    !     note:  activation fraction fluxes are defined as
    !     fluxn = [flux of activated aero. number into cloud [#/m^2/s]]
    !           / [aero. number conc. in updraft, just below cloudbase [#/m^3]]

    real(r8), allocatable :: coltend(:,:)       ! column tendency for diagnostic output
    real(r8), allocatable :: coltend_cw(:,:)    ! column tendency

#include "../../chemistry/yaml/cam_ndrop/f90_yaml/dropmixnuc_beg_yml.f90"

    !-------------------------------------------------------------------------------

    ! aerosol tendencies
    call physics_ptend_init(ptend, psetcols, 'ndrop_aero', lq=lq)

    !  Allocate / define local variables

    allocate( &
         nact(pver,ntot_amode),          &
         mact(pver,ntot_amode),          &
         raercol(pver,ncnst_tot,2),      &
         raercol_cw(pver,ncnst_tot,2),   &
         coltend(pcols,ncnst_tot),       &
         coltend_cw(pcols,ncnst_tot)          )

    dtinv = 1._r8/dtmicro

    !initialize variables to zero
    ndropmix(:,:) = 0._r8
    nsource(:,:) = 0._r8
    wtke(:,:)    = 0._r8
    factnum(:,:,:) = 0._r8

    !NOTE FOR C++ PORT: Get the cloud borne MMRs from AD in variable qcldbrn, do not port the code before END NOTE
    qcldbrn(:,:,:,:) = huge(qcldbrn) !store invalid values
    !END NOTE FOR C++ PORT

    ! overall_main_icol_loop
    do icol = 1, ncol

       nact(:,1:ntot_amode) = 0._r8
       mact(:,1:ntot_amode) = 0._r8
       cs(icol,:)  = pmid(icol,:)/(rair*temp(icol,:))        ! air density (kg/m3)
       dz(:)  = 1._r8/(cs(icol,:)*gravit*rpdel(icol,:)) ! layer thickness in m
       zn(:) = gravit*rpdel(icol,:)

       wtke(icol,:)     = max(wsub(icol,:),wmixmin)

       ! load number nucleated into qcld on cloud boundaries

       qcld(:)  = ncldwtr(icol,:)

       do kk = top_lev, pver-1
          zs(kk) = 1._r8/(zm(icol,kk) - zm(icol,kk+1))
          ekd(kk)   = min_max_bound(zkmin,zkmax,kvh(icol,kk+1))
          csbot(kk) = 2.0_r8*pint(icol,kk+1)/(rair*(temp(icol,kk) + temp(icol,kk+1)))
          csbot_cscen(kk) = csbot(kk)/cs(icol,kk)
       enddo
       zs(pver) = zs(pver-1)
       ekd(pver)   = 0._r8
       csbot(pver) = cs(icol,pver)
       csbot_cscen(pver) = 1.0_r8

       !  Initialize 1D (in space) versions of interstitial and cloud borne aerosol

       nsav = 1

       do imode = 1, ntot_amode
          mm = mam_idx(imode,0)
          raercol_cw(:,mm,nsav) = 0.0_r8
          raercol(:,mm,nsav)    = 0.0_r8
          raercol_cw(top_lev:pver,mm,nsav) = qqcw(mm)%fld(icol,top_lev:pver)
          num_idx = numptr_amode(imode)
          raercol(top_lev:pver,mm,nsav) = state_q(icol,top_lev:pver,num_idx)
          do lspec = 1, nspec_amode(imode)
             mm = mam_idx(imode,lspec)
             raercol_cw(top_lev:pver,mm,nsav) = qqcw(mm)%fld(icol,top_lev:pver)
             spc_idx=lmassptr_amode(lspec,imode)
             raercol(top_lev:pver,mm,nsav)    = state_q(icol,top_lev:pver,spc_idx)
          enddo
       enddo

       !  PART I:  changes of aerosol and cloud water from temporal changes in cloud fraction
       ! droplet nucleation/aerosol activation

       call update_from_newcld(cldn(icol,:),cldo(icol,:),dtinv,     &   ! in
            wtke(icol,:),temp(icol,:),cs(icol,:),state_q(icol,:,:),  &  !in
            qcld(:),raercol(:,:,nsav),raercol_cw(:,:,nsav), &  ! inout
            nsource(icol,:), factnum(icol,:,:))  ! inout

       !  PART II: changes in aerosol and cloud water from vertical profile of new cloud fraction

       call update_from_cldn_profile(cldn(icol,:),dtinv,wtke(icol,:),zs(:),dz(:),temp(icol,:), & ! in
            cs(icol,:),csbot_cscen(:),state_q(icol,:,:), & ! in
            raercol(:,:,nsav),raercol_cw(:,:,nsav), &  ! inout
            nsource(icol,:),qcld(:),factnum(icol,:,:),ekd(:),nact(:,:),mact(:,:))  !inout

       !  PART III:  perform explict integration of droplet/aerosol mixing using substepping

       nnew = 2

       call update_from_explmix(dtmicro,csbot,cldn(icol,:),zn,zs,ekd,   &  ! in
            nact,mact,qcld,raercol,raercol_cw,nsav,nnew)       ! inout
       ! droplet number

       ndropcol(icol) = 0._r8
       do kk = top_lev, pver
          ndropmix(icol,kk) = (qcld(kk) - ncldwtr(icol,kk))*dtinv - nsource(icol,kk)
          tendnd(icol,kk)   = (max(qcld(kk), 1.e-6_r8) - ncldwtr(icol,kk))*dtinv
          ndropcol(icol)   = ndropcol(icol) + ncldwtr(icol,kk)*pdel(icol,kk)
       enddo
       ndropcol(icol) = ndropcol(icol)/gravit


       raertend = 0._r8
       qqcwtend = 0._r8

       do imode = 1, ntot_amode
          do lspec = 0, nspec_amode(imode)

             mm   = mam_idx(imode,lspec)
             lptr = mam_cnst_idx(imode,lspec)

             qqcwtend(top_lev:pver) = (raercol_cw(top_lev:pver,mm,nnew) - qqcw(mm)%fld(icol,top_lev:pver))*dtinv
             qqcw(mm)%fld(icol,:) = 0.0_r8
             qqcw(mm)%fld(icol,top_lev:pver) = max(raercol_cw(top_lev:pver,mm,nnew),0.0_r8) ! update cloud-borne aerosol; HW: ensure non-negative

             if( lspec == 0 ) then
                num_idx = numptr_amode(imode)
                raertend(top_lev:pver) = (raercol(top_lev:pver,mm,nnew) - state_q(icol,top_lev:pver,num_idx))*dtinv
                qcldbrn_num(icol,top_lev:pver,imode) = qqcw(mm)%fld(icol,top_lev:pver)
             else
                spc_idx=lmassptr_amode(lspec,imode)
                raertend(top_lev:pver) = (raercol(top_lev:pver,mm,nnew) - state_q(icol,top_lev:pver,spc_idx))*dtinv
                ! Extract cloud borne MMRs from qqcw pointer
                qcldbrn(icol,top_lev:pver,lspec,imode) = qqcw(mm)%fld(icol,top_lev:pver)
             endif

             coltend(icol,mm)    = sum( pdel(icol,:)*raertend )/gravit
             coltend_cw(icol,mm) = sum( pdel(icol,:)*qqcwtend )/gravit

             ptend%q(icol,:,lptr) = 0.0_r8
             ptend%q(icol,top_lev:pver,lptr) = raertend(top_lev:pver)           ! set tendencies for interstitial aerosol

          enddo  ! lspec loop
       enddo   ! imode loop

    enddo  ! overall_main_i_loop
    ! end of main loop over i/longitude ....................................

    call outfld('NDROPCOL', ndropcol, pcols, lchnk)
    call outfld('NDROPSRC', nsource,  pcols, lchnk)
    call outfld('NDROPMIX', ndropmix, pcols, lchnk)
    call outfld('WTKE    ', wtke,     pcols, lchnk)

    !  Use interstitial and cloud-borne aerosol to compute output ccn fields.

    call ccncalc(state_q, temp, qcldbrn, qcldbrn_num, ncol, cs, ccn)

    do lsat = 1, psat
       call outfld(ccn_name(lsat), ccn(1,1,lsat), pcols, lchnk)
    enddo

    ! do column tendencies
    do imode = 1, ntot_amode
       do lspec = 0, nspec_amode(imode)
          mm = mam_idx(imode,lspec)
          call outfld(fieldname(mm),    coltend(:,mm),    pcols, lchnk)
          call outfld(fieldname_cw(mm), coltend_cw(:,mm), pcols, lchnk)
       enddo
    enddo

    deallocate( &
         nact,       &
         mact,       &
         raercol,    &
         raercol_cw, &
         coltend,    &
         coltend_cw       )

#include "../../chemistry/yaml/cam_ndrop/f90_yaml/dropmixnuc_end_yml.f90"

  end subroutine dropmixnuc

  !===============================================================================

  subroutine get_activate_frac(state_q_kload, cs_kload, cs_kk, wtke, tair, & !in
       fn, fm, fluxn, fluxm, flux_fullact)  ! out

    ! input arguments
    real(r8), intent(in) :: state_q_kload(:)        ! aerosol mmrs at level from which to load aerosol [kg/kg]
    real(r8), intent(in) :: cs_kload    ! air density at level from which to load aerosol [kg/m3]
    real(r8), intent(in) :: cs_kk       ! air density at actual vertical level [kg/m3]
    real(r8), intent(in) :: wtke        ! subgrid vertical velocity [m/s]
    real(r8), intent(in) :: tair        ! air temperature [K]

    ! output arguments
    real(r8), intent(out) :: fn(:)        ! number fraction of aerosols activated [fraction]
    real(r8), intent(out) :: fm(:)        ! mass fraction of aerosols activated [fraction]
    real(r8), intent(out) :: fluxn(:)     ! flux of activated aerosol number fraction into cloud [m/s]
    real(r8), intent(out) :: fluxm(:)     ! flux of activated aerosol mass fraction into cloud [m/s]
    real(r8), intent(out) :: flux_fullact ! flux of activated aerosol fraction assuming 100% activation [m/s]

    ! local arguments

    integer  :: phase       ! phase of aerosol
    integer  :: imode       ! mode index

    real(r8), parameter :: wmax = 10._r8    ! vertical velocity upper bound [m/s]

    real(r8) :: naermod(ntot_amode)  ! aerosol number concentration [#/m^3]
    real(r8) :: hygro(ntot_amode)    ! hygroscopicity of aerosol mode [dimensionless]
    real(r8) :: vaerosol(ntot_amode) ! aerosol volume conc [m^3/m^3]

    ! load aerosol properties, assuming external mixtures

    phase = 1 ! interstitial

    call loadaer( state_q_kload(:), &  ! in
         cs_kload, phase, & ! in
         naermod, vaerosol, hygro )   ! out

    ! Below is to avoid warning about not assigning value to intent(out)
    !  the assignment should have no affect because flux_fullact is intent(out)
    !  in activate_modal (and in that subroutine is initialized to zero anyway).

    flux_fullact=0._r8

    call activate_modal( wtke, wmax, tair, cs_kk, &   ! in
         naermod, ntot_amode, vaerosol, hygro, &  ! in
         fn, fm, fluxn, fluxm, flux_fullact )  ! out

  end subroutine get_activate_frac

  !===============================================================================

  subroutine update_from_newcld(cldn_col_in,cldo_col_in,dtinv, & ! in
       wtke_col_in,temp_col_in,cs_col_in,state_q_col_in, & ! in
       qcld,raercol_nsav,raercol_cw_nsav, &      ! inout
       nsource_col_out, factnum_col_out)              ! inout

    ! input arguments
    real(r8), intent(in) :: cldn_col_in(:)   ! cloud fraction [fraction]
    real(r8), intent(in) :: cldo_col_in(:)   ! cloud fraction on previous time step [fraction]
    real(r8), intent(in) :: dtinv     ! inverse time step for microphysics [s^{-1}]
    real(r8), intent(in) :: wtke_col_in(:)   ! subgrid vertical velocity [m/s]
    real(r8), intent(in) :: temp_col_in(:)   ! temperature [K]
    real(r8), intent(in) :: cs_col_in(:)     ! air density at actual level kk [kg/m^3]
    real(r8), intent(in) :: state_q_col_in(:,:) ! aerosol mmrs [kg/kg]

    real(r8), intent(inout) :: qcld(:)  ! cloud droplet number mixing ratio [#/kg]
    real(r8), intent(inout) :: nsource_col_out(:)   ! droplet number mixing ratio source tendency [#/kg/s]
    real(r8), intent(inout) :: raercol_nsav(:,:)   ! single column of saved aerosol mass, number mixing ratios [#/kg or kg/kg]
    real(r8), intent(inout) :: raercol_cw_nsav(:,:)  ! same as raercol_nsav but for cloud-borne phase [#/kg or kg/kg]
    real(r8), intent(inout) :: factnum_col_out(:,:)  ! activation fraction for aerosol number [fraction]


    !  local variables
    integer  :: imode           ! mode counter variable
    integer  :: lspec           ! species counter variable
    integer  :: mm              ! local array index for MAM number, species
    integer  :: kk              ! vertical level index
    integer  :: num_idx         ! number index
    integer  :: spc_idx         ! species index

    real(r8), parameter :: grow_cld_thresh = 0.01_r8   !  threshold cloud fraction growth [fraction]

    real(r8) :: delt_cld        ! new - old cloud fraction [fraction]
    real(r8) :: frac_delt_cld   ! fractional change in cloud fraction [fraction]
    real(r8) :: fm_delt_cld     ! fm change from fractional change in cloud fraction [fraction]
    real(r8) :: dact             ! cloud-borne aerosol tendency due to cloud frac tendency [#/kg or kg/kg]
    real(r8) :: fn(ntot_amode)              ! activation fraction for aerosol number [fraction]
    real(r8) :: fm(ntot_amode)              ! activation fraction for aerosol mass [fraction]
    real(r8) :: flux_fullact ! flux of activated aerosol fraction assuming 100% activation [m/s]
    real(r8) :: fluxn(ntot_amode)     ! flux of activated aerosol number fraction into cloud [m/s]
    real(r8) :: fluxm(ntot_amode)     ! flux of activated aerosol mass fraction into cloud [m/s]


    ! k-loop for growing/shrinking cloud calcs .............................
    do kk = top_lev, pver

       delt_cld = (cldn_col_in(kk) - cldo_col_in(kk))

       ! shrinking cloud ......................................................
       !    treat the reduction of cloud fraction from when cldn(i,k) < cldo(i,k)
       !    and also dissipate the portion of the cloud that will be regenerated

       if (cldn_col_in(kk) < cldo_col_in(kk)) then
          !  droplet loss in decaying cloud
          !++ sungsup
          nsource_col_out(kk) = nsource_col_out(kk) &
               + qcld(kk)*(cldn_col_in(kk) - cldo_col_in(kk))/cldo_col_in(kk)*dtinv
          qcld(kk) =  qcld(kk)*(1._r8 + (cldn_col_in(kk)-cldo_col_in(kk))/cldo_col_in(kk))
          !-- sungsup

          ! convert activated aerosol to interstitial in decaying cloud

          frac_delt_cld = (cldn_col_in(kk) - cldo_col_in(kk)) / cldo_col_in(kk)

          do imode = 1, ntot_amode
             mm = mam_idx(imode,0)
             dact   = raercol_cw_nsav(kk,mm)*frac_delt_cld
             raercol_cw_nsav(kk,mm) = raercol_cw_nsav(kk,mm) + dact   ! cloud-borne aerosol
             raercol_nsav(kk,mm)    = raercol_nsav(kk,mm) - dact
             do lspec = 1, nspec_amode(imode)
                mm = mam_idx(imode,lspec)
                dact    = raercol_cw_nsav(kk,mm)*frac_delt_cld
                raercol_cw_nsav(kk,mm) = raercol_cw_nsav(kk,mm) + dact  ! cloud-borne aerosol
                raercol_nsav(kk,mm)    = raercol_nsav(kk,mm) - dact
             enddo
          enddo
       endif   ! cldn(icol,kk) < cldo(icol,kk)

       ! growing cloud ......................................................
       !    treat the increase of cloud fraction from when cldn(i,k) > cldo(i,k)
       !    and also regenerate part of the cloud

       if (cldn_col_in(kk)-cldo_col_in(kk) > grow_cld_thresh) then

          call get_activate_frac(state_q_col_in(kk,:), cs_col_in(kk), cs_col_in(kk), & ! in
               wtke_col_in(kk), temp_col_in(kk), & ! in
               fn, fm, fluxn, fluxm, flux_fullact) ! out

          !  store for output activation fraction of aerosol
          factnum_col_out(kk,:) = fn

          do imode = 1, ntot_amode
             mm = mam_idx(imode,0)
             num_idx = numptr_amode(imode)
             dact = delt_cld*fn(imode)*state_q_col_in(kk,num_idx) ! interstitial only
             qcld(kk) = qcld(kk) + dact
             nsource_col_out(kk) = nsource_col_out(kk) + dact*dtinv
             raercol_cw_nsav(kk,mm) = raercol_cw_nsav(kk,mm) + dact  ! cloud-borne aerosol
             raercol_nsav(kk,mm)    = raercol_nsav(kk,mm) - dact
             fm_delt_cld = delt_cld * fm(imode)
             do lspec = 1, nspec_amode(imode)
                mm = mam_idx(imode,lspec)
                spc_idx=lmassptr_amode(lspec,imode)
                dact    = fm_delt_cld*state_q_col_in(kk,spc_idx) ! interstitial only
                raercol_cw_nsav(kk,mm) = raercol_cw_nsav(kk,mm) + dact  ! cloud-borne aerosol
                raercol_nsav(kk,mm)    = raercol_nsav(kk,mm) - dact
             enddo
          enddo
       endif    !  cldn(icol,kk)-cldo(icol,kk) > 0.01_r8
    enddo  ! end of k-loop for growing/shrinking cloud calcs ......................



  end subroutine update_from_newcld

  !===============================================================================

  subroutine update_from_cldn_profile(cldn_col_in,dtinv,wtke_col_in,zs,dz,  &  ! in
       temp_col_in, cs_col_in,csbot_cscen,state_q_col_in,  &  ! in
       raercol_nsav,raercol_cw_nsav,nsource_col, &  ! inout
       qcld,factnum_col,ekd,nact,mact)  ! inout

    ! input arguments
    real(r8), intent(in) :: cldn_col_in(:)   ! cloud fraction [fraction]
    real(r8), intent(in) :: dtinv      ! inverse time step for microphysics [s^{-1}]
    real(r8), intent(in) :: wtke_col_in(:)   ! subgrid vertical velocity [m/s]
    real(r8), intent(in) :: zs(:)            ! inverse of distance between levels [m^-1]
    real(r8), intent(in) :: dz(:)         ! geometric thickness of layers [m]
    real(r8), intent(in) :: temp_col_in(:)   ! temperature [K]
    real(r8), intent(in) :: cs_col_in(:)     ! air density [kg/m^3]
    real(r8), intent(in) :: csbot_cscen(:)   ! inverse normalized air density csbot(i)/cs(i,k) [dimensionless]
    real(r8), intent(in) :: state_q_col_in(:,:)    ! aerosol mmrs [kg/kg]

    real(r8), intent(inout) :: raercol_nsav(:,:)    ! single column of saved aerosol mass, number mixing ratios [#/kg or kg/kg]
    real(r8), intent(inout) :: raercol_cw_nsav(:,:) ! same as raercol but for cloud-borne phase [#/kg or kg/kg]
    real(r8), intent(inout) :: nsource_col(:)   ! droplet number mixing ratio source tendency [#/kg/s]
    real(r8), intent(inout) :: qcld(:)  ! cloud droplet number mixing ratio [#/kg]
    real(r8), intent(inout) :: factnum_col(:,:)  ! activation fraction for aerosol number [fraction]
    real(r8), intent(inout) :: ekd(:)     ! diffusivity for droplets [m^2/s]
    real(r8), intent(inout) :: nact(:,:)  ! fractional aero. number  activation rate [/s]
    real(r8), intent(inout) :: mact(:,:)  ! fractional aero. mass    activation rate [/s]

    ! local arguments

    integer :: kk         ! vertical level index
    integer :: kp1        ! bounded vertical level index + 1
    integer :: imode      ! mode counter variable
    integer :: lspec      ! species counter variable
    integer :: mm         ! local array index for MAM number, species

    real(r8), parameter :: cld_thresh = 0.01_r8   !  threshold cloud fraction [fraction]

    real(r8) :: delz_cld        ! vertical change in cloud raction [fraction]
    real(r8) :: crdz            ! conversion factor from flux to rate [m^{-1}]
    real(r8) :: fn(ntot_amode)              ! activation fraction for aerosol number [fraction]
    real(r8) :: fm(ntot_amode)              ! activation fraction for aerosol mass [fraction]
    real(r8) :: flux_fullact(pver) ! flux of activated aerosol fraction assuming 100% activation [m/s]
    real(r8) :: fluxn(ntot_amode)     ! flux of activated aerosol number fraction into cloud [m/s]
    real(r8) :: fluxm(ntot_amode)     ! flux of activated aerosol mass fraction into cloud [m/s]
    real(r8) :: fluxntot         ! flux of activated aerosol number into cloud [#/m^2/s]


    ! ......................................................................
    ! start of k-loop for calc of old cloud activation tendencies ..........
    !
    ! rce-comment
    !    changed this part of code to use current cloud fraction (cldn) exclusively
    !    consider case of cldo(:)=0, cldn(k)=1, cldn(k+1)=0
    !    previous code (which used cldo below here) would have no cloud-base activation
    !       into layer k.  however, activated particles in k mix out to k+1,
    !       so they are incorrectly depleted with no replacement

    do kk = top_lev, pver - 1

       kp1 = min0(kk+1, pver)

       if (cldn_col_in(kk) > cld_thresh) then

          if (cldn_col_in(kk) - cldn_col_in(kp1) > cld_thresh ) then

             ! cloud base

             ! rce-comments
             !   first, should probably have 1/zs(k) here rather than dz(i,k) because
             !      the turbulent flux is proportional to ekd(k)*zs(k),
             !      while the dz(i,k) is used to get flux divergences
             !      and mixing ratio tendency/change
             !   second and more importantly, using a single updraft velocity here
             !      means having monodisperse turbulent updraft and downdrafts.
             !      The sq2pi factor assumes a normal draft spectrum.
             !      The fluxn/fluxm from activate must be consistent with the
             !      fluxes calculated in explmix.
             ekd(kk) = wtke_col_in(kk)/zs(kk)
             ! rce-comment - use kp1 here as old-cloud activation involves
             !   aerosol from layer below

             call get_activate_frac(state_q_col_in(kp1,:),cs_col_in(kp1),  &   ! in
                  cs_col_in(kk), wtke_col_in(kk), temp_col_in(kk),  &   ! in
                  fn, fm, fluxn, fluxm, flux_fullact(kk) )  ! out

             !  store for output activation fraction of aerosol
             factnum_col(kk,:) = fn
             delz_cld = cldn_col_in(kk) - cldn_col_in(kp1)
             fluxntot = 0

             ! rce-comment 1
             !    flux of activated mass into layer k (in kg/m2/s)
             !       = "actmassflux" = dumc*fluxm*raercol(kp1,lmass)*csbot(k)
             !    source of activated mass (in kg/kg/s) = flux divergence
             !       = actmassflux/(cs(i,k)*dz(i,k))
             !    so need factor of csbot_cscen = csbot(k)/cs(i,k)
             !                   dum=1./(dz(i,k))
             crdz=csbot_cscen(kk)/(dz(kk))
             ! rce-comment 2
             !    code for k=pver was changed to use the following conceptual model
             !    in k=pver, there can be no cloud-base activation unless one considers
             !       a scenario such as the layer being partially cloudy,
             !       with clear air at bottom and cloudy air at top
             !    assume this scenario, and that the clear/cloudy portions mix with
             !       a timescale taumix_internal = dz(i,pver)/wtke_cen(i,pver)
             !    in the absence of other sources/sinks, qact (the activated particle
             !       mixratio) attains a steady state value given by
             !          qact_ss = fcloud*fact*qtot
             !       where fcloud is cloud fraction, fact is activation fraction,
             !       qtot=qact+qint, qint is interstitial particle mixratio
             !    the activation rate (from mixing within the layer) can now be
             !       written as
             !          d(qact)/dt = (qact_ss - qact)/taumix_internal
             !                     = qtot*(fcloud*fact*wtke/dz) - qact*(wtke/dz)
             !    note that (fcloud*fact*wtke/dz) is equal to the nact/mact
             !    also, d(qact)/dt can be negative.  in the code below
             !       it is forced to be >= 0
             !
             ! steve --
             !    you will likely want to change this.  i did not really understand
             !       what was previously being done in k=pver
             !    in the cam3_5_3 code, wtke(i,pver) appears to be equal to the
             !       droplet deposition velocity which is quite small
             !    in the cam3_5_37 version, wtke is done differently and is much
             !       larger in k=pver, so the activation is stronger there
             !

             do imode = 1, ntot_amode
                mm = mam_idx(imode,0)
                fluxn(imode) = fluxn(imode)*delz_cld
                fluxm(imode) = fluxm(imode)*delz_cld
                nact(kk,imode) = nact(kk,imode) + fluxn(imode)*crdz
                mact(kk,imode) = mact(kk,imode) + fluxm(imode)*crdz
                ! note that kp1 is used here
                fluxntot = fluxntot &
                     + fluxn(imode)*raercol_nsav(kp1,mm)*cs_col_in(kk)
             enddo
             nsource_col(kk) = nsource_col(kk) + fluxntot/(cs_col_in(kk)*dz(kk))

          endif  ! (cldn(icol,kk) - cldn(icol,kp1) > 0.01)

       else
          !  if cldn < 0.01_r8 at any level except kk=pver, deplete qcld, turn all raercol_cw to raercol, put appropriate tendency
          !  in nsource
          !  Note that if cldn(kk) >= 0.01_r8 but cldn(kk) - cldn(kp1)  <= 0.01, nothing is done.

          ! no cloud

          nsource_col(kk) = nsource_col(kk) - qcld(kk)*dtinv
          qcld(kk)      = 0

          ! convert activated aerosol to interstitial in decaying cloud

          do imode = 1, ntot_amode
             mm = mam_idx(imode,0)
             raercol_nsav(kk,mm)    = raercol_nsav(kk,mm) + raercol_cw_nsav(kk,mm)  ! cloud-borne aerosol
             raercol_cw_nsav(kk,mm) = 0._r8

             do lspec = 1, nspec_amode(imode)
                mm = mam_idx(imode,lspec)
                raercol_nsav(kk,mm)    = raercol_nsav(kk,mm) + raercol_cw_nsav(kk,mm) ! cloud-borne aerosol
                raercol_cw_nsav(kk,mm) = 0._r8
             enddo
          enddo

       endif  ! (cldn(icol,kk) > 0.01_r8) if-else structure

    enddo


  end subroutine update_from_cldn_profile

  !===============================================================================

  subroutine update_from_explmix(dtmicro,csbot,cldn_col,zn,zs,ekd,   &  ! in
       nact,mact,qcld,raercol,raercol_cw,nsav,nnew)  ! inout

    ! input arguments
    real(r8), intent(in) :: dtmicro     ! time step for microphysics [s]
    real(r8), intent(in) :: csbot(pver)       ! air density at bottom (interface) of layer [kg/m^3]
    real(r8), intent(in) :: cldn_col(:)   ! cloud fraction [fraction]
    real(r8), intent(in) :: zn(pver)   ! g/pdel for layer [m^2/kg]
    real(r8), intent(in) :: zs(:)            ! inverse of distance between levels [m^-1]
    real(r8), intent(in) :: ekd(:)     ! diffusivity for droplets [m^2/s]

    ! in/out arguments
    real(r8), intent(inout) :: nact(:,:)  ! fractional aero. number  activation rate [/s]
    real(r8), intent(inout) :: mact(:,:)  ! fractional aero. mass    activation rate [/s]
    real(r8), intent(inout) :: qcld(:)  ! cloud droplet number mixing ratio [#/kg]
    real(r8), intent(inout) :: raercol(:,:,:)    ! single column of saved aerosol mass, number mixing ratios [#/kg or kg/kg]
    real(r8), intent(inout) :: raercol_cw(:,:,:) ! same as raercol but for cloud-borne phase [#/kg or kg/kg]
    integer, intent(inout) :: nnew, nsav   ! indices for old, new time levels in substepping

    ! local arguments
    integer :: kk           ! vertical level index
    integer :: kp1          ! kk+1
    integer :: km1          ! kk-1
    integer  :: imode       ! mode counter variable
    integer  :: mm          ! local array index for MAM number, species
    integer  :: lspec       ! species counter variable
    integer  :: nsubmix, nsubmix_bnd  ! number of substeps and bound
    integer  :: ntemp   ! temporary index for substepping
    integer  :: isub       ! substep index

    real(r8), parameter :: overlap_cld_thresh = 1.e-10_r8   !  threshold cloud fraction to compute overlap [fraction]

    real(r8) :: ekk(0:pver)     ! density*diffusivity for droplets [kg/m/s]
    real(r8) :: dtmin      ! time step to determine subloop time step [s]
    real(r8) :: qncld(pver)     ! updated cloud droplet number mixing ratio [#/kg]
    real(r8) :: ekkp(pver)      ! zn*zs*density*diffusivity [/s]
    real(r8) :: ekkm(pver)      ! zn*zs*density*diffusivity   [/s]
    real(r8) :: overlapp(pver)  ! cloud overlap involving level kk+1 [fraction]
    real(r8) :: overlapm(pver)  ! cloud overlap involving level kk-1 [fraction]
    real(r8) :: source(pver)  !  source rate for activated number or species mass [/s]
    real(r8) :: tinv       ! inverse timescale of droplet diffusivity [/s]
    real(r8) :: dtt        ! timescale of droplet diffusivity [s]
    real(r8) :: dtmix    ! timescale for subloop [s]
    real(r8) :: tmpa             !  temporary aerosol tendency variable [/s]
    real(r8) :: srcn(pver)       ! droplet source rate [/s]

    ! load new droplets in layers above, below clouds

    dtmin     = dtmicro
    ekk(top_lev-1)    = 0.0_r8
    ekk(pver) = 0.0_r8
    do kk = top_lev, pver-1
       ! rce-comment -- ekd(k) is eddy-diffusivity at k/k+1 interface
       !   want ekk(k) = ekd(k) * (density at k/k+1 interface)
       !   so use pint(i,k+1) as pint is 1:pverp
       !           ekk(k)=ekd(k)*2.*pint(i,k)/(rair*(temp(i,k)+temp(i,k+1)))
       !           ekk(k)=ekd(k)*2.*pint(i,k+1)/(rair*(temp(i,k)+temp(i,k+1)))
       ekk(kk) = ekd(kk)*csbot(kk)
    enddo


    do kk = top_lev, pver
       kp1 = min(kk+1, pver)
       km1 = max(kk-1, top_lev)

       ! maximum overlap assumption
       if (cldn_col(kp1) > overlap_cld_thresh) then
          overlapp(kk) = min(cldn_col(kk)/cldn_col(kp1), 1._r8)
       else
          overlapp(kk) = 1._r8
       endif
       if (cldn_col(km1) > overlap_cld_thresh) then
          overlapm(kk) = min(cldn_col(kk)/cldn_col(km1), 1._r8)
       else
          overlapm(kk) = 1._r8
       endif

       ekkp(kk) = zn(kk)*ekk(kk)*zs(kk)
       ekkm(kk) = zn(kk)*ekk(kk-1)*zs(km1)
       tinv   = ekkp(kk) + ekkm(kk)

       ! rce-comment -- tinv is the sum of all first-order-loss-rates
       !    for the layer.  for most layers, the activation loss rate
       !    (for interstitial particles) is accounted for by the loss by
       !    turb-transfer to the layer above.
       !    k=pver is special, and the loss rate for activation within
       !    the layer must be added to tinv.  if not, the time step
       !    can be too big, and explmix can produce negative values.
       !    the negative values are reset to zero, resulting in an
       !    artificial source.

       if (tinv > 1.e-6_r8) then
          dtt   = 1._r8/tinv
          dtmin = min(dtmin, dtt)
       endif
    enddo

    dtmix   = 0.9_r8*dtmin
    nsubmix = dtmicro/dtmix + 1
    if (nsubmix > 100) then
       nsubmix_bnd = 100
    else
       nsubmix_bnd = nsubmix
    endif
    dtmix = dtmicro/nsubmix

    ! rce-comment
    !    the activation source(k) = mact(k,m)*raercol(kp1,lmass)
    !       should not exceed the rate of transfer of unactivated particles
    !       from kp1 to k which = ekkp(k)*raercol(kp1,lmass)
    !    however it might if things are not "just right" in subr activate
    !    the following is a safety measure to avoid negatives in explmix
    do kk = top_lev, pver-1
       do imode = 1, ntot_amode
          nact(kk,imode) = min( nact(kk,imode), ekkp(kk) )
          mact(kk,imode) = min( mact(kk,imode), ekkp(kk) )
       enddo
    enddo

    ! old_cloud_nsubmix_loop

    !  Note:  each pass in submix loop stores updated aerosol values at index nnew,
    !  current values at index nsav.  At the start of each pass, nnew values are
    !  copied to nsav.  However, this is accomplished by switching the values
    !  of nsav and nnew rather than a physical copying.  At end of loop nnew stores
    !  index of most recent updated values (either 1 or 2).


    do isub = 1, nsubmix
       qncld(:) = qcld(:)
       ! after first pass, switch nsav, nnew so that nsav is the recently updated aerosol
       if( isub > 1 ) then
          ntemp   = nsav
          nsav    = nnew
          nnew    = ntemp
       endif
       srcn(:) = 0.0_r8

       do imode = 1, ntot_amode
          mm = mam_idx(imode,0)

          ! update droplet source

          ! rce-comment- activation source in layer k involves particles from k+1
          !	       srcn(:)=srcn(:)+nact(:,m)*(raercol(:,mm,nsav))
          srcn(top_lev:pver-1) = srcn(top_lev:pver-1) + nact(top_lev:pver-1,imode)*(raercol(top_lev+1:pver,mm,nsav))

          ! rce-comment- new formulation for k=pver
          !              srcn(  pver  )=srcn(  pver  )+nact(  pver  ,m)*(raercol(  pver,mm,nsav))
          tmpa = raercol(pver,mm,nsav)*nact(pver,imode) &
               + raercol_cw(pver,mm,nsav)*nact(pver,imode)
          srcn(pver) = srcn(pver) + max(0.0_r8,tmpa)
       enddo
       call explmix(  qcld, &   ! out
            srcn, ekkp, ekkm, overlapp,  &   ! in
            overlapm, qncld,  &  ! in
            dtmix, .false. )     ! in
       ! update aerosol number

       ! rce-comment
       !    the interstitial particle mixratio is different in clear/cloudy portions
       !    of a layer, and generally higher in the clear portion.  (we have/had
       !    a method for diagnosing the the clear/cloudy mixratios.)  the activation
       !    source terms involve clear air (from below) moving into cloudy air (above).
       !    in theory, the clear-portion mixratio should be used when calculating
       !    source terms
       do imode = 1, ntot_amode
          mm = mam_idx(imode,0)
          ! rce-comment -   activation source in layer k involves particles from k+1
          !	              source(:)= nact(:,m)*(raercol(:,mm,nsav))
          source(top_lev:pver-1) = nact(top_lev:pver-1,imode)*(raercol(top_lev+1:pver,mm,nsav))
          ! rce-comment - new formulation for k=pver
          !               source(  pver  )= nact(  pver,  m)*(raercol(  pver,mm,nsav))
          tmpa = raercol(pver,mm,nsav)*nact(pver,imode) &
               + raercol_cw(pver,mm,nsav)*nact(pver,imode)
          source(pver) = max(0.0_r8, tmpa)

          call explmix( raercol_cw(:,mm,nnew), &  ! out
               source, ekkp, ekkm, overlapp, &   ! in
               overlapm, raercol_cw(:,mm,nsav),   &  ! in
               dtmix, .false. )  ! in
          call explmix( raercol(:,mm,nnew), &   ! out
               source, ekkp, ekkm, overlapp,  &    ! in
               overlapm, raercol(:,mm,nsav),  &  ! in
               dtmix, .true., &  ! in
               raercol_cw(:,mm,nsav))  ! optional in

          ! update aerosol species mass

          do lspec = 1, nspec_amode(imode)
             mm = mam_idx(imode,lspec)
             ! rce-comment -   activation source in layer k involves particles from k+1
             !	          source(:)= mact(:,m)*(raercol(:,mm,nsav))
             source(top_lev:pver-1) = mact(top_lev:pver-1,imode)*(raercol(top_lev+1:pver,mm,nsav))
             ! rce-comment- new formulation for k=pver
             !                 source(  pver  )= mact(  pver  ,m)*(raercol(  pver,mm,nsav))
             tmpa = raercol(pver,mm,nsav)*mact(pver,imode) &
                  + raercol_cw(pver,mm,nsav)*mact(pver,imode)
             source(pver) = max(0.0_r8, tmpa)

             call explmix( raercol_cw(:,mm,nnew), &  ! out
                  source, ekkp, ekkm, overlapp, &  ! in
                  overlapm, raercol_cw(:,mm,nsav),    &  ! in
                  dtmix, .false. ) ! in
             call explmix( raercol(:,mm,nnew), &  ! out
                  source, ekkp, ekkm, overlapp,  &  ! in
                  overlapm, raercol(:,mm,nsav),  &   ! in
                  dtmix, .true., & ! in
                  raercol_cw(:,mm,nsav))  ! optional in

          enddo  ! lspec loop
       enddo  !  imode loop

    enddo ! old_cloud_nsubmix_loop

    ! evaporate particles again if no cloud

    do kk = top_lev, pver
       if (cldn_col(kk) == 0._r8) then
          ! no cloud
          qcld(kk)=0._r8

          ! convert activated aerosol to interstitial in decaying cloud
          do imode = 1, ntot_amode
             mm = mam_idx(imode,0)
             raercol(kk,mm,nnew)    = raercol(kk,mm,nnew) + raercol_cw(kk,mm,nnew)
             raercol_cw(kk,mm,nnew) = 0._r8

             do lspec = 1, nspec_amode(imode)
                mm = mam_idx(imode,lspec)
                raercol(kk,mm,nnew)    = raercol(kk,mm,nnew) + raercol_cw(kk,mm,nnew)
                raercol_cw(kk,mm,nnew) = 0._r8
             enddo
          enddo
       endif
    enddo

  end subroutine update_from_explmix

  !===============================================================================

  subroutine explmix( qnew, &  ! out
       src, ekkp, ekkm, overlapp, overlapm, qold, dtmix, is_unact, &  ! in
       qactold )  ! optional in

    !  explicit integration of droplet/aerosol mixing
    !     with source due to activation/nucleation

    ! output arguments
    real(r8), intent(out) :: qnew(pver) ! number / mass mixing ratio to be updated [# or kg / kg]

    ! input arguments
    real(r8), intent(in) :: qold(pver) ! number / mass mixing ratio from previous time step [# or kg / kg]
    real(r8), intent(in) :: src(pver) ! source due to activation/nucleation [# or kg / (kg-s)]
    real(r8), intent(in) :: ekkp(pver) ! zn*zs*density*diffusivity (kg/m3 m2/s) at interface [/s]
    ! below layer k  (k,k+1 interface)
    real(r8), intent(in) :: ekkm(pver) ! zn*zs*density*diffusivity (kg/m3 m2/s) at interface [/s]
    ! above layer k  (k,k+1 interface)
    real(r8), intent(in) :: overlapp(pver) ! cloud overlap below [fraction]
    real(r8), intent(in) :: overlapm(pver) ! cloud overlap above [fraction]
    real(r8), intent(in) :: dtmix ! time step [s]
    logical, intent(in) :: is_unact ! true if this is an unactivated species
    real(r8), intent(in),optional :: qactold(pver) ! [# or kg / kg]
    ! number / mass mixing ratio of ACTIVATED species from previous step
    ! *** this should only be present
    !     if the current species is unactivated number/sfc/mass

    ! local arguments

    integer kk   !  vertical level index
    integer kp1  !  bounded vertical level index plus 1
    integer km1  !  bounded vertical level index minus 1

#include "../../chemistry/yaml/cam_ndrop/f90_yaml/explmix_beg_yml.f90"

    !     the qactold*(1-overlap) terms are resuspension of activated material
    do kk=top_lev,pver
       kp1=min(kk+1,pver)
       km1=max(kk-1,top_lev)

       if ( is_unact ) then
          qnew(kk) = qold(kk) + dtmix*( - src(kk) + ekkp(kk)*(qold(kp1) - qold(kk) +       &
               qactold(kp1)*(1.0_r8-overlapp(kk)))               &
               + ekkm(kk)*(qold(km1) - qold(kk) +     &
               qactold(km1)*(1.0_r8-overlapm(kk))) )
       else
          qnew(kk) = qold(kk) + dtmix*(src(kk) + ekkp(kk)*(overlapp(kk)*qold(kp1)-qold(kk)) +      &
               ekkm(kk)*(overlapm(kk)*qold(km1)-qold(kk)) )
       endif

       !        force to non-negative
       qnew(kk)=max(qnew(kk),0._r8)
    enddo

#include "../../chemistry/yaml/cam_ndrop/f90_yaml/explmix_end_yml.f90"

  end subroutine explmix

  !===============================================================================

  subroutine activate_modal(w_in, wmaxf, tair, rhoair,  &  ! in
       na, nmode, volume, hygro, &  ! in
       fn, fm, fluxn, fluxm, flux_fullact, smax_prescribed )  ! out

    !---------------------------------------------------------------------------------
    !Calculates number, surface, and mass fraction of aerosols activated as CCN
    !calculates flux of cloud droplets, surface area, and aerosol mass into cloud
    !assumes an internal mixture within each of up to nmode multiple aerosol modes
    !a gaussiam spectrum of updrafts can be treated.
    !
    !Units: SI (MKS)
    !
    !Reference: Abdul-Razzak and Ghan, A parameterization of aerosol activation.
    !      2. Multiple aerosol types. J. Geophys. Res., 105, 6837-6844.
    !---------------------------------------------------------------------------------

    !input
    real(r8), intent(in) :: w_in      ! vertical velocity [m/s]
    real(r8), intent(in) :: wmaxf     ! maximum updraft velocity for integration [m/s]
    real(r8), intent(in) :: tair      ! air temperature [K]
    real(r8), intent(in) :: rhoair    ! air density [kg/m3]
    real(r8), intent(in) :: na(:)     ! aerosol number concentration [#/m3]
    integer,  intent(in) :: nmode     ! number of aerosol modes
    real(r8), intent(in) :: volume(:) ! aerosol volume concentration [m3/m3]
    real(r8), intent(in) :: hygro(:)  ! hygroscopicity of aerosol mode [dimensionless]

    !output
    real(r8), intent(out) :: fn(:)        ! number fraction of aerosols activated [fraction]
    real(r8), intent(out) :: fm(:)        ! mass fraction of aerosols activated [fraction]
    real(r8), intent(out) :: fluxn(:)     ! flux of activated aerosol number fraction into cloud [m/s]
    real(r8), intent(out) :: fluxm(:)     ! flux of activated aerosol mass fraction into cloud [m/s]
    real(r8), intent(out) :: flux_fullact ! flux of activated aerosol fraction assuming 100% activation [m/s]
    !---------------------------------------------------------------------------------
    ! flux_fullact is used for consistency check -- this should match (ekd(k)*zs(k))
    ! also, fluxm/flux_fullact gives fraction of aerosol mass flux
    ! that is activated
    !---------------------------------------------------------------------------------

    !optional
    real(r8), optional :: smax_prescribed  ! prescribed max. supersaturation for secondary activation [fraction]

    !local
    integer imode  !  mode index

    real(r8), parameter :: p0 = 1013.25e2_r8    ! reference pressure [Pa]
    real(r8), parameter :: smallest_val = 1.e-39_r8
    real(r8), parameter :: smaller_val = 1.e-20_r8
    real(r8), parameter :: small_val = 1.e-10_r8
    real(r8) pres ! pressure [Pa]
    real(r8) diff0 ! vapor diffusivity [m2/s]
    real(r8) conduct0 ! thermal conductivity [J / (m-s-K)]
    real(r8) es ! saturation vapor pressure [Pa]
    real(r8) qs ! water vapor saturation specific humidity [kg/kg]
    real(r8) dqsdt ! change in qs with temperature  [(kg/kg)/T]
    real(r8) zeta, eta(nmode) ! [unitless]
    real(r8) alpha ! [/m]
    real(r8) gamma ![m3/kg]
    real(r8) beta ! [m2/s]
    real(r8) gthermfac ! thermodynamic function [m2/s]
    real(r8) :: amcube(nmode) ! cube of dry mode radius [m3]
    real(r8) smc(nmode) ! critical supersaturation for number mode radius [fraction]
    real(r8) :: lnsm(nmode) ! ln(critical supersaturation for activation) [unitless]

    real(r8) wnuc  ! nucleation w, but = w_in if wdiab == 0 [m/s]
    real(r8) alw ! [/s]
    real(r8) smax ! maximum supersaturation [fraction]
    real(r8) lnsmax ! ln(smax) [unitless]
    real(r8) arg_erf_n,arg_erf_m  ! [unitless]
    real(r8) etafactor1  ! [/ s^(3/2)]
    real(r8) etafactor2(nmode),etafactor2max ! [s^(3/2)]


    !initialize activated aerosols and their fluxes to zero for all the modes
    fn(:)=0._r8
    fm(:)=0._r8
    fluxn(:)=0._r8
    fluxm(:)=0._r8
    flux_fullact=0._r8

    !return if aerosol number is negligible in the accumulation mode
    if(na(1) < smaller_val)return

    !return if vertical velocity is 0 or negative
    if(w_in <= 0._r8)return

    if ( present( smax_prescribed ) ) then
       !return if max supersaturation is 0 or negative
       if (smax_prescribed <= 0.0_r8) return
    endif

    pres=rair*rhoair*tair !pressure
    !Obtain Saturation vapor pressure (es) and saturation specific humidity (qs)
    call qsat(tair, pres, es, qs) !es and qs are the outputs

    dqsdt=latvap/(rh2o*tair*tair)*qs
    alpha=gravit*(latvap/(cpair*rh2o*tair*tair)-1._r8/(rair*tair))
    gamma=(1+latvap/cpair*dqsdt)/(rhoair*qs)
    etafactor2max=1.e10_r8/(alpha*wmaxf)**1.5_r8 !this should make eta big if na is very small.

    diff0=0.211e-4_r8*(p0/pres)*(tair/t0)**1.94_r8
    conduct0=(5.69_r8+0.017_r8*(tair-t0))*4.186e2_r8*1.e-5_r8 !convert to J/m/s/deg
    gthermfac=1._r8/(rhoh2o/(diff0*rhoair*qs)                                    &
         +latvap*rhoh2o/(conduct0*tair)*(latvap/(rh2o*tair)-1._r8)) !gthermfac is same for all modes
    beta=2._r8*pi*rhoh2o*gthermfac*gamma
    wnuc = w_in
    alw=alpha*wnuc
    etafactor1=alw*sqrt(alw)
    zeta=twothird*sqrt(alw)*aten/sqrt(gthermfac)

    !Here compute smc, eta for all modes for maxsat calculation
    do imode=1,nmode
       if(volume(imode) > smallest_val .and. na(imode) > smallest_val)then
          !number mode radius (m)
          amcube(imode)=(3._r8*volume(imode)/(4._r8*pi*exp45logsig(imode)*na(imode)))  ! only if variable size dist
          !Growth coefficent Abdul-Razzak & Ghan 1998 eqn 16
          !should depend on mean radius of mode to account for gas kinetic effects
          !see Fountoukis and Nenes, JGR2005 and Meskhidze et al., JGR2006
          !for approriate size to use for effective diffusivity.

          etafactor2(imode)=1._r8/(na(imode)*beta*sqrt(gthermfac))
          if(hygro(imode) > small_val)then
             smc(imode)=2._r8*aten*sqrt(aten/(27._r8*hygro(imode)*amcube(imode))) ! only if variable size dist
          else
             smc(imode)=100._r8
          endif
       else
          smc(imode)=1._r8
          etafactor2(imode)=etafactor2max ! this should make eta big if na is very small.
       endif
       lnsm(imode)=log(smc(imode)) ! only if variable size dist
       eta(imode)=etafactor1*etafactor2(imode)
    enddo

    !Find maximum supersaturation
    !Use smax_prescribed if it is present; otherwise get smax from subr maxsat
    if ( present( smax_prescribed ) ) then
       smax = smax_prescribed
    else
       call maxsat(zeta,eta,nmode,smc,smax)
    endif
    lnsmax=log(smax)

    !Use maximum supersaturation to calculate aerosol activation output
    do imode=1,nmode
       arg_erf_n=twothird*(lnsm(imode)-lnsmax)/(sq2*alogsig(imode))
       fn(imode)=0.5_r8*(1._r8-erf(arg_erf_n)) !activated number
       arg_erf_m=arg_erf_n-1.5_r8*sq2*alogsig(imode)
       fm(imode)=0.5_r8*(1._r8-erf(arg_erf_m)) !activated mass
       fluxn(imode)=fn(imode)*w_in !activated aerosol number flux
       fluxm(imode)=fm(imode)*w_in !activated aerosol mass flux
    enddo
    flux_fullact = w_in

  end subroutine activate_modal

  !===============================================================================

  subroutine maxsat(zeta,eta,nmode,smc,smax)

    !      calculates maximum supersaturation for multiple
    !      competing aerosol modes.

    !      Abdul-Razzak and Ghan, A parameterization of aerosol activation.
    !      2. Multiple aerosol types. J. Geophys. Res., 105, 6837-6844.

    ! input arguments
    integer,  intent(in)  :: nmode ! number of modes
    real(r8), intent(in)  :: smc(nmode) ! critical supersaturation for number mode radius [fraction]
    real(r8), intent(in)  :: zeta ! [dimensionless]
    real(r8), intent(in)  :: eta(nmode) ! [dimensionless]

    ! output arguments
    real(r8), intent(out) :: smax ! maximum supersaturation [fraction]

    ! local arguments
    integer  :: m  ! mode index

    real(r8), parameter :: smaller_val = 1.e-20_r8

    real(r8) :: sum, g1, g2

    logical  :: weak_forcing  !  whether forcing is sufficiently weak or not

    weak_forcing = .true.
    do m=1, nmode
       if(zeta > 1.e5_r8*eta(m).or.smc(m)*smc(m) > 1.e5_r8*eta(m))then
          !weak forcing. essentially none activated
          smax=smaller_val
       else
          !significant activation of this mode. calc activation of all modes.
          weak_forcing = .false.
          exit
       endif
    enddo

    !if the forcing is weak, return
    if (weak_forcing) return

    sum=0
    do m=1,nmode
       if(eta(m) > smaller_val)then
          g1 = (zeta/eta(m)) * sqrt(zeta/eta(m))
          g2 = (smc(m)/sqrt(eta(m)+3._r8*zeta)) * sqrt(smc(m)/sqrt(eta(m)+3._r8*zeta))
          sum=sum+(f1(m)*g1+f2(m)*g2)/(smc(m)*smc(m))
       else
          sum=1.e20_r8
       endif
    enddo
    smax=1._r8/sqrt(sum)

  end subroutine maxsat

  !===============================================================================

  subroutine ccncalc(state_q, tair, qcldbrn, qcldbrn_num, ncol, cs, ccn)
    ! calculates number concentration of aerosols activated as CCN at
    ! supersaturation supersat.
    ! assumes an internal mixture of a multiple externally-mixed aerosol modes
    ! cgs units

    ! Ghan et al., Atmos. Res., 1993, 198-221.

    ! input arguments
    real(r8), intent(in)  :: state_q(:,:,:) ! aerosol mmrs [kg/kg]
    real(r8), intent(in)  :: tair(:,:)     ! air temperature [K]
    real(r8), intent(in)  :: qcldbrn(:,:,:,:), qcldbrn_num(:,:,:) ! cloud-borne aerosol mass / number  mixing ratios [kg/kg or #/kg]
    integer, intent(in)   :: ncol  ! number of columns
    real(r8), intent(in)  :: cs(pcols,pver)       ! air density [kg/m3]

    ! output arguments
    real(r8), intent(out) :: ccn(pcols,pver,psat) ! number conc of aerosols activated at supersat [#/m3]

    ! local
    real(r8) :: naerosol(ntot_amode,pcols) ! interstit+activated aerosol number conc [#/m3]
    real(r8) :: vaerosol(ntot_amode,pcols) ! interstit+activated aerosol volume conc [m3/m3]
    real(r8) :: amcube(pcols)  ! [m3]
    real(r8) :: amcubecoef(ntot_amode) ! [dimensionless]
    real(r8) :: argfactor(ntot_amode)  ! [dimensionless]
    real(r8) :: surften_coef  ! [m-K]
    real(r8) :: aparam(pcols) ! surface tension parameter  [m]
    real(r8) :: hygro(ntot_amode,pcols)  ! aerosol hygroscopicity [dimensionless]
    real(r8) :: sm(pcols)  ! critical supersaturation at mode radius [fraction]
    real(r8) :: arg_erf_ccn ! [dimensionless]
    real(r8) :: smcoef(pcols)  ! [m^(3/2)]
    integer  :: imode       ! mode index
    integer  :: kk          ! level index
    integer  :: icol        ! column index
    integer  :: lsat        ! level of supersaturation
    integer  :: phase       ! phase of aerosol

    !     mathematical constants
    real(r8), parameter ::  percent_to_fraction = 0.01_r8
    real(r8), parameter ::  per_m3_to_per_cm3 = 1.e-6_r8
    real(r8), parameter ::  smcoefcoef=2._r8/sqrt(27._r8)
    real(r8), parameter ::  nconc_thresh = 1.e-3_r8

    real(r8) super(psat) ! supersaturation [fraction]
    !-------------------------------------------------------------------------------

    phase=3 ! interstitial+cloudborne

    super(:)=supersat(:)*percent_to_fraction

    surften_coef=2._r8*mwh2o*surften/(r_universal*rhoh2o)

    ccn(:,:,:) = 0._r8

    !  C++ CODE PORT NOTE:  C++ utility functions could potentially replace the following
    amcubecoef(:)=3._r8/(4._r8*pi*exp45logsig(:))

    argfactor(:)=twothird/(sq2*alogsig(:))

    do kk=top_lev,pver

       do icol = 1, ncol
          call loadaer( state_q(icol,kk,:), &  ! in
               cs(icol,kk), phase, &  ! in
               naerosol(:,icol), vaerosol(:,icol), hygro(:,icol), &  ! out
               qcldbrn(icol,kk,:,:), qcldbrn_num(icol,kk,:) )  ! optional in
       enddo

       aparam(1:ncol) = surften_coef/tair(1:ncol,kk)
       smcoef(1:ncol)=smcoefcoef*aparam(1:ncol)*sqrt(aparam(1:ncol))

       do imode=1,ntot_amode

          where(naerosol(imode,:ncol) > nconc_thresh)
             amcube(:ncol)=amcubecoef(imode)*vaerosol(imode,:ncol)/naerosol(imode,:ncol)
             sm(:ncol)=smcoef(:ncol)/sqrt(hygro(imode,:ncol)*amcube(:ncol)) ! critical supersaturation
          elsewhere
             sm(:ncol)=1._r8 ! value shouldn't matter much since naerosol is small
          endwhere

          do lsat=1,psat
             do icol=1,ncol
                arg_erf_ccn=argfactor(imode)*log(sm(icol)/super(lsat))
                ccn(icol,kk,lsat)=ccn(icol,kk,lsat)+naerosol(imode,icol)*0.5_r8*(1._r8-erf(arg_erf_ccn))
             enddo
          enddo

       enddo

    enddo

    ccn(:ncol,:,:)=ccn(:ncol,:,:)*per_m3_to_per_cm3 ! convert from #/m3 to #/cm3

  end subroutine ccncalc

  !===============================================================================

  subroutine loadaer( &
       state_q, cs, phase, &   ! in
       naerosol, vaerosol, hygro, &  ! out
       qcldbrn1d, qcldbrn1d_num) ! optional in

    use modal_aero_data,   only: maxd_aspectype

    ! return aerosol number, volume concentrations, and bulk hygroscopicity at one specific column and level

    ! input arguments
    real(r8), intent(in) :: state_q(:)        ! aerosol mmrs [kg/kg]
    real(r8), intent(in) :: cs          ! air density [kg/m3]
    integer,  intent(in) :: phase       ! phase of aerosol: 1 for interstitial, 2 for cloud-borne, 3 for sum

    ! output arguments
    real(r8), intent(out) :: naerosol(ntot_amode)  ! number conc [#/m3]
    real(r8), intent(out) :: vaerosol(ntot_amode)  ! volume conc [m3/m3]
    real(r8), intent(out) :: hygro(ntot_amode)     ! bulk hygroscopicity of mode [dimensionless]

    ! optional input arguments
    real(r8), intent(in), optional  :: qcldbrn1d(:,:), qcldbrn1d_num(:) ! cloud-borne aerosol mass / number mixing ratios [kg/kg or #/kg]

    ! internal
    real(r8), parameter :: even_smaller_val = 1.0e-30_r8

    real(r8) :: vaerosolsum(ntot_amode)  ! sum to find volume conc [m3/kg]
    real(r8) :: hygrosum(ntot_amode)     ! sum to bulk hygroscopicity of mode [m3/kg]
    real(r8) :: qcldbrn_local(maxd_aspectype,ntot_amode)  ! local cloud-borne aerosol mass mixing ratios [kg/kg]
    real(r8) :: qcldbrn_num_local(ntot_amode) ! local cloud-borne aerosol number mixing ratios [#/kg]

    integer :: imode       ! mode index
    !-------------------------------------------------------------------------------

    !Currenly supports only phase 1 (interstitial) and 3 (interstitial+cldbrn)
    if (phase /= 1 .and. phase /=3) then
       write(iulog,*)'phase=',phase,' in loadaer'
       call endrun('phase error in loadaer')
    endif

    qcldbrn_local(:,:) = 0._r8
    if(present(qcldbrn1d)) qcldbrn_local = qcldbrn1d

    vaerosolsum(:) = 0._r8
    hygrosum(:)    = 0._r8

    !Sum over all species within imode to get bulk hygroscopicity and volume conc
    !phase == 1 is interstitial only.
    !phase == 3 is interstitial + cldborne
    !  Assumes iphase =1 or 3, so interstitial is always summed, added with cldbrn when present
    !  iphase = 2 would require alternate logic from following subroutine

    do imode=1,ntot_amode
       call get_aer_mmr_sum(imode, nspec_amode(imode), state_q(:), qcldbrn_local(:,imode), & !in
            vaerosolsum(imode), hygrosum(imode))    !inout
    enddo

    !  Finalize computation of bulk hygrospopicity and volume conc

    where (vaerosolsum(:ntot_amode) > even_smaller_val )
       hygro(:ntot_amode) = hygrosum(:ntot_amode) / vaerosolsum(:ntot_amode)
       vaerosol(:ntot_amode) = vaerosolsum(:ntot_amode) * cs
    elsewhere
       hygro(:ntot_amode) = 0.0_r8
       vaerosol(:ntot_amode) = 0.0_r8
    endwhere

    qcldbrn_num_local(:) = 0._r8
    if(present(qcldbrn1d_num)) qcldbrn_num_local = qcldbrn1d_num

    ! Compute aerosol number concentration
    do imode=1,ntot_amode
       call get_aer_num(imode, state_q(:), cs, vaerosol(imode), qcldbrn_num_local(imode), &!in
            naerosol(imode)) !out
    enddo

  end subroutine loadaer
  !===============================================================================

  subroutine get_aer_mmr_sum(imode, nspec, state_q, qcldbrn1d, & !in
       vaerosolsum, hygrosum)   !inout

    ! input arguments
    integer,  intent(in) :: imode       ! mode index
    integer,  intent(in) :: nspec       ! total # of species in mode imode
    real(r8), intent(in) :: state_q(:) ! interstitial aerosol mass mixing ratios [kg/kg]
    real(r8), intent(in) :: qcldbrn1d(:) ! cloud-borne aerosol mass mixing ratios [kg/kg]

    ! in/out arguments
    real(r8), intent(inout) :: vaerosolsum  ! sum to find volume conc [m3/kg]
    real(r8), intent(inout) :: hygrosum   ! sum to bulk hygroscopicity of mode [m3/kg]

    ! internal
    real(r8) :: density_sp ! density at species / mode indices [kg/m3]
    real(r8) :: hygro_sp   ! hygroscopicity at species / mode indices [dimensionless]
    real(r8) :: vol     !aerosol volume mixing ratio [m3/kg]

    integer  :: lspec, spc_idx, type_idx

    !Start to compute bulk volume conc / hygroscopicity by summing over species per mode.
    do lspec = 1, nspec
       type_idx = lspectype_amode(lspec,imode)
       density_sp  = specdens_amode(type_idx) !species density
       hygro_sp    = spechygro(type_idx)      !species hygroscopicity
       spc_idx   = lmassptr_amode(lspec,imode) !index of species in state_q array
       vol = max(state_q(spc_idx) + qcldbrn1d(lspec), 0._r8)/density_sp !volume = mmr/density
       vaerosolsum = vaerosolsum + vol        !bulk volume
       hygrosum    = hygrosum + vol*hygro_sp !bulk hygroscopicity
    enddo

  end subroutine get_aer_mmr_sum

  !===============================================================================

  subroutine get_aer_num(imode, state_q, cs, vaerosol, qcldbrn1d_num, &!in
       naerosol) !out

    use mam_support, only: min_max_bound

    ! input arguments
    integer,  intent(in) :: imode        ! mode index
    real(r8), intent(in) :: state_q(:) ! interstitial aerosol number mixing ratios [#/kg]
    real(r8), intent(in) :: cs           ! air density [kg/m3]
    real(r8), intent(in) :: vaerosol  ! volume conc [m3/m3]
    real(r8), intent(in) :: qcldbrn1d_num ! cloud-borne aerosol number mixing ratios [#/kg]

    !output arguments
    real(r8), intent(out) :: naerosol  ! number conc [#/m3]

    !internal
    integer  :: num_idx

    !convert number mixing ratios to number concentrations
    !Use bulk volume conc found previously to bound value

    num_idx = numptr_amode(imode)
    naerosol = (state_q(num_idx) + qcldbrn1d_num)*cs
    !adjust number so that dgnumlo < dgnum < dgnumhi
    naerosol = min_max_bound(vaerosol*voltonumbhi_amode(imode),vaerosol*voltonumblo_amode(imode),naerosol)

  end subroutine get_aer_num

  !===============================================================================

end module ndrop
