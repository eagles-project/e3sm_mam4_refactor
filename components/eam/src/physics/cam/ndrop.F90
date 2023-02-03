
module ndrop

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



implicit none
private
save

public ndrop_init, dropmixnuc, activate_modal

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
   state, ptend, dtmicro, pbuf, wsub, &
   cldn, cldo, tendnd, factnum)

   ! vertical diffusion and nucleation of cloud droplets
   ! assume cloud presence controlled by cloud fraction
   ! doesn't distinguish between warm, cold clouds

   use output_aerocom_aie , only: do_aerocom_ind3
   use modal_aero_data,   only: lmassptrcw_amode, numptrcw_amode, qqcw_get_field,maxd_aspectype

   ! arguments
   type(physics_state), target, intent(in)    :: state
   type(physics_ptend),         intent(out)   :: ptend
   real(r8),                    intent(in)    :: dtmicro     ! time step for microphysics (s)

   type(physics_buffer_desc), pointer :: pbuf(:)

   ! arguments
   real(r8), intent(in) :: wsub(pcols,pver)    ! subgrid vertical velocity
   real(r8), intent(in) :: cldn(pcols,pver)    ! cloud fraction
   real(r8), intent(in) :: cldo(pcols,pver)    ! cloud fraction on previous time step

   ! output arguments
   real(r8), intent(out) :: tendnd(pcols,pver) ! tendency in droplet number mixing ratio [#/kg/s]
   real(r8), intent(out) :: factnum(:,:,:)     ! activation fraction for aerosol number
   !--------------------Local storage-------------------------------------

   integer  :: lchnk               ! chunk identifier
   integer  :: ncol                ! number of columns
!  BJG not necessary   integer  :: loop_up_bnd
   real(r8), pointer :: ncldwtr(:,:) ! droplet number mixing ratio [#/kg]
   real(r8), pointer :: temp(:,:)    ! temperature [K]
! BJG not used   real(r8), pointer :: omega(:,:)   ! vertical velocity [Pa/s]
   real(r8), pointer :: pmid(:,:)    ! mid-level pressure [Pa]
   real(r8), pointer :: pint(:,:)    ! pressure at layer interfaces [Pa]
   real(r8), pointer :: pdel(:,:)    ! pressure thickess of layer [Pa]
   real(r8), pointer :: rpdel(:,:)   ! inverse of pressure thickess of layer [/Pa]
   real(r8), pointer :: zm(:,:)      ! geopotential height of level [m]
   real(r8), pointer :: state_q(:,:,:)        ! aerosol mmrs [kg/kg]

   real(r8), pointer :: kvh(:,:)     ! vertical diffusivity (m2/s)

   type(ptr2d_t), allocatable :: raer(:)     ! aerosol mass, number mixing ratios
   type(ptr2d_t), allocatable :: qqcw(:)
   real(r8) :: raertend(pver)  ! tendency of aerosol mass, number mixing ratios
   real(r8) :: qqcwtend(pver)  ! tendency of cloudborne aerosol mass, number mixing ratios


   real(r8), parameter :: zkmin = 0.01_r8, zkmax = 100._r8  ! min, max vertical diffusivity [m2/s]
   real(r8), parameter :: wmixmin = 0.1_r8        ! minimum turbulence vertical velocity [m/s]
   real(r8) :: sq2pi

! BJG i,k,l,m,n replaced by icol, kk, lspec/lsat,imode,isub which are defined below  integer  :: i, k, l, m, mm, n
   integer  :: mm
   integer  :: km1, kp1
   integer  :: nnew, nsav, ntemp
   integer  :: lptr
   integer  :: nsubmix, nsubmix_bnd
   integer, save :: count_submix(100)
   integer  :: phase ! phase of aerosol

   real(r8) :: arg
   real(r8) :: dtinv
   real(r8) :: dtmin, tinv, dtt

   real(r8) :: zs(pver) ! inverse of distance between levels (m)
   real(r8) :: qcld(pver) ! cloud droplet number mixing ratio (#/kg)
   real(r8) :: qncld(pver)     ! droplet number nucleated on cloud boundaries
   real(r8) :: srcn(pver)       ! droplet source rate (/s)
   real(r8) :: cs(pcols,pver)      ! air density (kg/m3)
   real(r8) :: csbot(pver)       ! air density at bottom (interface) of layer (kg/m3)
   real(r8) :: csbot_cscen(pver) ! csbot(i)/cs(i,k)
   real(r8) :: dz(pcols,pver)      ! geometric thickness of layers (m)

   real(r8) :: wtke(pcols,pver)     ! turbulent vertical velocity at base of layer k (m/s)
   real(r8) :: wtke_cen(pcols,pver) ! turbulent vertical velocity at center of layer k (m/s)
!  BJG  wmix, wmin not used or zero  real(r8) :: wbar, wmix, wmin, wmax
   real(r8) :: wbar, wmax

   real(r8) :: zn(pver)   ! g/pdel (m2/g) for layer
! BJG not used   real(r8) :: flxconv    ! convergence of flux into lowest layer

! BJG not used or zero   real(r8) :: wdiab           ! diabatic vertical velocity
   real(r8) :: ekd(pver)       ! diffusivity for droplets (m2/s)
   real(r8) :: ekk(0:pver)     ! density*diffusivity for droplets (kg/m3 m2/s)
   real(r8) :: ekkp(pver)      ! zn*zs*density*diffusivity
   real(r8) :: ekkm(pver)      ! zn*zs*density*diffusivity

   real(r8) :: dum, dumc
   real(r8) :: tmpa
   real(r8) :: dact
   real(r8) :: fluxntot         ! (#/cm2/s)
   real(r8) :: dtmix
!  BJG not used    real(r8) :: alogarg
   real(r8) :: overlapp(pver), overlapm(pver) ! cloud overlap

   real(r8) :: nsource(pcols,pver)            ! droplet number source (#/kg/s)
   real(r8) :: ndropmix(pcols,pver)           ! droplet number mixing (#/kg/s)
   real(r8) :: ndropcol(pcols)               ! column droplet number (#/m2)
   real(r8) :: cldo_tmp, cldn_tmp
! BJG not used   real(r8) :: tau_cld_regenerate
! BJG not currently nonzero   real(r8) :: taumix_internal_pver_inv ! 1/(internal mixing time scale for k=pver) (1/s)

   real(r8), allocatable :: nact(:,:)  ! fractional aero. number  activation rate (/s)
   real(r8), allocatable :: mact(:,:)  ! fractional aero. mass    activation rate (/s)

   real(r8), allocatable :: raercol(:,:,:)    ! single column of aerosol mass, number mixing ratios
   real(r8), allocatable :: raercol_cw(:,:,:) ! same as raercol but for cloud-borne phase


   real(r8) :: na(pcols), va(pcols), hy(pcols)
   real(r8), allocatable :: naermod(:)  ! (1/m3)
   real(r8), allocatable :: hygro(:)    ! hygroscopicity of aerosol mode
   real(r8), allocatable :: vaerosol(:) ! interstit+activated aerosol volume conc (cm3/cm3)

   real(r8) :: source(pver)

   real(r8), allocatable :: fn(:)              ! activation fraction for aerosol number
   real(r8), allocatable :: fm(:)              ! activation fraction for aerosol mass

   real(r8), allocatable :: fluxn(:)           ! number  activation fraction flux (cm/s)
   real(r8), allocatable :: fluxm(:)           ! mass    activation fraction flux (cm/s)
   real(r8)              :: flux_fullact(pver) ! 100%    activation fraction flux (cm/s)
   !     note:  activation fraction fluxes are defined as
   !     fluxn = [flux of activated aero. number into cloud (#/cm2/s)]
   !           / [aero. number conc. in updraft, just below cloudbase (#/cm3)]


   real(r8), allocatable :: coltend(:,:)       ! column tendency for diagnostic output
   real(r8), allocatable :: coltend_cw(:,:)    ! column tendency
   real(r8) :: ccn(pcols,pver,psat)    ! number conc of aerosols activated at supersat
   integer :: ccn3d_idx
   real(r8), pointer :: ccn3d(:, :)

!+++ AeroCOM IND3 output
   real(r8) :: ccn3col(pcols), ccn4col(pcols)
   real(r8) :: ccn3bl(pcols), ccn4bl(pcols)
   real(r8) :: zi2(pver+1), zm2(pver)
   integer  :: idx1000
   logical  :: zmflag

  integer :: icol        ! column index
  integer :: imode       ! mode index
  integer :: kk          ! level index
  integer :: lspec      ! species index for given mode
  integer :: lsat       !  level of supersaturation
  integer :: isub       ! substep index
  integer :: spc_idx, num_idx
  real(r8) :: qcldbrn(pcols,maxd_aspectype,pver,ntot_amode) ! ! cloud-borne aerosol mass mixing ratios [kg/kg]
  real(r8) :: qcldbrn_num(pcols,pver,ntot_amode) ! ! cloud-borne aerosol number mixing ratios [#/kg]
  real(r8), pointer :: fldcw(:,:)           !specie mmr/num (cloud borne)

   !-------------------------------------------------------------------------------

   sq2pi = sqrt(2._r8*pi)

   lchnk = state%lchnk
   ncol  = state%ncol

   ncldwtr  => state%q(:,:,numliq_idx)
   temp     => state%t
!  BJG   omega    => state%omega
   pmid     => state%pmid
   pint     => state%pint
   pdel     => state%pdel
   rpdel    => state%rpdel
   zm       => state%zm
    
   state_q  => state%q

   call pbuf_get_field(pbuf, kvh_idx, kvh)

!  BJG:  below probably always false?
   if(do_aerocom_ind3) then
       ccn3d_idx = pbuf_get_index('ccn3d')
       call pbuf_get_field(pbuf, ccn3d_idx, ccn3d)
   end if


   arg = 1.0_r8
   if (abs(0.8427_r8 - erf(arg))/0.8427_r8 > 0.001_r8) then
      write(iulog,*) 'erf(1.0) = ',ERF(arg)
      call endrun('dropmixnuc: Error function error')
   endif
   arg = 0.0_r8
   if (erf(arg) /= 0.0_r8) then
      write(iulog,*) 'erf(0.0) = ',erf(arg)
      write(iulog,*) 'dropmixnuc: Error function error'
      call endrun('dropmixnuc: Error function error')
   endif

   dtinv = 1._r8/dtmicro

   allocate( &
      nact(pver,ntot_amode),          &
      mact(pver,ntot_amode),          &
      raer(ncnst_tot),                &
      qqcw(ncnst_tot),                &
      raercol(pver,ncnst_tot,2),      &
      raercol_cw(pver,ncnst_tot,2),   &
      coltend(pcols,ncnst_tot),       &
      coltend_cw(pcols,ncnst_tot),    &
      naermod(ntot_amode),            &
      hygro(ntot_amode),              &
      vaerosol(ntot_amode),           &
      fn(ntot_amode),                 &
      fm(ntot_amode),                 &
      fluxn(ntot_amode),              &
      fluxm(ntot_amode)               )

   ! Init pointers to mode number and specie mass mixing ratios in
   ! intersitial and cloud borne phases.
   do imode = 1, ntot_amode
      mm = mam_idx(imode, 0)
      call rad_cnst_get_mode_num(0, imode, 'a', state, pbuf, raer(mm)%fld)
      call rad_cnst_get_mode_num(0, imode, 'c', state, pbuf, qqcw(mm)%fld)  ! cloud-borne aerosol
      do lspec = 1, nspec_amode(imode)
         mm = mam_idx(imode, lspec)
         call rad_cnst_get_aer_mmr(0, imode, lspec, 'a', state, pbuf, raer(mm)%fld)
         call rad_cnst_get_aer_mmr(0, imode, lspec, 'c', state, pbuf, qqcw(mm)%fld)  ! cloud-borne aerosol
      end do
   end do

   factnum = 0._r8
   wtke    = 0._r8

!  BJG:  below probably always true?
   if (prog_modal_aero) then
      ! aerosol tendencies
      call physics_ptend_init(ptend, state%psetcols, 'ndrop_aero', lq=lq)
   else
      ! no aerosol tendencies
      call physics_ptend_init(ptend, state%psetcols, 'ndrop')
   end if

   !initialize variables to zero
   ndropmix(:,:) = 0._r8
   nsource(:,:) = 0._r8

   ! overall_main_i_loop
   do icol = 1, ncol

      do kk = top_lev, pver-1
         zs(kk) = 1._r8/(zm(icol,kk) - zm(icol,kk+1))
      end do
      zs(pver) = zs(pver-1)

      ! load number nucleated into qcld on cloud boundaries

      do kk = top_lev, pver

         qcld(kk)  = ncldwtr(icol,kk)
         qncld(kk) = 0._r8
         srcn(kk)  = 0._r8
         cs(icol,kk)  = pmid(icol,kk)/(rair*temp(icol,kk))        ! air density (kg/m3)
         dz(icol,kk)  = 1._r8/(cs(icol,kk)*gravit*rpdel(icol,kk)) ! layer thickness in m

         do imode = 1, ntot_amode
            nact(kk,imode) = 0._r8
            mact(kk,imode) = 0._r8
         end do

         zn(kk) = gravit*rpdel(icol,kk)

         if (kk < pver) then
            ekd(kk)   = kvh(icol,kk+1)
            ekd(kk)   = max(ekd(kk), zkmin)
            ekd(kk)   = min(ekd(kk), zkmax)
            csbot(kk) = 2.0_r8*pint(icol,kk+1)/(rair*(temp(icol,kk) + temp(icol,kk+1)))
            csbot_cscen(kk) = csbot(kk)/cs(icol,kk)
         else
            ekd(kk)   = 0._r8
            csbot(kk) = cs(icol,kk)
            csbot_cscen(kk) = 1.0_r8
         end if

         ! rce-comment - define wtke at layer centers for new-cloud activation
         !    and at layer boundaries for old-cloud activation
         !++ag
         wtke_cen(icol,kk) = wsub(icol,kk)
         wtke(icol,kk)     = wsub(icol,kk)
         !--ag
         wtke_cen(icol,kk) = max(wtke_cen(icol,kk), wmixmin)
         wtke(icol,kk)     = max(wtke(icol,kk), wmixmin)

         nsource(icol,kk) = 0._r8

      end do

      nsav = 1
      nnew = 2
      do imode = 1, ntot_amode
         mm = mam_idx(imode,0)
         raercol_cw(:,mm,nsav) = 0.0_r8
         raercol(:,mm,nsav)    = 0.0_r8
         raercol_cw(top_lev:pver,mm,nsav) = qqcw(mm)%fld(icol,top_lev:pver)
         raercol(top_lev:pver,mm,nsav)    = raer(mm)%fld(icol,top_lev:pver)
         do lspec = 1, nspec_amode(imode)
            mm = mam_idx(imode,lspec)
            raercol_cw(top_lev:pver,mm,nsav) = qqcw(mm)%fld(icol,top_lev:pver)
            raercol(top_lev:pver,mm,nsav)    = raer(mm)%fld(icol,top_lev:pver)
         end do
      end do

      ! droplet nucleation/aerosol activation

      ! tau_cld_regenerate = time scale for regeneration of cloudy air
      !    by (horizontal) exchange with clear air
!  BJG      tau_cld_regenerate = 3600.0_r8 * 3.0_r8

      ! k-loop for growing/shrinking cloud calcs .............................
      ! grow_shrink_main_k_loop: &
      do kk = top_lev, pver

         ! shrinking cloud ......................................................
         !    treat the reduction of cloud fraction from when cldn(i,k) < cldo(i,k)
         !    and also dissipate the portion of the cloud that will be regenerated
         cldo_tmp = cldo(icol,kk)

! BJG         if(regen_fix) then
            cldn_tmp = cldn(icol,kk) !* exp( -dtmicro/tau_cld_regenerate )!HW: there is a bug here; turn off regeneration,01/10/2012
! BJG         else
! BJG            cldn_tmp = cldn(i,k) * exp( -dtmicro/tau_cld_regenerate )
! BJG         endif
         !    alternate formulation
         !    cldn_tmp = cldn(i,k) * max( 0.0_r8, (1.0_r8-dtmicro/tau_cld_regenerate) )

         if (cldn_tmp < cldo_tmp) then
            !  droplet loss in decaying cloud
            !++ sungsup
            nsource(icol,kk) = nsource(icol,kk) + qcld(kk)*(cldn_tmp - cldo_tmp)/cldo_tmp*dtinv
            qcld(kk)      = qcld(kk)*(1._r8 + (cldn_tmp - cldo_tmp)/cldo_tmp)
            !-- sungsup

            ! convert activated aerosol to interstitial in decaying cloud

            dumc = (cldn_tmp - cldo_tmp)/cldo_tmp
            do imode = 1, ntot_amode
               mm = mam_idx(imode,0)
               dact   = raercol_cw(kk,mm,nsav)*dumc
               raercol_cw(kk,mm,nsav) = raercol_cw(kk,mm,nsav) + dact   ! cloud-borne aerosol
               raercol(kk,mm,nsav)    = raercol(kk,mm,nsav) - dact
               do lspec = 1, nspec_amode(imode)
                  mm = mam_idx(imode,lspec)
                  dact    = raercol_cw(kk,mm,nsav)*dumc
                  raercol_cw(kk,mm,nsav) = raercol_cw(kk,mm,nsav) + dact  ! cloud-borne aerosol
                  raercol(kk,mm,nsav)    = raercol(kk,mm,nsav) - dact
               end do
            end do
         end if

         ! growing cloud ......................................................
         !    treat the increase of cloud fraction from when cldn(i,k) > cldo(i,k)
         !    and also regenerate part of the cloud
! BJG         if(regen_fix) then
            cldo_tmp = cldo(icol,kk)! HW turned off the regeneration growing
! BJG         else
! BJG            cldo_tmp = cldn_tmp
! BJG         endif
         cldn_tmp = cldn(icol,kk)

         if (cldn_tmp-cldo_tmp > 0.01_r8) then

            ! rce-comment - use wtke at layer centers for new-cloud activation
            wbar  = wtke_cen(icol,kk)
! BJG            wmix  = 0._r8
! BJG            wmin  = 0._r8
            wmax  = 10._r8
! BJG            wdiab = 0

            ! load aerosol properties, assuming external mixtures


            phase = 1 ! interstitial
            do imode = 1, ntot_amode


               call loadaer( &
                  state_q, icol, icol, kk, &
                  imode, nspec_amode(imode), cs, phase, na, va, &
                  hy )
               naermod(imode)  = na(icol)
               vaerosol(imode) = va(icol)
               hygro(imode)    = hy(icol)
            end do

            call activate_modal( &
               wbar, wmax,                       &
               temp(icol,kk), cs(icol,kk), naermod, ntot_amode, &
               vaerosol, hygro, fn, fm, fluxn,                      &
               fluxm,flux_fullact(kk))

            factnum(icol,kk,:) = fn

            dumc = (cldn_tmp - cldo_tmp)
            do imode = 1, ntot_amode
               mm = mam_idx(imode,0)
               dact   = dumc*fn(imode)*raer(mm)%fld(icol,kk) ! interstitial only
               qcld(kk) = qcld(kk) + dact
               nsource(icol,kk) = nsource(icol,kk) + dact*dtinv
               raercol_cw(kk,mm,nsav) = raercol_cw(kk,mm,nsav) + dact  ! cloud-borne aerosol
               raercol(kk,mm,nsav)    = raercol(kk,mm,nsav) - dact
               dum = dumc*fm(imode)
               do lspec = 1, nspec_amode(imode)
                  mm = mam_idx(imode,lspec)
                  dact    = dum*raer(mm)%fld(icol,kk) ! interstitial only
                  raercol_cw(kk,mm,nsav) = raercol_cw(kk,mm,nsav) + dact  ! cloud-borne aerosol
                  raercol(kk,mm,nsav)    = raercol(kk,mm,nsav) - dact
               enddo
            enddo
         endif

      enddo  ! grow_shrink_main_k_loop
      ! end of k-loop for growing/shrinking cloud calcs ......................

      ! ......................................................................
      ! start of k-loop for calc of old cloud activation tendencies ..........
      !
      ! rce-comment
      !    changed this part of code to use current cloud fraction (cldn) exclusively
      !    consider case of cldo(:)=0, cldn(k)=1, cldn(k+1)=0
      !    previous code (which used cldo below here) would have no cloud-base activation
      !       into layer k.  however, activated particles in k mix out to k+1,
      !       so they are incorrectly depleted with no replacement

      ! old_cloud_main_k_loop
! BJG      if(regen_fix) then
! BJG         loop_up_bnd = pver - 1
! BJG      else
! BJG         loop_up_bnd = pver
! BJG      endif
! BJG      do kk = top_lev, loop_up_bnd!pver
      do kk = top_lev, pver - 1
         kp1 = min0(kk+1, pver)
! BJG         taumix_internal_pver_inv = 0.0_r8

         if (cldn(icol,kk) > 0.01_r8) then

! BJG            wdiab = 0
!  BJG            wmix  = 0._r8                       ! single updraft
            wbar  = wtke(icol,kk)                   ! single updraft
! BJG            if (k == pver) wbar = wtke_cen(i,k) ! single updraft
            wmax  = 10._r8
! BJG            wmin  = 0._r8

! BJG            if (cldn(i,k) - cldn(i,kp1) > 0.01_r8 .or. k == pver) then
            if (cldn(icol,kk) - cldn(icol,kp1) > 0.01_r8 ) then

               ! cloud base

               ! ekd(k) = wtke(i,k)*dz(i,k)/sq2pi
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
               ekd(kk) = wbar/zs(kk)

! BJG               alogarg = max(1.e-20_r8, 1/cldn(icol,kk) - 1._r8)
!  BJG:  wmix = 0.    wmin    = wbar + wmix*0.25_r8*sq2pi*log(alogarg)
! BJG wmin not used               wmin    = wbar 
               phase   = 1   ! interstitial

               do imode = 1, ntot_amode
                  ! rce-comment - use kp1 here as old-cloud activation involves
                  !   aerosol from layer below
                  call loadaer( &
                     state_q, icol, icol, kp1,  &
                     imode, nspec_amode(imode), cs, phase, na, va,   &
                     hy )
                  naermod(imode)  = na(icol)
                  vaerosol(imode) = va(icol)
                  hygro(imode)    = hy(icol)
               end do

               call activate_modal( &
                  wbar, wmax,                       &
                  temp(icol,kk), cs(icol,kk), naermod, ntot_amode, &
                  vaerosol, hygro, fn, fm, fluxn,                      &
                  fluxm, flux_fullact(kk))

               factnum(icol,kk,:) = fn

! BJG               if (k < pver) then
                  dumc = cldn(icol,kk) - cldn(icol,kp1)
! BJG               else
! BJG                 if(regen_fix) then
! BJG                     dumc=0._r8
! BJG                 else
! BJG                    dumc = cldn(i,k)
! BJG                 endif
! BJG               endif

               fluxntot = 0

               ! rce-comment 1
               !    flux of activated mass into layer k (in kg/m2/s)
               !       = "actmassflux" = dumc*fluxm*raercol(kp1,lmass)*csbot(k)
               !    source of activated mass (in kg/kg/s) = flux divergence
               !       = actmassflux/(cs(i,k)*dz(i,k))
               !    so need factor of csbot_cscen = csbot(k)/cs(i,k)
               !                   dum=1./(dz(i,k))
               dum=csbot_cscen(kk)/(dz(icol,kk))

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
! BJG               if (k == pver) then
! BJG                  taumix_internal_pver_inv = flux_fullact(k)/dz(i,k)
! BJG               end if

               do imode = 1, ntot_amode
                  mm = mam_idx(imode,0)
                  fluxn(imode) = fluxn(imode)*dumc
                  fluxm(imode) = fluxm(imode)*dumc
                  nact(kk,imode) = nact(kk,imode) + fluxn(imode)*dum
                  mact(kk,imode) = mact(kk,imode) + fluxm(imode)*dum
! BJG                  if (k < pver) then
                     ! note that kp1 is used here
                     fluxntot = fluxntot &
                        + fluxn(imode)*raercol(kp1,mm,nsav)*cs(icol,kk)
!BJG                  else
! BJG                     tmpa = raercol(kp1,mm,nsav)*fluxn(m) &
! BJG                          + raercol_cw(kp1,mm,nsav)*(fluxn(m) &
! BJG                          - taumix_internal_pver_inv*dz(i,k))
! BJG                     fluxntot = fluxntot + max(0.0_r8, tmpa)*cs(i,k)
!BJG                  end if
               end do
               srcn(kk)      = srcn(kk) + fluxntot/(cs(icol,kk)*dz(icol,kk))
               nsource(icol,kk) = nsource(icol,kk) + fluxntot/(cs(icol,kk)*dz(icol,kk))

! BJG            endif  ! (cldn(i,k) - cldn(i,kp1) > 0.01 .or. k == pver)
            endif  ! (cldn(i,k) - cldn(i,kp1) > 0.01)

         else

            ! no cloud

            nsource(icol,kk) = nsource(icol,kk) - qcld(kk)*dtinv
            qcld(kk)      = 0

            ! convert activated aerosol to interstitial in decaying cloud

            do imode = 1, ntot_amode
               mm = mam_idx(imode,0)
               raercol(kk,mm,nsav)    = raercol(kk,mm,nsav) + raercol_cw(kk,mm,nsav)  ! cloud-borne aerosol
               raercol_cw(kk,mm,nsav) = 0._r8

               do lspec = 1, nspec_amode(imode)
                  mm = mam_idx(imode,lspec)
                  raercol(kk,mm,nsav)    = raercol(kk,mm,nsav) + raercol_cw(kk,mm,nsav) ! cloud-borne aerosol
                  raercol_cw(kk,mm,nsav) = 0._r8
               end do
            end do
         end if

      end do  ! old_cloud_main_k_loop

      ! switch nsav, nnew so that nnew is the updated aerosol
      ntemp = nsav
      nsav  = nnew
      nnew  = ntemp

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
      end do

      do kk = top_lev, pver
         km1     = max0(kk-1, top_lev)
         ekkp(kk) = zn(kk)*ekk(kk)*zs(kk)
         ekkm(kk) = zn(kk)*ekk(kk-1)*zs(km1)
         tinv    = ekkp(kk) + ekkm(kk)

         ! rce-comment -- tinv is the sum of all first-order-loss-rates
         !    for the layer.  for most layers, the activation loss rate
         !    (for interstitial particles) is accounted for by the loss by
         !    turb-transfer to the layer above.
         !    k=pver is special, and the loss rate for activation within
         !    the layer must be added to tinv.  if not, the time step
         !    can be too big, and explmix can produce negative values.
         !    the negative values are reset to zero, resulting in an
         !    artificial source.
! BJG always 0        if (kk == pver) tinv = tinv + taumix_internal_pver_inv
! BJG         if (kk == pver) tinv = tinv 

         if (tinv .gt. 1.e-6_r8) then
            dtt   = 1._r8/tinv
            dtmin = min(dtmin, dtt)
         end if
      end do

      dtmix   = 0.9_r8*dtmin
      nsubmix = dtmicro/dtmix + 1
      if (nsubmix > 100) then
         nsubmix_bnd = 100
      else
         nsubmix_bnd = nsubmix
      end if
      count_submix(nsubmix_bnd) = count_submix(nsubmix_bnd) + 1
      dtmix = dtmicro/nsubmix

      do kk = top_lev, pver
         kp1 = min(kk+1, pver)
         km1 = max(kk-1, top_lev)
         ! maximum overlap assumption
         if (cldn(icol,kp1) > 1.e-10_r8) then
            overlapp(kk) = min(cldn(icol,kk)/cldn(icol,kp1), 1._r8)
         else
            overlapp(kk) = 1._r8
         end if
         if (cldn(icol,km1) > 1.e-10_r8) then
            overlapm(kk) = min(cldn(icol,kk)/cldn(icol,km1), 1._r8)
         else
            overlapm(kk) = 1._r8
         end if
      end do


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
         end do
      end do


      ! old_cloud_nsubmix_loop
      do isub = 1, nsubmix
         qncld(:) = qcld(:)
         ! switch nsav, nnew so that nsav is the updated aerosol
         ntemp   = nsav
         nsav    = nnew
         nnew    = ntemp
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
! BJG always zero                 + raercol_cw(pver,mm,nsav)*(nact(pver,imode) - taumix_internal_pver_inv)
                 + raercol_cw(pver,mm,nsav)*nact(pver,imode)
            srcn(pver) = srcn(pver) + max(0.0_r8,tmpa)
         end do
         call explmix(  &
            qcld, srcn, ekkp, ekkm, overlapp,  &
            overlapm, qncld, pver, &
            dtmix, .false.)

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
! BJG always 0                 + raercol_cw(pver,mm,nsav)*(nact(pver,imode) - taumix_internal_pver_inv)
                 + raercol_cw(pver,mm,nsav)*nact(pver,imode)
            source(pver) = max(0.0_r8, tmpa)
!  BJG not used            flxconv = 0._r8

            call explmix( &
               raercol_cw(:,mm,nnew), source, ekkp, ekkm, overlapp, &
               overlapm, raercol_cw(:,mm,nsav), pver,   &
               dtmix, .false.)

            call explmix( &
               raercol(:,mm,nnew), source, ekkp, ekkm, overlapp,  &
               overlapm, raercol(:,mm,nsav), pver, &
               dtmix, .true., raercol_cw(:,mm,nsav))

            do lspec = 1, nspec_amode(imode)
               mm = mam_idx(imode,lspec)
               ! rce-comment -   activation source in layer k involves particles from k+1
               !	          source(:)= mact(:,m)*(raercol(:,mm,nsav))
               source(top_lev:pver-1) = mact(top_lev:pver-1,imode)*(raercol(top_lev+1:pver,mm,nsav))
               ! rce-comment- new formulation for k=pver
               !                 source(  pver  )= mact(  pver  ,m)*(raercol(  pver,mm,nsav))
               tmpa = raercol(pver,mm,nsav)*mact(pver,imode) &
!  BJG always 0                    + raercol_cw(pver,mm,nsav)*(mact(pver,imode) - taumix_internal_pver_inv)
                    + raercol_cw(pver,mm,nsav)*mact(pver,imode)
               source(pver) = max(0.0_r8, tmpa)
! BJG not used               flxconv = 0._r8

               call explmix( &
                  raercol_cw(:,mm,nnew), source, ekkp, ekkm, overlapp, &
                  overlapm, raercol_cw(:,mm,nsav), pver,   &
                  dtmix, .false.)

               call explmix( &
                  raercol(:,mm,nnew), source, ekkp, ekkm, overlapp,  &
                  overlapm, raercol(:,mm,nsav), pver, &
                  dtmix, .true., raercol_cw(:,mm,nsav))

            end do
         end do

      end do ! old_cloud_nsubmix_loop

      ! evaporate particles again if no cloud

      do kk = top_lev, pver
         if (cldn(icol,kk) == 0._r8) then
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
               end do
            end do
         end if
      end do

      ! droplet number

      ndropcol(icol) = 0._r8
      do kk = top_lev, pver
         ndropmix(icol,kk) = (qcld(kk) - ncldwtr(icol,kk))*dtinv - nsource(icol,kk)
         tendnd(icol,kk)   = (max(qcld(kk), 1.e-6_r8) - ncldwtr(icol,kk))*dtinv
         ndropcol(icol)   = ndropcol(icol) + ncldwtr(icol,kk)*pdel(icol,kk)
      end do
      ndropcol(icol) = ndropcol(icol)/gravit

!  BJG:  below probably always true?
      if (prog_modal_aero) then

         raertend = 0._r8
         qqcwtend = 0._r8

         do imode = 1, ntot_amode
            do lspec = 0, nspec_amode(imode)

               mm   = mam_idx(imode,lspec)
               lptr = mam_cnst_idx(imode,lspec)

               raertend(top_lev:pver) = (raercol(top_lev:pver,mm,nnew) - raer(mm)%fld(icol,top_lev:pver))*dtinv
               qqcwtend(top_lev:pver) = (raercol_cw(top_lev:pver,mm,nnew) - qqcw(mm)%fld(icol,top_lev:pver))*dtinv

               coltend(icol,mm)    = sum( pdel(icol,:)*raertend )/gravit
               coltend_cw(icol,mm) = sum( pdel(icol,:)*qqcwtend )/gravit

               ptend%q(icol,:,lptr) = 0.0_r8
               ptend%q(icol,top_lev:pver,lptr) = raertend(top_lev:pver)           ! set tendencies for interstitial aerosol
               qqcw(mm)%fld(icol,:) = 0.0_r8
               qqcw(mm)%fld(icol,top_lev:pver) = max(raercol_cw(top_lev:pver,mm,nnew),0.0_r8) ! update cloud-borne aerosol; HW: ensure non-negative
            end do
         end do

      end if

   end do  ! overall_main_i_loop
   ! end of main loop over i/longitude ....................................

   call outfld('NDROPCOL', ndropcol, pcols, lchnk)
   call outfld('NDROPSRC', nsource,  pcols, lchnk)
   call outfld('NDROPMIX', ndropmix, pcols, lchnk)
   call outfld('WTKE    ', wtke,     pcols, lchnk)

  !Note for C++ port: Get the cloud borne MMRs from AD in variable qcldbrn, do not port the code below

   ! Extract cloud borne MMRs from pbuf 

   qcldbrn(:,:,:,:) = huge(qcldbrn) !store invalid values
   do imode=1,ntot_amode
      do kk=top_lev,pver
         do lspec =1, nspec_amode(imode)
           spc_idx = lmassptrcw_amode(lspec,imode)
           fldcw => qqcw_get_field(pbuf, spc_idx, lchnk,.true.)
           qcldbrn(:,lspec,kk,imode) = fldcw(:,kk)
         enddo
         num_idx = numptrcw_amode(imode)
         fldcw => qqcw_get_field(pbuf,num_idx,lchnk,.true.)
         qcldbrn_num(:,kk,imode) = fldcw(:,kk)      
      enddo
   enddo

  !End note for C++ port

   call ccncalc(state_q, temp, qcldbrn, qcldbrn_num, ncol, cs, ccn)
   do lsat = 1, psat
      call outfld(ccn_name(lsat), ccn(1,1,lsat), pcols, lchnk)
   enddo


!  BJG:  below probably always false?
   if(do_aerocom_ind3) then
      ccn3d(:ncol, :) = ccn(:ncol, :, 4)
      ccn3col = 0.0_r8; ccn4col = 0.0_r8
      do icol=1, ncol
        do kk=1, pver
          ccn3col(icol) = ccn3col(icol) + ccn(icol,kk,3) * 1.0e6*   &
             pdel(icol,kk)/gravit/(pmid(icol,kk)/(temp(icol,kk)*rair))  !#/cm3 --> #/m2
          ccn4col(icol) = ccn4col(icol) + ccn(icol,kk,4) * 1.0e6*   &
             pdel(icol,kk)/gravit/(pmid(icol,kk)/(temp(icol,kk)*rair))  !#/cm3 --> #/m2
        enddo

! calculate CCN at 1km
        zi2 = 0.0
        zm2 = 0.0
        zmflag = .true.
        do kk=pver, 1, -1
          zi2(kk) = zi2(kk+1) + pdel(icol,kk)/gravit/(pmid(icol,kk)/(temp(icol,kk)*rair)) !
          zm2(kk) = (zi2(kk+1)+zi2(kk))/2._r8
          if(zm2(kk).gt.1000. .and. zmflag) then
            idx1000 = min(kk, pver-1)
            zmflag = .false.
          end if
        end do
        ccn3bl(icol) = (ccn(icol,idx1000,3)*(1000.-zm2(idx1000+1))+ccn(icol,idx1000+1,3) * (zm2(idx1000)-1000.)) &
                     /(zm2(idx1000)-zm2(idx1000+1)) * 1.0e6  ! #/cm3 -->#/m3
        ccn4bl(icol) = (ccn(icol,idx1000,4)*(1000.-zm2(idx1000+1))+ccn(icol,idx1000+1,4) * (zm2(idx1000)-1000.)) &
                     /(zm2(idx1000)-zm2(idx1000+1)) *1.0e6   ! #/cm3 -->#/m3
      enddo
      call outfld('colccn.1', ccn3col, pcols, lchnk)
      call outfld('colccn.3', ccn4col, pcols, lchnk)
      call outfld('ccn.1bl', ccn3bl, pcols, lchnk)
      call outfld('ccn.3bl', ccn4bl, pcols, lchnk)
   end if

!  BJG:  below probably always true?
   ! do column tendencies
   if (prog_modal_aero) then
      do imode = 1, ntot_amode
         do lspec = 0, nspec_amode(imode)
            mm = mam_idx(imode,lspec)
            call outfld(fieldname(mm),    coltend(:,mm),    pcols, lchnk)
            call outfld(fieldname_cw(mm), coltend_cw(:,mm), pcols, lchnk)
         end do
      end do
   end if

   deallocate( &
      nact,       &
      mact,       &
      raer,       &
      qqcw,       &
      raercol,    &
      raercol_cw, &
      coltend,    &
      coltend_cw, &
      naermod,    &
      hygro,      &
      vaerosol,   &
      fn,         &
      fm,         &
      fluxn,      &
      fluxm       )


end subroutine dropmixnuc

!===============================================================================

subroutine explmix( q, src, ekkp, ekkm, overlapp, overlapm, &
   qold, pver, dt, is_unact, qactold )

   !  explicit integration of droplet/aerosol mixing
   !     with source due to activation/nucleation


   integer, intent(in) :: pver ! number of levels
   real(r8), intent(out) :: q(pver) ! number / mass mixing ratio to be updated [# or kg / kg]
   real(r8), intent(in) :: qold(pver) ! number / mass mixing ratio from previous time step [# or kg / kg]
   real(r8), intent(in) :: src(pver) ! source due to activation/nucleation [# or kg / (kg-s)]
   real(r8), intent(in) :: ekkp(pver) ! zn*zs*density*diffusivity (kg/m3 m2/s) at interface [/s]
   ! below layer k  (k,k+1 interface)
   real(r8), intent(in) :: ekkm(pver) ! zn*zs*density*diffusivity (kg/m3 m2/s) at interface [/s]
   ! above layer k  (k,k+1 interface)
   real(r8), intent(in) :: overlapp(pver) ! cloud overlap below [fraction]
   real(r8), intent(in) :: overlapm(pver) ! cloud overlap above [fraction]
   real(r8), intent(in) :: dt ! time step [s]
   logical, intent(in) :: is_unact ! true if this is an unactivated species
   real(r8), intent(in),optional :: qactold(pver) ! [# or kg / kg]
   ! number / mass mixing ratio of ACTIVATED species from previous step
   ! *** this should only be present
   !     if the current species is unactivated number/sfc/mass

   integer k,kp1,km1

      !     the qactold*(1-overlap) terms are resuspension of activated material
      do k=top_lev,pver
         kp1=min(k+1,pver)
         km1=max(k-1,top_lev)

         if ( is_unact ) then
            q(k) = qold(k) + dt*( - src(k) + ekkp(k)*(qold(kp1) - qold(k) +       &
               qactold(kp1)*(1.0_r8-overlapp(k)))               &
               + ekkm(k)*(qold(km1) - qold(k) +     &
               qactold(km1)*(1.0_r8-overlapm(k))) )
         else
            q(k) = qold(k) + dt*(src(k) + ekkp(k)*(overlapp(k)*qold(kp1)-qold(k)) +      &
               ekkm(k)*(overlapm(k)*qold(km1)-qold(k)) )
         endif

         !        force to non-negative
         q(k)=max(q(k),0._r8)
      enddo

end subroutine explmix

!===============================================================================

subroutine activate_modal(wbar, wmaxf, tair, rhoair,  &
     na, nmode, volume, hygro, &
     fn, fm, fluxn, fluxm, flux_fullact, smax_prescribed )

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
  real(r8), intent(in) :: wbar      ! grid cell mean vertical velocity [m/s]
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
  real(r8), parameter :: p0 = 1013.25e2_r8    ! reference pressure [Pa]
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

  real(r8) wnuc  ! nucleation w, but = wbar if wdiab == 0 [m/s]
  real(r8) alw ! [/s]
  real(r8) smax ! maximum supersaturation [fraction]
  real(r8) lnsmax ! ln(smax) [unitless]
  real(r8) arg_erf_n,arg_erf_m  ! [unitless]
  real(r8) etafactor1  ! [/ s^(3/2)]
  real(r8) etafactor2(nmode),etafactor2max ! [s^(3/2)]
  integer imode

  !initialize activated aerosols and their fluxes to zero for all the modes
  fn(:)=0._r8
  fm(:)=0._r8
  fluxn(:)=0._r8
  fluxm(:)=0._r8
  flux_fullact=0._r8

  !return if aerosol number is negligible in the accumulation mode
  if(na(1) < 1.e-20_r8)return

  !return if mean vertical velocity is 0 or negative
  if(wbar <= 0._r8)return

  if ( present( smax_prescribed ) ) then
     !return if max supersaturation is 0 or negative
     if (smax_prescribed <= 0.0_r8) return
  end if

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
  wnuc = wbar
  alw=alpha*wnuc
  etafactor1=alw*sqrt(alw)
  zeta=twothird*sqrt(alw)*aten/sqrt(gthermfac)

  !Here compute smc, eta for all modes for maxsat calculation
  do imode=1,nmode
     if(volume(imode) > 1.e-39_r8 .and. na(imode) > 1.e-39_r8)then
        !number mode radius (m)
        amcube(imode)=(3._r8*volume(imode)/(4._r8*pi*exp45logsig(imode)*na(imode)))  ! only if variable size dist
        !Growth coefficent Abdul-Razzak & Ghan 1998 eqn 16
        !should depend on mean radius of mode to account for gas kinetic effects
        !see Fountoukis and Nenes, JGR2005 and Meskhidze et al., JGR2006
        !for approriate size to use for effective diffusivity.

        etafactor2(imode)=1._r8/(na(imode)*beta*sqrt(gthermfac))
        if(hygro(imode) > 1.e-10_r8)then
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
     fluxn(imode)=fn(imode)*wbar !activated aerosol number flux
     fluxm(imode)=fm(imode)*wbar !activated aerosol mass flux
  enddo
  flux_fullact = wbar

end subroutine activate_modal

!===============================================================================

subroutine maxsat(zeta,eta,nmode,smc,smax)

   !      calculates maximum supersaturation for multiple
   !      competing aerosol modes.

   !      Abdul-Razzak and Ghan, A parameterization of aerosol activation.
   !      2. Multiple aerosol types. J. Geophys. Res., 105, 6837-6844.

   integer,  intent(in)  :: nmode ! number of modes
   real(r8), intent(in)  :: smc(nmode) ! critical supersaturation for number mode radius [fraction]
   real(r8), intent(in)  :: zeta ! [dimensionless]
   real(r8), intent(in)  :: eta(nmode) ! [dimensionless]
   real(r8), intent(out) :: smax ! maximum supersaturation [fraction]
   integer  :: m  ! mode index
   real(r8) :: sum, g1, g2
   logical  :: weak_forcing  !  whether forcing is sufficiently weak or not


   weak_forcing = .true.
   do m=1, nmode
      if(zeta > 1.e5_r8*eta(m).or.smc(m)*smc(m) > 1.e5_r8*eta(m))then
         !weak forcing. essentially none activated
         smax=1.e-20_r8
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
      if(eta(m) > 1.e-20_r8)then
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
   real(r8), pointer, intent(in)  :: state_q(:,:,:) ! aerosol mmrs [kg/kg]     
   real(r8), pointer, intent(in)  :: tair(:,:)     ! air temperature [K]
   real(r8), intent(in)  :: qcldbrn(:,:,:,:), qcldbrn_num(:,:,:) ! cloud-borne aerosol mass / number  mixing ratios [kg/kg or #/kg]
   integer, intent(in)   :: ncol  ! number of columns
   real(r8), intent(in)  :: cs(pcols,pver)       ! air density [kg/m3]

  ! output arguments
   real(r8), intent(out) :: ccn(pcols,pver,psat) ! number conc of aerosols activated at supersat [#/m3]

   ! local


   real(r8) :: naerosol(pcols) ! interstit+activated aerosol number conc [#/m3]
   real(r8) :: vaerosol(pcols) ! interstit+activated aerosol volume conc [m3/m3]
   real(r8) :: amcube(pcols)  ! [m3]
   real(r8) :: amcubecoef(ntot_amode) ! [dimensionless]
   real(r8) :: argfactor(ntot_amode)  ! [dimensionless]
   real(r8) :: surften_coef  ! [m-K]
   real(r8) :: aparam(pcols) ! surface tension parameter  [m]
   real(r8) :: hygro(pcols)  ! aerosol hygroscopicity [dimensionless]
   real(r8) :: sm(pcols)  ! critical supersaturation at mode radius [fraction]
   real(r8) :: arg_erf_ccn ! [dimensionless] 
   real(r8) :: smcoef(pcols)  ! [m^(3/2)]
   integer lsat,imode,icol,kk
   integer phase ! phase of aerosol

   !     mathematical constants
   real(r8) percent_to_fraction, per_m3_to_per_cm3
   real(r8) smcoefcoef
   real(r8) super(psat) ! supersaturation [fraction]
   !-------------------------------------------------------------------------------

   phase=3 ! interstitial+cloudborne

   percent_to_fraction = 0.01_r8
   per_m3_to_per_cm3 = 1.e-6_r8
   super(:)=supersat(:)*percent_to_fraction
   smcoefcoef=2._r8/sqrt(27._r8)

   surften_coef=2._r8*mwh2o*surften/(r_universal*rhoh2o)

   ccn(:,:,:) = 0._r8

   do imode=1,ntot_amode

      amcubecoef(imode)=3._r8/(4._r8*pi*exp45logsig(imode))
      argfactor(imode)=twothird/(sq2*alogsig(imode))

      do kk=top_lev,pver

         call loadaer( &
            state_q, 1, ncol, kk, &
            imode, nspec_amode(imode), cs, phase, naerosol, vaerosol, &
            hygro, &
            qcldbrn(:,:,kk,imode), qcldbrn_num(:,kk,imode) )  ! optional in

         aparam(1:ncol) = surften_coef/tair(1:ncol,kk)
         smcoef(1:ncol)=smcoefcoef*aparam(1:ncol)*sqrt(aparam(1:ncol))

         where(naerosol(:ncol)>1.e-3_r8)
            amcube(:ncol)=amcubecoef(imode)*vaerosol(:ncol)/naerosol(:ncol)
            sm(:ncol)=smcoef(:ncol)/sqrt(hygro(:ncol)*amcube(:ncol)) ! critical supersaturation
         elsewhere
            sm(:ncol)=1._r8 ! value shouldn't matter much since naerosol is small
         endwhere

         do lsat=1,psat
            do icol=1,ncol
               arg_erf_ccn=argfactor(imode)*log(sm(icol)/super(lsat))
               ccn(icol,kk,lsat)=ccn(icol,kk,lsat)+naerosol(icol)*0.5_r8*(1._r8-erf(arg_erf_ccn))
            enddo
         enddo

      enddo

   enddo

   ccn(:ncol,:,:)=ccn(:ncol,:,:)*per_m3_to_per_cm3 ! convert from #/m3 to #/cm3

end subroutine ccncalc

!===============================================================================

subroutine loadaer( &
     state_q, istart, istop, kk, &
     imode, nspec, cs, phase, naerosol, &
     vaerosol, hygro, &
     qcldbrn1d, qcldbrn1d_num) ! optional in

  ! return aerosol number, volume concentrations, and bulk hygroscopicity

  ! input arguments
  real(r8), intent(in) :: state_q(:,:,:)        ! aerosol mmrs [kg/kg]

  integer,  intent(in) :: istart      ! start column index (1 <= istart <= istop <= pcols)
  integer,  intent(in) :: istop       ! stop column index
  integer,  intent(in) :: imode       ! mode index
  integer,  intent(in) :: nspec       ! total # of species in mode imode
  integer,  intent(in) :: kk          ! level index
  real(r8), intent(in) :: cs(:,:)     ! air density [kg/m3]
  integer,  intent(in) :: phase       ! phase of aerosol: 1 for interstitial, 2 for cloud-borne, 3 for sum

  ! output arguments
  real(r8), intent(out) :: naerosol(:)  ! number conc [#/m3]
  real(r8), intent(out) :: vaerosol(:)  ! volume conc [m3/m3]
  real(r8), intent(out) :: hygro(:)     ! bulk hygroscopicity of mode [dimensionless]

  ! optional input arguments
  real(r8), intent(in), optional  :: qcldbrn1d(:,:), qcldbrn1d_num(:) ! ! cloud-borne aerosol mass / number  mixing ratios [kg/kg or #/kg]

  ! internal

  real(r8) :: vaerosolsum(pcols)  ! sum to find volume conc [m3/kg]
  real(r8) :: hygrosum(pcols)     ! sum to bulk hygroscopicity of mode [m3/kg]
  real(r8) :: qcldbrn_local(pcols,nspec)  ! local cloud-borne aerosol mass mixing ratios [kg/kg]
  real(r8) :: qcldbrn_num_local(pcols) ! local cloud-borne aerosol number mixing ratios [#/kg]

  integer  :: icol, lspec, spc_idx
  !-------------------------------------------------------------------------------

  !Currenly supports only phase 1 (interstitial) and 3 (interstitial+cldbrn)
  if (phase /= 1 .and. phase /=3) then
     write(iulog,*)'phase=',phase,' in loadaer'
     call endrun('phase error in loadaer')
  endif

  qcldbrn_local(:,:) = 0._r8
  if(present(qcldbrn1d)) qcldbrn_local(:,:nspec) = qcldbrn1d(:,:nspec)

  vaerosolsum(:) = 0._r8
  hygrosum(:)    = 0._r8

  !Sum over all species within imode to get bulk hygroscopicity and volume conc
     !phase == 1 is interstitial only.
     !phase == 3 is interstitial + cldborne
!  Assumes iphase =1 or 3, so interstitial is always summed, added with cldbrn when present
!  iphase = 2 would require alternate logic from following subroutine

  call get_aer_mmr_sum(imode, nspec, istart, istop, state_q(:,kk,:), qcldbrn_local(:,:nspec), & !in
          vaerosolsum, hygrosum)    !inout

  !  Finalize computation of bulk hygrospopicity and volume conc
  do icol = istart, istop
     if (vaerosolsum(icol) > 1.0e-30_r8) then
        hygro(icol)    = hygrosum(icol)/(vaerosolsum(icol))
        vaerosol(icol) = vaerosolsum(icol)*cs(icol,kk)
     else
        hygro(icol)    = 0.0_r8
        vaerosol(icol) = 0.0_r8
     endif
  enddo

  qcldbrn_num_local(:) = 0._r8
  if(present(qcldbrn1d_num)) qcldbrn_num_local(:) = qcldbrn1d_num(:)

  ! Compute aerosol number concentration
  call get_aer_num(imode, istart, istop, state_q(:,kk,:), cs(:,kk), vaerosol, qcldbrn_num_local, &!in
          naerosol) !out

end subroutine loadaer

!===============================================================================

subroutine get_aer_mmr_sum(imode, nspec, istart, istop, state_q, qcldbrn1d, & !in
     vaerosolsum, hygrosum)   !inout

  !add these for direct access to mmr (in state_q array), density and hygroscopicity
  use modal_aero_data,   only: lspectype_amode, specdens_amode, spechygro, lmassptr_amode

  ! input arguments
  integer,  intent(in) :: imode       ! mode index
  integer,  intent(in) :: nspec       ! total # of species in mode imode
  integer,  intent(in) :: istart      ! start column index (1 <= istart <= istop <= pcols)
  integer,  intent(in) :: istop       ! stop column index
  real(r8), intent(in) :: state_q(:,:) ! interstitial aerosol mass mixing ratios [kg/kg]
  real(r8), intent(in) :: qcldbrn1d(:,:) ! cloud-borne aerosol mass mixing ratios [kg/kg]

  ! in/out arguments
  real(r8), intent(inout) :: vaerosolsum(:)  ! sum to find volume conc [m3/kg]
  real(r8), intent(inout) :: hygrosum(:)   ! sum to bulk hygroscopicity of mode [m3/kg]

  ! internal
  real(r8) :: density_sp ! density at species / mode indices [kg/m3]
  real(r8) :: hygro_sp   ! hygroscopicity at species / mode indices [dimensionless]
  real(r8) :: vol     !aerosol volume mixing ratio [m3/kg]

  integer  :: icol, lspec, spc_idx, type_idx

  !Start to compute bulk volume conc / hygroscopicity by summing over species per mode.
   do lspec = 1, nspec
      type_idx = lspectype_amode(lspec,imode)
      density_sp  = specdens_amode(type_idx) !species density
      hygro_sp    = spechygro(type_idx)      !species hygroscopicity
      spc_idx   = lmassptr_amode(lspec,imode) !index of species in state_q array
      do icol = istart, istop
         vol = max(state_q(icol,spc_idx) + qcldbrn1d(icol,lspec), 0._r8)/density_sp !volume = mmr/density
         vaerosolsum(icol) = vaerosolsum(icol) + vol        !bulk volume
         hygrosum(icol)    = hygrosum(icol) + vol*hygro_sp !bulk hygroscopicity
      enddo
   enddo

end subroutine get_aer_mmr_sum

!===============================================================================

subroutine get_aer_num(imode, istart, istop, state_q, cs, vaerosol, qcldbrn1d_num, &!in
           naerosol) !out

  use modal_aero_data, only:numptr_amode

  ! input arguments
  integer,  intent(in) :: imode        ! mode index
  integer,  intent(in) :: istart       ! start column index (1 <= istart <= istop <= pcols)
  integer,  intent(in) :: istop        ! stop column index
  real(r8), intent(in) :: state_q(:,:) ! interstitial aerosol number mixing ratios [#/kg]
  real(r8), intent(in) :: cs(:)        ! air density [kg/m3]
  real(r8), intent(in) :: vaerosol(:)  ! volume conc [m3/m3]
  real(r8), intent(in) :: qcldbrn1d_num(:) ! cloud-borne aerosol number mixing ratios [#/kg]

  !output arguments
  real(r8), intent(out) :: naerosol(:)  ! number conc [#/m3]

  !internal
  integer  :: icol, num_idx

  !convert number mixing ratios to number concentrations
  !Use bulk volume conc found previously to bound value

  num_idx = numptr_amode(imode)
  do icol = istart, istop
     naerosol(icol) = (state_q(icol,num_idx) + qcldbrn1d_num(icol))*cs(icol)
     !adjust number so that dgnumlo < dgnum < dgnumhi
     naerosol(icol) = max(naerosol(icol), vaerosol(icol)*voltonumbhi_amode(imode))
     naerosol(icol) = min(naerosol(icol), vaerosol(icol)*voltonumblo_amode(imode))
  enddo
end subroutine get_aer_num

!===============================================================================

end module ndrop
