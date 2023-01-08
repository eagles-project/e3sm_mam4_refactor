module wetdep

!----------------------------------------------------------------------- 
!
! Wet deposition routines for both aerosols and gas phase constituents.
! 
!-----------------------------------------------------------------------

use shr_kind_mod, only: r8 => shr_kind_r8
use ppgrid,       only: pcols, pver
use physconst,    only: gravit, rair, tmelt
use cam_abortutils, only: endrun
use mam_support,  only: min_max_bound

implicit none
save
private

public :: wetdepa_v2  ! scavenging codes for very soluble aerosols -- CAM5 version
public :: clddiag     ! calc of cloudy volume and rain mixing ratio
public :: faer_resusp_vs_fprec_evap_mpln

public :: wetdep_inputs_t
public :: wetdep_init
public :: wetdep_inputs_set
public :: wetdep_inputs_unset

type wetdep_inputs_t
   real(r8), pointer :: cldt(:,:) => null()  ! cloud fraction
   real(r8), pointer :: qme(:,:) => null()
   real(r8), pointer :: prain(:,:) => null()
   real(r8), pointer :: evapr(:,:) => null()
   real(r8), allocatable :: cldcu(:,:)  ! convective cloud fraction, currently empty
   real(r8), allocatable :: evapc(:,:)  ! Evaporation rate of convective precipitation
   real(r8), allocatable :: cmfdqr(:,:) ! convective production of rain
   real(r8), allocatable :: conicw(:,:) ! convective in-cloud water
   real(r8), allocatable :: totcond(:,:)! total condensate
   real(r8), allocatable :: cldv(:,:)   ! cloudy volume undergoing wet chem and scavenging
   real(r8), allocatable :: cldvcu(:,:) ! Convective precipitation area at the top interface of current layer
   real(r8), allocatable :: cldvst(:,:) ! Stratiform precipitation area at the top interface of current layer 
end type wetdep_inputs_t

integer :: cld_idx             = 0
integer :: qme_idx             = 0 
integer :: prain_idx           = 0 
integer :: nevapr_idx          = 0 
integer :: icwmrdp_idx         = 0 
integer :: icwmrsh_idx         = 0 
integer :: rprddp_idx          = 0 
integer :: rprdsh_idx          = 0 
integer :: sh_frac_idx         = 0 
integer :: dp_frac_idx         = 0 
integer :: nevapr_shcu_idx     = 0 
integer :: nevapr_dpcu_idx     = 0 
integer :: ixcldice, ixcldliq

real(r8), parameter :: rhoh2o = 1000._r8            ! density of water

! declare several small values that to avoid divided by zero used in the module
real(r8), parameter :: small_value_30 = 1.e-30_r8   
real(r8), parameter :: small_value_12 = 1.e-12_r8 
real(r8), parameter :: small_value_14 = 1.e-14_r8
real(r8), parameter :: small_value_36 = 1.e-36_r8  
real(r8), parameter :: small_value_5  = 1.e-5_r8        ! for cloud fraction
real(r8), parameter :: small_value_2  = 1.e-2_r8        ! for cloud fraction

!==============================================================================
contains
!==============================================================================

!==============================================================================
subroutine wetdep_init()
  use physics_buffer, only: pbuf_get_index
  use constituents,   only: cnst_get_ind
  use phys_control,   only: phys_getopts

  cld_idx             = pbuf_get_index('CLD')    
  qme_idx             = pbuf_get_index('QME')    
  prain_idx           = pbuf_get_index('PRAIN')  
  nevapr_idx          = pbuf_get_index('NEVAPR') 

  icwmrdp_idx         = pbuf_get_index('ICWMRDP') 
  rprddp_idx          = pbuf_get_index('RPRDDP')  
  icwmrsh_idx         = pbuf_get_index('ICWMRSH') 
  rprdsh_idx          = pbuf_get_index('RPRDSH')  
  sh_frac_idx         = pbuf_get_index('SH_FRAC' )
  dp_frac_idx         = pbuf_get_index('DP_FRAC') 
  nevapr_shcu_idx     = pbuf_get_index('NEVAPR_SHCU') 
  nevapr_dpcu_idx     = pbuf_get_index('NEVAPR_DPCU') 

  call cnst_get_ind('CLDICE', ixcldice)
  call cnst_get_ind('CLDLIQ', ixcldliq)

endsubroutine wetdep_init

!==============================================================================
subroutine wetdep_inputs_set( state, pbuf, inputs )
! -----------------------------------------------------------------------------
! gathers up the inputs needed for the wetdepa routines
! -----------------------------------------------------------------------------


  use physics_types,  only: physics_state
  use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx

  ! args

  type(physics_state),  intent(in )  :: state           !! physics state
  type(physics_buffer_desc), pointer :: pbuf(:)         !! physics buffer
  type(wetdep_inputs_t), intent(out) :: inputs          !! collection of wetdepa inputs

  ! local vars
  ! pointers to pbuf, will assign to variables later
  real(r8), pointer :: p_icwmrdp(:,:)    ! in cloud water mixing ratio, deep convection
  real(r8), pointer :: p_rprddp(:,:)     ! rain production, deep convection
  real(r8), pointer :: p_icwmrsh(:,:)    ! in cloud water mixing ratio, shallow convection
  real(r8), pointer :: p_rprdsh(:,:)     ! rain production, shallow convection
  real(r8), pointer :: p_sh_frac(:,:)    ! Shallow convective cloud fraction
  real(r8), pointer :: p_dp_frac(:,:)    ! Deep convective cloud fraction
  real(r8), pointer :: p_evapcsh(:,:)    ! Evaporation rate of shallow convective precipitation >=0.
  real(r8), pointer :: p_evapcdp(:,:)    ! Evaporation rate of deep    convective precipitation >=0.

  ! assign variables, some of them are in "state" and "inputs"
  ! in state
  real(r8) :: temperature(pcols,pver)        ! temperature [K]
  real(r8) :: q_liq(pcols,pver)    ! liquid water mixing ratio[kg/kg]
  real(r8) :: q_ice(pcols,pver)    ! ice water mixing ratio [kg/kg]
  real(r8) :: pmid(pcols,pver)     ! pressure at layer midpoints[Pa]
  real(r8) :: pdel(pcols,pver)     ! pressure difference across layers [Pa]
  ! data from pbuf pointers
  real(r8) :: icwmrdp(pcols,pver)     ! in cloud water mixing ratio, deep convection [kg/kg]
  real(r8) :: rprddp(pcols,pver)      ! rain production, deep convection [kg/kg/s]
  real(r8) :: icwmrsh(pcols,pver)     ! in cloud water mixing ratio, shallow convection [kg/kg]
  real(r8) :: rprdsh(pcols,pver)      ! rain production, shallow convection [kg/kg/s]
  real(r8) :: sh_frac(pcols,pver)     ! Shallow convective cloud fraction [fraction]
  real(r8) :: dp_frac(pcols,pver)     ! Deep convective cloud fraction [fraction]
  real(r8) :: evapcsh(pcols,pver)     ! Evaporation rate of shallow convective precipitation >=0. [kg/kg/s]
  real(r8) :: evapcdp(pcols,pver)     ! Evaporation rate of deep    convective precipitation >=0.[ kg/kg/s]
  ! in inputs
  real(r8) :: cldcu(pcols,pver)       ! cumulus cloud fraction [fraction]
  real(r8) :: cldt(pcols,pver)        ! total cloud fraction [fraction]
  real(r8) :: evapc(pcols,pver)       ! evaporation from convection (deep + shallow) [kg/kg/s]
  real(r8) :: evapr(pcols,pver)       ! evaporation from stratiform rain [kg/kg/s]
  real(r8) :: cmfdqr(pcols,pver)      ! dq/dt due to convective rainout [kg/kg/s]
  real(r8) :: qme(pcols,pver)         ! net condensation/evaporation of cloud water [kg/kg/s]
  real(r8) :: conicw(pcols,pver)      ! convective in-cloud water [kg/kg]
  real(r8) :: totcond(pcols,pver)     ! total condensate (ice+liq) [kg/kg]
  real(r8) :: prain(pcols,pver)       ! stratiform rain production rate [kg/kg/s]
  real(r8) :: cldv(pcols,pver)        ! total (rain+cloud) cloudy volume undergoing wet chem and scavenging [fraction]
  real(r8) :: cldvcu(pcols,pver)      ! Convective precipitation and cloudy volume,at the top interface of current layer [fraction]
  real(r8) :: cldvst(pcols,pver)      ! stratiform precipitation and cloudy volume,at the top interface of current layer [fraction]
  ! locally used
  real(r8) :: rainmr(pcols,pver)       ! mixing ratio of rain within cloud volume [kg/kg]
  real(r8) :: cldst(pcols,pver)        ! Stratiform cloud fraction [fraction]

  integer :: itim, ncol
  integer :: ierror

  ! **************************************************
  ! this part only extract/define variables from pbuf, state and inputs
  ! get variables from pbuf
  itim = pbuf_old_tim_idx()

  allocate (inputs%cldcu(pcols,pver), stat=ierror)
  if ( ierror /= 0 ) call endrun('WETDEP_INPUTS_SET error: allocation error cldcu')

  allocate (inputs%evapc(pcols,pver), stat=ierror)
  if ( ierror /= 0 ) call endrun('WETDEP_INPUTS_SET error: allocation error evapc')

  allocate (inputs%cmfdqr(pcols,pver), stat=ierror)
  if ( ierror /= 0 ) call endrun('WETDEP_INPUTS_SET error: allocation error cmfdqr')

  allocate (inputs%conicw(pcols,pver), stat=ierror)
  if ( ierror /= 0 ) call endrun('WETDEP_INPUTS_SET error: allocation error conicw')

  allocate (inputs%totcond(pcols,pver), stat=ierror)
  if ( ierror /= 0 ) call endrun('WETDEP_INPUTS_SET error: allocation error totcond')

  allocate (inputs%cldv(pcols,pver), stat=ierror)
  if ( ierror /= 0 ) call endrun('WETDEP_INPUTS_SET error: allocation error cldv')

  allocate (inputs%cldvcu(pcols,pver), stat=ierror)
  if ( ierror /= 0 ) call endrun('WETDEP_INPUTS_SET error: allocation error cldvcu')

  allocate (inputs%cldvst(pcols,pver), stat=ierror)
  if ( ierror /= 0 ) call endrun('WETDEP_INPUTS_SET error: allocation error cldvst')

  call pbuf_get_field(pbuf, cld_idx,         inputs%cldt, start=(/1,1,itim/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, qme_idx,         inputs%qme     )
  call pbuf_get_field(pbuf, prain_idx,       inputs%prain   )
  call pbuf_get_field(pbuf, nevapr_idx,      inputs%evapr   )
  call pbuf_get_field(pbuf, icwmrdp_idx,     p_icwmrdp )
  call pbuf_get_field(pbuf, icwmrsh_idx,     p_icwmrsh )
  call pbuf_get_field(pbuf, rprddp_idx,      p_rprddp  )
  call pbuf_get_field(pbuf, rprdsh_idx,      p_rprdsh  )
  call pbuf_get_field(pbuf, sh_frac_idx,     p_sh_frac )
  call pbuf_get_field(pbuf, dp_frac_idx,     p_dp_frac )
  call pbuf_get_field(pbuf, nevapr_shcu_idx, p_evapcsh )
  call pbuf_get_field(pbuf, nevapr_dpcu_idx, p_evapcdp )

  ! assign all variables from state and inputs
  ncol = state%ncol
  temperature = state%t
  q_liq = state%q(:,:,ixcldliq)
  q_ice = state%q(:,:,ixcldice)
  pmid = state%pmid
  pdel = state%pdel

  icwmrdp = p_icwmrdp
  icwmrsh = p_icwmrsh
  rprddp = p_rprddp
  rprdsh = p_rprdsh
  sh_frac = p_sh_frac
  dp_frac = p_dp_frac
  evapcsh = p_evapcsh
  evapcdp = p_evapcdp
  
  qme = inputs%qme     ! not used in this subroutine but will be used in wetdepa_v2
  cldt = inputs%cldt
  evapr = inputs%evapr
  prain = inputs%prain
  ! ************************************************** end

  ! ****************************************************
  ! this is the main calculation in this subroutine
  ! calculate some variables needed in wetdepa_v2
  cldcu(:ncol,:)  = dp_frac(:ncol,:) + sh_frac(:ncol,:)          !cumulus cloud fraction
  cldst(:ncol,:)  = cldt(:ncol,:) - cldcu(:ncol,:)               ! Stratiform cloud fraction
  evapc(:ncol,:)  = evapcsh(:ncol,:) + evapcdp(:ncol,:)          ! evaporation from convection (deep + shallow)
  cmfdqr(:ncol,:) = rprddp(:ncol,:)  + rprdsh(:ncol,:)           ! dq/dt due to convective rainout

  ! sum deep and shallow convection contributions
  conicw(:ncol,:) = (icwmrdp(:ncol,:)*dp_frac(:ncol,:) + icwmrsh(:ncol,:)*sh_frac(:ncol,:))/ &
                              max(small_value_2, sh_frac(:ncol,:) + dp_frac(:ncol,:))
  totcond(:ncol,:) = q_liq(:ncol,:) + q_ice(:ncol,:)

  ! calculate cldv, cldvcu and cldvst
  call clddiag( ncol, temperature, pmid,   pdel,   cmfdqr, evapc, & ! in
               cldt,  cldcu,       cldst,  evapr,  prain,         & ! in
               cldv,  cldvcu,      cldvst, rainmr                 ) ! out
  ! *********************************************** end

  ! **********************************************
  !  assign variables back in inputs
  inputs%cldcu = cldcu
  inputs%evapc = evapc
  inputs%cmfdqr = cmfdqr
  inputs%conicw = conicw
  inputs%totcond = totcond
  inputs%cldv = cldv
  inputs%cldvcu = cldvcu
  inputs%cldvst = cldvst
  inputs%qme = qme
  ! **********************************************end

end subroutine wetdep_inputs_set

!==============================================================================
subroutine wetdep_inputs_unset(inputs)
! -----------------------------------------------------------------------------
! deallocate storage assoicated with wetdep_inputs_t type variable
! -----------------------------------------------------------------------------

  ! args
  type(wetdep_inputs_t), intent(inout) :: inputs          !! collection of wetdepa inputs

  deallocate(inputs%cldcu)
  deallocate(inputs%evapc)
  deallocate(inputs%cmfdqr)
  deallocate(inputs%conicw)
  deallocate(inputs%totcond)
  deallocate(inputs%cldv)
  deallocate(inputs%cldvcu)
  deallocate(inputs%cldvst)

end subroutine wetdep_inputs_unset

!==============================================================================
subroutine clddiag(ncol, temperature, pmid, pdel, cmfdqr, evapc, & ! in
                   cldt, cldcu,       cldst,      evapr, prain,  & ! in
                   cldv, cldvcu,      cldvst,     rain           ) ! out
! ------------------------------------------------------------------------------------ 
! Estimate the cloudy volume which is occupied by rain or cloud water as
! the max between the local cloud amount or the
! sum above of (cloud*positive precip production)      sum total precip from above
!              ----------------------------------   x ------------------------
! sum above of     (positive precip           )        sum positive precip from above
! Author: P. Rasch
!         Sungsu Park. Mar.2010 
! ------------------------------------------------------------------------------------

   ! Input arguments:
   real(r8), intent(in) :: temperature(pcols,pver)        ! temperature [K]
   real(r8), intent(in) :: pmid(pcols,pver)     ! pressure at layer midpoints[Pa]
   real(r8), intent(in) :: pdel(pcols,pver)     ! pressure difference across layers [Pa]
   real(r8), intent(in) :: cmfdqr(pcols,pver)   ! dq/dt due to convective rainout [kg/kg/s]
   real(r8), intent(in) :: evapc(pcols,pver)    ! Evaporation rate of convective precipitation ( >= 0 ) [kg/kg/s]
   real(r8), intent(in) :: cldt(pcols,pver)     ! total cloud fraction [fraction, unitless]
   real(r8), intent(in) :: cldcu(pcols,pver)    ! Cumulus cloud fraction [fraction, unitless]
   real(r8), intent(in) :: cldst(pcols,pver)    ! Stratus cloud fraction [fraction, unitless]
   real(r8), intent(in) :: evapr(pcols,pver)    ! rate of evaporation of falling precipitation [kg/kg/s]
   real(r8), intent(in) :: prain(pcols,pver)    ! rate of conversion of condensate to precipitation [kg/kg/s]
   integer, intent(in) :: ncol

   ! Output arguments:
   real(r8), intent(out) :: cldv(pcols,pver)     ! fraction occupied by rain or cloud water [fraction, unitless]
   real(r8), intent(out) :: cldvcu(pcols,pver)   ! Convective precipitation volume [fraction, unitless]
   real(r8), intent(out) :: cldvst(pcols,pver)   ! Stratiform precipitation volume [fraction, unitless]
   real(r8), intent(out) :: rain(pcols,pver)     ! mixing ratio of rain [kg/kg]

   ! Local variables:
   real(r8) :: sumppr_all(pcols,pver)    ! precipitation rate in all vertical levels [kg/m2/s]
   real(r8) :: lprec(pcols,pver)         ! local production rate of precip [kg/m2/s]
   real(r8) :: sumppr_cu_all(pcols,pver) ! same as sumppr_all but for conv.precip. calculated but not used
   real(r8) :: lprec_cu(pcols,pver)      ! Local production rate of convective precip [kg/m2/s]
   real(r8) :: sumppr_st_all(pcols,pver) ! same as sumppr_all but for strat.precip. calculated but not used
   real(r8) :: lprec_st(pcols,pver)      ! Local production rate of stratiform precip [kg/m2/s]
   ! -----------------------------------------------------------------------

   !calculate local precipitation rate
   !FIXME: Possible bug: why there is no evapc in lprec calculation?
   call local_precip_production(ncol, pdel, prain+cmfdqr, evapr, lprec)
   call local_precip_production(ncol, pdel, cmfdqr, evapc, lprec_cu)
   call local_precip_production(ncol, pdel, prain, evapr, lprec_st)

   ! calculate cloud volume occupied by rain or cloud water
   call calculate_cloudy_volume(ncol, cldt, lprec, .true., cldv, sumppr_all) !total
   call calculate_cloudy_volume(ncol, cldcu, lprec_cu, .false., cldvcu, sumppr_cu_all) !convective
   call calculate_cloudy_volume(ncol, cldst, lprec_st, .false., cldvst, sumppr_st_all) !stratiform

   ! calculate rain mixing ratio
   call rain_mix_ratio(temperature, pmid, sumppr_all, ncol, rain)

end subroutine clddiag

!==============================================================================
subroutine local_precip_production(ncol, pdel, source_term, sink_term, lprec)
!----------------------------------------------------------------------------
! calculate local precipitation generation rate (kg/m2/s) from
! source (condensation) and sink (evaporation) terms
!----------------------------------------------------------------------------
    ! Input arguments:
    real(r8), intent(in) :: pdel(pcols,pver)        ! pressure difference across layers [Pa]
    real(r8), intent(in) :: source_term(pcols,pver) ! precipitation source term rate (condensation) [kg/kg/s]
    real(r8), intent(in) :: sink_term(pcols,pver)   ! precipitation sink term rate (evaporation) [kg/kg/s]
    integer,  intent(in) :: ncol

    ! Output arguments:
    real(r8), intent(out) :: lprec(pcols,pver)     ! local production rate of precip [kg/m2/s]

    ! Local variables:
    integer :: icol, kk

    !calculate local precipitation rate
    do icol=1,ncol
       do kk=1,pver
          lprec(icol,kk)  = (pdel(icol,kk)/gravit)*(source_term(icol,kk)-sink_term(icol,kk))
       enddo
    enddo

end subroutine local_precip_production

!==============================================================================
subroutine calculate_cloudy_volume(ncol, cld, lprec, is_tot_cld, cldv, sumppr_all)
! ------------------------------------------------------------------------------------
! Calculate cloudy volume which is occupied by rain or cloud water as
! the max between the local cloud amount or the
! sum above of (cloud*positive precip production)      sum total precip from above
!              ----------------------------------   x   ------------------------
! sum above of     (positive precip           )        sum positive precip from above
! ------------------------------------------------------------------------------------
   ! Input arguments:
   real(r8), intent(in) :: cld(pcols,pver)   ! cloud fraction [fraction, unitless]
   real(r8), intent(in) :: lprec(pcols,pver) ! local production rate of precip [kg/m2/s]
   integer, intent(in) :: ncol
   logical, intent(in) :: is_tot_cld

   ! Output arguments:
   real(r8), intent(out) :: cldv(pcols,pver)       ! fraction occupied by rain or cloud water [fraction, unitless]
   real(r8), intent(out) :: sumppr_all(pcols,pver) ! sum of precipitation rate above each layer, for calling rain_mix_ratio use [kg/m2/s]

   ! Local variables:
   integer :: icol,kk
   real(r8) :: sumppr(pcols)        ! precipitation rate [kg/m2/s]
   real(r8) :: sumpppr(pcols)       ! sum of positive precips from above
   real(r8) :: cldv1(pcols)         ! precip weighted cloud fraction from above [kg/m2/s]
   real(r8) :: lprecp               ! local production rate of precip [kg/m2/s] if positive


   ! initiate variables
   do icol=1,ncol
      sumppr(icol) = 0._r8
      cldv1(icol) = 0._r8
      sumpppr(icol) = small_value_36   ! not 0 because it will be divided
   end do

   do kk = 1,pver
      do icol = 1,ncol
         if (is_tot_cld) then
             cldv(icol,kk) = max( min(1._r8,cldv1(icol)/sumpppr(icol)) *sumppr(icol)/sumpppr(icol), cld(icol,kk))
         else
             ! For convective and stratiform precipitation volume at the top interface of each layer. 
             ! Neglect the current layer.
             cldv(icol,kk) = max( min(1._r8,cldv1(icol)/sumpppr(icol)) * (sumppr(icol)/sumpppr(icol)), 0._r8)
         endif
         lprecp = max(lprec(icol,kk), small_value_30)
         cldv1(icol) = cldv1(icol)  + cld(icol,kk)*lprecp
         sumppr(icol) = sumppr(icol) + lprec(icol,kk)
         sumppr_all(icol,kk) = sumppr(icol)      ! save all sumppr to callrain_mix_ratio
         sumpppr(icol) = sumpppr(icol) + lprecp
      enddo
   enddo

end subroutine calculate_cloudy_volume

!==============================================================================
subroutine rain_mix_ratio(temperature, pmid, sumppr, ncol, rain)
!-----------------------------------------------------------------------
! Purpose:
! calculate rain mixing ratio from precipitation rate above.
!
! extracted from clddiag subroutine
! for C++ portint, Shuaiqi Tang in 9/22/2022
!-----------------------------------------------------------------------

   ! Input arguments:
   real(r8), intent(in) :: temperature(pcols,pver)      ! temperature [K]
   real(r8), intent(in) :: pmid(pcols,pver)   ! pressure at layer midpoints [Pa]
   real(r8), intent(in) :: sumppr(pcols,pver) ! sum of precipitation rate above each layer [kg/m2/s]
   integer,  intent(in) :: ncol

   ! Output arguments:
   real(r8), intent(out) :: rain(pcols,pver) ! mixing ratio of rain [kg/kg]

   ! Local variables:
   integer  icol,kk

   ! constants used in fallspeed calculation; taken from findmcnew
   real(r8) :: convfw
   real(r8) :: rho     !air density
   real(r8) :: vfall

   ! define the constant convfw. taken from findmcnew, do not find the reference
   ! of the equation  -- by Shuaiqi when refactoring for C++
   convfw = 1.94_r8*2.13_r8*sqrt(rhoh2o*gravit*2.7e-4_r8)

   do kk = 1,pver
     do icol = 1,ncol
         rain(icol,kk) = 0._r8
         if(temperature(icol,kk) .gt. tmelt) then
            rho = pmid(icol,kk)/(rair*temperature(icol,kk))
            vfall = convfw/sqrt(rho)
            rain(icol,kk) = sumppr(icol,kk) / (rho*vfall)
            if (rain(icol,kk) .lt. small_value_14) rain(icol,kk) = 0._r8
         endif
      end do
   end do

end subroutine rain_mix_ratio

!==============================================================================
subroutine wetdepa_v2( ncol, deltat,  pdel, &
                       cmfdqr, evapc, dlf, conicw, &
                       precs, evaps, cwat, &
                       cldt, cldc, cldvcu, cldvst, &
                       sol_factb, sol_facti, sol_factic, &
                       mam_prevap_resusp_optcc, is_strat_cloudborne, scavcoef, f_act_conv, &
                       tracer, qqcw, &
                       fracis, scavt, iscavt, &
                       icscavt, isscavt, bcscavt, bsscavt, rcscavt, rsscavt )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! scavenging code for very soluble aerosols
      ! 
      ! Author: P. Rasch
      ! Modified by T. Bond 3/2003 to track different removals
      ! Sungsu Park. Mar.2010 : Impose consistencies with a few changes in physics.

      ! this section of code is for highly soluble aerosols,
      ! the assumption is that within the cloud that
      ! all the tracer is in the cloud water
      !
      ! for both convective and stratiform clouds,
      ! the fraction of cloud water converted to precip defines
      ! the amount of tracer which is pulled out.
      !-----------------------------------------------------------------------

      implicit none

      integer, intent(in) :: ncol

      real(r8), intent(in) ::&
         deltat,               &! time step [s]
         pdel(pcols,pver),     &! pressure thikness [Pa]
         cmfdqr(pcols,pver),   &! rate of production of convective precip [kg/kg/s]
         evapc(pcols,pver),    &! Evaporation rate of convective precipitation [kg/kg/s]
         dlf(pcols,pver),      &! Detrainment of convective condensate [kg/kg/s]
         conicw(pcols,pver),   &! convective cloud water [kg/kg]
         precs(pcols,pver),    &! rate of production of stratiform precip [kg/kg/s]
         evaps(pcols,pver),    &! rate of evaporation of precip [kg/kg/s]
         cwat(pcols,pver),     &! cloud water amount [kg/kg] 
         cldt(pcols,pver),     &! total cloud fraction [fraction]
         cldc(pcols,pver),     &! convective cloud fraction [fraction]
         cldvcu(pcols,pver),   &! Convective precipitation area at the top interface of each layer [fraction]
         cldvst(pcols,pver),   &! Stratiform precipitation area at the top interface of each layer [fraction]
         tracer(pcols,pver)     ! trace species [kg/kg]

      integer, intent(in) :: mam_prevap_resusp_optcc ! suspension options.
!     0 = no resuspension
!     1 = linear resuspension of aerosol mass or number following original mam
!     coding and history_aero_prevap_resusp = .false.
!     2 = same as 1 but history_aero_prevap_resusp = .true.
!     3 = same as 2 but with some added "xxx = max( 0, xxx)" lines
!   130 = non-linear resuspension of aerosol mass   based on scavenged aerosol mass
!   230 = non-linear resuspension of aerosol number based on raindrop number
!     (1,2,3 are not used in the current code)

      ! is_strat_cloudborne = .true. if tracer is stratiform-cloudborne aerosol; else .false. 
      logical,  intent(in) :: is_strat_cloudborne   
      real(r8), intent(in) :: scavcoef(pcols,pver) ! Dana and Hales coefficient [1/mm]
      ! f_act_conv = conv-cloud activation fraction when is_strat_cloudborne==.false.; else 0.0 
      real(r8), intent(in) :: f_act_conv(pcols,pver) ! [fraction]
      ! qqcw = strat-cloudborne aerosol corresponding to tracer when is_strat_cloudborne==.false.; else 0.0
      real(r8), intent(in) :: qqcw(pcols,pver)  ! [kg/kg]
      real(r8), intent(in) :: sol_factb   ! solubility factor (frac of aerosol scavenged below cloud) [fraction]
      real(r8), intent(in) :: sol_facti   ! solubility factor (frac of aerosol scavenged in cloud) [fraction]
      real(r8), intent(in) :: sol_factic(pcols,pver)  ! sol_facti for convective clouds [fraction]
         
      real(r8), intent(out) :: fracis(pcols,pver)  ! fraction of species not scavenged [fraction]
      real(r8), intent(out) :: scavt(pcols,pver)   ! scavenging tend [kg/kg/s]
      real(r8), intent(out) :: iscavt(pcols,pver)  ! incloud scavenging tends [kg/kg/s]
      real(r8), intent(out) :: icscavt(pcols,pver)  ! incloud, convective [kg/kg/s]
      real(r8), intent(out) :: isscavt(pcols,pver)  ! incloud, stratiform [kg/kg/s]
      real(r8), intent(out) :: bcscavt(pcols,pver)  ! below cloud, convective [kg/kg/s]
      real(r8), intent(out) :: bsscavt(pcols,pver)  ! below cloud, stratiform [kg/kg/s]
      real(r8), intent(out) :: rcscavt(pcols,pver)  ! resuspension, convective [kg/kg/s]
      real(r8), intent(out) :: rsscavt(pcols,pver)  ! resuspension, stratiform [kg/kg/s] 

      ! local variables
      integer  :: icol          ! column index
      integer  :: kk            ! z index
      real(r8) :: fracev_st     ! fraction of stratiform precip from above that is evaporating [fraction]
      real(r8) :: fracev_cu     ! Fraction of convective precip from above that is evaporating [fraction]
      real(r8) :: fracp         ! fraction of cloud water converted to precip [fraction]
      real(r8) :: precabc(pcols)        ! conv precip from above [kg/m2/s] 
      real(r8) :: precabs(pcols)        ! strat precip from above [kg/m2/s]
      real(r8) :: precabc_base(pcols)   ! conv precip at an effective cloud base for calculations in a particular layer [kg/m2/s]
      real(r8) :: precabs_base(pcols)   ! strat precip at an effective cloud base for calculations in a particular layer [kg/m2/s]
      real(r8) :: precnums_base(pcols)  ! stratiform precip number flux at the bottom of a particular layer [#/m2/s]
      real(r8) :: precnumc_base(pcols)  ! convective precip number flux at the bottom of a particular layer [#/m2/s]
      real(r8) :: scavabs(pcols)        ! stratiform scavenged tracer flux from above [kg/m2/s]
      real(r8) :: scavabc(pcols)        ! convective scavenged tracer flux from above [kg/m2/s]
      real(r8) :: precabx_tmp           ! temporary store precabc or precabs [kg/m2/s]
      real(r8) :: precabx_base_tmp      ! temporarily store precab*_base [kg/m2/s]
      real(r8) :: precnumx_base_tmp     ! temporarily store precnum*_base [#/m2/s]
      real(r8) :: scavabx_tmp   ! temporarily store scavab* [kg/m2/s]
      real(r8) :: resusp_c      ! aerosol mass re-suspension in a particular layer from convective rain [kg/m2/s]
      real(r8) :: resusp_s      ! aerosol mass re-suspension in a particular layer from stratiform rain [kg/m2/s]
      real(r8) :: srcc          ! tendency for convective rain scavenging [kg/kg/s]
      real(r8) :: srcs          ! tendency for stratiform rain scavenging [kg/kg/s]
      real(r8) :: srct          ! total scavenging tendency [kg/kg/s]
      real(r8) :: tracer_incu   ! in-cumulus tracer concentration [kg/kg]
      real(r8) :: tracer_mean   ! mean tracer concenration [kg/kg]
      real(r8) :: tracer_tmp    ! temporarily calculation of tracer [kg/kg]
      real(r8) :: fins          ! fraction of rem. rate by strat rain [fraction]
      real(r8) :: finc          ! fraction of rem. rate by conv. rain [fraction]
      real(r8) :: rat           ! ratio of amount available to amount removed [fraction]
      real(r8) :: arainx        ! precipitation and cloudy volume,at the top interface of current layer [fraction]
      real(r8) :: clddiff       ! cldt - cldc
      real(r8) :: scavt_ik      ! scavenging tend at current (icol,kk) [kg/kg/s]
      real(r8) :: iscavt_ik     ! incloud scavenging tends at current (icol,kk) [kg/kg/s]
      real(r8) :: icscavt_ik    ! incloud, convective scavenging tends at current (icol,kk) [kg/kg/s]
      real(r8) :: isscavt_ik    ! incloud, stratiform scavenging tends at current (icol,kk) [kg/kg/s]
      real(r8) :: bcscavt_ik    ! below cloud, convective scavenging tends at current (icol,kk) [kg/kg/s]
      real(r8) :: bsscavt_ik    ! below cloud, stratiform scavenging tends at current (icol,kk) [kg/kg/s]
      real(r8) :: rcscavt_ik    ! resuspension, convective tends at current (icol,kk) [kg/kg/s]
      real(r8) :: rsscavt_ik    ! resuspension, stratiform tends at current (icol,kk) [kg/kg/s]

      real(r8), parameter :: omsm = 1._r8-2*epsilon(1._r8) ! (1 - small number) used to prevent roundoff errors below zero      

      ! initiate variables
      precabs(:) = 0.0_r8
      precabc(:) = 0.0_r8
      scavabs(:) = 0.0_r8
      scavabc(:) = 0.0_r8

      precabs_base(:) = 0.0_r8
      precabc_base(:) = 0.0_r8
      precnums_base(:) = 0.0_r8
      precnumc_base(:) = 0.0_r8

main_k_loop: &
      do kk = 1,pver
main_i_loop: &
         do icol = 1,ncol
            ! ****************** Evaporation **************************
            ! stratiform
            call compute_evap_frac(                                     &
                        mam_prevap_resusp_optcc,                        & ! in
                        pdel(icol,kk), evaps(icol,kk), precabs(icol),   & ! in
                        fracev_st                                       ) ! out
            ! convective
            call compute_evap_frac(                                     &
                        mam_prevap_resusp_optcc,                        & ! in
                        pdel(icol,kk), evapc(icol,kk), precabc(icol),   & ! in
                        fracev_cu                                       ) ! out

           ! ****************** Scavenging **************************

            ! temporary saved tracer value 
            clddiff = cldt(icol,kk)-cldc(icol,kk)
            tracer_tmp = min(qqcw(icol,kk), tracer(icol,kk) * ( clddiff/max(small_value_2,(1._r8-clddiff)) ) )
            ! calculate in-cumulus and mean tracer values for wetdep_scavenging use
            tracer_incu = f_act_conv(icol,kk) * (tracer(icol,kk) + tracer_tmp)
            tracer_mean = tracer(icol,kk)*(1._r8-cldc(icol,kk)*f_act_conv(icol,kk)) - &
                          tracer_tmp*cldc(icol,kk)*f_act_conv(icol,kk)
            tracer_mean = max(0._r8,tracer_mean)

            ! now do the convective scavenging

            ! fracp: fraction of convective cloud water converted to rain
            ! Sungsu: Below new formula of 'fracp' is necessary since 'conicw'
            ! is a LWC/IWC that has already precipitated out, that is, 'conicw' does
            ! not contain precipitation at all !
            fracp = cmfdqr(icol,kk)*deltat / &
                    max(small_value_12,cldc(icol,kk)*conicw(icol,kk)+(cmfdqr(icol,kk)+dlf(icol,kk))*deltat)
            fracp = min_max_bound(0.0_r8, 1.0_r8, fracp) * cldc(icol,kk)
            call wetdep_scavenging(                                &
                        2,         is_strat_cloudborne,            & ! in, 2 for convective
                        deltat,    fracp,      precabc(icol),      & ! in
                        cldvcu(icol,kk),      scavcoef(icol,kk),   & ! in 
                        sol_factb,            sol_factic(icol,kk), & ! in
                        tracer_incu,          tracer_mean,         & ! in
                        srcc,                 finc                 ) ! out

            ! now do the stratiform scavenging

            ! fracp: fraction of convective cloud water converted to rain
            fracp = precs(icol,kk)*deltat/max(cwat(icol,kk)+precs(icol,kk)*deltat,small_value_12)
            fracp = min_max_bound(0.0_r8, 1.0_r8, fracp)
            call wetdep_scavenging(                                & 
                        1,         is_strat_cloudborne,            & ! in, 1 for stratiform
                        deltat,    fracp,     precabs(icol),       & ! in
                        cldvst(icol,kk),      scavcoef(icol,kk),   & ! in
                        sol_factb,            sol_facti,           & ! in
                        tracer(icol,kk),      tracer_mean,         & ! in
                        srcs,                 fins                 ) ! out


            ! make sure we dont take out more than is there
            ! ratio of amount available to amount removed
            rat = tracer(icol,kk)/max(deltat*(srcc+srcs),small_value_36)
            if (rat.lt.1._r8) then
               srcs = srcs*rat
               srcc = srcc*rat
            endif
            srct = (srcc+srcs)*omsm

            
            ! fraction that is not removed within the cloud
            ! (assumed to be interstitial, and subject to convective transport)
            fracp = deltat*srct/max(cldvst(icol,kk)*tracer(icol,kk),small_value_36)
            fracis(icol,kk) = 1._r8 - min_max_bound(0.0_r8, 1.0_r8, fracp)

            ! ****************** Resuspension **************************

! tend is all tracer removed by scavenging, plus all re-appearing from evaporation above
resusp_block_aa: &
            if ( mam_prevap_resusp_optcc >= 100) then

                ! for stratiform clouds
                arainx = max( cldvst(icol,min(kk+1,pver)), small_value_2 )    ! non-zero
                ! step 1 - do evaporation and resuspension
                call wetdep_resusp(                             &
                        1,              mam_prevap_resusp_optcc,& ! in
                        pdel(icol,kk),  evaps(icol,kk),         & ! in
                        precabs(icol),  precabs_base(icol),     & ! in
                        scavabs(icol),  precnums_base(icol),    & ! in
                        precabx_tmp,    precabx_base_tmp,       & ! out
                        scavabx_tmp,    precnumx_base_tmp,      & ! out
                        resusp_s                                ) ! out
                ! step 2 - do precip production and scavenging
                call wetdep_prevap(                             &
                        1,              mam_prevap_resusp_optcc,& ! in
                        pdel(icol,kk),  precs(icol,kk),         & ! in
                        srcs,           arainx,                 & ! in
                        precabx_tmp,    precabx_base_tmp,       & ! in
                        scavabx_tmp,    precnumx_base_tmp,      & ! in
                        precabs(icol),  precabs_base(icol),     & ! out
                        scavabs(icol),  precnums_base(icol)     ) ! out

                ! for convective clouds
                arainx = max( cldvcu(icol,min(kk+1,pver)), small_value_2 )     ! non-zero
                call wetdep_resusp(                             &
                        2,              mam_prevap_resusp_optcc,& ! in
                        pdel(icol,kk),  evapc(icol,kk),         & ! in
                        precabc(icol),  precabc_base(icol),     & ! in
                        scavabc(icol),  precnumc_base(icol),    & ! in
                        precabx_tmp,    precabx_base_tmp,       & ! out
                        scavabx_tmp,    precnumx_base_tmp,      & ! out
                        resusp_c                                ) ! out
                ! step 2 - do precip production and scavenging
                call wetdep_prevap(                             &
                        2,              mam_prevap_resusp_optcc,& ! in
                        pdel(icol,kk),  cmfdqr(icol,kk),        & ! in
                        srcc,           arainx,                 & ! in
                        precabx_tmp,    precabx_base_tmp,       & ! in
                        scavabx_tmp,    precnumx_base_tmp,      & ! in
                        precabc(icol),  precabc_base(icol),     & ! out
                        scavabc(icol),  precnumc_base(icol)     ) ! out

            else resusp_block_aa ! mam_prevap_resusp_optcc = 0, no resuspension

               resusp_c = fracev_cu*scavabc(icol)
               resusp_s = fracev_st*scavabs(icol)

            endif resusp_block_aa

            ! ****************** update scavengingfor output ***************

            call update_scavenging(                                     &
                   mam_prevap_resusp_optcc,            pdel(icol,kk),   & ! in
                   omsm,   srcc,   srcs,   srct,        fins,   finc,   & ! in
                   fracev_st,  fracev_cu,  resusp_c,    resusp_s,       & ! in
                   precs(icol,kk),         evaps(icol,kk),              & ! in
                   cmfdqr(icol,kk),        evapc(icol,kk),              & ! in
                   scavt_ik,   iscavt_ik,  icscavt_ik,  isscavt_ik,     & ! out
                   bcscavt_ik, bsscavt_ik, rcscavt_ik,  rsscavt_ik,     & ! out
                   scavabs(icol),          scavabc(icol),               & ! inout
                   precabc(icol),          precabs(icol)                ) ! inout

            scavt(icol,kk)   = scavt_ik
            iscavt(icol,kk)  = iscavt_ik
            icscavt(icol,kk) = icscavt_ik
            isscavt(icol,kk) = isscavt_ik
            bcscavt(icol,kk) = bcscavt_ik
            bsscavt(icol,kk) = bsscavt_ik
            rcscavt(icol,kk) = rcscavt_ik
            rsscavt(icol,kk) = rsscavt_ik

        enddo main_i_loop ! End of i = 1, ncol
    enddo main_k_loop ! End of k = 1, pver


   end subroutine wetdepa_v2

!==============================================================================
   subroutine compute_evap_frac(                        &
                        mam_prevap_resusp_optcc,        &
                        pdel_ik,   evap_ik,  precabx,  &
                        fracevx                        )       
! ------------------------------------------------------------------------------
! calculate the fraction of strat precip from above
!                 which evaporates within this layer
! ------------------------------------------------------------------------------
   integer, intent(in) :: mam_prevap_resusp_optcc       ! suspension options
   real(r8),intent(in) :: pdel_ik       ! pressure thikness at current column and level [Pa]
   real(r8),intent(in) :: evap_ik       ! evaporation in this layer [kg/kg/s]
   real(r8),intent(in) :: precabx       ! precipitation from above [kg/m2/s] 
   real(r8),intent(out) :: fracevx      ! fraction of evaporation [fraction]

   if (mam_prevap_resusp_optcc == 0) then
        fracevx = 0.0_r8
   else
        fracevx = evap_ik*pdel_ik/gravit/max(small_value_12,precabx)
        ! trap to ensure reasonable ratio bounds
        fracevx = min_max_bound(0._r8, 1._r8, fracevx)
   endif

   end subroutine compute_evap_frac

!==============================================================================
   subroutine wetdep_scavenging(                                &
                        is_st_cu,       is_strat_cloudborne,    &
                        deltat,         fracp,                  &
                        precabx,       cldv_ik, scavcoef_ik,    &
                        sol_factb,      sol_facti,              &
                        tracer_1,       tracer_2,               &
                        src,            fin                     )
! ------------------------------------------------------------------------------
! do scavenging for both convective and stratiform 
!
! assuming liquid clouds (no ice)
!
! set odds proportional to fraction of the grid box that is swept by the
! precipitation =precabc/rhoh20*(area of sphere projected on plane
!                                /volume of sphere)*deltat
! assume the radius of a raindrop is 1 e-3 m from Rogers and Yau,
! unless the fraction of the area that is cloud is less than odds, in which
! case use the cloud fraction (assumes precabs is in kg/m2/s)
! is really: precabs*3/4/1000./1e-3*deltat
! here I use .1 from Balkanski
! ------------------------------------------------------------------------------

   logical, intent(in) :: is_strat_cloudborne   ! if tracer is stratiform-cloudborne aerosol or not
   integer, intent(in) :: is_st_cu ! options for stratiform (1) or convective (2) clouds

   real(r8), intent(in) :: deltat       ! timestep [s]
   real(r8), intent(in) :: fracp        ! fraction of cloud water converted to precip [fraction]
   real(r8), intent(in) :: precabx     ! precip from above of the layer [kg/m2/s]
   real(r8), intent(in) :: cldv_ik      ! precipitation area at the top interface [fraction]
   real(r8), intent(in) :: scavcoef_ik  ! Dana and Hales coefficient [1/mm]
   real(r8), intent(in) :: sol_factb    ! solubility factor (frac of aerosol scavenged below cloud) [fraction]
   real(r8), intent(in) :: sol_facti    ! solubility factor (frac of aerosol scavenged in cloud) [fraction]
   real(r8), intent(in) :: tracer_1     ! tracer input for calculate src1 [kg/kg]
   real(r8), intent(in) :: tracer_2     ! tracer input for calculate src2 [kg/kg]
   real(r8), intent(out) :: src         ! total scavenging (incloud + belowcloud) [kg/kg/s]
   real(r8), intent(out) :: fin         ! fraction of incloud scavenging [fraction]

   real(r8) :: src1      ! incloud scavenging tendency [kg/kg/s]
   real(r8) :: src2      ! below-cloud scavenging tendency [kg/kg/s]
   real(r8) :: odds      ! limit on removal rate (proportional to prec) [fraction]

   ! calculate limitation of removal rate using Dana and Hales coefficient
   odds = precabx/max(cldv_ik,small_value_5) * scavcoef_ik*deltat
   odds = min_max_bound(0.0_r8, 1.0_r8, odds)

   if ( is_strat_cloudborne ) then
      if (is_st_cu == 2) then
         ! convective cloud does not affect strat-cloudborne aerosol
         src1 = 0._r8
      else
         ! strat in-cloud removal only affects strat-cloudborne aerosol
         ! in-cloud scavenging:
         src1 = sol_facti *fracp*tracer_1/deltat
      endif
      ! no below-cloud scavenging for strat-cloudborne aerosol
      src2 = 0._r8
  
   else
      if (is_st_cu == 2) then      ! convective
         src1 = sol_facti *fracp*tracer_1/deltat
      else      ! stratiform
         ! strat in-cloud removal only affects strat-cloudborne aerosol
         src1 = 0._r8
      endif
      src2 = sol_factb *cldv_ik*odds*tracer_2/deltat

   endif

   src = src1 + src2             ! total stratiform or convective scavenging
   fin=src1/(src + small_value_36)    ! fraction taken by incloud processes

   end subroutine wetdep_scavenging

!==============================================================================
   subroutine wetdep_resusp(                         &
                        is_st_cu,       mam_prevap_resusp_optcc,& ! in
                        pdel_ik,        evapx,                  & ! in
                        precabx_old,    precabx_base_old,       & ! in
                        scavabx_old,    precnumx_base_old,      & ! in
                        precabx_new,    precabx_base_new,       & ! out
                        scavabx_new,    precnumx_base_new,      & ! out
                        resusp_x                                ) ! out

! ------------------------------------------------------------------------------
! do precip production, resuspension and scavenging
! ------------------------------------------------------------------------------
   integer, intent(in) :: is_st_cu ! options for stratiform (1) or convective (2) clouds
                                        ! raindrop size distribution is
                                        ! different for different cloud:
                                        ! 1: assume marshall-palmer distribution
                                        ! 2: assume log-normal distribution
   integer, intent(in) :: mam_prevap_resusp_optcc       ! suspension options
   real(r8),intent(in) :: pdel_ik       ! pressure thikness at current column and level [Pa]
   real(r8),intent(in) :: evapx         ! evaporation at current layer [kg/kg/s]
   real(r8),intent(in) :: precabx_base_old ! input of precipitation at cloud base [kg/m2/s]
   real(r8),intent(in) :: precabx_old ! input of precipitation above this layer [kg/m2/s]
   real(r8),intent(in) :: scavabx_old ! input of scavenged tracer flux from above [kg/m2/s]
   real(r8),intent(in) :: precnumx_base_old ! input of precipitation number at cloud base [#/m2/s]
   real(r8),intent(out) :: precabx_base_new ! output of precipitation at cloud base [kg/m2/s]
   real(r8),intent(out) :: precabx_new ! output of precipitation above this layer [kg/m2/s]
   real(r8),intent(out) :: scavabx_new ! output of scavenged tracer flux from above [kg/m2/s]
   real(r8),intent(out) :: precnumx_base_new ! output of precipitation number at cloud base [#/m2/s]
   real(r8),intent(out) :: resusp_x    ! aerosol mass re-suspension in a particular layer [kg/m2/s]


   ! local variables
   real(r8)            :: tmpa ! temporary working variable

   ! initiate *_new in case they are not calculated
   scavabx_new = scavabx_old
   precnumx_base_new = precnumx_base_old
   precabx_base_new = precabx_base_old

   tmpa = max( 0.0_r8, evapx*pdel_ik/gravit )
   precabx_new = min_max_bound(0.0_r8,precabx_base_new, precabx_old-tmpa)

   if (precabx_new < small_value_30) then
      ! precip rate is essentially zero so do complete resuspension
      call wetdep_resusp_noprecip(                         &
                   is_st_cu,       mam_prevap_resusp_optcc,& ! in
                   precabx_old,    precabx_base_old,       & ! in
                   scavabx_old,    precnumx_base_old,      & ! in
                   precabx_new,    precabx_base_new,       & ! out
                   scavabx_new,    resusp_x                ) ! out

   elseif (evapx <= 0.0_r8) then
      ! no evap so no resuspension
      if ( mam_prevap_resusp_optcc <= 130) then
         scavabx_new = scavabx_old
      endif
      resusp_x = 0.0_r8

   else
      ! regular non-linear resuspension
      call wetdep_resusp_nonlinear(                          &
                     is_st_cu,       mam_prevap_resusp_optcc,& ! in
                     precabx_old,    precabx_base_old,       & ! in
                     scavabx_old,    precnumx_base_old,      & ! in
                     precabx_new,                            & ! in
                     scavabx_new,    resusp_x                ) ! out

   endif

   end subroutine wetdep_resusp


!==============================================================================
   subroutine wetdep_resusp_noprecip(                           &
                        is_st_cu,       mam_prevap_resusp_optcc,& ! in
                        precabx_old,    precabx_base_old,       & ! in
                        scavabx_old,    precnumx_base_old,      & ! in
                        precabx_new,    precabx_base_new,       & ! out
                        scavabx_new,    resusp_x                ) ! out

! ------------------------------------------------------------------------------
! do complete resuspension when precipitation rate is zero
! ------------------------------------------------------------------------------
   integer, intent(in) :: is_st_cu      ! options for stratiform (1) or convective (2) clouds
                                        ! raindrop size distribution is
                                        ! different for different cloud:
                                        ! 1: assume marshall-palmer distribution
                                        ! 2: assume log-normal distribution
   integer, intent(in) :: mam_prevap_resusp_optcc       ! suspension options
   real(r8),intent(in) :: precabx_base_old ! input of precipitation at cloud base [kg/m2/s]
   real(r8),intent(in) :: precabx_old ! input of precipitation above this layer [kg/m2/s]
   real(r8),intent(in) :: scavabx_old ! input of scavenged tracer flux from above [kg/m2/s]
   real(r8),intent(in) :: precnumx_base_old ! precipitation number at cloud base [#/m2/s]
   real(r8),intent(out) :: precabx_base_new ! output of precipitation at cloud base [kg/m2/s]
   real(r8),intent(out) :: precabx_new ! output of precipitation above this layer [kg/m2/s]
   real(r8),intent(out) :: scavabx_new ! output of scavenged tracer flux from above [kg/m2/s]
   real(r8),intent(out) :: resusp_x    ! aerosol mass re-suspension in a particular layer [kg/m2/s]


   ! local variables
   real(r8)     :: u_old,u_new  ! fraction of precabx and precabx_base
   real(r8)     :: x_old,x_new  ! fraction after calling function *_resusp_vs_fprec_evap_mpln

   if ( mam_prevap_resusp_optcc <= 130) then
        ! linear resuspension based on scavenged aerosol mass or number
        scavabx_new = 0.0_r8
        resusp_x = scavabx_old
   else
        if (precabx_base_old < small_value_30) then
            resusp_x = 0.0_r8
        else
            ! non-linear resuspension of aerosol number based on raindrop number
            u_old = min_max_bound(0.0_r8, 1.0_r8, precabx_old/precabx_base_old)
            x_old = 1.0_r8 - fprecn_resusp_vs_fprec_evap_mpln(1.0_r8-u_old, is_st_cu)
            x_old = min_max_bound(0.0_r8, 1.0_r8, x_old)
            x_new = 0.0_r8
            resusp_x = max( 0.0_r8, precnumx_base_old*(x_old - x_new) )
        endif
   endif
   ! setting both these precip rates to zero causes the resuspension
   ! calculations to start fresh if there is any more precip production
   precabx_new = 0.0_r8
   precabx_base_new = 0.0_r8

   end subroutine wetdep_resusp_noprecip

!==============================================================================
   subroutine wetdep_resusp_nonlinear(                          &
                        is_st_cu,       mam_prevap_resusp_optcc,& ! in
                        precabx_old,    precabx_base_old,       & ! in
                        scavabx_old,    precnumx_base_old,      & ! in
                        precabx_new,                            & ! in
                        scavabx_new,    resusp_x                ) ! out

! ------------------------------------------------------------------------------
! do nonlinear resuspension of aerosol mass or number
! ------------------------------------------------------------------------------
   integer, intent(in) :: is_st_cu      ! options for stratiform (1) or convective (2) clouds
                                        ! raindrop size distribution is
                                        ! different for different cloud:
                                        ! 1: assume marshall-palmer distribution
                                        ! 2: assume log-normal distribution
   integer, intent(in) :: mam_prevap_resusp_optcc       ! suspension options
   real(r8),intent(in) :: precabx_base_old ! input of precipitation at cloud base [kg/m2/s]
   real(r8),intent(in) :: precabx_old  ! input of precipitation above this layer [kg/m2/s]
   real(r8),intent(in) :: scavabx_old  ! input scavenged tracer flux from above [kg/m2/s]
   real(r8),intent(in) :: precnumx_base_old ! precipitation number at cloud base [#/m2/s]
   real(r8),intent(in) :: precabx_new  ! output of precipitation above this layer [kg/m2/s]
   real(r8),intent(out) :: scavabx_new ! output scavenged tracer flux from above [kg/m2/s]
   real(r8),intent(out) :: resusp_x    ! aerosol mass re-suspension in a particular layer [kg/m2/s]


   ! local variables
   real(r8)     :: u_old,u_new  ! fraction of precabx and precabx_base
   real(r8)     :: x_old,x_new  ! fraction after calling function *_resusp_vs_fprec_evap_mpln
   real(r8)     :: x_ratio      ! fraction of x_tmp/x_old

   u_old = min_max_bound(0.0_r8, 1.0_r8, precabx_old/precabx_base_old)
   if ( mam_prevap_resusp_optcc <= 130) then
        ! non-linear resuspension of aerosol mass
        x_old = 1.0_r8 - faer_resusp_vs_fprec_evap_mpln( 1.0_r8-u_old, is_st_cu)
   else
        ! non-linear resuspension of aerosol number based on raindrop number
        x_old = 1.0_r8 - fprecn_resusp_vs_fprec_evap_mpln(1.0_r8-u_old, is_st_cu)
   endif
   x_old = min_max_bound(0.0_r8, 1.0_r8, x_old)

   if (x_old < small_value_30) then
        x_new = 0.0_r8
        x_ratio = 0.0_r8
   else
        u_new = min_max_bound(0.0_r8, 1.0_r8, precabx_new/precabx_base_old)
        u_new = min( u_new, u_old )
        if ( mam_prevap_resusp_optcc <= 130) then
            ! non-linear resuspension of aerosol mass
            x_new = 1.0_r8 - faer_resusp_vs_fprec_evap_mpln(1.0_r8-u_new, is_st_cu)
        else
            ! non-linear resuspension of aerosol number based on raindrop number
            x_new = 1.0_r8 - fprecn_resusp_vs_fprec_evap_mpln(1.0_r8-u_new, is_st_cu)
        endif
        x_new = min_max_bound(0.0_r8, 1.0_r8, x_new)
        x_new = min( x_new, x_old )
        x_ratio = min_max_bound(0.0_r8, 1.0_r8, x_new/x_old)
   endif

   ! update aerosol resuspension
   if ( mam_prevap_resusp_optcc <= 130) then
        ! aerosol mass resuspension
        scavabx_new = max( 0.0_r8, scavabx_old * x_ratio )
        resusp_x = max( 0.0_r8, scavabx_old - scavabx_new )
   else
        ! number resuspension
        resusp_x = max( 0.0_r8, precnumx_base_old*(x_old - x_new) )
   endif


   end subroutine wetdep_resusp_nonlinear

!==============================================================================
   subroutine wetdep_prevap(                                    &
                        is_st_cu,       mam_prevap_resusp_optcc,& ! in
                        pdel_ik,        pprdx,  srcx,   arainx, & ! in
                        precabx_old,    precabx_base_old,       & ! in
                        scavabx_old,    precnumx_base_old,      & ! in
                        precabx_new,    precabx_base_new,       & ! out
                        scavabx_new,    precnumx_base_new       ) ! out

! ------------------------------------------------------------------------------
! do precip production and scavenging
! ------------------------------------------------------------------------------
   integer, intent(in) :: is_st_cu      ! options for stratiform (1) or convective (2) clouds
                                        ! raindrop size distribution is
                                        ! different for different cloud:
                                        ! 1: assume marshall-palmer distribution
                                        ! 2: assume log-normal distribution
   integer, intent(in) :: mam_prevap_resusp_optcc       ! suspension options
   real(r8),intent(in) :: pdel_ik       ! pressure thikness at current column and level [Pa]
   real(r8),intent(in) :: pprdx  ! precipitation generation rate [kg/kg/s]
   real(r8),intent(in) :: srcx   ! scavenging tendency [kg/kg/s]
   real(r8),intent(in) :: arainx ! precipitation and cloudy volume,at the top interface of current layer [fraction]
   real(r8),intent(in) :: precabx_base_old ! input of precipitation at cloud base [kg/m2/s]
   real(r8),intent(in) :: precabx_old ! input of precipitation above this layer [kg/m2/s]
   real(r8),intent(in) :: scavabx_old ! input scavenged tracer flux from above [kg/m2/s]
   real(r8),intent(in) :: precnumx_base_old ! input of rain number at cloud base [#/m2/s]
   real(r8),intent(out) :: precabx_base_new ! output of precipitation at cloud base [kg/m2/s]
   real(r8),intent(out) :: precabx_new ! output of precipitation above this layer [kg/m2/s]
   real(r8),intent(out) :: scavabx_new ! output scavenged tracer flux from above [kg/m2/s]
   real(r8),intent(out) :: precnumx_base_new ! output of rain number at cloud base [#/m2/s]

   ! local variables
   real(r8)            :: tmpa ! temporary working variable
 
   ! initiate *_new in case they are not calculated
   ! precabx_base_new and precabx_new are always calculated
   scavabx_new = scavabx_old
   precnumx_base_new = precnumx_base_old

   tmpa = max( 0.0_r8, pprdx*pdel_ik/gravit )
   precabx_base_new = max( 0.0_r8, precabx_base_old + tmpa )
   precabx_new = min_max_bound(0.0_r8, precabx_base_new,  precabx_old+tmpa)

   if ( mam_prevap_resusp_optcc <= 130) then
       ! aerosol mass scavenging
       tmpa = max( 0.0_r8, srcx*pdel_ik/gravit )
       scavabx_new = max( 0.0_r8, scavabx_old + tmpa )
   else
       ! raindrop number increase
       if (precabx_base_new < small_value_30) then
          precnumx_base_new = 0.0_r8
       elseif (precabx_base_new > precabx_base_old) then
          ! note - calc rainshaft number flux from rainshaft water flux,
          ! then multiply by rainshaft area to get grid-average number flux
          tmpa = arainx * flux_precnum_vs_flux_prec_mpln((precabx_base_new/arainx), is_st_cu)
          precnumx_base_new = max( 0.0_r8, tmpa )
       else
          precnumx_base_new = precnumx_base_old
       endif
   endif

   end subroutine wetdep_prevap

!==============================================================================
   subroutine update_scavenging(                                &
             mam_prevap_resusp_optcc,   pdel_ik,                & ! in
             omsm,   srcc,   srcs,      srct,    fins,   finc,  & ! in
             fracev_st, fracev_cu,      resusp_c,   resusp_s,   & ! in
             precs_ik,  evaps_ik,       cmfdqr_ik,  evapc_ik,   & ! in
             scavt_ik,  iscavt_ik,      icscavt_ik, isscavt_ik, & ! out
             bcscavt_ik,bsscavt_ik,     rcscavt_ik, rsscavt_ik, & ! out
             scavabs,   scavabc,        precabc,    precabs     ) ! inout
! ------------------------------------------------------------------------------
! update scavenging variables
! *_ik are variables at the grid (icol, kk)
! ------------------------------------------------------------------------------
   ! input variables
   integer, intent(in) :: mam_prevap_resusp_optcc       ! suspension options 
   real(r8),intent(in) :: pdel_ik       ! pressure thikness [Pa]
   real(r8),intent(in) :: omsm          ! 1 - (a small number), to prevent roundoff errors below zero
   real(r8),intent(in) :: srcc          ! tend for convective rain scavenging [kg/kg/s]
   real(r8),intent(in) :: srcs          ! tend for stratiform rain scavenging [kg/kg/s]
   real(r8),intent(in) :: srct          ! total scavenging tendency for conv+strat rain [kg/kg/s]
   real(r8),intent(in) :: fins          ! fraction of rem. rate by strat rain [fraction]
   real(r8),intent(in) :: finc          ! fraction of rem. rate by conv. rain [fraction]
   real(r8),intent(in) :: fracev_st     ! fraction of stratiform precip from above that is evaporating [fraction]
   real(r8),intent(in) :: fracev_cu     ! Fraction of convective precip from above that is evaporating [fraction]
   real(r8),intent(in) :: resusp_c      ! aerosol mass re-suspension in a particular layer from convective rain [kg/m2/s]
   real(r8),intent(in) :: resusp_s      ! aerosol mass re-suspension in a particular layer from stratiform rain [kg/m2/s]
   real(r8),intent(in) :: precs_ik      ! rate of production of stratiform precip [kg/kg/s]
   real(r8),intent(in) :: evaps_ik      ! rate of evaporation of precip [kg/kg/s]
   real(r8),intent(in) :: cmfdqr_ik     ! rate of production of convective precip [kg/kg/s]
   real(r8),intent(in) :: evapc_ik      ! Evaporation rate of convective precipitation [kg/kg/s]
   ! output variables
   real(r8), intent(out) :: scavt_ik    ! scavenging tend [kg/kg/s]
   real(r8), intent(out) :: iscavt_ik   ! incloud scavenging tends [kg/kg/s]
   real(r8), intent(out) :: icscavt_ik  ! incloud, convective [kg/kg/s]
   real(r8), intent(out) :: isscavt_ik  ! incloud, stratiform [kg/kg/s]
   real(r8), intent(out) :: bcscavt_ik  ! below cloud, convective [kg/kg/s]
   real(r8), intent(out) :: bsscavt_ik  ! below cloud, stratiform [kg/kg/s]
   real(r8), intent(out) :: rcscavt_ik  ! resuspension, convective [kg/kg/s]
   real(r8), intent(out) :: rsscavt_ik  ! resuspension, stratiform [kg/kg/s]
   real(r8), intent(inout) :: scavabs   ! stratiform scavenged tracer flux from above [kg/m2/s]
   real(r8), intent(inout) :: scavabc   ! convective scavenged tracer flux from above [kg/m2/s]
   real(r8), intent(inout) :: precabc   ! conv precip from above [kg/m2/s]
   real(r8), intent(inout) :: precabs   ! strat precip from above [kg/m2/s]

   if ( mam_prevap_resusp_optcc == 0) then
       scavt_ik = -srct + (fracev_st*scavabs + fracev_cu*scavabc)*gravit/pdel_ik
   else
       scavt_ik = -srct + (resusp_s+resusp_c)*gravit/pdel_ik
   endif

   iscavt_ik = -(srcc*finc + srcs*fins)*omsm
   icscavt_ik = -(srcc*finc) * omsm
   isscavt_ik = -(srcs*fins) * omsm

   if (mam_prevap_resusp_optcc == 0) then
          bcscavt_ik = -(srcc*(1-finc)) * omsm + fracev_cu*scavabc*gravit/pdel_ik
          bsscavt_ik = -(srcs*(1-fins)) * omsm + fracev_st*scavabs*gravit/pdel_ik
          rcscavt_ik = 0.0
          rsscavt_ik = 0.0
   else ! here mam_prevap_resusp_optcc == 130, 210, 230
          bcscavt_ik = -(srcc * (1-finc)) * omsm
          rcscavt_ik = resusp_c*gravit/pdel_ik
          bsscavt_ik = -(srcs * (1-fins)) * omsm
          rsscavt_ik = resusp_s*gravit/pdel_ik
   endif

   ! now keep track of scavenged mass and precip
   if (mam_prevap_resusp_optcc == 0) then
       scavabs = scavabs*(1-fracev_st) + srcs*pdel_ik/gravit
       scavabc = scavabc*(1-fracev_cu) + srcc*pdel_ik/gravit
       precabs = precabs + (precs_ik  - evaps_ik)*pdel_ik/gravit
       precabc = precabc + (cmfdqr_ik - evapc_ik)*pdel_ik/gravit
   endif

   end subroutine update_scavenging

!==============================================================================
      function flux_precnum_vs_flux_prec_mpln( flux_prec, jstrcnv )
! --------------------------------------------------------------------------------
!     flux_precnum_vs_flux_prec_mp = precipitation number flux at the cloud base [drops/m^2/s]
!     Options of assuming log-normal or marshall-palmer raindrop size distribution
! --------------------------------------------------------------------------------
      real(r8) :: flux_precnum_vs_flux_prec_mpln  ! [drops/m^2/s]
      real(r8), intent(in) :: flux_prec     ! [drops/m^2/s]
      integer,  intent(in) :: jstrcnv   ! current only two options: 1 for marshall-palmer distribution, 2 for log-normal distribution

      real(r8) :: a0, a1
      real(r8) :: x_var, y_var

      if (jstrcnv <= 1) then
        ! marshall-palmer distribution
        a0 =  1.0885896550304022E+01_r8
        a1 =  4.3660645528167907E-01_r8
      else
        ! log-normal distribution
        a0 =  9.9067806476181524E+00_r8
        a1 =  4.2690709912134056E-01_r8
      endif

      if (flux_prec >= small_value_36) then
         x_var = log( flux_prec )
         y_var = exp( a0 + a1*x_var )
      else
         y_var = 0.0_r8
      endif
      flux_precnum_vs_flux_prec_mpln = y_var

      return
      end function flux_precnum_vs_flux_prec_mpln


!==============================================================================
      function faer_resusp_vs_fprec_evap_mpln( fprec_evap, jstrcnv )
! --------------------------------------------------------------------------------
! corresponding fraction of precipitation-borne aerosol flux that is resuspended
! Options of assuming log-normal or marshall-palmer raindrop size distribution
! note that these fractions are relative to the cloud-base fluxes,
! and not to the layer immediately above fluxes
! --------------------------------------------------------------------------------

      real(r8) :: faer_resusp_vs_fprec_evap_mpln ! [fraction]
      real(r8), intent(in) :: fprec_evap ! [fraction]
      integer,  intent(in) :: jstrcnv   ! current only two options: 1 for marshall-palmer distribution, 2 for log-normal distribution

      real(r8) :: a01, a02, a03, a04, a05, a06, a07, a08, a09, x_lox_lin, y_lox_lin
      real(r8) :: x_var, y_var

      if (jstrcnv <= 1) then
        ! marshall-palmer distribution
        a01 =  8.6591133737322856E-02_r8
        a02 = -1.7389168499601941E+00_r8
        a03 =  2.7401882373663732E+01_r8
        a04 = -1.5861714653209464E+02_r8
        a05 =  5.1338179363011193E+02_r8
        a06 = -9.6835933124501412E+02_r8
        a07 =  1.0588489932213311E+03_r8
        a08 = -6.2184513459217271E+02_r8
        a09 =  1.5184126886039758E+02_r8
        x_lox_lin =  5.0000000000000003E-02_r8
        y_lox_lin =  2.5622471203221014E-03_r8
      else
        ! log-normal distribution
        a01 =  6.1944215103685640E-02_r8
        a02 = -2.0095166685965378E+00_r8
        a03 =  2.3882460251821236E+01_r8
        a04 = -1.2695611774753374E+02_r8
        a05 =  4.0086943562320101E+02_r8
        a06 = -7.4954272875943707E+02_r8
        a07 =  8.1701055892023624E+02_r8
        a08 = -4.7941894659538502E+02_r8
        a09 =  1.1710291076059025E+02_r8
        x_lox_lin =  1.0000000000000001E-01_r8
        y_lox_lin =  6.2227889828044350E-04_r8
      endif

      x_var = min_max_bound(0.0_r8, 1.0_r8, fprec_evap)
      if (x_var < x_lox_lin) then
         y_var = y_lox_lin * (x_var/x_lox_lin)
      else
         y_var = x_var*( a01 + x_var*( a02 + x_var*( a03 + x_var*( a04 + x_var*( a05 &
           + x_var*( a06 + x_var*( a07 + x_var*( a08 + x_var*a09 ))))))))
      endif
      faer_resusp_vs_fprec_evap_mpln = y_var

      return
      end function faer_resusp_vs_fprec_evap_mpln


!==============================================================================
      function fprecn_resusp_vs_fprec_evap_mpln( fprec_evap, jstrcnv )
! --------------------------------------------------------------------------------
! Rain number evaporation fraction
! Options of assuming log-normal or marshall-palmer raindrop size distribution
! note that these fractions are relative to the cloud-base fluxes,
! and not to the layer immediately above fluxes
! --------------------------------------------------------------------------------
      real(r8) :: fprecn_resusp_vs_fprec_evap_mpln  ! [fraction]
      real(r8), intent(in) :: fprec_evap     ! [fraction]
      integer,  intent(in) :: jstrcnv  ! current only two options: 1 for marshall-palmer distribution, 2 for log-normal distribution

      real(r8) :: a01, a02, a03, a04, a05, a06, a07, a08, a09, x_lox_lin, y_lox_lin
      real(r8) :: x_var, y_var

      if (jstrcnv <= 1) then
        !marshall-palmer distribution
        a01 =  4.5461070198414655E+00_r8
        a02 = -3.0381753620077529E+01_r8
        a03 =  1.7959619926085665E+02_r8
        a04 = -6.7152282193785618E+02_r8
        a05 =  1.5651931323557126E+03_r8
        a06 = -2.2743927701175126E+03_r8
        a07 =  2.0004645897056735E+03_r8
        a08 = -9.7351466279626209E+02_r8
        a09 =  2.0101198012962413E+02_r8
        x_lox_lin =  5.0000000000000003E-02_r8
        y_lox_lin =  1.7005858490684875E-01_r8
      else
        ! log-normal distribution
        a01 = -5.2335291116884175E-02_r8
        a02 =  2.7203158069178226E+00_r8
        a03 =  9.4730878152409375E+00_r8
        a04 = -5.0573187592544798E+01_r8
        a05 =  9.4732631441282862E+01_r8
        a06 = -8.8265926556465814E+01_r8
        a07 =  3.5247835268269142E+01_r8
        a08 =  1.5404586576716444E+00_r8
        a09 = -3.8228795492549068E+00_r8
        x_lox_lin =  1.0000000000000001E-01_r8
        y_lox_lin =  2.7247994766566485E-02_r8
      endif

      x_var = min_max_bound(0.0_r8, 1.0_r8, fprec_evap)
      if (x_var < x_lox_lin) then
         y_var = y_lox_lin * (x_var/x_lox_lin)
      else
         y_var = x_var*( a01 + x_var*( a02 + x_var*( a03 + x_var*( a04 + x_var*( a05 &
           + x_var*( a06 + x_var*( a07 + x_var*( a08 + x_var*a09 ))))))))
      endif
      fprecn_resusp_vs_fprec_evap_mpln = y_var

      return
      end function fprecn_resusp_vs_fprec_evap_mpln

!##############################################################################

end module wetdep
