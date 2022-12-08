module wetdep

!----------------------------------------------------------------------- 
!
! Wet deposition routines for both aerosols and gas phase constituents.
! 
!-----------------------------------------------------------------------

use shr_kind_mod, only: r8 => shr_kind_r8
use ppgrid,       only: pcols, pver
use physconst,    only: gravit, rair, tmelt
use phys_control, only: cam_physpkg_is
use cam_logfile,  only: iulog
use cam_abortutils, only: endrun
use spmd_utils,     only: masterproc
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

real(r8), parameter :: cmftau = 3600._r8
real(r8), parameter :: rhoh2o = 1000._r8            ! density of water
real(r8), parameter :: molwta = 28.97_r8            ! molecular weight dry air gm/mole

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


  use phys_control,   only: cam_physpkg_is
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
                              max(0.01_r8, sh_frac(:ncol,:) + dp_frac(:ncol,:))
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
   real(r8) sumppr_all(pcols,pver)    ! precipitation rate in all vertical levels [kg/m2/s]
   real(r8) lprec(pcols,pver)         ! local production rate of precip [kg/m2/s]
   real(r8) sumppr_cu_all(pcols,pver) ! same as sumppr_all but for conv.precip. calculated but not used
   real(r8) lprec_cu(pcols,pver)      ! Local production rate of convective precip [kg/m2/s]
   real(r8) sumppr_st_all(pcols,pver) ! same as sumppr_all but for strat.precip. calculated but not used
   real(r8) lprec_st(pcols,pver)      ! Local production rate of stratiform precip [kg/m2/s]
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
    integer  icol, kk

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
   integer  icol,kk
   real(r8) sumppr(pcols)        ! precipitation rate [kg/m2/s]
   real(r8) sumpppr(pcols)       ! sum of positive precips from above
   real(r8) cldv1(pcols)         ! precip weighted cloud fraction from above [kg/m2/s]
   real(r8) lprecp               ! local production rate of precip [kg/m2/s] if positive


   ! initiate variables
   do icol=1,ncol
      sumppr(icol) = 0._r8
      cldv1(icol) = 0._r8
      sumpppr(icol) = 1.e-36_r8   ! not 0 because it will be divided
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
         lprecp = max(lprec(icol,kk), 1.e-30_r8)
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
   ! of the equation  -- by Shuaiqi when refactoring
   convfw = 1.94_r8*2.13_r8*sqrt(rhoh2o*gravit*2.7e-4_r8)
   do kk = 1,pver
     do icol = 1,ncol
         rain(icol,kk) = 0._r8
         if(temperature(icol,kk) .gt. tmelt) then
            rho = pmid(icol,kk)/(rair*temperature(icol,kk))
            vfall = convfw/sqrt(rho)
            rain(icol,kk) = sumppr(icol,kk) / (rho*vfall)
            if (rain(icol,kk).lt.1.e-14_r8) rain(icol,kk) = 0._r8
         endif
      end do
   end do

end subroutine rain_mix_ratio

!==============================================================================

! REASTER 08/05/2015
! changed arguments
!    put them in a more logical order
!    optional arguments are now mandatory, and commented out:
!       all "if ( present(xx) )" tests
!       any code for ".not. present(xx)" cases
!    eliminated the sol_fact**_in - now just use sol_fact**

! old argument order
! ubroutine wetdepa_v2( t, p, q, pdel, &
!                       cldt, cldc, cmfdqr, evapc, conicw, precs, conds, &
!                       evaps, cwat, tracer, deltat, &
!                       scavt, iscavt, cldv, cldvcu, cldvst, dlf, fracis, sol_fact, ncol, &
!                       scavcoef, is_strat_cloudborne, rate1ord_cw2pr_st, qqcw, f_act_conv, &
!                       icscavt, isscavt, bcscavt, bsscavt, rcscavt, rsscavt, &  
!                       sol_facti_in, sol_factbi_in, sol_factii_in, &
!                       sol_factic_in, sol_factiic_in, resus_fix ) 
! note - p, q, cldv not needed

! new argument order
subroutine wetdepa_v2( ncol, deltat, &
                       t, p, q, pdel, &
                       cmfdqr, evapc, dlf, conicw, &
                       precs, conds, evaps, cwat, &
                       cldt, cldc, cldv, cldvcu, cldvst, &
                       sol_factb, sol_factbi, sol_facti, sol_factii, sol_factic, sol_factiic, &
                       mam_prevap_resusp_optcc, is_strat_cloudborne, scavcoef, rate1ord_cw2pr_st, f_act_conv, &
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
      !-----------------------------------------------------------------------

      use phys_control, only: phys_getopts

      implicit none

      integer, intent(in) :: ncol

      real(r8), intent(in) ::&
         deltat,               &! time step
         t(pcols,pver),        &! temperature
         p(pcols,pver),        &! pressure
         q(pcols,pver),        &! moisture
         pdel(pcols,pver),     &! pressure thikness
         cmfdqr(pcols,pver),   &! rate of production of convective precip
! Sungsu
         evapc(pcols,pver),    &! Evaporation rate of convective precipitation
         dlf(pcols,pver),      &! Detrainment of convective condensate [kg/kg/s]
! Sungsu
         conicw(pcols,pver),   &! convective cloud water
         precs(pcols,pver),    &! rate of production of stratiform precip
         conds(pcols,pver),    &! rate of production of condensate
         evaps(pcols,pver),    &! rate of evaporation of precip
         cwat(pcols,pver),     &! cloud water amount 
         cldt(pcols,pver),     &! total cloud fraction
         cldc(pcols,pver),     &! convective cloud fraction
         cldv(pcols,pver),     &! total cloud fraction
! Sungsu
         cldvcu(pcols,pver),   &! Convective precipitation area at the top interface of each layer
         cldvst(pcols,pver),   &! Stratiform precipitation area at the top interface of each layer
! Sungsu
         tracer(pcols,pver)     ! trace species
      ! If subroutine is called with just sol_fact:
            ! sol_fact is used for both in- and below-cloud scavenging
      ! If subroutine is called with optional argument sol_facti_in:
            ! sol_fact  is used for below cloud scavenging
            ! sol_facti is used for in cloud scavenging

      integer, intent(in) :: mam_prevap_resusp_optcc ! suspension options.
!     0 = no resuspension
!     1 = linear resuspension of aerosol mass or number following original mam
!     coding and history_aero_prevap_resusp = .false.
!     2 = same as 1 but history_aero_prevap_resusp = .true.
!     3 = same as 2 but with some added "xxx = max( 0, xxx)" lines
!   130 = non-linear resuspension of aerosol mass   based on scavenged aerosol mass
!   230 = non-linear resuspension of aerosol number based on raindrop number
!     (1,2,3 are not used in the current code)

!     logical, intent(in) :: resus_fix
      ! rce 2010/05/01
      ! is_strat_cloudborne = .true. if tracer is stratiform-cloudborne aerosol; else .false. 
      logical, intent(in) :: is_strat_cloudborne   
      real(r8), intent(in) :: scavcoef(pcols,pver) ! Dana and Hales coefficient (/mm) (0.1 if not MODAL_AERO)
      ! rate1ord_cw2pr_st = 1st order rate for strat cw to precip (1/s) 
      real(r8), intent(in) :: rate1ord_cw2pr_st(pcols,pver)
      ! qqcw = strat-cloudborne aerosol corresponding to tracer when is_strat_cloudborne==.false.; else 0.0 
      ! f_act_conv = conv-cloud activation fraction when is_strat_cloudborne==.false.; else 0.0 
      real(r8), intent(in) :: f_act_conv(pcols,pver)

      real(r8), intent(in) :: qqcw(pcols,pver)
      ! end rce 2010/05/01

!     real(r8), intent(in) :: sol_fact ! solubility factor (fraction of aer scavenged below & in, or just below or sol_facti is provided)
      real(r8), intent(in) :: sol_factb   ! solubility factor (frac of aerosol scavenged below cloud)
      real(r8), intent(in) :: sol_factbi  ! solubility factor (frac of aerosol scavenged below cloud by ice)
      real(r8), intent(in) :: sol_facti   ! solubility factor (frac of aerosol scavenged in cloud)
      real(r8), intent(in) :: sol_factii  ! solubility factor (frac of aerosol scavenged in cloud by ice)
      real(r8), intent(in) :: sol_factic(pcols,pver)  ! sol_facti for convective clouds
      real(r8), intent(in) :: sol_factiic ! sol_factii for convective clouds
         
      real(r8), intent(out) :: fracis(pcols,pver)  ! fraction of species not scavenged
      real(r8), intent(out) :: scavt(pcols,pver)   ! scavenging tend 
      real(r8), intent(out) :: iscavt(pcols,pver)  ! incloud scavenging tends

      real(r8), intent(out) :: icscavt(pcols,pver)  ! incloud, convective
      real(r8), intent(out) :: isscavt(pcols,pver)  ! incloud, stratiform
      real(r8), intent(out) :: bcscavt(pcols,pver)  ! below cloud, convective
      real(r8), intent(out) :: bsscavt(pcols,pver)  ! below cloud, stratiform
      real(r8), intent(out) :: rcscavt(pcols,pver)  ! resuspension, convective 
      real(r8), intent(out) :: rsscavt(pcols,pver)  ! resuspension, stratiform 

      ! local variables

      integer i                 ! x index
      integer k                 ! z index

      real(r8) adjfac               ! factor stolen from cmfmca
      real(r8) aqfrac               ! fraction of tracer in aqueous phase
      real(r8) cwatc                ! local convective total water amount 
      real(r8) cwats                ! local stratiform total water amount 
      real(r8) cwatp                ! local water amount falling from above precip
      real(r8) fracev(pcols)        ! fraction of precip from above that is evaporating
! Sungsu
      real(r8) fracev_cu(pcols)     ! Fraction of convective precip from above that is evaporating
! Sungsu
      real(r8) fracp                ! fraction of cloud water converted to precip
      real(r8) gafrac               ! fraction of tracer in gas phasea
      real(r8) hconst               ! henry's law solubility constant when equation is expressed
                                ! in terms of mixing ratios
      real(r8) mpla                 ! moles / liter H2O entering the layer from above
      real(r8) mplb                 ! moles / liter H2O leaving the layer below
      real(r8) omsm                 ! 1 - (a small number)
      real(r8) part                 !  partial pressure of tracer in atmospheres
      real(r8) patm                 ! total pressure in atmospheres
      real(r8) pdog                 ! work variable (pdel/gravit)
      real(r8) precabc(pcols)       ! conv precip from above (work array)
      real(r8) precabs(pcols)       ! strat precip from above (work array)
      real(r8) precbl               ! precip falling out of level (work array)
      real(r8) precmin              ! minimum convective precip causing scavenging
      real(r8) rat(pcols)           ! ratio of amount available to amount removed
      real(r8) scavab(pcols)        ! scavenged tracer flux from above (work array)
      real(r8) scavabc(pcols)       ! scavenged tracer flux from above (work array)
      real(r8) srcc                 ! tend for convective rain
      real(r8) srcs                 ! tend for stratiform rain
      real(r8) srct(pcols)          ! work variable
      real(r8) tracab(pcols)        ! column integrated tracer amount
!      real(r8) vfall                ! fall speed of precip
      real(r8) fins                 ! fraction of rem. rate by strat rain
      real(r8) finc                 ! fraction of rem. rate by conv. rain
      real(r8) srcs1                ! work variable
      real(r8) srcs2                ! work variable
      real(r8) tc                   ! temp in celcius
      real(r8) weight               ! fraction of condensate which is ice
      real(r8) cldmabs(pcols)       ! maximum cloud at or above this level
      real(r8) cldmabc(pcols)       ! maximum cloud at or above this level
      real(r8) odds                 ! limit on removal rate (proportional to prec)
      real(r8) dblchek(pcols)

    ! Jan.16.2009. Sungsu for wet scavenging below clouds.
    ! real(r8) cldovr_cu(pcols)     ! Convective precipitation area at the base of each layer
    ! real(r8) cldovr_st(pcols)     ! Stratiform precipitation area at the base of each layer

      real(r8) tracer_incu
      real(r8) tracer_mean

    ! End by Sungsu

!     real(r8) sol_facti,  sol_factb  ! in cloud and below cloud fraction of aerosol scavenged
!     real(r8) sol_factii, sol_factbi ! in cloud and below cloud fraction of aerosol scavenged by ice
!     real(r8) sol_factic(pcols,pver)             ! sol_facti for convective clouds
!     real(r8) sol_factiic            ! sol_factii for convective clouds
      ! sol_factic & solfact_iic added for MODAL_AERO.  
      ! For stratiform cloud, cloudborne aerosol is treated explicitly,
      !    and sol_facti is 1.0 for cloudborne, 0.0 for interstitial.
      ! For convective cloud, cloudborne aerosol is not treated explicitly,
      !    and sol_factic is 1.0 for both cloudborne and interstitial.

      integer  jstrcnv

      real(r8), parameter :: prec_smallaa = 1.0e-30_r8  ! 1e-30 kg/m2/s (or mm/s) = 3.2e-23 mm/yr
      real(r8), parameter :: x_smallaa = 1.0e-30_r8

      real(r8) arainx
      real(r8) evapx
      real(r8) pprdx
      real(r8) precabc_base(pcols)  ! conv precip at an effective cloud base for calculations in a particular layer
      real(r8) precabs_base(pcols)  ! strat precip at an effective cloud base for calculations in a particular layer
      real(r8) precabx_old, precabx_tmp, precabx_new
      real(r8) precabx_base_old, precabx_base_tmp, precabx_base_new
      real(r8) precnums_base(pcols)  ! stratiform precip number flux at the bottom of a particular layer
      real(r8) precnumc_base(pcols)  ! convective precip number flux at the bottom of a particular layer
      real(r8) precnumx_base_old, precnumx_base_tmp, precnumx_base_new
      real(r8) resusp_c          ! aerosol mass re-suspension in a particular layer from convective rain
      real(r8) resusp_s          ! aerosol mass re-suspension in a particular layer from stratiform rain

      real(r8) resusp_x
      real(r8) resusp_c_sv(pcols)
      real(r8) resusp_s_sv(pcols)
      real(r8) scavabx_old, scavabx_tmp, scavabx_new
      real(r8) srcx
      real(r8) tmpa, tmpb
      real(r8) u_old, u_tmp
      real(r8) x_old, x_tmp, x_ratio

      
#ifdef CRM_NZ
      ! crm_nz is used to disable warnings above the CRM with MMF
      integer, parameter :: crm_nz = CRM_NZ   
#else
      ! if not MMF, CRM_NZ is not defined, so set to pver to avoid build error
      integer, parameter :: crm_nz = pver
#endif
      logical use_MMF
      call phys_getopts( use_MMF_out = use_MMF)

! ------------------------------------------------------------------------
!      omsm = 1.-1.e-10          ! used to prevent roundoff errors below zero
      omsm = 1._r8-2*epsilon(1._r8) ! used to prevent roundoff errors below zero
      precmin =  0.1_r8/8.64e4_r8      ! set critical value to 0.1 mm/day in kg/m2/s

      adjfac = deltat/(max(deltat,cmftau)) ! adjustment factor from hack scheme

      ! assume 4 m/s fall speed currently (should be improved)
!      vfall = 4.
	
      ! default (if other sol_facts aren't in call, set all to required sol_fact
!     sol_facti = sol_fact
!     sol_factb = sol_fact
!     sol_factii = sol_fact
!     sol_factbi = sol_fact

!     if ( present(sol_facti_in) )  sol_facti = sol_facti_in
!     if ( present(sol_factii_in) )  sol_factii = sol_factii_in
!     if ( present(sol_factbi_in) )  sol_factbi = sol_factbi_in
!     sol_facti = sol_facti_in
!     sol_factii = sol_factii_in
!     sol_factbi = sol_factbi_in

!     sol_factic  = sol_facti
!     sol_factiic = sol_factii
!     if ( present(sol_factic_in ) )  sol_factic  = sol_factic_in
!     if ( present(sol_factiic_in) )  sol_factiic = sol_factiic_in
!     sol_factic  = sol_factic_in
!     sol_factiic = sol_factiic_in

      ! this section of code is for highly soluble aerosols,
      ! the assumption is that within the cloud that
      ! all the tracer is in the cloud water
      !
      ! for both convective and stratiform clouds, 
      ! the fraction of cloud water converted to precip defines
      ! the amount of tracer which is pulled out.

      do i = 1,pcols
         precabs(i) = 0
         precabc(i) = 0
         scavab(i) = 0
         scavabc(i) = 0
         tracab(i) = 0
         cldmabs(i) = 0
         cldmabc(i) = 0

         precabs_base(i) = 0.0_r8
         precabc_base(i) = 0.0_r8
         precnums_base(i) = 0.0_r8
         precnumc_base(i) = 0.0_r8
         resusp_c_sv(i) = 0.0_r8
         resusp_s_sv(i) = 0.0_r8

       ! Jan.16. Sungsu 
       ! I added below to compute vertically projected cumulus and stratus fractions from the top to the
       ! current model layer by assuming a simple independent maximum overlapping assumption for 
       ! each cloud.
       ! cldovr_cu(i) = 0._r8
       ! cldovr_st(i) = 0._r8
       ! End by Sungsu

      end do

main_k_loop: &
      do k = 1,pver
main_i_loop: &
         do i = 1,ncol
            tc     = t(i,k) - tmelt
            weight = max(0._r8,min(-tc*0.05_r8,1.0_r8)) ! fraction of condensate that is ice
            weight = 0._r8                                 ! assume no ice

            pdog = pdel(i,k)/gravit

            ! ****************** Evaporation **************************
            ! calculate the fraction of strat precip from above 
            !                 which evaporates within this layer
            fracev(i) = evaps(i,k)*pdel(i,k)/gravit &
                     /max(1.e-12_r8,precabs(i))

            ! trap to ensure reasonable ratio bounds
            fracev(i) = max(0._r8,min(1._r8,fracev(i)))

! Sungsu : Same as above but convective precipitation part
            fracev_cu(i) = evapc(i,k)*pdel(i,k)/gravit/max(1.e-12_r8,precabc(i))
            fracev_cu(i) = max(0._r8,min(1._r8,fracev_cu(i)))
! Sungsu

            if (mam_prevap_resusp_optcc == 0) then
               fracev(i) = 0.0_r8
               fracev_cu(i) = 0.0_r8
            endif

            ! ****************** Convection ***************************
            ! now do the convective scavenging

            ! set odds proportional to fraction of the grid box that is swept by the 
            ! precipitation =precabc/rhoh20*(area of sphere projected on plane
            !                                /volume of sphere)*deltat
            ! assume the radius of a raindrop is 1 e-3 m from Rogers and Yau,
            ! unless the fraction of the area that is cloud is less than odds, in which
            ! case use the cloud fraction (assumes precabs is in kg/m2/s)
            ! is really: precabs*3/4/1000./1e-3*deltat
            ! here I use .1 from Balkanski
            !
            ! use a local rate of convective rain production for incloud scav
            !odds=max(min(1._r8, &
            !     cmfdqr(i,k)*pdel(i,k)/gravit*0.1_r8*deltat),0._r8)
            !++mcb -- change cldc to cldt; change cldt to cldv (9/17/96)
            !            srcs1 =  cldt(i,k)*odds*tracer(i,k)*(1.-weight) &
            ! srcs1 =  cldv(i,k)*odds*tracer(i,k)*(1.-weight) &
            !srcs1 =  cldc(i,k)*odds*tracer(i,k)*(1.-weight) &
            !         /deltat 

            ! fraction of convective cloud water converted to rain
            ! Dec.29.2009 : Sungsu multiplied cldc(i,k) to conicw(i,k) below
            ! fracp = cmfdqr(i,k)*deltat/max(1.e-8_r8,conicw(i,k))
            ! fracp = cmfdqr(i,k)*deltat/max(1.e-8_r8,cldc(i,k)*conicw(i,k))
            ! Sungsu: Below new formula of 'fracp' is necessary since 'conicw' is a LWC/IWC
            !         that has already precipitated out, that is, 'conicw' does not contain
            !         precipitation at all ! 


    tracer_incu = f_act_conv(i,k)*(tracer(i,k)+&
        min(qqcw(i,k),tracer(i,k)*((cldt(i,k)-cldc(i,k))/max(0.01_r8,(1._r8-(cldt(i,k)-cldc(i,k)))))))
   tracer_mean = tracer(i,k)*(1._r8-cldc(i,k)*f_act_conv(i,k))-cldc(i,k)*f_act_conv(i,k)*&
           min(qqcw(i,k),tracer(i,k)*((cldt(i,k)-cldc(i,k))/max(0.01_r8,(1._r8-(cldt(i,k)-cldc(i,k))))))
    tracer_mean = max(0._r8,tracer_mean)

              fracp = cmfdqr(i,k)*deltat/max(1.e-12_r8,cldc(i,k)*conicw(i,k)+(cmfdqr(i,k)+dlf(i,k))*deltat)
            fracp = max(min(1._r8,fracp),0._r8)
            cldmabc(i) = cldvcu(i,k)

        call wetdep_scavenging(                                &
                        is_strat_cloudborne,    1,  &
                        deltat,         fracp,  cldc(i,k),           &
                        precabc(i),         cldmabc(i), scavcoef(i,k),       &
                        sol_factb,      sol_factbi,             &
                        sol_factic(i,k),      sol_factiic,             &
                        tracer_incu,       tracer_mean,    weight, &
                        srcc,            finc                     )

            ! ****************** Stratiform ***********************
            ! now do the stratiform scavenging
       fracp = precs(i,k)*deltat/max(cwat(i,k)+precs(i,k)*deltat,1.e-12_r8)
        fracp = max(0._r8,min(1._r8,fracp))
        cldmabs(i) = cldvst(i,k) ! Stratiform precipitation area at the top interface of current layer

        call wetdep_scavenging(                                &
                        is_strat_cloudborne,    2,  &
                        deltat,         fracp,  cldc(i,k),           &
                        precabs(i),         cldmabs(i), scavcoef(i,k),       &
                        sol_factb,      sol_factbi,             &
                        sol_facti,      sol_factii,             &
                        tracer(i,k),       tracer_mean,    weight, &
                        srcs,            fins                     )



            ! make sure we dont take out more than is there
            ! ratio of amount available to amount removed
            rat(i) = tracer(i,k)/max(deltat*(srcc+srcs),1.e-36_r8)
            if (rat(i).lt.1._r8) then
               srcs = srcs*rat(i)
               srcc = srcc*rat(i)
            endif
            srct(i) = (srcc+srcs)*omsm

            
            ! fraction that is not removed within the cloud
            ! (assumed to be interstitial, and subject to convective transport)
            fracp = deltat*srct(i)/max(cldmabs(i)*tracer(i,k),1.e-36_r8)
            fracis(i,k) = 1._r8 - min_max_bound(0.0_r8, 1.0_r8, fracp)

! tend is all tracer removed by scavenging, plus all re-appearing from evaporation above
! Sungsu added cumulus contribution in the below blocks
! mam_prevap_resusp_optcc values:
!     0 = no resuspension
!     1 = linear resuspension of aerosol mass or number following original mam coding
!     2 = same as 1 but resuspension tendencies are in rc/sscavt rather than combined with bc/sscavt
!     3 = same as 2 but with some added "xxx = max( 0, xxx)" lines
!   130 = non-linear resuspension of aerosol mass   based on scavenged aerosol mass
!   230 = non-linear resuspension of aerosol number based on raindrop number
!   (only 0, 130 and 230 are used in current code)
resusp_block_aa: &
            if ( mam_prevap_resusp_optcc >= 100) then
                ! force non-negative
                precabs_base(i) = max( 0.0_r8, precabs_base(i) )
                precabs(i)  = min_max_bound( 0.0_r8, precabs_base(i), precabs(i) )
                precabc_base(i) = max( 0.0_r8, precabc_base(i) )
                precabc(i)  = min_max_bound( 0.0_r8, precabc_base(i), precabc(i) )
                if ( mam_prevap_resusp_optcc <= 130) then
                   scavab(i) = max( 0.0_r8, scavab(i) )
                   scavabc(i) = max( 0.0_r8, scavabc(i) )
                else
                   precnums_base(i) = max( 0.0_r8, precnums_base(i) )
                   precnumc_base(i) = max( 0.0_r8, precnumc_base(i) )
                endif

                ! for stratiform clouds
                arainx = max( cldvst(i,min(k+1,pver)), 0.01_r8 )    ! non-zero
                ! step 1 - do evaporation and resuspension
                call wetdep_resusp(                             &
                        1,              mam_prevap_resusp_optcc,& ! in
                        pdel(i,k),      evaps(i,k),             & ! in
                        precabs(i),     precabs_base(i),        & ! in
                        scavab(i),      precnums_base(i),       & ! in
                        precabx_tmp,    precabx_base_tmp,       & ! out
                        scavabx_tmp,    precnumx_base_tmp,      & ! out
                        resusp_s                                ) ! out
                ! step 2 - do precip production and scavenging
                call wetdep_prevap(                             &
                        1,              mam_prevap_resusp_optcc,& ! in
                        pdel(i,k),      precs(i,k),             & ! in
                        srcs,           arainx,                 & ! in
                        precabx_tmp,    precabx_base_tmp,       & ! in
                        scavabx_tmp,    precnumx_base_tmp,      & ! in
                        precabs(i),     precabs_base(i),        & ! out
                        scavab(i),      precnums_base(i)        ) ! out

                ! for convective clouds
                arainx = max( cldvcu(i,min(k+1,pver)), 0.01_r8)     ! non-zero
                call wetdep_resusp(                             &
                        2,              mam_prevap_resusp_optcc,& ! in
                        pdel(i,k),      evapc(i,k),             & ! in
                        precabc(i),     precabc_base(i),        & ! in
                        scavabc(i),     precnumc_base(i),       & ! in
                        precabx_tmp,    precabx_base_tmp,       & ! out
                        scavabx_tmp,    precnumx_base_tmp,      & ! out
                        resusp_c                                ) ! out
                ! step 2 - do precip production and scavenging
                call wetdep_prevap(                             &
                        2,              mam_prevap_resusp_optcc,& ! in
                        pdel(i,k),      cmfdqr(i,k),            & ! in
                        srcc,           arainx,                 & ! in
                        precabx_tmp,    precabx_base_tmp,       & ! in
                        scavabx_tmp,    precnumx_base_tmp,      & ! in
                        precabc(i),     precabc_base(i),        & ! out
                        scavabc(i),     precnumc_base(i)        ) ! out

            else resusp_block_aa ! mam_prevap_resusp_optcc = 0, no resuspension

               resusp_c = fracev_cu(i)*scavabc(i)
               resusp_s = fracev(i)*scavab(i)

            end if resusp_block_aa

! --------------------------------------------------------------

           call update_scavenging( i,            k,             & ! in
                        mam_prevap_resusp_optcc,pdel,           & ! in
                        omsm,   srcc,   srcs,   srct,           & ! in
                        fins,   finc,   fracev, fracev_cu,      & ! in
                        resusp_c,       resusp_s,               & ! in
                        precs,  evaps,  cmfdqr, evapc,          & ! in
                        scavt,  iscavt, icscavt,isscavt,        & ! inout
                        bcscavt,bsscavt,rcscavt,rsscavt,        & ! inout
                        scavab, scavabc,precabc,precabs         ) ! inout

         end do main_i_loop ! End of i = 1, ncol

      end do main_k_loop ! End of k = 1, pver

   end subroutine wetdepa_v2

!==============================================================================
   subroutine wetdep_scavenging(                                &
                        is_strat_cloudborne,    is_conv_strat,  &
                        deltat,         fracp,  cldc_ik,           &
                        precab_i,         cldmab_i, scavcoef_ik,       &
                        sol_factb,      sol_factbi,             &
                        sol_facti,      sol_factii,             &
                        tracer_1,       tracer_mean,    weight, &
                        src,            fin                     )
! ------------------------------------------------------------------------------
! do the stratiform scavenging
! ------------------------------------------------------------------------------

   logical, intent(in) :: is_strat_cloudborne
   integer, intent(in) :: is_conv_strat ! 1: convective. 2: stratiform

   real(r8), intent(in) :: deltat       ! timestep
   real(r8), intent(in) :: fracp        ! fraction of cloud water converted to precip
   real(r8), intent(in) :: cldc_ik     ! convective cloud fraction
   real(r8), intent(in) :: precab_i       ! precip from above of the layer
   real(r8), intent(in) :: cldmab_i        ! precipitation area at the top interface
   real(r8), intent(in) :: scavcoef_ik ! Dana and Hales coefficient [1/mm]
   real(r8), intent(in) :: sol_factb   ! solubility factor (frac of aerosol scavenged below cloud)
   real(r8), intent(in) :: sol_factbi  ! solubility factor (frac of aerosol scavenged below cloud by ice)
   real(r8), intent(in) :: sol_facti   ! solubility factor (frac of aerosol scavenged in cloud)
   real(r8), intent(in) :: sol_factii  ! solubility factor (frac of aerosol scavenged in cloud by ice)
   real(r8), intent(in) :: tracer_1
   real(r8), intent(in) :: tracer_mean
   real(r8), intent(in) :: weight      ! fraction of condensate which is ice
   real(r8), intent(out) :: src
   real(r8), intent(out) :: fin

   real(r8), parameter :: small_value = 1.e-12_r8
   real(r8) :: src1      ! incloud scavenging tendency
   real(r8) :: src2      ! below-cloud scavenging tendency
   real(r8) :: odds       ! limit on removal rate (proportional to prec)

   ! calculate limitation of removal rate using Dana and Hales coefficient
   odds = precab_i/max(cldmab_i,1.e-5_r8)*scavcoef_ik*deltat
   odds = max(min(1._r8,odds),0._r8)

   if ( is_strat_cloudborne ) then
      if (is_conv_strat == 1) then
         ! convective cloud does not affect strat-cloudborne aerosol
         src1 = 0._r8
      else
         ! strat in-cloud removal only affects strat-cloudborne aerosol
         ! in-cloud scavenging:
         src1 = sol_facti *fracp*tracer_1/deltat*(1._r8-weight) & ! Liquid
              + sol_factii*fracp*tracer_1/deltat*(weight)          !Ice
      endif
      ! no below-cloud scavenging
      src2 = 0._r8
  
   else
      if (is_conv_strat == 1) then      ! convective
         src1 = sol_facti*cldc_ik*fracp*tracer_1*(1._r8-weight)/deltat &  ! Liquid
              + sol_factii*cldc_ik*fracp*tracer_1*(weight)/deltat              ! Ice
      else      ! stratiform
         ! strat in-cloud removal only affects strat-cloudborne aerosol
         src1 = 0._r8
      endif
      src2 = sol_factb *cldmab_i*odds*tracer_mean*(1._r8-weight)/deltat & ! Liquid
           + sol_factbi*cldmab_i*odds*tracer_mean*(weight)/deltat ! Ice

   endif

   src = src1 + src2             ! total stratiform or convective scavenging
   fin=src1/(src + 1.e-36_r8)    ! fraction taken by incloud processes

   end subroutine wetdep_scavenging

!==============================================================================
   subroutine wetdep_resusp(                         &
                        jstrcnv,        mam_prevap_resusp_optcc,& ! in
                        pdel_ik,        evapx,                  & ! in
                        precabx_old,    precabx_base_old,       & ! in
                        scavabx_old,    precnumx_base_old,      & ! in
                        precabx_new,    precabx_base_new,       & ! out
                        scavabx_new,    precnumx_base_new,      & ! out
                        resusp_x                                ) ! out

! ------------------------------------------------------------------------------
! do precip production, resuspension and scavenging
! ------------------------------------------------------------------------------
   integer, intent(in) :: jstrcnv       ! options for stratiform (1) or convective (2) clouds
                                        ! raindrop size distribution is
                                        ! different for different cloud:
                                        ! 1: assume marshall-palmer distribution
                                        ! 2: assume log-normal distribution
   integer, intent(in) :: mam_prevap_resusp_optcc       ! suspension options
   real(r8),intent(in) :: pdel_ik       ! pressure thikness at current column and level [Pa]
   real(r8),intent(in) :: evapx
   real(r8),intent(in) :: precabx_base_old
   real(r8),intent(in) :: precabx_old
   real(r8),intent(in) :: scavabx_old
   real(r8),intent(in) :: precnumx_base_old
   real(r8),intent(out) :: precabx_base_new
   real(r8),intent(out) :: precabx_new
   real(r8),intent(out) :: scavabx_new
   real(r8),intent(out) :: precnumx_base_new
   real(r8),intent(out) :: resusp_x


   ! local variables
   real(r8), parameter :: prec_smallaa = 1.0e-30_r8  ! 1e-30 kg/m2/s (or mm/s) = 3.2e-23 mm/yr
   real(r8)            :: tmpa ! temporary working variable

   ! initiate *_new in case they are not calculated
   scavabx_new = scavabx_old
   precnumx_base_new = precnumx_base_old
   precabx_base_new = precabx_base_old

   tmpa = max( 0.0_r8, evapx*pdel_ik/gravit )
   precabx_new = min_max_bound(0.0_r8,precabx_base_new, precabx_old-tmpa)

   if (precabx_new < prec_smallaa) then
      ! precip rate is essentially zero so do complete resuspension
      call wetdep_resusp_noprecip(                         &
                   jstrcnv,        mam_prevap_resusp_optcc,& ! in
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
                     jstrcnv,        mam_prevap_resusp_optcc,& ! in
                     precabx_old,    precabx_base_old,       & ! in
                     scavabx_old,    precnumx_base_old,      & ! in
                     precabx_new,                            & ! in
                     scavabx_new,    resusp_x                ) ! out

   endif

   end subroutine wetdep_resusp


!==============================================================================
   subroutine wetdep_resusp_noprecip(                           &
                        jstrcnv,        mam_prevap_resusp_optcc,& ! in
                        precabx_old,    precabx_base_old,       & ! in
                        scavabx_old,    precnumx_base_old,      & ! in
                        precabx_new,    precabx_base_new,       & ! out
                        scavabx_new,    resusp_x                ) ! out

! ------------------------------------------------------------------------------
! do complete resuspension when precipitation rate is zero
! ------------------------------------------------------------------------------
   integer, intent(in) :: jstrcnv       ! options for stratiform (1) or convective (2) clouds
                                        ! raindrop size distribution is
                                        ! different for different cloud:
                                        ! 1: assume marshall-palmer distribution
                                        ! 2: assume log-normal distribution
   integer, intent(in) :: mam_prevap_resusp_optcc       ! suspension options
   real(r8),intent(in) :: precabx_base_old
   real(r8),intent(in) :: precabx_old
   real(r8),intent(in) :: scavabx_old
   real(r8),intent(in) :: precnumx_base_old
   real(r8),intent(out) :: precabx_base_new
   real(r8),intent(out) :: precabx_new
   real(r8),intent(out) :: scavabx_new
   real(r8),intent(out) :: resusp_x


   ! local variables
   real(r8), parameter :: prec_smallaa = 1.0e-30_r8  ! 1e-30 kg/m2/s (or mm/s) = 3.2e-23 mm/yr
   real(r8)     :: u_old,u_new  ! fraction of precabx and precabx_base
   real(r8)     :: x_old,x_new  ! fraction after calling function *_resusp_vs_fprec_evap_mpln


   if ( mam_prevap_resusp_optcc <= 130) then
        ! linear resuspension based on scavenged aerosol mass or number
        scavabx_new = 0.0_r8
        resusp_x = scavabx_old
   else
        if (precabx_base_old < prec_smallaa) then
            resusp_x = 0.0_r8
        else
            ! non-linear resuspension of aerosol number based on raindrop number
            u_old = min_max_bound(0.0_r8, 1.0_r8, precabx_old/precabx_base_old)
            x_old = 1.0_r8 - fprecn_resusp_vs_fprec_evap_mpln(1.0_r8-u_old, jstrcnv )
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
                        jstrcnv,        mam_prevap_resusp_optcc,& ! in
                        precabx_old,    precabx_base_old,       & ! in
                        scavabx_old,    precnumx_base_old,      & ! in
                        precabx_new,                            & ! in
                        scavabx_new,    resusp_x                ) ! out

! ------------------------------------------------------------------------------
! do nonlinear resuspension of aerosol mass or number
! ------------------------------------------------------------------------------
   integer, intent(in) :: jstrcnv       ! options for stratiform (1) or convective (2) clouds
                                        ! raindrop size distribution is
                                        ! different for different cloud:
                                        ! 1: assume marshall-palmer distribution
                                        ! 2: assume log-normal distribution
   integer, intent(in) :: mam_prevap_resusp_optcc       ! suspension options
   real(r8),intent(in) :: precabx_base_old
   real(r8),intent(in) :: precabx_old
   real(r8),intent(in) :: scavabx_old
   real(r8),intent(in) :: precnumx_base_old
   real(r8),intent(in) :: precabx_new
   real(r8),intent(out) :: scavabx_new
   real(r8),intent(out) :: resusp_x


   ! local variables
   real(r8), parameter :: x_smallaa = 1.0e-30_r8     ! small value for x_old
   real(r8)     :: u_old,u_new  ! fraction of precabx and precabx_base
   real(r8)     :: x_old,x_new  ! fraction after calling function *_resusp_vs_fprec_evap_mpln
   real(r8)     :: x_ratio      ! fraction of x_tmp/x_old


   u_old = min_max_bound(0.0_r8, 1.0_r8, precabx_old/precabx_base_old)
   if ( mam_prevap_resusp_optcc <= 130) then
        ! non-linear resuspension of aerosol mass
        x_old = 1.0_r8 - faer_resusp_vs_fprec_evap_mpln( 1.0_r8-u_old, jstrcnv)
   else
        ! non-linear resuspension of aerosol number based on raindrop number
        x_old = 1.0_r8 - fprecn_resusp_vs_fprec_evap_mpln(1.0_r8-u_old, jstrcnv)
   endif
   x_old = min_max_bound(0.0_r8, 1.0_r8, x_old)

   if (x_old < x_smallaa) then
        x_new = 0.0_r8
        x_ratio = 0.0_r8
   else
        u_new = min_max_bound(0.0_r8, 1.0_r8, precabx_new/precabx_base_old)
        u_new = min( u_new, u_old )
        if ( mam_prevap_resusp_optcc <= 130) then
            ! non-linear resuspension of aerosol mass
            x_new = 1.0_r8 - faer_resusp_vs_fprec_evap_mpln(1.0_r8-u_new, jstrcnv )
        else
            ! non-linear resuspension of aerosol number based on raindrop number
            x_new = 1.0_r8 - fprecn_resusp_vs_fprec_evap_mpln(1.0_r8-u_new, jstrcnv )
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
                        jstrcnv,        mam_prevap_resusp_optcc,& ! in
                        pdel_ik,        pprdx,  srcx,   arainx, & ! in
                        precabx_old,    precabx_base_old,       & ! in
                        scavabx_old,    precnumx_base_old,      & ! in
                        precabx_new,    precabx_base_new,       & ! out
                        scavabx_new,    precnumx_base_new       ) ! out

! ------------------------------------------------------------------------------
! do precip production and scavenging
! ------------------------------------------------------------------------------
   integer, intent(in) :: jstrcnv       ! options for stratiform (1) or convective (2) clouds
                                        ! raindrop size distribution is
                                        ! different for different cloud:
                                        ! 1: assume marshall-palmer distribution
                                        ! 2: assume log-normal distribution
   integer, intent(in) :: mam_prevap_resusp_optcc       ! suspension options
   real(r8),intent(in) :: pdel_ik       ! pressure thikness at current column and level [Pa]
   real(r8),intent(in) :: pprdx
   real(r8),intent(in) :: srcx
   real(r8),intent(in) :: arainx ! precipitation and cloudy volume,at the top interface of current layer [fraction]
   real(r8),intent(in) :: precabx_base_old
   real(r8),intent(in) :: precabx_old
   real(r8),intent(in) :: scavabx_old
   real(r8),intent(in) :: precnumx_base_old
   real(r8),intent(out) :: precabx_base_new
   real(r8),intent(out) :: precabx_new
   real(r8),intent(out) :: scavabx_new
   real(r8),intent(out) :: precnumx_base_new

   ! local variables
   real(r8), parameter :: prec_smallaa = 1.0e-30_r8  ! 1e-30 kg/m2/s (or mm/s) = 3.2e-23 mm/yr
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
       if (precabx_base_new < prec_smallaa) then
          precnumx_base_new = 0.0_r8
       elseif (precabx_base_new > precabx_base_old) then
          ! note - calc rainshaft number flux from rainshaft water flux,
          ! then multiply by rainshaft area to get grid-average number flux
          tmpa = arainx * flux_precnum_vs_flux_prec_mpln((precabx_base_new/arainx), jstrcnv )
          precnumx_base_new = max( 0.0_r8, tmpa )
       else
          precnumx_base_new = precnumx_base_old
       endif
   endif

   end subroutine wetdep_prevap

!==============================================================================
   subroutine update_scavenging( icol,          kk,             & ! in
                        mam_prevap_resusp_optcc,pdel,           & ! in
                        omsm,   srcc,   srcs,   srct,           & ! in
                        fins,   finc,   fracev, fracev_cu,      & ! in 
                        resusp_c,       resusp_s,               & ! in
                        precs,  evaps,  cmfdqr, evapc,          & ! in
                        scavt,  iscavt, icscavt,isscavt,        & ! inout
                        bcscavt,bsscavt,rcscavt,rsscavt,        & ! inout
                        scavab, scavabc,precabc,precabs         ) ! inout
! ------------------------------------------------------------------------------
! update scavenging variables
! ------------------------------------------------------------------------------
   ! input variables
   integer, intent(in) :: icol                ! column index
   integer, intent(in) :: kk                  ! level index
   integer, intent(in) :: mam_prevap_resusp_optcc       ! suspension options 
   real(r8),intent(in) :: pdel(pcols,pver)    ! pressure thikness [Pa]
   real(r8),intent(in) :: omsm                ! 1 - (a small number), used to prevent roundoff errors below zero
   real(r8),intent(in) :: srcc                ! tend for convective rain [kg/kg/s]
   real(r8),intent(in) :: srcs                ! tend for stratiform rain [kg/kg/s]
   real(r8),intent(in) :: srct(pcols)         ! total tendency for conv+strat rain [kg/kg/s]
   real(r8),intent(in) :: fins                ! fraction of rem. rate by strat rain [fraction]
   real(r8),intent(in) :: finc                ! fraction of rem. rate by conv. rain [fraction]
   real(r8),intent(in) :: fracev(pcols)       ! fraction of precip from above that is evaporating [fraction]
   real(r8),intent(in) :: fracev_cu(pcols)    ! Fraction of convective precip from above that is evaporating [fraction]
   real(r8),intent(in) :: resusp_c            ! aerosol mass re-suspension in a particular layer from convective rain [kg/m2/s]
   real(r8),intent(in) :: resusp_s            ! aerosol mass re-suspension in a particular layer from stratiform rain [kg/m2/s]
   real(r8),intent(in) :: precs(pcols,pver)   ! rate of production of stratiform precip [kg/kg/s]
   real(r8),intent(in) :: evaps(pcols,pver)   ! rate of evaporation of precip [kg/kg/s]
   real(r8),intent(in) :: cmfdqr(pcols,pver)  ! rate of production of convective precip [kg/kg/s]
   real(r8),intent(in) :: evapc(pcols,pver)   ! Evaporation rate of convective precipitation [kg/kg/s]
   ! output variables
   real(r8), intent(inout) :: scavt(pcols,pver)    ! scavenging tend [kg/kg/s]
   real(r8), intent(inout) :: iscavt(pcols,pver)   ! incloud scavenging tends [kg/kg/s]
   real(r8), intent(inout) :: icscavt(pcols,pver)  ! incloud, convective [kg/kg/s]
   real(r8), intent(inout) :: isscavt(pcols,pver)  ! incloud, stratiform [kg/kg/s]
   real(r8), intent(inout) :: bcscavt(pcols,pver)  ! below cloud, convective [kg/kg/s]
   real(r8), intent(inout) :: bsscavt(pcols,pver)  ! below cloud, stratiform [kg/kg/s]
   real(r8), intent(inout) :: rcscavt(pcols,pver)  ! resuspension, convective [kg/kg/s]
   real(r8), intent(inout) :: rsscavt(pcols,pver)  ! resuspension, stratiform [kg/kg/s]
   real(r8), intent(inout) :: scavab(pcols)        ! scavenged tracer flux from above [kg/m2/s]
   real(r8), intent(inout) :: scavabc(pcols)       ! scavenged tracer flux from above [kg/m2/s]
   real(r8), intent(inout) :: precabc(pcols)       ! conv precip from above [kg/m2/s]
   real(r8), intent(inout) :: precabs(pcols)       ! strat precip from above [kg/m2/s]

   if ( mam_prevap_resusp_optcc == 0) then
       scavt(icol,kk) = -srct(icol) + &
            (fracev(icol)*scavab(icol)+fracev_cu(icol)*scavabc(icol))*gravit/pdel(icol,kk)
   else
       scavt(icol,kk) = -srct(icol) + (resusp_s+resusp_c)*gravit/pdel(icol,kk)
   endif

   iscavt(icol,kk) = -(srcc*finc + srcs*fins)*omsm
   icscavt(icol,kk) = -(srcc*finc) * omsm
   isscavt(icol,kk) = -(srcs*fins) * omsm

   if (mam_prevap_resusp_optcc == 0) then
          bcscavt(icol,kk) = -(srcc * (1-finc)) * omsm +  &
              fracev_cu(icol)*scavabc(icol)*gravit/pdel(icol,kk)

          bsscavt(icol,kk) = -(srcs * (1-fins)) * omsm +  &
              fracev(icol)*scavab(icol)*gravit/pdel(icol,kk)

          rcscavt(icol,kk) = 0.0
          rsscavt(icol,kk) = 0.0
   else ! here mam_prevap_resusp_optcc == 130, 210, 230
          bcscavt(icol,kk) = -(srcc * (1-finc)) * omsm
          rcscavt(icol,kk) = resusp_c*gravit/pdel(icol,kk)
          bsscavt(icol,kk) = -(srcs * (1-fins)) * omsm
          rsscavt(icol,kk) = resusp_s*gravit/pdel(icol,kk)
   endif

   ! now keep track of scavenged mass and precip
   if (mam_prevap_resusp_optcc == 0) then
       scavab(icol)  = scavab(icol) *(1-fracev   (icol)) + srcs*pdel(icol,kk)/gravit
       scavabc(icol) = scavabc(icol)*(1-fracev_cu(icol)) + srcc*pdel(icol,kk)/gravit
       precabs(icol) = precabs(icol) + (precs (icol,kk) - evaps(icol,kk))*pdel(icol,kk)/gravit
       precabc(icol) = precabc(icol) + (cmfdqr(icol,kk) - evapc(icol,kk))*pdel(icol,kk)/gravit
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

      if (flux_prec >= 1.0e-36_r8) then
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

      x_var = max( 0.0_r8, min( 1.0_r8, fprec_evap ) )
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

      x_var = max( 0.0_r8, min( 1.0_r8, fprec_evap ) )
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
