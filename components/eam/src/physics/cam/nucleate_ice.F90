module nucleate_ice

!-------------------------------------------------------------------------------
! Purpose:
!  A parameterization of ice nucleation.
!
!  *** This module is intended to be a "portable" code layer.  Ideally it should
!  *** not contain any use association of modules that belong to the model framework.
!
!
! Method:
!  The current method is based on Liu & Penner (2005) & Liu et al. (2007)
!  It related the ice nucleation with the aerosol number, temperature and the
!  updraft velocity. It includes homogeneous freezing of sulfate & immersion
!  freezing on mineral dust (soot disabled) in cirrus clouds, and 
!  Meyers et al. (1992) deposition nucleation in mixed-phase clouds
!
!  The effect of preexisting ice crystals on ice nucleation in cirrus clouds is included, 
!  and also consider the sub-grid variability of temperature in cirrus clouds,
!  following X. Shi et al. ACP (2014).
!
!  Ice nucleation in mixed-phase clouds now uses classical nucleation theory (CNT),
!  follows Y. Wang et al. ACP (2014), Hoose et al. (2010).
!
! Authors:
!  Xiaohong Liu, 01/2005, modifications by A. Gettelman 2009-2010
!  Xiangjun Shi & Xiaohong Liu, 01/2014.
!
!  With help from C. C. Chen and B. Eaton (2014)
!-------------------------------------------------------------------------------

use wv_saturation,  only: svp_water, svp_ice
use cam_logfile,    only: iulog
use phys_control,   only: cam_chempkg_is 

implicit none
private
save

integer, parameter :: r8 = selected_real_kind(12)

public :: nucleati_init, nucleati

logical  :: use_preexisting_ice
logical  :: use_hetfrz_classnuc
logical  :: use_nie_nucleate
logical  :: use_dem_nucleate

real(r8) :: pi
real(r8) :: mincld

! Subgrid scale factor on relative humidity (dimensionless)
real(r8) :: subgrid

real(r8), parameter :: rhoice = 0.5e3_r8   ! kg/m3, Wpice is not sensitive to rhoice 

real(r8) :: ci

!===============================================================================
contains
!===============================================================================

subroutine nucleati_init( &
   use_preexisting_ice_in, use_hetfrz_classnuc_in, &
   use_nie_nucleate_in, use_dem_nucleate_in, &
   iulog_in, pi_in, mincld_in, subgrid_in)

   logical,  intent(in) :: use_preexisting_ice_in
   logical,  intent(in) :: use_hetfrz_classnuc_in
   logical,  intent(in) :: use_nie_nucleate_in
   logical,  intent(in) :: use_dem_nucleate_in
   integer,  intent(in) :: iulog_in
   real(r8), intent(in) :: pi_in
   real(r8), intent(in) :: mincld_in
   real(r8), intent(in) :: subgrid_in

   use_preexisting_ice = use_preexisting_ice_in
   use_hetfrz_classnuc = use_hetfrz_classnuc_in
   use_nie_nucleate    = use_nie_nucleate_in
   use_dem_nucleate    = use_dem_nucleate_in
   iulog               = iulog_in
   pi                  = pi_in
   mincld              = mincld_in
   subgrid             = subgrid_in

   ci = rhoice*pi/6._r8

end subroutine nucleati_init

!===============================================================================

subroutine nucleati(  &
   wbar, tair, pmid, relhum, cldn,      &
   rhoair, so4_num, dst3_num,  &
   nuci, onihf, oniimm, onidep, onimey)

   !---------------------------------------------------------------
   ! Purpose:
   !  The parameterization of ice nucleation.
   !
   ! Method: The current method is based on Liu & Penner (2005)
   !  It related the ice nucleation with the aerosol number, temperature and the
   !  updraft velocity. It includes homogeneous freezing of sulfate, immersion
   !  freezing of soot, and Meyers et al. (1992) deposition nucleation
   !
   ! Authors: Xiaohong Liu, 01/2005, modifications by A. Gettelman 2009-2010
   !----------------------------------------------------------------

   ! Input Arguments
   real(r8), intent(in) :: wbar        ! grid cell mean vertical velocity [m/s]
   real(r8), intent(in) :: tair        ! temperature [K]
   real(r8), intent(in) :: pmid        ! pressure at layer midpoints [pa]
   real(r8), intent(in) :: relhum      ! relative humidity with respective to liquid [unitless]
   real(r8), intent(in) :: cldn        ! new value of cloud fraction    [fraction]
   real(r8), intent(in) :: rhoair      ! air density [kg/m3]
   real(r8), intent(in) :: so4_num     ! so4 aerosol number [#/cm^3]
   real(r8), intent(in) :: dst3_num     ! dust aerosol number [#/cm^3]

   ! Output Arguments
   real(r8), intent(out) :: nuci       ! ice number nucleated [#/kg]
   real(r8), intent(out) :: onihf      ! nucleated number from homogeneous freezing of so4 [#/kg]
   real(r8), intent(out) :: oniimm     ! nucleated number from immersion freezing [#/kg]
   real(r8), intent(out) :: onidep     ! nucleated number from deposition nucleation [#/kg]
   real(r8), intent(out) :: onimey     ! nucleated number from deposition nucleation  (meyers: mixed phase) [#/kg]

   ! Local workspace
   real(r8) :: nihf                      ! nucleated number from homogeneous freezing of so4
   real(r8) :: niimm                     ! nucleated number from immersion freezing
   real(r8) :: nidep                     ! nucleated number from deposition nucleation
   real(r8) :: nimey                     ! nucleated number from deposition nucleation (meyers)
   real(r8) :: n1, ni                    ! nucleated number
   real(r8) :: tc, regm                  ! work variable  
   real(r8) :: wbar1, wbar2

   !-------------------------------------------------------------------------------

   ! temp variables that depend on use_preexisting_ice
   wbar1 = wbar
   wbar2 = wbar

   ni = 0._r8
   tc = tair - 273.15_r8

   ! initialize
   niimm = 0._r8
   nidep = 0._r8
   nihf  = 0._r8


   if(so4_num >= 1.0e-10_r8 .and. dst3_num >= 1.0e-10_r8 .and. cldn > 0._r8) then

      if((tc <= -35.0_r8) .and. ((relhum*svp_water(tair)/svp_ice(tair)*subgrid) >= 1.2_r8)) then ! use higher RHi threshold

            call calculate_regm_nucleati(wbar1, dst3_num, regm)

            ! heterogeneous nucleation only
            if (tc > regm) then

               if(tc < -40._r8 .and. wbar1 > 1._r8) then ! exclude T<-40 & W>1m/s from hetero. nucleation

                  call hf(tc,wbar1,relhum,so4_num,nihf)
                  niimm=0._r8
                  nidep=0._r8   
                  n1=nihf
                  
               else

                  call hetero(tc,wbar2,dst3_num,niimm,nidep)
                  nihf=0._r8
                  n1=niimm+nidep

               endif

            ! homogeneous nucleation only
            else if (tc < regm-5._r8) then

               call hf(tc,wbar1,relhum,so4_num,nihf)
               niimm=0._r8
               nidep=0._r8
               n1=nihf
              
            ! transition between homogeneous and heterogeneous: interpolate in-between
            else

               if (tc < -40._r8 .and. wbar1 > 1._r8) then ! exclude T<-40 & W>1m/s from hetero. nucleation

                  call hf(tc, wbar1, relhum, so4_num, nihf)
                  niimm = 0._r8
                  nidep = 0._r8
                  n1 = nihf
                  
               else

                  call hf(regm-5._r8,wbar1,relhum,so4_num,nihf)
                  call hetero(regm,wbar2,dst3_num,niimm,nidep)

                  if (nihf <= (niimm+nidep)) then
                     n1 = nihf
                  else
                     n1=(niimm+nidep)*((niimm+nidep)/nihf)**((tc-regm)/5._r8)
                  endif

               end if
            end if

            ni = n1

         end if
      end if


   ! deposition/condensation nucleation in mixed clouds (-37<T<0C) (Meyers, 1992)
   ! this part is executed but is always replaced by 0, because CNT scheme takes over
   ! the calculation. use_hetfrz_classnuc is always true.

   nimey = 0._r8

   nuci=ni+nimey
   if(nuci > 9999._r8 .or. nuci < 0._r8) then
      nuci=0._r8
   endif
   
   nuci   = nuci*1.e+6_r8/rhoair    ! change unit from #/cm3 to #/kg
   onimey = nimey*1.e+6_r8/rhoair
   onidep = nidep*1.e+6_r8/rhoair
   oniimm = niimm*1.e+6_r8/rhoair
   onihf  = nihf*1.e+6_r8/rhoair

end subroutine nucleati

!===============================================================================

subroutine hetero(Temperature,w_vlc,Ns,Nis,Nid)

!-------------------------------------------------------------------------------
! Calculate number of ice crystals by heterogenous freezing (Nis) based on
! Eq. 4.7 in Liu & Penner (2005), Meteorol. Z.
!-------------------------------------------------------------------------------

    real(r8), intent(in)  :: Temperature     ! temperature [C]
    real(r8), intent(in)  :: w_vlc           ! vertical velocity [m/s]
    real(r8), intent(in)  :: Ns              ! aerosol concentrations [#/cm^3]
    real(r8), intent(out) :: Nis             ! ice number concentrations [#/cm^3]
    real(r8), intent(out) :: Nid             ! ice number concentrations [#/cm^3]  

!---------------------------------------------------------------------
! parameters
!---------------------------------------------------------------------

    real(r8), parameter :: A11 =  0.0263_r8
    real(r8), parameter :: A12 = -0.0185_r8
    real(r8), parameter :: A21 =  2.758_r8
    real(r8), parameter :: A22 =  1.3221_r8
    real(r8), parameter :: B11 = -0.008_r8
    real(r8), parameter :: B12 = -0.0468_r8
    real(r8), parameter :: B21 = -0.2667_r8
    real(r8), parameter :: B22 = -1.4588_r8

!---------------------------------------------------------------------
! local variables
!---------------------------------------------------------------------
   
   real(r8) lnNs, lnw, B_coef, C_coef

   lnNs = log(Ns)
   lnw  = log(w_vlc)

!  ice from immersion nucleation (cm^-3)

   B_coef = (A11 + B11*lnNs) * lnw + (A12 + B12*lnNs) 
   C_coef = A21 + B21*lnNs

   Nis = exp(A22) * Ns**B22 * exp(B_coef*Temperature) * w_vlc**C_coef
   Nis = min(Nis,Ns)

   Nid = 0.0_r8    ! don't include deposition nucleation for cirrus clouds when T<-37C

end subroutine hetero

!===============================================================================

subroutine hf(Temperature, w_vlc, RH, Na, Ni)

!-------------------------------------------------------------------------------
! Calculate number of ice crystals by homogeneous freezing (Ni) based on
! Liu & Penner (2005), Meteorol. Z.
!-------------------------------------------------------------------------------

      real(r8), intent(in)  :: Temperature     ! temperature [C]
      real(r8), intent(in)  :: w_vlc           ! vertical velocity [m/s]
      real(r8), intent(in)  :: RH              ! unitless relative humidity
      real(r8), intent(in)  :: Na              ! aerosol number concentrations [#/cm^3]
      real(r8), intent(out) :: Ni              ! ice number concentrations [#/cm^3]

!---------------------------------------------------------------------
! parameters
!---------------------------------------------------------------------

      real(r8), parameter :: A1_fast  =  0.0231_r8
      real(r8), parameter :: A21_fast = -1.6387_r8  !(T>-64 deg)
      real(r8), parameter :: A22_fast = -6.045_r8   !(T<=-64 deg)
      real(r8), parameter :: B1_fast  = -0.008_r8
      real(r8), parameter :: B21_fast = -0.042_r8   !(T>-64 deg)
      real(r8), parameter :: B22_fast = -0.112_r8   !(T<=-64 deg)
      real(r8), parameter :: C1_fast  =  0.0739_r8
      real(r8), parameter :: C2_fast  =  1.2372_r8

      real(r8), parameter :: A1_slow  = -0.3949_r8
      real(r8), parameter :: A2_slow  =  1.282_r8
      real(r8), parameter :: B1_slow  = -0.0156_r8
      real(r8), parameter :: B2_slow  =  0.0111_r8
      real(r8), parameter :: B3_slow  =  0.0217_r8
      real(r8), parameter :: C1_slow  =  0.120_r8
      real(r8), parameter :: C2_slow  =  2.312_r8

!---------------------------------------------------------------------
! local variables
!---------------------------------------------------------------------
      real(r8) A2_fast, B2_fast, B4_slow
      real(r8) lnw, regm, RHw 

      lnw = log(w_vlc)
      Ni = 0.0_r8

      call calculate_RHw_hf(Temperature, lnw, RHw)

      if((Temperature <= -37.0_r8) .and. ((RH*subgrid) >= RHw)) then

        regm = 6.07_r8*lnw-55.0_r8

        if(Temperature >= regm) then    ! fast-growth regime

          if(Temperature > -64.0_r8) then
            A2_fast=A21_fast
            B2_fast=B21_fast
          else
            A2_fast=A22_fast
            B2_fast=B22_fast
          endif

          call calculate_Ni_hf(A1_fast, B1_fast, C1_fast, A2_fast, B2_fast, C2_fast, Temperature, lnw, Na, Ni)

        else       ! slow-growth regime
 
          B4_slow = B2_slow + B3_slow*lnw   

          call calculate_Ni_hf(A1_slow, B1_slow, C1_slow, A2_slow, B4_slow, C2_slow, Temperature, lnw, Na, Ni)

        endif

      endif

end subroutine hf

!===============================================================================

subroutine calculate_regm_nucleati(w_vlc, Na, regm)

!-------------------------------------------------------------------------------
! Calculate temperature regime for ice nucleation based on
! Eq. 4.5 in Liu & Penner (2005), Meteorol. Z.
!-------------------------------------------------------------------------------

   real(r8), intent(in)  :: w_vlc            ! vertical velocity [m/s]
   real(r8), intent(in)  :: Na               ! aerosol number concentration [#/cm^3]
   real(r8), intent(out) :: regm             ! threshold temperature [C] 

   real(r8) A_coef, B_coef, lnNa

   lnNa   = log(Na)

   A_coef = -1.4938_r8 * lnNa + 12.884_r8
   B_coef = -10.41_r8  * lnNa - 67.69_r8

   regm = A_coef * log(w_vlc) + B_coef
end subroutine calculate_regm_nucleati


subroutine calculate_RHw_hf(Temperature, lnw, RHw)

!-------------------------------------------------------------------------------
! Calculate threshold relative humidity with respective to water (RHw) based on
! Eq. 3.1 in Liu & Penner (2005), Meteorol. Z.
!-------------------------------------------------------------------------------

   real(r8), intent(in)  :: Temperature     ! temperature [C]
   real(r8), intent(in)  :: lnw             ! ln of vertical velocity
   real(r8), intent(out) :: RHw             ! relative humidity threshold 

   real(r8) A_coef, B_coef, C_coef
 
   A_coef = 6.0e-4_r8*lnw + 6.6e-3_r8
   B_coef = 6.0e-2_r8*lnw + 1.052_r8
   C_coef = 1.68_r8  *lnw + 129.35_r8
   
   RHw = (A_coef*Temperature*Temperature + B_coef*Temperature + C_coef)*0.01_r8

end subroutine calculate_RHw_hf


subroutine calculate_Ni_hf(A1, B1, C1, A2, B2, C2, Temperature, lnw, Na, Ni)
   
!-------------------------------------------------------------------------------
! Calculate number of ice crystals (Ni) based on
! Eq. 3.3 in Liu & Penner (2005), Meteorol. Z.
!-------------------------------------------------------------------------------

   real(r8), intent(in)  :: A1, B1, C1     ! Coefficients        
   real(r8), intent(in)  :: A2, B2, C2     ! Coefficients
   real(r8), intent(in)  :: Temperature    ! temperature [C]
   real(r8), intent(in)  :: lnw            ! ln of vertical velocity
   real(r8), intent(in)  :: Na             ! aerosol number concentrations [#/cm^3]
   real(r8), intent(out) :: Ni             ! ice number concentrations [#/cm^3]
   
   real(r8) k1, k2

   k1 = exp(A2 + B2*Temperature + C2*lnw)   
   k2 = A1 + B1*Temperature + C1*lnw

   Ni = k1*Na**(k2)
   Ni = min(Ni,Na)

end subroutine calculate_Ni_hf

end module nucleate_ice

