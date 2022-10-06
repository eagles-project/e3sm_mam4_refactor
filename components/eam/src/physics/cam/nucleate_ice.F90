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

real(r8), parameter :: Shet   = 1.3_r8     ! het freezing threshold
real(r8), parameter :: rhoice = 0.5e3_r8   ! kg/m3, Wpice is not sensitive to rhoice
real(r8), parameter :: minweff= 0.001_r8   ! m/s
real(r8), parameter :: gamma1=1.0_r8 
real(r8), parameter :: gamma2=1.0_r8 
real(r8), parameter :: gamma3=2.0_r8 
real(r8), parameter :: gamma4=6.0_r8 

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
   qc, qi, ni_in, rhoair,               &
   so4_num, dst_num, soot_num,          &
   dst1_sfc_to_num, dst3_sfc_to_num,    &
   nuci, onihf, oniimm, onidep, onimey, &
   wpice, weff, fhom,                   &
   dst1_num,dst2_num,dst3_num,dst4_num, &
   organic_num, clim_modal_aero )

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
   real(r8), intent(in) :: wbar        ! grid cell mean vertical velocity (m/s)
   real(r8), intent(in) :: tair        ! temperature (K)
   real(r8), intent(in) :: pmid        ! pressure at layer midpoints (pa)
   real(r8), intent(in) :: relhum      ! relative humidity with respective to liquid
   real(r8), intent(in) :: cldn        ! new value of cloud fraction    (fraction)
   real(r8), intent(in) :: qc          ! liquid water mixing ratio (kg/kg)
   real(r8), intent(in) :: qi          ! grid-mean preexisting cloud ice mass mixing ratio (kg/kg)
   real(r8), intent(in) :: ni_in       ! grid-mean preexisting cloud ice number conc (#/kg) 
   real(r8), intent(in) :: rhoair      ! air density (kg/m3)
   real(r8), intent(in) :: so4_num     ! so4 aerosol number (#/cm^3)
   real(r8), intent(in) :: dst_num     ! total dust aerosol number (#/cm^3)
   real(r8), intent(in) :: soot_num    ! soot (hydrophilic) aerosol number (#/cm^3)
   real(r8), intent(in) :: dst1_num     ! dust aerosol number (#/cm^3)
   real(r8), intent(in) :: dst2_num     ! dust aerosol number (#/cm^3)
   real(r8), intent(in) :: dst3_num     ! dust aerosol number (#/cm^3)
   real(r8), intent(in) :: dst4_num     ! dust aerosol number (#/cm^3)
   real(r8), intent(in) :: dst1_sfc_to_num
   real(r8), intent(in) :: dst3_sfc_to_num
   real(r8), intent(in) :: organic_num  !organic aerosol number (primary carbon) (#/cm^3)
   logical,  intent(in) :: clim_modal_aero !whether MAM is used or not

   ! Output Arguments
   real(r8), intent(out) :: nuci       ! ice number nucleated (#/kg)
   real(r8), intent(out) :: onihf      ! nucleated number from homogeneous freezing of so4
   real(r8), intent(out) :: oniimm     ! nucleated number from immersion freezing
   real(r8), intent(out) :: onidep     ! nucleated number from deposition nucleation
   real(r8), intent(out) :: onimey     ! nucleated number from deposition nucleation  (meyers: mixed phase)
   real(r8), intent(out) :: wpice      ! diagnosed Vertical velocity Reduction caused by preexisting ice (m/s), at Shom
   real(r8), intent(out) :: weff       ! effective Vertical velocity for ice nucleation (m/s); weff=wbar-wpice
   real(r8), intent(out) :: fhom       ! how much fraction of cloud can reach Shom

   ! Local workspace
   real(r8) :: nihf                      ! nucleated number from homogeneous freezing of so4
   real(r8) :: niimm                     ! nucleated number from immersion freezing
   real(r8) :: nidep                     ! nucleated number from deposition nucleation
   real(r8) :: nimey                     ! nucleated number from deposition nucleation (meyers)
   real(r8) :: n1, ni                    ! nucleated number
   real(r8) :: tc, A, B, regm            ! work variable
   real(r8) :: esl, esi, deles           ! work variable

   !(09/29/2014)DeMott for mixed-phase cloud ice nucleation 
   real(r8)  :: na500, na500_1                            ! aerosol number with D>500 nm (#/cm^3)
   real(r8)  :: na500stp                                  ! aerosol number with D>500 nm (#/cm^3) at STP
   real(r8)  :: nimeystp                                  ! nucleated number from ice nucleation (meyers) at STP
   real(r8)  :: ad, bd   
   real(r8)  :: wbar1, wbar2

   ! Niemand et al. for mixed-phase cloud immersion ice nucleation (surface area based, dust)
   real(r8)  :: an
   real(r8)  :: ns_dust_imm    ! dust active site surface densities from AIDA experiments (m-2)

   ! used in SUBROUTINE Vpreice
   real(r8) :: Ni_preice        ! cloud ice number conc (1/m3)   
   real(r8) :: lami,Ri_preice   ! mean cloud ice radius (m)
   real(r8) :: Shom             ! initial ice saturation ratio; if <1, use hom threshold Si
   real(r8) :: detaT,RHimean    ! temperature standard deviation, mean cloudy RHi
   real(r8) :: wpicehet   ! diagnosed Vertical velocity Reduction caused by preexisting ice (m/s), at shet

   real(r8) :: weffhet    ! effective Vertical velocity for ice nucleation (m/s)  weff=wbar-wpicehet 
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
   deles = 0._r8
   esi   = 0._r8

   if(so4_num >= 1.0e-10_r8 .and. (soot_num+dst3_num) >= 1.0e-10_r8 .and. cldn > 0._r8) then

#ifdef USE_XLIU_MOD
!++ Mod from Xiaohong is the following two line conditional.
!   It changes answers so needs climate validation.
      if ((relhum*svp_water(tair)/svp_ice(tair)*subgrid).ge.1.2_r8) then
         if ( ((tc.le.0.0_r8).and.(tc.ge.-37.0_r8).and.(qc.lt.1.e-12_r8)).or.(tc.le.-37.0_r8)) then
#else
      if((tc.le.-35.0_r8) .and. ((relhum*svp_water(tair)/svp_ice(tair)*subgrid).ge.1.2_r8)) then ! use higher RHi threshold
#endif

            A = -1.4938_r8 * log(soot_num+dst3_num) + 12.884_r8
            B = -10.41_r8  * log(soot_num+dst3_num) - 67.69_r8

            regm = A * log(wbar1) + B

            ! heterogeneous nucleation only
            if (tc .gt. regm) then

               if(tc.lt.-40._r8 .and. wbar1.gt.1._r8) then ! exclude T<-40 & W>1m/s from hetero. nucleation

                  call hf(tc,wbar1,relhum,so4_num,nihf)
                  niimm=0._r8
                  nidep=0._r8   
                  n1=nihf
                  
               else

                  call hetero(tc,wbar2,soot_num+dst3_num,niimm,nidep)
                  nihf=0._r8
                  n1=niimm+nidep

               endif

            ! homogeneous nucleation only
            else if (tc.lt.regm-5._r8) then

               call hf(tc,wbar1,relhum,so4_num,nihf)
               niimm=0._r8
               nidep=0._r8
               n1=nihf
              
            ! transition between homogeneous and heterogeneous: interpolate in-between
            else

               if (tc.lt.-40._r8 .and. wbar1.gt.1._r8) then ! exclude T<-40 & W>1m/s from hetero. nucleation

                  call hf(tc, wbar1, relhum, so4_num, nihf)
                  niimm = 0._r8
                  nidep = 0._r8
                  n1 = nihf
                  
               else

                  call hf(regm-5._r8,wbar1,relhum,so4_num,nihf)
                  call hetero(regm,wbar2,soot_num+dst3_num,niimm,nidep)

                  if (nihf .le. (niimm+nidep)) then
                     n1 = nihf
                  else
                     n1=(niimm+nidep)*((niimm+nidep)/nihf)**((tc-regm)/5._r8)
                  endif

               end if
            end if

            ni = n1

         end if
      end if
#ifdef USE_XLIU_MOD
   end if
#endif

   ! deposition/condensation nucleation in mixed clouds (-37<T<0C) (Meyers, 1992)
   if(tc.lt.0._r8 .and. tc.gt.-37._r8 .and. qc.gt.1.e-12_r8) then
      if (use_dem_nucleate) then          ! use DeMott et al.         
         !++iceMP
         nimeystp = 1.e-3_r8 * 3.0_r8 * (na500stp**ad) * exp(bd)               ! cm^-3 
         nimey=nimeystp*273.15_r8*pmid/(101325_r8*tair)
      else if (use_nie_nucleate) then
         ns_dust_imm = exp(an)          ! m^-2
         nimey = dst1_num * (1._r8 - exp(-1.0_r8 * dst1_sfc_to_num * ns_dust_imm))&   ! cm^-3
               + dst3_num * (1._r8 - exp(-1.0_r8 * dst3_sfc_to_num * ns_dust_imm))

      else
         !--iceMP
         esl = svp_water(tair)     ! over water in mixed clouds
         esi = svp_ice(tair)     ! over ice
         deles = (esl - esi)
         nimey=1.e-3_r8*exp(12.96_r8*deles/esi - 0.639_r8) 
      endif 
   else
      nimey=0._r8
   endif

   if (use_hetfrz_classnuc) nimey = 0._r8

   nuci=ni+nimey
   if(nuci.gt.9999._r8.or.nuci.lt.0._r8) then
      write(iulog, *) 'Warning: incorrect ice nucleation number (nuci reset =0)'
      write(iulog, *) ni, tair, relhum, wbar, nihf, niimm, nidep,deles,esi,dst_num,so4_num
      nuci=0._r8
   endif

   nuci   = nuci*1.e+6_r8/rhoair    ! change unit from #/cm3 to #/kg
   onimey = nimey*1.e+6_r8/rhoair
   onidep = nidep*1.e+6_r8/rhoair
   oniimm = niimm*1.e+6_r8/rhoair
   onihf  = nihf*1.e+6_r8/rhoair

end subroutine nucleati

!===============================================================================

subroutine hetero(T,ww,Ns,Nis,Nid)

    real(r8), intent(in)  :: T, ww, Ns
    real(r8), intent(out) :: Nis, Nid

    real(r8) A11,A12,A21,A22,B11,B12,B21,B22
    real(r8) B,C

!---------------------------------------------------------------------
! parameters

      A11 = 0.0263_r8
      A12 = -0.0185_r8
      A21 = 2.758_r8
      A22 = 1.3221_r8
      B11 = -0.008_r8
      B12 = -0.0468_r8
      B21 = -0.2667_r8
      B22 = -1.4588_r8

!     ice from immersion nucleation (cm^-3)

      B = (A11+B11*log(Ns)) * log(ww) + (A12+B12*log(Ns))
      C =  A21+B21*log(Ns)

      Nis = exp(A22) * Ns**B22 * exp(B*T) * ww**C
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

      if((Temperature.le.-37.0_r8) .and. ((RH*subgrid).ge.RHw)) then

        regm = 6.07_r8*lnw-55.0_r8

        if(Temperature.ge.regm) then    ! fast-growth regime

          if(Temperature.gt.-64.0_r8) then
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

end subroutine


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

end subroutine

end module nucleate_ice

