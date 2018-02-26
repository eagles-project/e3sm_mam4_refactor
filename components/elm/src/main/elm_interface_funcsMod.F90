module elm_interface_funcsMod
!=================================================================================================
! CLM Theraml-Hydrology (TH) & BioGeoChemistry (BGC) Interface: Modules
! created: 8/25/2015
! updated: June-2017
! update: 5/13/2019 (all are IMPLICITLY column-based coupling,
!        esp. after ELM v2 data-structure modification @3/12/2019, commit 1bf22e32d)
!=================================================================================================

#include "shr_assert.h"


  ! MODULE: elm_interface_funcsMod
  !--------------------------------------------------------------------------------------
  ! DESCRIPTION:
  ! Coupling of CLM with any specific Soil BGC module Consists of 3 STEPS:
  ! STEP-1:   elm vars             -> elm_interface_data (i.e. elm_interface_dataType)  ; pass elm vars to elm_interface_data
  ! STEP-2:   elm_interface_data   -> soil bgc module -> elm_interface_data
  !      2.1: elm_interface_data   -> soil bgc module
  !      2.2: run soil bgc module
  !      2.3: soil bgc module      -> elm_interface_data
  ! STEP-3:   elm_interface_data   -> elm vars
  !--------------------------------------------------------------------------------------


  !--------------------------------------------------------------------------------------
  ! !USES:
  ! elm g/l/c/p constants
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use shr_kind_mod          , only : r8 => shr_kind_r8

  ! global variables
  use GridcellType          , only : grc_pp
  use LandunitType          , only : lun_pp
  use ColumnType            , only : col_pp 
  use VegetationType        , only : veg_pp
  use ColumnDataType        , only : col_es, col_ef
  use ColumnDataType        , only : col_ws, col_wf
  use ColumnDataType        , only : col_cs, col_cf
  use ColumnDataType        , only : col_ns, col_nf
  use ColumnDataType        , only : col_ps, col_pf

  ! (dummy) variable definitions
  use decompMod             , only : bounds_type

  use atm2lndType           , only : atm2lnd_type
  use SoilStateType         , only : soilstate_type
  use SoilHydrologyType     , only : soilhydrology_type

  use CNStateType           , only : cnstate_type

  use CH4Mod                , only : ch4_type

  use PhotosynthesisType    , only : photosyns_type
  use cropType              , only : crop_type
  use CanopyStateType       , only : canopystate_type

  use SoilWaterRetentionCurveMod    , only : soil_water_retention_curve_type

  use elm_interface_dataType        , only : elm_interface_data_type
  use elm_interface_thType          , only : elm_interface_th_datatype
  use elm_interface_bgcType         , only : elm_interface_bgc_datatype

  ! most used constants in this module
  use elm_varpar            , only : nlevsoi, nlevsno, nlevgrnd, nlevdecomp_full
  use elm_varpar            , only : ndecomp_pools, ndecomp_cascade_transitions
  use elm_varpar            , only : max_patch_per_col
  use elm_varcon            , only : denh2o, denice, tfrz, dzsoi_decomp
  use ColumnDataType        , only : column_water_state, column_water_flux
  use ColumnDataType        , only : column_energy_state, column_energy_flux
  use ColumnDataType        , only : column_carbon_state, column_carbon_flux
  use ColumnDataType        , only : column_nitrogen_state, column_nitrogen_flux
  use ColumnDataType        , only : column_phosphorus_state, column_phosphorus_flux

  use landunit_varcon       , only : istsoil, istcrop

  ! misc.
  use abortutils            , only : endrun
  use elm_varctl            , only : nu_com

  use timeinfoMod
  !--------------------------------------------------------------------------------------

  implicit none

  save

  private    ! By default everything is private

  ! LOCAL VARIABLES:


  !--------------------------------------------------------------------------------------
  ! (1) GENERIC SUBROUTINES: used by any specific soil BGC/TH module
  ! pass elm variables to elm_interface_data
  public    :: get_elm_data                 !! STEP-1: elm vars -> elm_interface_data

  !! pass elm variables to elm_interface_data, called by get_elm_data
  private   :: get_elm_soil_property        !! STEP-1.1: soil properties
  private   :: get_elm_soil_th_state        !! STEP-1.2: thermohydrology (TH) state vars
  private   :: get_elm_soil_th_flux         !! STEP-1.3: thermohydrology (TH) flux vars
  private   :: get_elm_bgc_state            !! STEP-1.4: state vars
  private   :: get_elm_bgc_flux             !! STEP-1.5: flux vars

  ! STEP-3.x: elm_interface_data -> elm vars
  ! update elm variables from elm_interface_data,
  ! e.g., called in 'update_bgc_data_elm2elm' and 'update_bgc_data_pf2elm'
  ! specific bgc-module (e.g., PFLOTRAN) requires certain combination of these subroutines
  private   :: update_bgc_state_decomp
  private   :: update_bgc_state_smin
  private   :: update_bgc_flux_decomp_sourcesink
  private   :: update_bgc_flux_decomp_cascade
  private   :: update_bgc_flux_smin
  private   :: update_bgc_flux_nitdenit
  private   :: update_bgc_flux_gas_pf
  private   :: update_soil_moisture
  private   :: update_soil_temperature

  !--------------------------------------------------------------------------------------
  ! (2) SPECIFIC SUBROUTINES: used by a specific soil BGC module
  ! (2.1) Specific Subroutines for running elm-bgc (CN or BGC) through interface
  ! if (use_elm_interface .and. use_elm_bgc)
  public    :: elm_bgc_run              ! STEP-2:   elm_interface_data  -> elm-bgc module -> elm_interface_data    ; called in elm_driver
  private   :: elm_bgc_get_data         ! STEP-2.1: elm_interface_data  -> elm-bgc module                          ; called in elm_bgc_run
                                        ! STEP-2.2: eun elm-bgc module                                             ; see SoilLittDecompAlloc in SoilLittDecompMod
  private   :: elm_bgc_update_data      ! STEP-2.3: elm-bgc module-> elm_interface_data                            ; called in elm_bgc_run
  public    :: update_bgc_data_elm2elm  ! STEP-3:   elm_interface_data  -> elm vars                                ; called in elm_driver

  ! (2.2) Specific Subroutines for CLM-PFLOTRAN Coupling: update elm variables from pflotran
  ! if (use_elm_interface .and. use_pflotran)
  public    :: update_bgc_data_pf2elm   ! STEP-3:   elm_interface_data  -> elm vars                                ; called in elm_driver
                                        ! STEP-2:   see 'elm_pf_run' in elm_interface_pflotranMod

  public    :: update_th_data_pf2elm
  !--------------------------------------------------------------------------------------

contains


!--------------------------------------------------------------------------------------
  subroutine get_elm_data(elm_idata,                              &
           bounds, num_soilc, filter_soilc,                       &
           num_soilp, filter_soilp,                               &
           atm2lnd_vars, soilstate_vars,                          &
           cnstate_vars,                                          &
           ch4_vars                                               &
           )

    implicit none

    ! !ARGUMENTS:
    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                     , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                     , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                     , intent(in)    :: filter_soilp(:)   ! filter for soil patches
    type(atm2lnd_type)          , intent(in)    :: atm2lnd_vars
    type(soilstate_type)        , intent(in)    :: soilstate_vars

    type(cnstate_type)          , intent(in)    :: cnstate_vars
    type(ch4_type)              , intent(in)    :: ch4_vars

    type(elm_interface_data_type), intent(inout) :: elm_idata

    ! LOCAL

    character(len=256) :: subname = "get_elm_data"
    !-----------------------------------------------------------------------

     associate ( &
      elm_idata_th  => elm_idata%th,  &
      elm_idata_bgc => elm_idata%bgc  &
     )

    call get_elm_soil_property(elm_idata,                   &
                    bounds, num_soilc, filter_soilc,        &
                    soilstate_vars, cnstate_vars)

    call get_elm_soil_th_state(elm_idata_th,                &
                    bounds, num_soilc, filter_soilc,        &
                    atm2lnd_vars, soilstate_vars,           &
                    col_ws, col_es)                           ! global vars

    call get_elm_soil_th_flux(elm_idata_th,                 &
                    bounds, num_soilc, filter_soilc,        &
                    col_wf, col_ef)                           ! global vars

    call get_elm_bgc_state(elm_idata_bgc,                   &
                    bounds, num_soilc, filter_soilc,        &
                    atm2lnd_vars, soilstate_vars,           &
                    col_cs, col_ns, col_ps,                 & ! global vars
                    ch4_vars)

    call get_elm_bgc_flux(elm_idata_bgc,                    &
                    bounds, num_soilc, filter_soilc,        &
                    cnstate_vars,                           &
                    col_cf, col_nf, col_pf,                 & ! global vars
                    ch4_vars)

    end associate
  end subroutine get_elm_data
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
  subroutine get_elm_soil_property(elm_idata,               &
                        bounds, num_soilc, filter_soilc,    &
                        soilstate_vars, cnstate_vars)

    !
    ! !DESCRIPTION:
    ! get soil column physical and biogeochemical properties
    !
    ! !USES:
    use CNDecompCascadeConType, only : decomp_cascade_con
    use elm_varctl            , only : spinup_state
    use elm_varcon            , only : tkice, thk_bedrock

    ! !ARGUMENTS:

    implicit none

    ! LOCAL VARAIBLES:

    type(bounds_type)        , intent(in) :: bounds           ! bounds
    integer                  , intent(in) :: num_soilc        ! number of column soil points in column filter
    integer                  , intent(in) :: filter_soilc(:)  ! column filter for soil points
    type(soilstate_type)     , intent(in) :: soilstate_vars
    type(cnstate_type)       , intent(in) :: cnstate_vars

    type(elm_interface_data_type), intent(inout) :: elm_idata

    integer  :: fc, g, l, c, j, k      ! indices
    integer  :: gcount, cellcount

    character(len= 32) :: subname = 'get_elm_soil_property' ! subroutine name

    associate ( &
         ! Assign local pointer to derived subtypes components (column-level)
         z                  => col_pp%z                                             , & !  [real(r8) (:,:)]  layer depth (m)
         dz                 => col_pp%dz                                            , & !  [real(r8) (:,:)]  layer thickness depth (m)
         zi                 => col_pp%zi                                            , & !
         snl                => col_pp%snl                                           , & !

         bd                 => soilstate_vars%bd_col                                , & !
         bsw                => soilstate_vars%bsw_col                               , & !  [real(r8) (:,:)]  Clapp and Hornberger "b" (nlevgrnd)
         hksat              => soilstate_vars%hksat_col                             , & !  [real(r8) (:,:)]  hydraulic conductivity at saturation (mm H2O /s) (nlevgrnd)
         sucsat             => soilstate_vars%sucsat_col                            , & !  [real(r8) (:,:)]  minimum soil suction (mm) (nlevgrnd)
         watsat             => soilstate_vars%watsat_col                            , & !  [real(r8) (:,:)]  volumetric soil water at saturation (porosity) (nlevgrnd)
         watfc              => soilstate_vars%watfc_col                             , & !  [real(r8) (:,:)]  volumetric soil water at saturation (porosity) (nlevgrnd)
         watmin             => soilstate_vars%watmin_col                            , & !  [real(r8) (:,:)]  col minimum volumetric soil water (nlevsoi)
         sucmin             => soilstate_vars%sucmin_col                            , & !  [real(r8) (:,:)]  col minimum allowable soil liquid suction pressure (mm) [Note: sucmin is a negative value, while sucsat is a positive quantity]
         !
         tkmg               => soilstate_vars%tkmg_col                              , & !  [real(r8) (:,:)]  ! col thermal conductivity, soil minerals  [W/m-K] (nlevgrnd)
         tkdry              => soilstate_vars%tkdry_col                             , & !  [real(r8) (:,:)]  ! col thermal conductivity, dry soil [W/m-K] (nlevgrnd)
         tksatu             => soilstate_vars%tksatu_col                            , & !  [real(r8) (:,:)]  ! col thermal conductivity, saturated soil [W/m-K] (nlevgrnd)
         hcsoil             => soilstate_vars%csol_col                              , & !  [real(r8) (:,:)]  ! col heat capacity, soil solids (J/m**3/Kelvin) (nlevgrnd)
         !
         cellorg            => soilstate_vars%cellorg_col                           , & !  Input:  [real(r8) (:,:)  ]  column 3D org (kg/m3 organic matter) (nlevgrnd)
         !
         !porosity           => soilstate_vars%porosity_col                         , &
         eff_porosity       => soilstate_vars%eff_porosity_col                      , &
         !
         rootfr             => soilstate_vars%rootfr_col                            , & ! pft-level effective fraction of roots in each soil layer
         !
         initial_cn_ratio   => decomp_cascade_con%initial_cn_ratio                  , &
         initial_cp_ratio   => decomp_cascade_con%initial_cp_ratio                  , &
         !
         decomp_pool_name   => decomp_cascade_con%decomp_pool_name_history          , &
         floating_cn_ratio  => decomp_cascade_con%floating_cn_ratio_decomp_pools    , &
         floating_cp_ratio  => decomp_cascade_con%floating_cp_ratio_decomp_pools    , &
         decomp_k_pools     => decomp_cascade_con%decomp_k_pools                    , &
         adfactor_kd_pools  => decomp_cascade_con%spinup_factor                     , &
         rf_decomp_cascade       => cnstate_vars%rf_decomp_cascade_col              , &
         pathfrac_decomp_cascade => cnstate_vars%pathfrac_decomp_cascade_col          &

         )

!-------------------------------------------------------------------------------------
    ! constants:
    elm_idata%bgc%ndecomp_pools          = ndecomp_pools
    elm_idata%bgc%decomp_pool_name(:)    = decomp_pool_name(:)
    elm_idata%bgc%floating_cn_ratio(:)   = floating_cn_ratio(:)
    elm_idata%bgc%floating_cp_ratio(:)   = floating_cp_ratio(:)

    elm_idata%bgc%initial_cn_ratio(:)    = initial_cn_ratio(:)
    elm_idata%bgc%initial_cp_ratio(:)    = initial_cp_ratio(:)
    elm_idata%bgc%decomp_k_pools(:)      = decomp_k_pools(1:ndecomp_pools)
    if (spinup_state == 1) then
        elm_idata%bgc%adfactor_kd_pools(:) = adfactor_kd_pools(1:ndecomp_pools)
    else
        elm_idata%bgc%adfactor_kd_pools(:) = 1.0_r8
    end if

    do fc = 1, num_soilc
        c = filter_soilc(fc)

        elm_idata%z(c,:)             = z(c,:)
        elm_idata%zi(c,:)            = zi(c,:)
        elm_idata%dz(c,:)            = dz(c,:)
        elm_idata%snl(c)             = snl(c)
        elm_idata%bd(c,:)            = bd(c,:)
        elm_idata%bsw(c,:)           = bsw(c,:)
        elm_idata%hksat(c,:)         = hksat(c,:)
        elm_idata%sucsat(c,:)        = sucsat(c,:)
        elm_idata%watsat(c,:)        = watsat(c,:)
        elm_idata%watfc(c,:)         = watfc(c,:)
        elm_idata%watmin(c,:)        = watmin(c,:)
        elm_idata%sucmin(c,:)        = sucmin(c,:)

        elm_idata%porosity(c,:)      = watsat(c,:)
        elm_idata%eff_porosity(c,:)  = eff_porosity(c,:)

        elm_idata%cellorg(c,:)       = cellorg(c,:)

        elm_idata%rootfr(c,:)        = rootfr(c,:)

        elm_idata%tkwet(c,:)         = tksatu(c,:)
        elm_idata%tkdry(c,:)         = tkdry(c,:)
        elm_idata%tkfrz(c,:)         = tkice
        elm_idata%csol(c,:)          = hcsoil(c,:)

        !
        do k = 1, ndecomp_cascade_transitions
            elm_idata%bgc%rf_decomp_cascade(c,:,k)           = rf_decomp_cascade(c,:,k)
            elm_idata%bgc%pathfrac_decomp_cascade(c,:,k)     = pathfrac_decomp_cascade(c,:,k)
        end do

    end do

  end associate
  end subroutine get_elm_soil_property
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
  subroutine get_elm_soil_th_state(elm_idata_th,            &
                       bounds, num_soilc, filter_soilc,     &
                       atm2lnd_vars, soilstate_vars,        &
                       waterstate_vars, energystate_vars)
  !
  ! !DESCRIPTION:
  !  get soil temperature/saturation from CLM to soil BGC module
  !
  ! !USES:
    use shr_const_mod       , only : SHR_CONST_G


  ! !ARGUMENTS:
    implicit none

    type(bounds_type)        , intent(in) :: bounds           ! bounds
    integer                  , intent(in) :: num_soilc        ! number of column soil points in column filter
    integer                  , intent(in) :: filter_soilc(:)  ! column filter for soil points
    type(atm2lnd_type)       , intent(in) :: atm2lnd_vars
    type(soilstate_type)     , intent(in) :: soilstate_vars

    ! the following currently at column-level, but may be re-defined as grid-level if available (2019-10-28, fmy)
    type(column_water_state) , intent(in) :: waterstate_vars
    type(column_energy_state), intent(in) :: energystate_vars

    type(elm_interface_th_datatype)  , intent(inout) :: elm_idata_th

  ! !LOCAL VARIABLES:
    integer  :: fc, g, c, j         ! indices

  !EOP
  !-----------------------------------------------------------------------
    associate ( &
      soilpsi               => waterstate_vars%soilp                , & ! [real(r8) (:,:) ]  soil water pressure (Pa)
      frac_sno_eff          => waterstate_vars%frac_sno_eff         , & !
      frac_h2osfc           => waterstate_vars%frac_h2osfc          , & !
      h2osoi_vol            => waterstate_vars%h2osoi_vol           , & ! [real(r8) (:,:)] volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevgrnd)
      h2osoi_liq            => waterstate_vars%h2osoi_liq           , & ! [real(r8) (:,:)] liquid water (kg/m2) (-nlevsno+1:nlevgrnd)
      h2osoi_ice            => waterstate_vars%h2osoi_ice           , & ! [real(r8) (:,:)] ice lens (kg/m2) (-nlevsno+1:nlevgrnd)
      h2osfc                => waterstate_vars%h2osfc               , & ! [real(r8) (:)] surface water (mmH2O)
      !
      t_soisno              => energystate_vars%t_soisno            , & ! [real(r8) (:,:)] snow-soil temperature (Kelvin) (-nlevsno+1:nlevgrnd)
      t_grnd                => energystate_vars%t_grnd              , & ! [real(r8) (:)] ground (snow/soil1/surfwater-mixed) temperature (Kelvin)
      t_h2osfc              => energystate_vars%t_h2osfc            , & ! [real(r8) (:)] surface water temperature (Kelvin)
      t_nearsurf            => energystate_vars%t_nearsurf          , & ! [real(r8) (:)] near surface air temperature (Kelvin)
      !
      forc_pbot             => atm2lnd_vars%forc_pbot_not_downscaled_grc  & ! atmospheric pressure (Pa)
      )


    !--------------------------------------------------------------------------------------
    do fc = 1,num_soilc
        c = filter_soilc(fc)
        g = col_pp%gridcell(c)

        elm_idata_th%forc_pbot(c)            = forc_pbot(g)

        elm_idata_th%frac_sno_eff(c)         = frac_sno_eff(c)
        elm_idata_th%frac_h2osfc(c)          = frac_h2osfc(c)
        elm_idata_th%h2osfc(c)               = h2osfc(c)

        elm_idata_th%t_grnd(c)               = t_grnd(c)
        elm_idata_th%t_h2osfc(c)             = t_h2osfc(c)
        elm_idata_th%t_nearsurf(c)           = t_nearsurf(c)

        do j = -nlevsno+1,nlevgrnd
            if(j>=1) then
                elm_idata_th%soilpsi(c,j)        = soilpsi(c,j)
                elm_idata_th%h2osoi_vol(c,j)     = h2osoi_vol(c,j)
            endif

            elm_idata_th%h2osoi_liq(c,j)         = h2osoi_liq(c,j)
            elm_idata_th%h2osoi_ice(c,j)         = h2osoi_ice(c,j)
            elm_idata_th%t_soisno(c,j)           = t_soisno(c,j)
        end do

    end do

    end associate
  end subroutine get_elm_soil_th_state
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
  subroutine get_elm_soil_th_flux(elm_idata_th,             &
                       bounds, num_soilc, filter_soilc,     &
                       waterflux_vars, energyflux_vars)
  !
  ! !DESCRIPTION:
  !  get soil thermal-hydrological fluxes from CLM to soil TH module
  !
  ! !USES:
    use shr_const_mod       , only : SHR_CONST_G


  ! !ARGUMENTS:
    implicit none

    type(bounds_type)        , intent(in) :: bounds           ! bounds
    integer                  , intent(in) :: num_soilc        ! number of column soil points in column filter
    integer                  , intent(in) :: filter_soilc(:)  ! column filter for soil points

    ! the following currently at column-level, but may be re-defined as grid-level if available (2019-10-28, fmy)
    type(column_water_flux)  , intent(in) :: waterflux_vars
    type(column_energy_flux) , intent(in) :: energyflux_vars

    type(elm_interface_th_datatype)       , intent(inout) :: elm_idata_th

  ! !LOCAL VARIABLES:
    integer  :: fc, c, j         ! indices

  !EOP
  !-----------------------------------------------------------------------
    associate ( &
      qflx_top_soil     => waterflux_vars%qflx_top_soil         , & ! [real(:)] ! net liq. water input into top of soil column (mmH2O/s)
      qflx_evap_soil    => waterflux_vars%qflx_ev_soil          , & ! [real(:)] ! col soil surface evaporation (mm H2O/s) (+ = to atm)
      qflx_evap_h2osfc  => waterflux_vars%qflx_ev_h2osfc        , & ! [real(:)] ! col water surface evaporation (mm H2O/s) (+ = to atm)
      qflx_evap_snow    => waterflux_vars%qflx_ev_snow          , & ! [real(:)] ! col snow surface evaporation (mm H2O/s) (+ = to atm)
      qflx_subl_snow    => waterflux_vars%qflx_sub_snow         , & ! [real(:)] ! col snow sublimation (mm H2O/s) (+ = to atm)
      qflx_tran_veg     => waterflux_vars%qflx_tran_veg         , & ! [real(:)] ! col plant transpiration (mm H2O/s) (+ = to atm)
      qflx_rootsoil     => waterflux_vars%qflx_rootsoi          , & ! [real(:,:)] ! col vertically-resolved (nlevsoi) root and soil water exchange [mm H2O/s] [+ into root]
      !
      htvp              => energyflux_vars%htvp                 , & ! [real(:) ! latent heat of vapor of water (or sublimation) [j/kg]
      eflx_bot          => energyflux_vars%eflx_bot             , & ! [real(:) ! col heat flux from beneath the soil or ice column (W/m**2)
      eflx_soil_grnd    => energyflux_vars%eflx_soil_grnd       , & ! [real(:) ! col soil (ground) heat flux (W/m**2) [+ = into ground]
      eflx_fgr0_snow    => energyflux_vars%eflx_fgr0_snow       , & ! [real(:) ! col ground heat flux from snow bottom to first soil layer (W/m**2) [+ = into soil]
      eflx_fgr0_h2osfc  => energyflux_vars%eflx_fgr0_h2osfc     , & ! [real(:) ! col ground heat flux from surface water bottom to first soil layer (W/m**2) [+ = into soil]
      eflx_fgr0_soil    => energyflux_vars%eflx_fgr0_soil       , & ! [real(:) ! col ground heat flux from near-surface air to first soil layer (W/m**2) [+ = into soil]
      eflx_rnet_soil    => energyflux_vars%eflx_rnet_soil         & ! [real(:) ! net radiation flux between soil layer 1 and above-air, excluding SH and LE (i.e. radiation form only ) (W/m2) [+ = into soil]

    )

    ! a few notes:
    !   - 'qflx_evap_soil' appears for total soil surface, esp. bare soil; 'qflx_ev_soil/snow/h2osfc' are actually applied for in soil water modules
    !   - 'qflx_ev_snow' vs. 'qflx_sub_snow': the former is for total evap from both solid/liq., the latter is from solid snow pack (normally shall be same)
    !                        there is another variable 'qlfx_evap_grnd', which are those from liq. water when snow
    !--------------------------------------------------------------------------------------
!
    do fc = 1,num_soilc
        c = filter_soilc(fc)

        elm_idata_th%qflx_top_soil(c)        = qflx_top_soil(c)
        elm_idata_th%qflx_evap_soil(c)       = qflx_evap_soil(c)
        elm_idata_th%qflx_evap_h2osfc(c)     = qflx_evap_h2osfc(c)
        elm_idata_th%qflx_evap_snow(c)       = qflx_evap_snow(c)
        elm_idata_th%qflx_subl_snow(c)       = qflx_subl_snow(c)
        elm_idata_th%qflx_tran_veg(c)        = qflx_tran_veg(c)

        do j = 1,nlevgrnd
            if (j<=nlevsoi) then
               elm_idata_th%qflx_rootsoil(c,j) = qflx_rootsoil(c,j)
            else
               elm_idata_th%qflx_rootsoil(c,j) = 0._r8
            end if
        end do

        elm_idata_th%htvp(c)                 = htvp(c)
        elm_idata_th%eflx_bot(c)             = eflx_bot(c)
        elm_idata_th%eflx_soil_grnd(c)       = eflx_soil_grnd(c)
        elm_idata_th%eflx_fgr0_snow(c)       = eflx_fgr0_snow(c)
        elm_idata_th%eflx_fgr0_h2osfc(c)     = eflx_fgr0_h2osfc(c)
        elm_idata_th%eflx_fgr0_soil(c)       = eflx_fgr0_soil(c)
        elm_idata_th%eflx_rnet_soil(c)       = eflx_rnet_soil(c)

    end do

    end associate
  end subroutine get_elm_soil_th_flux
!--------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------
  subroutine get_elm_bgc_state(elm_bgc_data,                    &
                        bounds, num_soilc, filter_soilc,        &
                        atm2lnd_vars, soilstate_vars,           &
                        carbonstate_vars, nitrogenstate_vars,   &
                        phosphorusstate_vars,                   &
                        ch4_vars)

    ! get elm bgc state variables
    implicit none

    type(bounds_type)           , intent(in) :: bounds
    integer                     , intent(in) :: num_soilc         ! number of soil columns in filter
    integer                     , intent(in) :: filter_soilc(:)   ! filter for soil columns

    type(atm2lnd_type)             , intent(in) :: atm2lnd_vars
    type(soilstate_type)           , intent(in) :: soilstate_vars

    ! the following currently at column-level, but may be re-defined as grid-level if available (2019-10-28, fmy)
    type(column_carbon_state)      , intent(in) :: carbonstate_vars
    type(column_nitrogen_state)    , intent(in) :: nitrogenstate_vars
    type(column_phosphorus_state)  , intent(in) :: phosphorusstate_vars
    type(ch4_type)                 , intent(in) :: ch4_vars          ! not yet used, but will be.

    type(elm_interface_bgc_datatype), intent(inout) :: elm_bgc_data

    character(len=256) :: subname = "get_elm_bgc_state"

    ! Local variables
    integer  :: fc, g, c, j, k

    !------------------------------------------------------------------------------------------
    !
    associate ( &
       decomp_cpools_vr=> carbonstate_vars%decomp_cpools_vr     , & ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
       decomp_npools_vr=> nitrogenstate_vars%decomp_npools_vr   , & ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
       decomp_ppools_vr=> phosphorusstate_vars%decomp_ppools_vr , & ! [real(r8) (:,:,:) ! col (gP/m3) vertically-resolved decomposing (litter, cwd, soil) P pools
       smin_no3_vr     => nitrogenstate_vars%smin_no3_vr        , & ! (gN/m3) vertically-resolved soil mineral NO3
       smin_nh4_vr     => nitrogenstate_vars%smin_nh4_vr        , & ! (gN/m3) vertically-resolved soil mineral NH4
       smin_nh4sorb_vr => nitrogenstate_vars%smin_nh4sorb_vr    , & ! (gN/m3) vertically-resolved soil mineral NH4 absorbed
       !
       solutionp_vr    => phosphorusstate_vars%solutionp_vr     , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil solution P
       labilep_vr      => phosphorusstate_vars%labilep_vr       , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil labile mineral P
       secondp_vr      => phosphorusstate_vars%secondp_vr       , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil secondary mineralP
       sminp_vr        => phosphorusstate_vars%sminp_vr         , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil mineral P = solutionp + labilep + secondp
       occlp_vr        => phosphorusstate_vars%occlp_vr         , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil occluded mineral P
       primp_vr        => phosphorusstate_vars%primp_vr         , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil primary mineral P
       !
       forc_pco2             => atm2lnd_vars%forc_pco2_grc      , & ! partial pressure co2 (Pa)
       forc_pch4             => atm2lnd_vars%forc_pch4_grc      , & ! partial pressure ch4 (Pa)
       !
       o2stress_sat          => ch4_vars%o2stress_sat_col           , & ! Input:  [real(r8) (:,:)  ]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
       o2stress_unsat        => ch4_vars%o2stress_unsat_col         , & ! Input:  [real(r8) (:,:)  ]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
       finundated            => ch4_vars%finundated_col             , & ! Input:  [real(r8) (:)    ]  fractional inundated area (excluding dedicated wetland columns)
       o2_decomp_depth_unsat => ch4_vars%o2_decomp_depth_unsat_col  , & ! Input:  [real(r8) (:,:)  ]  O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
       conc_o2_unsat         => ch4_vars%conc_o2_unsat_col          , & ! Input:  [real(r8) (:,:)  ]  O2 conc in each soil layer (mol/m3) (nlevsoi)
       o2_decomp_depth_sat   => ch4_vars%o2_decomp_depth_sat_col    , & ! Input:  [real(r8) (:,:)  ]  O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
       conc_o2_sat           => ch4_vars%conc_o2_sat_col              & ! Input:  [real(r8) (:,:)  ]  O2 conc in each soil layer (mol/m3) (nlevsoi)
    )
!


    do fc = 1, num_soilc
        c = filter_soilc(fc)
        g = col_pp%gridcell(c)

        elm_bgc_data%forc_pco2(c)                   = forc_pco2(g)
        elm_bgc_data%forc_pch4(c)                   = forc_pch4(g)

        do k = 1, ndecomp_pools
            elm_bgc_data%decomp_cpools_vr(c,:,k)    = decomp_cpools_vr(c,:,k)
            elm_bgc_data%decomp_npools_vr(c,:,k)    = decomp_npools_vr(c,:,k)
            elm_bgc_data%decomp_ppools_vr(c,:,k)    = decomp_ppools_vr(c,:,k)

        end do

        elm_bgc_data%smin_no3_vr(c,:)           = smin_no3_vr(c,:)
        elm_bgc_data%smin_nh4_vr(c,:)           = smin_nh4_vr(c,:)
        !elm_bgc_data%smin_nh4sorb_vr(c,:)       = smin_nh4sorb_vr(c,:)
        elm_bgc_data%smin_nh4sorb_vr(c,:)       = 0._r8 ! currently not yet available from CLM bgc

        elm_bgc_data%solutionp_vr(c,:)          = solutionp_vr(c,:)
        elm_bgc_data%labilep_vr(c,:)            = labilep_vr(c,:)
        elm_bgc_data%secondp_vr(c,:)            = secondp_vr(c,:)
        elm_bgc_data%sminp_vr(c,:)              = solutionp_vr(c,:) + labilep_vr(c,:) + secondp_vr(c,:)
        elm_bgc_data%occlp_vr(c,:)              = occlp_vr(c,:)
        elm_bgc_data%primp_vr(c,:)              = primp_vr(c,:)

        elm_bgc_data%finundated(c)              = finundated(c)
        elm_bgc_data%o2stress_unsat(c,:)        = o2stress_unsat(c,:)
        elm_bgc_data%o2stress_sat(c,:)          = o2stress_sat(c,:)
        elm_bgc_data%o2_decomp_depth_unsat(c,:) = o2_decomp_depth_unsat(c,:)
        elm_bgc_data%conc_o2_unsat(c,:)         = conc_o2_unsat(c,:)
        elm_bgc_data%o2_decomp_depth_sat(c,:)   = o2_decomp_depth_sat(c,:)
        elm_bgc_data%conc_o2_sat(c,:)           = conc_o2_sat(c,:)

    end do

!-----------------------------------------------------------------------------

    end associate
  end subroutine get_elm_bgc_state
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
  subroutine get_elm_bgc_flux(elm_bgc_data,                      &
                        bounds, num_soilc, filter_soilc,         &
                        cnstate_vars, carbonflux_vars,           &
                        nitrogenflux_vars, phosphorusflux_vars,  &
                        ch4_vars)

  !
  ! !DESCRIPTION:
  ! get elm bgc flux variables: external inputs to bgc state variables (pools)
  !
  ! !USES:
    use elm_varctl            , only : spinup_state
    use CNDecompCascadeConType, only : decomp_cascade_con

  ! !ARGUMENTS:
    implicit none

    type(bounds_type) , intent(in) :: bounds
    integer           , intent(in) :: num_soilc       ! number of soil columns in filter
    integer           , intent(in) :: filter_soilc(:) ! filter for soil columns

    type(cnstate_type)                  , intent(in) :: cnstate_vars

    ! the following currently at column-level, but may be re-defined as grid-level if available (2019-10-28, fmy)
    type(column_carbon_flux)            , intent(in) :: carbonflux_vars
    type(column_nitrogen_flux)          , intent(in) :: nitrogenflux_vars
    type(column_phosphorus_flux)        , intent(in) :: phosphorusflux_vars

    type(ch4_type)                      , intent(in) :: ch4_vars          ! not yet used, but will be.

    type(elm_interface_bgc_datatype)   , intent(inout) :: elm_bgc_data

    character(len=256) :: subname = "get_elm_bgc_flux"


 ! !LOCAL VARIABLES:
    integer  :: fc, c, g, j, k                         ! do loop indices

    real(r8) :: dtime                               ! land model time step (sec)
    integer  :: year, mon, day, sec

    ! ratios of NH4:NO3 in N deposition and fertilization (temporarily set here, will be as inputs)
    real(r8) :: r_nh4_no3_dep(bounds%begc:bounds%endc)
    real(r8) :: r_nh4_no3_fert(bounds%begc:bounds%endc)
    real(r8) :: fnh4_dep, fnh4_fert

    ! C/N source/sink rates as inputs for soil bgc

    !
    !---------------------------------------------------------------------------
    !
    associate ( &
      ! plant litering and removal + SOM/LIT vertical transport
      externalc_to_decomp_cpools_vr    => carbonflux_vars%externalc_to_decomp_cpools        , &
      externaln_to_decomp_npools_vr    => nitrogenflux_vars%externaln_to_decomp_npools      , &
      externalp_to_decomp_ppools_vr    => phosphorusflux_vars%externalp_to_decomp_ppools    , &
      t_scalar                         => carbonflux_vars%t_scalar                          , & ! Output: [real(r8) (:,:)   ]  soil temperature scalar for decomp
      w_scalar                         => carbonflux_vars%w_scalar                          , & ! Output: [real(r8) (:,:)   ]  soil water scalar for decomp
      o_scalar                         => carbonflux_vars%o_scalar                          , & ! Output: [real(r8) (:,:)   ]  fraction by which decomposition is limited by anoxia
      ! inorg. nitrogen source
      ndep_to_sminn                    => nitrogenflux_vars%ndep_to_sminn                   , &
      nfix_to_sminn                    => nitrogenflux_vars%nfix_to_sminn                   , &
      fert_to_sminn                    => nitrogenflux_vars%fert_to_sminn                   , &
      soyfixn_to_sminn                 => nitrogenflux_vars%soyfixn_to_sminn                , &
      supplement_to_sminn_vr           => nitrogenflux_vars%supplement_to_sminn_vr          , &
      !
      nfixation_prof                   => cnstate_vars%nfixation_prof_col                   , &
      ndep_prof                        => cnstate_vars%ndep_prof_col                        , &
      !
      decomp_k_scalar                  => cnstate_vars%scalaravg_col                        , &
      !
      no3_net_transport_vr             => nitrogenflux_vars%no3_net_transport_vr            , &
      nh4_net_transport_vr             => nitrogenflux_vars%nh4_net_transport_vr            , &
      plant_ndemand_vr                 => nitrogenflux_vars%plant_ndemand_vr                , &
      plant_ndemand                    => nitrogenflux_vars%plant_ndemand                   , &
      plant_pdemand                    => phosphorusflux_vars%plant_pdemand                 , &
      plant_pdemand_vr                 => phosphorusflux_vars%plant_pdemand_vr              , &
      pdep_to_sminp                    => phosphorusflux_vars%pdep_to_sminp                 , &
      ! assume pdep_prof = ndep_prof
      fert_p_to_sminp                  => phosphorusflux_vars%fert_p_to_sminp               , &
      supplement_to_sminp_vr           => phosphorusflux_vars%supplement_to_sminp_vr        , &
      sminp_net_transport_vr           => phosphorusflux_vars%sminp_net_transport_vr          &
    )

    !
    year=year_curr
    mon=mod_curr
    day=day_curr
    sec=secs_curr

!
    r_nh4_no3_dep(:)  = 1.0_r8      ! temporarily assuming half of N dep is in NH4 and another half in NO3
    r_nh4_no3_fert(:) = 1.0_r8      ! temporarily assiming half of N fertilization is in NH4 and another half in NO3

!
    do fc = 1,num_soilc
        c = filter_soilc(fc)

        elm_bgc_data%t_scalar(c,:) = t_scalar(c,:)
        elm_bgc_data%w_scalar(c,:) = w_scalar(c,:)
        elm_bgc_data%o_scalar(c,:) = o_scalar(c,:)


        elm_bgc_data%plant_ndemand(c)   = plant_ndemand(c)
        elm_bgc_data%plant_pdemand(c)   = plant_pdemand(c)

        fnh4_dep  = max(0._r8, min(1.0_r8, 1._r8/(r_nh4_no3_dep(c)+1._r8)))
        fnh4_fert = max(0._r8, min(1.0_r8, 1._r8/(r_nh4_no3_fert(c)+1._r8)))

        do k = 1, ndecomp_pools
            elm_bgc_data%externalc_to_decomp_cpools(c,:,k)  = externalc_to_decomp_cpools_vr(c,:,k)
            elm_bgc_data%externaln_to_decomp_npools(c,:,k)  = externaln_to_decomp_npools_vr(c,:,k)
            elm_bgc_data%externalp_to_decomp_ppools(c,:,k)  = externalp_to_decomp_ppools_vr(c,:,k)
       end do

       ! the following is for CTC ad-spinup.
       ! There is a 'time' control here, so MUST be called each time-step,
       ! and then better put the code here rather than in 'get_elm_bgc_state'
       if (spinup_state == 1 .and. year >= 40 .and. nu_com .eq. 'RD') then
            elm_bgc_data%sitefactor_kd_vr(c,:) = decomp_k_scalar(c,:)
       else
            elm_bgc_data%sitefactor_kd_vr(c,:) = 1.0_r8
       end if



       elm_bgc_data%externaln_to_nh4(c,:)          =   fnh4_dep*ndep_to_sminn(c) * ndep_prof(c,:) +  &
                                                       fnh4_fert*fert_to_sminn(c) * ndep_prof(c,:) + &
                                                       fnh4_fert*supplement_to_sminn_vr(c,:) +       &
                                                       nfix_to_sminn(c) * nfixation_prof(c,:) +      &
                                                       soyfixn_to_sminn(c) * nfixation_prof(c,:)

       elm_bgc_data%externaln_to_no3(c,:)          =   (1._r8-fnh4_dep)*ndep_to_sminn(c) * ndep_prof(c, :) +  &
                                                       (1._r8-fnh4_fert)*fert_to_sminn(c) * ndep_prof(c, :) + &
                                                       (1._r8-fnh4_fert)*supplement_to_sminn_vr(c,:)



       elm_bgc_data%externalp_to_primp(c,:)        =   pdep_to_sminp(c)*ndep_prof(c, :)
       elm_bgc_data%externalp_to_labilep(c,:)      =   fert_p_to_sminp(c)*ndep_prof(c, :)
       elm_bgc_data%externalp_to_solutionp(c,:)    =   supplement_to_sminp_vr(c,:)

       elm_bgc_data%no3_net_transport_vr(c,:)      = no3_net_transport_vr(c,:)
       elm_bgc_data%nh4_net_transport_vr(c,:)      = nh4_net_transport_vr(c,:)
       elm_bgc_data%sminp_net_transport_vr(c,:)    = sminp_net_transport_vr(c,:)  !from solutionp

       elm_bgc_data%plant_ndemand_vr(c,:)          = plant_ndemand_vr(c,:)
       elm_bgc_data%plant_pdemand_vr(c,:)          = plant_pdemand_vr(c,:)

    end do ! fc = 1,num_soilc

    end associate
  end subroutine get_elm_bgc_flux
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
  subroutine update_soil_moisture(elm_idata_th,     &
           bounds, num_soilc, filter_soilc,         &
           soilstate_vars, soilhydrology_vars,      &
           waterstate_vars, waterflux_vars)

  !
  ! !DESCRIPTION:
  !
  !
  ! !USES:

  ! !ARGUMENTS:
    implicit none

    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_soilc        ! number of column soil points in column filter
    integer, intent(in) :: filter_soilc(:)  ! column filter for soil points
    type(soilstate_type),     intent(inout) :: soilstate_vars
    type(soilhydrology_type), intent(inout) :: soilhydrology_vars

    ! the following currently at column-level, but may be re-defined as grid-level if available (2019-10-28, fmy)
    type(column_water_state), intent(inout) :: waterstate_vars
    type(column_water_flux),  intent(inout) :: waterflux_vars

    type(elm_interface_th_datatype), intent(in) :: elm_idata_th

  ! !LOCAL VARIABLES:
    integer  :: fc, c, j, g, gcount      ! indices

  !EOP
  !-----------------------------------------------------------------------
    associate ( &
         soilpsi    =>  waterstate_vars%soilp           , &
         !
         h2osoi_liq =>  waterstate_vars%h2osoi_liq      , &
         h2osoi_ice =>  waterstate_vars%h2osoi_ice      , &
         h2osoi_vol =>  waterstate_vars%h2osoi_vol      , &
         !
         qflx_drain_perched     => waterflux_vars%qflx_drain_perched      , & ! Output: [real(r8) (:)   ]  sub-surface runoff from perched zwt (mm H2O /s)
         qflx_rsub_sat          => waterflux_vars%qflx_rsub_sat           , & ! Output: [real(r8) (:)   ]  soil saturation excess flow (exfiltraion) [mm h2o/s]
         qflx_drain             => waterflux_vars%qflx_drain              , & ! Output: [real(r8) (:)   ]  sub-surface runoff/drainage at bottom (mm H2O /s)
         qflx_lateral           => waterflux_vars%qflx_lateral            , & ! Output: [real(r8) (:)   ]  sub-surface runoff/drainage laterally (mm H2O /s)
         qflx_surf              => waterflux_vars%qflx_surf               , & ! Output: [real(r8) (:)   ]  soil surface runoff (mm H2O /s)
         qflx_h2osfc_surf       => waterflux_vars%qflx_h2osfc_surf        , & ! Output: [real(r8) (:)   ]  surface (pond) water runoff (mm/s)
         qflx_infl              => waterflux_vars%qflx_infl               , & ! Output: [real(r8) (:)   ]  infiltration (mm H2O /s)
         qflx_qrgwl             => waterflux_vars%qflx_qrgwl              , & ! Output: [real(r8) (:)   ]  qflx_surf at glaciers, wetlands, lakes
         qflx_runoff            => waterflux_vars%qflx_runoff               & ! Output: [real(r8) (:)   ]  total water losses to currents from column (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
    )

    ! states
    do fc = 1,num_soilc
        c = filter_soilc(fc)

        soilpsi(c,:)    =  elm_idata_th%soilpsi(c,:)
        ! currently 'col_ws%soilp' (in Pa) NOT really used elsewhere in ELM,
        ! rather the following 'soilstate_vars%soilpsi_col' (but in MPa)
        soilstate_vars%soilpsi_col(c,:) = elm_idata_th%soilpsi(c,:)*1.e-6_r8

        h2osoi_liq(c,:) =  elm_idata_th%h2osoi_liq(c,:)
        h2osoi_ice(c,:) =  elm_idata_th%h2osoi_ice(c,:)
        h2osoi_vol(c,:) =  elm_idata_th%h2osoi_vol(c,:)
    end do

    ! fluxes
    do fc = 1,num_soilc
        c = filter_soilc(fc)

        ! NOTES from pflotran coupling
        ! (1) 'qflx_drain_perched' setting to zero because included in 'qflx_drain'
        !    (if needed, it requires re-calculation from vertical-drain/lateral-flow in enclosed/local saturation zone)
        ! (2) 'qflx_drain' is the NET water flow in/out of soil column bottom interface.
        ! (3) 'qflx_surf' is the NET water flow in/out of soil column top interface.
        ! (4) 'qflx_lateral' is the NET water flow in/out of soil column sides interface.

        qflx_drain (c)        = elm_idata_th%qflx_drain(c)
        qflx_rsub_sat (c)     = elm_idata_th%qflx_exfl(c)
        qflx_infl (c)         = elm_idata_th%qflx_infl(c)
        qflx_surf (c)         = elm_idata_th%qflx_surf(c)
        qflx_h2osfc_surf(c)   = elm_idata_th%qflx_h2osfc(c)
        qflx_lateral(c)       = elm_idata_th%qflx_lateral(c)

        qflx_qrgwl(c)         = 0._r8
        qflx_drain_perched(c) = 0._r8

        ! add amount of ET adjusted by PFLOTRAN into 'qflx_surf' so that counted correctly in balance-checking
        ! NOTE: this is a work-around, especially when surface module not coupled with pflotran.
        qflx_surf(c) = qflx_surf(c) - elm_idata_th%qflx_et_reduced(c)

        !summary of all water loss to currents
        qflx_runoff(c) = qflx_drain(c) + qflx_surf(c)  + qflx_h2osfc_surf(c) + qflx_qrgwl(c) + qflx_drain_perched(c)


    end do


    end associate
  end subroutine update_soil_moisture
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
  subroutine update_soil_temperature(elm_idata_th,     &
           bounds, num_soilc, filter_soilc,            &
           energystate_vars)

  !
  ! !DESCRIPTION:
  !
  !
  ! !USES:

  ! !ARGUMENTS:
    implicit none

    type(bounds_type)                   , intent(in)    :: bounds
    integer                             , intent(in)    :: num_soilc        ! number of column soil points in column filter
    integer                             , intent(in)    :: filter_soilc(:)  ! column filter for soil points

    ! the following currently at column-level, but may be re-defined as grid-level if available (2019-10-28, fmy)
    type(column_energy_state)           , intent(inout) :: energystate_vars

    type(elm_interface_th_datatype)     , intent(in)    :: elm_idata_th

  ! !LOCAL VARIABLES:
    integer  :: fc, c, j, g, gcount      ! indices

  !EOP
  !-----------------------------------------------------------------------
    associate ( &
         t_soisno   => energystate_vars%t_soisno   & ! snow-soil temperature (Kelvin)
      )

      do fc = 1,num_soilc
         c = filter_soilc(fc)
         t_soisno(c,:)   = elm_idata_th%t_soisno(c,:)
      end do

    end associate
  end subroutine update_soil_temperature
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
  subroutine update_th_data_pf2elm(elm_idata_th,           &
           bounds, num_soilc, filter_soilc,                &
           soilstate_vars, soilhydrology_vars)

    ! USES
    use elm_varctl          , only : use_pflotran, pf_tmode, pf_hmode

    implicit none

    ! !ARGUMENTS:
    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                     , intent(in)    :: filter_soilc(:)   ! filter for soil columns

    type(soilstate_type)        , intent(inout) :: soilstate_vars
    type(soilhydrology_type)    , intent(inout) :: soilhydrology_vars

    type(elm_interface_th_datatype), intent(in) :: elm_idata_th

    !-----------------------------------------------------------------------

    character(len=256) :: subname = "update_th_data_pf2elm"

    if (pf_tmode) then
        call update_soil_temperature(elm_idata_th,      &
                   bounds, num_soilc, filter_soilc,     &
                   col_es)                                ! global vars
    end if

    if (pf_hmode) then
        call update_soil_moisture(elm_idata_th,         &
                   bounds, num_soilc, filter_soilc,     &
                   soilstate_vars, soilhydrology_vars,  &
                   col_ws, col_wf)                        ! global vars
    end if

  end subroutine update_th_data_pf2elm
!--------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------
  subroutine update_bgc_state_decomp(elm_bgc_data,  &
           bounds, num_soilc, filter_soilc,         &
           carbonstate_vars, nitrogenstate_vars,    &
           phosphorusstate_vars                     &
           )

    implicit none

    type(bounds_type)                   , intent(in)    :: bounds
    integer                             , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                             , intent(in)    :: filter_soilc(:)   ! filter for soil columns

    ! the following currently at column-level, but may be re-defined as grid-level if available (2019-10-28, fmy)
    type(column_carbon_state)              , intent(inout) :: carbonstate_vars
    type(column_nitrogen_state)            , intent(inout) :: nitrogenstate_vars
    type(column_phosphorus_state)          , intent(inout) :: phosphorusstate_vars

    type(elm_interface_bgc_datatype)    , intent(in)    :: elm_bgc_data

    character(len=256) :: subname = "update_soil_bgc_state"

    integer  :: fc,c,j,k

!------------------------------------------------------------------------------------
    !
    associate ( &
       decomp_cpools_vr             => carbonstate_vars%decomp_cpools_vr           , &
       decomp_npools_vr             => nitrogenstate_vars%decomp_npools_vr         , &
       decomp_ppools_vr             => phosphorusstate_vars%decomp_ppools_vr         &
       )
! ------------------------------------------------------------------------
!
       do fc = 1, num_soilc
            c = filter_soilc(fc)
            do k = 1, ndecomp_pools
                decomp_cpools_vr(c,:,k) = elm_bgc_data%decomp_cpools_vr(c,:,k)
                decomp_npools_vr(c,:,k) = elm_bgc_data%decomp_npools_vr(c,:,k)
                decomp_ppools_vr(c,:,k) = elm_bgc_data%decomp_ppools_vr(c,:,k)
            end do
       end do

    end associate
  end subroutine update_bgc_state_decomp
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
  subroutine update_bgc_state_smin(elm_bgc_data,      &
           bounds, num_soilc, filter_soilc,             &
           nitrogenstate_vars, phosphorusstate_vars)

    use CNDecompCascadeConType, only : decomp_cascade_con

    implicit none

    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                     , intent(in)    :: filter_soilc(:)   ! filter for soil columns

    ! the following currently at column-level, but may be re-defined as grid-level if available (2019-10-28, fmy)
    type(column_nitrogen_state)    , intent(inout) :: nitrogenstate_vars
    type(column_phosphorus_state)  , intent(inout) :: phosphorusstate_vars

    type(elm_interface_bgc_datatype), intent(in) :: elm_bgc_data

    character(len=256) :: subname = "update_bgc_state_smin"

    integer  :: fc,c,j

!------------------------------------------------------------------------------------
    !
    associate ( &
      sminn_vr           => nitrogenstate_vars%sminn_vr           , &
      smin_no3_vr        => nitrogenstate_vars%smin_no3_vr        , &
      smin_nh4_vr        => nitrogenstate_vars%smin_nh4_vr        , &
      smin_nh4sorb_vr    => nitrogenstate_vars%smin_nh4sorb_vr    , &

      solutionp_vr       => phosphorusstate_vars%solutionp_vr     , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil solution P
      labilep_vr         => phosphorusstate_vars%labilep_vr       , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil labile mineral P
      secondp_vr         => phosphorusstate_vars%secondp_vr       , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil secondary mineralP
      sminp_vr           => phosphorusstate_vars%sminp_vr         , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil mineral P = solutionp + labilep + second
      occlp_vr           => phosphorusstate_vars%occlp_vr         , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil occluded mineral P
      primp_vr           => phosphorusstate_vars%primp_vr           & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil primary mineral P

      )
! ------------------------------------------------------------------------
!
      do fc = 1, num_soilc
            c = filter_soilc(fc)
            smin_no3_vr(c,:)        = elm_bgc_data%smin_no3_vr(c,:)
            smin_nh4_vr(c,:)        = elm_bgc_data%smin_nh4_vr(c,:)
            smin_nh4sorb_vr(c,:)    = elm_bgc_data%smin_nh4sorb_vr(c,:)
            sminn_vr(c,:)           = elm_bgc_data%sminn_vr(c,:)

            solutionp_vr(c,:)       = elm_bgc_data%solutionp_vr(c,:)
            labilep_vr(c,:)         = elm_bgc_data%labilep_vr(c,:)
            secondp_vr(c,:)         = elm_bgc_data%secondp_vr(c,:)
            sminp_vr(c,:)           = elm_bgc_data%sminp_vr(c,:)
            occlp_vr(c,:)           = elm_bgc_data%occlp_vr(c,:)
            primp_vr(c,:)           = elm_bgc_data%primp_vr(c,:)
      end do

    end associate
  end subroutine update_bgc_state_smin
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
  subroutine update_bgc_flux_decomp_sourcesink(elm_bgc_data,       &
           bounds, num_soilc, filter_soilc,                          &
           carbonflux_vars, nitrogenflux_vars,                       &
           phosphorusflux_vars)

    use CNDecompCascadeConType, only : decomp_cascade_con

    implicit none

    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                     , intent(in)    :: filter_soilc(:)   ! filter for soil columns

    ! the following currently at column-level, but may be re-defined as grid-level if available (2019-10-28, fmy)
    type(column_carbon_flux)    , intent(inout) :: carbonflux_vars
    type(column_nitrogen_flux)  , intent(inout) :: nitrogenflux_vars
    type(column_phosphorus_flux), intent(inout) :: phosphorusflux_vars

    type(elm_interface_bgc_datatype), intent(in):: elm_bgc_data

    integer :: fc, c, j, k
    character(len=256) :: subname = "update_bgc_flux_decomp_sourcesink"

    associate ( &
      decomp_cpools_sourcesink_vr  => carbonflux_vars%decomp_cpools_sourcesink    , &
      decomp_npools_sourcesink_vr  => nitrogenflux_vars%decomp_npools_sourcesink  , &
      decomp_ppools_sourcesink_vr  => phosphorusflux_vars%decomp_ppools_sourcesink  &
      )

      do fc = 1, num_soilc
        c = filter_soilc(fc)
            do k = 1, ndecomp_pools
                decomp_cpools_sourcesink_vr(c,:,k) = elm_bgc_data%decomp_cpools_sourcesink(c,:,k)
                decomp_npools_sourcesink_vr(c,:,k) = elm_bgc_data%decomp_npools_sourcesink(c,:,k)
                decomp_ppools_sourcesink_vr(c,:,k) = elm_bgc_data%decomp_ppools_sourcesink(c,:,k)
            end do
      end do
    end associate
  end subroutine update_bgc_flux_decomp_sourcesink
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
    subroutine update_bgc_flux_decomp_cascade(elm_bgc_data, &
           bounds, num_soilc, filter_soilc,                 &
           carbonflux_vars, nitrogenflux_vars,              &
           phosphorusflux_vars)

    use CNDecompCascadeConType, only : decomp_cascade_con

    implicit none

    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                     , intent(in)    :: filter_soilc(:)   ! filter for soil columns

    ! the following currently at column-level, but may be re-defined as grid-level if available (2019-10-28, fmy)
    type(column_carbon_flux)    , intent(inout) :: carbonflux_vars
    type(column_nitrogen_flux)  , intent(inout) :: nitrogenflux_vars
    type(column_phosphorus_flux), intent(inout) :: phosphorusflux_vars

    type(elm_interface_bgc_datatype), intent(in):: elm_bgc_data

    integer :: fc, c, j, k
    character(len=256) :: subname = "update_bgc_flux_decomp_cascade"

    associate ( &
     phr_vr                       => carbonflux_vars%phr_vr                                 , & ! Output: [real(r8) (:,:)   ]  potential HR (gC/m3/s)
     fphr                         => carbonflux_vars%fphr                                   , & ! Output: [real(r8) (:,:)   ]  fraction of potential SOM + LITTER heterotrophic
     !
     decomp_cascade_hr_vr         => carbonflux_vars%decomp_cascade_hr_vr                   , &
     !
     decomp_cascade_ctransfer_vr  => carbonflux_vars%decomp_cascade_ctransfer_vr            , &
     decomp_cascade_ntransfer_vr  => nitrogenflux_vars%decomp_cascade_ntransfer_vr          , &
     decomp_cascade_ptransfer_vr  => phosphorusflux_vars%decomp_cascade_ptransfer_vr        , &
     !
     decomp_cascade_sminn_flux_vr => nitrogenflux_vars%decomp_cascade_sminn_flux_vr         , & ! Output: [real(r8) (:,:,:) ]  vert-res mineral N flux for transition along decomposition cascade (gN/m3/s)
     decomp_cascade_sminp_flux_vr => phosphorusflux_vars%decomp_cascade_sminp_flux_vr       , & ! Output: [real(r8) (:,:,:) ]  vert-res mineral P flux for transition along decomposition cascade (gP/m3/s)
     !
     sminn_to_denit_decomp_cascade_vr => nitrogenflux_vars%sminn_to_denit_decomp_cascade_vr   & ! Output: [real(r8) (:,:,:) ]
     )

    do fc = 1, num_soilc
        c = filter_soilc(fc)
            phr_vr(c,:) = elm_bgc_data%phr_vr(c,:)
            fphr(c,:)   = elm_bgc_data%fphr(c,:)

            do k = 1, ndecomp_cascade_transitions
                decomp_cascade_hr_vr(c,:,k)         = elm_bgc_data%decomp_cascade_hr_vr(c,:,k)
                decomp_cascade_ctransfer_vr(c,:,k)  = elm_bgc_data%decomp_cascade_ctransfer_vr(c,:,k)
                decomp_cascade_ntransfer_vr(c,:,k)  = elm_bgc_data%decomp_cascade_ntransfer_vr(c,:,k)
                decomp_cascade_ptransfer_vr(c,:,k)  = elm_bgc_data%decomp_cascade_ptransfer_vr(c,:,k)

                decomp_cascade_sminn_flux_vr(c,:,k)     = elm_bgc_data%decomp_cascade_sminn_flux_vr(c,:,k)
                decomp_cascade_sminp_flux_vr(c,:,k)     = elm_bgc_data%decomp_cascade_sminp_flux_vr(c,:,k)
                sminn_to_denit_decomp_cascade_vr(c,:,k) = elm_bgc_data%sminn_to_denit_decomp_cascade_vr(c,:,k)
            end do
    end do
    end associate
    end subroutine update_bgc_flux_decomp_cascade
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
  subroutine update_bgc_flux_smin(elm_bgc_data,     &
           bounds, num_soilc, filter_soilc,         &
           cnstate_vars,                            &
           nitrogenflux_vars, phosphorusflux_vars)

    implicit none

    type(bounds_type)                   , intent(in)    :: bounds
    integer                             , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                             , intent(in)    :: filter_soilc(:)   ! filter for soil columns

    type(cnstate_type)                  , intent(inout) :: cnstate_vars

    ! the following currently at column-level, but may be re-defined as grid-level if available (2019-10-28, fmy)
    type(column_nitrogen_flux)          , intent(inout) :: nitrogenflux_vars
    type(column_phosphorus_flux)        , intent(inout) :: phosphorusflux_vars

    type(elm_interface_bgc_datatype)    , intent(in)    :: elm_bgc_data

    integer :: fc, c, j
    character(len=256) :: subname = "update_bgc_flux_smin"

    associate ( &
     fpg                          => cnstate_vars%fpg_col                        , & ! Output: [real(r8) (:)   ]  fraction of potential gpp (no units)
     fpi                          => cnstate_vars%fpi_col                        , & ! Output: [real(r8) (:)   ]  fraction of potential immobilization (no units)
     fpi_vr                       => cnstate_vars%fpi_vr_col                     , & ! Output: [real(r8) (:,:) ]  fraction of potential immobilization (no units)
     fpg_p                        => cnstate_vars%fpg_p_col                      , & ! Output: [real(r8) (:)   ]  fraction of potential gpp (no units)
     fpi_p                        => cnstate_vars%fpi_p_col                      , & ! Output: [real(r8) (:)   ]  fraction of potential immobilization (no units)
     fpi_p_vr                     => cnstate_vars%fpi_p_vr_col                   , & ! Output: [real(r8) (:,:) ]  fraction of potential immobilization (no units)
     !
     potential_immob              => nitrogenflux_vars%potential_immob           , & ! Output: [real(r8) (:)   ]
     actual_immob                 => nitrogenflux_vars%actual_immob              , & ! Output: [real(r8) (:)   ]
     sminn_to_plant               => nitrogenflux_vars%sminn_to_plant            , & ! Output: [real(r8) (:)   ]
     !
     sminn_to_plant_vr            => nitrogenflux_vars%sminn_to_plant_vr         , &
     smin_no3_to_plant_vr         => nitrogenflux_vars%smin_no3_to_plant_vr      , &
     smin_nh4_to_plant_vr         => nitrogenflux_vars%smin_nh4_to_plant_vr      , &
     potential_immob_vr           => nitrogenflux_vars%potential_immob_vr        , &
     actual_immob_vr              => nitrogenflux_vars%actual_immob_vr           , &
     actual_immob_no3_vr          => nitrogenflux_vars%actual_immob_no3_vr       , & ! Output: [real(r8) (:,:) ]
     actual_immob_nh4_vr          => nitrogenflux_vars%actual_immob_nh4_vr       , & ! Output: [real(r8) (:,:) ]
     gross_nmin_vr                => nitrogenflux_vars%gross_nmin_vr             , &
     net_nmin_vr                  => nitrogenflux_vars%net_nmin_vr               , & ! Output: [real(r8) (:,:)   ]
     !
     sminn_to_denit_excess_vr     => nitrogenflux_vars%sminn_to_denit_excess_vr  , & ! Output: [real(r8) (:,:) ]
     supplement_to_sminn_vr       => nitrogenflux_vars%supplement_to_sminn_vr    , & ! Output: [real(r8) (:,:) ]
     !
     no3_net_transport_vr         => nitrogenflux_vars%no3_net_transport_vr      , & ! Output: updated from PF, if coupled
     nh4_net_transport_vr         => nitrogenflux_vars%nh4_net_transport_vr      , & ! Output: updated from PF, if coupled
     !
     potential_immob_p            => phosphorusflux_vars%potential_immob_p       , & ! Output: [real(r8) (:)   ]
     actual_immob_p               => phosphorusflux_vars%actual_immob_p          , & ! Output: [real(r8) (:)   ]
     sminp_to_plant               => phosphorusflux_vars%sminp_to_plant          , & ! Output: [real(r8) (:)   ]
     !
     supplement_to_sminp_vr       => phosphorusflux_vars%supplement_to_sminp_vr  , & ! Output: [real(r8) (:,:) ]
     !
     sminp_to_plant_vr            => phosphorusflux_vars%sminp_to_plant_vr       , &
     potential_immob_p_vr         => phosphorusflux_vars%potential_immob_p_vr    , & ! Output: [real(r8) (:,:)   ]
     actual_immob_p_vr            => phosphorusflux_vars%actual_immob_p_vr       , &
     gross_pmin_vr                => phosphorusflux_vars%gross_pmin_vr           , & ! Output: [real(r8) (:,:)   ]
     net_pmin_vr                  => phosphorusflux_vars%net_pmin_vr               & ! Output: [real(r8) (:,:)   ]
     )

     do fc = 1, num_soilc
        c = filter_soilc(fc)

        fpg(c)              = elm_bgc_data%fpg(c)
        fpi(c)              = elm_bgc_data%fpi(c)
        fpg_p(c)            = elm_bgc_data%fpg_p(c)
        fpi_p(c)            = elm_bgc_data%fpi_p(c)

        potential_immob(c)  = elm_bgc_data%potential_immob(c)
        actual_immob(c)     = elm_bgc_data%actual_immob(c)
        sminn_to_plant(c)   = elm_bgc_data%sminn_to_plant(c)

        potential_immob_p(c)  = elm_bgc_data%potential_immob_p(c)
        actual_immob_p(c)     = elm_bgc_data%actual_immob_p(c)
        sminp_to_plant(c)     = elm_bgc_data%sminp_to_plant(c)

        fpi_vr(c,:)                 = elm_bgc_data%fpi_vr(c,:)
        fpi_p_vr(c,:)               = elm_bgc_data%fpi_p_vr(c,:)

        sminn_to_plant_vr(c,:)      = elm_bgc_data%sminn_to_plant_vr(c,:)
        smin_no3_to_plant_vr(c,:)   = elm_bgc_data%smin_no3_to_plant_vr(c,:)
        smin_nh4_to_plant_vr(c,:)   = elm_bgc_data%smin_nh4_to_plant_vr(c,:)

        potential_immob_vr(c,:)     = elm_bgc_data%potential_immob_vr(c,:)
        actual_immob_vr(c,:)        = elm_bgc_data%actual_immob_vr(c,:)
        actual_immob_no3_vr(c,:)    = elm_bgc_data%actual_immob_no3_vr(c,:)
        actual_immob_nh4_vr(c,:)    = elm_bgc_data%actual_immob_nh4_vr(c,:)
        gross_nmin_vr(c,:)          = elm_bgc_data%gross_nmin_vr(c,:)
        net_nmin_vr(c,:)            = elm_bgc_data%net_nmin_vr(c,:)

        sminn_to_denit_excess_vr(c,:)=elm_bgc_data%sminn_to_denit_excess_vr(c,:)
        supplement_to_sminn_vr(c,:) = elm_bgc_data%supplement_to_sminn_vr(c,:)

        no3_net_transport_vr(c,:)   = elm_bgc_data%no3_net_transport_vr(c,:)
        nh4_net_transport_vr(c,:)   = elm_bgc_data%nh4_net_transport_vr(c,:)

        supplement_to_sminp_vr(c,:) = elm_bgc_data%supplement_to_sminp_vr(c,:)

        sminp_to_plant_vr(c,:)      = elm_bgc_data%sminp_to_plant_vr(c,:)
        potential_immob_p_vr(c,:)   = elm_bgc_data%potential_immob_p_vr(c,:)
        actual_immob_p_vr(c,:)      = elm_bgc_data%actual_immob_p_vr(c,:)
        gross_pmin_vr(c,:)          = elm_bgc_data%gross_pmin_vr(c,:)
        net_pmin_vr(c,:)            = elm_bgc_data%net_pmin_vr(c,:)     !NOT available in PF

      end do
    end associate
  end subroutine update_bgc_flux_smin
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
  subroutine update_bgc_flux_nitdenit(elm_bgc_data,    &
           bounds, num_soilc, filter_soilc,            &
           nitrogenflux_vars)

    implicit none

    type(bounds_type)                   , intent(in)    :: bounds
    integer                             , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                             , intent(in)    :: filter_soilc(:)   ! filter for soil columns

    ! the following currently at column-level, but may be re-defined as grid-level if available (2019-10-28, fmy)
    type(column_nitrogen_flux)          , intent(inout) :: nitrogenflux_vars

    type(elm_interface_bgc_datatype)    , intent(in)    :: elm_bgc_data

    integer :: fc, c, j
    character(len=256) :: subname = "update_bgc_flux_nitdenit"

    associate ( &
         pot_f_nit_vr                 => nitrogenflux_vars%pot_f_nit_vr                  , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) potential soil nitrification flux
         pot_f_denit_vr               => nitrogenflux_vars%pot_f_denit_vr                , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) potential soil denitrification flux
         f_nit_vr                     => nitrogenflux_vars%f_nit_vr                      , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) soil nitrification flux
         f_denit_vr                   => nitrogenflux_vars%f_denit_vr                    , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) soil denitrification flux
         n2_n2o_ratio_denit_vr        => nitrogenflux_vars%n2_n2o_ratio_denit_vr         , & ! Output: [real(r8) (:,:) ]  ratio of N2 to N2O production by denitrification [gN/gN]
         f_n2o_denit_vr               => nitrogenflux_vars%f_n2o_denit_vr                , & ! Output: [real(r8) (:,:) ]  flux of N2O from denitrification [gN/m3/s]
         f_n2o_nit_vr                 => nitrogenflux_vars%f_n2o_nit_vr                    & ! Output: [real(r8) (:,:) ]  flux of N2O from nitrification [gN/m3/s]
    )

    do fc = 1, num_soilc
        c = filter_soilc(fc)
        pot_f_nit_vr(c,:)           = elm_bgc_data%pot_f_nit_vr(c,:)
        pot_f_denit_vr(c,:)         = elm_bgc_data%pot_f_denit_vr(c,:)
        f_nit_vr(c,:)               = elm_bgc_data%f_nit_vr(c,:)
        f_denit_vr(c,:)             = elm_bgc_data%f_denit_vr(c,:)
        n2_n2o_ratio_denit_vr(c,:)  = elm_bgc_data%n2_n2o_ratio_denit_vr(c,:)
        f_n2o_denit_vr(c,:)         = elm_bgc_data%f_n2o_denit_vr(c,:)
        f_n2o_nit_vr(c,:)           = elm_bgc_data%f_n2o_nit_vr(c,:)
    end do
    end associate
  end subroutine update_bgc_flux_nitdenit
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
  subroutine update_bgc_flux_gas_pf(elm_bgc_data,  &
     bounds, num_soilc, filter_soilc,              &
     carbonflux_vars, nitrogenflux_vars)

     ! PFLOTRAN gas fluxes
     implicit none

     type(bounds_type)                  , intent(in)    :: bounds
     integer                            , intent(in)    :: num_soilc       ! number of soil columns in filter
     integer                            , intent(in)    :: filter_soilc(:) ! filter for soil columns

    ! the following currently at column-level, but may be re-defined as grid-level if available (2019-10-28, fmy)
     type(column_carbon_flux)           , intent(inout) :: carbonflux_vars
     type(column_nitrogen_flux)         , intent(inout) :: nitrogenflux_vars

     type(elm_interface_bgc_datatype)   , intent(in)    :: elm_bgc_data

     !character(len=256) :: subname = "get_pf_bgc_gaslosses"

     integer  :: fc, c

!------------------------------------------------------------------------------------
    associate ( &
     hr_vr                        => carbonflux_vars%hr_vr              , &
     f_co2_soil_vr                => carbonflux_vars%f_co2_soil_vr      , &
     f_n2o_soil_vr                => nitrogenflux_vars%f_n2o_soil_vr    , &
     f_n2_soil_vr                 => nitrogenflux_vars%f_n2_soil_vr     , &
     f_ngas_decomp_vr             => nitrogenflux_vars%f_ngas_decomp_vr , &
     f_ngas_nitri_vr              => nitrogenflux_vars%f_ngas_nitri_vr  , &
     f_ngas_denit_vr              => nitrogenflux_vars%f_ngas_denit_vr    &
     )
! ------------------------------------------------------------------------
    do fc = 1,num_soilc
        c = filter_soilc(fc)
        f_co2_soil_vr(c,:)         = elm_bgc_data%f_co2_soil_vr(c,:)
        f_n2_soil_vr(c,:)          = elm_bgc_data%f_n2_soil_vr(c,:)
        f_n2o_soil_vr(c,:)         = elm_bgc_data%f_n2o_soil_vr(c,:)

        hr_vr(c,:)                 = elm_bgc_data%hr_vr(c,:)
        f_ngas_decomp_vr(c,:)      = elm_bgc_data%f_ngas_decomp_vr(c,:)
        f_ngas_nitri_vr(c,:)       = elm_bgc_data%f_ngas_nitri_vr(c,:)
        f_ngas_denit_vr(c,:)       = elm_bgc_data%f_ngas_denit_vr(c,:)

     enddo ! do c = begc, endc
!
    end associate
  end subroutine update_bgc_flux_gas_pf
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
  subroutine update_bgc_data_pf2elm(elm_bgc_data, bounds,         &
           num_soilc, filter_soilc,                               &
           num_soilp, filter_soilp,                               &
           cnstate_vars,                                          &
           ch4_vars)
    ! USES
    use elm_varctl          , only : use_pflotran, pf_cmode

    implicit none

    ! !ARGUMENTS:
    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                     , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                     , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                     , intent(in)    :: filter_soilp(:)   ! filter for soil patches

    type(cnstate_type)          , intent(inout) :: cnstate_vars
    type(ch4_type)              , intent(inout) :: ch4_vars

    type(elm_interface_bgc_datatype), intent(in):: elm_bgc_data

    !-----------------------------------------------------------------------

    character(len=256) :: subname = "update_bgc_data_pf2elm"


    if (pf_cmode) then
        ! bgc_state_decomp is updated in CLM
        ! by passing bgc_flux_decomp_sourcesink into SoilLittVertTransp
        call update_bgc_flux_decomp_sourcesink(elm_bgc_data,    &
                    bounds, num_soilc, filter_soilc,            &
                    col_cf, col_nf, col_pf)                       ! global vars

        call update_bgc_state_smin(elm_bgc_data,                &
                    bounds, num_soilc, filter_soilc,            &
                    col_ns, col_ps)                               ! global vars

        call update_bgc_flux_smin(elm_bgc_data,                 &
                    bounds, num_soilc, filter_soilc,            &
                    cnstate_vars,                               &
                    col_nf, col_pf)                               ! global vars

        call update_bgc_flux_gas_pf(elm_bgc_data,               &
                    bounds, num_soilc, filter_soilc,            &
                    col_cf, col_nf)                               ! global vars

    end if

  end subroutine update_bgc_data_pf2elm

!--------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------
! BEG of CLM-bgc through interface
!--------------------------------------------------------------------------------------
  ! !INTERFACE:
  subroutine elm_bgc_run(elm_interface_data, bounds,        &
                num_soilc, filter_soilc,                    &
                num_soilp, filter_soilp,                    &
                canopystate_vars, soilstate_vars,           &
                cnstate_vars, ch4_vars)

    ! USES:
    use SoilLittDecompMod          , only: SoilLittDecompAlloc

    ! ARGUMENTS:
    type(bounds_type)                   , intent(in)    :: bounds
    integer                             , intent(in)    :: num_soilc          ! number of soil columns in filter
    integer                             , intent(in)    :: filter_soilc(:)    ! filter for soil columns
    integer                             , intent(in)    :: num_soilp          ! number of soil patches in filter
    integer                             , intent(in)    :: filter_soilp(:)    ! filter for soil patches
    type(canopystate_type)              , intent(inout) :: canopystate_vars
    type(soilstate_type)                , intent(inout) :: soilstate_vars
    type(cnstate_type)                  , intent(inout) :: cnstate_vars
    type(ch4_type)                      , intent(inout) :: ch4_vars

    type(elm_interface_data_type)       , intent(inout) :: elm_interface_data
    real(r8) :: dt

    dt = dtime_mod
    !-------------------------------------------------------------
    ! STEP-2: (i) pass data from elm_bgc_data to SoilLittDecompAlloc
    call elm_bgc_get_data(elm_interface_data, bounds,       &
                num_soilc, filter_soilc,                    &
                canopystate_vars, soilstate_vars,           &
                cnstate_vars, ch4_vars,                     &
                col_es, col_ws,                             &
                col_cs, col_cf,                             &
                col_ns, col_nf,                             &
                col_ps, col_pf)

    ! STEP-2: (ii) run SoilLittDecompAlloc
    call SoilLittDecompAlloc (bounds, num_soilc, filter_soilc,    &
               num_soilp, filter_soilp,                     &
               canopystate_vars, soilstate_vars,            &
               cnstate_vars, ch4_vars, dt )

    ! STEP-2: (iii) update elm_bgc_data from SoilLittDecompAlloc
    call elm_bgc_update_data(elm_interface_data%bgc, bounds, &
                num_soilc, filter_soilc,                     &
                cnstate_vars,                                &
                col_cf, col_nf, col_pf)

  end subroutine elm_bgc_run
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
  ! !INTERFACE:
  ! pass data from elm_bgc_data to elm original data-types that used by SoilLittDecompAlloc
  subroutine elm_bgc_get_data(elm_interface_data,       &
            bounds, num_soilc, filter_soilc,            &
            canopystate_vars, soilstate_vars,           &
            cnstate_vars, ch4_vars,                     &
            energystate_vars, waterstate_vars,          &
            carbonstate_vars, carbonflux_vars,          &
            nitrogenstate_vars, nitrogenflux_vars,      &
            phosphorusstate_vars, phosphorusflux_vars)

    ! USES:

    ! ARGUMENTS:
    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: num_soilc          ! number of soil columns in filter
    integer                     , intent(in)    :: filter_soilc(:)    ! filter for soil columns
    type(canopystate_type)      , intent(inout) :: canopystate_vars
    type(soilstate_type)        , intent(inout) :: soilstate_vars
    type(cnstate_type)          , intent(inout) :: cnstate_vars
    type(ch4_type)              , intent(inout) :: ch4_vars

    ! the following currently at column-level, but may be re-defined as grid-level if available (2019-10-28, fmy)
    type(column_energy_state)     , intent(inout) :: energystate_vars
    type(column_water_state)      , intent(inout) :: waterstate_vars
    type(column_carbon_state)     , intent(inout) :: carbonstate_vars
    type(column_carbon_flux)      , intent(inout) :: carbonflux_vars
    type(column_nitrogen_state)   , intent(inout) :: nitrogenstate_vars
    type(column_nitrogen_flux)    , intent(inout) :: nitrogenflux_vars
    type(column_phosphorus_state) , intent(inout) :: phosphorusstate_vars
    type(column_phosphorus_flux)  , intent(inout) :: phosphorusflux_vars

    type(elm_interface_data_type)  , intent(in) :: elm_interface_data

    ! LOCAL VARIABLES:
    integer :: fc, c, j, k
    !-----------------------------------------------------------------------

    associate(&
        decomp_cpools_vr        => carbonstate_vars%decomp_cpools_vr        , & ! Input:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
        decomp_npools_vr        => nitrogenstate_vars%decomp_npools_vr      , & ! Input:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
        decomp_ppools_vr        => phosphorusstate_vars%decomp_ppools_vr    , & ! Input:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) P pools
        !
        smin_no3_vr             => nitrogenstate_vars%smin_no3_vr           , & ! [real(r8) (:,:) (gN/m3) vertically-resolved soil mineral NO3
        smin_nh4_vr             => nitrogenstate_vars%smin_nh4_vr           , & ! [real(r8) (:,:) (gN/m3) vertically-resolved soil mineral NH4
        smin_nh4sorb_vr         => nitrogenstate_vars%smin_nh4sorb_vr       , & ! [real(r8) (:,:) (gN/m3) vertically-resolved soil mineral NH4 absorbed
        !
        solutionp_vr            => phosphorusstate_vars%solutionp_vr        , & ! [real(r8) (:,:) (gP/m3) vertically-resolved soil solution P
        labilep_vr              => phosphorusstate_vars%labilep_vr          , & ! [real(r8) (:,:) (gP/m3) vertically-resolved soil labile mineral P
        secondp_vr              => phosphorusstate_vars%secondp_vr          , & ! [real(r8) (:,:) (gP/m3) vertically-resolved soil secondary mineralP
        sminp_vr                => phosphorusstate_vars%sminp_vr            , & ! [real(r8) (:,:) (gP/m3) vertically-resolved soil mineral P = solutionp + labilep + secondp
        occlp_vr                => phosphorusstate_vars%occlp_vr            , & ! [real(r8) (:,:) (gP/m3) vertically-resolved soil occluded mineral P
        primp_vr                => phosphorusstate_vars%primp_vr            , & ! [real(r8) (:,:) (gP/m3) vertically-resolved soil primary mineral P
        !
        plant_ndemand           => nitrogenflux_vars%plant_ndemand          , &
        plant_pdemand           => phosphorusflux_vars%plant_pdemand        , &
        !
        !alt_indx               => canopystate_vars%alt_indx                , & ! Input:  [integer  (:)     ]  current depth of thaw
        !
        watsat                  => soilstate_vars%watsat_col                , & ! Input:  [real(r8) (:,:)  ]  volumetric soil water at saturation (porosity) (nlevgrnd)
        bd                      => soilstate_vars%bd_col                    , & ! Input:  [real(r8) (:,:)  ]  bulk density of dry soil material [kg/m3]
        watfc                   => soilstate_vars%watfc_col                 , & ! Input:  [real(r8) (:,:)  ]  volumetric soil water at field capacity (nlevsoi)
        bsw                     => soilstate_vars%bsw_col                   , & ! Input:  [real(r8) (:,:)  ]  Clapp and Hornberger "b" (nlevgrnd)
        cellorg                 => soilstate_vars%cellorg_col               , & ! Input:  [real(r8) (:,:)  ]  column 3D org (kg/m3 organic matter) (nlevgrnd)
        sucsat                  => soilstate_vars%sucsat_col                , & ! Input:  [real(r8) (:,:)  ]  minimum soil suction (mm)
        !
        soilpsi                 => waterstate_vars%soilp                    , & ! Input:  [real(r8) (:,:)  ]  soil water potential in each soil layer (Pa)
        h2osoi_vol              => waterstate_vars%h2osoi_vol               , & ! Input:  [real(r8) (:,:)  ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevgrnd)
        h2osoi_liq              => waterstate_vars%h2osoi_liq               , & ! Input:  [real(r8) (:,:)  ]  liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)
        !
        t_soisno                => energystate_vars%t_soisno                , & ! Input:  [real(r8) (:,:)  ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
        !
        o2_decomp_depth_unsat   => ch4_vars%o2_decomp_depth_unsat_col       , & ! Input:  [real(r8) (:,:)  ]  O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
        o2_decomp_depth_sat     => ch4_vars%o2_decomp_depth_sat_col         , & ! Input:  [real(r8) (:,:)  ]  O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
        conc_o2_unsat           => ch4_vars%conc_o2_unsat_col               , & ! Input:  [real(r8) (:,:)  ]  O2 conc in each soil layer (mol/m3) (nlevsoi)
        conc_o2_sat             => ch4_vars%conc_o2_sat_col                 , & ! Input:  [real(r8) (:,:)  ]  O2 conc in each soil layer (mol/m3) (nlevsoi)
        o2stress_unsat          => ch4_vars%o2stress_unsat_col              , & ! Input:  [real(r8) (:,:)   ]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
        o2stress_sat            => ch4_vars%o2stress_sat_col                , & ! Input:  [real(r8) (:,:)   ]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
        !
        finundated              => ch4_vars%finundated_col                    & ! Input:  [real(r8) (:)     ]  fractional inundated area (excluding dedicated wetland columns)
    )

    ! soil properties & thermohydrology

    do fc = 1, num_soilc
        c = filter_soilc(fc)

        plant_ndemand(c)            = elm_interface_data%bgc%plant_ndemand(c)
        plant_pdemand(c)            = elm_interface_data%bgc%plant_pdemand(c)

        finundated(c)               = elm_interface_data%bgc%finundated(c)

        bd(c,:)                     = elm_interface_data%bd(c,:)
        watsat(c,:)                 = elm_interface_data%watsat(c,:)
        bsw(c,:)                    = elm_interface_data%bsw(c,:)
        sucsat(c,:)                 = elm_interface_data%sucsat(c,:)
        watfc(c,:)                  = elm_interface_data%watfc(c,:)
        cellorg(c,:)                = elm_interface_data%cellorg(c,:)

        soilpsi(c,:)                = elm_interface_data%th%soilpsi(c,:)
        h2osoi_vol(c,:)             = elm_interface_data%th%h2osoi_vol(c,:)
        h2osoi_liq(c,:)             = elm_interface_data%th%h2osoi_liq(c,:)

        t_soisno(c,:)               = elm_interface_data%th%t_soisno(c,:)

        o2stress_unsat(c,:)         = elm_interface_data%bgc%o2stress_unsat(c,:)
        o2stress_sat(c,:)           = elm_interface_data%bgc%o2stress_sat(c,:)
        o2_decomp_depth_unsat(c,:)  = elm_interface_data%bgc%o2_decomp_depth_unsat(c,:)
        conc_o2_unsat(c,:)          = elm_interface_data%bgc%conc_o2_unsat(c,:)
        o2_decomp_depth_sat(c,:)    = elm_interface_data%bgc%o2_decomp_depth_sat(c,:)
        conc_o2_sat(c,:)            = elm_interface_data%bgc%conc_o2_sat(c,:)

    end do

    !state variables
    do fc = 1, num_soilc
        c = filter_soilc(fc)
            do k = 1, ndecomp_pools
                decomp_cpools_vr(c,:,k) = elm_interface_data%bgc%decomp_cpools_vr(c,:,k)
                decomp_npools_vr(c,:,k) = elm_interface_data%bgc%decomp_npools_vr(c,:,k)
                decomp_ppools_vr(c,:,k) = elm_interface_data%bgc%decomp_ppools_vr(c,:,k)
            end do

            smin_no3_vr(c,:)        = elm_interface_data%bgc%smin_no3_vr(c,:)
            smin_nh4_vr(c,:)        = elm_interface_data%bgc%smin_nh4_vr(c,:)
            smin_nh4sorb_vr(c,:)    = elm_interface_data%bgc%smin_nh4sorb_vr(c,:)

            solutionp_vr(c,:)       = elm_interface_data%bgc%solutionp_vr(c,:)
            labilep_vr(c,:)         = elm_interface_data%bgc%labilep_vr(c,:)
            secondp_vr(c,:)         = elm_interface_data%bgc%secondp_vr(c,:)
            sminp_vr(c,:)           = elm_interface_data%bgc%sminp_vr(c,:)
            occlp_vr(c,:)           = elm_interface_data%bgc%occlp_vr(c,:)
            primp_vr(c,:)           = elm_interface_data%bgc%primp_vr(c,:)
    end do

    end associate
  end subroutine elm_bgc_get_data
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
  ! !INTERFACE:
  ! pass data from elm_bgc to elm_bgc_data
  subroutine elm_bgc_update_data(elm_bgc_data, bounds,  &
            num_soilc, filter_soilc,                    &
            cnstate_vars, carbonflux_vars,              &
            nitrogenflux_vars, phosphorusflux_vars)

    ! USES:

    ! ARGUMENTS:
    type(bounds_type)                   , intent(in)    :: bounds
    integer                             , intent(in)    :: num_soilc          ! number of soil columns in filter
    integer                             , intent(in)    :: filter_soilc(:)    ! filter for soil columns

    type(cnstate_type)                  , intent(in)    :: cnstate_vars

    ! the following currently at column-level, but may be re-defined as grid-level if available (2019-10-28, fmy)
    type(column_carbon_flux)            , intent(in)    :: carbonflux_vars
    type(column_nitrogen_flux)          , intent(in)    :: nitrogenflux_vars
    type(column_phosphorus_flux)        , intent(in)    :: phosphorusflux_vars

    type(elm_interface_bgc_datatype)    , intent(inout) :: elm_bgc_data

    ! LOCAL VARIABLES:
    integer :: fc, c, j, k

    !-----------------------------------------------------------------------

    associate(                                                                             &
         fpg                          => cnstate_vars%fpg_col                            , & ! Output: [real(r8) (:)   ]  fraction of potential gpp (no units)
         fpi                          => cnstate_vars%fpi_col                            , & ! Output: [real(r8) (:)   ]  fraction of potential immobilization (no units)
         fpi_vr                       => cnstate_vars%fpi_vr_col                         , & ! Output: [real(r8) (:,:) ]  fraction of potential immobilization (no units)
         fpg_p                        => cnstate_vars%fpg_p_col                          , & ! Output: [real(r8) (:)   ]  fraction of potential gpp (no units)
         fpi_p                        => cnstate_vars%fpi_p_col                          , & ! Output: [real(r8) (:)   ]  fraction of potential immobilization (no units)
         fpi_p_vr                     => cnstate_vars%fpi_p_vr_col                       , & ! Output: [real(r8) (:,:) ]  fraction of potential immobilization (no units)
         !
         potential_immob              => nitrogenflux_vars%potential_immob                          , & ! Output: [real(r8) (:)   ]
         actual_immob                 => nitrogenflux_vars%actual_immob                             , & ! Output: [real(r8) (:)   ]
         sminn_to_plant               => nitrogenflux_vars%sminn_to_plant                           , & ! Output: [real(r8) (:)   ]
         sminn_to_denit_excess_vr     => nitrogenflux_vars%sminn_to_denit_excess_vr                 , & ! Output: [real(r8) (:,:) ]
         pot_f_nit_vr                 => nitrogenflux_vars%pot_f_nit_vr                             , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) potential soil nitrification flux
         pot_f_denit_vr               => nitrogenflux_vars%pot_f_denit_vr                           , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) potential soil denitrification flux
         f_nit_vr                     => nitrogenflux_vars%f_nit_vr                                 , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) soil nitrification flux
         f_denit_vr                   => nitrogenflux_vars%f_denit_vr                               , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) soil denitrification flux
         actual_immob_no3_vr          => nitrogenflux_vars%actual_immob_no3_vr                      , & ! Output: [real(r8) (:,:) ]
         actual_immob_nh4_vr          => nitrogenflux_vars%actual_immob_nh4_vr                      , & ! Output: [real(r8) (:,:) ]
         smin_no3_to_plant_vr         => nitrogenflux_vars%smin_no3_to_plant_vr                     , & ! Output: [real(r8) (:,:) ]
         smin_nh4_to_plant_vr         => nitrogenflux_vars%smin_nh4_to_plant_vr                     , & ! Output: [real(r8) (:,:) ]
         n2_n2o_ratio_denit_vr        => nitrogenflux_vars%n2_n2o_ratio_denit_vr                    , & ! Output: [real(r8) (:,:) ]  ratio of N2 to N2O production by denitrification [gN/gN]
         f_n2o_denit_vr               => nitrogenflux_vars%f_n2o_denit_vr                           , & ! Output: [real(r8) (:,:) ]  flux of N2O from denitrification [gN/m3/s]
         f_n2o_nit_vr                 => nitrogenflux_vars%f_n2o_nit_vr                             , & ! Output: [real(r8) (:,:) ]  flux of N2O from nitrification [gN/m3/s]
         supplement_to_sminn_vr       => nitrogenflux_vars%supplement_to_sminn_vr                   , & ! Output: [real(r8) (:,:) ]
         sminn_to_plant_vr            => nitrogenflux_vars%sminn_to_plant_vr                        , & ! Output: [real(r8) (:,:) ]
         actual_immob_vr              => nitrogenflux_vars%actual_immob_vr                          , & ! Output: [real(r8) (:,:) ]
         !
         potential_immob_p            => phosphorusflux_vars%potential_immob_p                        , & ! Output: [real(r8) (:)   ]
         actual_immob_p               => phosphorusflux_vars%actual_immob_p                           , & ! Output: [real(r8) (:)   ]
         sminp_to_plant               => phosphorusflux_vars%sminp_to_plant                           , & ! Output: [real(r8) (:)   ]
         supplement_to_sminp_vr       => phosphorusflux_vars%supplement_to_sminp_vr                   , & ! Output: [real(r8) (:,:) ]
         sminp_to_plant_vr            => phosphorusflux_vars%sminp_to_plant_vr                        , & ! Output: [real(r8) (:,:) ]
         actual_immob_p_vr            => phosphorusflux_vars%actual_immob_p_vr                        , & ! Output: [real(r8) (:,:) ]
         !
         decomp_cascade_ntransfer_vr      => nitrogenflux_vars%decomp_cascade_ntransfer_vr       , & ! Output: [real(r8) (:,:,:) ]  vert-res transfer of N from donor to receiver pool along decomp. cascade (gN/m3/s)
         decomp_cascade_sminn_flux_vr     => nitrogenflux_vars%decomp_cascade_sminn_flux_vr      , & ! Output: [real(r8) (:,:,:) ]  vert-res mineral N flux for transition along decomposition cascade (gN/m3/s)
         potential_immob_vr               => nitrogenflux_vars%potential_immob_vr                , & ! Output: [real(r8) (:,:)   ]
         sminn_to_denit_decomp_cascade_vr => nitrogenflux_vars%sminn_to_denit_decomp_cascade_vr  , & ! Output: [real(r8) (:,:,:) ]
         gross_nmin_vr                    => nitrogenflux_vars%gross_nmin_vr                     , & ! Output: [real(r8) (:,:)   ]
         net_nmin_vr                      => nitrogenflux_vars%net_nmin_vr                       , & ! Output: [real(r8) (:,:)   ]
         gross_nmin                       => nitrogenflux_vars%gross_nmin                        , & ! Output: [real(r8) (:)     ]  gross rate of N mineralization (gN/m2/s)
         net_nmin                         => nitrogenflux_vars%net_nmin                          , & ! Output: [real(r8) (:)     ]  net rate of N mineralization (gN/m2/s)
         ! add phosphorus
         decomp_cascade_ptransfer_vr      => phosphorusflux_vars%decomp_cascade_ptransfer_vr     , & ! Output: [real(r8) (:,:,:) ]  vert-res transfer of P from donor to receiver pool along decomp. cascade (gP/m3/s)
         decomp_cascade_sminp_flux_vr     => phosphorusflux_vars%decomp_cascade_sminp_flux_vr    , & ! Output: [real(r8) (:,:,:) ]  vert-res mineral P flux for transition along decomposition cascade (gP/m3/s)
         potential_immob_p_vr             => phosphorusflux_vars%potential_immob_p_vr            , & ! Output: [real(r8) (:,:)   ]
         gross_pmin_vr                    => phosphorusflux_vars%gross_pmin_vr                   , & ! Output: [real(r8) (:,:)   ]
         net_pmin_vr                      => phosphorusflux_vars%net_pmin_vr                     , & ! Output: [real(r8) (:,:)   ]
         gross_pmin                       => phosphorusflux_vars%gross_pmin                      , & ! Output: [real(r8) (:)     ]  gross rate of P mineralization (gP/m2/s)
         net_pmin                         => phosphorusflux_vars%net_pmin                        , & ! Output: [real(r8) (:)     ]  net rate of P mineralization (gP/m2/s)
         !
         decomp_cascade_hr_vr             => carbonflux_vars%decomp_cascade_hr_vr             , & ! Output: [real(r8) (:,:,:) ]  vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
         decomp_cascade_ctransfer_vr      => carbonflux_vars%decomp_cascade_ctransfer_vr      , & ! Output: [real(r8) (:,:,:) ]  vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
         phr_vr                           => carbonflux_vars%phr_vr                           , & ! Output: [real(r8) (:,:)   ]  potential HR (gC/m3/s)
         fphr                             => carbonflux_vars%fphr                               & ! Output: [real(r8) (:,:)   ]  fraction of potential SOM + LITTER heterotrophic
         )

    !---------------------------------------------------------------------------
        do fc = 1,num_soilc
            c = filter_soilc(fc)
            elm_bgc_data%fpg(c)                                   = fpg(c)
            elm_bgc_data%fpi(c)                                   = fpi(c)
            elm_bgc_data%fpi_vr(c,:)                              = fpi_vr(c,:)
            elm_bgc_data%fpg_p(c)                                 = fpg_p(c)
            elm_bgc_data%fpi_p(c)                                 = fpi_p(c)
            elm_bgc_data%fpi_p_vr(c,:)                            = fpi_p_vr(c,:)
            elm_bgc_data%potential_immob(c)                       = potential_immob(c)
            elm_bgc_data%actual_immob(c)                          = actual_immob(c)
            elm_bgc_data%sminn_to_plant(c)                        = sminn_to_plant(c)
            elm_bgc_data%sminn_to_denit_excess_vr(c,:)            = sminn_to_denit_excess_vr(c,:)
            elm_bgc_data%pot_f_nit_vr(c,:)                        = pot_f_nit_vr(c,:)
            elm_bgc_data%pot_f_denit_vr(c,:)                      = pot_f_denit_vr(c,:)
            elm_bgc_data%f_nit_vr(c,:)                            = f_nit_vr(c,:)
            elm_bgc_data%f_denit_vr(c,:)                          = f_denit_vr(c,:)
            elm_bgc_data%actual_immob_no3_vr(c,:)                 = actual_immob_no3_vr(c,:)
            elm_bgc_data%actual_immob_nh4_vr(c,:)                 = actual_immob_nh4_vr(c,:)
            elm_bgc_data%smin_no3_to_plant_vr(c,:)                = smin_no3_to_plant_vr(c,:)
            elm_bgc_data%smin_nh4_to_plant_vr(c,:)                = smin_nh4_to_plant_vr(c,:)
            elm_bgc_data%n2_n2o_ratio_denit_vr(c,:)               = n2_n2o_ratio_denit_vr(c,:)
            elm_bgc_data%f_n2o_denit_vr(c,:)                      = f_n2o_denit_vr(c,:)
            elm_bgc_data%f_n2o_nit_vr(c,:)                        = f_n2o_nit_vr(c,:)
            elm_bgc_data%supplement_to_sminn_vr(c,:)              = supplement_to_sminn_vr(c,:)
            elm_bgc_data%sminn_to_plant_vr(c,:)                   = sminn_to_plant_vr(c,:)
            elm_bgc_data%potential_immob_vr(c,:)                  = potential_immob_vr(c,:)
            elm_bgc_data%actual_immob_vr(c,:)                     = actual_immob_vr(c,:)
            elm_bgc_data%potential_immob_p(c)                     = potential_immob_p(c)
            elm_bgc_data%actual_immob_p(c)                        = actual_immob_p(c)
            elm_bgc_data%sminp_to_plant(c)                        = sminp_to_plant(c)
            elm_bgc_data%supplement_to_sminp_vr(c,:)              = supplement_to_sminp_vr(c,:)
            elm_bgc_data%sminp_to_plant_vr(c,:)                   = sminp_to_plant_vr(c,:)
            elm_bgc_data%potential_immob_p_vr(c,:)                = potential_immob_p_vr(c,:)
            elm_bgc_data%actual_immob_p_vr(c,:)                   = actual_immob_p_vr(c,:)

            elm_bgc_data%decomp_cascade_ntransfer_vr(c,:,:)       = decomp_cascade_ntransfer_vr(c,:,:)
            elm_bgc_data%decomp_cascade_sminn_flux_vr(c,:,:)      = decomp_cascade_sminn_flux_vr(c,:,:)
            elm_bgc_data%sminn_to_denit_decomp_cascade_vr(c,:,:)  = sminn_to_denit_decomp_cascade_vr(c,:,:)
            elm_bgc_data%gross_nmin_vr(c,:)                       = gross_nmin_vr(c,:)
            elm_bgc_data%net_nmin_vr(c,:)                         = net_nmin_vr(c,:)

            ! phosphorus
            elm_bgc_data%decomp_cascade_ptransfer_vr(c,:,:)       = decomp_cascade_ptransfer_vr(c,:,:)
            elm_bgc_data%decomp_cascade_sminp_flux_vr(c,:,:)      = decomp_cascade_sminp_flux_vr(c,:,:)
            elm_bgc_data%potential_immob_p_vr(c,:)                = potential_immob_p_vr(c,:)
            elm_bgc_data%gross_pmin_vr(c,:)                       = gross_pmin_vr(c,:)
            elm_bgc_data%net_pmin_vr(c,:)                         = net_pmin_vr(c,:)

            elm_bgc_data%decomp_cascade_hr_vr(c,:,:)              = decomp_cascade_hr_vr(c,:,:)
            elm_bgc_data%decomp_cascade_ctransfer_vr(c,:,:)       = decomp_cascade_ctransfer_vr(c,:,:)

            elm_bgc_data%phr_vr(c,:)                              = phr_vr(c,:)
            elm_bgc_data%fphr(c,:)                                = fphr(c,:)
        end do


    end associate
  end subroutine elm_bgc_update_data
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
  subroutine update_bgc_data_elm2elm(elm_bgc_data,bounds,         &
           num_soilc, filter_soilc,                               &
           num_soilp, filter_soilp,                               &
           cnstate_vars,                                          &
           ch4_vars)
    ! USES

    implicit none

    ! !ARGUMENTS:
    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                     , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                     , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                     , intent(in)    :: filter_soilp(:)   ! filter for soil patches

    type(cnstate_type)          , intent(inout) :: cnstate_vars
    type(ch4_type)              , intent(inout) :: ch4_vars

    type(elm_interface_bgc_datatype), intent(in):: elm_bgc_data

    !-----------------------------------------------------------------------

    character(len=256) :: subname = "update_bgc_data_elm2elm"

    ! bgc_state_decomp is updated in CLM
    ! by passing bgc_flux_decomp_sourcesink into SoilLittVertTransp
    call update_bgc_flux_decomp_cascade(elm_bgc_data,   &
                bounds, num_soilc, filter_soilc,        &
                col_cf, col_nf, col_pf)

    call update_bgc_flux_smin(elm_bgc_data,             &
                bounds, num_soilc, filter_soilc,        &
                cnstate_vars,                           &
                col_nf, col_pf)

    call update_bgc_flux_nitdenit(elm_bgc_data,         &
                bounds, num_soilc, filter_soilc,        &
                col_nf)

  end subroutine update_bgc_data_elm2elm
!--------------------------------------------------------------------------------------
! END of CLM-bgc through interface
!--------------------------------------------------------------------------------------


end module elm_interface_funcsMod

