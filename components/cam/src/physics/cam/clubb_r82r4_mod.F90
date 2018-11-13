module clubb_r82r4_mod

implicit none

public ::  &

clubb_r82r4_core

private ! Default Scope

contains

subroutine clubb_r82r4_core &
             ( hydromet_dim, dtime, fcor, sfc_elevation, &          ! intent(in)
               thlm_forcing, rtm_forcing, um_forcing, vm_forcing, & ! intent(in)
               sclrm_forcing, edsclrm_forcing, wprtp_forcing, &     ! intent(in)
               wpthlp_forcing, rtp2_forcing, thlp2_forcing, &       ! intent(in)
               rtpthlp_forcing, wm_zm, wm_zt, &                     ! intent(in)
               wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, &         ! intent(in)
               wpsclrp_sfc, wpedsclrp_sfc, &                        ! intent(in)
               p_in_Pa, rho_zm, rho_in, exner, &                       ! intent(in)
               rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &             ! intent(in)
               invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, hydromet, &   ! intent(in)
               rfrzm, radf, &                                       ! intent(in)
               wphydrometp, wp2hmp, rtphmp_zt, thlphmp_zt, &        ! intent(in)
               host_dx, host_dy, &                                  ! intent(in)
               um_in, vm_in, upwp_in, vpwp_in, up2_in, vp2_in, &    ! intent(in)
               thlm_in, rtm_in, rvm_in, wprtp_in, wpthlp_in, &      ! intent(in)
               wp2_in, wp3_in, rtp2_in, thlp2_in, rtpthlp_in, &     ! intent(in)
               sclrm,   &                                           ! intent(in) 
               sclrp2, sclrprtp, sclrpthlp, &                       ! intent(in)
               wpsclrp, edsclr_in,  &                                 ! intent(in)
               rcm_out, wprcp_out, cloud_frac_out, ice_supersat_frac, &         ! intent(in)
               rcm_in_layer_out, cloud_cover_out, &                         ! intent(in)
               khzm_out, khzt_out, &                                        ! intent(in)
               qclvar_out, thlprcp_out, &                               ! intent(in)
               dtime_r4, fcor_r4, sfc_elevation_r4, &                           ! intent(inout)
               thlm_forcing_r4, rtm_forcing_r4, um_forcing_r4, vm_forcing_r4, & ! intent(inout)
               sclrm_forcing_r4, edsclrm_forcing_r4, wprtp_forcing_r4, &        ! intent(inout)
               wpthlp_forcing_r4, rtp2_forcing_r4, thlp2_forcing_r4, &          ! intent(inout)
               rtpthlp_forcing_r4, wm_zm_r4, wm_zt_r4, &                           ! intent(inout)
               wpthlp_sfc_r4, wprtp_sfc_r4, upwp_sfc_r4, vpwp_sfc_r4, &         ! intent(inout)
               wpsclrp_sfc_r4, wpedsclrp_sfc_r4, &                              ! intent(inout)
               p_in_Pa_r4, rho_zm_r4, rho_in_r4, exner_r4, &                       ! intent(inout)
               rho_ds_zm_r4, rho_ds_zt_r4, invrs_rho_ds_zm_r4, &                ! intent(inout)
               invrs_rho_ds_zt_r4, thv_ds_zm_r4, thv_ds_zt_r4, hydromet_r4, &   ! intent(inout)
               rfrzm_r4, radf_r4, &                                             ! intent(inout)
               wphydrometp_r4, wp2hmp_r4, rtphmp_zt_r4, thlphmp_zt_r4, &        ! intent(inout)
               host_dx_r4, host_dy_r4, &                                  ! intent(in)
               um_in_r4, vm_in_r4, upwp_in_r4, vpwp_in_r4, up2_in_r4, vp2_in_r4, &          ! intent(inout)
               thlm_in_r4, rtm_in_r4, rvm_in_r4, wprtp_in_r4, wpthlp_in_r4, &               ! intent(inout)
               wp2_in_r4, wp3_in_r4, rtp2_in_r4, thlp2_in_r4, rtpthlp_in_r4, &              ! intent(inout)
               sclrm_r4,   &
               sclrp2_r4, sclrprtp_r4, sclrpthlp_r4, &                             ! intent(inout)
               wpsclrp_r4, edsclr_in_r4, &                                         ! intent(inout)
               rcm_out_r4, wprcp_out_r4, cloud_frac_out_r4, ice_supersat_frac_r4, &         ! intent(inout)
               rcm_in_layer_out_r4, cloud_cover_out_r4, &                               ! intent(inout)
               khzm_out_r4, khzt_out_r4, &                                              ! intent(inout)
               qclvar_out_r4, thlprcp_out_r4 &                                         ! intent(inout)
                        )

   use shr_kind_mod,              only: r8=>shr_kind_r8
   use clubb_precision,           only: core_rknd                     ! Constant(s)
   use parameter_indices,         only: nparams
   use ppgrid,                    only: pver, pverp, pcols
   use pdf_parameter_module,      only: pdf_parameter ! Type
   use parameters_model,          only: sclr_dim, edsclr_dim  ! Variable(s)

   implicit none
   integer :: hydromet_dim                       ! The hydromet array in SAM-CLUBB is currently 0 elements

   real(r8), intent(in) :: host_dx
   real(r8), intent(in) :: host_dy
   real(r8), intent(in) :: dtime                            ! CLUBB time step                              [s]   
   real(r8), intent(in) :: fcor                             ! Coriolis forcing                              [s^-1]
   real(r8), intent(in) :: sfc_elevation                    ! Elevation of ground                           [m AMSL]
   real(r8), intent(in) :: thlm_forcing(pverp)              ! theta_l forcing (thermodynamic levels)        [K/s]
   real(r8), intent(in) :: rtm_forcing(pverp)               ! r_t forcing (thermodynamic levels)            [(kg/kg)/s]                              
   real(r8), intent(in) :: um_forcing(pverp)                ! u wind forcing (thermodynamic levels)         [m/s/s]
   real(r8), intent(in) :: vm_forcing(pverp)                ! v wind forcing (thermodynamic levels)         [m/s/s]
   real(r8), intent(in) :: sclrm_forcing(pverp,sclr_dim)    ! Passive scalar forcing                        [{units vary}/s]
   real(r8), intent(in) :: edsclrm_forcing(pverp,edsclr_dim)! Eddy passive scalar forcing                   [{units vary}/s]
   real(r8), intent(in) :: wprtp_forcing(pverp)
   real(r8), intent(in) :: wpthlp_forcing(pverp)
   real(r8), intent(in) :: rtp2_forcing(pverp)
   real(r8), intent(in) :: thlp2_forcing(pverp)
   real(r8), intent(in) :: rtpthlp_forcing(pverp)
   real(r8), intent(in) :: wm_zm(pverp)                     ! w mean wind component on momentum levels      [m/s]
   real(r8), intent(in) :: wm_zt(pverp)                     ! w mean wind component on thermo. levels       [m/s]
   real(r8), intent(in) :: wpthlp_sfc                       ! w' theta_l' at surface                        [(m K)/s]
   real(r8), intent(in) :: wprtp_sfc                        ! w' r_t' at surface                            [(kg m)/( kg s)]
   real(r8), intent(in) :: upwp_sfc                         ! u'w' at surface                               [m^2/s^2]
   real(r8), intent(in) :: vpwp_sfc                         ! v'w' at surface                               [m^2/s^2]   
   real(r8), intent(in) :: wpsclrp_sfc(sclr_dim)            ! Scalar flux at surface                        [{units vary} m/s]
   real(r8), intent(in) :: wpedsclrp_sfc(edsclr_dim)        ! Eddy-scalar flux at surface                   [{units vary} m/s]
   real(r8), intent(in) :: p_in_Pa(pverp)                   ! Air pressure (thermodynamic levels)           [Pa]
   real(r8), intent(in) :: exner(pverp)                     ! Exner function (thermodynamic levels)         [-]
   real(r8), intent(in) :: rho_ds_zm(pverp)                 ! Dry, static density on momentum levels        [kg/m^3]
   real(r8), intent(in) :: rho_ds_zt(pverp)                 ! Dry, static density on thermodynamic levels   [kg/m^3]
   real(r8), intent(in) :: rho_zm(pverp)                 ! Dry, static density on momentum levels        [kg/m^3]
   real(r8), intent(in) :: rho_in(pverp)                 ! Dry, static density on thermodynamic levels   [kg/m^3]
   real(r8), intent(in) :: invrs_rho_ds_zm(pverp)           ! Inv. dry, static density on momentum levels   [m^3/kg]
   real(r8), intent(in) :: invrs_rho_ds_zt(pverp)           ! Inv. dry, static density on thermo. levels    [m^3/kg]
   real(r8), intent(in) :: thv_ds_zm(pverp)                 ! Dry, base-state theta_v on momentum levels    [K]
   real(r8), intent(in) :: thv_ds_zt(pverp)                 ! Dry, base-state theta_v on thermo. levels     [K]
   real(r8), intent(in) :: hydromet(pverp,hydromet_dim)
   real(r8), intent(in) :: rfrzm(pverp)
   real(r8), intent(in) :: radf(pverp)
   real(r8), intent(in) :: wphydrometp(pverp,hydromet_dim)
   real(r8), intent(in) :: wp2hmp(pverp,hydromet_dim)
   real(r8), intent(in) :: rtphmp_zt(pverp,hydromet_dim)
   real(r8), intent(in) :: thlphmp_zt (pverp,hydromet_dim)
   real(r8), intent(in) :: um_in(pverp)                     ! meridional wind                              [m/s]
   real(r8), intent(in) :: vm_in(pverp)                     ! zonal wind        
   real(r8), intent(in) :: up2_in(pverp)                    ! meridional wind variance                     [m^2/s^2]
   real(r8), intent(in) :: vp2_in(pverp)                    ! zonal wind variance                          [m^2/s^2]
   real(r8), intent(in) :: upwp_in(pverp)                   ! meridional wind flux                         [m^2/s^2]
   real(r8), intent(in) :: vpwp_in(pverp)                   ! zonal wind flux                              [m^2/s^2]
   real(r8), intent(in) :: thlm_in(pverp)                   ! liquid water potential temperature (thetal)  [K]
   real(r8), intent(in) :: rtm_in(pverp)                    ! total water mixing ratio                     [kg/kg]
   real(r8), intent(in) :: rvm_in(pverp)                    !  water vapor mixing ratio                     [kg/kg]
   real(r8), intent(in) :: wprtp_in(pverp)                  ! turbulent flux of total water                [kg/kg m/s]
   real(r8), intent(in) :: wpthlp_in(pverp)                 ! turbulent flux of thetal                     [K m/s]
   real(r8), intent(in) :: wp2_in(pverp)                    ! vertical velocity variance (CLUBB)           [m^2/s^2]
   real(r8), intent(in) :: wp3_in(pverp)                    ! third moment vertical velocity               [m^3/s^3]
   real(r8), intent(in) :: rtp2_in(pverp)                   ! total water variance                         [kg^2/k^2]
   real(r8), intent(in) :: thlp2_in(pverp)                  ! thetal variance                              [K^2]
   real(r8), intent(in) :: rtpthlp_in(pverp)                ! covariance of thetal and qt                  [kg/kg K]
   real(r8), intent(in) :: sclrm(pverp,sclr_dim)            ! Passive scalar mean (thermo. levels)          [units vary]
   real(r8), intent(in) :: sclrp2(pverp,sclr_dim)           ! sclr'^2 (momentum levels)                     [{units vary}^2]
   real(r8), intent(in) :: sclrprtp(pverp,sclr_dim)         ! sclr'rt' (momentum levels)                    [{units vary} (kg/kg)]
   real(r8), intent(in) :: sclrpthlp(pverp,sclr_dim)        ! sclr'thlp' (momentum levels)                  [{units vary} (K)]
   real(r8), intent(in) :: wpsclrp(pverp,sclr_dim)          ! w'sclr' (momentum levels)                     [{units vary} m/s]
   real(r8), intent(in) :: edsclr_in(pverp,edsclr_dim)      ! Scalars to be diffused through CLUBB         [units vary]
   real(r8), intent(in) :: rcm_out(pverp)                   ! CLUBB output of liquid water mixing ratio     [kg/kg]
   real(r8), intent(in) :: wprcp_out(pverp)                 ! CLUBB output of flux of liquid water          [kg/kg m/s]
   real(r8), intent(in) :: cloud_frac_out(pverp)            ! CLUBB output of cloud fraction                [fraction]
   real(r8), intent(in) :: ice_supersat_frac(pverp)
   real(r8), intent(in) :: rcm_in_layer_out(pverp)          ! CLUBB output of in-cloud liq. wat. mix. ratio [kg/kg]
   real(r8), intent(in) :: cloud_cover_out(pverp)           ! CLUBB output of in-cloud cloud fraction       [fraction]
   real(r8), intent(in) :: khzm_out(pverp)                  ! eddy diffusivity on momentum grids            [m^2/s]
   real(r8), intent(in) :: khzt_out(pverp)                  ! eddy diffusivity on thermo grids              [m^2/s]
   real(r8), intent(in) :: qclvar_out(pverp)                ! cloud water variance                          [kg^2/kg^2]
   real(r8), intent(in) :: thlprcp_out(pverp)


   real( kind = core_rknd ), intent(inout) :: dtime_r4 ! time-step size 
   real( kind = core_rknd ), intent(inout) :: fcor_r4            ! Coriolis forcing             [s^-1]
   real( kind = core_rknd ), intent(inout) :: sfc_elevation_r4     ! Elevation of ground level    [m AMSL]

   real( kind = core_rknd ), intent(inout), dimension(pverp) ::  &
      thlm_forcing_r4,    & ! theta_l forcing (thermodynamic levels)    [K/s]
      rtm_forcing_r4,     & ! r_t forcing (thermodynamic levels)        [(kg/kg)/s]
      um_forcing_r4,      & ! u wind forcing (thermodynamic levels)     [m/s/s]
      vm_forcing_r4,      & ! v wind forcing (thermodynamic levels)     [m/s/s]
      wprtp_forcing_r4,   & ! <w'r_t'> forcing (momentum levels)    [m*K/s^2]
      wpthlp_forcing_r4,  & ! <w'th_l'> forcing (momentum levels)   [m*(kg/kg)/s^2]
      rtp2_forcing_r4,    & ! <r_t'^2> forcing (momentum levels)    [(kg/kg)^2/s]
      thlp2_forcing_r4,   & ! <th_l'^2> forcing (momentum levels)   [K^2/s]
      rtpthlp_forcing_r4, & ! <r_t'th_l'> forcing (momentum levels) [K*(kg/kg)/s]
      wm_zm_r4,           & ! w mean wind component on momentum levels  [m/s]
      wm_zt_r4,           & ! w mean wind component on thermo. levels   [m/s]
      p_in_Pa_r4,         & ! Air pressure (thermodynamic levels)       [Pa]
      rho_zm_r4,          & ! Air density on momentum levels            [kg/m^3]
      rho_in_r4,          & ! Air density on thermodynamic levels       [kg/m^3]
      exner_r4,           & ! Exner function (thermodynamic levels)     [-]
      rho_ds_zm_r4,       & ! Dry, static density on momentum levels    [kg/m^3]
      rho_ds_zt_r4,       & ! Dry, static density on thermo. levels     [kg/m^3]
      invrs_rho_ds_zm_r4, & ! Inv. dry, static density @ momentum levs. [m^3/kg]
      invrs_rho_ds_zt_r4, & ! Inv. dry, static density @ thermo. levs.  [m^3/kg]
      thv_ds_zm_r4,       & ! Dry, base-state theta_v on momentum levs. [K]
      thv_ds_zt_r4,       & ! Dry, base-state theta_v on thermo. levs.  [K]
      rfrzm_r4              ! Total ice-phase water mixing ratio        [kg/kg]


    real( kind = core_rknd ), intent(inout), dimension(pverp,hydromet_dim) :: hydromet_r4        ! Collection of hydrometeors                [units vary]
    real( kind = core_rknd ), intent(inout), dimension(pverp) :: radf_r4        ! Buoyancy production at the CL top due to LW radiative cooling [m^2/s^3]

    real( kind = core_rknd ), intent(inout), dimension(pverp,hydromet_dim) :: &
      wphydrometp_r4, & ! Covariance of w and a hydrometeor      [(m/s) <hm units>]
      wp2hmp_r4,      & ! Third-order moment:  < w'^2 hm' >    [(m/s)^2 <hm units>]
      rtphmp_zt_r4,   & ! Covariance of rt and hm (on t-levs.) [(kg/kg) <hm units>]
      thlphmp_zt_r4     ! Covariance of thl and hm (on t-levs.)      [K <hm units>]

    real( kind = core_rknd ), intent(inout) ::  &
      wpthlp_sfc_r4,   & ! w' theta_l' at surface   [(m K)/s]
      wprtp_sfc_r4,    & ! w' r_t' at surface       [(kg m)/( kg s)]
      upwp_sfc_r4,     & ! u'w' at surface          [m^2/s^2]
      vpwp_sfc_r4        ! v'w' at surface          [m^2/s^2]

    ! Passive scalar variables
    real( kind = core_rknd ), intent(inout), dimension(pverp,sclr_dim) :: &
      sclrm_forcing_r4    ! Passive scalar forcing         [{units vary}/s]

    real( kind = core_rknd ), intent(inout), dimension(sclr_dim) ::  &
      wpsclrp_sfc_r4      ! Scalar flux at surface         [{units vary} m/s]

    ! Eddy passive scalar variables
    real( kind = core_rknd ), intent(inout), dimension(pverp,edsclr_dim) :: &
      edsclrm_forcing_r4  ! Eddy passive scalar forcing    [{units vary}/s]
                          ! edsclr_dim  is defined at the beginning of the club_intr.F90, 
                          ! also it is setup in the setup_clubb_core() subroutine. 

    real( kind = core_rknd ), intent(inout), dimension(edsclr_dim) ::  &
      wpedsclrp_sfc_r4    ! Eddy-Scalar flux at surface    [{units vary} m/s]

    !!! Input/Output Variables
    ! These are prognostic or are planned to be in the future
    real( kind = core_rknd ), intent(inout), dimension(pverp) ::  &
      um_in_r4,      & ! u mean wind component (thermodynamic levels)   [m/s]
      upwp_in_r4,    & ! u'w' (momentum levels)                         [m^2/s^2]
      vm_in_r4,      & ! v mean wind component (thermodynamic levels)   [m/s]
      vpwp_in_r4,    & ! v'w' (momentum levels)                         [m^2/s^2]
      up2_in_r4,     & ! u'^2 (momentum levels)                         [m^2/s^2]
      vp2_in_r4,     & ! v'^2 (momentum levels)                         [m^2/s^2]
      rtm_in_r4,     & ! total water mixing ratio, r_t (thermo. levels) [kg/kg]
      rvm_in_r4,     & ! water vapor mixing ratio, r_t (thermo. levels) [kg/kg]
      wprtp_in_r4,   & ! w' r_t' (momentum levels)                      [(kg/kg) m/s]
      thlm_in_r4,    & ! liq. water pot. temp., th_l (thermo. levels)   [K]
      wpthlp_in_r4,  & ! w' th_l' (momentum levels)                     [(m/s) K]
      rtp2_in_r4,    & ! r_t'^2 (momentum levels)                       [(kg/kg)^2]
      thlp2_in_r4,   & ! th_l'^2 (momentum levels)                      [K^2]
      rtpthlp_in_r4, & ! r_t' th_l' (momentum levels)                   [(kg/kg) K]
      wp2_in_r4,     & ! w'^2 (momentum levels)                         [m^2/s^2]
      wp3_in_r4        ! w'^3 (thermodynamic levels)                    [m^3/s^3]

    ! Passive scalar variables
    real( kind = core_rknd ), intent(inout), dimension(pverp,sclr_dim) :: &
      sclrm_r4,     & ! Passive scalar mean (thermo. levels) [units vary]
      wpsclrp_r4,   & ! w'sclr' (momentum levels)            [{units vary} m/s]
      sclrp2_r4,    & ! sclr'^2 (momentum levels)            [{units vary}^2]
      sclrprtp_r4,  & ! sclr'rt' (momentum levels)           [{units vary} (kg/kg)]
      sclrpthlp_r4    ! sclr'thl' (momentum levels)          [{units vary} K]

    ! Eddy passive scalar variable
    real( kind = core_rknd ), intent(inout), dimension(pverp,edsclr_dim) :: &
      edsclr_in_r4   ! Eddy passive scalar mean (thermo. levels)   [units vary]

    ! Host model horizontal grid spacing, if part of host model.
    real( kind = core_rknd ), intent(inout) :: &
      host_dx_r4,  & ! East-West horizontal grid spacing     [m]
      host_dy_r4     ! North-South horizontal grid spacing   [m]

    real( kind = core_rknd ), intent(inout), dimension(pverp) ::  &
      rcm_out_r4,          & ! cloud water mixing ratio, r_c (thermo. levels)  [kg/kg]
      rcm_in_layer_out_r4, & ! rcm in cloud layer                              [kg/kg]
      cloud_cover_out_r4     ! cloud cover    

    real( kind = core_rknd ), intent(inout), dimension(pverp) :: &
      khzt_out_r4, &       ! eddy diffusivity on thermo levels
      khzm_out_r4, &       ! eddy diffusivity on momentum levels
      thlprcp_out_r4

    real( kind = core_rknd), intent(inout), dimension(pverp) :: &
      qclvar_out_r4        ! cloud water varian

    ! Variables that need to be output for use in host models
    real( kind = core_rknd ), intent(inout), dimension(pverp) ::  &
      wprcp_out_r4,            & ! w'r_c' (momentum levels)                  [(kg/kg) m/s]
      cloud_frac_out_r4,       & ! cloud fraction (thermodynamic levels)     [-]
      ice_supersat_frac_r4   ! ice cloud fraction (thermodynamic levels) [-]


    !!! SZhang and HWan, reduce precision for the CLUBB call, convert the r8 variables to r4 precision
      dtime_r4             = real(dtime, kind = core_rknd)
      fcor_r4              = real(fcor, kind = core_rknd)
      sfc_elevation_r4     = real(sfc_elevation, kind = core_rknd)
      thlm_forcing_r4      = real(thlm_forcing, kind = core_rknd)
      rtm_forcing_r4       = real(rtm_forcing, kind = core_rknd)
      um_forcing_r4        = real(um_forcing, kind = core_rknd)
      vm_forcing_r4        = real(vm_forcing, kind = core_rknd)
      sclrm_forcing_r4     = real(sclrm_forcing, kind = core_rknd)
      edsclrm_forcing_r4   = real(edsclrm_forcing, kind = core_rknd)
      wprtp_forcing_r4     = real(wprtp_forcing, kind = core_rknd)
      wpthlp_forcing_r4    = real(wpthlp_forcing, kind = core_rknd)
      rtp2_forcing_r4      = real(rtp2_forcing, kind = core_rknd)
      thlp2_forcing_r4     = real(thlp2_forcing, kind = core_rknd)
      rtpthlp_forcing_r4   = real(rtpthlp_forcing, kind = core_rknd)
      wm_zm_r4             = real(wm_zm, kind = core_rknd)
      wm_zt_r4             = real(wm_zt, kind = core_rknd)
      wpthlp_sfc_r4        = real(wpthlp_sfc, kind = core_rknd)
      wprtp_sfc_r4         = real(wprtp_sfc, kind = core_rknd)
      upwp_sfc_r4          = real(upwp_sfc, kind = core_rknd)
      vpwp_sfc_r4          = real(vpwp_sfc, kind = core_rknd)
      wpsclrp_sfc_r4       = real(wpsclrp_sfc, kind = core_rknd)
      wpedsclrp_sfc_r4     = real(wpedsclrp_sfc, kind = core_rknd)
      p_in_Pa_r4           = real(p_in_Pa, kind = core_rknd)
      rho_zm_r4            = real(rho_zm, kind = core_rknd)
      rho_in_r4            = real(rho_in, kind = core_rknd)
      exner_r4             = real(exner, kind = core_rknd)
      rho_ds_zm_r4         = real(rho_ds_zm, kind = core_rknd)
      rho_ds_zt_r4         = real(rho_ds_zt, kind = core_rknd)
      invrs_rho_ds_zm_r4   = real(invrs_rho_ds_zm, kind = core_rknd)
      invrs_rho_ds_zt_r4   = real(invrs_rho_ds_zt, kind = core_rknd)
      thv_ds_zm_r4         = real(thv_ds_zm, kind = core_rknd)
      thv_ds_zt_r4         = real(thv_ds_zt, kind = core_rknd)
      hydromet_r4          = real(hydromet, kind = core_rknd)
      rfrzm_r4             = real(rfrzm, kind = core_rknd)
      radf_r4              = real(radf, kind = core_rknd)
      wphydrometp_r4       = real(wphydrometp, kind = core_rknd)
      wp2hmp_r4            = real(wp2hmp, kind = core_rknd)
      rtphmp_zt_r4         = real(rtphmp_zt, kind = core_rknd)
      thlphmp_zt_r4        = real(thlphmp_zt, kind = core_rknd)
      host_dx_r4           = real(host_dx, kind = core_rknd)
      host_dy_r4           = real(host_dy, kind = core_rknd)
      um_in_r4             = real(um_in, kind = core_rknd)
      vm_in_r4             = real(vm_in, kind = core_rknd)
      upwp_in_r4           = real(upwp_in, kind = core_rknd)
      vpwp_in_r4           = real(vpwp_in, kind = core_rknd)
      up2_in_r4            = real(up2_in, kind = core_rknd)
      vp2_in_r4            = real(vp2_in, kind = core_rknd)
      thlm_in_r4           = real(thlm_in, kind = core_rknd)
      rtm_in_r4            = real(rtm_in, kind = core_rknd)
      rvm_in_r4            = real(rvm_in, kind = core_rknd)
      wprtp_in_r4          = real(wprtp_in, kind = core_rknd)
      wpthlp_in_r4         = real(wpthlp_in, kind = core_rknd)
      wp2_in_r4            = real(wp2_in, kind = core_rknd)
      wp3_in_r4            = real(wp3_in, kind = core_rknd)
      rtp2_in_r4           = real(rtp2_in, kind = core_rknd)
      thlp2_in_r4          = real(thlp2_in, kind = core_rknd)
      rtpthlp_in_r4        = real(rtpthlp_in, kind = core_rknd)
      sclrm_r4             = real(sclrm, kind = core_rknd)
      sclrp2_r4            = real(sclrp2, kind = core_rknd)
      sclrprtp_r4          = real(sclrprtp, kind = core_rknd)
      sclrpthlp_r4         = real(sclrpthlp, kind = core_rknd)
      wpsclrp_r4           = real(wpsclrp, kind = core_rknd)
      edsclr_in_r4         = real(edsclr_in, kind = core_rknd)
      rcm_out_r4           = real(rcm_out, kind = core_rknd)
      wprcp_out_r4         = real(wprcp_out, kind = core_rknd)
      cloud_frac_out_r4    = real(cloud_frac_out, kind = core_rknd)
      khzm_out_r4          = real(khzm_out, kind = core_rknd)
      khzt_out_r4          = real(khzt_out, kind = core_rknd)

end subroutine clubb_r82r4_core


end module clubb_r82r4_mod
