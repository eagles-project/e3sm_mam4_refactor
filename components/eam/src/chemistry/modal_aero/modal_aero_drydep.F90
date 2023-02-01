!===============================================================================
module modal_aero_drydep

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use constituents,   only: pcnst, cnst_name
  use modal_aero_data,only: cnst_name_cw
  use ppgrid,         only: pcols, pver, pverp
  use modal_aero_data,only: ntot_amode
  use aerodep_flx,    only: aerodep_flx_prescribed
  use physconst,      only: gravit, rair, rhoh2o
  use camsrfexch,     only: cam_in_t, cam_out_t
  use physics_types,  only: physics_state, physics_ptend, physics_ptend_init
  use physics_buffer, only: physics_buffer_desc
  use physics_buffer, only: pbuf_get_field, pbuf_get_index, pbuf_set_field

  use cam_history,    only: outfld

  implicit none
  private

  public :: aero_model_drydep     ! aerosol dry deposition and sediment

contains
  
  !=============================================================================
  subroutine aero_model_drydep  ( state, pbuf, obklen, ustar, cam_in, dt, cam_out, ptend )

    use aero_model,        only: drydep_lq, dgnumwet_idx, nmodes, wetdens_ap_idx
    use drydep_mod,        only: calcram
    use modal_aero_data,   only: qqcw_get_field
    use modal_aero_data,   only: cnst_name_cw
    use modal_aero_data,   only: alnsg_amode
    use modal_aero_data,   only: sigmag_amode
    use modal_aero_data,   only: nspec_amode
    use modal_aero_data,   only: numptr_amode
    use modal_aero_data,   only: numptrcw_amode
    use modal_aero_data,   only: lmassptr_amode
    use modal_aero_data,   only: lmassptrcw_amode
    use modal_aero_deposition, only: set_srf_drydep

  ! args 
    type(physics_state),target,intent(in) :: state     ! Physics state variables
    real(r8),               intent(in)    :: obklen(:) ! Obukhov length [m]
    real(r8),               intent(in)    :: ustar(:)  ! sfc friction velocity [m/s]
    type(cam_in_t), target, intent(in)    :: cam_in    ! import state
    real(r8),               intent(in)    :: dt        ! time step [s]
    type(cam_out_t),        intent(inout) :: cam_out   ! export state
    type(physics_ptend),    intent(out)   :: ptend     ! indivdual parameterization tendencies
    type(physics_buffer_desc),    pointer :: pbuf(:)

  ! local vars
    real(r8), pointer :: landfrac(:) ! land fraction
    real(r8), pointer :: icefrac(:)  ! ice fraction
    real(r8), pointer :: ocnfrac(:)  ! ocean fraction
    real(r8), pointer :: fvin(:)     !
    real(r8), pointer :: ram1in(:)   ! for dry dep velocities from land model for progseasalts

    real(r8), pointer :: tair(:,:)   ! air temperture [k]
    real(r8), pointer :: pmid(:,:)   ! air pressure at layer midpoint [Pa]
    real(r8), pointer :: pint(:,:)   ! air pressure at layer interface [Pa]
    real(r8), pointer :: pdel(:,:)   ! layer thickness [Pa]

    real(r8) :: fv(pcols)            ! for dry dep velocities, from land modified over ocean & ice
    real(r8) :: ram1(pcols)          ! for dry dep velocities, from land modified over ocean & ice

    integer :: lchnk   ! chunk identifier
    integer :: ncol    ! number of active atmospheric columns
    integer :: lspec   ! index for aerosol number / chem-mass / water-mass
    integer :: imode   ! aerosol mode index
    integer :: icnst   ! tracer index
    integer :: imnt    ! moment of the aerosol size distribution. 0 = number; 3 = volume
    integer :: jvlc    ! index for last dimension of vlc_xxx arrays

    real(r8) :: rho(pcols,pver)      ! air density in kg/m3
    real(r8) :: sflx(pcols)          ! deposition flux
    real(r8) :: dqdt_tmp(pcols,pver) ! temporary array to hold tendency for 1 species

    real(r8) :: rad_drop(pcols,pver)
    real(r8) :: dens_drop(pcols,pver)
    real(r8) :: sg_drop(pcols,pver)
    real(r8) :: rad_aer(pcols,pver)
    real(r8) :: dens_aer(pcols,pver)
    real(r8) :: sg_aer(pcols,pver)

    real(r8) :: vlc_dry(pcols,pver,4)     ! dep velocity
    real(r8) :: vlc_grv(pcols,pver,4)     ! dep velocity
    real(r8)::  vlc_trb(pcols,4)          ! dep velocity

    real(r8) :: aerdepdryis(pcols,pcnst)  ! aerosol dry deposition (interstitial)
    real(r8) :: aerdepdrycw(pcols,pcnst)  ! aerosol dry deposition (cloud water)

    real(r8), pointer :: qq(:,:)            ! mixing ratio of a single tracer [kg/kg] or [1/kg]
    real(r8), pointer :: dgncur_awet(:,:,:)
    real(r8), pointer :: wetdens(:,:,:)

    !---------------------------------------------------------------------------
    ! Retrieve input variables; initialize output (i.e., ptend).
    !---------------------------------------------------------------------------
    landfrac => cam_in%landfrac(:)
    icefrac  => cam_in%icefrac(:)
    ocnfrac  => cam_in%ocnfrac(:)
    fvin     => cam_in%fv(:)
    ram1in   => cam_in%ram1(:)

    lchnk = state%lchnk
    ncol  = state%ncol

    tair => state%t(:,:)
    pmid => state%pmid(:,:)
    pint => state%pint(:,:)
    pdel => state%pdel(:,:)

    rho(:ncol,:)=  pmid(:ncol,:)/(rair*tair(:ncol,:))

    call pbuf_get_field(pbuf, dgnumwet_idx,   dgncur_awet, start=(/1,1,1/), kount=(/pcols,pver,nmodes/) ) 
    call pbuf_get_field(pbuf, wetdens_ap_idx, wetdens,     start=(/1,1,1/), kount=(/pcols,pver,nmodes/) ) 

    call physics_ptend_init(ptend, state%psetcols, 'aero_model_drydep_ma', lq=drydep_lq)

    !---------------------------------------------------------------------------
    ! For turbulent dry deposition: calculate ram and fv over ocean and sea ice; 
    ! copy values over land
    !---------------------------------------------------------------------------
    call calcram( ncol,         &! in: state%ncol
                  landfrac,     &! in: cam_in%landfrac
                  icefrac,      &! in: cam_in%icefrac
                  ocnfrac,      &! in: cam_in%ocnfrac
                  obklen,       &! in: calculated in tphysac
                  ustar,        &! in: calculated in tphysac
                  tair(:,pver), &! in. note: bottom level only
                  pmid(:,pver), &! in. note: bottom level only
                  pdel(:,pver), &! in. note: bottom level only
                  ram1in,       &! in: cam_in%ram1
                  fvin,         &! in: cam_in%fv
                  ram1,         &! out: aerodynamical resistance (s/m)
                  fv            &! out: friction velocity
                  )

    call outfld( 'airFV', fv(:), pcols, lchnk )
    call outfld( 'RAM1', ram1(:), pcols, lchnk )
 
    !======================
    ! cloud-borne aerosols
    !---------------------------------------------------------------------------------------
    ! Calculate gravitational settling and dry deposition velocities for cloud droplets 
    ! (and hence the cloud-borne aerosols therein).
    ! There is one set of velocities for number mixing ratios of all aerosol modes
    ! and one set of velocities for all mass mixing ratios of all modes.
    !---------------------------------------------------------------------------------------
    ! *** mean drop radius should eventually be computed from ndrop and qcldwtr

    rad_drop(:,:) = 5.0e-6_r8
    dens_drop(:,:) = rhoh2o
    sg_drop(:,:) = 1.46_r8

    jvlc = 3 ; imnt = 0  ! cloud-borne aerosol number
    call modal_aero_depvel_part( ncol, lchnk, tair, pmid, ram1, fv,  &! in
                                 rad_drop, dens_drop, sg_drop, imnt, &! in
                                 vlc_dry(:,:,jvlc),                  &! out
                                 vlc_trb(:,  jvlc),                  &! out
                                 vlc_grv(:,:,jvlc)                   )! out

    jvlc = 4 ; imnt = 3  ! cloud-borne aerosol volume/mass
    call modal_aero_depvel_part( ncol, lchnk, tair, pmid, ram1, fv,  &! in
                                 rad_drop, dens_drop, sg_drop, imnt, &! in
                                 vlc_dry(:,:,jvlc),                  &! out
                                 vlc_trb(:,  jvlc),                  &! out
                                 vlc_grv(:,:,jvlc)                   )! out

    !----------------------------------------------------------------------------------
    ! Loop over all modes and all aerosol tracers (number + mass species).
    ! Calculate the drydep-induced tendencies, then update the mixing ratios.
    !----------------------------------------------------------------------------------
    do imode = 1, ntot_amode         ! main loop over aerosol modes
    do lspec = 0, nspec_amode(imode) ! loop over number + constituents

       if (lspec == 0) then   ! number
           icnst = numptrcw_amode(imode) ; jvlc = 3
       else ! aerosol mass
           icnst = lmassptrcw_amode(lspec,imode) ; jvlc = 4
       endif

       qq => qqcw_get_field(pbuf,icnst,lchnk)
       call sedimentation_solver_for_1_tracer( ncol, dt, vlc_dry(:,:,jvlc), qq,    &! in
                                               gravit, rho, tair, pint, pmid, pdel,&! in
                                               dqdt_tmp, sflx                      )! out

       aerdepdrycw(:ncol,icnst) = sflx(:ncol)

       ! Update mixing ratios here. Recall that mixing ratios of cloud-borne aerosols
       ! are stored in pbuf, not as part of the state variable

       qq(1:ncol,:) = qq(1:ncol,:) + dqdt_tmp(1:ncol,:) * dt

       ! Get and save diagnostics

       call drydep_diags_for_1_tracer( lchnk, ncol, trim(cnst_name_cw(icnst)), &! in
                                       vlc_dry(:,:,jvlc), vlc_trb(:,jvlc),     &! in
                                       vlc_grv(:,:,jvlc), sflx                 )! in

    enddo ! loop over number + constituents
    enddo ! imode = 1, ntot_amode

    !====================
    ! interstial aerosol
    !====================
    do imode = 1, ntot_amode   ! main loop over aerosol modes

       !-----------------------------------------------------------------
       ! Calculate gravitational settling and dry deposition velocities for 
       ! interstitial aerosol particles in a single lognormal mode. Note:
       !  One set of velocities for number mixing ratio of the mode;
       !  One set of velocities for all mass mixing ratios of the mode.
       !-----------------------------------------------------------------
       ! rad_aer = volume mean wet radius (m)
       ! dgncur_awet = geometric mean wet diameter for number distribution (m)

       rad_aer(1:ncol,:) = 0.5_r8*dgncur_awet(1:ncol,:,imode)   &
                           *exp(1.5_r8*(alnsg_amode(imode)**2))

       ! dens_aer(1:ncol,:) = wet density (kg/m3)
       dens_aer(1:ncol,:) = wetdens(1:ncol,:,imode)
       sg_aer(1:ncol,:) = sigmag_amode(imode)

       jvlc = 1  ; imnt = 0  ! interstitial aerosol number
       call modal_aero_depvel_part( ncol, lchnk, tair, pmid, ram1, fv, &! in
                                    rad_aer, dens_aer, sg_aer, imnt,   &! in
                                    vlc_dry(:,:,jvlc),                 &! out
                                    vlc_trb(:,  jvlc),                 &! out
                                    vlc_grv(:,:,jvlc)                  )! out

       jvlc = 2  ; imnt = 3  ! interstitial aerosol volume/mass
       call modal_aero_depvel_part( ncol, lchnk, tair, pmid, ram1, fv, &! in
                                    rad_aer, dens_aer, sg_aer, imnt,   &! in
                                    vlc_dry(:,:,jvlc),                 &! out
                                    vlc_trb(:,  jvlc),                 &! out
                                    vlc_grv(:,:,jvlc)                  )! out

       !-----------------------------------------------------------
       ! Loop over number + mass species of the mode. 
       ! Calculate drydep-induced tendencies; save to ptend.
       !-----------------------------------------------------------
       do lspec = 0, nspec_amode(imode)

          if (lspec == 0) then   ! number
             icnst = numptr_amode(imode) ; jvlc = 1
          else ! aerosol mass
             icnst = lmassptr_amode(lspec,imode) ; jvlc = 2
          endif

          qq => state%q(:,:,icnst)
          call sedimentation_solver_for_1_tracer( ncol, dt, vlc_dry(:,:,jvlc), qq,    &! in
                                                  gravit, rho, tair, pint, pmid, pdel,&! in
                                                  dqdt_tmp, sflx                      )! out

          aerdepdryis(:ncol,icnst) = sflx(:ncol)
          ptend%lq(icnst) = .TRUE.
          ptend%q(:ncol,:,icnst) = dqdt_tmp(:ncol,:)

          ! Get and save diagnostics

          call drydep_diags_for_1_tracer( lchnk, ncol, trim(cnst_name(icnst)), &! in
                                          vlc_dry(:,:,jvlc), vlc_trb(:,jvlc),  &! in
                                          vlc_grv(:,:,jvlc), sflx,             &! in
                                          ptend%q(:,:,icnst)                   )! in

          call outfld( trim(cnst_name(icnst))//'DDV', vlc_dry(:ncol,:,jvlc), pcols, lchnk )

      enddo ! lspec = 1, nspec_amode(m)
    enddo   ! imode = 1, ntot_amode

    !=====================================================================
    ! if the user has specified prescribed aerosol dep fluxes then 
    ! do not set cam_out dep fluxes according to the prognostic aerosols
    !=====================================================================
    if (.not.aerodep_flx_prescribed()) then
       call set_srf_drydep(aerdepdryis, aerdepdrycw, cam_out)
    endif

  end subroutine aero_model_drydep

  !-----------------------------------------------------------------------
  ! Numerically solve the sedimentation equation for 1 tracer
  !-----------------------------------------------------------------------
  subroutine sedimentation_solver_for_1_tracer( ncol, dt, sed_vel, qq_in,            &! in
                                                gravit, rho, tair, pint, pmid, pdel, &! in
                                                dqdt_sed, sflx                       )! out

    use shr_kind_mod,      only: r8 => shr_kind_r8
    use ppgrid,            only: pcols, pver, pverp
   !use dust_sediment_mod, only: dust_sediment_tend
    use dust_sediment_mod, only: getflx 

    integer , intent(in) :: ncol
    real(r8), intent(in) :: dt
    real(r8), intent(in) :: gravit
    real(r8), intent(in) :: rho(pcols,pver)           ! air density in kg/m3
    real(r8), intent(in) :: tair(pcols,pver)          ! air temperature
    real(r8), intent(in) :: pint(pcols,pverp)         ! air pressure at layer interfaces 
    real(r8), intent(in) :: pmid(pcols,pver)          ! air pressure at layer midpoints
    real(r8), intent(in) :: pdel(pcols,pver)          ! pressure layer thickness
    real(r8), intent(in) :: sed_vel(pcols,pver)       ! dep velocity
    real(r8), intent(in) :: qq_in(pcols,pver)

    real(r8), intent(out) :: dqdt_sed(pcols,pver) ! tendency
    real(r8), intent(out) :: sflx(pcols)          ! deposition flux at the Earth's surface

    real (r8), parameter :: mxsedfac = 0.99_r8    ! maximum sedimentation flux factor

    real(r8) :: pvmzaer(pcols,pverp)     ! sedimentation velocity in Pa (positive = down)
    real(r8) :: dtmassflux(pcols,pverp)  ! dt * mass fluxes at layer interfaces (positive = down)

    integer  :: ii,kk

    !---------------------------------------------------------------------------------------
    ! Set sedimentation velocity to zero at the top interface of the model domain.
    ! (This was referred to as the "pvprogseasalts method" in the code before refactoring.)

    pvmzaer(:ncol,1)=0._r8

    ! Assume the sedimentation velocities passed in are velocities
    ! at the bottom interface of each model layer.

    pvmzaer(:ncol,2:pverp) = sed_vel(:ncol,:)

    ! Convert from velocity to (gravitiy * mass fluxes of the air);
    ! units: convert from meters/sec to pascals/sec.
    ! (This was referred to as "Phil's method" in the code before refactoring.)
 
    pvmzaer(:ncol,2:pverp) = pvmzaer(:ncol,2:pverp) * rho(:ncol,:)*gravit

    !------------------------------------------------------
    ! Calculate mass flux * dt at each layer interface
    !------------------------------------------------------
    call getflx(ncol, pint, qq_in, pvmzaer, dt, dtmassflux)

    ! Filter out any negative fluxes from the getflx routine

    do kk = 2,pver
       dtmassflux(:ncol,kk) = max(0._r8, dtmassflux(:ncol,kk))
    end do

    ! Upper and lower boundaries

    do ii = 1,ncol
       dtmassflux(ii,1)     = 0                                         ! no flux at model top 
       dtmassflux(ii,pverp) = qq_in(ii,pver) * pvmzaer(ii,pverp) * dt   ! surface flux by upstream scheme
    end do

    ! Limit the flux out of the bottom of each column:
    ! apply mxsedfac to prevent generating very small negative mixing ratio.
    ! *** Should we include the flux through the top interface, to accommodate thin surface layers?

    do kk = 1,pver
       do ii = 1,ncol
          dtmassflux(ii,kk+1) = min( dtmassflux(ii,kk+1), mxsedfac * qq_in(ii,kk) * pdel(ii,kk))
       end do
    end do

    !-----------------------------------------------------------------------
    ! Calculate the mixing ratio tendencies resulting from flux divergence
    !-----------------------------------------------------------------------
    do kk = 1,pver
       do ii = 1,ncol
          dqdt_sed(ii,kk)  = (dtmassflux(ii,kk) - dtmassflux(ii,kk+1)) / (dt * pdel(ii,kk))
       end do
    end do

    !-----------------------------------------------------------------------
    ! Convert flux out the bottom to mass units Pa -> kg/m2/s
    !-----------------------------------------------------------------------
    sflx(:ncol) = dtmassflux(:ncol,pverp) / (dt*gravit)

  end subroutine sedimentation_solver_for_1_tracer


  subroutine drydep_diags_for_1_tracer( lchnk, ncol, cnst_name_in, vlc_dry, vlc_trb, vlc_grv, sflx, dqdt_sed )

    integer, intent(in) :: lchnk  ! chunk index
    integer, intent(in) :: ncol   ! # of active columns 

    character(len=*), intent(in) :: cnst_name_in  ! tracer name

    real(r8),intent(in) :: vlc_dry(pcols,pver)
    real(r8),intent(in) :: vlc_trb(pcols)
    real(r8),intent(in) :: vlc_grv(pcols,pver)
    real(r8),intent(in) ::    sflx(pcols)

    real(r8),intent(in),optional :: dqdt_sed(pcols,pver)

    real(r8) :: dep_trb(pcols)       !kg/m2/s
    real(r8) :: dep_grv(pcols)       !kg/m2/s (total of grav and trb)

    integer :: ii

    ! apportion dry deposition into turb and gravitational settling for tapes

    do ii=1,ncol
       if (vlc_dry(ii,pver) .ne. 0._r8) then
          dep_trb(ii)=sflx(ii)*vlc_trb(ii)/vlc_dry(ii,pver)
          dep_grv(ii)=sflx(ii)*vlc_grv(ii,pver)/vlc_dry(ii,pver)
       endif
    enddo

    ! send diagnostics to output

    call outfld( cnst_name_in//'DDF', sflx,     pcols, lchnk)
    call outfld( cnst_name_in//'TBF', dep_trb,  pcols, lchnk)
    call outfld( cnst_name_in//'GVF', dep_grv,  pcols, lchnk)

    if (present(dqdt_sed)) &
    call outfld( cnst_name_in//'DTQ', dqdt_sed, pcols, lchnk)

  end subroutine drydep_diags_for_1_tracer
  !=============================================================================

  !==========================================================================================
  ! Calculate deposition velocities caused by turbulent dry deposition and
  ! gravitational settling of aerosol particles

  ! Reference: 
  !  L. Zhang, S. Gong, J. Padro, and L. Barrie:
  !  A size-seggregated particle dry deposition scheme for an atmospheric aerosol module
  !  Atmospheric Environment, 35, 549-560, 2001.
  !
  ! Authors: X. Liu
  !==========================================================================================
  subroutine modal_aero_depvel_part( ncol, lchnk, tair, pmid, ram1, fv,            &! in
                                     radius_part, density_part, sig_part, moment,  &! in
                                     vlc_dry, vlc_trb, vlc_grv                     )! out

    use physconst,     only: pi,boltz, gravit, rair
    use mo_drydep,     only: n_land_type, fraction_landuse
    use ieee_arithmetic, only: ieee_is_nan
    use phys_control,    only: phys_getopts

    implicit none

    integer,  intent(in) :: ncol
    integer,  intent(in) :: lchnk

    real(r8), intent(in) :: tair(pcols,pver)    ! air temperature [K]
    real(r8), intent(in) :: pmid(pcols,pver)    ! air pressure [Pa]
    real(r8), intent(in) :: fv(pcols)           ! friction velocity [m/s]
    real(r8), intent(in) :: ram1(pcols)         ! aerodynamical resistance [s/m]

    real(r8), intent(in) :: radius_part(pcols,pver)    ! mean (volume or number) particle radius [m]
    real(r8), intent(in) :: density_part(pcols,pver)   ! density of particle material [kg/m3]
    real(r8), intent(in) :: sig_part(pcols,pver)       ! geometric standard deviation of particle size distribution
    integer,  intent(in) :: moment                     ! moment of size distribution (0 for number, 2 for surface area, 3 for volume)

    real(r8), intent(out) :: vlc_trb(pcols)         ! turbulent dry deposition velocity [m/s]
    real(r8), intent(out) :: vlc_grv(pcols,pver)    ! gravitational deposition velocity [m/s]
    real(r8), intent(out) :: vlc_dry(pcols,pver)    ! total     dry deposition velocity [m/s]
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------
    ! Local Variables

    integer  :: ii, kk                    ! grid column and layer indices
    real(r8) :: rho                       ! air density [kg/m**3]
    real(r8) :: vsc_dyn_atm(pcols,pver)   ! [kg m-1 s-1] Dynamic viscosity of air
    real(r8) :: vsc_knm_atm(pcols,pver)   ! [m2 s-1] Kinematic viscosity of atmosphere
    real(r8) :: mfp_atm(pcols,pver)       ! [m] Mean free path of air
    real(r8) :: slp_crc(pcols,pver)       ! [frc] Slip correction factor
    real(r8) :: radius_moment(pcols,pver) ! median radius [m] for moment

    real(r8) :: shm_nbr       ! [frc] Schmidt number
    real(r8) :: stk_nbr       ! [frc] Stokes number
    real(r8) :: dff_aer       ! [m2 s-1] Brownian diffusivity of particle
    real(r8) :: rss_trb       ! [s m-1] Resistance to turbulent deposition
    real(r8) :: rss_lmn       ! [s m-1] Quasi-laminar layer resistance
    real(r8) :: brownian      ! collection efficiency for Browning diffusion
    real(r8) :: impaction     ! collection efficiency for impaction
    real(r8) :: interception  ! collection efficiency for interception
    real(r8) :: stickfrac     ! fraction of particles sticking to surface
    real(r8) :: lnsig         ! ln(sig_part)
    real(r8) :: dispersion    ! accounts for influence of size dist dispersion on bulk settling velocity
                              ! assuming radius_part is number mode radius * exp(1.5 ln(sigma))

    integer  :: lt              ! land type index
    real(r8) :: lnd_frc         ! land type fraction [unitless]
    real(r8) :: vlc_trb_ontype  ! turbulent dry dep. velocity on single land type [m/s]
    real(r8) :: vlc_trb_wgtsum  ! turbulent dry dep. velocity averaged over land types [m/s]
    real(r8) :: vlc_dry_wgtsum  ! total     dry dep. velocity averaged over land types [m/s]

    logical  :: use_MMF

    ! constants

    real(r8) gamma(11)      ! exponent of schmidt number
    data gamma/0.56e+00_r8,  0.54e+00_r8,  0.54e+00_r8,  0.56e+00_r8,  0.56e+00_r8, &        
               0.56e+00_r8,  0.50e+00_r8,  0.54e+00_r8,  0.54e+00_r8,  0.54e+00_r8, &
               0.54e+00_r8/
    save gamma

    real(r8) alpha(11)      ! parameter for impaction
    data alpha/1.50e+00_r8,   1.20e+00_r8,  1.20e+00_r8,  0.80e+00_r8,  1.00e+00_r8, &
               0.80e+00_r8, 100.00e+00_r8, 50.00e+00_r8,  2.00e+00_r8,  1.20e+00_r8, &
              50.00e+00_r8/
    save alpha

    real(r8) radius_collector(11) ! radius (m) of surface collectors
    data radius_collector/10.00e-03_r8,  3.50e-03_r8,  3.50e-03_r8,  5.10e-03_r8,  2.00e-03_r8, &
                           5.00e-03_r8, -1.00e+00_r8, -1.00e+00_r8, 10.00e-03_r8,  3.50e-03_r8, &
                          -1.00e+00_r8/
    save radius_collector

    integer  :: iwet(11) ! flag for wet surface = 1, otherwise = -1
    data iwet/-1,  -1,   -1,   -1,   -1,  &
              -1,   1,   -1,    1,   -1,  &
              -1/
    save iwet


    ! use a maximum radius of 50 microns when calculating deposition velocity

    real(r8),parameter :: radiaus_max = 50.0e-6_r8

    call phys_getopts( use_MMF_out = use_MMF )

    !------------------------------------------------------------------------
    ! Calculate particle velocity of gravitational settling
    !------------------------------------------------------------------------
    do kk=1,pver
       do ii=1,ncol

          lnsig = log(sig_part(ii,kk))
          radius_moment(ii,kk) = min(radiaus_max,radius_part(ii,kk))*exp((float(moment)-1.5_r8)*lnsig*lnsig)
          dispersion = exp(2._r8*lnsig*lnsig)

          rho=pmid(ii,kk)/rair/tair(ii,kk)

          ! Quasi-laminar layer resistance: call rss_lmn_get
          ! Size-independent thermokinetic properties

          vsc_dyn_atm(ii,kk) = 1.72e-5_r8 * ( (tair(ii,kk)/273.0_r8)**1.5_r8) * 393.0_r8 &
                                             /(tair(ii,kk)+120.0_r8)  ![kg m-1 s-1] RoY94 p. 102
          mfp_atm(ii,kk) = 2.0_r8 * vsc_dyn_atm(ii,kk) / &   ![m] SeP97 p. 455
                          ( pmid(ii,kk)*sqrt(8.0_r8/(pi*rair*tair(ii,kk))) )
          vsc_knm_atm(ii,kk) = vsc_dyn_atm(ii,kk) / rho ![m2 s-1] Kinematic viscosity of air

          slp_crc(ii,kk) = 1.0_r8 + mfp_atm(ii,kk) * &
                           (1.257_r8+0.4_r8*exp(-1.1_r8*radius_moment(ii,kk)/(mfp_atm(ii,kk)))) / &
                           radius_moment(ii,kk)   ![frc] Slip correction factor SeP97 p. 464

          vlc_grv(ii,kk) = (4.0_r8/18.0_r8) * radius_moment(ii,kk)*radius_moment(ii,kk)*density_part(ii,kk)* &
                           gravit*slp_crc(ii,kk) / vsc_dyn_atm(ii,kk) ![m s-1] Stokes' settling velocity SeP97 p. 466
          vlc_grv(ii,kk) = vlc_grv(ii,kk) * dispersion

          ! in the MMF NaN's were occurring here but the root cause was not
          ! identified, so this check was added to work around the issue
          if (use_MMF) then
             if ( ieee_is_nan(vlc_grv(ii,kk)) ) vlc_grv(ii,kk) = 0.0_r8 
          end if

       enddo
    enddo
    vlc_dry(:ncol,:)=vlc_grv(:ncol,:)

    !------------------------------------------------------------------------------------
    ! Calculate particle velocity of turbulent dry deposition.
    ! This process is assumed to only occur in the lowest model layer.
    !------------------------------------------------------------------------------------
    kk = pver
    do ii=1,ncol

       dff_aer = boltz * tair(ii,kk) * slp_crc(ii,kk) / &    ![m2 s-1]
                 (6.0_r8*pi*vsc_dyn_atm(ii,kk)*radius_moment(ii,kk)) !SeP97 p.474

       shm_nbr = vsc_knm_atm(ii,kk) / dff_aer                        ![frc] SeP97 p.972

       vlc_trb_wgtsum = 0._r8
       vlc_dry_wgtsum = 0._r8

       do lt = 1,n_land_type

          lnd_frc = fraction_landuse(ii,lt,lchnk)

          if ( lnd_frc /= 0._r8 ) then
             brownian = shm_nbr**(-gamma(lt))

             if (radius_collector(lt) > 0.0_r8) then ! vegetated surface
                stk_nbr = vlc_grv(ii,kk) * fv(ii) / (gravit*radius_collector(lt))
                interception = 2.0_r8*(radius_moment(ii,kk)/radius_collector(lt))**2.0_r8

             else ! non-vegetated surface
                stk_nbr = vlc_grv(ii,kk) * fv(ii) * fv(ii) / (gravit*vsc_knm_atm(ii,kk))  ![frc] SeP97 p.965
                interception = 0.0_r8
             endif
             impaction = (stk_nbr/(alpha(lt)+stk_nbr))**2.0_r8   

             if (iwet(lt) > 0) then
                stickfrac = 1.0_r8
             else
                stickfrac = exp(-sqrt(stk_nbr))
                if (stickfrac < 1.0e-10_r8) stickfrac = 1.0e-10_r8
             endif
             rss_lmn = 1.0_r8 / (3.0_r8 * fv(ii) * stickfrac * (brownian+interception+impaction))
             rss_trb = ram1(ii) + rss_lmn + ram1(ii)*rss_lmn*vlc_grv(ii,kk)

             vlc_trb_ontype = 1.0_r8 / rss_trb
             vlc_trb_wgtsum = vlc_trb_wgtsum + lnd_frc*( vlc_trb_ontype )
             vlc_dry_wgtsum = vlc_dry_wgtsum + lnd_frc*( vlc_trb_ontype + vlc_grv(ii,kk) )
          endif
       enddo  ! n_land_type

       vlc_trb(ii)    = vlc_trb_wgtsum
       vlc_dry(ii,kk) = vlc_dry_wgtsum

    enddo !ncol

  end subroutine modal_aero_depvel_part

end module modal_aero_drydep
