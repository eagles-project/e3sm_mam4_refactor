!===============================================================================
! Calculations of 
!  - changes in interstitial aerosol mixing ratios caused by
!    gravitational settling and turbulent dry deposition of aerosol particles,
!  - changes in cloud-borne aerosol mixing ratios caused by 
!    gravitational settling and turbulent dry deposition of cloud droplets.
!===============================================================================
module modal_aero_drydep

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use constituents,   only: pcnst, cnst_name
  use modal_aero_data,only: cnst_name_cw
  use ppgrid,         only: pcols, pver, pverp
  use modal_aero_data,only: ntot_amode
  use aerodep_flx,    only: aerodep_flx_prescribed
  use physconst,      only: gravit, rair, rhoh2o, pi, boltz
  use camsrfexch,     only: cam_in_t, cam_out_t
  use physics_types,  only: physics_state, physics_ptend, physics_ptend_init
  use physics_buffer, only: physics_buffer_desc
  use physics_buffer, only: pbuf_get_field, pbuf_get_index, pbuf_set_field

  use cam_history,    only: outfld

  implicit none
  private

  public :: aero_model_drydep

contains
  
  !=============================================================================
  ! Main subroutine of aerosol dry deposition parameterization.
  ! Also serves as the interface routine called by EAM's physics driver.
  !=============================================================================
  subroutine aero_model_drydep  ( state, pbuf, obklen, ustar, cam_in, dt, cam_out, ptend )

    use aero_model,        only: drydep_lq, dgnumwet_idx, nmodes, wetdens_ap_idx
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

    ! Arguments

    type(physics_state),target,intent(in) :: state     ! Physics state variables
    real(r8),               intent(in)    :: obklen(:) ! Obukhov length [m]
    real(r8),               intent(in)    :: ustar(:)  ! sfc friction velocity [m/s]
    type(cam_in_t), target, intent(in)    :: cam_in    ! import state
    real(r8),               intent(in)    :: dt        ! time step [s]
    type(cam_out_t),        intent(inout) :: cam_out   ! export state
    type(physics_ptend),    intent(out)   :: ptend     ! indivdual parameterization tendencies
    type(physics_buffer_desc), pointer    :: pbuf(:)

    ! Local variables

    real(r8), pointer :: landfrac(:) ! land fraction [unitless]
    real(r8), pointer :: icefrac(:)  ! ice fraction [unitless]
    real(r8), pointer :: ocnfrac(:)  ! ocean fraction [unitless]
    real(r8), pointer :: fricvelin(:)! friction velocity from land model [m/s]
    real(r8), pointer :: ram1in(:)   ! aerodynamical resistance from land model [s/m]

    real(r8), pointer :: tair(:,:)   ! air temperture [k]
    real(r8), pointer :: pmid(:,:)   ! air pressure at layer midpoint [Pa]
    real(r8), pointer :: pint(:,:)   ! air pressure at layer interface [Pa]
    real(r8), pointer :: pdel(:,:)   ! layer thickness [Pa]

    real(r8) :: fricvel(pcols)     ! friction velocity used in the calculaiton of turbulent dry deposition velocity [m/s]
    real(r8) ::    ram1(pcols)     ! aerodynamical resistance used in the calculaiton of turbulent dry deposition velocity [s/m]

    integer :: lchnk   ! chunk identifier
    integer :: ncol    ! number of active atmospheric columns
    integer :: lspec   ! index for aerosol number / chem-mass / water-mass
    integer :: imode   ! aerosol mode index
    integer :: icnst   ! tracer index
    integer :: imnt    ! moment of the aerosol size distribution. 0 = number; 3 = volume
    integer :: jvlc    ! index for last dimension of vlc_xxx arrays

    real(r8) :: rho(pcols,pver)      ! air density [kg/m3]
    real(r8) :: sflx(pcols)          ! surface deposition flux of a single species [kg/m2/s] or [1/m2/s]
    real(r8) :: dqdt_tmp(pcols,pver) ! temporary array to hold tendency for 1 species, [kg/kg/s] or [1/kg/s]

    real(r8) ::  rad_drop(pcols,pver)  ! cloud droplet radius [m]
    real(r8) :: dens_drop(pcols,pver)  ! cloud droplet density [kg/m3]
    real(r8) ::   sg_drop(pcols,pver)  ! assumed geometric standard deviation of droplet size distribution

    real(r8) ::   rad_aer(pcols,pver)  ! volume mean wet radius of interstitial aerosols [m]
    real(r8) ::  dens_aer(pcols,pver)  ! wet density of interstitial aerosols [kg/m3]
    real(r8) ::    sg_aer(pcols,pver)  ! assumed geometric standard deviation of particle size distribution

    real(r8) :: vlc_dry(pcols,pver,4)     ! dep velocity [m/s]
    real(r8) :: vlc_grv(pcols,pver,4)     ! dep velocity [m/s]
    real(r8)::  vlc_trb(pcols,4)          ! dep velocity [m/s]

    real(r8) :: aerdepdryis(pcols,pcnst)  ! surface deposition flux of interstitial aerosols, [kg/m2/s] or [1/m2/s]
    real(r8) :: aerdepdrycw(pcols,pcnst)  ! surface deposition flux of cloud-borne  aerosols, [kg/m2/s] or [1/m2/s]

    ! The pointers below are used for retrieving information from pbuf

    real(r8), pointer :: qq(:,:)            ! mixing ratio of a single tracer [kg/kg] or [1/kg]
    real(r8), pointer :: dgncur_awet(:,:,:) ! geometric mean wet diameter for number distribution [m]
    real(r8), pointer :: wetdens(:,:,:)     ! wet density of interstitial aerosol [kg/m3]

    !---------------------------------------------------------------------------
    ! Retrieve input variables; initialize output (i.e., ptend).
    !---------------------------------------------------------------------------
    landfrac  => cam_in%landfrac(:)
    icefrac   => cam_in%icefrac(:)
    ocnfrac   => cam_in%ocnfrac(:)
    fricvelin => cam_in%fv(:)
    ram1in    => cam_in%ram1(:)

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

    !--------------------------------------------------------------------------------
    ! For turbulent dry deposition: calculate ram and fricvel over ocean and sea ice; 
    ! copy values over land
    !--------------------------------------------------------------------------------
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
                  fricvelin,    &! in: cam_in%fv
                  ram1,         &! out: aerodynamical resistance (s/m)
                  fricvel       &! out: bulk friction velocity of a grid cell
                  )

    call outfld( 'airFV', fricvel(:), pcols, lchnk )
    call outfld( 'RAM1',     ram1(:), pcols, lchnk )
 
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
    call modal_aero_depvel_part( ncol, lchnk, tair, pmid, ram1, fricvel,  &! in
                                 rad_drop, dens_drop, sg_drop, imnt,      &! in
                                 vlc_dry(:,:,jvlc),                       &! out
                                 vlc_trb(:,  jvlc),                       &! out
                                 vlc_grv(:,:,jvlc)                        )! out

    jvlc = 4 ; imnt = 3  ! cloud-borne aerosol volume/mass
    call modal_aero_depvel_part( ncol, lchnk, tair, pmid, ram1, fricvel, &! in
                                 rad_drop, dens_drop, sg_drop, imnt,     &! in
                                 vlc_dry(:,:,jvlc),                      &! out
                                 vlc_trb(:,  jvlc),                      &! out
                                 vlc_grv(:,:,jvlc)                       )! out

    !----------------------------------------------------------------------------------
    ! Loop over all modes and all aerosol tracers (number + mass species).
    ! Calculate the drydep-induced tendencies, then update the mixing ratios.
    !----------------------------------------------------------------------------------
    do imode = 1, ntot_amode         ! loop over aerosol modes
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

    !=====================
    ! interstial aerosols
    !=====================
    do imode = 1, ntot_amode   ! loop over aerosol modes

       !-----------------------------------------------------------------
       ! Calculate gravitational settling and dry deposition velocities for 
       ! interstitial aerosol particles in a single lognormal mode. Note:
       !  One set of velocities for number mixing ratio of the mode;
       !  One set of velocities for all mass mixing ratios of the mode.
       !-----------------------------------------------------------------
        rad_aer(1:ncol,:) = 0.5_r8*dgncur_awet(1:ncol,:,imode) *exp(1.5_r8*(alnsg_amode(imode)**2))
       dens_aer(1:ncol,:) = wetdens(1:ncol,:,imode)
         sg_aer(1:ncol,:) = sigmag_amode(imode)

       jvlc = 1  ; imnt = 0  ! interstitial aerosol number
       call modal_aero_depvel_part( ncol, lchnk, tair, pmid, ram1, fricvel, &! in
                                    rad_aer, dens_aer, sg_aer, imnt,        &! in
                                    vlc_dry(:,:,jvlc),                      &! out
                                    vlc_trb(:,  jvlc),                      &! out
                                    vlc_grv(:,:,jvlc)                       )! out

       jvlc = 2  ; imnt = 3  ! interstitial aerosol volume/mass
       call modal_aero_depvel_part( ncol, lchnk, tair, pmid, ram1, fricvel, &! in
                                    rad_aer, dens_aer, sg_aer, imnt,        &! in
                                    vlc_dry(:,:,jvlc),                      &! out
                                    vlc_trb(:,  jvlc),                      &! out
                                    vlc_grv(:,:,jvlc)                       )! out

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
          ptend%lq(icnst) = .true.
          ptend%q(:ncol,:,icnst) = dqdt_tmp(:ncol,:)

          ! Get and save diagnostics

          call drydep_diags_for_1_tracer( lchnk, ncol, trim(cnst_name(icnst)), &! in
                                          vlc_dry(:,:,jvlc), vlc_trb(:,jvlc),  &! in
                                          vlc_grv(:,:,jvlc), sflx,             &! in
                                          ptend%q(:,:,icnst)                   )! in

          call outfld( trim(cnst_name(icnst))//'DDV', vlc_dry(:ncol,:,jvlc), pcols, lchnk )

       enddo ! lspec = 1, nspec_amode(m)
    enddo    ! imode = 1, ntot_amode

    !=====================================================================
    ! Unless the user has specified prescribed aerosol dep fluxes,
    ! copy the fluxes calculated here to cam_out to be passed to other 
    ! components of the Earth System Model.
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
    use mo_spitfire_transport, only: getflx 

    integer , intent(in) :: ncol
    real(r8), intent(in) :: dt
    real(r8), intent(in) :: gravit                    ! gravitational acceleration [m/s2]
    real(r8), intent(in) :: rho(pcols,pver)           ! air density [kg/m3]
    real(r8), intent(in) :: tair(pcols,pver)          ! air temperature [K]
    real(r8), intent(in) :: pint(pcols,pverp)         ! air pressure at layer interfaces [Pa]
    real(r8), intent(in) :: pmid(pcols,pver)          ! air pressure at layer midpoints  [Pa]
    real(r8), intent(in) :: pdel(pcols,pver)          ! pressure layer thickness [Pa]
    real(r8), intent(in) :: sed_vel(pcols,pver)       ! deposition velocity [m/s]
    real(r8), intent(in) :: qq_in(pcols,pver)         ! tracer mixing ratio, [kg/kg] or [1/kg]

    real(r8), intent(out) :: dqdt_sed(pcols,pver) ! tracer mixing ratio tendency [kg/kg/s] or [1/kg/s]
    real(r8), intent(out) :: sflx(pcols)          ! deposition flux at the Earth's surface [kg/m2/s] or [1/m2/s]

    real (r8), parameter :: mxsedfac = 0.99_r8    ! maximum sedimentation flux factor

    real(r8) :: pvmzaer(pcols,pverp)     ! sedimentation velocity in Pa (positive = down)
    real(r8) :: dtmassflux(pcols,pverp)  ! dt * mass fluxes at layer interfaces (positive = down)

    integer  :: ii,kk

    !---------------------------------------------------------------------------------------
    ! Set sedimentation velocity to zero at the top interface of the model domain.

    pvmzaer(:ncol,1)=0._r8

    ! Assume the sedimentation velocities passed in are velocities
    ! at the bottom interface of each model layer, like an upwind scheme.

    pvmzaer(:ncol,2:pverp) = sed_vel(:ncol,:)

    ! Convert velocity from height coordinate to pressure coordinate; 
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
    enddo

    ! Set values for the upper and lower boundaries

    do ii = 1,ncol
       dtmassflux(ii,1)     = 0                                         ! no flux at model top 
       dtmassflux(ii,pverp) = qq_in(ii,pver) * pvmzaer(ii,pverp) * dt   ! surface flux by upwind scheme
    enddo

    ! Limit the flux out of the bottom of each column:
    ! apply mxsedfac to prevent generating very small negative mixing ratio.
    ! *** Should we include the flux through the top interface, to accommodate thin surface layers?

    do kk = 1,pver
       do ii = 1,ncol
          dtmassflux(ii,kk+1) = min( dtmassflux(ii,kk+1), mxsedfac * qq_in(ii,kk) * pdel(ii,kk))
       enddo
    enddo

    !-----------------------------------------------------------------------
    ! Calculate the mixing ratio tendencies resulting from flux divergence
    !-----------------------------------------------------------------------
    do kk = 1,pver
       do ii = 1,ncol
          dqdt_sed(ii,kk)  = (dtmassflux(ii,kk) - dtmassflux(ii,kk+1)) / (dt * pdel(ii,kk))
       enddo
    enddo

    !-----------------------------------------------------------------------
    ! Convert flux out the bottom to mass units [kg/m2/s]
    !-----------------------------------------------------------------------
    sflx(:ncol) = dtmassflux(:ncol,pverp) / (dt*gravit)

  end subroutine sedimentation_solver_for_1_tracer


  !=================================================================================
  ! Calculate some diagnostics for output. This does not affect time integration.
  !=================================================================================
  subroutine drydep_diags_for_1_tracer( lchnk, ncol, cnst_name_in, vlc_dry, vlc_trb, vlc_grv, sflx, dqdt_sed )

    integer, intent(in) :: lchnk  ! chunk index
    integer, intent(in) :: ncol   ! # of active columns 

    character(len=*), intent(in) :: cnst_name_in  ! tracer name

    real(r8),intent(in) :: vlc_trb(pcols)       ! deposition velocity of turbulent dry deposition [m/s]
    real(r8),intent(in) :: vlc_grv(pcols,pver)  ! deposition velocity of gravitational settling [m/s]
    real(r8),intent(in) :: vlc_dry(pcols,pver)  ! deposition velocity of both mechanisms combined  [m/s]
    real(r8),intent(in) ::    sflx(pcols)       ! total deposition flux at the surface for one species [kg/m2/s] or [1/m2/s]

    real(r8),intent(in),optional :: dqdt_sed(pcols,pver)

    real(r8) :: dep_trb(pcols)       ! turbulent dry deposition portion of sflx [kg/m2/s] or [1/m2/s]
    real(r8) :: dep_grv(pcols)       ! gravitational settling   portion of slfx [kg/m2/s] or [1/m2/s]
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
  ! History: 
  !  - Original version by X. Liu.
  !  - Calculations for gravitational and turbulent dry deposition separated into 
  !    different subroutines by Hui Wan, 2023.
  !==========================================================================================
  subroutine modal_aero_depvel_part( ncol, lchnk, tair, pmid, ram1, fricvel,       &! in
                                     radius_part, density_part, sig_part, moment,  &! in
                                     vlc_dry, vlc_trb, vlc_grv                     )! out


    integer,  intent(in) :: moment       ! moment of size distribution (0 for number, 2 for surface area, 3 for volume)
    integer,  intent(in) :: ncol         ! # of grid columns to do calculations for
    integer,  intent(in) :: lchnk        ! chunk index 

    real(r8), intent(in) :: tair(pcols,pver)    ! air temperature [K]
    real(r8), intent(in) :: pmid(pcols,pver)    ! air pressure [Pa]

    real(r8), intent(in) :: ram1(pcols)     ! aerodynamical resistance [s/m]
    real(r8), intent(in) :: fricvel(pcols)  ! sfc friction velocity from land model [m/s]

    real(r8), intent(in) :: radius_part(pcols,pver)    ! mean (volume or number) particle radius [m]
    real(r8), intent(in) :: density_part(pcols,pver)   ! density of particle material [kg/m3]
    real(r8), intent(in) :: sig_part(pcols,pver)       ! geometric standard deviation of particle size distribution

    real(r8), intent(out) :: vlc_grv(pcols,pver)    ! gravitational deposition velocity [m/s]
    real(r8), intent(out) :: vlc_dry(pcols,pver)    ! total dry deposition velocity [m/s]
    real(r8), intent(out) :: vlc_trb(pcols)         ! turbulent dry deposition velocity [m/s]

    ! use a maximum radius of 50 microns when calculating deposition velocity

    real(r8),parameter :: radius_max = 50.0e-6_r8

    !------------------------------------------------------------------------------------
    ! Calculate deposition velocity of gravitational settling in all grid layers
    !------------------------------------------------------------------------------------
    call modal_aero_gravit_settling_velocity( moment, ncol, pcols, pver, radius_max,           &! in
                                              tair, pmid, radius_part, density_part, sig_part, &! in
                                              vlc_grv                                          )! out

    ! vlc_dry is just the gravitational settling velocity for now.
    vlc_dry(:ncol,:)=vlc_grv(:ncol,:)

    !------------------------------------------------------------------------------------
    ! For the lowest model layer:
    !  - Calculate turbulent dry deposition velocity, vlc_trb.
    !  - Add vlc_trb to vlc_grv to give the total deposition velocity, vlc_dry.
    !------------------------------------------------------------------------------------
    call modal_aero_turb_drydep_velocity( moment, ncol, pcols, lchnk, radius_max,    &! in
                                          tair(:,pver), pmid(:,pver),                &! in
                                          radius_part(:,pver), density_part(:,pver), &! in
                                          sig_part(:,pver),                          &! in
                                          fricvel(:), ram1(:), vlc_grv(:,pver),      &! in
                                          vlc_trb(:), vlc_dry(:,pver)                )! out

  end subroutine modal_aero_depvel_part
  !---------------------------------------------------------------------------------
  ! !DESCRIPTION: 
  !  
  ! Calc aerodynamic resistance over oceans and sea ice from Seinfeld and Pandis, p.963.
  !  
  ! Author: Natalie Mahowald
  ! Code refactor: Hui Wan, 2023
  !---------------------------------------------------------------------------------
  subroutine calcram(ncol,landfrac,icefrac,ocnfrac, &
                     obklen,ustar,                  &
                     tair, pmid, pdel,              &
                     ram1_in, fv_in,                &
                     ram1_out,fv_out )

     implicit none
     integer, intent(in) :: ncol

     real(r8), intent(in) :: landfrac(pcols)  ! land fraction [unitless]
     real(r8), intent(in) :: icefrac(pcols)   ! ice fraction  [unitless]
     real(r8), intent(in) :: ocnfrac(pcols)   ! ocean fraction [unitless]

     real(r8), intent(in) :: obklen(pcols)    ! Obukhov length [m]
     real(r8), intent(in) :: ustar(pcols)     ! Surface friction velocity [m/s]

     real(r8), intent(in) :: tair(pcols)      ! air temperature [K]
     real(r8), intent(in) :: pmid(pcols)      ! air pressure [Pa]
     real(r8), intent(in) :: pdel(pcols)      ! layer pressure thickness [Pa]

     real(r8),intent(in) :: ram1_in(pcols)     ! aerodynamical resistance [s/m]
     real(r8),intent(in) :: fv_in(pcols)       ! sfc friction velocity from land model [m/s]

     real(r8),intent(out) :: ram1_out(pcols)   ! aerodynamical resistance [s/m]
     real(r8),intent(out) :: fv_out(pcols)     ! sfc friction velocity of the current grid cell [m/s]

     real(r8), parameter :: zzocen = 0.0001_r8   ! Ocean aerodynamic roughness length
     real(r8), parameter :: zzsice = 0.0400_r8   ! Sea ice aerodynamic roughness length
     real(r8), parameter :: xkar   = 0.4_r8      ! Von Karman constant

     ! local variables

     real(r8) :: zz,psi,psi0,nu,nu0,temp
     integer :: ii

     real(r8), parameter :: lndfrc_threshold = 0.000000001_r8  ! fraction, unitless

     !---------------------------------------------------------------------------
     ! Friction velocity:
     !  - If the grid cell has a land fraction larger than a threshold (~zero),
     !    then use cam_in%fv.
     !  - Otherwise, use the ustar calculated in the atmosphere.
     !---------------------------------------------------------------------------
     where( landfrac(:ncol) > lndfrc_threshold )
       fv_out(:ncol) = fv_in(:ncol)
     elsewhere
       fv_out(:ncol) = ustar(:ncol)
     endwhere

     ! fvitt -- fv == 0 causes a floating point exception in
     ! dry dep of sea salts and dust

     where ( fv_out(:ncol) == 0._r8 )
        fv_out(:ncol) = 1.e-12_r8
     endwhere

     !-------------------------------------------------------------------
     ! Aerodynamic resistence
     !-------------------------------------------------------------------
     do ii=1,ncol

        ! If the grid cell has a land fraction larger than a threshold (~zero),
        ! simply use cam_in%ram1

        if (landfrac(ii) > lndfrc_threshold) then
           ram1_out(ii)=ram1_in(ii)

        else
        ! If the grid cell has a land fraction smaller than the threshold,
        ! calculate aerodynamic resistence


           ! calculate psi, psi0, temp

           zz=pdel(ii)*rair*tair(ii)/pmid(ii)/gravit/2.0_r8   !use half the layer height like Ganzefeld and Lelieveld, 1995
           if(obklen(ii).eq.0) then
              psi=0._r8
              psi0=0._r8
           else
              psi=min(max(zz/obklen(ii),-1.0_r8),1.0_r8)
              psi0=min(max(zzocen/obklen(ii),-1.0_r8),1.0_r8)
           endif

           temp=zz/zzocen

           ! special treatment for ice-dominant cells

           if(icefrac(ii) > 0.5_r8) then 
              if(obklen(ii).gt.0) then 
                 psi0=min(max(zzsice/obklen(ii),-1.0_r8),1.0_r8)
              else
                 psi0=0.0_r8
              endif
              temp=zz/zzsice
           endif

           ! calculate aerodynamic resistence

           if(psi> 0._r8) then
              ram1_out(ii) =1/xkar/ustar(ii)*(log(temp)+4.7_r8*(psi-psi0))
           else
              nu=(1.00_r8-15.000_r8*psi)**(.25_r8)
              nu0=(1.000_r8-15.000_r8*psi0)**(.25_r8)
              if(ustar(ii).ne.0._r8) then
                 ram1_out(ii) =1/xkar/ustar(ii)*(log(temp) &
                      +log(((nu0**2+1.00_r8)*(nu0+1.0_r8)**2)/((nu**2+1.0_r8)*(nu+1.00_r8)**2)) &
                      +2.0_r8*(atan(nu)-atan(nu0)))
              else
                 ram1_out(ii) =0._r8
              endif
           endif

        endif  ! if grid cell has land fraction > threshold or not

     enddo ! loop over grid columns

  end subroutine calcram

  !==========================================================================
  ! Calculate particle velocity of gravitational settling
  !==========================================================================
  subroutine modal_aero_gravit_settling_velocity( moment, ncol, pcols, nver, radius_max,           &! in
                                                  tair, pmid, radius_part, density_part, sig_part, &! in
                                                  vlc_grv                                          )! out

    ! Arguments

    integer,  intent(in) :: moment       ! moment of size distribution (0 for number, 2 for surface area, 3 for volume)
    integer,  intent(in) :: ncol         ! # of grid columns to do calculations for
    integer,  intent(in) :: pcols        ! dimension size (# of columns) 
    integer,  intent(in) :: nver         ! dimension size (# of model layers) 

    real(r8), intent(in) :: radius_max        ! upper bound of radius used for the calculation of deposition velocity [m]

    real(r8), intent(in) :: tair(pcols,nver)    ! air temperature [K]
    real(r8), intent(in) :: pmid(pcols,nver)    ! air pressure [Pa]

    real(r8), intent(in) :: radius_part(pcols,nver)    ! mean (volume or number) particle radius [m]
    real(r8), intent(in) :: density_part(pcols,nver)   ! density of particle material [kg/m3]
    real(r8), intent(in) :: sig_part(pcols,nver)       ! geometric standard deviation of particle size distribution

    real(r8), intent(out) :: vlc_grv(pcols,nver)    ! gravitational deposition velocity [m/s]

    ! Local Variables

    integer  :: ii, kk        ! grid column and layer indices
    real(r8) :: vsc_dyn_atm   ! Dynamic viscosity of air [kg m-1 s-1]
    real(r8) :: slp_crc       ! Slip correction factor [unitless]
    real(r8) :: radius_moment ! median radius for moment [m]

    !-----------
    do kk=1,nver
    do ii=1,ncol

       vsc_dyn_atm = air_dynamic_viscosity( tair(ii,kk) )

       radius_moment = radius_for_moment( moment,sig_part(ii,kk),radius_part(ii,kk),radius_max )

       slp_crc = slip_correction_factor( vsc_dyn_atm, pmid(ii,kk), tair(ii,kk), rair, pi, radius_moment ) 

       vlc_grv(ii,kk) = gravit_settling_velocity( radius_moment, density_part(ii,kk),      &
                                                  gravit, slp_crc, vsc_dyn_atm, sig_part(ii,kk) )
    enddo
    enddo
    !----

  end subroutine modal_aero_gravit_settling_velocity

  !------------------------------------------------------------------------------------
  ! Calculate particle velocity of turbulent dry deposition.
  ! This process is assumed to only occur in the lowest model layer.
  !------------------------------------------------------------------------------------
  subroutine modal_aero_turb_drydep_velocity( moment, ncol, pcols, lchnk, radius_max,  &! in
                                              tair, pmid,                              &! in
                                              radius_part, density_part, sig_part,     &! in
                                              fricvel, ram1, vlc_grv,                  &! in
                                              vlc_trb, vlc_dry                         )! out

    use mo_drydep,     only: n_land_type, fraction_landuse

    implicit none

    integer,  intent(in) :: moment           ! moment of size distribution (0 for number, 2 for surface area, 3 for volume)
    integer,  intent(in) :: ncol             ! # of grid columns to do calculations for
    integer,  intent(in) :: pcols            ! dimension size (# of grid columns)
    integer,  intent(in) :: lchnk            ! index of grid chunk

    real(r8), intent(in) :: radius_max       ! upper bound of radius used for the calculation of deposition velocity [m]

    real(r8), intent(in) :: tair(pcols)      ! air temperature [K]
    real(r8), intent(in) :: pmid(pcols)      ! air pressure [Pa]

    real(r8), intent(in) ::  radius_part(pcols)   ! mean (volume or number) particle radius [m]
    real(r8), intent(in) :: density_part(pcols)   ! density of particle material [kg/m3]
    real(r8), intent(in) ::     sig_part(pcols)   ! geometric standard deviation of particle size distribution

    real(r8), intent(in) :: fricvel(pcols)   ! friction velocity [m/s]
    real(r8), intent(in) :: ram1(pcols)      ! aerodynamical resistance [s/m]

    real(r8), intent(in)  :: vlc_grv(pcols)  ! gravitational deposition velocity [m/s]
    real(r8), intent(out) :: vlc_trb(pcols)  ! turbulent dry deposition velocity [m/s]
    real(r8), intent(out) :: vlc_dry(pcols)  ! total     dry deposition velocity [m/s]

    !------------------------------------------------------------------------
    ! Local Variables

    integer  :: ii            ! grid column index
    real(r8) :: vsc_dyn_atm   ! Dynamic viscosity of air [kg m-1 s-1]
    real(r8) :: vsc_knm_atm   ! Kinematic viscosity of atmosphere [m2 s-1]
    real(r8) :: slp_crc       ! Slip correction factor [unitless]
    real(r8) :: radius_moment ! median radius for moment [m]

    real(r8) :: shm_nbr       ! Schmidt number [unitless]
    real(r8) :: stk_nbr       ! Stokes number [unitless]
    real(r8) :: rss_trb       ! Resistance to turbulent deposition [s m-1]
    real(r8) :: rss_lmn       ! Quasi-laminar layer resistance [s m-1]
    real(r8) :: brownian      ! collection efficiency for Browning diffusion
    real(r8) :: impaction     ! collection efficiency for impaction
    real(r8) :: interception  ! collection efficiency for interception
    real(r8) :: stickfrac     ! fraction of particles sticking to surface

    integer  :: lt              ! land type index
    real(r8) :: lnd_frc         ! land type fraction [unitless]
    real(r8) :: vlc_trb_ontype  ! turbulent dry dep. velocity on single land type [m/s]
    real(r8) :: vlc_trb_wgtsum  ! turbulent dry dep. velocity averaged over land types [m/s]
    real(r8) :: vlc_dry_wgtsum  ! total     dry dep. velocity averaged over land types [m/s]

    ! constants

    real(r8),parameter :: eps0 = 3.0_r8 ! empirical parameter $\varepsilon_0$ in Eq. (5) of Zhang L. et al. (2001)
    real(r8),parameter :: beta = 2.0_r8 ! empirical parameter $\beta$ in Eq. (7c) of Zhang L. et al. (2001)
    real(r8),parameter :: stickfrac_lowerbnd  = 1.0e-10_r8  ! lower bound of stick fraction

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

    !---------
    do ii=1,ncol

       ! Calculate size-INdependent thermokinetic properties of the air

       vsc_dyn_atm = air_dynamic_viscosity( tair(ii) )
       vsc_knm_atm = air_kinematic_viscosity( tair(ii), pmid(ii), rair )

       ! Calculate the mean radius and Schmidt number of the moment

       radius_moment = radius_for_moment( moment,sig_part(ii),radius_part(ii),radius_max )
       shm_nbr = schmidt_number( tair(ii), pmid(ii), radius_moment, rair, pi, boltz, vsc_dyn_atm, vsc_knm_atm )

       ! Initialize deposition velocities averages over different land surface types

       vlc_trb_wgtsum = 0._r8
       vlc_dry_wgtsum = 0._r8

       ! Loop over different land surface types. Calculate deposition velocities of 
       ! those different surface types. The overall deposition velocity of a grid cell 
       ! is the area-weighted average of those land-type-specific velocities.

       do lt = 1,n_land_type

          lnd_frc = fraction_landuse(ii,lt,lchnk)

          if ( lnd_frc /= 0._r8 ) then

             !----------------------------------------------------------------------
             ! Collection efficiency of deposition mechanism 1 - Brownian diffusion
             !----------------------------------------------------------------------
             brownian = shm_nbr**(-gamma(lt))

             !----------------------------------------------------------------------
             ! Collection efficiency of deposition mechanism 2 - interception 
             !----------------------------------------------------------------------
             if (radius_collector(lt) > 0.0_r8) then ! vegetated surface
                interception = 2.0_r8*(radius_moment/radius_collector(lt))**2.0_r8
             else ! non-vegetated surface
                interception = 0.0_r8
             endif

             !----------------------------------------------------------------------
             ! Collection efficiency of deposition mechanism 3 - impaction 
             !----------------------------------------------------------------------
             if (radius_collector(lt) > 0.0_r8) then ! vegetated surface
                stk_nbr = vlc_grv(ii) * fricvel(ii) / (gravit*radius_collector(lt))
             else ! non-vegetated surface
                stk_nbr = vlc_grv(ii) * fricvel(ii) * fricvel(ii) / (gravit*vsc_knm_atm)  ! SeP97 p.965
             endif

             impaction = (stk_nbr/(alpha(lt)+stk_nbr))**beta   ! Eq. (7c) of Zhang L. et al.  (2001)

             !-----------------------------------------------------
             ! Stick fraction, Eq. (10) of Zhang L. et al.  (2001)
             !-----------------------------------------------------
             if (iwet(lt) > 0) then
                stickfrac = 1.0_r8
             else
                stickfrac = max( stickfrac_lowerbnd, exp(-sqrt(stk_nbr)) )
             endif

             !----------------------------------------------------------------------------------
             ! Using the numbers calculated above, compute the quasi-laminar layer resistance
             ! following Zhang L. et al. (2001), Eq. (5)
             !----------------------------------------------------------------------------------
             rss_lmn = 1.0_r8 / (eps0 * fricvel(ii) * stickfrac * (brownian+interception+impaction))

             !--------------------------------------------------------------------
             ! Total resistence and deposition velocity of turbulent deposition,
             ! see Eq. (21) of Zender et al. (2003)
             !--------------------------------------------------------------------
             rss_trb = ram1(ii) + rss_lmn + ram1(ii)*rss_lmn*vlc_grv(ii)
             vlc_trb_ontype = 1.0_r8 / rss_trb

             !--------------------------------------------------------------------------------
             ! Contributions to the single-value bulk deposition velocities of the grid cell
             !--------------------------------------------------------------------------------
             vlc_trb_wgtsum = vlc_trb_wgtsum + lnd_frc*( vlc_trb_ontype )
             vlc_dry_wgtsum = vlc_dry_wgtsum + lnd_frc*( vlc_trb_ontype + vlc_grv(ii) )
          endif

       enddo  ! lt=1,n_land_type

       vlc_trb(ii) = vlc_trb_wgtsum
       vlc_dry(ii) = vlc_dry_wgtsum

    enddo ! ii=1,ncol
    !-------

  end subroutine modal_aero_turb_drydep_velocity

  !==============================================================================
  ! Calculate the radius for a moment of a lognormal size distribution
  !==============================================================================
  real(r8) function radius_for_moment( moment,sig_part,radius_part,radius_max)

    integer, intent(in) :: moment        ! moment of the particle size distribution. 0 = number; 3 = volume
    real(r8),intent(in) :: sig_part      ! geometric standard deviation of particle size distribution
    real(r8),intent(in) :: radius_part   ! mean (volume or number) particle radius [m]
    real(r8),intent(in) :: radius_max    ! developer-specified upper bound of mean radius [m]

    real(r8) :: lnsig
    
    lnsig = log(sig_part)
    radius_for_moment = min(radius_max,radius_part)*exp((float(moment)-1.5_r8)*lnsig*lnsig)

  end function radius_for_moment

  !==========================================================================
  ! Calculate dynamic viscosity of air, unit [kg m-1 s-1]. See RoY94 p. 102
  !==========================================================================
  real(r8) function air_dynamic_viscosity( temp )

    real(r8),intent(in) :: temp   ! air temperature [K]

    air_dynamic_viscosity = 1.72e-5_r8 * ( (temp/273.0_r8)**1.5_r8) * 393.0_r8 &
                                          /(temp+120.0_r8)

  end function air_dynamic_viscosity

  !==========================================================================
  ! Calculate kinematic viscosity of air, unit [m2 s-1]
  !==========================================================================
  real(r8) function air_kinematic_viscosity( temp, pres, rair )

    real(r8),intent(in) :: temp   ! air temperature [K]
    real(r8),intent(in) :: pres   ! air pressure [Pa]
    real(r8),intent(in) :: rair   ! gas constant of air [J/K/kg]

    real(r8) :: vsc_dyn_atm  ! dynamic viscosity of air [kg m-1 s-1]
    real(r8) :: rho          ! density of air [kg/m3]

    vsc_dyn_atm = air_dynamic_viscosity( temp)
    rho = pres/rair/temp

    air_kinematic_viscosity = vsc_dyn_atm/rho

  end function air_kinematic_viscosity

  !======================================================
  ! Slip correction factor [unitless]. 
  ! See, e.g., SeP97 p. 464 and Zhang L. et al. (2001),
  ! DOI: 10.1016/S1352-2310(00)00326-5, Eq. (3).
  !======================================================
  real(r8) function slip_correction_factor( dyn_visc, pres, temp, rair, pi, particle_radius ) 

    real(r8),intent(in) :: dyn_visc         ! dynamic viscosity of air [kg m-1 s-1]
    real(r8),intent(in) :: pres             ! air pressure [Pa]
    real(r8),intent(in) :: temp             ! air temperature [K]
    real(r8),intent(in) :: rair             ! gas constant of air [J/K/kg]
    real(r8),intent(in) :: pi               ! constant pi = 3.14159....
    real(r8),intent(in) :: particle_radius  ! particle radius [m]

    real(r8) :: mean_free_path  ! [m]

    mean_free_path = 2.0_r8*dyn_visc/( pres*sqrt( 8.0_r8/(pi*rair*temp) ) )

    slip_correction_factor = 1.0_r8 + mean_free_path * &
                             ( 1.257_r8+0.4_r8*exp(-1.1_r8*particle_radius/mean_free_path) ) / &
                             particle_radius

  end function slip_correction_factor

  !====================================================================
  ! Calculate the Schmidt number of air [unitless], see SeP97 p.972
  !====================================================================
  real(r8) function schmidt_number( temp, pres, radius, rair, pi, boltz, vsc_dyn_atm, vsc_knm_atm )

    real(r8),intent(in) :: temp             ! air temperature [K]
    real(r8),intent(in) :: pres             ! air pressure [Pa]
    real(r8),intent(in) :: radius           ! particle radius [m]
    real(r8),intent(in) :: rair             ! gas constant of air [J/K/kg]
    real(r8),intent(in) :: pi               ! constant pi = 3.14159....
    real(r8),intent(in) :: boltz            ! Boltzman's constant [J/K/molecule] 
    real(r8),intent(in) :: vsc_dyn_atm      ! dynamic viscosity of air [kg m-1 s-1]
    real(r8),intent(in) :: vsc_knm_atm      ! kinematic viscosity of air [m2 s-1]

    real(r8) :: slp_crc   ! slip correction factor [unitless]
    real(r8) :: dff_aer   ! Brownian diffusivity of particle [m2/s], see SeP97 p.474

    slp_crc = slip_correction_factor( vsc_dyn_atm, pres, temp, rair, pi, radius )
    dff_aer = boltz * temp * slp_crc / (6.0_r8*pi*vsc_dyn_atm*radius) 

    schmidt_number = vsc_knm_atm / dff_aer

  end function schmidt_number

  !=======================================================================================
  ! Calculate the bulk gravitational settling velocity [m/s-1] 
  !  - using the terminal velocity of sphere falling in a fluid based on Stokes's law and 
  !  - taking into account the influces of size distribution.
  !=======================================================================================
  real(r8) function gravit_settling_velocity( particle_radius, particle_density, gravit,       &
                                              slip_correction, dynamic_viscosity, particle_sig )

    real(r8),intent(in) :: particle_radius   ! [m]
    real(r8),intent(in) :: particle_density  ! [kg/m3]
    real(r8),intent(in) :: gravit            ! [kg/m2/s], gravitational acceleration
    real(r8),intent(in) :: slip_correction   ! [unitless]
    real(r8),intent(in) :: dynamic_viscosity ! [kg/m/s], dynamic viscosity of air
    real(r8),intent(in) :: particle_sig      ! geometric standard deviation of particle size distribution

    real(r8) :: lnsig         ! ln(particle_sig)
    real(r8) :: dispersion    ! accounts for influence of size dist dispersion on bulk settling velocity
                              ! assuming radius_part is number mode radius * exp(1.5 ln(sigma))

    ! Calculate terminal velocity following, e.g., 
    !  -  Seinfeld and Pandis (1997),  p. 466
    !  - Zhang L. et al. (2001), DOI: 10.1016/S1352-2310(00)00326-5, Eq. 2.

    gravit_settling_velocity = (4.0_r8/18.0_r8) * particle_radius*particle_radius* &
                               particle_density*gravit*slip_correction&
                               /dynamic_viscosity

    ! Account for size distribution (i.e., we are calculating the bulk velocity
    ! for a particle population instead of a single particle).

    lnsig = log(particle_sig)
    dispersion = exp(2._r8*lnsig*lnsig)

    gravit_settling_velocity = gravit_settling_velocity * dispersion

    !---------

  end function gravit_settling_velocity

end module modal_aero_drydep
