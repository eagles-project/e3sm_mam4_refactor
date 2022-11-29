
module modal_aero_convproc
!---------------------------------------------------------------------------------
! Purpose:
!
! CAM interface to aerosol/trace-gas convective cloud processing scheme
!
! currently these routines assume stratiform and convective clouds only interact 
! through the detrainment of convective cloudborne material into stratiform clouds
!
! thus the stratiform-cloudborne aerosols (in the qqcw array) are not processed 
! by the convective up/downdrafts, but are affected by the detrainment
!
! Author: R. C. Easter
!
!---------------------------------------------------------------------------------

   use shr_kind_mod,  only: r8=>shr_kind_r8
   use ppgrid,        only: pver, pcols, pverp
   use cam_history,   only: outfld, addfld, horiz_only, add_default
   use cam_logfile,   only: iulog
   use cam_abortutils,only: endrun
   use spmd_utils,    only: masterproc

   use physconst,     only: gravit, rair, rhoh2o,             &
                            spec_class_aerosol, spec_class_gas
   use physics_types, only: physics_state, physics_ptend
   use phys_control,  only: phys_getopts
   use constituents,  only: pcnst, cnst_name

   use wetdep, only: faer_resusp_vs_fprec_evap_mpln
   use ndrop,  only: activate_modal

   use mam_support,     only: min_max_bound
   use modal_aero_data, only: ntot_amode, nspec_amode,  &
                        lmassptr_amode, lmassptrcw_amode, lspectype_amode, &
                        numptr_amode,  numptrcw_amode, &
                        specdens_amode, spechygro, &
                        voltonumblo_amode, voltonumbhi_amode, &
                        mmtoo_prevap_resusp

   implicit none

   save
   private   ! Make default type private to the module

! Public methods
   public :: &
      ma_convproc_init,             &!
      ma_convproc_intr               !


! module data
   real(r8) :: hund_ovr_g  ! = 100.0_r8/gravit, calculated once and used frequently
   real(r8), parameter :: mbsth = 1.e-15 ! threshold below which we treat the mass fluxes as zero [mb/s]
   logical  :: convproc_do_gas, convproc_do_aer

!=========================================================================================
  contains



!=========================================================================================
subroutine ma_convproc_init

!----------------------------------------
! Purpose:  declare output fields, initialize variables needed by convection
!----------------------------------------

  implicit none

  logical :: history_aerosol      ! Output the MAM aerosol tendencies

! 
! Add history fields
!
    call phys_getopts( history_aerosol_out=history_aerosol, &
        convproc_do_aer_out = convproc_do_aer, &
        convproc_do_gas_out = convproc_do_gas  )

    call addfld(      'SH_MFUP_MAX', horiz_only, 'A', 'kg/m2', &
                      'Shallow conv. column-max updraft mass flux' )
    call addfld(      'SH_WCLDBASE', horiz_only, 'A', 'm/s', &
                      'Shallow conv. cloudbase vertical velocity' )
    call addfld(      'SH_KCLDBASE', horiz_only, 'A', '1', &
                      'Shallow conv. cloudbase level index' )

    call addfld(      'DP_MFUP_MAX', horiz_only, 'A', 'kg/m2', &
                      'Deep conv. column-max updraft mass flux' )
    call addfld(      'DP_WCLDBASE', horiz_only, 'A', 'm/s', &
                      'Deep conv. cloudbase vertical velocity' )
    call addfld(      'DP_KCLDBASE', horiz_only, 'A', '1', &
                      'Deep conv. cloudbase level index' )

    if ( history_aerosol .and. convproc_do_aer ) then
       call add_default( 'SH_MFUP_MAX', 1, ' ' )
       call add_default( 'SH_WCLDBASE', 1, ' ' )
       call add_default( 'SH_KCLDBASE', 1, ' ' )
       call add_default( 'DP_MFUP_MAX', 1, ' ' )
       call add_default( 'DP_WCLDBASE', 1, ' ' )
       call add_default( 'DP_KCLDBASE', 1, ' ' )
    endif

!
! Print control variable settings
!
   if ( masterproc ) then 
           write(*,'(a,l12)')     'ma_convproc_init - convproc_do_aer               = ', &
                 convproc_do_aer
           write(*,'(a,l12)')     'ma_convproc_init - convproc_do_gas               = ', &
                 convproc_do_gas
   endif

   return
end subroutine ma_convproc_init



!=========================================================================================
subroutine ma_convproc_intr( state, ztodt,                          & ! in
                           dp_frac, icwmrdp, rprddp, evapcdp,       & ! in
                           sh_frac, icwmrsh, rprdsh, evapcsh,       & ! in
                           dlf, dlfsh, cmfmcsh, sh_e_ed_ratio,      & ! in
                           nsrflx_mzaer2cnvpr, qsrflx_mzaer2cnvpr,  & ! in
                           mu, md, du, eu,ed, dp,                   & ! in
                           jt, maxg, ideep, lengath, species_class, & ! in
                           ptend, aerdepwetis                       ) ! inout
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Convective cloud processing (transport, activation/resuspension,
!    wet removal) of aerosols and trace gases.
!    (Currently no aqueous chemistry and no trace-gas wet removal)
! Does aerosols    when convproc_do_aer is .true.
! Does trace gases when convproc_do_gas is .true.
!
! Does deep and shallow convection
! Uses mass fluxes, cloud water, precip production from the
!    convective cloud routines
! 
! Author: R. Easter
! 
!-----------------------------------------------------------------------

 
! Arguments
   type(physics_state), intent(in )   :: state          ! Physics state variables
   type(physics_ptend), intent(inout) :: ptend          ! indivdual parameterization tendencies
   real(r8), intent(in)    :: ztodt                     ! 2 delta t (model time step, not sure why it is "2" delta t) [s]
   real(r8), intent(in)    :: dp_frac(pcols,pver)       ! Deep conv cloud frac [fraction]
   real(r8), intent(in)    :: icwmrdp(pcols,pver)       ! Deep conv cloud condensate (in cloud) [kg/kg]
   real(r8), intent(in)    :: rprddp(pcols,pver)        ! Deep conv precip production (grid avg) [kg/kg/s]
   real(r8), intent(in)    :: evapcdp(pcols,pver)       ! Deep conv precip evaporation (grid avg) [kg/kg/s]
   real(r8), intent(in)    :: sh_frac(pcols,pver)       ! Shal conv cloud frac [fraction]
   real(r8), intent(in)    :: icwmrsh(pcols,pver)       ! Shal conv cloud condensate (in cloud) [kg/kg]
   real(r8), intent(in)    :: rprdsh(pcols,pver)        ! Shal conv precip production (grid avg) [kg/kg/s]
   real(r8), intent(in)    :: evapcsh(pcols,pver)       ! Shal conv precip evaporation (grid avg) [kg/kg/s]
   real(r8), intent(in)    :: dlf(pcols,pver)           ! Tot  conv cldwtr detrainment (grid avg) [kg/kg/s]
   real(r8), intent(in)    :: dlfsh(pcols,pver)         ! Shal conv cldwtr detrainment (grid avg) [kg/kg/s]
   real(r8), intent(in)    :: cmfmcsh(pcols,pverp)      ! Shal conv mass flux [kg/m2/s]
   real(r8), intent(in)    :: sh_e_ed_ratio(pcols,pver) ! shallow conv [ent/(ent+det)] ratio [fraction]
   integer,  intent(in)    :: nsrflx_mzaer2cnvpr
   real(r8), intent(in)    :: qsrflx_mzaer2cnvpr(pcols,pcnst,nsrflx_mzaer2cnvpr)
   real(r8), intent(inout) :: aerdepwetis(pcols,pcnst)  ! aerosol wet deposition (interstitial) [kg/m2/s]

                                                ! mu, md, ..., ideep, lengath are all deep conv variables
                                                ! *** AND ARE GATHERED ***
   real(r8), intent(in)    :: mu(pcols,pver)    ! Updraft mass flux (positive) [mb/s]
   real(r8), intent(in)    :: md(pcols,pver)    ! Downdraft mass flux (negative) [mb/s]
                                                ! eu, ed, du are "d(massflux)/dp" and are all positive
   real(r8), intent(in)    :: eu(pcols,pver)    ! Mass entrain rate into updraft [1/s]
   real(r8), intent(in)    :: ed(pcols,pver)    ! Mass entrain rate into downdraft [1/s]
   real(r8), intent(in)    :: du(pcols,pver)    ! Mass detrain rate from updraft [1/s]
   real(r8), intent(in)    :: dp(pcols,pver)    ! Delta pressure between interfaces [mb]
   integer,  intent(in)    :: jt(pcols)         ! Index of cloud top (updraft top) for each column in w grid
   integer,  intent(in)    :: maxg(pcols)       ! Index of cloud base (level of maximum moist static energy) for each column in w grid
   integer,  intent(in)    :: ideep(pcols)      ! Gathering array
   integer,  intent(in)    :: lengath           ! Gathered min lon indices over which to operate
   integer,  intent(in)    :: species_class(:)  ! species index


! Local variables
   integer  :: ncol             ! total column number. from state%ncol
   integer  :: n, ll, l,lc, lchnk  ! indices
   logical  :: dotend(pcnst)    ! if do tendency

   real(r8) :: dlfdp(pcols,pver)                ! Deep (Total-Shallow) conv cldwtr detrainment (grid avg) [kg/kg/s] 
   real(r8) :: dqdt(pcols,pver,pcnst)           ! time tendency of q [kg/kg/s]
   real(r8) :: qnew(pcols,pver,pcnst)           ! tracer mixing ratio from state%q [kg/kg]
                                                ! qnew is updated through the processes in this subroutine but does not update into state
   real(r8) :: sflxic(pcols,pcnst)
   real(r8) :: sflxid(pcols,pcnst)
   real(r8) :: sflxec(pcols,pcnst)
   real(r8) :: sflxed(pcols,pcnst)

   integer, parameter :: nsrflx = 6             ! last dimension of qsrflx
   real(r8) :: qsrflx(pcols,pcnst,nsrflx)       ! process-specific column tracer tendencies
                            !  1 = activation   of interstial to  conv-cloudborne
                            !  2 = resuspension of conv-cloudborne to interstital
                            !  3 = aqueous chemistry (not implemented yet, so  zero)
                            !  4 = wet removal
                            !  5 = actual precip-evap resuspension (what actually is applied to a species)
                            !  6 = pseudo precip-evap resuspension (for history file)


!
! Initialize
!
  ncol  = state%ncol
  lchnk = state%lchnk
  qnew(1:ncol,:,:) = state%q(1:ncol,:,:)
  dotend(:) = ptend%lq(:)
  dqdt(:,:,:) = ptend%q(:,:,:)
  hund_ovr_g = 100.0_r8/gravit
!  used with zm_conv mass fluxes and delta-p. This is also used in other
!  subroutines in this file since it is declared at the beginning.
!     for mu = [mbar/s],   mu*hund_ovr_g = [kg/m2/s]
!     for dp = [mbar] and q = [kg/kg],   q*dp*hund_ovr_g = [kg/m2]

!
! prepare for processing
!
  qsrflx(:,:,:) = 0.0_r8
  call update_qnew_ptend(                                         &
                         dotend,                    .false.,      &  ! in
                         ncol,     species_class,   dqdt,         &  ! in
                         qsrflx,   ztodt,                         &  ! in
                         ptend,    qnew                           )  ! inout

  ! calculate variables for output
   sflxic(:,:) = 0.0_r8
   sflxid(:,:) = 0.0_r8
   sflxec(:,:) = 0.0_r8
   sflxed(:,:) = 0.0_r8
   do l = 1, pcnst
      if ( (species_class(l) == spec_class_aerosol) .and. ptend%lq(l) ) then
         sflxec(1:ncol,l) = qsrflx_mzaer2cnvpr(1:ncol,l,1)
         sflxed(1:ncol,l) = qsrflx_mzaer2cnvpr(1:ncol,l,2)
      endif
   enddo

  if (convproc_do_aer .or. convproc_do_gas) then
     !
     ! do deep conv processing
     !
     dqdt(:,:,:) = 0.0_r8
     qsrflx(:,:,:) = 0.0_r8
     dlfdp(1:ncol,:) = max( (dlf(1:ncol,:) - dlfsh(1:ncol,:)), 0.0_r8 )
     call ma_convproc_dp_intr(                    &
        state, ztodt,                             & ! in
        dp_frac, icwmrdp, rprddp, evapcdp, dlfdp, & ! in
        mu, md, du, eu, ed, dp,                   & ! in
        jt, maxg, ideep, lengath, qnew,           & ! in
        nsrflx,  species_class,                   & ! in
        dqdt, qsrflx,                             & ! inout
        dotend                                    ) ! out

     ! apply deep conv processing tendency and prepare for shallow conv processing
     call update_qnew_ptend(                       &
            dotend,                  .true.,       &  ! in
            ncol,   species_class,   dqdt,         &  ! in
            qsrflx, ztodt,                         &  ! in
            ptend,  qnew                           )  ! inout

     ! update variables for output
     do l = 1, pcnst
        if ( .not. dotend(l) ) cycle
        if ((species_class(l) == spec_class_aerosol) .or. &
            (species_class(l) == spec_class_gas    )) then
           ! these used for history file wetdep diagnostics
           sflxic(1:ncol,l) = sflxic(1:ncol,l) + qsrflx(1:ncol,l,4)
           sflxid(1:ncol,l) = sflxid(1:ncol,l) + qsrflx(1:ncol,l,4)
           sflxec(1:ncol,l) = sflxec(1:ncol,l) + qsrflx(1:ncol,l,6)
           sflxed(1:ncol,l) = sflxed(1:ncol,l) + qsrflx(1:ncol,l,6)
        endif

        if (species_class(l) == spec_class_aerosol) then
           ! this used for surface coupling
           aerdepwetis(1:ncol,l) = aerdepwetis(1:ncol,l) &
                + qsrflx(1:ncol,l,4) + qsrflx(1:ncol,l,5)
        endif
     enddo

     !
     ! do shallow conv processing
     !
     dqdt(:,:,:) = 0.0_r8
     qsrflx(:,:,:) = 0.0_r8
     call ma_convproc_sh_intr(                    &
        state, ztodt,                             & ! in
        sh_frac, icwmrsh, rprdsh, evapcsh, dlfsh, & ! in
        cmfmcsh, sh_e_ed_ratio,   qnew,           & ! in
        nsrflx,  species_class,                   & ! in
        dqdt,    qsrflx,                          & ! inout
        dotend                                    ) ! out

     ! apply shallow conv processing tendency
     call update_qnew_ptend(                         &
              dotend,                  .true.,       &  ! in
              ncol,   species_class,   dqdt,         &  ! in
              qsrflx, ztodt,                         &  ! in
              ptend,  qnew                           )  ! inout

     ! update variables for output
     do l = 1, pcnst
        if ( .not. dotend(l) ) cycle
        if ((species_class(l) == spec_class_aerosol) .or. &
            (species_class(l) == spec_class_gas    )) then
           sflxic(1:ncol,l) = sflxic(1:ncol,l) + qsrflx(1:ncol,l,4)
           sflxec(1:ncol,l) = sflxec(1:ncol,l) + qsrflx(1:ncol,l,6)
        endif

        if (species_class(l) == spec_class_aerosol) then
           ! this used for surface coupling
           aerdepwetis(1:ncol,l) = aerdepwetis(1:ncol,l) &
                + qsrflx(1:ncol,l,4) + qsrflx(1:ncol,l,5)
        endif
     enddo


  endif ! (convproc_do_aer  .or. convproc_do_gas) then


! output wet deposition fields to history
!    I = in-cloud removal;     E = precip-evap resuspension
!    C = convective (total);   D = deep convective
! note that the precip-evap resuspension includes that resulting from
!    below-cloud removal, calculated in mz_aero_wet_intr
  if (convproc_do_aer) then
     do n = 1, ntot_amode
     do ll = 0, nspec_amode(n)

        call assign_la_lc( n,   ll,   l,   lc   )

        call outfld( trim(cnst_name(l))//'SFWET', aerdepwetis(:,l), pcols, lchnk)
        call outfld( trim(cnst_name(l))//'SFSIC', sflxic(:,l), pcols, lchnk )
        call outfld( trim(cnst_name(l))//'SFSEC', sflxec(:,l), pcols, lchnk )
        call outfld( trim(cnst_name(l))//'SFSID', sflxid(:,l), pcols, lchnk )
        call outfld( trim(cnst_name(l))//'SFSED', sflxed(:,l), pcols, lchnk )
     enddo ! ll
     enddo ! n
  endif

end subroutine ma_convproc_intr

!=========================================================================================
subroutine update_qnew_ptend(                                         &
                           dotend, is_update_ptend,                   &  ! in
                           ncol,   species_class,   dqdt,             &  ! in
                           qsrflx, ztodt,                             &  ! in
                           ptend,  qnew                               )  ! inout
! ---------------------------------------------------------------------------------------
! update qnew, ptend (%q and %lq)
! ---------------------------------------------------------------------------------------  
use physics_types, only: physics_ptend
use constituents,  only: pcnst

   ! Arguments
   type(physics_ptend), intent(inout) :: ptend          ! indivdual parameterization tendencies. ptend%q [kg/kg/s] and ptend%lq [logical] will be updated
   logical,  intent(in)    :: dotend(pcnst)             ! if do tendency
   logical,  intent(in)    :: is_update_ptend           ! if add dqdt onto ptend%q
   integer,  intent(in)    :: ncol                      ! index
   integer,  intent(in)    :: species_class(:)          ! species index
   real(r8), intent(in)    :: dqdt(pcols,pver,pcnst)    ! time tendency of tracer [kg/kg/s]
   real(r8), intent(in)    :: qsrflx(:,:,:)             ! process-specific column tracer tendencies. see ma_convproc_tend for detail info [kg/m2/s]
   real(r8), intent(in)    :: ztodt                     ! 2 delta t (model time step, not sure why it is "2" delta t) [s]
   real(r8), intent(inout) :: qnew(pcols,pver,pcnst)    ! Tracer array including moisture [kg/kg]

   ! Local variables
   integer  :: ll                         ! index
   real(r8) :: qtmp(pcols,pver,pcnst)     ! temporary q [kg/kg]


   do ll = 1, pcnst
     if ( .not. dotend(ll) ) cycle

     ! calc new q (after ma_convproc_sh_intr)
     qtmp(1:ncol,:,ll) = qnew(1:ncol,:,ll) + ztodt*dqdt(1:ncol,:,ll)
     qnew(1:ncol,:,ll) = max( 0.0_r8, qtmp(1:ncol,:,ll) )

     ! add dqdt onto ptend%q and set ptend%lq
     if ( is_update_ptend ) then
        ptend%q(1:ncol,:,ll) = ptend%q(1:ncol,:,ll) + dqdt(1:ncol,:,ll)
        ptend%lq(ll) = .true.
     endif

   enddo ! ll
end subroutine update_qnew_ptend

!=========================================================================================
subroutine ma_convproc_dp_intr(                &
     state,  dt,                               & ! in
     dp_frac, icwmrdp, rprddp, evapcdp, dlfdp, & ! in
     mu, md, du, eu, ed, dp,                   & ! in
     jt, maxg, ideep, lengath,   qnew,         & ! in
     nsrflx,  species_class,                   & ! in
     dqdt, qsrflx,                             & ! inout
     dotend                                    ) ! out
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Convective cloud processing (transport, activation/resuspension,
!    wet removal) of aerosols and trace gases.
!    (Currently no aqueous chemistry and no trace-gas wet removal)
! Does aerosols    when convproc_do_aer is .true.
! Does trace gases when convproc_do_gas is .true.
!
! This routine does deep convection
! Uses mass fluxes, cloud water, precip production from the
!    convective cloud routines
! 
! Author: R. Easter
! 
!-----------------------------------------------------------------------

 
! Arguments
   type(physics_state), intent(in ) :: state          ! Physics state variables

   real(r8), intent(in)    :: dt                         ! delta t (model time increment) [s]

   real(r8), intent(in)    :: qnew(pcols,pver,pcnst)     ! tracer mixing ratio including water vapor [kg/kg]
   real(r8), intent(inout) :: dqdt(pcols,pver,pcnst)     ! time tendency of q [kg/kg/s]
   logical,  intent(out)   :: dotend(pcnst)              ! if do tendency
   integer,  intent(in)    :: nsrflx                     ! last dimension of qsrflx
   real(r8), intent(inout) :: qsrflx(pcols,pcnst,nsrflx) ! process-specific column tracer tendencies (see ma_convproc_intr for more information)

   real(r8), intent(in)    :: dp_frac(pcols,pver) ! Deep conv cloud fraction [0-1]
   real(r8), intent(in)    :: icwmrdp(pcols,pver) ! Deep conv cloud condensate (in cloud) [kg/kg]
   real(r8), intent(in)    :: rprddp(pcols,pver)  ! Deep conv precip production (grid avg) [kg/kg/s]
   real(r8), intent(in)    :: evapcdp(pcols,pver) ! Deep conv precip evaporation (grid avg) [kg/kg/s]
   real(r8), intent(in)    :: dlfdp(pcols,pver)   ! Deep conv cldwtr detrainment (grid avg) [kg/kg/s]

                                               ! mu, md, ..., ideep, lengath are all deep conv variables
   real(r8), intent(in)    :: mu(pcols,pver)   ! Updraft mass flux (positive) [mb/s]
   real(r8), intent(in)    :: md(pcols,pver)   ! Downdraft mass flux (negative) [mb/s]
   real(r8), intent(in)    :: du(pcols,pver)   ! Mass detrain rate from updraft [1/s]
   real(r8), intent(in)    :: eu(pcols,pver)   ! Mass entrain rate into updraft [1/s]
   real(r8), intent(in)    :: ed(pcols,pver)   ! Mass entrain rate into downdraft [1/s]
                           ! eu, ed, du are "d(massflux)/dp" and are all positive
   real(r8), intent(in)    :: dp(pcols,pver)   ! Delta pressure between interfaces [mb]

   integer,  intent(in)    :: jt(pcols)         ! Index of cloud top for each column
   integer,  intent(in)    :: maxg(pcols)       ! Index of cloud top for each column
   integer,  intent(in)    :: ideep(pcols)      ! Gathering array
   integer,  intent(in)    :: lengath           ! Gathered min lon indices over which to operate
   integer,  intent(in)    :: species_class(:)  ! species index

! Local variables
   integer :: ii

   real(r8) :: dpdry(pcols,pver)     ! layer delta-p-dry [mb]
   ! diagnostics variables to write out. 
   ! xx_kcldbase is filled with index kk. maybe better define as integer
   ! keep it in C++ refactoring for BFB comparison (by Shuaiqi Tang, 2022)
   real(r8) :: xx_mfup_max(pcols), xx_wcldbase(pcols), xx_kcldbase(pcols)

!
! Initialize
!


! initialize dpdry (units=mb), which is used for tracers of dry mixing ratio type
   dpdry = 0._r8
   do ii = 1, lengath
      dpdry(ii,:) = state%pdeldry(ideep(ii),:)/100._r8
   end do

! turn on/off calculations for aerosols and trace gases
   call assign_dotend( species_class, dotend)

!
! do ma_convproc_tend call
!
! question/issue - when computing first-order removal rate for convective cloud water,
!    should dlf be included as is done in wetdepa?
! detrainment does not change the in-cloud (= in updraft) cldwtr mixing ratio
! when you have detrainment, the updraft air mass flux is decreasing with height,
!    and the cldwtr flux may be decreasing also, 
!    but the in-cloud cldwtr mixing ratio is not changed by detrainment itself
! this suggests that wetdepa is incorrect, and dlf should not be included
!
! if dlf should be included, then you want to calculate
!    rprddp / (dp_frac*icwmrdp + dt*(rprddp + dlfdp)]
! so need to pass both rprddp and dlfdp to ma_convproc_tend
!

   call ma_convproc_tend(                                            &
                     'deep',                                         &
                     state%lchnk,      pcnst,            dt,         &
                     state%t,    state%pmid,             qnew,       &   
                     mu,         md,         du,         eu,         &   
                     ed,         dp,         dpdry,      jt,         &   
                     maxg,       ideep,      1,          lengath,    &       
                     dp_frac,    icwmrdp,    rprddp,     evapcdp,    &
                     dqdt,                                           & ! out
                     dotend,     nsrflx,                             &
                     qsrflx,                                         & ! out
                     species_class,                                  &
                     xx_mfup_max, xx_wcldbase, xx_kcldbase           ) ! out


    ! output diagnostics fields
    call outfld( 'DP_MFUP_MAX', xx_mfup_max, pcols, state%lchnk )
    call outfld( 'DP_WCLDBASE', xx_wcldbase, pcols, state%lchnk )
    call outfld( 'DP_KCLDBASE', xx_kcldbase, pcols, state%lchnk )

end subroutine ma_convproc_dp_intr

!=========================================================================================
subroutine ma_convproc_sh_intr(                 &
     state, dt,                                 & ! in
     sh_frac, icwmrsh, rprdsh, evapcsh, dlfsh,  & ! in
     cmfmcsh, sh_e_ed_ratio,    qnew,           & ! in
     nsrflx,  species_class,                    & ! in
     dqdt,    qsrflx,                           & ! inout
     dotend                                     ) ! out
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Convective cloud processing (transport, activation/resuspension,
!    wet removal) of aerosols and trace gases.
!    (Currently no aqueous chemistry and no trace-gas wet removal)
! Does aerosols    when convproc_do_aer is .true.
! Does trace gases when convproc_do_gas is .true.
!
! This routine does shallow convection
! Uses mass fluxes, cloud water, precip production from the
!    convective cloud routines
! 
! Author: R. Easter
! 
!-----------------------------------------------------------------------

! Arguments
   type(physics_state), intent(in ) :: state          ! Physics state variables
   real(r8), intent(in)    :: dt                      ! delta t (model time increment) [s]
   real(r8), intent(in)    :: qnew(pcols,pver,pcnst)  ! tracer mixing ratio (TMR) including water vapor [kg/kg]
   real(r8), intent(inout) :: dqdt(pcols,pver,pcnst)  ! time tendency of TMR [kg/kg/s]
   logical,  intent(out)   :: dotend(pcnst)           ! flag if do tendency
   integer,  intent(in)    :: nsrflx                  ! last dimension of qsrflx
   real(r8), intent(inout) :: qsrflx(pcols,pcnst,nsrflx)  ! process-specific column tracer tendencies  (see ma_convproc_intr for more information)

   real(r8), intent(in)    :: sh_frac(pcols,pver)       ! Shallow conv cloud frac [0-1]
   real(r8), intent(in)    :: icwmrsh(pcols,pver)       ! Shallow conv cloud condensate (in cloud) [kg/kg]
   real(r8), intent(in)    :: rprdsh(pcols,pver)        ! Shallow conv precip production (grid avg) [kg/kg/s]
   real(r8), intent(in)    :: evapcsh(pcols,pver)       ! Shallow conv precip evaporation (grid avg) [kg/kg/s]
   real(r8), intent(in)    :: dlfsh(pcols,pver)         ! Shallow conv cldwtr detrainment (grid avg) [kg/kg/s]
   real(r8), intent(in)    :: cmfmcsh(pcols,pverp)      ! Shallow conv mass flux [kg/m2/s]
   real(r8), intent(in)    :: sh_e_ed_ratio(pcols,pver) ! shallow conv [ent/(ent+det)] ratio
   integer,  intent(in)    :: species_class(:)          ! species index

! Local variables
   integer :: icol, ncol

   real(r8) :: dpdry(pcols,pver)     ! layer delta-p-dry [mb]
   real(r8) :: xx_mfup_max(pcols), xx_wcldbase(pcols), xx_kcldbase(pcols)  ! output of ma_convproc_tend, may not used

! variables that mimic the zm-deep counterparts
                                               ! mu, md, ..., ideep, lengath are all deep conv variables
   real(r8)  :: mu(pcols,pver)   ! Updraft mass flux (positive) [mb/s]
   real(r8)  :: md(pcols,pver)   ! Downdraft mass flux (negative) [mb/s]
        ! eu, ed, du are "d(massflux)/dp" and are all positive
   real(r8)  :: du(pcols,pver)   ! Mass detrain rate from updraft [1/s]
   real(r8)  :: eu(pcols,pver)   ! Mass entrain rate into updraft [1/s]
   real(r8)  :: ed(pcols,pver)   ! Mass entrain rate into downdraft [1/s]
   real(r8)  :: dp(pcols,pver)   ! Delta pressure between interfaces [mb]
   integer   :: jt(pcols)         ! Index of cloud top for each column
   integer   :: maxg(pcols)       ! Index of cloud bot for each column
   integer   :: ideep(pcols)      ! Gathering array [index]
   integer   :: lengath           ! Gathered min lon indices over which to operate


!
! Initialize
!

   ncol  = state%ncol

   ! md and ed are assumed zero in shallow convection in ma_convproc_tend
   md(:,:) = 0.0_r8
   ed(:,:) = 0.0_r8

! these dp and dpdry have units of mb
   dpdry(1:ncol,:) = state%pdeldry(1:ncol,:)/100._r8
   dp(   1:ncol,:) = state%pdel(   1:ncol,:)/100._r8

   ideep(:) = -1
   do icol = 1, ncol
      ideep(icol) = icol
   enddo

! mimic variables counterparts as in zm-deep
   call mimic_deep_counterparts( ncol,                   & ! in
                        dpdry,  cmfmcsh,  sh_e_ed_ratio, & ! in
                        jt, maxg, mu, eu, du, lengath    ) ! out


! turn on/off calculations for aerosols and trace gases
   call assign_dotend( species_class, dotend)

!
! do ma_convproc_tend call
!
! question/issue - when computing first-order removal rate for convective cloud water,
!    should dlf be included as is done in wetdepa?
! detrainment does not change the in-cloud (= in updraft) cldwtr mixing ratio
! when you have detrainment, the updraft air mass flux is decreasing with height,
!    and the cldwtr flux may be decreasing also, 
!    but the cldwtr mixing ratio does not change
! this suggests that wetdepa is incorrect, and dlf should not be included
!
! if dlf should be included, then you want to calculate
!    rprddp / (dp_frac*icwmrdp + dt*(rprddp + dlfdp)]
! so need to pass both rprddp and dlfdp to ma_convproc_tend
!

   call ma_convproc_tend(                                            &
                     'uwsh',                                         &
                     state%lchnk,      pcnst,            dt,         &
                     state%t,    state%pmid,             qnew,       &   
                     mu,         md,         du,         eu,         &   
                     ed,         dp,         dpdry,      jt,         &   
                     maxg,       ideep,      1,          lengath,    &       
                     sh_frac,    icwmrsh,    rprdsh,     evapcsh,    &
                     dqdt,                                           & ! out
                     dotend,     nsrflx,                             &
                     qsrflx,                                         & ! out
                     species_class,                                  &
                     xx_mfup_max, xx_wcldbase, xx_kcldbase           ) ! out

    ! output diagnostics fields
    call outfld( 'SH_MFUP_MAX', xx_mfup_max, pcols, state%lchnk )
    call outfld( 'SH_WCLDBASE', xx_wcldbase, pcols, state%lchnk )
    call outfld( 'SH_KCLDBASE', xx_kcldbase, pcols, state%lchnk )

end subroutine ma_convproc_sh_intr

!=========================================================================================
subroutine assign_dotend( species_class,  & ! in
                          dotend          ) ! out
!---------------------------------------------------------------------
! assign do-tendency flag from species_class, convproc_do_aer and convproc_do_gas.
! convproc_do_aer and convproc_do_gas are assigned in the beginning of the  module
!---------------------------------------------------------------------

use constituents,   only: pcnst

   integer,  intent(in)    :: species_class(:)
   logical,  intent(out)   :: dotend(pcnst)

   integer  :: ll

! turn on/off calculations for aerosols and trace gases
   do ll = 1, pcnst
      if (species_class(ll) == spec_class_aerosol .and. convproc_do_aer) then
         dotend(ll) = .true.
      elseif (species_class(ll) == spec_class_gas .and. convproc_do_gas) then
         dotend(ll) = .true.
      else
         dotend(ll) = .false.
      endif
   enddo

end subroutine assign_dotend

!=========================================================================================
subroutine mimic_deep_counterparts( ncol,                & ! in
                        dpdry,  cmfmcsh,  sh_e_ed_ratio, & ! in
                        jt, maxg, mu, eu, du, lengath    ) ! out
!-----------------------------------------------------------------------------
! create mass flux, entrainment, detrainment, and delta-p arrays
! with the same units as the zm-deep
!-----------------------------------------------------------------------------

   integer,  intent(in)    :: ncol              ! total number of column
   real(r8), intent(in)    :: dpdry(pcols,pver) ! layer delta-p-dry [mb]
   real(r8), intent(in)    :: cmfmcsh(pcols,pverp) ! Shallow conv mass flux [kg/m2/s]
   real(r8), intent(in)    :: sh_e_ed_ratio(pcols,pver)  ! shallow conv [ent/(ent+det)] ratio

   integer,  intent(out)   :: jt(pcols)        ! Index of cloud top for each column
   integer,  intent(out)   :: maxg(pcols)      ! Index of cloud bot for each column
   integer,  intent(out)   :: lengath          ! min lon indices over which to operate
   real(r8), intent(out)   :: mu(pcols,pver)   ! Updraft mass flux (positive) [mb/s]
   real(r8), intent(out)   :: eu(pcols,pver)   ! Mass entrain rate into updraft [1/s]
   real(r8), intent(out)   :: du(pcols,pver)   ! Mass detrain rate from updraft [1/s]


   integer      :: icol                 ! index of column
   integer      :: tot_conv_layer       ! total layers of convection in this column
   integer      :: maxg_init            ! initial maxg at icol


   lengath = ncol   ! mimic lengath in zm-deep

   do icol = 1, ncol
        ! load updraft mass flux from cmfmcsh
        call load_updraft_massflux(             &
                icol,   cmfmcsh,                & ! in
                tot_conv_layer,   jt, maxg, mu  ) ! out

        if (tot_conv_layer <= 0) cycle  ! current column has no convection

        ! extend below-cloud source region downwards
        call extend_belowcloud_downward(           &
                        icol,     dpdry,           & ! in
                        maxg,     mu,              & ! inout
                        maxg_init                  ) ! out

! calc ent / detrainment, using the [ent/(ent+det)] ratio from uw scheme
!    which is equal to [fer_out/(fer_out+fdr_out)]  (see uwshcu.F90)
! note that the ratio is set to -1.0 (invalid) when both fer and fdr are very small
!    and the ratio values are often strange (??) at topmost layer\
        call calculate_ent_det(                 &
                icol,   maxg,  maxg_init,  jt,  & ! in
                dpdry,  mu,    sh_e_ed_ratio,   & ! in
                eu,     du                      ) ! out

   enddo ! icol

end subroutine mimic_deep_counterparts


!=========================================================================================
subroutine load_updraft_massflux(               &
                icol,   cmfmcsh,                & ! in
                tot_conv_layer,   jt, maxg, mu  ) ! out
!-----------------------------------------------------------------------------
! load updraft mass flux from cmfmcsh
!-----------------------------------------------------------------------------

   real(r8), intent(in)    :: cmfmcsh(pcols,pverp) ! Shallow conv mass flux [kg/m2/s]
   integer,  intent(in)    :: icol             ! index of column
   integer,  intent(out)   :: tot_conv_layer   ! total layers of convection in this column
   integer,  intent(out)   :: jt(pcols)        ! Index of cloud top for each column
   integer,  intent(out)   :: maxg(pcols)      ! Index of cloud bot for each column
   real(r8), intent(out)   :: mu(pcols,pver)   ! Updraft mass flux (positive) [mb/s]

   integer              :: kk
   real(r8), parameter :: small_massflux = 1.0e-7_r8  ! if mass-flux < 1e-7 kg/m2/s ~= 1e-7 m/s ~= 1 cm/day, treat as zero


   ! initiate variables
   tot_conv_layer = 0 ! total layers of convection in this column
   mu(:,:) = 0.0_r8
   jt(:) = -1
   maxg(:) = -1

   do kk = 2, pver
       if (cmfmcsh(icol,kk) >= small_massflux) then
            ! mu has units of mb/s
            mu(icol,kk) = cmfmcsh(icol,kk) / hund_ovr_g
            tot_conv_layer = tot_conv_layer + 1
            if (tot_conv_layer == 1) jt(icol) = kk - 1
            maxg(icol) = kk
       endif
   enddo ! kk

end subroutine load_updraft_massflux

!=========================================================================================
subroutine extend_belowcloud_downward(             &
                        icol,     dpdry,           & ! in
                        maxg,     mu,              & ! inout
                        maxg_init                  ) ! out
!-----------------------------------------------------------------------------
! extend below-cloud source region downwards
!-----------------------------------------------------------------------------
   integer,  intent(in)    :: icol              ! index of column
   real(r8), intent(in)    :: dpdry(pcols,pver) ! layer delta-p-dry [mb]
   integer,  intent(inout) :: maxg(pcols)       ! Index of cloud bot for each column
   real(r8), intent(inout) :: mu(pcols,pver)    ! Updraft mass flux (positive) [mb/s]
   integer,  intent(out)   :: maxg_init         ! initial maxg at icol

   integer              :: kk, k_cldbot, k_4_below_cldbot         ! vertical index
   real(r8)             :: dp_sum               ! sum of dpdry for weighting purpose
   integer              :: maxg_minval

   ! initiate variables
   maxg_minval = pver*2  ! this variable seems not used
   maxg_minval = min( maxg_minval, maxg(icol) )

   k_cldbot = maxg(icol)          ! cloud bot level
   k_4_below_cldbot = min( k_cldbot+4, pver )  ! 4 levels below cloud bottom
   if (k_4_below_cldbot > k_cldbot) then  ! make sure cloud bot is not the bottom model level
      dp_sum = sum( dpdry(icol,k_cldbot:k_4_below_cldbot) )
      do kk = k_cldbot+1, k_4_below_cldbot
         ! extend mass flux below cloud to k_4_below_cldbot
         mu(icol,kk) = mu(icol,k_cldbot)*sum( dpdry(icol,kk:k_4_below_cldbot) )/dp_sum
      enddo ! kk
      maxg(icol) = k_4_below_cldbot
   endif

   ! assign initial maxg value at icol for calculate_ent_det use
   maxg_init = k_cldbot

end subroutine extend_belowcloud_downward

!=========================================================================================
subroutine calculate_ent_det(                   &
                icol,   maxg,  maxg_init,  jt,  & ! in
                dpdry,  mu,    sh_e_ed_ratio,   & ! in
                eu,     du                      ) ! out
!-----------------------------------------------------------------------------
! calc ent / detrainment, using the [ent/(ent+det)] ratio from uw scheme
!    which is equal to [fer_out/(fer_out+fdr_out)]  (see uwshcu.F90)
!
! note that the ratio is set to -1.0 (invalid) when both fer and fdr are very
! small and the ratio values are often strange (??) at topmost layer
!-----------------------------------------------------------------------------

   integer,  intent(in)    :: icol              ! index of column
   integer,  intent(in)    :: maxg(pcols)       ! Index of cloud bot for each column
   integer,  intent(in)    :: jt(pcols)         ! Index of cloud top for each column
   integer,  intent(in)    :: maxg_init         ! initial maxg at icol
   real(r8), intent(in)    :: dpdry(pcols,pver) ! layer delta-p-dry [mb]
   real(r8), intent(in)    :: mu(pcols,pver)    ! Updraft mass flux (positive) [mb/s]
   real(r8), intent(in)    :: sh_e_ed_ratio(pcols,pver)  ! shallow conv [ent/(ent+det)] ratio
   real(r8), intent(out)   :: eu(pcols,pver)   ! Mass entrain rate into updraft [1/s]
   real(r8), intent(out)   :: du(pcols,pver)   ! Mass detrain rate from updraft [1/s]

   integer              :: kk
   real(r8)             :: tmp_ratio            ! sh_e_ed_ratio(icol,kk) [fraction]
   real(r8)             :: tmp_mu_rate          ! ent/det rate from mu/dpdry [1/s]

   ! initiate variables
   du(:,:) = 0.0_r8
   eu(:,:) = 0.0_r8

   do kk = jt(icol), maxg(icol)
      if (kk < pver) then
         tmp_mu_rate = (mu(icol,kk) - mu(icol,kk+1))/dpdry(icol,kk)
      else
         tmp_mu_rate = mu(icol,kk)/dpdry(icol,kk)
      endif
      tmp_ratio = sh_e_ed_ratio(icol,kk)

      if (tmp_ratio < -1.0e-5_r8) then ! do ent only or det only
         if (tmp_mu_rate >= 0.0_r8) then
            eu(icol,kk) = tmp_mu_rate  ! net entrainment
         else
            du(icol,kk) = -tmp_mu_rate ! net detrainment
         endif

      else   ! do both ent and det
         if (tmp_mu_rate >= 0.0_r8) then
            ! net entrainment
            if (kk >= maxg_init .or. tmp_ratio < 0.0_r8) then
               ! layers at/below initial maxg (cloud base), or sh_e_ed_ratio is invalid
               eu(icol,kk) = tmp_mu_rate
            else
               tmp_ratio = max( tmp_ratio, 0.571_r8 ) ! not sure why 0.571 is used
               eu(icol,kk) = tmp_mu_rate*(tmp_ratio/(2.0_r8*tmp_ratio - 1.0_r8))
               du(icol,kk) = eu(icol,kk) - tmp_mu_rate
            endif
         else
            ! net detrainment
            if (kk <= jt(icol) .or. tmp_ratio < 0.0_r8) then
               ! layers at/above jt (cloud top), or sh_e_ed_ratio is invalid
               du(icol,kk) = - tmp_mu_rate
            else
               tmp_ratio = min( tmp_ratio, 0.429_r8 )
               du(icol,kk) = - tmp_mu_rate*(1.0_r8 - tmp_ratio)/(1.0_r8 - 2.0_r8*tmp_ratio)
               eu(icol,kk) = du(icol,kk) + tmp_mu_rate
            endif
         endif
      endif
   enddo ! kk

end subroutine calculate_ent_det

!=========================================================================================
subroutine ma_convproc_tend(                                         &
                     convtype,                                       & ! in
                     lchnk,      ncnst,      dt,                     & ! in
                     t,          pmid,       q,                      & ! in
                     mu,         md,         du,         eu,         & ! in
                     ed,         dp,         dpdry,      jt,         & ! in
                     mx,         ideep,      il1g,       il2g,       & ! in  
                     cldfrac,    icwmr,      rprd,       evapc,      & ! in
                     dqdt,                                           & ! out
                     doconvproc, nsrflx,                             & ! in
                     qsrflx,                                         & ! out
                     species_class,                                  & ! in
                     xx_mfup_max, xx_wcldbase, xx_kcldbase           ) ! out

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Convective transport of trace species.
! The trace species need not be conservative, and source/sink terms for
!    activation, resuspension, aqueous chemistry and gas uptake, and
!    wet removal are all applied.
! Currently this works with the ZM deep convection, but we should be able
!    to adapt it for both Hack and McCaa shallow convection
!
!
! Compare to subr convproc which does conservative trace species.
!
! A distinction between "moist" and "dry" mixing ratios is not currently made.
! (P. Rasch comment:  Note that we are still assuming that the tracers are 
!  in a moist mixing ratio this will change soon)

! 
! Method: 
! Computes tracer mixing ratios in updraft and downdraft "cells" in a
! Lagrangian manner, with source/sinks applied in the updraft other.
! Then computes grid-cell-mean tendencies by considering
!    updraft and downdraft fluxes across layer boundaries
!    environment subsidence/lifting fluxes across layer boundaries
!    sources and sinks in the updraft
!    resuspension of activated species in the grid-cell as a whole
!
! Note1:  A better estimate or calculation of either the updraft velocity
!         or fractional area is needed.
! Note2:  If updraft area is a small fraction of over cloud area, 
!         then aqueous chemistry is underestimated.  These are both
!         research areas.
! 
! Authors: O. Seland and R. Easter, based on convtran by P. Rasch
! 
!-----------------------------------------------------------------------

   implicit none

!-----------------------------------------------------------------------
!
! Input arguments
!
   character(len=*), intent(in) :: convtype  ! identifies the type of
                                             ! convection ("deep", "shcu")
   integer,  intent(in) :: lchnk             ! chunk identifier
   integer,  intent(in) :: ncnst             ! number of tracers to transport
   real(r8), intent(in) :: dt                ! Model timestep [s]
   real(r8), intent(in) :: t(pcols,pver)     ! Temperature [K]
   real(r8), intent(in) :: pmid(pcols,pver)  ! Pressure at model levels [Pa]
   real(r8), intent(in) :: q(pcols,pver,ncnst) ! Tracer array including moisture [kg/kg]

   real(r8), intent(in) :: mu(pcols,pver)    ! Updraft mass flux (positive) [mb/s]
   real(r8), intent(in) :: md(pcols,pver)    ! Downdraft mass flux (negative) [mb/s]
   real(r8), intent(in) :: du(pcols,pver)    ! Mass detrain rate from updraft [1/s]
   real(r8), intent(in) :: eu(pcols,pver)    ! Mass entrain rate into updraft [1/s]
   real(r8), intent(in) :: ed(pcols,pver)    ! Mass entrain rate into downdraft [1/s]
! *** note1 - mu, md, eu, ed, du, dp, dpdry are GATHERED ARRAYS ***
! *** note2 - mu and md units are (mb/s), which is used in the zm_conv code
!           - eventually these should be changed to (kg/m2/s)
! *** note3 - eu, ed, du are "d(massflux)/dp" (with dp units = mb), and are all >= 0

   real(r8), intent(in) :: dp(pcols,pver)    ! Delta pressure between interfaces [mb]
   real(r8), intent(in) :: dpdry(pcols,pver) ! Delta dry-pressure [mb]
   integer,  intent(in) :: jt(pcols)         ! Index of cloud top for each column
   integer,  intent(in) :: mx(pcols)         ! Index of cloud bottom for each column
   integer,  intent(in) :: ideep(pcols)      ! Gathering array indices 
   integer,  intent(in) :: il1g              ! Gathered min lon indices over which to operate
   integer,  intent(in) :: il2g              ! Gathered max lon indices over which to operate
! *** note4 -- for il1g <= i <= il2g,  icol = ideep(i) is the "normal" chunk column index

   real(r8), intent(in) :: cldfrac(pcols,pver)  ! Convective cloud fractional area [fraction]
   real(r8), intent(in) :: icwmr(pcols,pver)    ! Convective cloud water from zhang [kg/kg]
   real(r8), intent(in) :: rprd(pcols,pver)     ! Convective precipitation formation rate [kg/kg/s]
   real(r8), intent(in) :: evapc(pcols,pver)    ! Convective precipitation evaporation rate [kg/kg/s]

   real(r8), intent(out):: dqdt(pcols,pver,ncnst)  ! Tracer tendency array [kg/kg/s]
   logical,  intent(in) :: doconvproc(ncnst) ! flag for doing convective transport
   integer,  intent(in) :: nsrflx            ! last dimension of qsrflx
   real(r8), intent(out):: qsrflx(pcols,pcnst,nsrflx)
                              ! process-specific column tracer tendencies [kg/m2/s]
                              !  1 = activation   of interstial to conv-cloudborne
                              !  2 = resuspension of conv-cloudborne to interstital
                              !  3 = aqueous chemistry (not implemented yet, so zero)
                              !  4 = wet removal
                              !  5 = actual precip-evap resuspension (what actually is applied to a species)
                              !  6 = pseudo precip-evap resuspension (for history file) 
   integer,  intent(in) :: species_class(:)  ! specify what kind of species it is. defined at physconst.F90
                                                ! undefined  = 0
                                                ! cldphysics = 1
                                                ! aerosol    = 2
                                                ! gas        = 3
                                                ! other      = 4
   real(r8), intent(out):: xx_mfup_max(pcols)   ! diagnostic field of column maximum updraft mass flux [mb/s]
   real(r8), intent(out):: xx_wcldbase(pcols)
   real(r8), intent(out):: xx_kcldbase(pcols)


!--------------------------Local Variables------------------------------

! cloudborne aerosol, so the arrays are dimensioned with pcnst_extd = pcnst*2
   integer, parameter :: pcnst_extd = pcnst*2

   integer :: lun             ! unit number for diagnostic output
   integer :: ii, icol        ! Work index
   integer :: iconvtype       ! 1=deep, 2=uw shallow
   integer :: iflux_method    ! 1=as in convtran (deep), 2=simpler
   integer :: jtsub           ! Work index
   integer :: kk              ! Work index
   integer :: icnst           ! work index
   integer :: kbot            ! Cloud-flux bottom layer for current i (=mx(i))
   integer :: kbot_prevap     ! Lowest layer for doing resuspension from evaporating precip 
   integer :: ktop            ! Cloud-flux top    layer for current i (=jt(i))
                              ! Layers between kbot,ktop have mass fluxes
                              !    but not all have cloud water, because the
                              !    updraft starts below the cloud base
   integer :: la, lc          ! Work index
   integer :: imode, ispec    ! Work index
   integer :: ntsub           ! 

   logical  doconvproc_extd(pcnst_extd) ! flag for doing convective transport

   real(r8) aqfrac(pcnst_extd)       ! aqueous fraction of constituent in updraft [fraction]
   real(r8) cldfrac_i(pver)          ! cldfrac at current i (with adjustments) [fraction]
   real(r8) chat(pcnst_extd,pverp)   ! mix ratio in env at interfaces  [kg/kg]
   real(r8) cond(pcnst_extd,pverp)   ! mix ratio in downdraft at interfaces [kg/kg]
   real(r8) const(pcnst_extd,pver)   ! gathered tracer array [kg/kg]
   real(r8) conu(pcnst_extd,pverp)   ! mix ratio in updraft at interfaces [kg/kg]

   real(r8) dcondt(pcnst_extd,pver)  ! grid-average TMR tendency for current column [kg/kg/s]
   real(r8) dcondt_prevap(pcnst_extd,pver) ! portion of dcondt from precip evaporation [kg/kg/s]
   real(r8) dcondt_prevap_hist(pcnst_extd,pver) ! similar but used for history output [kg/kg/s]
   real(r8) dcondt_resusp(pcnst_extd,pver) ! portion of dcondt from resuspension [kg/kg/s]
   real(r8) dcondt_wetdep(pcnst_extd,pver) ! portion of dcondt from wet deposition [kg/kg/s]
   real(r8) dconudt_activa(pcnst_extd,pverp) ! d(conu)/dt by activation [kg/kg/s]
   real(r8) dconudt_wetdep(pcnst_extd,pverp) ! d(conu)/dt by wet removal [kg/kg/s]

   real(r8) sumactiva(pcnst_extd)    ! sum (over layers) of dp*dconudt_activa [kg/kg/s]
   real(r8) sumaqchem(pcnst_extd)    ! sum (over layers) of dp*dconudt_aqchem [kg/kg/s]
   real(r8) sumprevap(pcnst_extd)    ! sum (over layers) of dp*dcondt_prevap [kg/kg/s]
   real(r8) sumprevap_hist(pcnst_extd) ! sum (over layers) of dp*dcondt_prevap_hist [kg/kg/s]
   real(r8) sumresusp(pcnst_extd)    ! sum (over layers) of dp*dcondt_resusp [kg/kg/s]
   real(r8) sumwetdep(pcnst_extd)    ! sum (over layers) of dp*dconudt_wetdep [kg/kg/s]

   real(r8) dddp(pver)           ! dd(i,k)*dp(i,k) at current i [mb/s]
   real(r8) dp_i(pver)           ! dp(i,k) at current i [mb]
   real(r8) dpdry_i(pver)        ! dpdry(i,k) at current i [mb]
   real(r8) fa_u_dp              ! fa_u * dp at the current level [mb]
   real(r8) dudp(pver)           ! du(i,k)*dp(i,k) at current i [mb/s]
   real(r8) dz                   ! working layer thickness [m] 
   real(r8) eddp(pver)           ! ed(i,k)*dp(i,k) at current i [mb/s]
   real(r8) eudp(pver)           ! eu(i,k)*dp(i,k) at current i [mb/s]
   real(r8) fa_u(pver)           ! fractional area of in the updraft [fraction]
   real(r8) md_i(pverp)          ! md(i,k) at current i (note pverp dimension) [mb/s]
   real(r8) mu_i(pverp)          ! mu(i,k) at current i (note pverp dimension) [mb/s]
   ! md_i, mu_i, are all "dry" mass fluxes
   real(r8) q_i(pver,pcnst)      ! q(i,k,m) at current i [kg/kg]
   real(r8) rhoair_i(pver)       ! air density at current i [kg/m3]
   real(r8) zmagl(pver)          ! working height above surface [m]


!-----------------------------------------------------------------------
!

   lun = iulog


   if (convtype == 'deep') then
      iconvtype = 1
      iflux_method = 1
   else if (convtype == 'uwsh') then
      iconvtype = 2
      iflux_method = 2
   else
      call endrun( '*** ma_convproc_tend -- convtype is not |deep| or |uwsh|' )
   end if

   ! initiate output variables
   qsrflx(:,:,:) = 0.0_r8
   dqdt(:,:,:) = 0.0_r8
   xx_mfup_max(:) = 0.0_r8
   xx_wcldbase(:) = 0.0_r8
   xx_kcldbase(:) = 0.0_r8

! set doconvproc_extd (extended array) values
! inititialize aqfrac to 1.0 for activated aerosol species, 0.0 otherwise
   call set_cloudborne_vars(                            &
                                ncnst,  doconvproc,     & ! in
                                aqfrac, doconvproc_extd ) ! out

! Loop ever each column that has convection
! *** i is index to gathered arrays; ideep(i) is index to "normal" chunk arrays
i_loop_main_aa: &
   do ii = il1g, il2g
      icol = ideep(ii)

      ! skip some columns
      if ( (jt(ii) <= 0) .and. (mx(ii) <= 0) .and. (iconvtype /= 1) ) then
          ! shallow conv case with jt,mx <= 0, which means there is no shallow conv
          ! in this column -- skip this column
          cycle i_loop_main_aa
      elseif ( (jt(ii) < 1) .or. (mx(ii) > pver) .or. (jt(ii) > mx(ii)) ) then
          ! invalid cloudtop and cloudbase indices -- skip this column
          write(lun,9010) 'illegal jt, mx', convtype, lchnk, icol, ii, jt(ii), mx(ii)
          cycle i_loop_main_aa
      elseif (jt(ii) == mx(ii)) then
          ! cloudtop = cloudbase (1 layer cloud) -- skip this column
          write(lun,9010) 'jt == mx', convtype, lchnk, icol, ii, jt(ii), mx(ii)
          cycle i_loop_main_aa
      endif
9010  format( '*** ma_convproc_tend error -- ', a, 5x, 'convtype = ', a /   &
              '*** lchnk, icol, il, jt, mx = ', 5(1x,i10) )


      ! Load some variables in current column for further subroutine use
      do kk = 1, pver
         dp_i(kk) = dp(ii,kk)
         dpdry_i(kk) = dpdry(ii,kk)
         cldfrac_i(kk) = cldfrac(icol,kk)
         rhoair_i(kk) = pmid(icol,kk)/(rair*t(icol,kk))
      enddo
!  load tracer mixing ratio array, which will be updated at the end of each
!  jtsub interation
      q_i(1:pver,1:pcnst) = q(icol,1:pver,1:pcnst)
      ktop = jt(ii)
      kbot = mx(ii)

      ! calculate dry mass fluxes at cloud layer
      call compute_massflux(                            &
                        ii,     icol,   ktop,   kbot,   & ! in
                        dpdry_i,du,     eu,     ed,     & ! in
                        mu_i,   md_i,                   & ! out
                        xx_mfup_max                     ) ! inout

      ! compute entraintment*dp and detraintment*dp and calculate ntsub
      call compute_ent_det_dp(                                  &
                        ii,     ktop,   kbot,   dt,     dpdry_i,& ! in
                        mu_i,   md_i,   du,     eu,     ed,     & ! in
                        ntsub,  eudp,   dudp,   eddp,   dddp    ) ! out

      ! calculate height of layer interface above ground
      call compute_midlev_height( dpdry_i, rhoair_i, & ! in
                                  zmagl              ) ! out


jtsub_loop_main_aa: &
      do jtsub = 1, ntsub

        ! initialize some tracer mixing ratio arrays
        call initialize_tmr_array(     ncnst,                  & ! in
                       iconvtype,      doconvproc_extd,  q_i,  & ! in
                       const,          chat,    conu,    cond  ) ! out

        ! Compute updraft mixing ratios from cloudbase to cloudtop
        ! ---------------------------------------------------------------------------
        ! there is a bug here for dz, which is not a vertically-varying profile. 
        ! see the subroutine compute_updraft_mixing_ratio for the bug information
        ! by Shuaiqi Tang, 2022
        ! ---------------------------------------------------------------------------
        dz = dpdry_i(1)*hund_ovr_g/rhoair_i(1)

        call compute_updraft_mixing_ratio(                              &
                doconvproc_extd, icol,  ktop,   kbot,   iconvtype,      & ! in
                dt,     dp_i,   dpdry_i,        cldfrac,                & ! in
                rhoair_i,       zmagl,  dz,     mu_i,   eudp,           & ! in
                const,  t,      aqfrac, icwmr,          rprd,           & ! in
                fa_u,   dconudt_wetdep,         dconudt_activa,         & ! out
                conu,           xx_wcldbase,    xx_kcldbase             ) ! inout


        ! Compute downdraft mixing ratios from cloudtop to cloudbase
        call compute_downdraft_mixing_ratio(                    &
                        doconvproc_extd,       ktop, kbot,      & ! in
                        md_i,           eddp,        const,     & ! in
                        cond                                    ) ! inout

        ! Now computes fluxes and tendencies
        ! NOTE:  The approach used in convtran applies to inert tracers and
        !        must be modified to include source and sink terms
        call initialize_dcondt(                                 &
                doconvproc_extd, iflux_method, ktop, kbot,      & ! in
                dpdry_i, fa_u,      mu_i,       md_i,           & ! in
                chat,    const,     conu,       cond,           & ! in
                dconudt_activa,     dconudt_wetdep,             & ! in
                dudp,    dddp,      eudp,       eddp,           & ! in
                dcondt                                          ) ! out

        ! compute dcondt_wetdep for next subroutine
        dcondt_wetdep(:,:) = 0.0
        do kk = ktop, kbot
           fa_u_dp = fa_u(kk)*dpdry_i(kk) ! simply cancelling dpdry_i causes BFB test fail
           do icnst = 2, pcnst_extd
              if (doconvproc_extd(icnst)) then
                 dcondt_wetdep(icnst,kk) = fa_u_dp*dconudt_wetdep(icnst,kk)/dpdry_i(kk)
              endif
           enddo
        enddo

        ! calculate effects of precipitation evaporation
        call ma_precpevap_convproc(                             & 
               dcondt_wetdep,                                   & ! in
               rprd,   evapc,           dpdry_i,                & ! in
               icol,   ktop,            pcnst_extd,             & ! in
               doconvproc_extd,         species_class,          & ! in
               dcondt_prevap,          dcondt_prevap_hist,      & ! out
               dcondt                                           & ! inout
               ) 


! make adjustments to dcondt for activated & unactivated aerosol species
!    pairs to account any (or total) resuspension of convective-cloudborne aerosol

! usually the updraft ( & downdraft) start ( & end ) at kbot=pver, but sometimes kbot < pver
! transport, activation, resuspension, and wet removal only occur between kbot >= k >= ktop
! resuspension from evaporating precip can occur at k > kbot when kbot < pver
! in the first version of this routine, the precp evap resusp tendencies for k > kbot were ignored,
!    but that is now fixed
! this was a minor bug with quite minor affects on the aerosol,
!    because convective precip evap is (or used to be) much less than stratiform precip evap )
!      kbot_prevap = kbot
! apply this minor fix when doing resuspend to coarse mode
        kbot_prevap = pver

        call ma_resuspend_convproc( dcondt,                             & ! inout
                                    ktop, kbot_prevap, pcnst_extd,      & ! in
                                    dcondt_resusp                       ) ! out

        ! calculate new column-tendency variables
        call compute_tendency_resusp_evap(                              &
                doconvproc_extd, ktop,          kbot_prevap,   dpdry_i, & ! in
                dcondt_resusp,  dcondt_prevap,  dcondt_prevap_hist,     & ! in
                dconudt_activa, dconudt_wetdep, fa_u,                   & ! in
                sumactiva,      sumaqchem,      sumwetdep,              & ! out
                sumresusp,      sumprevap,      sumprevap_hist          ) ! out

        ! update tendencies
        call update_tendency_final(                             &
                        ktop,   kbot_prevap,    ntsub,  jtsub,  & ! in
                        ncnst,  nsrflx,         icol,           & ! in
                        dt,     dcondt,         doconvproc,     & ! in
                        sumactiva,              sumaqchem,      & ! inout
                        sumwetdep,              sumresusp,      & ! inout
                        sumprevap,              sumprevap_hist, & ! inout
                        dqdt,   q_i,            qsrflx          ) ! inout

      enddo jtsub_loop_main_aa  ! of the main "do jtsub = 1, ntsub" loop
   enddo i_loop_main_aa  ! of the main "do i = il1g, il2g" loop

   return
   end subroutine ma_convproc_tend

!====================================================================================
   subroutine set_cloudborne_vars(                      &
                                ncnst,  doconvproc,     & ! in
                                aqfrac, doconvproc_extd ) ! out
!-----------------------------------------------------------------------
! set cloudborne aerosol related variables:
! doconvproc_extd: extended array for both activated and unactivated aerosols
! aqfrac: set as 1.0 for activated aerosols and 0.0 otherwise
!-----------------------------------------------------------------------

   ! cloudborne aerosol, so the arrays are dimensioned with pcnst_extd = pcnst*2
   integer,  parameter  :: pcnst_extd = pcnst*2
   integer,  intent(in) :: ncnst                ! number of tracers to transport
   logical,  intent(in) :: doconvproc(ncnst)    ! flag for doing convective transport


   logical,  intent(out) :: doconvproc_extd(pcnst_extd)    ! flag for doing convective transport
   real(r8), intent(out) :: aqfrac(pcnst_extd)  ! aqueous fraction of constituent in updraft [fraction]

   ! local variables
   integer  :: la, lc, imode, ispec   ! indices


   doconvproc_extd(:) = .false.
   doconvproc_extd(2:ncnst) = doconvproc(2:ncnst)
   aqfrac(:) = 0.0_r8
   do imode = 1, ntot_amode
      do ispec = 0, nspec_amode(imode)
         call assign_la_lc(imode, ispec, la, lc)
         if ( doconvproc(la) ) then
            doconvproc_extd(lc) = .true.
            aqfrac(lc) = 1.0_r8
         end if
      enddo
   enddo 

   end subroutine set_cloudborne_vars

!====================================================================================
   subroutine compute_massflux(                         &
                        ii,     icol,   ktop,   kbot,   & ! in
                        dpdry_i,du,     eu,     ed,     & ! in
                        mu_i,   md_i,                   & ! out
                        xx_mfup_max                     ) ! inout
!-----------------------------------------------------------------------
! compute dry mass fluxes
! This is approximate because the updraft air is has different temp and qv than
! the grid mean, but the whole convective parameterization is highly approximate
! values below a threshold and at "top of cloudtop", "base of cloudbase" are set as zero
!-----------------------------------------------------------------------

   integer,  intent(in)  :: ii                   ! index to gathered arrays
   integer,  intent(in)  :: icol                 ! column index
   integer,  intent(in)  :: ktop                 ! top level index
   integer,  intent(in)  :: kbot                 ! bottom level index
   real(r8), intent(in)  :: dpdry_i(pver)        ! dp [mb]
   real(r8), intent(in)  :: du(pcols,pver)       ! Mass detrain rate from updraft [1/s]
   real(r8), intent(in)  :: eu(pcols,pver)       ! Mass entrain rate into updraft [1/s]
   real(r8), intent(in)  :: ed(pcols,pver)       ! Mass entrain rate into downdraft [1/s]
   real(r8), intent(out) :: mu_i(pverp)          ! mu at current i (note pverp dimension, see ma_convproc_tend) [mb/s]
   real(r8), intent(out) :: md_i(pverp)          ! md at current i (note pverp dimension) [mb/s]
   real(r8), intent(inout) :: xx_mfup_max(pcols) ! diagnostic field of column maximum updraft mass flux [mb/s]


   ! local variables
   integer  :: kk            ! index
   real(r8) :: mu_x(pverp)   ! updraft mass flux calculated for all interface layers [mb/s]
   real(r8) :: md_x(pverp)   ! downdraft mass flux calculated for all interface layers [mb/s]

! first calculate updraft and downdraft mass fluxes for all layers
   mu_x(:) = 0.0
   md_x(:) = 0.0
   ! (eu-du) = d(mu)/dp -- integrate upwards, multiplying by dpdry
   do kk = pver, 1, -1
         mu_x(kk) = mu_x(kk+1) + (eu(ii,kk)-du(ii,kk))*dpdry_i(kk)
         xx_mfup_max(icol) = max( xx_mfup_max(icol), mu_x(kk) )
   enddo
   ! (ed) = d(md)/dp -- integrate downwards, multiplying by dpdry
   do kk = 2, pver
         md_x(kk) = md_x(kk-1) - ed(ii,kk-1)*dpdry_i(kk-1)
   enddo

! then load mass fluxes only over cloud layers
! excluding "top of cloudtop", "base of cloudbase"
   mu_i(:) = 0.0
   md_i(:) = 0.0
   do kk = ktop+1, kbot
      ! load mass fluxes at cloud layers
      mu_i(kk) = mu_x(kk)
      md_i(kk) = md_x(kk)
      ! zero out values below threshold
      if (mu_i(kk) <= mbsth) mu_i(kk) = 0.0
      if (md_i(kk) >= -mbsth) md_i(kk) = 0.0
   enddo

   end subroutine compute_massflux

!====================================================================================
   subroutine compute_ent_det_dp(                               &
                        ii,     ktop,   kbot,   dt,     dpdry_i,& ! in
                        mu_i,   md_i,   du,     eu,     ed,     & ! in
                        ntsub,  eudp,   dudp,   eddp,   dddp    ) ! out
!-----------------------------------------------------------------------
! calculate mass flux change (from entrainment or detrainment) in the current dp
! also get number of time substeps
!-----------------------------------------------------------------------

   implicit none

   integer,  intent(in)  :: ii                   ! index to gathered arrays
   integer,  intent(in)  :: ktop                 ! top level index
   integer,  intent(in)  :: kbot                 ! bottom level index
   real(r8), intent(in)  :: dt                   ! delta t (model time increment) [s]
   real(r8), intent(in)  :: dpdry_i(pver)        ! dp [mb]
   real(r8), intent(in)  :: mu_i(pverp)          ! mu at current i (note pverp dimension, see ma_convproc_tend) [mb/s]
   real(r8), intent(in)  :: md_i(pverp)          ! md at current i (note pverp dimension) [mb/s]
   real(r8), intent(in)  :: du(pcols,pver)       ! Mass detrain rate from updraft [1/s]
   real(r8), intent(in)  :: eu(pcols,pver)       ! Mass entrain rate into updraft [1/s]
   real(r8), intent(in)  :: ed(pcols,pver)       ! Mass entrain rate into downdraft [1/s]

   integer,  intent(out) :: ntsub               ! number of sub timesteps
   real(r8), intent(out) :: eudp(pver)          ! eu(i,k)*dp(i,k) at current i [mb/s]
   real(r8), intent(out) :: dudp(pver)          ! du(i,k)*dp(i,k) at current i [mb/s]
   real(r8), intent(out) :: eddp(pver)          ! ed(i,k)*dp(i,k) at current i [mb/s]
   real(r8), intent(out) :: dddp(pver)          ! dd(i,k)*dp(i,k) at current i [mb/s]

   ! local variables
   integer  :: kk            ! index
   real(r8) :: courantmax    ! maximum value of courant number [unitless]


   ! initiate variables
   eudp(:) = 0.0
   dudp(:) = 0.0
   eddp(:) = 0.0
   dddp(:) = 0.0
   courantmax = 0.0
   ntsub = 1

   !  Compute updraft and downdraft "entrainment*dp" from eu and ed
   !  Compute "detrainment*dp" from mass conservation (total is mass flux
   !  difference between the top an bottom interface of this layer)
   do kk = ktop, kbot
       if ((mu_i(kk) > 0) .or. (mu_i(kk+1) > 0)) then
            if (du(ii,kk) <= 0.0) then
               eudp(kk) = mu_i(kk) - mu_i(kk+1)
            else
               eudp(kk) = max( eu(ii,kk)*dpdry_i(kk), 0.0_r8 )
               dudp(kk) = (mu_i(kk+1) + eudp(kk)) - mu_i(kk)
               if (dudp(kk) < 1.0e-12*eudp(kk)) then
                  eudp(kk) = mu_i(kk) - mu_i(kk+1)
                  dudp(kk) = 0.0
               endif
            endif
       endif
       if ((md_i(kk) < 0) .or. (md_i(kk+1) < 0)) then
            eddp(kk) = max( ed(ii,kk)*dpdry_i(kk), 0.0_r8 )
            dddp(kk) = (md_i(kk+1) + eddp(kk)) - md_i(kk)
            if (dddp(kk) < 1.0e-12*eddp(kk)) then
               eddp(kk) = md_i(kk) - md_i(kk+1)
               dddp(kk) = 0.0
            endif
       endif
       ! get courantmax to calculate ntsub
       courantmax = max( courantmax, ( mu_i(kk+1)+eudp(kk)-md_i(kk)+eddp(kk))*dt/dpdry_i(kk) )
   enddo ! kk

   ! number of time substeps needed to maintain "courant number" <= 1
   if (courantmax > (1.0_r8 + 1.0e-6_r8)) then
       ntsub = 1 + int( courantmax )
   endif

   end subroutine compute_ent_det_dp

!====================================================================================
   subroutine compute_midlev_height( dpdry_i, rhoair_i, & ! in
                                     zmagl              ) ! out
!-----------------------------------------------------------------------
! compute height above surface for middle of level kk
!-----------------------------------------------------------------------

   real(r8), intent(in)  :: dpdry_i(pver)       ! dp [mb]
   real(r8), intent(in)  :: rhoair_i(pver)      ! air density [kg/m3]
   real(r8), intent(out) :: zmagl(pver)         ! height above surface at middle level [m]   
   
   ! local variables
   integer  :: kk            ! index
   real(r8) :: dz            ! layer thickness [m]

   ! at surface (level pver)
   dz = dpdry_i(pver)*hund_ovr_g/rhoair_i(pver)
   zmagl(pver) = 0.5*dz
   ! other levels
   do kk = pver-1, 1, -1
       zmagl(kk) = zmagl(kk+1) + 0.5*dz         ! add half layer below
       dz = dpdry_i(kk)*hund_ovr_g/rhoair_i(kk) ! update layer thickness at level kk
       zmagl(kk) = zmagl(kk) + 0.5*dz           ! add half layer in this level
   enddo

   end subroutine compute_midlev_height

!====================================================================================
   subroutine initialize_tmr_array(     ncnst,                  &
                        iconvtype,      doconvproc_extd,  q_i,  & ! in
                        const,          chat,      conu,  cond  ) ! out
!-----------------------------------------------------------------------
! initialize tracer mixing ratio arrays (const, chat, conu, cond)
! chat, conu and cond are at interfaces; interpolation needed
! Note: for deep convection, some values between the two layers 
! differ significantly, use geometric averaging under certain conditions
!-----------------------------------------------------------------------

   integer,  parameter  :: pcnst_extd = pcnst*2

   integer, intent(in)  :: ncnst                     ! number of tracers to transport
   integer, intent(in)  :: iconvtype                 ! 1=deep, 2=uw shallow
   logical, intent(in)  :: doconvproc_extd(pcnst_extd)    ! flag for doing convective transport
   real(r8), intent(in) :: q_i(pver,pcnst)           ! q(icol,kk,icnst) at current icol

   real(r8), intent(out) :: const(pcnst_extd,pver)   ! gathered tracer array [kg/kg]
   real(r8), intent(out) :: chat(pcnst_extd,pverp)   ! mix ratio in env at interfaces [kg/kg]
   real(r8), intent(out) :: conu(pcnst_extd,pverp)   ! mix ratio in updraft at interfaces [kg/kg]
   real(r8), intent(out) :: cond(pcnst_extd,pverp)   ! mix ratio in downdraft at interfaces [kg/kg]

   ! local variables
   integer  :: icnst            ! index
   integer  :: kk               ! vertical index
   integer  :: km1              ! index of kk-1
   real(r8) :: max_con          ! max const concentration at level kk and kk-1 [kg/kg]
   real(r8) :: min_con          ! min const concentration at level kk and kk-1 [kg/kg]
   real(r8) :: c_dif_rel        ! relative difference between level kk and kk-1 [unitless]
   real(r8) :: c_above          ! const at the above (kk-1 level) [kg/kg]
   real(r8) :: c_below          ! const at the below (kk level) [kg/kg]

   real(r8), parameter :: small_con = 1.e-36        ! threshold of constitute as zero [kg/kg]
   real(r8), parameter :: small_rel = 1.0e-6_r8         ! small value for relative comparison


   ! initiate variables
    const(:,:) = 0.0
    chat(:,:) = 0.0 
    conu(:,:) = 0.0
    cond(:,:) = 0.0   

    do icnst = 2, ncnst
        if ( doconvproc_extd(icnst) ) then
            ! Gather up the constituent
            do kk = 1,pver
                const(icnst,kk) = q_i(kk,icnst)
            enddo

           ! Interpolate environment tracer values to interfaces
           do kk = 1,pver
                km1 = max(1,kk-1)

                ! get relative difference between the two levels
                min_con = min(const(icnst,km1),const(icnst,kk))
                max_con = max(const(icnst,km1),const(icnst,kk))
                if (min_con < 0) then
                   c_dif_rel = 0.
                else
                   c_dif_rel = abs(const(icnst,kk)-const(icnst,km1))/max(max_con,small_con)
                endif

                ! If the two layers differ significantly use a geometric averaging procedure
                ! But only do that for deep convection.  For shallow, use the simple
                ! averaging which is used in subr cmfmca
                if (iconvtype /= 1) then ! simple averaging for non-deep convection
                    chat(icnst,kk) = 0.5* (const(icnst,kk)+const(icnst,km1))
                elseif (c_dif_rel > small_rel) then ! deep convection using geometric averaging
                    c_above = max(const(icnst,km1),max_con*1.e-12_r8)
                    c_below = max(const(icnst,kk),max_con*1.e-12_r8)
                    chat(icnst,kk) = log(c_above/c_below)/(c_above-c_below)*c_above*c_below
                else             ! Small diff, so just arithmetic mean
                    chat(icnst,kk) = 0.5* (const(icnst,kk)+const(icnst,km1))
                endif

                ! Set provisional up and down draft values, and tendencies
                conu(icnst,kk) = chat(icnst,kk)
                cond(icnst,kk) = chat(icnst,kk)
            enddo ! kk

            ! Values at surface inferface == values in lowest layer
            chat(icnst,pver+1) = const(icnst,pver)
            conu(icnst,pver+1) = const(icnst,pver)
            cond(icnst,pver+1) = const(icnst,pver)

        endif ! doconvproc_extd(icnst)
   enddo ! icnst

   end subroutine initialize_tmr_array

!====================================================================================
   subroutine compute_updraft_mixing_ratio(                             &
                doconvproc_extd, icol,  ktop,   kbot,   iconvtype,      & ! in
                dt,     dp_i,   dpdry_i,        cldfrac,                & ! in
                rhoair_i,       zmagl,  dz,     mu_i,   eudp,           & ! in
                const,  t,      aqfrac, icwmr,          rprd,           & ! in
                fa_u,   dconudt_wetdep,         dconudt_activa,         & ! out
                conu,           xx_wcldbase,    xx_kcldbase             ) ! inout 
!-----------------------------------------------------------------------
! Compute updraft mixing ratios from cloudbase to cloudtop
! No special treatment is needed at k=pver because arrays are dimensioned 1:pver+1
! A time-split approach is used.  First, entrainment is applied to produce
!    an initial conu(m,k) from conu(m,k+1).  Next, chemistry/physics are
!    applied to the initial conu(m,k) to produce a final conu(m,k).
!    Detrainment from the updraft uses this final conu(m,k).
! Note that different time-split approaches would give somewhat different results
!-----------------------------------------------------------------------

   ! cloudborne aerosol, so the arrays are dimensioned with pcnst_extd = pcnst*2
   integer,  parameter  :: pcnst_extd = pcnst*2
   logical, intent(in)  :: doconvproc_extd(pcnst_extd) ! flag for doing convective transport
   integer, intent(in)  :: icol                 ! column index
   integer, intent(in)  :: ktop                 ! top level index
   integer, intent(in)  :: kbot                 ! bottom level index
   integer, intent(in)  :: iconvtype            ! 1=deep, 2=uw shallow
   real(r8),intent(in)  :: dt                   ! Model timestep [s]
   real(r8),intent(in)  :: dp_i(pver)           ! dp [mb]
   real(r8),intent(in)  :: dpdry_i(pver)        ! dp [mb]
   real(r8),intent(in)  :: cldfrac(pcols, pver) ! cloud fraction [fraction]
   real(r8),intent(in)  :: rhoair_i(pver)       ! air density at current i [kg/m3]
   real(r8),intent(in)  :: zmagl(pver)          ! height above surface [m]
   real(r8),intent(in)  :: dz                   ! layer thickness [m]
   real(r8),intent(in)  :: mu_i(pverp)          ! mu at current i (note pverp dimension, see ma_convproc_tend) [mb/s]
   real(r8),intent(in)  :: eudp(pver)           ! eu(i,k)*dp(i,k) at current i [mb/s]
   real(r8),intent(in)  :: const(pcnst_extd,pver)   ! gathered tracer array [kg/kg]

   real(r8), intent(in) :: t(pcols,pver)        ! Temperature [K]
   real(r8), intent(in) :: aqfrac(pcnst_extd)   ! aqueous fraction of constituent in updraft [fraction]
   real(r8), intent(in) :: icwmr(pcols,pver)    ! Convective cloud water from zm scheme [kg/kg]
   real(r8), intent(in) :: rprd(pcols,pver)     ! Convective precipitation formation rate [kg/kg/s]
   real(r8), intent(out) :: dconudt_wetdep(pcnst_extd,pverp) ! d(conu)/dt by wet removal[kg/kg/s]
   real(r8), intent(out) :: dconudt_activa(pcnst_extd,pverp) ! d(conu)/dt by activation [kg/kg/s]
   real(r8), intent(out) :: fa_u(pver)           ! fractional area of in the updraft
   real(r8), intent(inout) :: conu(pcnst_extd,pverp)   ! mix ratio in updraft at interfaces [kg/kg]
   real(r8), intent(inout) :: xx_wcldbase(pcols) ! w at first cloudy layer [m/s]
   real(r8), intent(inout) :: xx_kcldbase(pcols) ! level of cloud base

   ! local variables
   integer      :: icnst        ! index of pcnst_extd
   integer      :: kk           ! index of vertical level
   integer      :: kp1          ! index of kk+1
   integer      :: kactcnt      ! Counter for no. of levels having activation
   integer      :: kactfirst    ! Lowest layer with activation (= cloudbase)
   real(r8)     :: dt_u(pver)           ! lagrangian transport time in the updraft [s]
   real(r8)     :: mu_p_eudp(pver)      ! = mu_i(kp1) + eudp(k) [mb/s]
   real(r8)     :: cldfrac_i(pver)      ! cldfrac at current icol (with adjustments) [fraction]
   real(r8)     :: wup(pver)            ! mean updraft vertical velocity at current level updraft [m/s]
   real(r8)     :: f_ent                ! fraction of updraft massflux that was entrained across this layer == eudp/mu_p_eudp [fraction]

   real(r8), parameter :: cldfrac_cut = 0.005_r8        ! cutoff value of cloud fraction to remove zero cloud [fraction]

    ! initiate variables
    kactcnt = 0 ; kactfirst = 1
    wup(:) = 0.0
    dconudt_wetdep(:,:) = 0.0
    dconudt_activa(:,:) = 0.0

    ! cldfrac = conv cloud fractional area.  This could represent anvil cirrus area,
    ! and may not useful for aqueous chem and wet removal calculations
    ! adjust to remove zero clouds
    do kk = kbot, ktop, -1
        cldfrac_i(kk) = cldfrac(icol,kk)
        cldfrac_i(kk) = max( cldfrac_i(kk), cldfrac_cut )
    enddo

    do kk = kbot, ktop, -1
        kp1 = kk+1

        !Initialized so that it has a value if the following "if" check yeilds
        !.false.
        fa_u(kk) = 0.0_r8

        ! mu_p_eudp(k) = updraft massflux at k, without detrainment between kp1,k
        mu_p_eudp(kk) = mu_i(kp1) + eudp(kk)
        if (mu_p_eudp(kk) > mbsth) then
        ! if (mu_p_eudp(k) <= mbsth) the updraft mass flux is negligible 
        ! at base and top of current layer,
        ! so current layer is a "gap" between two unconnected updrafts,
        ! so essentially skip all the updraft calculations for this layer

            ! First apply changes from entrainment
            f_ent = eudp(kk)/mu_p_eudp(kk)
            f_ent = min_max_bound(0.0_r8, 1.0_r8, f_ent)
            do icnst = 2, pcnst_extd
               if (doconvproc_extd(icnst)) then
                  conu(icnst,kk) = (1.0_r8 - f_ent)*conu(icnst,kp1) + f_ent*const(icnst,kk)
               endif
            enddo

            ! estimate updraft velocity (wup)
            call compute_wup(                           &
                        iconvtype,      kk,    mu_i,    & ! in
                        cldfrac_i,   rhoair_i, zmagl,   & ! in
                        wup                             ) ! inout
            ! compute lagrangian transport time (dt_u)
            ! -------------------------------------------------------------------
            ! There is a bug here that dz is not vertically varying but fixed as
            ! the thickness of the lowest level calculated previously in the
            ! subroutine ma_convproc_tend. keep it here for C++ refacoring but
            ! may need to fix later
            ! found by Shuaiqi Tang, 2022
            ! -------------------------------------------------------------------
            dt_u(kk) = dz/wup(kk)
            dt_u(kk) = min( dt_u(kk), dt )

! Now apply transformation and removal changes

            ! aerosol activation - method 2
            call compute_activation_tend(                       &
                        icol,           kk,             f_ent,  & ! in

                        cldfrac_i,      rhoair_i,       mu_i,   & ! in
                        dt_u,   wup,    icwmr,          t,      & ! in
                        kactcnt,        kactfirst,              & ! inout
                        conu,           dconudt_activa,         & ! inout
                        xx_wcldbase,    xx_kcldbase             ) ! inout

            ! wet removal
            call compute_wetdep_tend(                              &
                doconvproc_extd,        icol,   kk,     dt,        & ! in
                dt_u,   dp_i,           cldfrac_i,      mu_p_eudp, & ! in
                aqfrac,                 icwmr,          rprd,      & ! in
                conu,                   dconudt_wetdep             ) ! inout

            ! compute updraft fractional area; for update fluxes use
            ! *** these must obey  dt_u(k)*mu_p_eudp(k) = dpdry_i(k)*fa_u(k)
            fa_u(kk) = dt_u(kk)*(mu_p_eudp(kk)/dpdry_i(kk))

       endif    ! "(mu_p_eudp(k) > mbsth)"
   enddo  ! "kk = kbot, ktop, -1"

   end subroutine compute_updraft_mixing_ratio

!====================================================================================
   subroutine compute_wup(                              &
                        iconvtype,      kk,     mu_i,   & ! in
                        cldfrac_i,   rhoair_i, zmagl,   & ! in
                        wup                             ) ! inout
!-----------------------------------------------------------------------
! estimate updraft velocity (wup)
! do it differently for deep and shallow convection
!-----------------------------------------------------------------------

   integer,  intent(in) :: iconvtype            ! 1=deep, 2=uw shallow
   integer,  intent(in) :: kk                   ! vertical level index
   real(r8), intent(in) :: mu_i(pverp)          ! mu at current i (note pverp dimension) [mb/s]
   real(r8), intent(in) :: cldfrac_i(pver)      ! cldfrac at current icol (with adjustments) [fraction]
   real(r8), intent(in) :: rhoair_i(pver)       ! air density at current i [kg/m3]
   real(r8), intent(in) :: zmagl(pver)          ! height above surface [m]
   real(r8), intent(inout)      :: wup(pver)    ! mean updraft vertical velocity at current level updraft [m/s]

   integer      :: kp1                  ! kk + 1
   real(r8)     :: zkm                  ! height above surface [km]
   real(r8), parameter  :: w_peak = 4.0_r8      ! pre-defined peak updraft [m/s]

! shallow - wup = (mup in kg/m2/s) / [rhoair * (updraft area)]
    if (iconvtype /= 1) then
       wup(kk) = (mu_i(kp1) + mu_i(kk))*0.5_r8*hund_ovr_g &
                      / (rhoair_i(kk) * (cldfrac_i(kk)*0.5_r8))
       wup(kk) = max( 0.1_r8, wup(kk) )

! deep - the above method overestimates updraft area and underestimate wup
!    the following is based lemone and zipser (j atmos sci, 1980, p. 2455)
!    peak updraft (= 4 m/s) is sort of a "grand median" from their GATE data
!       and Thunderstorm Project data which they also show
!    the vertical profile shape is a crude fit to their median updraft profile
   else
       zkm = zmagl(kk)*1.0e-3
       if (zkm .ge. 1.0) then
           wup(kk) = w_peak*((zkm/w_peak)**0.21_r8)
       else
           wup(kk) = 2.9897_r8*(zkm**0.5_r8)
       endif
       wup(kk) = min_max_bound(0.1_r8, w_peak, wup(kk))
   endif

   end subroutine compute_wup

!====================================================================================
   subroutine compute_activation_tend(                          &
                        icol,           kk,             f_ent,  & ! in
                        cldfrac_i,      rhoair_i,       mu_i,   & ! in
                        dt_u,   wup,    icwmr,          t,      & ! in
                        kactcnt,        kactfirst,              & ! inout
                        conu,           dconudt_activa,         & ! inout
                        xx_wcldbase,    xx_kcldbase             ) ! inout
!-----------------------------------------------------------------------
! aerosol activation - method 2
!    when kactcnt=1 (first/lowest layer with cloud water)
!       apply "primary" activatation to the entire updraft
!    when kactcnt>1
!       apply secondary activatation to the entire updraft
!       do this for all levels above cloud base (even if completely glaciated)
!          (this is something for sensitivity testing)
!-----------------------------------------------------------------------

   ! cloudborne aerosol, so the arrays are dimensioned with pcnst_extd = pcnst*2
   integer,  parameter  :: pcnst_extd = pcnst*2
   integer,  intent(in) :: icol                 ! column index
   integer,  intent(in) :: kk                   ! vertical level index
   real(r8), intent(in) :: f_ent                ! fraction of updraft massflux that was entrained across this layer == eudp/mu_p_eudp [fraction]
   real(r8), intent(in) :: cldfrac_i(pver)      ! cldfrac at current icol (with adjustments) [fraction]
   real(r8), intent(in) :: rhoair_i(pver)       ! air density at current i [kg/m3]
   real(r8), intent(in) :: mu_i(pverp)          ! mu at current i (note pverp dimension) [mb/s]
   real(r8), intent(in) :: dt_u(pver)           ! lagrangian transport time in the updraft at current level [s]
   real(r8), intent(in) :: wup(pver)            ! mean updraft vertical velocity at current level updraft [m/s]
   real(r8), intent(in) :: t(pcols,pver)        ! Temperature [K]
   real(r8), intent(in) :: icwmr(pcols,pver)    ! Convective cloud water from zm scheme [kg/kg]
   integer,  intent(inout) :: kactcnt           !  Counter for no. of levels having activation
   integer,  intent(inout) :: kactfirst         ! Lowest layer with activation (= cloudbase)
   real(r8), intent(inout) :: conu(pcnst_extd,pverp)   ! mix ratio in updraft at interfaces [kg/kg]
   real(r8), intent(inout) :: dconudt_activa(pcnst_extd,pverp) ! d(conu)/dt by activation [kg/kg/s]
   real(r8), intent(inout) :: xx_wcldbase(pcols) ! w at first cloudy layer [m/s]
   real(r8), intent(inout) :: xx_kcldbase(pcols) ! level of cloud base
                                ! xx_kcldbase is filled with index kk. maybe better define as integer

   logical      :: do_act_this_lev      ! flag for doing activation at current level
   integer      :: kp1                  ! kk + 1
   real(r8), parameter :: clw_cut = 1.0e-6      ! cutoff value of cloud water for doing updraft [kg/kg]
                ! Skip levels where icwmr(icol,k) <= clw_cut (=1.0e-6) to
                ! eliminate occasional very small icwmr values from the ZM module

   ! initiate variables
   do_act_this_lev = .false.
   kp1 = kk+1

   if (kactcnt <= 0) then
       if (icwmr(icol,kk) > clw_cut) then
            do_act_this_lev = .true.
            kactcnt = 1
            kactfirst = kk
            ! diagnostic fields
            ! xx_wcldbase = w at first cloudy layer, estimated from mu and cldfrac
            xx_wcldbase(icol) = (mu_i(kp1) + mu_i(kk))*0.5_r8*hund_ovr_g &
                         / (rhoair_i(kk) * (cldfrac_i(kk)*0.5_r8))
            xx_kcldbase(icol) = kk
        endif
   else
        do_act_this_lev = .true.
        kactcnt = kactcnt + 1
   endif

   if ( do_act_this_lev ) then
        call ma_activate_convproc(                                 &
                     conu(:,kk), dconudt_activa(:,kk),             & ! inout
                     f_ent,      dt_u(kk),            wup(kk),     & ! in
                     t(icol,kk), rhoair_i(kk),                     & ! in
                     pcnst_extd, kk,                  kactfirst    ) ! in
   endif

   end subroutine compute_activation_tend

!=========================================================================================
   subroutine ma_activate_convproc(             &
              conu,       dconudt,              & ! inout
              f_ent,      dt_u,      wup,       & ! in
              tair,       rhoair,               & ! in
              pcnst_extd, kk,        kactfirst  ) ! in
!-----------------------------------------------------------------------
!
! Purpose:
! Calculate activation of aerosol species in convective updraft
! for a single column and level
!
! Method:
! conu(l)    = Updraft TMR (tracer mixing ratio) at k/k-1 interface
! f_ent      = Fraction of the "before-detrainment" updraft massflux at
!              k/k-1 interface" resulting from entrainment of level k air
!              (where k is the current level in subr ma_convproc_tend)
!
! On entry to this routine, the conu(l) represents the updraft TMR
! after entrainment, but before chemistry/physics and detrainment.
!
! This routine applies aerosol activation to the conu tracer mixing ratios,
! then adjusts the conu so that on exit,
!   conu(la) = conu_incoming(la) - conu(la)*f_act(la)
!   conu(lc) = conu_incoming(lc) + conu(la)*f_act(la)
! where
!   la, lc   = indices for an unactivated/activated aerosol component pair
!   f_act    = fraction of conu(la) that is activated.  The f_act are
!              calculated with the Razzak-Ghan activation parameterization.
!              The f_act differ for each mode, and for number/surface/mass.
!
! At cloud base (k==kactfirst), primary activation is done using the
! "standard" code in subr activate do diagnose maximum supersaturation.
! Above cloud base, secondary activation is done using a
! prescribed supersaturation.
!
! *** The updraft velocity used for activation calculations is rather
!     uncertain and needs more work.  However, an updraft of 1-3 m/s
!     will activate essentially all of accumulation and coarse mode particles.
!
! Author: R. Easter
!
!-----------------------------------------------------------------------

   implicit none

!-----------------------------------------------------------------------
! arguments  (note:  TMR = tracer mixing ratio)
   integer, intent(in)     :: pcnst_extd        ! indices
   ! conu = tracer mixing ratios in updraft at top of this (current) level. The
   ! conu are changed by activation
   real(r8), intent(inout) :: conu(pcnst_extd)    ! [#/kg or kg/kg]
   real(r8), intent(inout) :: dconudt(pcnst_extd) ! TMR tendencies due to activation [#/kg/s or kg/kg/s]
   real(r8), intent(in)    :: f_ent     ! fraction of updraft massflux that was entrained across this layer == eudp/mu_p_eudp [fraction]
   real(r8), intent(in)    :: dt_u      ! lagrangian transport time in the updraft at current level [s]
   real(r8), intent(in)    :: wup       ! mean updraft vertical velocity at current level updraft [m/s]
   real(r8), intent(in)    :: tair      ! Temperature [K]
   real(r8), intent(in)    :: rhoair    ! air density [kg/m3]
   integer,  intent(in)    :: kk        ! level index
   integer,  intent(in)    :: kactfirst ! k at cloud base

!-----------------------------------------------------------------------
! local variables
   integer  :: la, lc, imode, ispec   ! indices

   real(r8) :: dt_u_inv               ! 1.0/dt_u [1/s]
   real(r8) :: fluxm(ntot_amode)      ! output of activate_modal not used in this subroutine. to understand this, see subr activate_modal
   real(r8) :: fluxn(ntot_amode)      ! output of activate_modal not used in this subroutine. to understand this, see subr activate_modal
   real(r8) :: flux_fullact           ! output of activate_modal not used in this subroutine. to understand this, see subr activate_modal
   real(r8) :: fm(ntot_amode)         ! mass fraction of aerosols activated [fraction]
   real(r8) :: fn(ntot_amode)         ! number fraction of aerosols activated [fraction]
   real(r8) :: hygro(ntot_amode)      ! current hygroscopicity for int+act [unitless]
   real(r8) :: naerosol(ntot_amode)   ! interstitial+activated number conc [#/m3]
   real(r8) :: sigw                   ! standard deviation of updraft velocity [m/s]
   real(r8) :: act_frac               ! activation fraction [fraction]
   real(r8) :: vaerosol(ntot_amode)   ! int+act volume [m3/m3]
   real(r8) :: wbar                   ! mean updraft velocity [m/s]
   real(r8) :: wdiab                  ! diabatic vertical velocity [m/s]
   real(r8) :: wminf, wmaxf           ! limits for integration over updraft spectrum [m/s]

!  activate_smaxmax = the uniform or peak supersat value (as 0-1 fraction =
!  percent*0.01)
   real(r8), parameter :: activate_smaxmax = 0.003_r8

!-----------------------------------------------------------------------

! check f_ent > 0
   if (f_ent <= 0.0_r8) return

! calculate aerosol (or a+cw) volume, number and hygroscopicity
   call aer_vol_num_hygro( pcnst_extd, conu,     rhoair,    & ! in
                           vaerosol,   naerosol, hygro      ) ! out


! call Razzak-Ghan activation routine with single updraft
   wbar = max( wup, 0.5_r8 )  ! force wbar >= 0.5 m/s for now
   sigw = 0.0_r8
   wdiab = 0.0_r8
   wminf = wbar
   wmaxf = wbar

   if (kk == kactfirst) then
! at cloud base - do primary activation
      call activate_modal(                                             &
         wbar, wmaxf, tair, rhoair,                                    & ! in
         naerosol, ntot_amode, vaerosol, hygro,                        & ! in
         fn, fm, fluxn, fluxm, flux_fullact                            ) ! out
   else
! above cloud base - do secondary activation with prescribed supersat
! that is constant with height
      call activate_modal(                                             &
         wbar, wmaxf, tair, rhoair,                                    & ! in
         naerosol, ntot_amode, vaerosol, hygro,                        & ! in
         fn, fm, fluxn, fluxm, flux_fullact,                           & ! out
         activate_smaxmax                                              ) ! optional in
   endif

! apply the activation fractions to the updraft aerosol mixing ratios
   dt_u_inv = 1.0_r8/dt_u

   do imode = 1, ntot_amode
        ! for aerosol number
        call assign_la_lc(imode, 0, la, lc)
        act_frac = fn(imode)
        call update_conu_from_act_frac ( conu,       dconudt,          & ! inout
                   pcnst_extd,  la,    lc,  act_frac,    dt_u_inv      ) ! in
        ! for aerosol mass
        do ispec = 1, nspec_amode(imode)
            call assign_la_lc(imode, ispec, la, lc)
            act_frac = fm(imode)
            call update_conu_from_act_frac ( conu,      dconudt,          & ! inout
                   pcnst_extd,  la,    lc,   act_frac,    dt_u_inv        ) ! in
        enddo  ! ispec
   enddo   ! imode

   return
   end subroutine ma_activate_convproc

!========================================================================================
   subroutine aer_vol_num_hygro( pcnst_extd, conu,     rhoair,    & ! in
                                 vaerosol,   naerosol, hygro      ) ! out
!-----------------------------------------------------------------------
! calculate aerosol volume, number and hygroscopicity
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! arguments:
   implicit none
   integer, intent(in)     :: pcnst_extd             ! TMR index
   ! conu = tracer mixing ratios in updraft at top of this (current) level. The
   ! conu are changed by activation
   real(r8), intent(in)    :: conu(pcnst_extd)       ! TMR [#/kg or kg/kg]
   real(r8), intent(in)    :: rhoair                 ! air density [kg/m3]
   real(r8), intent(out)   :: vaerosol(ntot_amode)   ! int+act volume [m3/m3]
   real(r8), intent(out)   :: naerosol(ntot_amode)   ! interstitial+activated number conc [#/m3]
   real(r8), intent(out)   :: hygro(ntot_amode)      ! current hygroscopicity for int+act [unitless]

! local variables
   integer  :: imode, ispec          ! index
   real(r8) :: tmp_vol_spec          ! aerosol volume for a single species [m3/kg]
   real(r8) :: tmp_vol               ! aerosol volume [m3/kg]
   real(r8) :: tmp_num               ! aerosol number [#/kg]
   real(r8) :: tmp_hygro             ! aerosol hygroscopicity * volume [m3/kg]
   real(r8) :: n_min, n_max          ! min and max bound of naerosol

   real(r8), parameter :: small_vol = 1.0e-35_r8   ! a small value that variables smaller than it are considered as zero for aerosol volume [m3/kg]

!-----------------------------------------------------------------------

   do imode = 1, ntot_amode

      ! compute aerosol (or a+cw) volume and hygroscopicity
      tmp_vol = 0.0_r8
      tmp_hygro = 0.0_r8
      do ispec = 1, nspec_amode(imode)
         tmp_vol_spec = max( conu(lmassptr_amode(ispec,imode)), 0.0_r8 ) / specdens_amode(lspectype_amode(ispec,imode))    ! mass divided by density
         tmp_vol = tmp_vol + tmp_vol_spec ! total aerosol volume
         tmp_hygro = tmp_hygro + tmp_vol_spec * spechygro(lspectype_amode(ispec,imode))   ! volume*hygro suming up for all species
      enddo
      vaerosol(imode) = tmp_vol * rhoair   ! change volume from m3/kgair to m3/m3air
      if (tmp_vol < small_vol) then
         hygro(imode) = 0.2_r8
      else
         hygro(imode) = tmp_hygro/tmp_vol
      endif

      ! computer a (or a+cw) number and bound it
      tmp_num = max( conu(numptr_amode(imode)), 0.0_r8 )
      n_min = vaerosol(imode)*voltonumbhi_amode(imode)
      n_max = vaerosol(imode)*voltonumblo_amode(imode)
      naerosol(imode) = min_max_bound(n_min, n_max, tmp_num*rhoair)

   enddo  ! imode

   return
   end subroutine aer_vol_num_hygro

!========================================================================================
   subroutine update_conu_from_act_frac ( conu,      dconudt,          & ! inout
                                        pcnst_extd,  la,    lc,        & ! in
                                        act_frac,    dt_u_inv          ) ! in
!-----------------------------------------------------------------------
! update conu and dconudt from activation fraction
!-----------------------------------------------------------------------
! arguments:
   integer, intent(in)     :: pcnst_extd
   real(r8), intent(inout) :: conu(pcnst_extd)    ! TMR concentration [#/kg or kg/kg]
   real(r8), intent(inout) :: dconudt(pcnst_extd) ! TMR tendencies due to activation [#/kg/s or kg/kg/s]
   real(r8), intent(in)    :: act_frac            ! activation fraction [fraction]
   real(r8), intent(in)    :: dt_u_inv            ! 1.0/dt_u  [1/s]
   integer, intent(in)     :: la                  ! indices for interstitial aerosols
   integer, intent(in)     :: lc                  ! indices for in-cloud water aerosols
! local variable
   real(r8)     :: delact   ! change of aerosols due to activation [#/kg or kg/kg]

   delact = min_max_bound(0.0_r8, conu(la), conu(la)*act_frac)
   ! update conu in interstitial and in-cloud condition
   conu(la) = conu(la) - delact
   conu(lc) = conu(lc) + delact
   ! update dconu/dt
   dconudt(la) = -delact*dt_u_inv
   dconudt(lc) =  delact*dt_u_inv

   return
   end subroutine update_conu_from_act_frac



!====================================================================================
   subroutine compute_wetdep_tend(                              &
                doconvproc_extd,        icol,   kk,     dt,     & ! in
                dt_u,   dp_i,   cldfrac_i,      mu_p_eudp,      & ! in
                aqfrac,         icwmr,          rprd,           & ! in
                conu,           dconudt_wetdep                  ) ! inout
!-----------------------------------------------------------------------
! compute tendency from wet deposition
! 
!    rprd               = precip formation as a grid-cell average (kgW/kgA/s)
!    icwmr              = cloud water MR within updraft area (kgW/kgA)
!    fupdr              = updraft fractional area (--)
!    A = rprd/fupdr     = precip formation rate within updraft area (kgW/kgA/s)
!    clw_preloss = cloud water MR before loss to precip
!                = icwmr + dt*(rprd/fupdr)
!    B = A/clw_preloss  = (rprd/fupdr)/(icwmr + dt*rprd/fupdr)
!                       = rprd/(fupdr*icwmr + dt*rprd)
!                       = first-order removal rate (1/s)
!    C = dp/(mup/fupdr) = updraft air residence time in the layer (s)
!
!    fraction removed = (1.0 - exp(-cdt)) where
!                 cdt = B*C = (fupdr*dp/mup)*[rprd/(fupdr*icwmr + dt*rprd)]
!
!    Note1:  *** cdt is now sensitive to fupdr, which we do not really know,
!                and is not the same as the convective cloud fraction
!    Note2:  dt is appropriate in the above cdt expression, not dtsub
!
!    Apply wet removal at levels where
!       icwmr(icol,k) > clw_cut  AND  rprd(icol,k) > 0.0
!    as wet removal occurs in both liquid and ice clouds
!-----------------------------------------------------------------------

   ! cloudborne aerosol, so the arrays are dimensioned with pcnst_extd = pcnst*2
   integer,  parameter  :: pcnst_extd = pcnst*2
   logical,  intent(in) :: doconvproc_extd(pcnst_extd) ! flag for doing convective transport
   integer,  intent(in) :: icol                 ! column index
   integer,  intent(in) :: kk                   ! vertical level index
   real(r8), intent(in) :: dt                   ! Model timestep [s]
   real(r8), intent(in) :: dt_u(pver)           ! lagrangian transport time in the updraft[s]
   real(r8), intent(in) :: dp_i(pver)            ! dp [mb]
   real(r8), intent(in) :: cldfrac_i(pver)      ! cldfrac at current icol (with adjustments) [fraction]
   real(r8), intent(in) :: mu_p_eudp(pver)      ! = mu_i(kp1) + eudp(k) [mb/s]
   real(r8), intent(in) :: aqfrac(pcnst_extd)   ! aqueous fraction of constituent in updraft [fraction]
   real(r8), intent(in) :: icwmr(pcols,pver)    ! Convective cloud water from zm scheme [kg/kg]
   real(r8), intent(in) :: rprd(pcols,pver)     ! Convective precipitation formation rate [kg/kg/s]
   real(r8), intent(inout) :: conu(pcnst_extd,pverp)   ! mix ratio in updraft at interfaces [kg/kg]
   real(r8), intent(inout) :: dconudt_wetdep(pcnst_extd,pverp) ! d(conu)/dt by wet removal[kg/kg/s]


   integer      :: icnst        ! index of pcnst_extd
   real(r8)     :: cdt          ! (in-updraft first order wet removal rate) * dt [unitless]
   real(r8)     :: expcdtm1     ! exp(cdt) - 1 [unitless]
   real(r8)     :: half_cld     ! 0.5 * cldfrac [fraction]

   real(r8), parameter :: clw_cut = 1.0e-6      ! cutoff value of cloud water for doing updraft [kg/kg]
                ! Skip levels where icwmr(icol,k) <= clw_cut (=1.0e-6) to
                ! eliminate occasional very small icwmr values from the ZM module

   cdt = 0.0_r8
   if ((icwmr(icol,kk) > clw_cut) .and. (rprd(icol,kk) > 0.0)) then
         half_cld = 0.5_r8*cldfrac_i(kk)
         cdt = (half_cld*dp_i(kk)/mu_p_eudp(kk)) * rprd(icol,kk) / &
                      (half_cld*icwmr(icol,kk) + dt*rprd(icol,kk))
   endif

   if (cdt > 0.0_r8) then
      expcdtm1 = exp(-cdt) - 1.0
      do icnst = 2, pcnst_extd
         if (doconvproc_extd(icnst)) then
            dconudt_wetdep(icnst,kk) = conu(icnst,kk)*aqfrac(icnst)*expcdtm1
            conu(icnst,kk) = conu(icnst,kk) + dconudt_wetdep(icnst,kk)
            dconudt_wetdep(icnst,kk) = dconudt_wetdep(icnst,kk) / dt_u(kk)
         endif
      enddo
   endif

   end subroutine compute_wetdep_tend

!====================================================================================
   subroutine compute_downdraft_mixing_ratio(              &
                doconvproc_extd,         ktop,  kbot,      & ! in
                md_i,           eddp,           const,     & ! in
                cond                                       ) ! inout
!-----------------------------------------------------------------------
! Compute downdraft mixing ratios from cloudtop to cloudbase
! No special treatment is needed at k=2
! No transformation or removal is applied in the downdraft
!-----------------------------------------------------------------------

   ! cloudborne aerosol, so the arrays are dimensioned with pcnst_extd = pcnst*2
   integer, parameter   :: pcnst_extd = pcnst*2
   logical, intent(in)  :: doconvproc_extd(pcnst_extd) ! flag for doing convective transport
   integer, intent(in)  :: ktop                     ! top level index
   integer, intent(in)  :: kbot                     ! bottom level index
   real(r8),intent(in)  :: md_i(pverp)              ! md at current i (note pverp dimension) [mb/s]
   real(r8),intent(in)  :: eddp(pver)               ! ed(i,k)*dp(i,k) at current i [mb/s]
   real(r8),intent(in)  :: const(pcnst_extd,pver)   ! gathered tracer array [kg/kg]
   real(r8),intent(inout) :: cond(pcnst_extd,pverp) ! mix ratio in downdraft at interfaces [kg/kg]

   integer      :: kk           ! vertical index
   integer      :: kp1          ! index of kk+1
   integer      :: icnst        ! index of pcnst_extd
   real(r8)     :: md_m_eddp    ! downdraft massflux at kp1 [mb/s]

      do kk = ktop, kbot
         kp1 = kk + 1
! md_m_eddp = downdraft massflux at kp1, without detrainment between k,kp1
         md_m_eddp = md_i(kk) - eddp(kk)
         if (md_m_eddp < -mbsth) then
            do icnst = 2, pcnst_extd
               if (doconvproc_extd(icnst)) then
                  cond(icnst,kp1) = ( md_i(kk)*cond(icnst,kk)             &
                                - eddp(kk)*const(icnst,kk) ) / md_m_eddp
               endif
            enddo
         endif
      enddo ! kk

   end subroutine compute_downdraft_mixing_ratio

!====================================================================================
   subroutine initialize_dcondt (                               &
                doconvproc_extd, iflux_method, ktop, kbot,      & ! in
                dpdry_i, fa_u,      mu_i,       md_i,           & ! in
                chat,    const,     conu,       cond,           & ! in
                dconudt_activa,     dconudt_wetdep,             & ! in
                dudp,    dddp,      eudp,       eddp,           & ! in
                dcondt                                          ) ! out
!-----------------------------------------------------------------------
! initialize dondt and update with aerosol activation and wetdeposition
! will update later with dcondt_prevap and dcondt_resusp
! NOTE:  The approach used in convtran applies to inert tracers and
!        must be modified to include source and sink terms
!-----------------------------------------------------------------------

   ! cloudborne aerosol, so the arrays are dimensioned with pcnst_extd = pcnst*2
   integer, parameter :: pcnst_extd = pcnst*2
   logical, intent(in)  :: doconvproc_extd(pcnst_extd) ! flag for doing convective transport
   integer, intent(in)  :: iflux_method             ! 1=as in convtran (deep), 2=uwsh
   integer, intent(in)  :: ktop                     ! top level index
   integer, intent(in)  :: kbot                     ! bottom level index
   real(r8),intent(in)  :: dpdry_i(pver)            ! dp [mb]
   real(r8),intent(in)  :: fa_u(pver)               ! fractional area of the updraft [fraction]
   real(r8),intent(in)  :: mu_i(pverp)              ! mu at current i (note pverp dimension, see ma_convproc_tend) [mb/s] 
   real(r8),intent(in)  :: md_i(pverp)              ! md at current i (note pverp dimension) [mb/s]
   real(r8),intent(in)  :: chat(pcnst_extd,pverp)   ! mix ratio in env at interfaces [kg/kg]
   real(r8),intent(in)  :: const(pcnst_extd,pver)   ! gathered tracer array [kg/kg]
   real(r8),intent(in)  :: conu(pcnst_extd,pverp)   ! mix ratio in updraft at interfaces [kg/kg]
   real(r8),intent(in)  :: cond(pcnst_extd,pverp)   ! mix ratio in downdraft at interfaces [kg/kg]
   real(r8),intent(in)  :: dconudt_activa(pcnst_extd,pverp) ! d(conu)/dt by activation [kg/kg/s]
   real(r8),intent(in)  :: dconudt_wetdep(pcnst_extd,pverp) ! d(conu)/dt by wet removal[kg/kg/s]
   real(r8),intent(in)  :: dudp(pver)           ! du(i,k)*dp(i,k) at current i [mb/s]
   real(r8),intent(in)  :: dddp(pver)           ! dd(i,k)*dp(i,k) at current i [mb/s]
   real(r8),intent(in)  :: eudp(pver)           ! eu(i,k)*dp(i,k) at current i [mb/s]
   real(r8),intent(in)  :: eddp(pver)           ! ed(i,k)*dp(i,k) at current i [mb/s]
   real(r8),intent(out)    :: dcondt(pcnst_extd,pver)  ! grid-average TMR tendency for current column  [kg/kg/s]


   integer      :: kk           ! vertical index
   integer      :: kp1, km1, kp1x, km1x ! index of kk+1 and kk-1. x: with bounds
   integer      :: icnst        ! index of pcnst_extd
   real(r8)     :: fa_u_dp      ! fa_u * dp at the current level [mb]
   real(r8)     :: fluxin, fluxout, netflux   ! in, out and net flux [kg/kg/s * mb]
   real(r8)     :: tmpa         ! working variable of mass flux [mb/s]
   real(r8)     :: netsrce      ! working variable of flux source [kg/kg/s * mb]

   ! initialize variables
   dcondt(:,:) = 0.0

   ! loop from ktop to kbot
   do kk = ktop, kbot
      kp1 = kk+1
      km1 = kk-1
      kp1x = min( kp1, pver )
      km1x = max( km1, 1 )
      fa_u_dp = fa_u(kk)*dpdry_i(kk)
      do icnst = 2, pcnst_extd
         if (doconvproc_extd(icnst)) then
            ! First compute fluxes using environment subsidence/lifting and
            ! entrainment/detrainment into up/downdrafts,
            ! to provide an additional mass balance check
            ! (this could be deleted after the code is well tested)
            fluxin  = mu_i(kk )*min(chat(icnst,kk ),const(icnst,km1x))     &
                    - md_i(kp1)*min(chat(icnst,kp1),const(icnst,kp1x))     &
                    + dudp(kk)*conu(icnst,kk) + dddp(kk)*cond(icnst,kp1)
            fluxout = mu_i(kp1)*min(chat(icnst,kp1),const(icnst,kk))       &
                    - md_i(kk )*min(chat(icnst,kk ),const(icnst,kk))       &
                    + (eudp(kk) + eddp(kk)) * const(icnst,kk)

            netflux = fluxin - fluxout

            ! Now compute fluxes as in convtran, and also source/sink terms
            ! (version 3 limit fluxes outside convection to mass in appropriate layer
            ! (these limiters are probably only safe for positive definite quantitities
            ! (it assumes that mu and md already satify a courant number limit of 1)
            if (iflux_method /= 2) then
               fluxin  =     mu_i(kp1)*conu(icnst,kp1)                          &
                           + mu_i(kk )*min(chat(icnst,kk ),const(icnst,km1x))   &
                         - ( md_i(kk )*cond(icnst,kk)                           &
                           + md_i(kp1)*min(chat(icnst,kp1),const(icnst,kp1x)) )
               fluxout =     mu_i(kk )*conu(icnst,kk)                           &
                           + mu_i(kp1)*min(chat(icnst,kp1),const(icnst,kk ))    &
                         - ( md_i(kp1)*cond(icnst,kp1)                          &
                           + md_i(kk )*min(chat(icnst,kk ),const(icnst,kk )) )
            else
               fluxin  =     mu_i(kp1)*conu(icnst,kp1) - ( md_i(kk )*cond(icnst,kk ) )
               fluxout =     mu_i(kk )*conu(icnst,kk ) - ( md_i(kp1)*cond(icnst,kp1) )

               ! new method -- simple upstream method for the env subsidence
               ! tmpa = net env mass flux (positive up) at top of layer k
               tmpa = -( mu_i(kk  ) + md_i(kk  ) )
               if (tmpa <= 0.0_r8) then
                  fluxin  = fluxin  - tmpa*const(icnst,km1x)
               else
                  fluxout = fluxout + tmpa*const(icnst,kk   )
               endif
               ! tmpa = net env mass flux (positive up) at base of layer k
               tmpa = -( mu_i(kp1) + md_i(kp1) )
               if (tmpa >= 0.0_r8) then
                  fluxin  = fluxin  + tmpa*const(icnst,kp1x)
               else
                  fluxout = fluxout - tmpa*const(icnst,kk   )
               endif
            endif
            netflux = fluxin - fluxout
            netsrce = fa_u_dp*(dconudt_activa(icnst,kk) + dconudt_wetdep(icnst,kk))

            dcondt(icnst,kk) = (netflux+netsrce)/dpdry_i(kk)

         endif   ! "doconvproc_extd"
      enddo      ! "icnst = 2,pcnst_extd"
   enddo ! "kk = ktop, kbot"

   end subroutine initialize_dcondt

!=========================================================================================
   subroutine ma_precpevap_convproc(                           &
              dcondt_wetdep,                                   & ! in
              rprd,    evapc,         dpdry_i,                 & ! in
              icol,    ktop,          pcnst_extd,              & ! in
              doconvproc_extd,        species_class,           & ! in
              dcondt_prevap,          dcondt_prevap_hist,      & ! out
              dcondt                                           & ! inout
              ) 
!-----------------------------------------------------------------------
!
! Purpose:
! Calculate resuspension of wet-removed aerosol species resulting precip evaporation
!
!     for aerosol mass   species, do non-linear resuspension to coarse mode
!     for aerosol number species, all the resuspension is done in wetdepa_v2, so do nothing here
!
! Author: R. Easter
!
!-----------------------------------------------------------------------


   implicit none

!-----------------------------------------------------------------------
! arguments
! (note:  TMR = tracer mixing ratio)
   integer,  intent(in)    :: pcnst_extd        ! aerosol dimension [unitless]
   real(r8), intent(inout) :: dcondt(pcnst_extd,pver)
                              ! overall TMR tendency from convection [kg/kg/s]
   real(r8), intent(in)    :: dcondt_wetdep(pcnst_extd,pver)
                              ! portion of TMR tendency due to wet removal [kg/kg/s]
   real(r8), intent(out)   :: dcondt_prevap(pcnst_extd,pver)
                              ! portion of TMR tendency due to precip evaporation [kg/kg/s]
                              ! (actually, due to the adjustments made here)
                              ! (on entry, this is 0.0)
   real(r8), intent(out)   :: dcondt_prevap_hist(pcnst_extd,pver) 
                              ! this determines what goes into the history 
                              !    precip-evap SFSEC variables
                              ! currently, the SFSEC resuspension are attributed
                              !    to the species that got scavenged,
                              !    WHICH IS NOT the species that actually
                              !    receives the resuspension
                              !    when modal_aero_wetdep_resusp_opt > 0
                              ! so when scavenged so4_c1 is resuspended as so4_a1, 
                              !    this resuspension column-tendency shows
                              !    up in so4_c1SFSES
                              ! this is done to allow better tracking of the
                              !    resuspension in the mass-budget post-processing scripts

   real(r8), intent(in)    :: rprd(pcols,pver)  ! conv precip production  rate (gathered) [kg/kg/s]
   real(r8), intent(in)    :: evapc(pcols,pver)  ! conv precip evaporation rate (gathered) [kg/kg/s]
   real(r8), intent(in)    :: dpdry_i(pver) ! pressure thickness of level [mb]

   integer,  intent(in)    :: icol  ! normal (ungathered) i index for current column
   integer,  intent(in)    :: ktop  ! index of top cloud level for current column

   logical,  intent(in)    :: doconvproc_extd(pcnst_extd)  ! indicates which species to process
   integer,  intent(in)    :: species_class(:)  ! specify what kind of species it is. defined at physconst.F90
                                                ! undefined  = 0
                                                ! cldphysics = 1
                                                ! aerosol    = 2
                                                ! gas        = 3
                                                ! other      = 4

!-----------------------------------------------------------------------
! local variables
   integer  :: kk
   real(r8) :: pr_flux               ! precip flux at base of current layer [(kg/kg/s)*mb]
   real(r8) :: pr_flux_tmp           ! precip flux at base of current layer, after adjustment of resuspension in step 1 [(kg/kg/s)*mb]
   real(r8) :: pr_flux_base          ! precip flux at an effective cloud base for calculations in a particular layer
   real(r8) :: wd_flux(pcnst_extd)   ! tracer wet deposition flux at base of current layer [(kg/kg/s)*mb]
   real(r8) :: x_ratio               ! ratio of adjusted and old fraction of precipitation-borne aerosol flux that is NOT resuspended, calculated in step 1 and used in step 2 (see below)

!-----------------------------------------------------------------------

!
! *** note use of non-standard units
!
! precip
!    dpdry_i is mb
!    rprd and evapc are kgwtr/kgair/s
!    pr_flux = dpdry_i(kk)*rprd is mb*kgwtr/kgair/s
!
! precip-borne aerosol
!    dcondt_wetdep is kgaero/kgair/s
!    wd_flux = tmpdp*dcondt_wetdep is mb*kgaero/kgair/s
!    dcondt_prevap = del_wd_flux_evap/dpdry_i is kgaero/kgair/s
! so this works ok too
!
! *** dilip switched from tmpdg (or dpdry_i) to tmpdpg = tmpdp/gravit
! that is incorrect, but probably does not matter
!    for aerosol, wd_flux units do not matter
!        only important thing is that tmpdp (or tmpdpg) is used
!        consistently when going from dcondt to wd_flux then to dcondt

   ! initiate variables that are integrated in vertical
   pr_flux = 0.0_r8
   pr_flux_base = 0.0_r8
   wd_flux(:) = 0.0_r8
   dcondt_prevap(:,:) = 0.0
   dcondt_prevap_hist(:,:) = 0.0

   do kk = ktop, pver
! step 1 - precip evaporation and aerosol resuspension
      call ma_precpevap(                                        &
                            icol,     kk,    dpdry_i,   evapc,  &  ! in
                            pr_flux,                            &  ! in
                            pr_flux_base,                       &  ! inout
                            pr_flux_tmp,     x_ratio            &  ! out
                        )

! step 2 - precip production and aerosol scavenging
      call  ma_precpprod(                                         &
              rprd,   dpdry_i,  icol,    kk,    pcnst_extd,       & ! in
              doconvproc_extd,  x_ratio,        species_class,    & ! in
              pr_flux,        pr_flux_tmp,       pr_flux_base,    & ! inout
              wd_flux,        dcondt_wetdep,                      & ! inout
              dcondt,         dcondt_prevap,  dcondt_prevap_hist  & ! inout
              )

   enddo ! kk

   return
   end subroutine ma_precpevap_convproc

!=========================================================================================
   subroutine ma_precpevap(                                     &
                            icol,     kk,    dpdry_i,  evapc,   & ! in
                            pr_flux,                            & ! in
                            pr_flux_base,                       & ! inout
                            pr_flux_tmp,     x_ratio            ) ! out
!------------------------------------------
! step 1 in ma_precpevap_convproc: aerosol resuspension from precipitation evaporation
!------------------------------------------

   integer,  intent(in)      :: icol  ! normal (ungathered) i index for current column
   integer,  intent(in)      :: kk    ! vertical level
   real(r8), intent(in)      :: evapc(pcols,pver) ! conv precipitataion evaporation rate [kg/kg/s]
   real(r8), intent(in)      :: dpdry_i(pver)     ! pressure thickness of level [mb]
   real(r8), intent(in)      :: pr_flux           ! precip flux at base of current layer [(kg/kg/s)*mb]
   real(r8), intent(inout)   :: pr_flux_base      ! precip flux at an effective cloud base for calculations in a particular layer
   real(r8), intent(out)     :: pr_flux_tmp       ! precip flux at base of current layer, after adjustment in step 1 [(kg/kg/s)*mb]
   real(r8), intent(out)     :: x_ratio           ! ratio of adjusted and old fraction of precipitation-borne aerosol flux that is NOT resuspended, used in step 2

   ! local variables
   real(r8) :: ev_flux_local ! local precip flux from evaporation [(kg/kg/s)*mb]
   real(r8) :: pr_ratio_old, pr_ratio_tmp  ! ratio of pr_flux and pr_flux_base, before and after adjustment in step 1
   real(r8) :: frac_aer_resusp_old, frac_aer_resusp_tmp  ! fraction of precipitation-borne aerosol flux that is NOT resuspended, before and after adjustment in step 1

   real(r8), parameter :: small_value = 1.0e-30_r8   ! a small value that variables smaller than it are considered as zero


   ! adjust pr_flux due to local evaporation
   ev_flux_local = max( 0.0_r8, evapc(icol,kk)*dpdry_i(kk) )
   pr_flux_tmp = min_max_bound(0.0_r8, pr_flux_base, pr_flux-ev_flux_local)

   x_ratio = 0.0_r8
   if (pr_flux_base < small_value) then
         pr_flux_base = 0.0_r8    ! this will start things fresh at the next layer
         pr_flux_tmp  = 0.0_r8
         return
   endif

   ! calculate fraction of resuspension
   pr_ratio_old = pr_flux/pr_flux_base
   pr_ratio_old = min_max_bound(0.0_r8, 1.0_r8, pr_ratio_old)
   frac_aer_resusp_old = 1.0_r8 - faer_resusp_vs_fprec_evap_mpln(1.0_r8-pr_ratio_old, 2)  ! 2: log-normal distribution
   frac_aer_resusp_old = min_max_bound(0.0_r8, 1.0_r8, frac_aer_resusp_old)

   pr_ratio_tmp = min_max_bound(0.0_r8, 1.0_r8, pr_flux_tmp/pr_flux_base)
   pr_ratio_tmp = min( pr_ratio_tmp, pr_ratio_old )
   frac_aer_resusp_tmp = 1.0_r8 - faer_resusp_vs_fprec_evap_mpln(1.0_r8-pr_ratio_tmp, 2)
   frac_aer_resusp_tmp = min_max_bound(0.0_r8, 1.0_r8, frac_aer_resusp_tmp)
   frac_aer_resusp_tmp = min( frac_aer_resusp_tmp, frac_aer_resusp_old )

   !compute x_ratio
   if (frac_aer_resusp_tmp > small_value) then
          x_ratio = frac_aer_resusp_tmp/frac_aer_resusp_old
   else   ! this will start things fresh at the next layer
          pr_flux_base = 0.0_r8
          pr_flux_tmp  = 0.0_r8
   endif

   return
   end subroutine ma_precpevap

!=========================================================================================
   subroutine ma_precpprod(                                       &
              rprd,   dpdry_i,  icol,   kk,     pcnst_extd,       &  ! in
              doconvproc_extd,  x_ratio,        species_class,    &  ! in
              pr_flux,        pr_flux_tmp,       pr_flux_base,    &  ! inout
              wd_flux,    dcondt_wetdep,                          &  ! inout
              dcondt,    dcondt_prevap,  dcondt_prevap_hist       &  ! inout                            
              ) 
!------------------------------------------
! step 2 in ma_precpevap_convproc: aerosol scavenging from precipitation production
!------------------------------------------

   integer,  intent(in)    :: icol  ! normal (ungathered) i index for current column
   integer,  intent(in)    :: kk
   integer,  intent(in)    :: pcnst_extd        ! aerosol dimension [unitless]
   real(r8), intent(in)    :: rprd(pcols,pver) ! conv precip production  rate (at a certain level) [kg/kg/s]
   real(r8), intent(in)    :: dcondt_wetdep(pcnst_extd,pver) ! portion of TMR tendency due to wet removal [kg/kg/s]
   real(r8), intent(in)    :: dpdry_i(pver)      ! pressure thickness of level [mb]
   logical,  intent(in)    :: doconvproc_extd(pcnst_extd)  ! indicates which species to process
   real(r8), intent(in)    :: x_ratio    ! ratio of adjusted and old fraction of precipitation-borne aerosol flux that is NOT resuspended, calculated in step 1
   integer,  intent(in)    :: species_class(:)

   real(r8), intent(inout)    :: pr_flux   ! precip flux at base of current layer [(kg/kg/s)*mb]
   real(r8), intent(inout)    :: pr_flux_tmp   ! precip flux at base of current layer, after adjustment in step 1 [(kg/kg/s)*mb]
   real(r8), intent(inout)    :: pr_flux_base   ! precip flux at an effective cloud base for calculations in a particular layer
   real(r8), intent(inout)    :: wd_flux(pcnst_extd)   ! tracer wet deposition flux at base of current layer [(kg/kg/s)*mb]
   real(r8), intent(inout)    :: dcondt(pcnst_extd,pver)  ! overall TMR tendency from convection at a certain layer [kg/kg/s]
   real(r8), intent(inout)    :: dcondt_prevap(pcnst_extd,pver)  ! portion of TMR tendency due to precip evaporation [kg/kg/s]
   real(r8), intent(inout)    :: dcondt_prevap_hist(pcnst_extd,pver)   ! dcondt_prevap_hist at a certain layer [kg/kg/s]

   ! local variables
   real(r8) :: pr_flux_local            ! local precip flux [(kg/kg/s)*mb]
   real(r8) :: wd_flux_local            ! local wet deposition flux [(kg/kg/s)*mb]
   real(r8) :: wd_flux_tmp              ! updated wet deposition flux [(kg/kg/s)*mb]
   real(r8) :: del_wd_flux_evap         ! change to wet deposition flux from evaporation [(kg/kg/s)*mb]
   real(r8) :: dcondt_wdflux            ! dcondt due to wet deposition flux change [kg/kg/s]
   integer  :: icnst, icnst2, mmtoo


      pr_flux_local = max( 0.0_r8, rprd(icol,kk)*dpdry_i(kk) )
      pr_flux_base = max( 0.0_r8, pr_flux_base + pr_flux_local )
      pr_flux = min_max_bound(0.0_r8, pr_flux_base, pr_flux_tmp+pr_flux_local)

      do icnst = 2, pcnst_extd
         if ( .not. doconvproc_extd(icnst) ) cycle

        ! wet deposition flux from the aerosol resuspension
        ! wd_flux_tmp (updated) = (wd_flux coming into the layer) - (resuspension ! decrement)
         wd_flux_tmp = max( 0.0_r8, wd_flux(icnst) * x_ratio )
         del_wd_flux_evap = max( 0.0_r8, wd_flux(icnst) - wd_flux_tmp )

        ! wet deposition flux from the aerosol scavenging
        ! wd_flux (updated) = (wd_flux after resuspension) - (scavenging ! increment)
         wd_flux_local  = max( 0.0_r8, -dcondt_wetdep(icnst,kk)*dpdry_i(kk) )
         wd_flux(icnst) = max( 0.0_r8, wd_flux_tmp + wd_flux_local )

         dcondt_wdflux = del_wd_flux_evap/dpdry_i(kk)

         icnst2 = mod( icnst-1, pcnst ) + 1 ! for interstitial icnst2=icnst;  for activated icnst2=icnst-pcnst

         ! not sure what this mean exactly. Only do it for aerosol mass species
         ! (mmtoo>0).  mmtoo<=0 represents aerosol number species
         mmtoo = mmtoo_prevap_resusp(icnst2)
         if ( species_class(icnst2) == spec_class_aerosol ) then
            if (mmtoo > 0) then
               ! add the precip-evap (resuspension) to the history-tendency of the current species
               dcondt_prevap_hist(icnst,kk) = dcondt_prevap_hist(icnst,kk) + dcondt_wdflux
               ! add the precip-evap (resuspension) to the actual tendencies of appropriate coarse-mode species
               dcondt_prevap(mmtoo,kk) = dcondt_prevap(mmtoo,kk) + dcondt_wdflux
               dcondt(mmtoo,kk) = dcondt(mmtoo,kk) + dcondt_wdflux
            endif

         else ! ( species_class(m2) /= spec_class_aerosol )
            ! do this for trace gases (although currently modal_aero_convproc does not treat trace gases)
            dcondt_prevap_hist(icnst,kk) = dcondt_prevap_hist(icnst,kk) + dcondt_wdflux
            dcondt_prevap(icnst,kk) = dcondt_prevap(icnst,kk) + dcondt_wdflux
            dcondt(icnst,kk) = dcondt(icnst,kk) + dcondt_wdflux
         endif

      enddo ! m = 2, pcnst_extd

    end subroutine ma_precpprod

!=========================================================================================
   subroutine ma_resuspend_convproc(   dcondt,          & ! inout
              ktop,  kbot_prevap,  pcnst_extd,          & ! in
              dcondt_resusp                             ) ! out

!-----------------------------------------------------------------------
!
! Purpose:
! Calculate resuspension of activated aerosol species resulting from both
!    detrainment from updraft and downdraft into environment
!    subsidence and lifting of environment, which may move air from
!       levels with large-scale cloud to levels with no large-scale cloud
!
! Method:
! Three possible approaches were considered:
!
! 1. Ad-hoc #1 approach.  At each level, adjust dcondt for the activated
!    and unactivated portions of a particular aerosol species so that the
!    ratio of dcondt (activated/unactivate) is equal to the ratio of the
!    mixing ratios before convection.
!    THIS WAS IMPLEMENTED IN MIRAGE2
!
! 2. Ad-hoc #2 approach.  At each level, adjust dcondt for the activated
!    and unactivated portions of a particular aerosol species so that the
!    change to the activated portion is minimized (zero if possible).  The
!    would minimize effects of convection on the large-scale cloud.
!    THIS IS CURRENTLY IMPLEMENTED IN CAM5 where we assume that convective
!    clouds have no impact on the stratiform-cloudborne aerosol
!
! 3. Mechanistic approach that treats the details of interactions between
!    the large-scale and convective clouds.  (Something for the future.)
!
! Author: R. Easter
!
! C++ porting: only method #2 is implemented.
!-----------------------------------------------------------------------

   implicit none

!-----------------------------------------------------------------------
! arguments
! (note:  TMR = tracer mixing ratio)
   integer,  intent(in)    :: pcnst_extd ! cloudborne aerosol dimension [unitless]
   real(r8), intent(inout) :: dcondt(pcnst_extd,pver)
                              ! overall TMR tendency from convection [#/kg/s or kg/kg/s]
   real(r8), intent(out)   :: dcondt_resusp(pcnst_extd,pver)
                              ! portion of TMR tendency due to resuspension [#/kg/s or kg/kg/s]
                              ! (actually, due to the adjustments made here)

   integer,  intent(in)    :: ktop, kbot_prevap ! indices of top and bottom cloud levels


!-----------------------------------------------------------------------
! local variables
   integer  :: imode, ispec, la, lc      ! indices
!-----------------------------------------------------------------------

   dcondt_resusp(:,:) = 0.0

   do imode = 1, ntot_amode
      do ispec = 0, nspec_amode(imode)
         call assign_la_lc(imode, ispec, la, lc)
         call tmr_tendency(pcnst_extd, dcondt, dcondt_resusp, la, lc, ktop, kbot_prevap)
      enddo   ! "ll = -1, nspec_amode(n)"
   enddo      ! "n = 1, ntot_amode"


   return
   end subroutine ma_resuspend_convproc

!=========================================================================================
   subroutine tmr_tendency(pcnst_extd, dcondt, dcondt_resusp, &
                        la, lc, ktop, kbot_prevap)
!-----------------------------------------------------------------------
! calculate tendency of TMR
!-----------------------------------------------------------------------

   implicit none

   ! arguments (note:  TMR = tracer mixing ratio)
   integer,  intent(in)    :: pcnst_extd ! cloudborne aerosol dimension [unitless]
   real(r8), intent(inout) :: dcondt(pcnst_extd,pver) ! overall TMR tendency from convection [#/kg/s or kg/kg/s]
   real(r8), intent(inout) :: dcondt_resusp(pcnst_extd,pver) ! portion of TMR tendency due to resuspension [#/kg/s or kg/kg/s]
   integer, intent(in)     :: la, lc            ! indices
   integer,  intent(in)    :: ktop, kbot_prevap ! indices of top and bottom cloud levels

   ! local variables
   integer  :: kk   ! indices of vertical levels

   ! only apply adjustments to dcondt for pairs of unactivated (la) and
   ! activated (lc) aerosol species
   if ((la > 0) .and. (la <= pcnst_extd) .and. (lc > 0) .and. (lc <= pcnst_extd)) then
       do kk = ktop, kbot_prevap
            ! cam5 approach
            dcondt(la,kk) = dcondt(la,kk) + dcondt(lc,kk)
            dcondt_resusp(la,kk) = dcondt(lc,kk)
            dcondt_resusp(lc,kk) =  - dcondt(lc,kk)
            dcondt(lc,kk) = 0.0
       enddo
   endif

   return
   end subroutine tmr_tendency


!=========================================================================================
   subroutine compute_tendency_resusp_evap(                             &
                doconvproc_extd, ktop,          kbot_prevap,  dpdry_i,  & ! in
                dcondt_resusp,  dcondt_prevap,  dcondt_prevap_hist,     & ! in
                dconudt_activa, dconudt_wetdep, fa_u,                   & ! in
                sumactiva,      sumaqchem,      sumwetdep,              & ! out
                sumresusp,      sumprevap,      sumprevap_hist          ) ! out
!-----------------------------------------------------------------------
! calculate column-tendency of resuspension and evaporation
!-----------------------------------------------------------------------

   ! cloudborne aerosol, so the arrays are dimensioned with pcnst_extd = pcnst*2
   integer, parameter :: pcnst_extd = pcnst*2
   logical, intent(in)  :: doconvproc_extd(pcnst_extd) ! flag for doing convective transport
   integer, intent(in)  :: ktop                     ! top level index
   integer, intent(in)  :: kbot_prevap              ! bottom level index, for resuspension and evaporation only
   real(r8),intent(in)  :: dpdry_i(pver)            ! dp [mb]
   real(r8),intent(in)  :: fa_u(pver)               ! fractional area of the updraft [fraction]
   real(r8),intent(in)  :: dconudt_activa(pcnst_extd,pverp) ! d(conu)/dt by activation [kg/kg/s]
   real(r8),intent(in)  :: dconudt_wetdep(pcnst_extd,pverp) ! d(conu)/dt by wet removal[kg/kg/s]
   real(r8),intent(in)  :: dcondt_resusp(pcnst_extd,pver)  ! portion of TMR tendency due to resuspension [kg/kg/s]
   real(r8),intent(in)  :: dcondt_prevap(pcnst_extd,pver)  ! portion of TMR tendency due to precip evaporation [kg/kg/s]
   real(r8),intent(in)  :: dcondt_prevap_hist(pcnst_extd,pver) ! portion of TMR tendency due to precip evaporation, goes into the history [kg/kg/s]

   real(r8),intent(out)  :: sumactiva(pcnst_extd)    ! sum (over layers) of dp*dconudt_activa [kg/kg/s * mb]
   real(r8),intent(out)  :: sumaqchem(pcnst_extd)    ! sum (over layers) of dp*dconudt_aqchem [kg/kg/s * mb]
   real(r8),intent(out)  :: sumwetdep(pcnst_extd)    ! sum (over layers) of dp*dconudt_wetdep [kg/kg/s * mb]
   real(r8),intent(out)  :: sumresusp(pcnst_extd)    ! sum (over layers) of dp*dconudt_resusp [kg/kg/s * mb]
   real(r8),intent(out)  :: sumprevap(pcnst_extd)    ! sum (over layers) of dp*dconudt_prevap [kg/kg/s * mb]
   real(r8),intent(out)  :: sumprevap_hist(pcnst_extd)    ! sum (over layers) of dp*dconudt_prevap_hist [kg/kg/s * mb]

   integer      :: kk           ! vertical index
   integer      :: icnst        ! index of pcnst_extd
   real(r8)     :: dconudt_aqchem(pcnst_extd,pverp)  ! d(conu)/dt by aqueous chem [kg/kg/s]


   ! initialize variables
   sumactiva(:) = 0.0
   sumaqchem(:) = 0.0
   sumwetdep(:) = 0.0
   sumresusp(:) = 0.0
   sumprevap(:) = 0.0
   sumprevap_hist(:) = 0.0
   dconudt_aqchem(:,:) = 0.0    ! aqueous chemistry is ignored in current code
 
   do icnst = 2, pcnst_extd
      if (doconvproc_extd(icnst)) then
         ! should go to kk=pver for dcondt_prevap, and this should be safe for other sums
         do kk = ktop, kbot_prevap
            sumactiva(icnst) = sumactiva(icnst) + dconudt_activa(icnst,kk)*dpdry_i(kk)*fa_u(kk)
            sumaqchem(icnst) = sumaqchem(icnst) + dconudt_aqchem(icnst,kk)*dpdry_i(kk)*fa_u(kk)
            sumwetdep(icnst) = sumwetdep(icnst) + dconudt_wetdep(icnst,kk)*dpdry_i(kk)*fa_u(kk)

            sumresusp(icnst) = sumresusp(icnst) + dcondt_resusp(icnst,kk)*dpdry_i(kk)
            sumprevap(icnst) = sumprevap(icnst) + dcondt_prevap(icnst,kk)*dpdry_i(kk)
            sumprevap_hist(icnst) = sumprevap_hist(icnst) + dcondt_prevap_hist(icnst,kk)*dpdry_i(kk)
         enddo
      endif
   enddo
   end subroutine compute_tendency_resusp_evap

!====================================================================================
   subroutine update_tendency_final(                            &
                        ktop,   kbot_prevap,    ntsub,  jtsub,  & ! in
                        ncnst,  nsrflx,         icol,           & ! in
                        dt,     dcondt,         doconvproc,     & ! in
                        sumactiva,              sumaqchem,      & ! inout
                        sumwetdep,              sumresusp,      & ! inout
                        sumprevap,              sumprevap_hist, & ! inout
                        dqdt,   q_i,            qsrflx          ) ! inout
!-----------------------------------------------------------------------
! update tendencies to final output of ma_convproc_tend
!
! note that the ma_convproc_tend does not apply convective cloud processing
!    to the stratiform-cloudborne aerosol
! within this routine, cloudborne aerosols are convective-cloudborne
!
! before tendencies (dcondt, which is loaded into dqdt) are returned,
!    the convective-cloudborne aerosol tendencies must be combined
!    with the interstitial tendencies
! ma_resuspend_convproc has already done this for the dcondt
!
! the individual process column tendencies (sumwetdep, sumprevap, ...)
!    are just diagnostic fields that can be written to history
! tendencies for interstitial and convective-cloudborne aerosol could
!    both be passed back and output, if desired
! currently, however, the interstitial and convective-cloudborne tendencies
!    are combined (in the next code block) before being passed back (in qsrflx)
!-----------------------------------------------------------------------

   integer, parameter   :: pcnst_extd = pcnst*2
   integer, intent(in)  :: ktop              ! top level index
   integer, intent(in)  :: kbot_prevap       ! bottom level index, for resuspension and evaporation only
   integer, intent(in)  :: ntsub             ! number of sub timesteps
   integer, intent(in)  :: jtsub             ! index of sub timesteps from the outer loop
   integer, intent(in)  :: ncnst             ! number of tracers to transport
   integer, intent(in)  :: nsrflx            ! last dimension of qsrflx
   integer, intent(in)  :: icol              ! column index
   real(r8),intent(in)  :: dt                ! delta t (model time increment) [s]
   real(r8),intent(in)  :: dcondt(pcnst_extd,pver)  ! grid-average TMR tendency for current column  [kg/kg/s]
   logical, intent(in)  :: doconvproc(ncnst) ! flag for doing convective transport
   real(r8), intent(inout) :: sumactiva(pcnst_extd)            ! sum (over layers) of dp*dconudt_activa [kg/kg/s * mb]
   real(r8), intent(inout) :: sumaqchem(pcnst_extd)            ! sum (over layers) of dp*dconudt_aqchem [kg/kg/s * mb]
   real(r8), intent(inout) :: sumwetdep(pcnst_extd)            ! sum (over layers) of dp*dconudt_wetdep [kg/kg/s * mb]
   real(r8), intent(inout) :: sumresusp(pcnst_extd)     ! sum (over layers) of dp*dcondt_resusp [kg/kg/s * mb]
   real(r8), intent(inout) :: sumprevap(pcnst_extd)     ! sum (over layers) of dp*dcondt_prevap [kg/kg/s * mb]
   real(r8), intent(inout) :: sumprevap_hist(pcnst_extd)! sum (over layers) of dp*dcondt_prevap_hist [kg/kg/s * mb]
   real(r8), intent(inout) :: dqdt(pcols,pver,ncnst)    ! Tracer tendency array
   real(r8), intent(inout) :: q_i(pver,pcnst)           ! q(icol,kk,icnst) at current icol
   real(r8), intent(inout) :: qsrflx(pcols,pcnst,nsrflx)
                              ! process-specific column tracer tendencies
                              !  1 = activation   of interstial to conv-cloudborne
                              !  2 = resuspension of conv-cloudborne to interstital
                              !  3 = aqueous chemistry (not implemented yet, so zero)
                              !  4 = wet removal
                              !  5 = actual precip-evap resuspension (what actually is applied to a species)
                              !  6 = pseudo precip-evap resuspension (for history file)

   ! local variables
   integer  :: la, lc, imode, ispec             ! indices
   integer  :: icnst                            ! index
   integer  :: kk                               ! vertical index
   real(r8) :: dtsub                            ! delta t of sub timestep (dt/ntsub) [s]
   real(r8) :: xinv_ntsub                       ! inverse of ntsub (1.0/ntsub)
   real(r8) :: dqdt_i(pver,pcnst)               ! dqdt(icol,kk,icnst) at current icol
   real(r8) :: qsrflx_i(pcnst,nsrflx)           ! qsrflx(i,m,n) at current i

   ! initiate variables
   qsrflx_i(:,:) = 0.0
   dqdt_i(:,:) = 0.0
   xinv_ntsub = 1.0_r8/ntsub
   dtsub = dt*xinv_ntsub


   do imode = 1, ntot_amode
        do ispec = 0, nspec_amode(imode)
             call assign_la_lc(imode, ispec, la, lc)
             if (doconvproc(la)) then
                   sumactiva(la) = sumactiva(la) + sumactiva(lc)
                   sumresusp(la) = sumresusp(la) + sumresusp(lc)
                   sumaqchem(la) = sumaqchem(la) + sumaqchem(lc)
                   sumwetdep(la) = sumwetdep(la) + sumwetdep(lc)
                   sumprevap(la) = sumprevap(la) + sumprevap(lc)
                   sumprevap_hist(la) = sumprevap_hist(la) + sumprevap_hist(lc)
             endif
        enddo ! ispec
   enddo ! imode

   ! scatter overall tendency back to full array
   do icnst = 2, ncnst
        if (doconvproc(icnst)) then
            ! scatter overall dqdt tendency back
            do kk = ktop, kbot_prevap  ! should go to k=pver because of prevap
               dqdt_i(kk,icnst) = dcondt(icnst,kk)
               dqdt(icol,kk,icnst) = dqdt(icol,kk,icnst) + dqdt_i(kk,icnst)*xinv_ntsub
               ! update the q_i for the next interation of the jtsub loop
               if (jtsub < ntsub) then
                  q_i(kk,icnst) = max( (q_i(kk,icnst) + dqdt_i(kk,icnst)*dtsub), 0.0_r8 )
               endif
            enddo

            ! scatter column burden tendencies for various processes to qsrflx
            qsrflx_i(icnst,1) = sumactiva(icnst)*hund_ovr_g
            qsrflx_i(icnst,2) = sumresusp(icnst)*hund_ovr_g
            qsrflx_i(icnst,3) = sumaqchem(icnst)*hund_ovr_g
            qsrflx_i(icnst,4) = sumwetdep(icnst)*hund_ovr_g
            qsrflx_i(icnst,5) = sumprevap(icnst)*hund_ovr_g
            qsrflx_i(icnst,6) = sumprevap_hist(icnst)*hund_ovr_g
            qsrflx(icol,icnst,1:6) = qsrflx(icol,icnst,1:6) + qsrflx_i(icnst,1:6)*xinv_ntsub
       endif
   enddo ! icnst

   end subroutine update_tendency_final


!=========================================================================================
   subroutine assign_la_lc( imode,      ispec,          & ! in
                            la,         lc              ) ! out
!-----------------------------------------------------------------------
! get the index of interstital (la) and cloudborne (lc) aerosols
! from mode index and species index
!-----------------------------------------------------------------------

   integer, intent(in)     :: imode            ! index of MAM4 modes
   integer, intent(in)     :: ispec            ! index of species, in which:
                                               ! 0 = number concentration
                                               ! other = mass concentration
   integer, intent(out)    :: la               ! index of interstitial aerosol
   integer, intent(out)    :: lc               ! index of cloudborne aerosol (la + pcnst)

   if (ispec == 0) then
      la = numptr_amode(imode)
      lc = numptrcw_amode(imode) + pcnst
   else
      la = lmassptr_amode(ispec,imode)
      lc = lmassptrcw_amode(ispec,imode) + pcnst
   endif

   end subroutine assign_la_lc

!=========================================================================================



end module modal_aero_convproc
