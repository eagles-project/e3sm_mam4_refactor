
!--------------------------------------------------------------------------------
! Call aerosol microphysics processes for a single (cloudy or clear) subarea
! (with indices = lchnk,ii,kk,jsubarea)
!
! qgas3, qaer3, qaercw3, qnum3, qnumcw3 are the current incoming TMRs
! qgas_cur, qaer_cur, qaercw_cur, qnum_cur, qnumcw_cur are the updated outgoing TMRs
!
! In a clear subarea, calculate
!  - gas-aerosol exchange (condensation/evaporation)
!  - growth from smaller to larger modes (renaming) due to condensation
!  - new particle nucleation
!  - coagulation
!  - transfer of particles from hydrophobic modes to hydrophilic modes (aging)
!    due to condensation and coagulation
!
! In a cloudy subarea,
!  - when do_cond = false, this routine only calculate changes involving
!    growth from smaller to larger modes (renaming) following cloud chemistry
!    so gas TMRs are not changed
!  - when do_cond = true, this routine also calculates changes involving
!    gas-aerosol exchange (condensation/evaporation)
!  - transfer of particles from hydrophobic modes to hydrophilic modes (aging)
!       due to condensation
! Currently, in a cloudy subarea, this routine does not do
!  - new particle nucleation - because h2so4 gas conc. should be very low in cloudy air
!  - coagulation - because cloud-borne aerosol would need to be included
!--------------------------------------------------------------------------------
subroutine mam_amicphys_1subarea(&
                                 do_cond_sub,            do_rename_sub,      &
                                 do_newnuc_sub,          do_coag_sub,        &
                                 nstep,      lchnk,      ii,      kk,        &
                                 latndx,     lonndx,     lund,               &
                                 loffset,    deltat,                         &
                                 jsubarea,               nsubarea,           &
                                 iscldy_subarea,         afracsub,           &
                                 temp,       pmid,       pdel,               &
                                 zmid,       pblh,       relhum,             &
                                 dgn_a,      dgn_awet,   wetdens,            &
                                 qgas1,   qgas3,   qgas_cur,   qgas_delaa,   &
                                          qnum3,   qnum_cur,   qnum_delaa,   &
                                 qaer2,   qaer3,   qaer_cur,   qaer_delaa,   &
                                          qwtr3,   qwtr_cur,                 &
                                          qnumcw3, qnumcw_cur, qnumcw_delaa, &
                                 qaercw2, qaercw3, qaercw_cur, qaercw_delaa, &
                                 misc_vars_aa_sub                            )

   use modal_aero_amicphys_control, only: r8, ntot_poaspec=>npoa, ntot_soaspec=>nsoa &
                                        , max_mode, max_gas, max_aer &
                                        , misc_vars_aa_type &
                                        , gaexch_h2so4_uptake_optaa, ngas, igas_h2so4, igas_nh3 &
                                        , nait, nacc, max_agepair, n_agepair
   use modal_aero_amicphys_diags,   only: nqtendaa, nqqcwtendaa &
                                        , iqtend_cond, iqtend_rnam, iqtend_nnuc, iqtend_coag, iqqcwtend_rnam

   use modal_aero_data,   only: n_mode=>ntot_amode
   use physconst,         only: r_universal
   use modal_aero_coag,   only: mam_coag_1subarea
   use modal_aero_rename, only: mam_rename_1subarea

   logical,  intent(in)    :: do_cond_sub, do_rename_sub    ! true if the aero. microp. process is calculated in this subarea
   logical,  intent(in)    :: do_newnuc_sub, do_coag_sub    ! true if the aero.  microp. process is calculated in this subarea
   logical,  intent(in)    :: iscldy_subarea        ! true if sub-area is cloudy
   integer,  intent(in)    :: lchnk                 ! chunk identifier
   integer,  intent(in)    :: nstep                 ! model time-step number
   integer,  intent(in)    :: ii, kk                ! column and level indices
   integer,  intent(in)    :: latndx, lonndx        ! lat and lon indices
   integer,  intent(in)    :: lund                  ! logical unit for diagnostic output
   integer,  intent(in)    :: loffset
   integer,  intent(in)    :: jsubarea, nsubarea    ! sub-area index, number of sub-areas

   real(r8), intent(in)    :: afracsub              ! fractional area of subarea [unitless, 0-1]
   real(r8), intent(in)    :: deltat                ! time step [s]

   real(r8), intent(in)    :: temp                  ! air temperature at model levels [K]
   real(r8), intent(in)    :: pmid                  ! air pressure at layer center [Pa]
   real(r8), intent(in)    :: pdel                  ! pressure thickness of layer [Pa]
   real(r8), intent(in)    :: zmid                  ! altitude (above ground) at layer center [m]
   real(r8), intent(in)    :: pblh                  ! planetary boundary layer depth [m]
   real(r8), intent(in)    :: relhum                ! relative humidity [unitless, 0-1]

   real(r8), intent(inout) :: dgn_a   (max_mode)    ! dry geo. mean diameter [m] of number distribution
   real(r8), intent(inout) :: dgn_awet(max_mode)    ! wet geo. mean diameter [m] of number distribution
   real(r8), intent(inout) :: wetdens (max_mode)    ! interstitial aerosol wet density [kg/m3]

   ! Subare mixing ratios qXXXN (X=gas,aer,wat,num; N=1:4):
   !
   !    XXX=gas - gas species [kmol/kmol]
   !    XXX=aer - aerosol mass species (excluding water) [kmol/kmol]
   !    XXX=wat - aerosol water [kmol/kmol]
   !    XXX=num - aerosol number [#/kmol]
   !
   !    N=1 - before gas-phase chemistry
   !    N=2 - before cloud chemistry
   !    N=3 - current incoming values (before gas-aerosol exchange, newnuc, coag)
   !    N=_cur - updated outgoing values (after  gas-aerosol exchange, newnuc, coag)
   !
   real(r8), intent(in   ), dimension(max_gas) :: qgas1, qgas3  ! [kmol/kmol]
   real(r8), intent(inout), dimension(max_gas) :: qgas_cur      ! [kmol/kmol]

   real(r8), intent(in   ), dimension(max_mode) :: qnum3        ! [#/kmol]
   real(r8), intent(inout), dimension(max_mode) :: qnum_cur     ! [#/kmol]

   real(r8), intent(in   ), dimension(max_aer,max_mode) :: qaer2, qaer3  ! [kmol/kmol]
   real(r8), intent(inout), dimension(max_aer,max_mode) :: qaer_cur      ! [kmol/kmol]

   real(r8), intent(in   ), dimension(max_mode) :: qnumcw3      ! [#/kmol]
   real(r8), intent(inout), dimension(max_mode) :: qnumcw_cur   ! [#/kmol]

   real(r8), intent(in   ), dimension(max_aer,max_mode) :: qaercw2, qaercw3  ! [kmol/kmol]
   real(r8), intent(inout), dimension(max_aer,max_mode) :: qaercw_cur        ! [kmol/kmol]

   real(r8), intent(in   ), dimension(max_mode) :: qwtr3       ! [kmol/kmol]
   real(r8), intent(inout), dimension(max_mode) :: qwtr_cur    ! [kmol/kmol]

   ! qXXX_delaa are TMR changes (increments, not tendencies) of different microphysics processes.
   ! These are diagnostics sent to history output; they do not directly affect time integration.

   real(r8), intent(inout), dimension(max_gas, nqtendaa)             :: qgas_delaa    ! [kmol/kmol]
   real(r8), intent(inout), dimension(max_mode,nqtendaa)             :: qnum_delaa    ! [   #/kmol]
   real(r8), intent(inout), dimension(max_aer,max_mode,nqtendaa)     :: qaer_delaa    ! [kmol/kmol]
   real(r8), intent(inout), dimension(max_mode,nqqcwtendaa)          :: qnumcw_delaa  ! [   #/kmol]
   real(r8), intent(inout), dimension(max_aer,max_mode,nqqcwtendaa ) :: qaercw_delaa  ! [kmol/kmol]

   type ( misc_vars_aa_type ), intent(inout) :: misc_vars_aa_sub

   ! Local variables

   integer,parameter ::  ntsubstep = 1
   integer  :: jtsubstep     ! sub-timestep index
   real(r8) :: dtsubstep     ! length time substep [s]

   ! if dest_mode_of_mode(imode) >  0, then mode imode gets renamed into mode dest_mode_of_mode(imode)
   ! if dest_mode_of_mode(imode) <= 0, then mode imode does not have renaming
   integer :: dest_mode_of_mode(max_mode)

   ! Tmp variables of mixing ratios used for diagnosing increments or for process coupling

   real(r8), dimension( 1:max_gas ) :: qgas_sv1, qgas_avg  ! [kmol/kmol]

   real(r8), dimension( 1:max_mode ) :: qnum_sv1    ! [#/kmol]
   real(r8), dimension( 1:max_mode ) :: qnumcw_sv1  ! [#/kmol]

   real(r8), dimension( 1:max_aer, 1:max_mode ) :: qaer_sv1     ! [kmol/kmol]
   real(r8), dimension( 1:max_aer, 1:max_mode ) :: qaercw_sv1   ! [kmol/kmol]

   ! Mixing ratio increments of sub-timesteps used for process coupling

   real(r8), dimension( 1:max_mode )            :: qnum_delsub_cond, qnum_delsub_coag   ! [   #/kmol]
   real(r8), dimension( 1:max_aer, 1:max_mode ) :: qaer_delsub_cond, qaer_delsub_coag   ! [kmol/kmol]
   real(r8), dimension( 1:max_aer, 1:max_mode ) :: qaer_delsub_grow4rnam                ! [kmol/kmol]
   real(r8), dimension( 1:max_aer, 1:max_mode ) :: qaercw_delsub_grow4rnam              ! [kmol/kmol]
   real(r8), dimension( 1:max_aer, 1:max_agepair ) :: qaer_delsub_coag_in               ! [kmol/kmol]

   real(r8) :: del_h2so4_gasprod   ! [kmol/kmol]
   real(r8) :: del_h2so4_aeruptk   ! [kmol/kmol]

   ! Tendencies and coefficients usef for process coupling

   real(r8) :: qgas_netprod_otrproc(max_gas)
             ! qgas_netprod_otrproc = gas net production rate from other processes
             !    such as gas-phase chemistry and emissions [kmol/kmol/s]
             ! this allows the condensation (gasaerexch) routine to apply production and condensation loss
             !    together, which is more accurate numerically
             ! NOTE - must be >= zero, as numerical method can fail when it is negative
             ! NOTE - currently only the values for h2so4 and nh3 should be non-zero

   real(r8) :: uptkaer(max_gas,max_mode) ! uptake rates of individual gas species and aerosol modes [1/s]
   real(r8) :: uptkrate_h2so4            ! total uptake rate of H2SO4 (summed over all modes) [1/s]
   real(r8) :: dnclusterdt_substep       ! cluster nucleation rate [#/m3/s]

   ! Miscellaneous

   real(r8):: aircon                ! air molar density [kmol/m3]
   integer :: igas                  ! gas index
   logical :: do_aging_in_subarea

   !---------------------------------------------------------------------------------------
   ! Calculate air molar density [kmol/m3] to be passed on to individual parameterizations
   !---------------------------------------------------------------------------------------
   aircon = pmid/(r_universal*temp) 

   !----------------------------------------------------------
   ! Initializ mixing ratios with the before-amicphys values
   !----------------------------------------------------------
   qgas_cur = qgas3
   qaer_cur = qaer3
   qnum_cur = qnum3
   qwtr_cur = qwtr3

   if (iscldy_subarea) then
      qnumcw_cur = qnumcw3
      qaercw_cur = qaercw3
   end if

   !---------------------------------------------------------------------
   ! Diagnose net production rate of H2SO4 gas production
   ! cause by other processes (e.g., gas chemistry and cloud chemistry) 
   !---------------------------------------------------------------------
   qgas_netprod_otrproc(:) = 0.0_r8

   ! If gaexch_h2so4_uptake_optaa == 2, then
   !  - if qgas increases from pre-gaschem to post-cldchem,
   !    start from the pre-gaschem mix-ratio and add in the production during the integration
   !  - if it decreases,  start from post-cldchem mix-ratio
   ! *** currently just do this for h2so4 (and nh3 if considered in model)

   if ( (do_cond_sub) .and. (gaexch_h2so4_uptake_optaa==2) ) then
      do igas = 1, ngas
         if ((igas == igas_h2so4) .or. (igas == igas_nh3)) then
            qgas_netprod_otrproc(igas) = (qgas3(igas) - qgas1(igas))/deltat
            if ( qgas_netprod_otrproc(igas) >= 0.0_r8 ) then
               qgas_cur(igas) = qgas1(igas)
            else
               qgas_netprod_otrproc(igas) = 0.0_r8
            end if
         end if
      end do ! igas
   end if
   del_h2so4_gasprod = max( qgas3(igas_h2so4)-qgas1(igas_h2so4), 0.0_r8 )/ntsubstep

   !-----------------------------------
   ! Initialize increment diagnostics
   !-----------------------------------
   qgas_delaa(:,:)   = 0._r8
   qnum_delaa(:,:)   = 0._r8
   qaer_delaa(:,:,:) = 0._r8

   qnumcw_delaa(:,:)   = 0._r8 
   qaercw_delaa(:,:,:) = 0._r8

   misc_vars_aa_sub%ncluster_tend_nnuc_1grid = 0._r8

   !***********************************
   ! loop over multiple time sub-steps
   !***********************************
   dtsubstep = deltat/ntsubstep

   jtsubstep_loop: do jtsubstep = 1, ntsubstep

      !======================
      ! Gas-aerosol exchange
      !======================
      uptkrate_h2so4 = 0.0_r8

      if ( do_cond_sub ) then

         qgas_sv1 = qgas_cur
         qnum_sv1 = qnum_cur
         qaer_sv1 = qaer_cur

         call mam_gasaerexch_1subarea(                    &
           nstep,             lchnk,                      &
           ii,                kk,               jsubarea, &
           jtsubstep,         ntsubstep,                  &
           latndx,            lonndx,           lund,     &
           dtsubstep,                                     &
           temp,              pmid,             aircon,   &
           n_mode,                                        &
           qgas_cur,          qgas_avg,                   &
           qgas_netprod_otrproc,                          &
           qaer_cur,                                      &
           qnum_cur,                                      &
           qwtr_cur,                                      &
           dgn_a,             dgn_awet,         wetdens,  &
           uptkaer,           uptkrate_h2so4              )

         qgas_delaa(:,iqtend_cond) = qgas_delaa(:,iqtend_cond) + (qgas_cur - (qgas_sv1 + qgas_netprod_otrproc*dtsubstep)) 

         qnum_delsub_cond = qnum_cur - qnum_sv1
         qaer_delsub_cond = qaer_cur - qaer_sv1

         del_h2so4_aeruptk = qgas_cur(igas_h2so4) &
                           - (qgas_sv1(igas_h2so4) + qgas_netprod_otrproc(igas_h2so4)*dtsubstep)

      else ! do_cond_sub

        qgas_avg(1:ngas) = qgas_cur(1:ngas)
        qaer_delsub_cond = 0.0_r8
        qnum_delsub_cond = 0.0_r8
        del_h2so4_aeruptk = 0.0_r8

      end if !do_cond_sub

      !====================================
      ! Renaming after "continuous growth"
      !====================================
      if ( do_rename_sub ) then

         dest_mode_of_mode(:) = 0
         dest_mode_of_mode(nait) = nacc

         !---------------------------------------------------------
         ! Calculate changes in aerosol mass mixing ratios due to 
         !  - gas condensation/evaporation
         !  - cloud chemistry (if the subarea is cloudy)
         !---------------------------------------------------------
         qaer_delsub_grow4rnam  = qaer_delsub_cond

         if (iscldy_subarea) then 
            qaer_delsub_grow4rnam   = (qaer3 - qaer2)/ntsubstep  + qaer_delsub_grow4rnam
            qaercw_delsub_grow4rnam = (qaercw3 - qaercw2)/ntsubstep
         end if

         !----------
         ! Renaming
         !----------
         qnum_sv1 = qnum_cur
         qaer_sv1 = qaer_cur

         qnumcw_sv1 = qnumcw_cur
         qaercw_sv1 = qaercw_cur

         call mam_rename_1subarea(                          &
            iscldy_subarea,                                 &
            dest_mode_of_mode, n_mode,                      &
            qnum_cur,   qaer_cur,    qaer_delsub_grow4rnam, &
            qnumcw_cur, qaercw_cur, qaercw_delsub_grow4rnam )

         !------------------------
         ! Accumulate increments
         !------------------------
         qnum_delaa  (:,iqtend_rnam) = qnum_delaa  (:,iqtend_rnam) + (qnum_cur - qnum_sv1) 
         qaer_delaa(:,:,iqtend_rnam) = qaer_delaa(:,:,iqtend_rnam) + (qaer_cur - qaer_sv1) 

         if (iscldy_subarea) then
            qnumcw_delaa  (:,iqqcwtend_rnam) = qnumcw_delaa  (:,iqqcwtend_rnam) + (qnumcw_cur - qnumcw_sv1) 
            qaercw_delaa(:,:,iqqcwtend_rnam) = qaercw_delaa(:,:,iqqcwtend_rnam) + (qaercw_cur - qaercw_sv1)
         end if

      end if !do_rename_sub

      !====================================
      ! New particle formation (nucleation)
      !====================================
      if ( do_newnuc_sub ) then

         qgas_sv1 = qgas_cur
         qnum_sv1 = qnum_cur
         qaer_sv1 = qaer_cur

         call mam_newnuc_1subarea(                                     &
            nstep,             lchnk,                                  &
            ii,                kk,               jsubarea,             &
            latndx,            lonndx,           lund,                 &
            dtsubstep,                                                 &
            temp,              pmid,             aircon,               &
            zmid,              pblh,             relhum,               &
            uptkrate_h2so4,   del_h2so4_gasprod, del_h2so4_aeruptk,    &
            n_mode,                                                    &
            qgas_cur,          qgas_avg,                               &
            qnum_cur,                                                  &
            qaer_cur,                                                  &
            qwtr_cur,                                                  &
            dnclusterdt_substep                                        )

         qgas_delaa(:,iqtend_nnuc) = qgas_delaa(:,iqtend_nnuc) + (qgas_cur - qgas_sv1)
         qnum_delaa(:,iqtend_nnuc) = qnum_delaa(:,iqtend_nnuc) + (qnum_cur - qnum_sv1)
         qaer_delaa(:,:,iqtend_nnuc) = qaer_delaa(:,:,iqtend_nnuc) + (qaer_cur - qaer_sv1)

         misc_vars_aa_sub%ncluster_tend_nnuc_1grid = &
         misc_vars_aa_sub%ncluster_tend_nnuc_1grid + dnclusterdt_substep*(dtsubstep/deltat) 

      end if

      !====================================
      ! Coagulation
      !====================================
      if ( do_coag_sub ) then

         qnum_sv1 = qnum_cur
         qaer_sv1 = qaer_cur

         call mam_coag_1subarea( &
            lchnk,      ii,      kk,        &
            dtsubstep,                                &! in
            temp,      pmid,     aircon,              &! in
            dgn_a,     dgn_awet, wetdens,             &! in
            qnum_cur,  qaer_cur, qaer_delsub_coag_in  )! inout, inout, out

         qnum_delsub_coag = qnum_cur - qnum_sv1
         qaer_delsub_coag = qaer_cur - qaer_sv1

         qnum_delaa  (:,iqtend_coag) = qnum_delaa  (:,iqtend_coag) + qnum_delsub_coag 
         qaer_delaa(:,:,iqtend_coag) = qaer_delaa(:,:,iqtend_coag) + qaer_delsub_coag 

      else
         qaer_delsub_coag_in = 0.0_r8
         qaer_delsub_coag = 0.0_r8
         qnum_delsub_coag = 0.0_r8
      end if

      !====================================
      ! primary carbon aging
      !====================================
      do_aging_in_subarea = ( n_agepair>0 ) .and. &
                            ( (.not.iscldy_subarea).or.(iscldy_subarea.and.do_cond_sub) )

      if (do_aging_in_subarea) then

         call mam_pcarbon_aging_1subarea(                          &
            dgn_a,             n_mode,                             &! input
            qnum_cur,          qnum_delsub_cond, qnum_delsub_coag, &! in-outs
            qaer_cur,          qaer_delsub_cond, qaer_delsub_coag, &! in-outs
            qaer_delsub_coag_in)                                    ! in-outs

      end if

      ! The following block has to be placed here (after both condensation and aging)
      ! as both can change the values of qnum_delsub_cond and qaer_delsub_cond.

      if ( do_cond_sub ) then
         qnum_delaa  (:,iqtend_cond) = qnum_delaa  (:,iqtend_cond) + qnum_delsub_cond 
         qaer_delaa(:,:,iqtend_cond) = qaer_delaa(:,:,iqtend_cond) + qaer_delsub_cond 
      end if

   end do jtsubstep_loop
   !***********************************

end subroutine mam_amicphys_1subarea
!--------------------------------------------------------------------------------
