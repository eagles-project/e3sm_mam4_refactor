!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
      subroutine mam_amicphys_1subarea(        &
         do_cond,                do_rename,          &
         do_newnuc,              do_coag,            &
         nstep,      lchnk,      i,        k,        &
         latndx,     lonndx,     lund,               &
         loffset,    deltat,                         &
         jsub,                   nsubarea,           &
         iscldy_subarea,         afracsub,           &
         temp,       pmid,       pdel,               &
         zmid,       pblh,       relhum,             &
         dgn_a,      dgn_awet,   wetdens,            &
         qgas1,      qgas3,      qgas4,              &
         qgas_delaa,                                 &
                qnum3,      qnum4,   qnum_delaa,     &
         qaer2, qaer3,      qaer4,   qaer_delaa,     &
         qwtr3,      qwtr4,                          &
         qnumcw3,    qnumcw4,                        &
         qnumcw_delaa,                               &
         qaercw2,    qaercw3,    qaercw4,            &
         qaercw_delaa,                               &
         misc_vars_aa_sub                            )

  use modal_aero_amicphys_control
!
! CLEAR SUBAREA ----
! calculates changes to gas and aerosol sub-area TMRs (tracer mixing ratios)
!    for a single clear sub-area (with indices = lchnk,i,k,jsub)
! qgas3, qaer3, qnum3 are the current incoming TMRs
! qgas4, qaer4, qnum4 are the updated outgoing TMRs
!
! this routine calculates changes involving
!    gas-aerosol exchange (condensation/evaporation)
!    growth from smaller to larger modes (renaming) due to condensation
!    new particle nucleation
!    coagulation
!    transfer of particles from hydrophobic modes to hydrophilic modes (aging)
!       due to condensation and coagulation
!
! CLOUDY SUBAREA
!
! calculates changes to gas and aerosol sub-area TMRs (tracer mixing ratios)
!    for a single cloudy sub-area (with indices = lchnk,i,k,jsub)
! qgas3, qaer3, qaercw3, qnum3, qnumcw3 are the current incoming TMRs
! qgas4, qaer4, qaercw4, qnum4, qnumcw4 are the updated outgoing TMRs
!
! when do_cond = false, this routine only calculates changes involving
!    growth from smaller to larger modes (renaming) following cloud chemistry
!    so gas TMRs are not changed
! when do_cond = true, this routine also calculates changes involving
!    gas-aerosol exchange (condensation/evaporation)
!    transfer of particles from hydrophobic modes to hydrophilic modes (aging)
!       due to condensation
! currently this routine does not do
!    new particle nucleation - because h2so4 gas conc. should be very low in cloudy air
!    coagulation - because cloud-borne aerosol would need to be included
!
      use physconst, only:  r_universal
      use modal_aero_coag, only: mam_coag_1subarea
      use modal_aero_rename, only: mam_rename_1subarea

      logical,  intent(in)    :: do_cond, do_rename, do_newnuc, do_coag
      logical,  intent(in)    :: iscldy_subarea        ! true if sub-area is cloudy
      integer,  intent(in)    :: lchnk                 ! chunk identifier
      integer,  intent(in)    :: nstep                 ! model time-step number
      integer,  intent(in)    :: i, k                  ! column and level indices
      integer,  intent(in)    :: latndx, lonndx        ! lat and lon indices
      integer,  intent(in)    :: lund                  ! logical unit for diagnostic output
      integer,  intent(in)    :: loffset
      integer,  intent(in)    :: jsub, nsubarea        ! sub-area index, number of sub-areas

      real(r8), intent(in)    :: afracsub              ! fractional area of sub-area (0-1)
      real(r8), intent(in)    :: deltat                ! time step (s)

      real(r8), intent(in)    :: temp                  ! temperature at model levels (K)
      real(r8), intent(in)    :: pmid                  ! pressure at layer center (Pa)
      real(r8), intent(in)    :: pdel                  ! pressure thickness of layer (Pa)
      real(r8), intent(in)    :: zmid                  ! altitude (above ground) at layer center (m)
      real(r8), intent(in)    :: pblh                  ! planetary boundary layer depth (m)
      real(r8), intent(in)    :: relhum                ! relative humidity (0-1)

      real(r8), intent(inout) :: dgn_a(max_mode)
      real(r8), intent(inout) :: dgn_awet(max_mode)
                                    ! dry & wet geo. mean dia. (m) of number distrib.
      real(r8), intent(inout) :: wetdens(max_mode)
                                    ! interstitial aerosol wet density (kg/m3)
                                    ! dry & wet geo. mean dia. (m) of number distrib.

! qXXXN (X=gas,aer,wat,num; N=1:4) are sub-area mixing ratios
!    XXX=gas - gas species
!    XXX=aer - aerosol mass  species (excluding water)
!    XXX=wat - aerosol water
!    XXX=num - aerosol number
!    N=1 - before gas-phase chemistry
!    N=2 - before cloud chemistry
!    N=3 - current incoming values (before gas-aerosol exchange, newnuc, coag)
!    N=4 - updated outgoing values (after  gas-aerosol exchange, newnuc, coag)
!
! qXXX_delaa are TMR changes (not tendencies)
!    for different processes, which are used to produce history output
! for a clear sub-area, the processes are condensation/evaporation (and associated aging),
!    renaming, coagulation, and nucleation
      real(r8), intent(in   ), dimension( 1:max_gas ) :: qgas1, qgas3
      real(r8), intent(inout), dimension( 1:max_gas ) :: qgas4
      real(r8), intent(inout), dimension( 1:max_gas, 1:nqtendaa ) :: qgas_delaa

      real(r8), intent(in   ), dimension( 1:max_mode ) ::  qnum3
      real(r8), intent(inout), dimension( 1:max_mode ) ::  qnum4
      real(r8), intent(inout), dimension( 1:max_mode, 1:nqtendaa ) :: qnum_delaa

      real(r8), intent(in   ), dimension( 1:max_aer, 1:max_mode ) :: qaer2, qaer3
      real(r8), intent(inout), dimension( 1:max_aer, 1:max_mode ) ::  qaer4
      real(r8), intent(inout), dimension( 1:max_aer, 1:max_mode, 1:nqtendaa ) :: qaer_delaa

      real(r8), intent(in   ), dimension( 1:max_mode ) :: qwtr3
      real(r8), intent(inout), dimension( 1:max_mode ) :: qwtr4

      real(r8), intent(in   ), dimension( 1:max_mode ) :: qnumcw3
      real(r8), intent(inout), dimension( 1:max_mode ) :: qnumcw4
      real(r8), intent(inout), dimension( 1:max_mode, 1:nqqcwtendaa ) :: qnumcw_delaa

      real(r8), intent(in   ), dimension( 1:max_aer, 1:max_mode ) :: qaercw2, qaercw3
      real(r8), intent(inout), dimension( 1:max_aer, 1:max_mode ) :: qaercw4
      real(r8), intent(inout), dimension( 1:max_aer, 1:max_mode, 1:nqqcwtendaa ) :: qaercw_delaa

      type ( misc_vars_aa_type ), intent(inout) :: misc_vars_aa_sub

! local
      integer, parameter :: ntot_poaspec = npoa
      integer, parameter :: ntot_soaspec = nsoa

      integer :: iaer, igas, ip
      integer :: jtsubstep
      integer :: ll
! if dest_mode_of_mode(n) >  0, then mode n gets renamed into mode dest_mode_of_mode(n)
! if dest_mode_of_mode(n) <= 0, then mode n does not have renaming
      integer :: dest_mode_of_mode(max_mode)
      integer :: n, ntsubstep
      integer :: n_mode
      integer :: ntot_soamode

      logical, parameter :: flag_pcarbon_opoa_frac_zero   = .true.
      logical, parameter :: flag_nh4_lt_2so4_each_step  = .false.

      logical :: skip_soamode(max_mode)   ! true if this mode does not have soa

      real(r8), dimension( 1:max_gas ) :: qgas_cur, qgas_sv1, qgas_avg
      real(r8), dimension( 1:max_gas ) :: qgas_del_cond, qgas_del_nnuc, qgas_netprod_otrproc
                ! qgas_netprod_otrproc = gas net production rate from other processes
                !    such as gas-phase chemistry and emissions (mol/mol/s)
                ! this allows the condensation (gasaerexch) routine to apply production and condensation loss
                !    together, which is more accurate numerically
                ! NOTE - must be >= zero, as numerical method can fail when it is negative
                ! NOTE - currently only the values for h2so4 and nh3 should be non-zero

! qxxx_del_yyyy    are mix-ratio changes over full time step (deltat)
! qxxx_delsub_yyyy are mix-ratio changes over time sub-step (dtsubstep)
      real(r8), dimension( 1:max_mode ) :: qnum_cur, qnum_sv1
      real(r8), dimension( 1:max_mode ) :: qnum_del_cond, qnum_del_rnam, qnum_del_nnuc, qnum_del_coag, &
         qnum_delsub_cond, qnum_delsub_coag

      real(r8), dimension( 1:max_mode ) :: qnumcw_cur, qnumcw_sv1
      real(r8), dimension( 1:max_mode ) :: qnumcw_del_rnam

      real(r8), dimension( 1:max_aer, 1:max_mode ) :: qaer_cur, qaer_sv1
      real(r8), dimension( 1:max_aer, 1:max_agepair ) :: qaer_delsub_coag_in
      real(r8), dimension( 1:max_aer, 1:max_mode ) :: qaer_del_cond, qaer_del_rnam, qaer_del_nnuc, qaer_del_coag, &
         qaer_delsub_grow4rnam, &
         qaer_delsub_cond, qaer_delsub_coag

      real(r8), dimension( 1:max_aer, 1:max_mode ) :: qaercw_cur, qaercw_sv1
      real(r8), dimension( 1:max_aer, 1:max_mode ) :: qaercw_del_rnam, qaercw_delsub_grow4rnam

      real(r8), dimension( 1:max_mode ) :: qwtr_cur

      real(r8) :: aircon                           ! air molar density (kmol/m3)
      real(r8) :: del_h2so4_gasprod
      real(r8) :: del_h2so4_aeruptk
      real(r8) :: dnclusterdt, dnclusterdt_substep
      real(r8) :: dtsubstep                        ! time sub-step
      real(r8) :: gas_diffus(max_gas)              ! gas diffusivity at current temp and pres (m2/s)
      real(r8) :: gas_freepath(max_gas)            ! gas mean free path at current temp and pres (m)

      real(r8) :: tmpa, tmpb, tmpc, tmpd, tmpe, tmpf
      real(r8) :: tmp_relhum
      real(r8) :: uptkaer(max_gas,max_mode)
      real(r8) :: uptkrate_h2so4

      logical :: do_aging_in_subarea

      aircon = pmid/(r_universal*temp) ! air molar density (kmol/m3)
      n_mode = ntot_amode

      qgas_cur = qgas3
      qaer_cur = qaer3
      qnum_cur = qnum3
      qwtr_cur = qwtr3

      if (iscldy_subarea) then
         qnumcw_cur = qnumcw3
         qaercw_cur = qaercw3
      end if

      qgas_netprod_otrproc(:) = 0.0_r8
      if ( ( do_cond                         ) .and. &
           ( gaexch_h2so4_uptake_optaa == 2 ) ) then
         do igas = 1, ngas
            if ((igas == igas_h2so4) .or. (igas == igas_nh3)) then
! if gaexch_h2so4_uptake_optaa == 2, then
!    if qgas increases from pre-gaschem to post-cldchem,
!       start from the pre-gaschem mix-ratio and add in the production
!       during the integration
!    if it decreases,
!       start from post-cldchem mix-ratio
! *** currently just do this for h2so4 and nh3
               qgas_netprod_otrproc(igas) = (qgas3(igas) - qgas1(igas))/deltat
               if ( qgas_netprod_otrproc(igas) >= 0.0_r8 ) then
                  qgas_cur(igas) = qgas1(igas)
               else
                  qgas_netprod_otrproc(igas) = 0.0_r8
               end if
            end if
         end do ! igas
      end if

      qgas_del_cond = 0.0_r8
      qgas_del_nnuc = 0.0_r8

      qaer_del_cond = 0.0_r8
      qaer_del_rnam = 0.0_r8
      qaer_del_nnuc = 0.0_r8
      qaer_del_coag = 0.0_r8
      qaer_delsub_cond = 0.0_r8

     if (iscldy_subarea) then
      qaercw_del_rnam = 0.0_r8
     else ! clear
      qaer_delsub_coag_in = 0.0_r8
      qaer_delsub_coag = 0.0_r8
     end if

      qnum_del_cond = 0.0_r8
      qnum_del_rnam = 0.0_r8
      qnum_del_nnuc = 0.0_r8
      qnum_del_coag = 0.0_r8
      qnum_delsub_cond = 0.0_r8

     if (iscldy_subarea) then
      qnumcw_del_rnam = 0.0_r8
     else
      qnum_delsub_coag = 0.0_r8
     end if

      dnclusterdt = 0.0_r8

      ntsubstep = 1
      dtsubstep = deltat
      if (ntsubstep > 1) dtsubstep = deltat/ntsubstep

      del_h2so4_gasprod = max( qgas3(igas_h2so4)-qgas1(igas_h2so4), 0.0_r8 )/ntsubstep

!
!
! loop over multiple time sub-steps
!
!
jtsubstep_loop: &
      do jtsubstep = 1, ntsubstep


!
!
! gas-aerosol exchange
!
!
      uptkrate_h2so4 = 0.0_r8

do_cond_if_block10: &
      if ( do_cond ) then

      qgas_sv1 = qgas_cur
      qnum_sv1 = qnum_cur
      qaer_sv1 = qaer_cur

         call mam_gasaerexch_1subarea(                                &
           nstep,             lchnk,                                  &
           i,                 k,                jsub,                 &
           jtsubstep,         ntsubstep,                              &
           latndx,            lonndx,           lund,                 &
           dtsubstep,                                                 &
           temp,              pmid,             aircon,               &
           n_mode,                                                    &
           qgas_cur,          qgas_avg,                               &
           qgas_netprod_otrproc,                                      &
           qaer_cur,                                                  &
           qnum_cur,                                                  &
           qwtr_cur,                                                  &
           dgn_a,             dgn_awet,         wetdens,              &
           uptkaer,           uptkrate_h2so4                          )


      if (newnuc_h2so4_conc_optaa == 11) then
         qgas_avg(igas_h2so4) = 0.5_r8*(qgas_sv1(igas_h2so4) + qgas_cur(igas_h2so4))
      else if (newnuc_h2so4_conc_optaa == 12) then
         qgas_avg(igas_h2so4) = qgas_cur(igas_h2so4)
      end if

      qgas_del_cond = qgas_del_cond + (qgas_cur - (qgas_sv1 + qgas_netprod_otrproc*dtsubstep))
      qnum_delsub_cond = qnum_cur - qnum_sv1
      qaer_delsub_cond = qaer_cur - qaer_sv1
! qaer_delsub_grow4rnam = change in qaer_del_cond during latest condensation calculations
      qaer_delsub_grow4rnam = qaer_cur - qaer_sv1

      del_h2so4_aeruptk = qgas_cur(igas_h2so4) &
                       - (qgas_sv1(igas_h2so4) + qgas_netprod_otrproc(igas_h2so4)*dtsubstep)

      else ! do_cond_if_block10

      qgas_avg(1:ngas) = qgas_cur(1:ngas)
      qaer_delsub_grow4rnam(:,:) = 0.0_r8

      del_h2so4_aeruptk = 0.0_r8

      end if do_cond_if_block10


!
!
! renaming after "continuous growth"
!
!
do_rename_if_block30: &
      if ( do_rename ) then

      dest_mode_of_mode(:) = 0
      dest_mode_of_mode(nait) = nacc

      if (iscldy_subarea) then 
! qaer_delsub_grow4rnam   = change in qaer from cloud chemistry and gas condensation
! qaercw_delsub_grow4rnam = change in qaercw from cloud chemistry
      qaer_delsub_grow4rnam   = (qaer3 - qaer2)/ntsubstep  + qaer_delsub_grow4rnam
      qaercw_delsub_grow4rnam = (qaercw3 - qaercw2)/ntsubstep
      end if

      qnum_sv1 = qnum_cur
      qaer_sv1 = qaer_cur

      if (iscldy_subarea) then 
      qnumcw_sv1 = qnumcw_cur
      qaercw_sv1 = qaercw_cur
      end if

      if (.not.iscldy_subarea) then 
      call mam_rename_1subarea(                                    &
         iscldy_subarea,                                           &
         dest_mode_of_mode,                                        &
         n_mode,                                                   &
         qnum_cur,                                                 &
         qaer_cur,          qaer_delsub_grow4rnam                  )

      else
      call mam_rename_1subarea(                                      &
         iscldy_subarea,                                             &
         dest_mode_of_mode,                                          &
         n_mode,                                                     &
         qnum_cur,                                                   &
         qaer_cur,          qaer_delsub_grow4rnam,                   &
         qnumcw_cur,                                                 &
         qaercw_cur,        qaercw_delsub_grow4rnam                  )

      end if

      qnum_del_rnam = qnum_del_rnam + (qnum_cur - qnum_sv1)
      qaer_del_rnam = qaer_del_rnam + (qaer_cur - qaer_sv1)

      if (iscldy_subarea) then
      qnumcw_del_rnam = qnumcw_del_rnam + (qnumcw_cur - qnumcw_sv1)
      qaercw_del_rnam = qaercw_del_rnam + (qaercw_cur - qaercw_sv1)
      end if

      end if do_rename_if_block30



   if (.not.iscldy_subarea) then

!
!
! new particle formation (nucleation)
!
!
      do_newnuc_if_block50: &
      if ( do_newnuc ) then

      qgas_sv1 = qgas_cur
      qnum_sv1 = qnum_cur
      qaer_sv1 = qaer_cur

      call mam_newnuc_1subarea(                                     &
         nstep,             lchnk,                                  &
         i,                 k,                jsub,                 &
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

      qgas_del_nnuc = qgas_del_nnuc + (qgas_cur - qgas_sv1)
      qnum_del_nnuc = qnum_del_nnuc + (qnum_cur - qnum_sv1)
      qaer_del_nnuc = qaer_del_nnuc + (qaer_cur - qaer_sv1)
      dnclusterdt = dnclusterdt + dnclusterdt_substep*(dtsubstep/deltat)

      end if do_newnuc_if_block50


!
!
! coagulation part
!
!
      if ( do_coag ) then

      qnum_sv1 = qnum_cur
      qaer_sv1 = qaer_cur

      call mam_coag_1subarea(                                       &
         dtsubstep,                                                 &! in
         temp,              pmid,             aircon,               &! in
         dgn_a,             dgn_awet,         wetdens,              &! in
         qnum_cur,          qaer_cur,         qaer_delsub_coag_in   )! inout, inout, out

      qnum_delsub_coag = qnum_cur - qnum_sv1
      qaer_delsub_coag = qaer_cur - qaer_sv1

      end if

   end if ! .not.iscldy_subarea

!
!
! primary carbon aging
!
!

      do_aging_in_subarea = ( n_agepair>0 ) .and. &
                            ( (.not.iscldy_subarea).or.(iscldy_subarea.and.do_cond) )

      if (do_aging_in_subarea) then

         if (iscldy_subarea) then
            qaer_delsub_coag_in = 0.0_r8
            qaer_delsub_coag = 0.0_r8
            qnum_delsub_coag = 0.0_r8
         end if

         call mam_pcarbon_aging_1subarea(                              &
            dgn_a,             n_mode,                                 &  ! input
            qnum_cur,          qnum_delsub_cond, qnum_delsub_coag,     &  ! in-outs
            qaer_cur,          qaer_delsub_cond, qaer_delsub_coag,     &  ! in-outs
            qaer_delsub_coag_in)                                          ! in-outs

      end if

! accumulate sub-step q-dels

      if ( do_coag .and. (.not.iscldy_subarea) ) then
         qnum_del_coag = qnum_del_coag + qnum_delsub_coag
         qaer_del_coag = qaer_del_coag + qaer_delsub_coag
      end if

      if ( do_cond ) then
         qnum_del_cond = qnum_del_cond + qnum_delsub_cond
         qaer_del_cond = qaer_del_cond + qaer_delsub_cond
      end if

end do jtsubstep_loop


!
!
! final mix ratios
!
!
      qgas4 = qgas_cur
      qaer4 = qaer_cur
      qnum4 = qnum_cur
      qwtr4 = qwtr_cur

     if (iscldy_subarea) then
      qnumcw4 = qnumcw_cur
      qaercw4 = qaercw_cur
     end if

! final mix ratio changes

      qgas_delaa(:,iqtend_cond) = qgas_del_cond(:)
      qgas_delaa(:,iqtend_rnam) = 0.0_r8
      qgas_delaa(:,iqtend_coag) = 0.0_r8

     if (iscldy_subarea) then
        qgas_delaa(:,iqtend_nnuc) = 0.0_r8
        qnum_delaa(:,iqtend_nnuc) = 0.0_r8
        qnum_delaa(:,iqtend_coag) = 0.0_r8
        qaer_delaa(:,:,iqtend_nnuc) = 0.0_r8
        qaer_delaa(:,:,iqtend_coag) = 0.0_r8
      qnumcw_delaa(:,iqqcwtend_rnam) = qnumcw_del_rnam(:)
      qaercw_delaa(:,:,iqqcwtend_rnam) = qaercw_del_rnam(:,:)
     else
        qgas_delaa(:,iqtend_nnuc) = qgas_del_nnuc(:)
        qnum_delaa(:,iqtend_nnuc) = qnum_del_nnuc(:)
        qnum_delaa(:,iqtend_coag) = qnum_del_coag(:)
        qaer_delaa(:,:,iqtend_nnuc) = qaer_del_nnuc(:,:)
        qaer_delaa(:,:,iqtend_coag) = qaer_del_coag(:,:)
     end if

      qnum_delaa(:,iqtend_cond) = qnum_del_cond(:)
      qnum_delaa(:,iqtend_rnam) = qnum_del_rnam(:)

      qaer_delaa(:,:,iqtend_cond) = qaer_del_cond(:,:)
      qaer_delaa(:,:,iqtend_rnam) = qaer_del_rnam(:,:)

      misc_vars_aa_sub%ncluster_tend_nnuc_1grid = dnclusterdt

      return
      end subroutine mam_amicphys_1subarea

