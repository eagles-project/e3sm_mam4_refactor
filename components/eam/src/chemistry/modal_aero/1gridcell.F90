
!--------------------------------------------------------------------------------
subroutine mam_amicphys_1gridcell(          &
         do_cond,            do_rename,           &
         do_newnuc,          do_coag,             &
         nstep,    lchnk,    i,        k,         &
         latndx,   lonndx,   lund,                &
         loffset,  deltat,                        &
         nsubarea,  ncldy_subarea,                &
         iscldy_subarea,     afracsub,            &
         temp,     pmid,     pdel,                &
         zmid,     pblh,     relhumsub,           &
         dgn_a,    dgn_awet, wetdens,             &
         qsub1,                                   &
         qsub2, qqcwsub2,                         &
         qsub3, qqcwsub3, qaerwatsub3,            &
         qsub4, qqcwsub4, qaerwatsub4,            &
         qsub_tendaa, qqcwsub_tendaa,             &
         misc_vars_aa                             )
!
! calculates changes to gas and aerosol sub-area TMRs (tracer mixing ratios)
!    for the current grid cell (with indices = lchnk,i,k)
! qsub3 and qqcwsub3 are the incoming current TMRs
! qsub4 and qqcwsub4 are the outgoing updated TMRs
!
  use modal_aero_amicphys_control

      logical,  intent(in)    :: do_cond, do_rename, do_newnuc, do_coag
      logical,  intent(in)    :: iscldy_subarea(maxsubarea)

      integer,  intent(in)    :: nstep                 ! model time-step number
      integer,  intent(in)    :: lchnk                 ! chunk identifier
      integer,  intent(in)    :: i, k                  ! column and level indices
      integer,  intent(in)    :: latndx, lonndx        ! lat and lon indices
      integer,  intent(in)    :: lund                  ! logical unit for diagnostic output
      integer,  intent(in)    :: loffset
      integer,  intent(in)    :: nsubarea, ncldy_subarea

      real(r8), intent(in)    :: deltat                ! time step (s)
      real(r8), intent(in)    :: afracsub(maxsubarea)   ! sub-area fractional area (0-1)

      real(r8), intent(in)    :: temp                  ! temperature at model levels (K)
      real(r8), intent(in)    :: pmid                  ! pressure at layer center (Pa)
      real(r8), intent(in)    :: pdel                  ! pressure thickness of layer (Pa)
      real(r8), intent(in)    :: zmid                  ! altitude (above ground) at layer center (m)
      real(r8), intent(in)    :: pblh                  ! planetary boundary layer depth (m)
      real(r8), intent(in)    :: relhumsub(maxsubarea) ! sub-area relative humidity (0-1)
      real(r8), intent(inout) :: dgn_a(max_mode)       ! dry geo. mean dia. (m) of number distrib.
      real(r8), intent(inout) :: dgn_awet(max_mode)    ! wet geo. mean dia. (m) of number distrib.
      real(r8), intent(inout) :: wetdens(max_mode)     ! interstitial aerosol wet density (kg/m3)

! qsubN and qqcwsubN (N=1:4) are tracer mixing ratios (TMRs, mol/mol or #/kmol) in sub-areas
!    currently there are just clear and cloudy sub-areas
!    the N=1:4 have same meanings as for qgcmN
!    N=1 - before gas-phase chemistry
!    N=2 - before cloud chemistry
!    N=3 - incoming values (before gas-aerosol exchange, newnuc, coag)
!    N=4 - outgoing values (after  gas-aerosol exchange, newnuc, coag)
      real(r8), intent(in   ), dimension( 1:gas_pcnst, 1:maxsubarea ) :: qsub1, qsub2, qsub3, qqcwsub2, qqcwsub3
      real(r8), intent(inout), dimension( 1:gas_pcnst, 1:maxsubarea ) :: qsub4, qqcwsub4
      real(r8), intent(inout), dimension( 1:ntot_amode_extd, 1:maxsubarea ) :: qaerwatsub3, qaerwatsub4   ! aerosol water mixing ratios (mol/mol)

! qsub_tendaa and qqcwsub_tendaa are TMR tendencies
!    for different processes, which are used to produce history output
! the processes are condensation/evaporation (and associated aging),
!    renaming, coagulation, and nucleation
      real(r8), intent(inout), dimension( 1:gas_pcnst, 1:nqtendaa, 1:maxsubarea ) ::  qsub_tendaa
      real(r8), intent(inout), dimension( 1:gas_pcnst, 1:nqqcwtendaa, 1:maxsubarea ) ::  qqcwsub_tendaa
      type ( misc_vars_aa_type ), intent(inout) :: misc_vars_aa

! local
      integer :: iaer, igas
      integer :: jsub
      integer :: icnst
      integer :: imode
      logical :: do_cond_sub, do_rename_sub, do_newnuc_sub, do_coag_sub
      logical :: do_map_gas_sub

      real(r8), dimension( 1:max_gas ) :: qgas1, qgas3, qgas4
      real(r8), dimension( 1:max_mode ) :: qnum3, qnum4, qnumcw3, qnumcw4
      real(r8), dimension( 1:max_aer, 1:max_mode ) :: qaer2, qaer3, qaer4, qaercw2, qaercw3, qaercw4
      real(r8), dimension( 1:max_mode ) :: qwtr3, qwtr4

      real(r8), dimension( 1:max_gas, 1:nqtendaa ) ::  qgas_delaa
      real(r8), dimension( 1:max_mode, 1:nqtendaa ) :: qnum_delaa
      real(r8), dimension( 1:max_mode, 1:nqqcwtendaa ) :: qnumcw_delaa
      real(r8), dimension( 1:max_aer, 1:max_mode, 1:nqtendaa ) ::  qaer_delaa
      real(r8), dimension( 1:max_aer, 1:max_mode, 1:nqqcwtendaa ) :: qaercw_delaa

    ! real(r8) :: tmpa, tmpb, tmpc, tmpd, tmpe, tmpf, tmpn

      type ( misc_vars_aa_type ), dimension(nsubarea) :: misc_vars_aa_sub


   ! the q--4 values will be equal to q--3 values unless they get changed
   qsub4(:,1:nsubarea) = qsub3(:,1:nsubarea)
   qqcwsub4(:,1:nsubarea) = qqcwsub3(:,1:nsubarea)
   qaerwatsub4(:,1:nsubarea) = qaerwatsub3(:,1:nsubarea)

   qsub_tendaa(:,:,1:nsubarea) = 0.0_r8
   qqcwsub_tendaa(:,:,1:nsubarea) = 0.0_r8

main_jsubarea_loop: &
   do jsub = 1, nsubarea

      if ( iscldy_subarea(jsub) .eqv. .true. ) then
         do_cond_sub   = do_cond
         do_rename_sub = do_rename
         do_newnuc_sub = .false.
         do_coag_sub   = .false.
         if (mdo_gaexch_cldy_subarea <= 0) do_cond_sub = .false.
      else
         do_cond_sub   = do_cond
         do_rename_sub = do_rename
         do_newnuc_sub = do_newnuc
         do_coag_sub   = do_coag
      end if

      ! for cldy subarea, only do gases if doing gaexch
      do_map_gas_sub = do_cond_sub .or. do_newnuc_sub

! map incoming sub-area mix-ratios to gas/aer/num arrays

      !--------------------
      ! Gases
      !--------------------
      qgas1(:) = 0.0_r8
      qgas3(:) = 0.0_r8
      if ( do_map_gas_sub .eqv. .true. ) then
         do igas = 1, ngas
          icnst = lmap_gas(igas)
          qgas1(igas) = qsub1(icnst,jsub)*fcvt_gas(igas)
          qgas3(igas) = qsub3(icnst,jsub)*fcvt_gas(igas)
         end do
      end if

      !-----------------------------------------
      ! interstitial aerosol mass and number
      !-----------------------------------------
      qnum3(:)   = 0.0_r8
      do imode = 1, ntot_amode
         icnst = lmap_num(imode)
         qnum3(imode) = qsub3(icnst,jsub)*fcvt_num
      end do ! imode

      qaer2(:,:) = 0.0_r8
      qaer3(:,:) = 0.0_r8
      do imode = 1, ntot_amode
         do iaer = 1, naer
            icnst = lmap_aer(iaer,imode)
            if (icnst > 0) then
               qaer2(iaer,imode) = qsub2(icnst,jsub)*fcvt_aer(iaer)
               qaer3(iaer,imode) = qsub3(icnst,jsub)*fcvt_aer(iaer)
            end if
         end do
      end do ! imode

      !-----------------------------------------
      ! aerosol water
      !-----------------------------------------
      qwtr3(:)   = 0.0_r8
      do imode = 1, ntot_amode
         qwtr3(imode) = qaerwatsub3(imode,jsub)*fcvt_wtr
      end do ! n

      !-----------------------------------------
      ! cloud-borne aerosol mass and number
      !-----------------------------------------
      ! only do cloud-borne for cloudy

      if ( iscldy_subarea(jsub) .eqv. .true. ) then

      qnumcw3(:)   = 0.0_r8
      do imode = 1, ntot_amode
         icnst = lmap_numcw(imode)
         qnumcw3(imode) = qqcwsub3(icnst,jsub)*fcvt_num
      end do ! imode

      qaercw2(:,:) = 0.0_r8
      qaercw3(:,:) = 0.0_r8
      do imode = 1, ntot_amode
         do iaer = 1, naer
            icnst = lmap_aercw(iaer,imode)
            if (icnst > 0) then
               qaercw2(iaer,imode) = qqcwsub2(icnst,jsub)*fcvt_aer(iaer)
               qaercw3(iaer,imode) = qqcwsub3(icnst,jsub)*fcvt_aer(iaer)
            end if
         end do
      end do ! imode

      end if

      !----------------------------------------------------
      ! Calculate aerosol microphysics for a single subarea
      !----------------------------------------------------
      call mam_amicphys_1subarea(                    &
         do_cond_sub,            do_rename_sub,      &
         do_newnuc_sub,          do_coag_sub,        &
         nstep,      lchnk,      i,        k,        &
         latndx,     lonndx,     lund,               &
         loffset,    deltat,                         &
         jsub,                   nsubarea,           &
         iscldy_subarea(jsub),   afracsub(jsub),     &
         temp,       pmid,       pdel,               &
         zmid,       pblh,       relhumsub(jsub),    &
         dgn_a,      dgn_awet,   wetdens,            &
         qgas1,      qgas3,      qgas4,              &
         qgas_delaa,                                 &
         qnum3,      qnum4,                          &
         qnum_delaa,                                 &
         qaer2,      qaer3,      qaer4,              &
         qaer_delaa,                                 &
         qwtr3,      qwtr4,                          &
         qnumcw3,    qnumcw4,                        &
         qnumcw_delaa,                               &
         qaercw2,    qaercw3,    qaercw4,            &
         qaercw_delaa,                               &
         misc_vars_aa_sub(jsub)                      )

      !----------------------------------------------------
      if ((nsubarea == 1) .or. (iscldy_subarea(jsub) .eqv. .false.)) then
         misc_vars_aa%ncluster_tend_nnuc_1grid = misc_vars_aa%ncluster_tend_nnuc_1grid &
                                               + misc_vars_aa_sub(jsub)%ncluster_tend_nnuc_1grid*afracsub(jsub)
      end if

     ! map gas/aer/num arrays (mix-ratio and del=change) back to sub-area arrays

      !----------------------------------------------------
      ! Gases
      !----------------------------------------------------
      if ( do_map_gas_sub .eqv. .true. ) then
       do igas = 1, ngas
          icnst = lmap_gas(igas)
          qsub4(icnst,jsub ) = qgas4(igas)/fcvt_gas(igas)
          qsub_tendaa(icnst,:,jsub) = qgas_delaa(igas,:)/(fcvt_gas(igas)*deltat)
       end do
      end if

      !----------------------------------------------------
      ! Aerosol water
      !----------------------------------------------------
      do imode = 1, ntot_amode
         qaerwatsub4(imode,jsub) = qwtr4(imode)/fcvt_wtr
      end do ! imode

      !----------------------------------------------------
      ! Interstitial aerosol number
      !----------------------------------------------------
      do imode = 1, ntot_amode
         icnst = lmap_num(imode)
         qsub4(icnst,jsub) = qnum4(imode)/fcvt_num
         qsub_tendaa(icnst,:,jsub) = qnum_delaa(imode,:)/(fcvt_num*deltat)
      end do ! imode

      !----------------------------------------------------
      ! Interstitial aerosol mass
      !----------------------------------------------------
      do imode = 1, ntot_amode
      do iaer = 1, naer
          icnst = lmap_aer(iaer,imode)
          if (icnst > 0) then
             qsub4(icnst,jsub) = qaer4(iaer,imode)/fcvt_aer(iaer)
             qsub_tendaa(icnst,:,jsub) = qaer_delaa(iaer,imode,:)/(fcvt_aer(iaer)*deltat)
          end if
      end do ! iaer
      end do ! imode

      if ( iscldy_subarea(jsub) ) then

       !----------------------------------------------------
       ! Cloud-borne aerosol number
       !----------------------------------------------------
       do imode = 1, ntot_amode
          icnst = lmap_numcw(imode)
          qqcwsub4(icnst,jsub) = qnumcw4(imode)/fcvt_num
          qqcwsub_tendaa(icnst,:,jsub) = qnumcw_delaa(imode,:)/(fcvt_num*deltat)
       end do ! imode

       !----------------------------------------------------
       ! Cloud-borne aerosol mass
       !----------------------------------------------------
       do imode = 1, ntot_amode
       do iaer = 1, naer
          icnst = lmap_aercw(iaer,imode)
          if (icnst > 0) then
             qqcwsub4(icnst,jsub) = qaercw4(iaer,imode)/fcvt_aer(iaer)
             qqcwsub_tendaa(icnst,:,jsub) = qaercw_delaa(iaer,imode,:)/(fcvt_aer(iaer)*deltat)
          end if
       end do ! iaer
       end do ! imode

      end if

   end do main_jsubarea_loop

end subroutine mam_amicphys_1gridcell
