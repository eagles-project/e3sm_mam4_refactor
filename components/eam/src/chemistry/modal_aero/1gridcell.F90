
!--------------------------------------------------------------------------------
! Calculates changes to gas and aerosol subarea TMRs (tracer mixing ratios)
!    for all active subareas in the current grid cell (with indices = lchnk,i,k)
! qsub3 and qqcwsub3 are the incoming current TMRs
! qsub4 and qqcwsub4 are the outgoing updated TMRs
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

      real(r8), dimension( 1:max_gas ) :: qgas1_mam, qgas3_mam, qgas4_mam
      real(r8), dimension( 1:max_mode ) :: qnum3_mam, qnum4_mam, qnumcw3_mam, qnumcw4_mam
      real(r8), dimension( 1:max_aer, 1:max_mode ) :: qaer2_mam, qaer3_mam, qaer4_mam, qaercw2_mam, qaercw3_mam, qaercw4_mam
      real(r8), dimension( 1:max_mode ) :: qwtr3_mam, qwtr4_mam

      real(r8), dimension( 1:max_gas, 1:nqtendaa ) ::  qgas_delaa_mam
      real(r8), dimension( 1:max_mode, 1:nqtendaa ) :: qnum_delaa_mam
      real(r8), dimension( 1:max_mode, 1:nqqcwtendaa ) :: qnumcw_delaa_mam
      real(r8), dimension( 1:max_aer, 1:max_mode, 1:nqtendaa ) ::  qaer_delaa_mam
      real(r8), dimension( 1:max_aer, 1:max_mode, 1:nqqcwtendaa ) :: qaercw_delaa_mam

      type ( misc_vars_aa_type ), dimension(nsubarea) :: misc_vars_aa_sub


   !--------------------------------------------------------------------------------
   ! the q--4 values will be equal to q--3 values unless they get changed
   !--------------------------------------------------------------------------------
   qsub4(:,1:nsubarea) = qsub3(:,1:nsubarea)
   qqcwsub4(:,1:nsubarea) = qqcwsub3(:,1:nsubarea)
   qaerwatsub4(:,1:nsubarea) = qaerwatsub3(:,1:nsubarea)

   qsub_tendaa(:,:,1:nsubarea) = 0.0_r8
   qqcwsub_tendaa(:,:,1:nsubarea) = 0.0_r8

   !=================================
   ! Loop over cloudy/clear subareas
   !=================================
   main_jsubarea_loop: &
   do jsub = 1, nsubarea

      !--------------------------------------------------------------------------------
      ! Set switches. These are used in the subroutines called below
      !--------------------------------------------------------------------------------
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

      !--------------------------------------------------------------------------------
      ! Map incoming (host model's) subarea mixing ratios to the aerosol
      ! microphysics package's internal gas/aer/num arrays
      !--------------------------------------------------------------------------------
      call map_info_from_host_to_mam( do_map_gas_sub, iscldy_subarea(jsub),        &! in
                                      qsub1(:,jsub), qsub2(:,jsub), qsub3(:,jsub), &! in
                                      qgas1_mam, qgas3_mam,                        &! out
                                      qnum3_mam, qaer2_mam, qaer3_mam,             &! out
                                      qqcwsub2(:,jsub), qqcwsub3(:,jsub),          &! in
                                      qnumcw3_mam, qaercw2_mam, qaercw3_mam,       &! out
                                      qaerwatsub3(:,jsub),                         &! in
                                      qwtr3_mam                                    )! out


      !--------------------------------------------------------------------------------
      ! Calculate aerosol microphysics for a single subarea
      !--------------------------------------------------------------------------------
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
         qgas1_mam,      qgas3_mam,    qgas4_mam,       qgas_delaa_mam,  &
         qnum3_mam,      qnum4_mam,    qnum_delaa_mam,                   &
         qaer2_mam,      qaer3_mam,    qaer4_mam,       qaer_delaa_mam,  &
         qwtr3_mam,      qwtr4_mam,                                      &
         qnumcw3_mam,    qnumcw4_mam,  qnumcw_delaa_mam,                 &
         qaercw2_mam,    qaercw3_mam,  qaercw4_mam,     qaercw_delaa_mam,&
         misc_vars_aa_sub(jsub)                                          )

      !-------
      if ((nsubarea == 1) .or. (iscldy_subarea(jsub) .eqv. .false.)) then
         misc_vars_aa%ncluster_tend_nnuc_1grid = misc_vars_aa%ncluster_tend_nnuc_1grid &
                                               + misc_vars_aa_sub(jsub)%ncluster_tend_nnuc_1grid*afracsub(jsub)
      end if

      !--------------------------------------------------------------------------------
      ! Map the aerosol microphysics package's internal gas/aer/num arrays (mix-ratios 
      ! and increments) back to (host model's) subarea arrays
      !--------------------------------------------------------------------------------
      call map_info_from_mam_to_host( do_map_gas_sub, iscldy_subarea(jsub), deltat, &! in
                                      qgas4_mam, qgas_delaa_mam,                    &! in
                                      qnum4_mam, qnum_delaa_mam,                    &! in
                                      qaer4_mam, qaer_delaa_mam,                    &! in
                                      qsub4(:,jsub), qsub_tendaa(:,:,jsub),         &! inout
                                      qnumcw4_mam, qnumcw_delaa_mam,                &! in
                                      qaercw4_mam, qaercw_delaa_mam,                &! in
                                      qqcwsub4(:,jsub), qqcwsub_tendaa(:,:,jsub),   &! inout
                                      qwtr4_mam,                                    &! in
                                      qaerwatsub4(:,jsub)                           )! inout

   !=================================
   ! Done with the subarea
   !=================================
   end do main_jsubarea_loop

end subroutine mam_amicphys_1gridcell

!--------------------------------------------------------------------------------
! Purpose: 
! Map mixing ratios from arrays that use the host model's indexing to
! MAM's amicphys internal arrays.
! Also reconcile possible differences in molecular weights.
!--------------------------------------------------------------------------------
subroutine map_info_from_host_to_mam( do_map_gas, do_map_cldbrn,               &! in
                                      q1_host,       q2_host,     q3_host,     &! in
                                      qgas1_mam,     qgas3_mam,                &! out
                                      qnum3_mam,     qaer2_mam,   qaer3_mam,   &! out
                                      qqcw2_host,    qqcw3_host,               &! in
                                      qnumcw3_mam,   qaercw2_mam, qaercw3_mam, &! out
                                      qaerwat3_host,                           &! in
                                      qwtr3_mam                                )! out

   use modal_aero_amicphys_control,only: r8
   use modal_aero_amicphys_control,only: ngas, lmap_gas, fcvt_gas, &
                                         ntot_amode, naer, lmap_num, lmap_aer, fcvt_num, fcvt_aer, &
                                         lmap_numcw, lmap_aercw, &
                                         fcvt_wtr

   logical,intent(in) :: do_map_gas
   logical,intent(in) :: do_map_cldbrn

   real(r8),intent(in)  :: q1_host(:), q2_host(:), q3_host(:)

   real(r8),intent(out) :: qgas1_mam(:), qgas3_mam(:)
   real(r8),intent(out) :: qnum3_mam(:)
   real(r8),intent(out) :: qaer2_mam(:,:), qaer3_mam(:,:)

   real(r8),intent(in)  :: qqcw2_host(:),  qqcw3_host(:)

   real(r8),intent(out) :: qnumcw3_mam(:)
   real(r8),intent(out) :: qaercw2_mam(:,:), qaercw3_mam(:,:)

   real(r8),intent(in)  :: qaerwat3_host(:)
   real(r8),intent(out) :: qwtr3_mam(:)

   integer :: igas, iaer, imode, icnst

   !--------------------
   ! Gases
   !--------------------
   if ( do_map_gas ) then
      do igas = 1, ngas
       icnst = lmap_gas(igas)
       qgas1_mam(igas) = q1_host(icnst)*fcvt_gas(igas)
       qgas3_mam(igas) = q3_host(icnst)*fcvt_gas(igas)
      end do
   end if

   !-----------------------------------------
   ! Interstitial aerosol mass and number
   !-----------------------------------------
   do imode = 1, ntot_amode
      icnst = lmap_num(imode)
      qnum3_mam(imode) = q3_host(icnst)*fcvt_num
   end do ! imode

   do imode = 1, ntot_amode
      do iaer = 1, naer
         icnst = lmap_aer(iaer,imode)
         if (icnst > 0) then
            qaer2_mam(iaer,imode) = q2_host(icnst)*fcvt_aer(iaer)
            qaer3_mam(iaer,imode) = q3_host(icnst)*fcvt_aer(iaer)
         end if
      end do
   end do ! imode

   !-----------------------------------------
   ! Aerosol water
   !-----------------------------------------
   do imode = 1, ntot_amode
      qwtr3_mam(imode) = qaerwat3_host(imode)*fcvt_wtr
   end do ! imode 

   !-----------------------------------------
   ! Cloud-borne aerosol mass and number
   !-----------------------------------------
   if ( do_map_cldbrn ) then

      do imode = 1, ntot_amode
         icnst = lmap_numcw(imode)
         qnumcw3_mam(imode) = qqcw3_host(icnst)*fcvt_num
      end do ! imode

      do imode = 1, ntot_amode
      do iaer = 1, naer
         icnst = lmap_aercw(iaer,imode)
         if (icnst > 0) then
            qaercw2_mam(iaer,imode) = qqcw2_host(icnst)*fcvt_aer(iaer)
            qaercw3_mam(iaer,imode) = qqcw3_host(icnst)*fcvt_aer(iaer)
         end if
      end do ! iaer
      end do ! imode

   end if

end subroutine map_info_from_host_to_mam

!--------------------------------------------------------------------------------
! Purpose: 
! 1. Map mixing ratios from MAM's amicphys internal arrays to arrays
!    that use the host model's indexing. Also reconcile possible
!    differences in molecular weights.
! 2. Do the same for budget terms and convert increments to tendencies.
!--------------------------------------------------------------------------------
subroutine map_info_from_mam_to_host( do_map_gas, do_map_cldbrn, deltat, &! in
                                      qgas_mam,   qgas_delaa_mam,        &! in
                                      qnum_mam,   qnum_delaa_mam,        &! in
                                      qaer_mam,   qaer_delaa_mam,        &! in
                                      q_host,     q_tendaa_host,         &! inout
                                      qnumcw_mam, qnumcw_delaa_mam,      &! in
                                      qaercw_mam, qaercw_delaa_mam,      &! in
                                      qqcw_host,  qqcw_tendaa_host,      &! inout
                                      qwtr_mam,                          &! in
                                      qaerwat_host                       )! inout

    use modal_aero_amicphys_control,only: r8
    use modal_aero_amicphys_control,only: ngas, lmap_gas, fcvt_gas, &
                                          ntot_amode, naer, lmap_num, lmap_aer, fcvt_num, fcvt_aer, &
                                          lmap_numcw, lmap_aercw, &
                                          fcvt_wtr 
   
    logical,intent(in) :: do_map_cldbrn
    logical,intent(in) :: do_map_gas
    real(r8),intent(in):: deltat

    real(r8),intent(in)    :: qgas_mam(:),   qgas_delaa_mam(:,:)
    real(r8),intent(in)    :: qnum_mam(:),   qnum_delaa_mam(:,:)
    real(r8),intent(in)    :: qaer_mam(:,:), qaer_delaa_mam(:,:,:)
    real(r8),intent(inout) :: q_host(:),     q_tendaa_host(:,:)

    real(r8),intent(in)    :: qnumcw_mam(:),   qnumcw_delaa_mam(:,:)
    real(r8),intent(in)    :: qaercw_mam(:,:), qaercw_delaa_mam(:,:,:)
    real(r8),intent(inout) :: qqcw_host(:),    qqcw_tendaa_host(:,:)

    real(r8),intent(in)    :: qwtr_mam(:)
    real(r8),intent(inout) :: qaerwat_host(:)

    integer :: igas, iaer, imode, icnst

    !----------------------------------------------------
    ! Gases
    !----------------------------------------------------
    if ( do_map_gas ) then
     do igas = 1, ngas
        icnst = lmap_gas(igas)
        q_host(icnst) = qgas_mam(igas)/fcvt_gas(igas)
        q_tendaa_host(icnst,:) = qgas_delaa_mam(igas,:)/(fcvt_gas(igas)*deltat)
     end do
    end if

    !----------------------------------------------------
    ! Aerosol water
    !----------------------------------------------------
    do imode = 1, ntot_amode
       qaerwat_host(imode) = qwtr_mam(imode)/fcvt_wtr
    end do ! imode

    !----------------------------------------------------
    ! Interstitial aerosol number
    !----------------------------------------------------
    do imode = 1, ntot_amode
       icnst = lmap_num(imode)
       q_host(icnst) = qnum_mam(imode)/fcvt_num
       q_tendaa_host(icnst,:) = qnum_delaa_mam(imode,:)/(fcvt_num*deltat)
    end do ! imode

    !----------------------------------------------------
    ! Interstitial aerosol mass
    !----------------------------------------------------
    do imode = 1, ntot_amode
    do iaer = 1, naer
        icnst = lmap_aer(iaer,imode)
        if (icnst > 0) then
           q_host(icnst) = qaer_mam(iaer,imode)/fcvt_aer(iaer)
           q_tendaa_host(icnst,:) = qaer_delaa_mam(iaer,imode,:)/(fcvt_aer(iaer)*deltat)
        end if
    end do ! iaer
    end do ! imode


    if ( do_map_cldbrn ) then
     !----------------------------------------------------
     ! Cloud-borne aerosol number
     !----------------------------------------------------
     do imode = 1, ntot_amode
        icnst = lmap_numcw(imode)
        qqcw_host(icnst) = qnumcw_mam(imode)/fcvt_num
        qqcw_tendaa_host(icnst,:) = qnumcw_delaa_mam(imode,:)/(fcvt_num*deltat)
     end do ! imode

     !----------------------------------------------------
     ! Cloud-borne aerosol mass
     !----------------------------------------------------
     do imode = 1, ntot_amode
     do iaer = 1, naer
        icnst = lmap_aercw(iaer,imode)
        if (icnst > 0) then
           qqcw_host(icnst) = qaercw_mam(iaer,imode)/fcvt_aer(iaer)
           qqcw_tendaa_host(icnst,:) = qaercw_delaa_mam(iaer,imode,:)/(fcvt_aer(iaer)*deltat)
        end if
     end do ! iaer
     end do ! imode

    end if

end subroutine map_info_from_mam_to_host


