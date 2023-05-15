module mo_imp_sol
  use shr_kind_mod, only : r8 => shr_kind_r8
  use chem_mods, only : clscnt4, gas_pcnst, clsmap
  use cam_logfile, only : iulog

!==================================================================================
! module variables
  implicit none
  private
  public :: imp_slv_inti, imp_sol
  save
  real(r8), parameter :: rel_err = 1.e-3_r8
  real(r8), parameter :: high_rel_err = 1.e-4_r8
  !-----------------------------------------------------------------------
  ! Newton-Raphson iteration limits
  !-----------------------------------------------------------------------
  integer, parameter :: itermax = 11
  integer, parameter :: cut_limit = 5
  real(r8) :: small
  real(r8) :: epsilon(clscnt4)
  logical :: factor(itermax)
  integer :: ox_ndx
  integer :: o1d_ndx = -1
  integer :: h2o_ndx = -1
  integer :: oh_ndx, ho2_ndx, ch3o2_ndx, po2_ndx, ch3co3_ndx
  integer :: c2h5o2_ndx, isopo2_ndx, macro2_ndx, mco3_ndx, c3h7o2_ndx
  integer :: ro2_ndx, xo2_ndx, no_ndx, no2_ndx, no3_ndx, n2o5_ndx
  integer :: c2h4_ndx, c3h6_ndx, isop_ndx, mvk_ndx, c10h16_ndx
  integer :: ox_p1_ndx, ox_p2_ndx, ox_p3_ndx, ox_p4_ndx, ox_p5_ndx
  integer :: ox_p6_ndx, ox_p7_ndx, ox_p8_ndx, ox_p9_ndx, ox_p10_ndx
  integer :: ox_p11_ndx
  integer :: ox_l1_ndx, ox_l2_ndx, ox_l3_ndx, ox_l4_ndx, ox_l5_ndx
  integer :: ox_l6_ndx, ox_l7_ndx, ox_l8_ndx, ox_l9_ndx, usr4_ndx
  integer :: usr16_ndx, usr17_ndx, r63_ndx,c2o3_ndx,ole_ndx
  integer :: tolo2_ndx, terpo2_ndx, alko2_ndx, eneo2_ndx, eo2_ndx, meko2_ndx
  integer :: ox_p17_ndx,ox_p12_ndx,ox_p13_ndx,ox_p14_ndx,ox_p15_ndx,ox_p16_ndx
  integer :: lt_cnt
  ! for xnox ozone chemistry diagnostics
  integer :: o3a_ndx, xno2_ndx, no2xno3_ndx, xno2no3_ndx, xno3_ndx, o1da_ndx, xno_ndx
  integer :: usr4a_ndx, usr16a_ndx, usr16b_ndx, usr17b_ndx
contains

!==================================================================================
  subroutine imp_slv_inti
    !-----------------------------------------------------------------------
    ! ... Initialize the implict solver
    !-----------------------------------------------------------------------
    use mo_chem_utls, only : get_spc_ndx, get_rxt_ndx
    use cam_abortutils, only : endrun
    use cam_history, only : addfld, add_default
    use ppgrid, only : pver
    use mo_tracname, only : solsym
    implicit none
    !-----------------------------------------------------------------------
    ! ... Local variables
    !-----------------------------------------------------------------------
    integer :: m
    real(r8) :: eps(gas_pcnst)
    integer :: wrk(27)
    integer :: i,j
    ! small = 1.e6_r8 * tiny( small )
    small = 1.e-40_r8
    factor(:) = .true.
    eps(:) = rel_err
    ox_ndx = get_spc_ndx( 'OX' )
    h2o_ndx = get_spc_ndx( 'H2O' )
    if( ox_ndx < 1 ) then
       ox_ndx = get_spc_ndx( 'O3' )
       o1d_ndx = get_spc_ndx( 'O1D' )
       o1da_ndx = get_spc_ndx( 'O1DA' )
    end if
    if( ox_ndx > 0 ) then
       eps(ox_ndx) = high_rel_err
    end if
    do m = 1,clscnt4
       epsilon(m) = eps(clsmap(m,4))
    end do
    has_o3_chem: if( ox_ndx > 0 ) then
       ox_p1_ndx = get_rxt_ndx( 'ox_p1' )
       ox_p2_ndx = get_rxt_ndx( 'ox_p2' )
       ox_p3_ndx = get_rxt_ndx( 'ox_p3' )
       ox_p4_ndx = get_rxt_ndx( 'ox_p4' )
       ox_p5_ndx = get_rxt_ndx( 'ox_p5' )
       ox_p6_ndx = get_rxt_ndx( 'ox_p6' )
       ox_p7_ndx = get_rxt_ndx( 'ox_p7' )
       ox_p8_ndx = get_rxt_ndx( 'ox_p8' )
       ox_p9_ndx = get_rxt_ndx( 'ox_p9' )
       ox_p10_ndx = get_rxt_ndx( 'ox_p10' )
       ox_p11_ndx = get_rxt_ndx( 'ox_p11' )
       ox_p12_ndx = get_rxt_ndx( 'ox_p12' )
       ox_p13_ndx = get_rxt_ndx( 'ox_p13' )
       ox_p14_ndx = get_rxt_ndx( 'ox_p14' )
       ox_p15_ndx = get_rxt_ndx( 'ox_p15' )
       ox_p16_ndx = get_rxt_ndx( 'ox_p16' )
       ox_p17_ndx = get_rxt_ndx( 'ox_p17' )
       wrk(1:17) = (/ ox_p1_ndx, ox_p2_ndx, ox_p3_ndx, ox_p4_ndx, ox_p5_ndx, &
            ox_p6_ndx, ox_p7_ndx, ox_p8_ndx, ox_p9_ndx, ox_p10_ndx, ox_p11_ndx, &
            ox_p12_ndx, ox_p13_ndx, ox_p14_ndx, ox_p15_ndx, ox_p16_ndx, ox_p17_ndx /)
       r63_ndx = get_rxt_ndx( 'r63' )
       wrk(1:4) = (/ ox_p1_ndx, ox_p2_ndx, ox_p3_ndx, r63_ndx/)
    endif has_o3_chem
    do i = 1,clscnt4
       j = clsmap(i,4)
       call addfld( trim(solsym(j))//'_CHMP', (/ 'lev' /), 'I', '/cm3/s', 'chemical production rate' )
       call addfld( trim(solsym(j))//'_CHML', (/ 'lev' /), 'I', '/cm3/s', 'chemical loss rate' )
    enddo
  end subroutine imp_slv_inti
!==================================================================================
  subroutine imp_sol( base_sol,                 & ! inout
                      reaction_rates, het_rates, extfrc, delt, & ! in
                      xhnm, ncol, lchnk, ltrop ) ! in
    !-----------------------------------------------------------------------
    ! ... imp_sol advances the volumetric mixing ratio
    ! forward one time step via the fully implicit euler scheme.
    ! this source is meant for small l1 cache machines such as
    ! the intel pentium and itanium cpus
    !-----------------------------------------------------------------------
    use chem_mods,      only : rxntot, extcnt, nzcnt, permute, cls_rxt_cnt
    use mo_tracname,    only : solsym
    use ppgrid,         only : pver
    use mo_lin_matrix,  only : linmat
!    use mo_nln_matrix,  only : nlnmat
!    use mo_lu_factor,   only : lu_fac
!    use mo_lu_solve,    only : lu_slv
!    use mo_prod_loss,   only : imp_prod_loss
    use mo_indprd,      only : indprd
    use time_manager,   only : get_nstep
    use cam_history,    only : outfld

    implicit none
    !-----------------------------------------------------------------------
    ! ... dummy args
    !-----------------------------------------------------------------------
    integer,  intent(in) :: ncol ! columns in chunck
    integer,  intent(in) :: lchnk ! chunk id
    real(r8), intent(in) :: delt ! time step (s)
    real(r8), intent(in) :: reaction_rates(ncol,pver,max(1,rxntot)), & ! rxt rates (1/cm^3/s)
                            extfrc(ncol,pver,max(1,extcnt)), & ! external in-situ forcing (1/cm^3/s)
                            het_rates(ncol,pver,max(1,gas_pcnst)) ! washout rates (1/s)
    real(r8), intent(in) :: xhnm(ncol,pver)
    integer,  intent(in) :: ltrop(ncol) ! chemistry troposphere boundary (index)
    real(r8), intent(inout) :: base_sol(ncol,pver,gas_pcnst) ! species mixing ratios (vmr)

    !-----------------------------------------------------------------------
    ! ... local variables
    !-----------------------------------------------------------------------
    integer :: nr_iter, &
               lev, icol, &
               jj, kk, mm 
    integer :: fail_cnt, cut_cnt, stp_con_cnt
    integer :: nstep
    real(r8) :: interval_done, dt, dti
    real(r8) :: max_delta(max(1,clscnt4))
    real(r8) :: sys_jac(max(1,nzcnt))
    real(r8) :: lin_jac(max(1,nzcnt))
    real(r8), dimension(max(1,clscnt4)) :: &
                         solution, &
                         forcing, &
                         iter_invariant, &
                         prod, loss
    real(r8) :: lrxt(max(1,rxntot))
    real(r8) :: lsol(max(1,gas_pcnst))
    real(r8) :: lhet(max(1,gas_pcnst))
    real(r8), dimension(ncol,pver,max(1,clscnt4)) :: ind_prd
    real(r8), dimension(ncol,pver,max(1,clscnt4)) :: prod_out, loss_out
    logical :: frc_mask
    logical :: converged(max(1,clscnt4))
    logical :: convergence      ! all converged(:) are true

    ! initiate variables
    prod_out(:,:,:) = 0._r8
    loss_out(:,:,:) = 0._r8
    solution(:) = 0._r8
    !-----------------------------------------------------------------------
    ! ... class independent forcing
    !-----------------------------------------------------------------------
    call indprd( 4, ind_prd, clscnt4, base_sol, extfrc, &
            reaction_rates, ncol )
    level_loop : do lev = 1,pver
       column_loop : do icol = 1,ncol
          if (lev <= ltrop(icol)) cycle column_loop
          !-----------------------------------------------------------------------
          ! ... transfer from base to local work arrays
          !-----------------------------------------------------------------------
          do mm = 1,rxntot
             lrxt(mm) = reaction_rates(icol,lev,mm)
          enddo
          if( gas_pcnst > 0 ) then
             do mm = 1,gas_pcnst
                lhet(mm) = het_rates(icol,lev,mm)
             enddo
          endif
          !-----------------------------------------------------------------------
          ! ... time step loop
          !-----------------------------------------------------------------------
          dt = delt
          cut_cnt = 0
          fail_cnt = 0
          stp_con_cnt = 0
          interval_done = 0._r8
          time_step_loop : do
             dti = 1._r8 / dt
             !-----------------------------------------------------------------------
             ! ... transfer from base to local work arrays
             !-----------------------------------------------------------------------
             do mm = 1,gas_pcnst
                lsol(mm) = base_sol(icol,lev,mm)
             enddo
             !-----------------------------------------------------------------------
             ! ... transfer from base to class array
             !-----------------------------------------------------------------------
             do kk = 1,clscnt4
                jj = clsmap(kk,4)
                mm = permute(kk,4)
                solution(mm) = lsol(jj)
             enddo
             !-----------------------------------------------------------------------
             ! ... set the iteration invariant part of the function f(y)
             !-----------------------------------------------------------------------
             do mm = 1,clscnt4
                   iter_invariant(mm) = dti * solution(mm) + ind_prd(icol,lev,mm)
             enddo
             !-----------------------------------------------------------------------
             ! ... the linear component
             !-----------------------------------------------------------------------
             call linmat( lin_jac,           & ! out
                          lsol, lrxt, lhet   ) ! in
             !=======================================================================
             ! the newton-raphson iteration for f(y) = 0
             !=======================================================================

             call newton_raphson_iter( dti, lin_jac, lrxt, lhet, & ! in
                                  iter_invariant,                & ! in
                                  lsol,   solution,              & ! inout
                                  converged, convergence,        & ! out
                                  prod, loss, max_delta          ) ! out

             !-----------------------------------------------------------------------
             ! ... check for newton-raphson convergence
             !-----------------------------------------------------------------------
             if( .not. convergence ) then
                !-----------------------------------------------------------------------
                ! ... non-convergence
                !-----------------------------------------------------------------------
                fail_cnt = fail_cnt + 1
                nstep = get_nstep()
                write(iulog,'('' imp_sol: Time step '',1p,e21.13,'' failed to converge @ (lchnk,lev,col,nstep) = '',4i6)') &
                     dt,lchnk,lev,icol,nstep
                stp_con_cnt = 0
                if( cut_cnt < cut_limit ) then
                   cut_cnt = cut_cnt + 1
                   if( cut_cnt < cut_limit ) then
                      dt = .5_r8 * dt
                   else
                      dt = .1_r8 * dt
                   endif
                   cycle time_step_loop
                else
                   write(iulog,'('' imp_sol: Failed to converge @ (lchnk,lev,col,nstep,dt,time) = '',4i6,1p,2e21.13)') &
                        lchnk,lev,icol,nstep,dt,interval_done+dt
                   do mm = 1,clscnt4
                      if( .not. converged(mm) ) then
                         write(iulog,'(1x,a8,1x,1pe10.3)') solsym(clsmap(mm,4)), max_delta(mm)
                      endif
                   enddo
                endif
             endif ! if( .not. convergence )
             !-----------------------------------------------------------------------
             ! ... check for interval done
             !-----------------------------------------------------------------------
             interval_done = interval_done + dt
             if( abs( delt - interval_done ) <= .0001_r8 ) then
                if( fail_cnt > 0 ) then
                   write(iulog,*) 'imp_sol : @ (lchnk,lev,col) = ',lchnk,lev,icol,' failed ',fail_cnt,' times'
                endif
                exit time_step_loop
             else
                !-----------------------------------------------------------------------
                ! ... transfer latest solution back to base array
                !-----------------------------------------------------------------------
                if( convergence ) then
                   stp_con_cnt = stp_con_cnt + 1
                endif
                do mm = 1,gas_pcnst
                   base_sol(icol,lev,mm) = lsol(mm)
                enddo
                if( stp_con_cnt >= 2 ) then
                   dt = 2._r8*dt
                   stp_con_cnt = 0
                endif
                dt = min( dt,delt-interval_done )
             endif
          enddo time_step_loop
          !-----------------------------------------------------------------------
          ! ... Transfer latest solution back to base array
          !     and calculate Prod/Loss history buffers
          !-----------------------------------------------------------------------
          cls_loop: do kk = 1,clscnt4
             jj = clsmap(kk,4)
             mm = permute(kk,4)
             ! ... Transfer latest solution back to base array
             base_sol(icol,lev,jj) = solution(mm)
             ! ... Prod/Loss history buffers...
             prod_out(icol,lev,kk) = prod(mm) + ind_prd(icol,lev,mm)
             loss_out(icol,lev,kk) = loss(mm)
          enddo cls_loop

       enddo column_loop
    enddo level_loop

    ! diagnose variables
    do icol = 1,clscnt4
       jj = clsmap(icol,4)
       prod_out(:,:,icol) = prod_out(:,:,icol)*xhnm
       loss_out(:,:,icol) = loss_out(:,:,icol)*xhnm
       call outfld( trim(solsym(jj))//'_CHMP', prod_out(:,:,icol), ncol, lchnk )
       call outfld( trim(solsym(jj))//'_CHML', loss_out(:,:,icol), ncol, lchnk )
    enddo
  end subroutine imp_sol

!==================================================================================
  subroutine newton_raphson_iter( dti, lin_jac, lrxt, lhet, & ! in
                        iter_invariant,                & ! in
                        lsol,   solution,               & ! inout
                        converged, convergence,         & ! out
                        prod, loss, max_delta           ) ! out
!-----------------------------------------------------
! the newton-raphson iteration for f(y) = 0
!-----------------------------------------------------
    use chem_mods,      only : rxntot, nzcnt, permute
    use mo_nln_matrix,  only : nlnmat
    use mo_lu_factor,   only : lu_fac
    use mo_lu_solve,    only : lu_slv
    use mo_prod_loss,   only : imp_prod_loss

    real(r8),intent(in) :: lin_jac(:)
    real(r8),intent(in) :: dti
    real(r8),intent(in) :: lrxt(max(1,rxntot))
    real(r8),intent(in) :: lhet(max(1,gas_pcnst))
    real(r8),intent(in) :: iter_invariant(max(1,clscnt4))

    real(r8),intent(inout) :: solution(max(1,clscnt4))
    real(r8),intent(inout) :: lsol(max(1,gas_pcnst))

    logical, intent(out) :: converged(max(1,clscnt4))
    logical, intent(out) :: convergence      ! all converged(:) are true
    real(r8),intent(out) :: prod(max(1,clscnt4))
    real(r8),intent(out) :: loss(max(1,clscnt4))
    real(r8),intent(out) :: max_delta(max(1,clscnt4))

    ! ... local variables
    integer :: nr_iter
    integer :: jj,kk,mm
    logical :: frc_mask

    real(r8) :: sys_jac(max(1,nzcnt))
    real(r8) :: forcing(max(1,clscnt4))



    iter_loop : do nr_iter = 1,itermax
         !-----------------------------------------------------------------------
         ! ... the non-linear component
         !-----------------------------------------------------------------------
         if( factor(nr_iter) ) then
            call nlnmat( sys_jac,     & ! out
                         lin_jac, dti ) ! in
            !-----------------------------------------------------------------------
            ! ... factor the "system" matrix
            !-----------------------------------------------------------------------
            call lu_fac( sys_jac )
         endif
         !-----------------------------------------------------------------------
         ! ... form f(y)
         !-----------------------------------------------------------------------
         call imp_prod_loss( prod, loss,  & ! out
                         lsol, lrxt, lhet ) ! in
         do mm = 1,clscnt4
            forcing(mm) = solution(mm)*dti - (iter_invariant(mm) + prod(mm) - loss(mm))
         enddo
         !-----------------------------------------------------------------------
         ! ... solve for the mixing ratio at t(n+1)
         !-----------------------------------------------------------------------
         call lu_slv( sys_jac, forcing )
         do mm = 1,clscnt4
            solution(mm) = solution(mm) + forcing(mm)
         enddo
         !-----------------------------------------------------------------------
         ! ... convergence measures
         !-----------------------------------------------------------------------
         if( nr_iter > 1 ) then
            do kk = 1,clscnt4
               mm = permute(kk,4)
               if( abs(solution(mm)) > 1.e-20_r8 ) then
                  max_delta(kk) = abs( forcing(mm)/solution(mm) )
               else
                  max_delta(kk) = 0._r8
               endif
            enddo
         endif
         !-----------------------------------------------------------------------
         ! ... limit iterate
         !-----------------------------------------------------------------------
         where( solution(:) < 0._r8 )
            solution(:) = 0._r8
         endwhere
         !-----------------------------------------------------------------------
         ! ... transfer latest solution back to work array
         !-----------------------------------------------------------------------
         do kk = 1,clscnt4
            jj = clsmap(kk,4)
            mm = permute(kk,4)
            lsol(jj) = solution(mm)
         enddo
         !-----------------------------------------------------------------------
         ! ... check for convergence
         !-----------------------------------------------------------------------
         converged(:) = .true.
         if( nr_iter > 1 ) then
            do kk = 1,clscnt4
               mm = permute(kk,4)
               frc_mask = abs( forcing(mm) ) > small
               if( frc_mask ) then
                  converged(kk) = abs(forcing(mm)) <= epsilon(kk)*abs(solution(mm))
               else
                  converged(kk) = .true.
               endif
            enddo
            convergence = all( converged(:) )
            if( convergence ) then
               exit
            endif
         endif
      enddo iter_loop

  end subroutine newton_raphson_iter
!==================================================================================
end module mo_imp_sol
