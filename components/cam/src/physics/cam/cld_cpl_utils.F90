module cld_cpl_utils

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use cam_abortutils, only: endrun

  implicit none
  public

  integer,parameter :: DRIB_TQ   = 11
  integer,parameter :: DRIB_UV   = 12
  integer,parameter :: DRIB_TQUV = 13

  integer,parameter :: FORC_TQ   = 21
  integer,parameter :: FORC_UV   = 22
  integer,parameter :: FORC_TQUV = 23

contains

   subroutine set_state_and_tendencies( state, pbuf, cld_cpl_opt, ztodt, p0, rair, cpair, latvap, &
                                        ptend_dribble, thlm_forcing, rtm_forcing, um_forcing, vm_forcing )

   use physics_types,    only: physics_state, physics_ptend, physics_ptend_init
   use physics_buffer,   only: physics_buffer_desc, pbuf_get_field
   use physics_buffer,   only: pbuf_get_index
   use constituents,     only: pcnst, cnst_get_ind
   use ppgrid,           only: pcols, pver

   type(physics_state), intent(inout) :: state
   type(physics_buffer_desc), pointer :: pbuf(:)
   integer, intent(in) :: cld_cpl_opt                    ! scheme identifier 
   real(r8),intent(in) :: ztodt, p0, rair, cpair, latvap

   type(physics_ptend),intent(out) :: ptend_dribble

   real(r8),intent(out) :: thlm_forcing(pcols,pver)
   real(r8),intent(out) ::  rlm_forcing(pcols,pver)
   real(r8),intent(out) ::   um_forcing(pcols,pver)
   real(r8),intent(out) ::   vm_forcing(pcols,pver)

   ! local variables

   integer :: ifld
   integer :: ncol
   integer :: ixcldliq, ixcldice, ixnumliq, ixnumice, ixq
   logical :: lq(pcnst)

    real(r8), pointer, dimension(:,:) :: s_after_macmic
    real(r8), pointer, dimension(:,:) :: t_after_macmic
    real(r8), pointer, dimension(:,:) :: q_after_macmic
    real(r8), pointer, dimension(:,:) :: ql_after_macmic
    real(r8), pointer, dimension(:,:) :: qi_after_macmic
    real(r8), pointer, dimension(:,:) :: nl_after_macmic
    real(r8), pointer, dimension(:,:) :: ni_after_macmic

    real(r8), pointer, dimension(:,:) :: thlm_after_macmic
    real(r8), pointer, dimension(:,:) ::  rtm_after_macmic
    real(r8), pointer, dimension(:,:) ::  utm_after_macmic
    real(r8), pointer, dimension(:,:) ::  vtm_after_macmic

    real(r8) :: thlm_current(pcols,pver)
    real(r8) ::  rlm_current(pcols,pver)
    real(r8) ::   um_current(pcols,pver)
    real(r8) ::   vm_current(pcols,pver)

   !-----------------------------------------------
   ncol = state%ncol

   call cnst_get_ind('Q',      ixq)
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)
   call cnst_get_ind('NUMLIQ', ixnumliq)
   call cnst_get_ind('NUMICE', ixnumice)

   select case (cld_cpl_opt)
   !=============================================
   ! Dribbling for s, T, q, and cloud condensate 
   !=============================================
   case(DRIB_TQ)

      !=====================================================================
      ! thlm forcing and rtm forcing to be passed to CLUBB are set to zero
      !=====================================================================
      thlm_forcing(:,:) = 0._r8
       rtm_forcing(:,:) = 0._r8
        um_forcing(:,:) = 0._r8
        vm_forcing(:,:) = 0._r8

      !=====================================================================================================
      ! s, T, q, and cloud condensate calculate tendencies to be dirbbled and revert state to old snapshot
      !=====================================================================================================
      ifld = pbuf_get_index( 'S_After_MACMIC'); call pbuf_get_field(pbuf, ifld,  s_after_macmic )
      ifld = pbuf_get_index( 'T_After_MACMIC'); call pbuf_get_field(pbuf, ifld,  t_after_macmic )
      ifld = pbuf_get_index( 'Q_After_MACMIC'); call pbuf_get_field(pbuf, ifld,  q_after_macmic )
      ifld = pbuf_get_index('QL_After_MACMIC'); call pbuf_get_field(pbuf, ifld, ql_after_macmic )
      ifld = pbuf_get_index('QI_After_MACMIC'); call pbuf_get_field(pbuf, ifld, qi_after_macmic )
      ifld = pbuf_get_index('NL_After_MACMIC'); call pbuf_get_field(pbuf, ifld, nl_after_macmic )
      ifld = pbuf_get_index('NI_After_MACMIC'); call pbuf_get_field(pbuf, ifld, ni_after_macmic )

      ! Calculate tendencies to be dribbled, and save in ptend_dribble

      lq(:)        = .FALSE.
      lq(ixq)      = .TRUE.
      lq(ixcldliq) = .TRUE.
      lq(ixcldice) = .TRUE.
      lq(ixnumliq) = .TRUE.
      lq(ixnumice) = .TRUE.

      call physics_ptend_init(ptend_dribble, state%psetcols, 'macmic_dribble_tend', ls= .true., lq=lq)

      ptend_dribble%s(:ncol,:pver)          = (state%s(:ncol,:pver)          -   s_after_macmic(:ncol,:pver))  / ztodt
      ptend_dribble%q(:ncol,:pver,ixq)      = (state%q(:ncol,:pver,ixq)      -   q_after_macmic(:ncol,:pver))  / ztodt
      ptend_dribble%q(:ncol,:pver,ixcldliq) = (state%q(:ncol,:pver,ixcldliq) -  ql_after_macmic(:ncol,:pver))  / ztodt
      ptend_dribble%q(:ncol,:pver,ixcldice) = (state%q(:ncol,:pver,ixcldice) -  qi_after_macmic(:ncol,:pver))  / ztodt
      ptend_dribble%q(:ncol,:pver,ixnumliq) = (state%q(:ncol,:pver,ixnumliq) -  nl_after_macmic(:ncol,:pver))  / ztodt
      ptend_dribble%q(:ncol,:pver,ixnumice) = (state%q(:ncol,:pver,ixnumice) -  ni_after_macmic(:ncol,:pver))  / ztodt

      ! Reset state back to the snapshot at old time step  

      state%s(:ncol,:pver)          =  s_after_macmic(:ncol,:pver)
      state%t(:ncol,:pver)          =  t_after_macmic(:ncol,:pver)
      state%q(:ncol,:pver,ixq)      =  q_after_macmic(:ncol,:pver)
      state%q(:ncol,:pver,ixcldliq) = ql_after_macmic(:ncol,:pver)
      state%q(:ncol,:pver,ixcldice) = qi_after_macmic(:ncol,:pver)
      state%q(:ncol,:pver,ixnumliq) = nl_after_macmic(:ncol,:pver)
      state%q(:ncol,:pver,ixnumice) = ni_after_macmic(:ncol,:pver)

   !================
   ! Forcing method 
   !================
   case (FORC_TQ) 

      !==================================================================================================================
      ! thlm and rtm in CLUBB correspond to s, T, q, ql in host model.
      ! Need to calculate tendencies of thlm and rtm resulting from rest of the model,
      ! and revert s, T, q, ql to an old state (the corresponding old thlm and rtm will be calculated in clubb_tend_cam.
      !==================================================================================================================
      ! thlm: calculate tendency (thlm_forcing), theta_l = T* (p0/p)**(Rair/Cpair) - (Lv/Cpair)*ql

      ifld = pbuf_get_index('thlm_After_MACMIC'); call pbuf_get_field(pbuf, ifld, thlm_after_macmic )

      thlm_current(:ncol,:pver) = state%t(:ncol,:pver) * ( p0/state%pmid(:ncol,:per) )**(rair/cpair) &
                                  - (latvap/cpair)*state%q(:ncol,:pver,ixcldliq)

      thlm_forcing(:ncol,:pver) = ( thlm_current(:ncol,:pver) - thlm_after_macmic(:ncol,:pver))/ztodt

      !------------------------------------------------------------
      ! rtm: calculate tendency (rtm_forcing), rtm:  rt = qv + ql

      ifld = pbuf_get_index('rtm_After_MACMIC'); call pbuf_get_field(pbuf, ifld, rtm_after_macmic )

      rtm_current(:ncol,:pver) = state%q(:ncol,pver,ixq) + state%q(:ncol,:pver,ixcldliq) 

      rtm_forcing(:ncol,:pver) = ( rtm_new(:ncol,:pver) - rtm_after_macmic(:ncol,:pver))/ztodt

      !-----------------------------------------------
      ! revert s, T, q, ql in "state" to old values
      !-----------------------------------------------
      ifld = pbuf_get_index( 'S_After_MACMIC'); call pbuf_get_field(pbuf, ifld,  s_after_macmic )
      ifld = pbuf_get_index( 'T_After_MACMIC'); call pbuf_get_field(pbuf, ifld,  t_after_macmic )
      ifld = pbuf_get_index( 'Q_After_MACMIC'); call pbuf_get_field(pbuf, ifld,  q_after_macmic )
      ifld = pbuf_get_index('QL_After_MACMIC'); call pbuf_get_field(pbuf, ifld, ql_after_macmic )

      state%s(:ncol,:pver)          =  s_after_macmic(:ncol,:pver)
      state%t(:ncol,:pver)          =  t_after_macmic(:ncol,:pver)
      state%q(:ncol,:pver,ixq)      =  q_after_macmic(:ncol,:pver)
      state%q(:ncol,:pver,ixcldliq) = ql_after_macmic(:ncol,:pver)

      !==================================================================================
      ! CLUBB does not deal with qi, ni, nl, so these will have to be dribbled for now
      !==================================================================================
      ifld = pbuf_get_index('QI_After_MACMIC'); call pbuf_get_field(pbuf, ifld, qi_after_macmic )
      ifld = pbuf_get_index('NL_After_MACMIC'); call pbuf_get_field(pbuf, ifld, nl_after_macmic )
      ifld = pbuf_get_index('NI_After_MACMIC'); call pbuf_get_field(pbuf, ifld, ni_after_macmic )

      ! Calculate tendencies to be dribbled, and save in ptend_dribble

      lq(:)        = .FALSE.
      lq(ixcldice) = .TRUE.
      lq(ixnumliq) = .TRUE.
      lq(ixnumice) = .TRUE.

      call physics_ptend_init(ptend_dribble, state%psetcols, 'macmic_dribble_tend', ls= .false., lq=lq)

      ptend_dribble%q(:ncol,:pver,ixcldice) = (state%q(:ncol,:pver,ixcldice) -  qi_after_macmic(:ncol,:pver))  / ztodt
      ptend_dribble%q(:ncol,:pver,ixnumliq) = (state%q(:ncol,:pver,ixnumliq) -  nl_after_macmic(:ncol,:pver))  / ztodt
      ptend_dribble%q(:ncol,:pver,ixnumice) = (state%q(:ncol,:pver,ixnumice) -  ni_after_macmic(:ncol,:pver))  / ztodt

      ! Reset state back to the snapshot at old time step  

      state%q(:ncol,:pver,ixcldice) = qi_after_macmic(:ncol,:pver)
      state%q(:ncol,:pver,ixnumliq) = nl_after_macmic(:ncol,:pver)
      state%q(:ncol,:pver,ixnumice) = ni_after_macmic(:ncol,:pver)

      !==================================================================================
      ! u and v: 
      ! Use isolated sequential splitting between macmic subcycles and rest of model.
      ! No dribbling. Forcing terms for CLUBB are set to zero because the state has
      ! already been updated.
      !==================================================================================
      um_forcing(:,:) = 0._r8
      vm_forcing(:,:) = 0._r8

   case default
      call endrun('set_state_and_tendencies: choice of cld_cpl_opt not yet supported.') 
   end select 


   !!if (cld_cpl_opt == FORC_UV) then ! forcing method for u and v

   !!   ifld = pbuf_get_index('utm_After_MACMIC'); call pbuf_get_field(pbuf, ifld, um_after_macmic )
   !!   um_forcing(:ncol,:pver) = ( um_new(:ncol,:pver) - um_after_macmic(:ncol,:pver))/ztodt

   !!   ifld = pbuf_get_index('utm_After_MACMIC'); call pbuf_get_field(pbuf, ifld, vm_after_macmic )
   !!   vm_forcing(:ncol,:pver) = ( vm_new(:ncol,:pver) - vm_after_macmic(:ncol,:pver))/ztodt

   !!end if

   !!if (cld_cpl_opt == DRIB_UV .or. cld_cpl_opt == FORC_UV) then ! use dribbling or forcing method for u and v

   !!   state%u(:ncol,:pver) =  u_after_macmic(:ncol,:pver)
   !!   state%v(:ncol,:pver) =  v_after_macmic(:ncol,:pver)

   !!end if

   end subroutine set_state_and_tendencies

   !---------------------------------------------------------------------------------------------------
   subroutine save_state_snapshot_to_pbuf( state, pbuf, cld_cpl_opt, p0, rair, cpair, latvap )

   use physics_types,  only: physics_state
   use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_get_index
   use ppgrid,         only: pcols, pver
   use constituents,   only: cnst_get_ind

   ! Arguments

   type(physics_state), intent(in)    :: state
   type(physics_buffer_desc), pointer :: pbuf(:)

   integer, intent(in) :: cld_cpl_opt
   real(r8),intent(in) :: p0, rair, cpair, latvap   ! various constants needed for calculating theta_l

   ! Local variables

   real(r8), pointer, dimension(:,:) :: ptr2d
   integer :: ifld, ncol
   integer :: ixcldliq, ixcldice, ixnumliq, ixnumice, ixq

   real(r8) :: thlm_current(pcols,pver)
   real(r8) ::  rlm_current(pcols,pver)

   char(len=20) :: varname

   !---------------------
   ncol = state%ncol

   ! Pbuf variable indices of qv and ql will be used multiple times, so get them already now.

   call cnst_get_ind('Q', ixq)
   call cnst_get_ind('CLDLIQ', ixcldliq)

   ! Save s, T, and cloud condensate 

   varname = 'S_After_MACMIC'; call pbuf_get_field(pbuf, pbuf_get_index(trim(varname)), ptr2d)
   ptr2d(:ncol,:pver) = state%s(:ncol,:pver)

   varname = 'T_After_MACMIC'; call pbuf_get_field(pbuf, pbuf_get_index(trim(varname)), ptr2d)
   ptr2d(:ncol,:pver) = state%t(:ncol,:pver)

   varname = 'Q_After_MACMIC'; call pbuf_get_field(pbuf, pbuf_get_index(trim(varname)), ptr2d)
   ptr2d(:ncol,:pver) = state%q(:ncol,:pver,ixq)

   varname = 'QL_After_MACMIC'; call pbuf_get_field(pbuf, pbuf_get_index(trim(varname)), ptr2d)
   ptr2d(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq)

   varname = 'QI_After_MACMIC'; call pbuf_get_field(pbuf, pbuf_get_index(trim(varname)), ptr2d)
   call cnst_get_ind('CLDICE', ixcldice)
   ptr2d(:ncol,:pver) = state%q(:ncol,:pver,ixcldice)

   varname = 'NL_After_MACMIC'; call pbuf_get_field(pbuf, pbuf_get_index(trim(varname)), ptr2d)
   call cnst_get_ind('NUMLIQ', ixnumliq)
   ptr2d(:ncol,:pver) = state%q(:ncol,:pver,ixnumliq)

   varname = 'NI_After_MACMIC'; call pbuf_get_field(pbuf, pbuf_get_index(trim(varname)), ptr2d)
   call cnst_get_ind('NUMICE', ixnumice)
   ptr2d(:ncol,:pver) = state%q(:ncol,:pver,ixnumice)

   ! Save theta_l and rt.
   ! Because T, q, and ql likely have changed during microphsics calculation,
   ! we need to re-diagnose theta_l and rt instead of saving the values
   ! that came out of the last call of CLUBB.

   if (cld_cpl_opt == FORC_TQ .or. cld_cpl_opt == FORC_TQUV) then

      ! theta_l = T* (p0/p)**(Rair/Cpair) - (Lv/Cpair)*ql

      thlm_current(:ncol,:pver) = state%t(:ncol,:pver) * ( p0/state%pmid(:ncol,:per) )**(rair/cpair) &
                                  - (latvap/cpair)*state%q(:ncol,:pver,ixcldliq)

      varname = 'THLM_After_MACMIC'; call pbuf_get_field(pbuf, pbuf_get_index(trim(varname)), ptr2d)
      ptr2d(:ncol,:pver) = thlm_current(:ncol,:pver)

      ! rt = qv + ql

      rtm_current(:ncol,:pver) = state%q(:ncol,pver,ixq) + state%q(:ncol,:pver,ixcldliq)

      varname = 'RTM_After_MACMIC'; call pbuf_get_field(pbuf, pbuf_get_index(trim(varname)), ptr2d)
      ptr2d(:ncol,:pver) = rtm_current(:ncol,:pver)

   end if

   ! Save u and v

   if ( cld_cpl_opt == DRIB_UV .or. cld_cpl_opt == DRIB_TQUV .or. &
        cld_cpl_opt == FORC_UV .or. cld_cpl_opt == FORC_TQUV      ) then

      varname = 'U_After_MACMIC'; call pbuf_get_field(pbuf, pbuf_get_index(trim(varname)), ptr2d)
      ptr2d(:ncol,:pver) = state%u(:ncol,:pver)

      varname = 'V_After_MACMIC'; call pbuf_get_field(pbuf, pbuf_get_index(trim(varname)), ptr2d)
      ptr2d(:ncol,:pver) = state%v(:ncol,:pver)

   end if

   end subroutine save_state_snapshot_to_pbuf


end module cld_cpl_utils
