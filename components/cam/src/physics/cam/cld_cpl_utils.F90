module cld_cpl_utils

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use cam_abortutils, only: endrun

  implicit none
  public

  integer,parameter :: DRIB_TQ   = 11
  integer,parameter :: DRIB_UV   = 12
  integer,parameter :: DRIB_TQUV = 13

  integer,parameter :: FORC_TQ_dTdt   = 211
  integer,parameter :: FORC_TQ_dsdt   = 212
  integer,parameter :: FORC_UV   = 22
  integer,parameter :: FORC_TQUV = 23

contains

   subroutine cld_cpl_register( cld_cpl_opt )

     use physics_buffer, only: pbuf_add_field, dtype_r8
     use ppgrid,         only: pcols, pver

     integer, intent(in) :: cld_cpl_opt

     integer :: idxtmp   ! pbuf component index
 
     if (cld_cpl_opt > 0) then
        call pbuf_add_field(    'S_AFT_MACMIC', 'global', dtype_r8, (/pcols,pver/), idxtmp)
        call pbuf_add_field(    'T_AFT_MACMIC', 'global', dtype_r8, (/pcols,pver/), idxtmp)
        call pbuf_add_field(    'Q_AFT_MACMIC', 'global', dtype_r8, (/pcols,pver/), idxtmp)
        call pbuf_add_field(   'QL_AFT_MACMIC', 'global', dtype_r8, (/pcols,pver/), idxtmp)
        call pbuf_add_field(   'QI_AFT_MACMIC', 'global', dtype_r8, (/pcols,pver/), idxtmp)
        call pbuf_add_field(   'NL_AFT_MACMIC', 'global', dtype_r8, (/pcols,pver/), idxtmp)
        call pbuf_add_field(   'NI_AFT_MACMIC', 'global', dtype_r8, (/pcols,pver/), idxtmp)
        call pbuf_add_field(   'UM_AFT_MACMIC', 'global', dtype_r8, (/pcols,pver/), idxtmp)
        call pbuf_add_field(   'VM_AFT_MACMIC', 'global', dtype_r8, (/pcols,pver/), idxtmp)
        call pbuf_add_field( 'PMID_AFT_MACMIC', 'global', dtype_r8, (/pcols,pver/), idxtmp)
        call pbuf_add_field( 'THLM_AFT_MACMIC', 'global', dtype_r8, (/pcols,pver/), idxtmp)
     end if

   end subroutine cld_cpl_register

   !---------------------------------------------------------------------------------------------------------
   subroutine set_state_and_tendencies( state, pbuf, cld_cpl_opt, ztodt, p0, rair, cpair, latvap, tend, &
                                        ptend_dribble, thlm_forcing, rtm_forcing, um_forcing, vm_forcing )

   use physics_types,    only: physics_state, physics_ptend, physics_ptend_init, physics_tend
   use physics_buffer,   only: physics_buffer_desc, pbuf_get_field
   use physics_buffer,   only: pbuf_get_index
   use constituents,     only: pcnst, cnst_get_ind
   use ppgrid,           only: pcols, pver

   type(physics_state), intent(inout) :: state
   type(physics_buffer_desc), pointer :: pbuf(:)
   integer, intent(in) :: cld_cpl_opt                    ! scheme identifier 
   real(r8),intent(in) :: ztodt, p0, rair, cpair, latvap

   type(physics_tend ),intent(inout) :: tend
   type(physics_ptend),intent(out)   :: ptend_dribble

   real(r8),intent(out) :: thlm_forcing(pcols,pver)
   real(r8),intent(out) ::  rtm_forcing(pcols,pver)
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
    real(r8), pointer, dimension(:,:) ::   um_after_macmic
    real(r8), pointer, dimension(:,:) ::   vm_after_macmic
    real(r8), pointer, dimension(:,:) :: pmid_after_macmic

    real(r8) :: thlm_current(pcols,pver)
    real(r8) ::   um_current(pcols,pver)
    real(r8) ::   vm_current(pcols,pver)

    real(r8) ::   dTdt(pcols,pver)
    real(r8) ::  dqldt(pcols,pver)
    real(r8) ::   dqdt(pcols,pver)

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

      !--------------------------------------------------------------------
      ! Thlm forcing and rtm forcing to be passed to CLUBB are set to zero
      !--------------------------------------------------------------------
      thlm_forcing(:,:) = 0._r8
       rtm_forcing(:,:) = 0._r8
        um_forcing(:,:) = 0._r8
        vm_forcing(:,:) = 0._r8

      !---------------------------------------------------------------------------------------------------------
      ! For s, T, q, and cloud condensate, calculate tendencies to be dribbled and revert state to old snapshot
      !---------------------------------------------------------------------------------------------------------
      ifld = pbuf_get_index( 'S_AFT_MACMIC'); call pbuf_get_field(pbuf, ifld,  s_after_macmic )
      ifld = pbuf_get_index( 'Q_AFT_MACMIC'); call pbuf_get_field(pbuf, ifld,  q_after_macmic )
      ifld = pbuf_get_index('QL_AFT_MACMIC'); call pbuf_get_field(pbuf, ifld, ql_after_macmic )
      ifld = pbuf_get_index('QI_AFT_MACMIC'); call pbuf_get_field(pbuf, ifld, qi_after_macmic )
      ifld = pbuf_get_index('NL_AFT_MACMIC'); call pbuf_get_field(pbuf, ifld, nl_after_macmic )
      ifld = pbuf_get_index('NI_AFT_MACMIC'); call pbuf_get_field(pbuf, ifld, ni_after_macmic )

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

      ! Reset state back to an old snapshot

      state%s(:ncol,:pver)          =  s_after_macmic(:ncol,:pver)
      state%q(:ncol,:pver,ixq)      =  q_after_macmic(:ncol,:pver)
      state%q(:ncol,:pver,ixcldliq) = ql_after_macmic(:ncol,:pver)
      state%q(:ncol,:pver,ixcldice) = qi_after_macmic(:ncol,:pver)
      state%q(:ncol,:pver,ixnumliq) = nl_after_macmic(:ncol,:pver)
      state%q(:ncol,:pver,ixnumice) = ni_after_macmic(:ncol,:pver)

   !================
   ! Forcing method 
   !================
  !case (FORC_TQ_dTdt,FORC_TQ_dsdt) 
   case (211,212,213,214,215,216) 

      !==================================================================================================================
      ! Thlm and rtm in CLUBB correspond to s (T), q, ql in host model.
      ! Need to calculate tendencies of thlm and rtm resulting from rest of the model,
      ! and revert s, q, ql to an old state (the corresponding old thlm and rtm will be calculated in clubb_tend_cam.
      !==================================================================================================================
      ifld = pbuf_get_index( 'T_AFT_MACMIC'); call pbuf_get_field(pbuf, ifld,  t_after_macmic )
      ifld = pbuf_get_index( 'Q_AFT_MACMIC'); call pbuf_get_field(pbuf, ifld,  q_after_macmic )
      ifld = pbuf_get_index('QL_AFT_MACMIC'); call pbuf_get_field(pbuf, ifld, ql_after_macmic )

      !---------------------------------------
      ! Tendency/forcing calculation for rtm 
      !---------------------------------------
       dqdt(:ncol,:pver) = ( state%q(:ncol,:pver,ixq)      -  q_after_macmic(:ncol,:pver) )/ztodt
      dqldt(:ncol,:pver) = ( state%q(:ncol,:pver,ixcldliq) - ql_after_macmic(:ncol,:pver) )/ztodt

      ! rtm = q + ql

      rtm_forcing(:ncol,:pver) = dqldt(:ncol,:pver) + dqdt(:ncol,:pver)

      if (cld_cpl_opt==216) then
      rtm_forcing(:ncol,:pver) = 0._r8
      end if
      !---------------------------------------
      ! Tendency/forcing calculation for thlm 
      !---------------------------------------
      ! thlm = T* (p0/p)**(Rair/Cpair) - (Lv/Cpair)*ql
      !--------------------------------------------------------------------------------------------
      ! Subtract s or T tendency from variable "tend".
      ! (This is needed because "call physics_update()" in tphysbc after "call clubb_tend_cam" 
      ! has an actual argument "tend", meaning that this "call physics_update()" not only updates 
      ! the model state but also accumulate tendencies in "tend". When the forcing method is used,
      ! the out-of-mac-mic tendencies will be included in the ptend returned by 
      ! "call clubb_tend_cam". To avoid double-counting, we need to subtract the
      ! out-of-mac-mic tendencies here.)
      ! Then revert s to an old snapshot.
      !--------------------------------------------------------------------------------------------

      select case(cld_cpl_opt)
      case(211) 

         dTdt(:ncol,:pver) = ( state%t(:ncol,:pver) - t_after_macmic(:ncol,:pver) )/ztodt

         thlm_forcing(:ncol,:pver) =   dTdt(:ncol,:pver) * ( p0/state%pmid(:ncol,:pver) )**(rair/cpair) &
                                    - dqldt(:ncol,:pver)* (latvap/cpair)

         tend%dtdt(:ncol,:pver) = tend%dtdt(:ncol,:pver) - dTdt(:ncol,:pver) 
           state%s(:ncol,:pver) =   state%s(:ncol,:pver) - cpair*dTdt(:ncol,:pver)*ztodt

      case(212)

         dTdt(:ncol,:pver) = ( state%t(:ncol,:pver) - t_after_macmic(:ncol,:pver) )/ztodt

         ifld = pbuf_get_index('PMID_AFT_MACMIC'); call pbuf_get_field(pbuf, ifld, pmid_after_macmic )
         thlm_forcing(:ncol,:pver) =   dTdt(:ncol,:pver) * ( p0/pmid_after_macmic(:ncol,:pver) )**(rair/cpair) &
                                    - dqldt(:ncol,:pver)* (latvap/cpair)

         tend%dtdt(:ncol,:pver) = tend%dtdt(:ncol,:pver) - dTdt(:ncol,:pver) 
           state%s(:ncol,:pver) =   state%s(:ncol,:pver) - cpair*dTdt(:ncol,:pver)*ztodt

      case(213,216)

         dTdt(:ncol,:pver) = ( state%t(:ncol,:pver) - t_after_macmic(:ncol,:pver) )/ztodt

         ifld = pbuf_get_index('PMID_AFT_MACMIC'); call pbuf_get_field(pbuf, ifld, pmid_after_macmic )
         thlm_forcing(:ncol,:pver) =   dTdt(:ncol,:pver) * ( p0/pmid_after_macmic(:ncol,:pver) )**(rair/cpair) !&
                                  ! - dqldt(:ncol,:pver)* (latvap/cpair)

         tend%dtdt(:ncol,:pver) = tend%dtdt(:ncol,:pver) - dTdt(:ncol,:pver) 
           state%s(:ncol,:pver) =   state%s(:ncol,:pver) - cpair*dTdt(:ncol,:pver)*ztodt

      case(214)

         dTdt(:ncol,:pver) = ( state%t(:ncol,:pver) - t_after_macmic(:ncol,:pver) )/ztodt

         ifld = pbuf_get_index('PMID_AFT_MACMIC'); call pbuf_get_field(pbuf, ifld, pmid_after_macmic )
         thlm_forcing(:ncol,:pver) =   dTdt(:ncol,:pver) * ( p0/pmid_after_macmic(:ncol,:pver) )**(rair/cpair) &
                                    - dqldt(:ncol,:pver)* (latvap/cpair)


         ifld = pbuf_get_index( 'S_AFT_MACMIC'); call pbuf_get_field(pbuf, ifld,  s_after_macmic )
         tend%dtdt(:ncol,:pver) = tend%dtdt(:ncol,:pver) - (state%s(:ncol,:pver)-s_after_macmic(:ncol,:pver))/ztodt/cpair
           state%s(:ncol,:pver) = s_after_macmic(:ncol,:pver)

      case(215)

         ifld = pbuf_get_index('THLM_AFT_MACMIC'); call pbuf_get_field(pbuf, ifld, thlm_after_macmic )

         thlm_current(:ncol,:pver) =  state%t(:ncol,:pver) * ( p0/state%pmid(:ncol,:pver) )**(rair/cpair) &
                                    - state%q(:ncol,:pver,ixcldliq)* (latvap/cpair)

         thlm_forcing(:ncol,:pver) = ( thlm_current(:ncol,:pver) - thlm_after_macmic(:ncol,:pver) )/ztodt

         ifld = pbuf_get_index( 'S_AFT_MACMIC'); call pbuf_get_field(pbuf, ifld,  s_after_macmic )
         tend%dtdt(:ncol,:pver) = tend%dtdt(:ncol,:pver) - (state%s(:ncol,:pver)-s_after_macmic(:ncol,:pver))/ztodt/cpair
           state%s(:ncol,:pver) = s_after_macmic(:ncol,:pver)

      case default
         call endrun('set_state_and_tendencies: unrecognized value for cld_cpl_opt')
      end select

      !==================================================================================
      ! u and v: 
      !==================================================================================
      ! Use isolated sequential splitting between macmic subcycles and rest of model.
      ! No dribbling. Forcing terms for CLUBB are set to zero because the state has
      ! already been updated.
      !-----------------------------------------------------------------------------------
      um_forcing(:,:) = 0._r8
      vm_forcing(:,:) = 0._r8

      !----------------------------------------------------------
      ! PLACE HOLDER: if applying forcing method to u, v, then 
      ! - calculate um_forcing, vm_forcing
      ! - subtract from tend%dudt and tend%dvdt
      ! - revert state%u and state%v to old values
      !----------------------------------------------------------

      !==================================================================================
      ! CLUBB does not deal with qi, ni, nl, so these will have to be dribbled for now.
      ! Also, s and qv are reverted back to old values, we need to re-diagnose T and geopotential height,
      ! for which we set ls = .true. with zero s tendency in ptend_dribble.
      !==================================================================================
      ! Calculate condensate tendencies to be dribbled, and save in ptend_dribble

      lq(:)        = .FALSE.
      lq(ixcldice) = .TRUE.
      lq(ixnumliq) = .TRUE.
      lq(ixnumice) = .TRUE.

      if (cld_cpl_opt==216) then
         lq(ixq) = .true.
         lq(ixcldliq) = .true.
      end if

      call physics_ptend_init(ptend_dribble, state%psetcols, 'macmic_dribble_tend', ls=.true., lq=lq)

      ptend_dribble%s(:ncol,:pver) = 0._r8

      if (cld_cpl_opt==216) then
      ptend_dribble%q(:ncol,:pver,ixq)      = (state%q(:ncol,:pver,ixq)      -   q_after_macmic(:ncol,:pver))  / ztodt
      ptend_dribble%q(:ncol,:pver,ixcldliq) = (state%q(:ncol,:pver,ixcldliq) -  ql_after_macmic(:ncol,:pver))  / ztodt
      end if

      ifld = pbuf_get_index('QI_AFT_MACMIC'); call pbuf_get_field(pbuf, ifld, qi_after_macmic )
      ptend_dribble%q(:ncol,:pver,ixcldice) = (state%q(:ncol,:pver,ixcldice) -  qi_after_macmic(:ncol,:pver))  / ztodt

      ifld = pbuf_get_index('NL_AFT_MACMIC'); call pbuf_get_field(pbuf, ifld, nl_after_macmic )
      ptend_dribble%q(:ncol,:pver,ixnumliq) = (state%q(:ncol,:pver,ixnumliq) -  nl_after_macmic(:ncol,:pver))  / ztodt

      ifld = pbuf_get_index('NI_AFT_MACMIC'); call pbuf_get_field(pbuf, ifld, ni_after_macmic )
      ptend_dribble%q(:ncol,:pver,ixnumice) = (state%q(:ncol,:pver,ixnumice) -  ni_after_macmic(:ncol,:pver))  / ztodt

      ! Reset state back to the snapshot at old time step  

      state%q(:ncol,:pver,ixcldice) = qi_after_macmic(:ncol,:pver)
      state%q(:ncol,:pver,ixnumliq) = nl_after_macmic(:ncol,:pver)
      state%q(:ncol,:pver,ixnumice) = ni_after_macmic(:ncol,:pver)

      !-----------------------------------------------
      ! Revert q, ql in "state" to old values
      !-----------------------------------------------
      state%q(:ncol,:pver,ixq)      =  q_after_macmic(:ncol,:pver)
      state%q(:ncol,:pver,ixcldliq) = ql_after_macmic(:ncol,:pver)

   case default
      call endrun('set_state_and_tendencies: choice of cld_cpl_opt not yet supported.') 
   end select 


   !!if (cld_cpl_opt == FORC_UV) then ! forcing method for u and v

   !!   ifld = pbuf_get_index('UM_AFT_MACMIC'); call pbuf_get_field(pbuf, ifld, um_after_macmic )
   !!   um_forcing(:ncol,:pver) = ( um_new(:ncol,:pver) - um_after_macmic(:ncol,:pver))/ztodt

   !!   ifld = pbuf_get_index('VM_AFT_MACMIC'); call pbuf_get_field(pbuf, ifld, vm_after_macmic )
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

   real(r8) :: thlm(pcols,pver)

   character(len=20) :: varname

   !---------------------
   ncol = state%ncol

   ! Pbuf variable indices of qv and ql will be used multiple times, so get them already now.

   call cnst_get_ind('Q', ixq)
   call cnst_get_ind('CLDLIQ', ixcldliq)

   ! Save s, T, and cloud condensate 

   varname = 'S_AFT_MACMIC'; call pbuf_get_field(pbuf, pbuf_get_index(trim(varname)), ptr2d)
   ptr2d(:ncol,:pver) = state%s(:ncol,:pver)

   varname = 'T_AFT_MACMIC'; call pbuf_get_field(pbuf, pbuf_get_index(trim(varname)), ptr2d)
   ptr2d(:ncol,:pver) = state%t(:ncol,:pver)

   varname = 'Q_AFT_MACMIC'; call pbuf_get_field(pbuf, pbuf_get_index(trim(varname)), ptr2d)
   ptr2d(:ncol,:pver) = state%q(:ncol,:pver,ixq)

   varname = 'QL_AFT_MACMIC'; call pbuf_get_field(pbuf, pbuf_get_index(trim(varname)), ptr2d)
   ptr2d(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq)

   varname = 'QI_AFT_MACMIC'; call pbuf_get_field(pbuf, pbuf_get_index(trim(varname)), ptr2d)
   call cnst_get_ind('CLDICE', ixcldice)
   ptr2d(:ncol,:pver) = state%q(:ncol,:pver,ixcldice)

   varname = 'NL_AFT_MACMIC'; call pbuf_get_field(pbuf, pbuf_get_index(trim(varname)), ptr2d)
   call cnst_get_ind('NUMLIQ', ixnumliq)
   ptr2d(:ncol,:pver) = state%q(:ncol,:pver,ixnumliq)

   varname = 'NI_AFT_MACMIC'; call pbuf_get_field(pbuf, pbuf_get_index(trim(varname)), ptr2d)
   call cnst_get_ind('NUMICE', ixnumice)
   ptr2d(:ncol,:pver) = state%q(:ncol,:pver,ixnumice)

   varname = 'PMID_AFT_MACMIC'; call pbuf_get_field(pbuf, pbuf_get_index(trim(varname)), ptr2d)
   ptr2d(:ncol,:pver) = state%pmid(:ncol,:pver)

   ! Save thlm = T* (p0/p)**(Rair/Cpair) - (Lv/Cpair)*ql

   thlm(:ncol,:pver) =  state%t(:ncol,:pver) * ( p0/state%pmid(:ncol,:pver) )**(rair/cpair) &
                      - state%q(:ncol,:pver,ixcldliq)* (latvap/cpair)

   varname = 'THLM_AFT_MACMIC'; call pbuf_get_field(pbuf, pbuf_get_index(trim(varname)), ptr2d)
   ptr2d(:ncol,:pver) = thlm(:ncol,:pver)


   ! Save u and v

   if ( cld_cpl_opt == DRIB_UV .or. cld_cpl_opt == DRIB_TQUV .or. &
        cld_cpl_opt == FORC_UV .or. cld_cpl_opt == FORC_TQUV      ) then

      varname = 'UM_AFT_MACMIC'; call pbuf_get_field(pbuf, pbuf_get_index(trim(varname)), ptr2d)
      ptr2d(:ncol,:pver) = state%u(:ncol,:pver)

      varname = 'VM_AFT_MACMIC'; call pbuf_get_field(pbuf, pbuf_get_index(trim(varname)), ptr2d)
      ptr2d(:ncol,:pver) = state%v(:ncol,:pver)

   end if

   end subroutine save_state_snapshot_to_pbuf
   !---------------------------------------------------------------------------------------------------------


end module cld_cpl_utils
