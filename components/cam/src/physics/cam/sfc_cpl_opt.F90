module sfc_cpl_opt

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use physconst,     only: gravit
  use ppgrid,        only: pcols, pver
  use constituents,  only: pcnst

  implicit none
  public

contains

  subroutine cflx_tend(state, cam_in, ptend)

    use physics_types,   only: physics_state, physics_ptend, physics_ptend_init
    use camsrfexch,      only: cam_in_t 

    implicit none

    type(physics_state), intent(in)     :: state                ! Physics state variables
    type(cam_in_t),      intent(in)     :: cam_in
    type(physics_ptend), intent(out)    :: ptend                ! Individual parameterization tendencies

    logical  :: lq(pcnst)
    integer  :: ncol, m
   !real(r8) :: tmp1(pcols)

    ncol = state%ncol

   !----------
   !lq(:) = .TRUE.
   !call physics_ptend_init(ptend, state%psetcols, 'cflx_tend', lq=lq)

   !rztodt                 = 1._r8/ztodt
   !ptend%q(:ncol,:pver,:) = state%q(:ncol,:pver,:)
   !tmp1(:ncol)            = ztodt * gravit * state%rpdel(:ncol,pver)

   !do m = 2, pcnst
   !  ptend%q(:ncol,pver,m) = ptend%q(:ncol,pver,m) + tmp1(:ncol) * cam_in%cflx(:ncol,m)
   !enddo

   !ptend%q(:ncol,:pver,:) = (ptend%q(:ncol,:pver,:) - state%q(:ncol,:pver,:)) * rztodt
   !----------

    !----------
    lq(:) = .TRUE.
    call physics_ptend_init(ptend, state%psetcols, 'cflx_tend', lq=lq)
    
    do m = 2, pcnst
       ptend%q(:ncol,pver,m) = gravit * state%rpdel(:ncol,pver)* cam_in%cflx(:ncol,m)
    enddo
    !----------


  end subroutine cflx_tend

end module sfc_cpl_opt
