module misc_diagnostics

use shr_kind_mod,   only: r8 => shr_kind_r8

implicit none
public

integer, parameter :: iCAPE_NEW_PCL_NEW_ENV   = 0
integer, parameter :: iCAPE_NEW_PCL_FIXED_ENV = 1
integer, parameter :: idCAPEe                 = 2
integer, parameter :: idCAPEp                 = 3

contains


!------------------------------------------------
! saturation specific humidity wrt ice.
! The calculation in this subroutine follows 
! what is used in the model for 
!  - history output
!  - ice nucleation parameterization
! 
subroutine qsat_ice( ncol, pver, tair, pair, qsati )

  use wv_saturation, only: qsat_water, svp_ice

  integer, intent(in)  :: ncol, pver

  real(r8),intent(in)  :: tair(ncol,pver)
  real(r8),intent(in)  :: pair(ncol,pver)

  real(r8),intent(out) :: qsati(ncol,pver)

  real(r8) :: qsatw(ncol,pver)
  real(r8) ::   esl(ncol,pver)
  real(r8) ::   esi(ncol,pver)

  integer :: i,k

  call qsat_water( tair, pair, esl, qsatw )

  do i=1,ncol
  do k=1,pver
     esi(i,k)=svp_ice(tair(i,k))
  end do
  end do

  qsati = qsatw*esi/esl

end subroutine qsat_ice
!----------------------------------

!-----------------------------------------------------------
! supersaturation wrt water (liquid) given as mixing ratio
!
subroutine supersat_q_water( ncol, pver, tair, pair, qv, qssatw )

  use wv_saturation, only: qsat_water

  integer, intent(in)  :: ncol, pver

  real(r8),intent(in)  :: tair(ncol,pver)
  real(r8),intent(in)  :: pair(ncol,pver)
  real(r8),intent(in)  ::   qv(ncol,pver)

  real(r8),intent(out) :: qssatw(ncol,pver)

  real(r8) :: qsatw(ncol,pver)
  real(r8) ::   esl(ncol,pver) !not used

  call qsat_water( tair, pair, esl, qsatw )
  qssatw = qv-qsatw

end subroutine supersat_q_water
!----------------------------------

!-----------------------------------------------------------
! supersaturation wrt ice given as mixing ratio
!
subroutine supersat_q_ice( ncol, pver, tair, pair, qv, qssati )

  use wv_saturation, only: qsat_water, svp_ice

  integer, intent(in)  :: ncol, pver

  real(r8),intent(in)  :: tair(ncol,pver)
  real(r8),intent(in)  :: pair(ncol,pver)
  real(r8),intent(in)  ::   qv(ncol,pver)

  real(r8),intent(out) :: qssati(ncol,pver)

  real(r8) :: qsatw(ncol,pver)
  real(r8) ::   esl(ncol,pver)
  real(r8) ::   esi(ncol,pver)

  integer :: i,k

  call qsat_water( tair, pair, esl, qsatw )

  do i=1,ncol
  do k=1,pver
     esi(i,k)=svp_ice(tair(i,k))
  end do
  end do

  qssati = qv - qsatw*esi/esl


end subroutine supersat_q_ice
!----------------------------------

subroutine relhum_water_percent( ncol, pver, tair, pair, qv,  rhw_percent )

  use wv_saturation, only: qsat_water

  integer, intent(in)  :: ncol, pver

  real(r8),intent(in)  :: tair(ncol,pver)
  real(r8),intent(in)  :: pair(ncol,pver)
  real(r8),intent(in)  ::   qv(ncol,pver)

  real(r8),intent(out) :: rhw_percent(ncol,pver)

  real(r8) :: qsatw(ncol,pver)
  real(r8) ::   esl(ncol,pver) !not used

  call qsat_water( tair, pair, esl, qsatw )
  rhw_percent = qv/qsatw*100._r8

end subroutine relhum_water_percent


!----------------------------------
subroutine relhum_ice_percent( ncol, pver, tair, pair, qv,  rhi_percent )

  use wv_saturation, only: qsat_water, svp_ice

  integer, intent(in)  :: ncol, pver

  real(r8),intent(in)  :: tair(ncol,pver)
  real(r8),intent(in)  :: pair(ncol,pver)
  real(r8),intent(in)  ::   qv(ncol,pver)

  real(r8),intent(out) :: rhi_percent(ncol,pver)

  real(r8) :: qsatw(ncol,pver)
  real(r8) ::   esl(ncol,pver)
  real(r8) ::   esi(ncol,pver)

  integer :: i,k

  call qsat_water( tair, pair, esl, qsatw )

  do i=1,ncol
  do k=1,pver
     esi(i,k)=svp_ice(tair(i,k))
  end do
  end do

  rhi_percent = 100._r8* qv/qsatw* esl/esi

end subroutine relhum_ice_percent

subroutine compute_cape_diags( state, pbuf, pcols, pver, iopt, out1d )
!-------------------------------------------------------------------------------------------
! Purpose: 
! - CAPE, the convecitve available potential energy
! - CAPEp, CAPE assuming fixed environment within one time step but evolving parcle properties
! - dCAPEe, the change in convecitve available potential energy
!   caused by environment change (i.e. assuming fixed parcel property)
! - dCAPEp, the change in convecitve available potential energy caused by parcel change.
!
! History: first version by Hui Wan and Xiaoliang Song, 2021
!-------------------------------------------------------------------------------------------

  use physics_types,  only: physics_state
  use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
  use physconst,      only: cpair, gravit, rair, latvap
  use zm_conv,        only: buoyan_dilute, limcnv

  type(physics_state),intent(in),target:: state
  type(physics_buffer_desc),pointer    :: pbuf(:)
  integer,                  intent(in) :: pver
  integer,                  intent(in) :: pcols
  integer,                  intent(in) :: iopt

  real(r8),                 intent(out) :: out1d(pcols)

  ! local variables used for providing the same input to two calls of subroutine buoyan_dilute

  real(r8) :: pmid_in_hPa(pcols,pver)
  real(r8) :: pint_in_hPa(pcols,pver+1)

  real(r8) ::   zs(pcols)
  real(r8) :: pblt(pcols)
  real(r8) :: zmid_above_sealevel(pcols,pver)

  real(r8),pointer :: tpert(:), pblh(:)

  integer :: idx, kk, lchnk, ncol, msg

  ! variables that distinguish different calls of buoyan_dilute 

  real(r8),pointer ::    qv_new(:,:)  ! new qv   from current state
  real(r8),pointer ::  temp_new(:,:)  ! new temp from current state

  real(r8),pointer ::    qv_old(:,:)  ! old qv   from pbuf
  real(r8),pointer ::  temp_old(:,:)  ! old temp from pbuf

  logical :: l_find_lnch_lvl    ! whether or not to let buoyan_dilute find new launching level

  integer  ::  maxi_new(pcols)  ! index of launching level in new environment,  calculated by buoyan_dilute
  real(r8) ::  q_mx_new(pcols)  ! new specific humidity at new launching level, calculated by buoyan_dilute
  real(r8) ::  t_mx_new(pcols)  ! new temperature       at new launching level, calculated by buoyan_dilute

  integer, pointer :: maxi_old(:)  ! old launching level from pbuf
  real(r8),pointer :: q_mx_old(:)  ! old qv   at launching level from pbuf
  real(r8),pointer :: t_mx_old(:)  ! old temp at launching level from pbuf

  ! variables returned by buoyan_dilute but not needed here

  real(r8) ::   ztp(pcols,pver) ! parcel temperatures.
  real(r8) :: zqstp(pcols,pver) ! grid slice of parcel temp. saturation mixing ratio.
  real(r8) ::   ztl(pcols)      ! parcel temperature at lcl.
  integer  ::  zlcl(pcols)      ! base level index of deep cumulus convection.
  integer  ::  zlel(pcols)      ! index of highest theoretical convective plume.
  integer  ::  zlon(pcols)      ! index of onset level for deep convection.

  real(r8) ::  cape_new_pcl_new_env(pcols) ! cape in new environment (assuming new launching level and parcel properties) 
  real(r8) ::  cape_new_pcl_fixed_env(pcols) ! cape in fixed environment (assuming new launching level and parcel properties) 
  real(r8) ::  cape_new_pcl_old_env(pcols) ! cape in old environment (assuming new launching level and parcel properties) 
  real(r8) ::  cape_old_pcl_new_env(pcols) ! cape in new environment (assuming old launching level and parcel properties) 

  !----------------------------------------------------------------------- 
  ncol  = state%ncol
  lchnk = state%lchnk

  !---------------------------------------------------------------
  ! Temperature and specific humidity in new and old environments
  !---------------------------------------------------------------
  qv_new   => state%q(:,:,1)
  temp_new => state%t

  select case(iopt)

  case(iCAPE_NEW_PCL_FIXED_ENV)
    idx = pbuf_get_index('Q_fixed_4CAPE')  ; call pbuf_get_field( pbuf, idx,   qv_old )
    idx = pbuf_get_index('T_fixed_4CAPE')  ; call pbuf_get_field( pbuf, idx, temp_old )

  case(idCAPEe)
    idx = pbuf_get_index('Q_old_4CAPE')    ; call pbuf_get_field( pbuf, idx,   qv_old )
    idx = pbuf_get_index('T_old_4CAPE')    ; call pbuf_get_field( pbuf, idx, temp_old )

  case(idCAPEp)
    idx = pbuf_get_index('Q_mx_old_4CAPE')  ; call pbuf_get_field( pbuf, idx, q_mx_old )
    idx = pbuf_get_index('T_mx_old_4CAPE')  ; call pbuf_get_field( pbuf, idx, t_mx_old )
    idx = pbuf_get_index('maxi_old_4CAPE')  ; call pbuf_get_field( pbuf, idx, maxi_old )

  end select

  !-----------------------------------
  ! Pressure (in the new environment)
  !-----------------------------------
  pmid_in_hPa(1:ncol,:) = state%pmid(1:ncol,:) * 0.01_r8
  pint_in_hPa(1:ncol,:) = state%pint(1:ncol,:) * 0.01_r8

  !-----------------------------------
  ! Some time-independent quantities 
  !-----------------------------------
  msg = limcnv - 1  ! limcnv is the top interface level limit for convection

  idx = pbuf_get_index('tpert') ; call pbuf_get_field( pbuf, idx, tpert )

  ! Surface elevation (m) is needed to calculate height above sea level (m) 
  ! Note that zm (and zi) stored in state are height above surface. 
  ! The layer midpoint height provided to buoyan_dilute is height above sea level. 

  zs(1:ncol) = state%phis(1:ncol)/gravit

  !------------------------------------------
  ! Height above sea level at layer midpoints
  !------------------------------------------
  do kk = 1,pver
     zmid_above_sealevel(1:ncol,kk) = state%zm(1:ncol,kk)+zs(1:ncol)
  end do

  !--------------------------
  ! layer index for PBL top
  !--------------------------
  idx = pbuf_get_index('pblh')  ; call pbuf_get_field( pbuf, idx, pblh )

  pblt(:) = pver
  do kk = pver-1, msg+1, -1
     where( abs(zmid_above_sealevel(:ncol,kk)-zs(:ncol)-pblh(:ncol))   &
            < (state%zi(:ncol,kk)-state%zi(:ncol,kk+1))*0.5_r8       ) & 
     pblt(:ncol) = kk
  end do

  !------------------------------------------------------------------------
  ! Calculate CAPE using the new state; also return launching level index
  ! and T, qv values at (new) launching level
  !------------------------------------------------------------------------
  l_find_lnch_lvl = .true.
  call buoyan_dilute(lchnk ,ncol,                      &! in
                     qv_new, temp_new,                 &! in  !!
                     pmid_in_hPa, zmid_above_sealevel, &! in
                     pint_in_hPa,                      &! in
                     ztp, zqstp, ztl,                  &! out
                     latvap,                           &! in
                     cape_new_pcl_new_env,             &! out !!
                     pblt,                             &! in
                     zlcl, zlel, zlon,                 &! out
                     maxi_new,                         &! out !!
                     rair, gravit, cpair, msg, tpert,  &! in
                     l_find_lnch_lvl,                  &! in  !!
                     q_mx_new, t_mx_new                )! out !!

  select case (iopt)
  case(iCAPE_NEW_PCL_NEW_ENV)
  !---------------------------------------------------------------------
  ! Output is new CAPE. Copy to out1d and we are done.
  !---------------------------------------------------------------------
     out1d(:ncol) = cape_new_pcl_new_env(:ncol)
 
  case(iCAPE_NEW_PCL_FIXED_ENV) 
  !---------------------------------------------------------------------
  ! Calculate CAPE using 
  !  - a fixed old state (T, qv profiles)
  !  - newly diagnosed launching level and parcel T, qv
  !---------------------------------------------------------------------
    l_find_lnch_lvl = .false.
    call buoyan_dilute(lchnk ,ncol,                      &! in
                       qv_old, temp_old,                 &! in  !!!
                       pmid_in_hPa, zmid_above_sealevel, &! in
                       pint_in_hPa,                      &! in
                       ztp, zqstp, ztl,                  &! out
                       latvap,                           &! in
                       cape_new_pcl_fixed_env,           &! out !!!
                       pblt,                             &! in
                       zlcl, zlel, zlon,                 &! out
                       maxi_new,                         &! in  !!!
                       rair, gravit, cpair, msg, tpert,  &! in
                       l_find_lnch_lvl,                  &! in  !!!
                       q_mx_new, t_mx_new                )! in  !!!

    ! Result to be passed to calling routine is CAPEp (CAPE with fixed env but evolving parcle property)
    out1d = cape_new_pcl_fixed_env
  
    ! No need to anything in pbuf
 
  case(idCAPEe) 
  !---------------------------------------------------------------------
  ! Calculate CAPE using 
  !  - an old state that evolves within a time step 
  !  - newly diagnosed launching level and parcel T, qv
  !---------------------------------------------------------------------
    l_find_lnch_lvl = .false.
    call buoyan_dilute(lchnk ,ncol,                      &! in
                       qv_old, temp_old,                 &! in  !!!
                       pmid_in_hPa, zmid_above_sealevel, &! in
                       pint_in_hPa,                      &! in
                       ztp, zqstp, ztl,                  &! out
                       latvap,                           &! in
                       cape_new_pcl_old_env,             &! out !!!
                       pblt,                             &! in
                       zlcl, zlel, zlon,                 &! out
                       maxi_new,                         &! in  !!!
                       rair, gravit, cpair, msg, tpert,  &! in
                       l_find_lnch_lvl,                  &! in  !!!
                       q_mx_new, t_mx_new                )! in  !!!
 
    ! Result to be passed to calling routine is dCAPEe (CAPE difference caused by environment change).

    out1d(:ncol) = cape_new_pcl_new_env(:ncol) - cape_new_pcl_old_env(:ncol)

    ! Also update the "old" temperature and specific humidity values in pbuf for next call

    temp_old(:ncol,:) = temp_new(:ncol,:)
      qv_old(:ncol,:) =   qv_new(:ncol,:)

  case(idCAPEp) 
  !---------------------------------------------------------------------
  ! Calculate CAPE using 
  !  - new state (T, qv profiles)
  !  - old launching level and parcel T, qv
  !---------------------------------------------------------------------
    l_find_lnch_lvl = .false.
    call buoyan_dilute(lchnk ,ncol,                      &! in
                       qv_new, temp_new,                 &! in  !!!
                       pmid_in_hPa, zmid_above_sealevel, &! in
                       pint_in_hPa,                      &! in
                       ztp, zqstp, ztl,                  &! out
                       latvap,                           &! in
                       cape_old_pcl_new_env,             &! out !!!
                       pblt,                             &! in
                       zlcl, zlel, zlon,                 &! out
                       maxi_old,                         &! in  !!!
                       rair, gravit, cpair, msg, tpert,  &! in
                       l_find_lnch_lvl,                  &! in  !!!
                       q_mx_old, t_mx_old                )! in  !!!

    ! Result to be passed to calling routine is dCAPEp (CAPE difference caused by parcel change).
    out1d(:ncol) = cape_new_pcl_new_env(:ncol) - cape_old_pcl_new_env(:ncol)

    ! Also update "old" parcel properties in pbuf for next call
    t_mx_old(:ncol) = t_mx_new(:ncol)
    q_mx_old(:ncol) = q_mx_new(:ncol)
    maxi_old(:ncol) = maxi_new(:ncol)

  end select

 end subroutine compute_cape_diags
!---------------------------

end module misc_diagnostics
