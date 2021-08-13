module misc_diagnostics

use shr_kind_mod,   only: r8 => shr_kind_r8

implicit none
public

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

subroutine compute_cape( state, pbuf, pcols, pver, cape )
!----------------------------------------------------------------------
! Purpose: compute CAPE (convecitve available potential energy)
!          using subroutine buoyan_dilute from the ZM deep convection
!          parameterization (module file zm_conv.F90)
! History: first version by Hui Wan, 2021-05
!----------------------------------------------------------------------

  use physics_types,  only: physics_state
  use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
  use physconst,      only: cpair, gravit, rair, latvap
  use zm_conv,        only: buoyan_dilute, limcnv

  type(physics_state),intent(in),target:: state
  type(physics_buffer_desc),pointer    :: pbuf(:)
  integer,                  intent(in) :: pver
  integer,                  intent(in) :: pcols
  real(r8),                intent(out) :: cape(pcols)

  ! local variables used for providing input to subroutine buoyan_dilute

  real(r8) :: pmid_in_hPa(pcols,pver)
  real(r8) :: pint_in_hPa(pcols,pver+1)

  real(r8),pointer ::    qv(:,:)
  real(r8),pointer ::  temp(:,:)

  real(r8) ::   zs(pcols)
  real(r8) :: pblt(pcols)
  real(r8) :: zmid_above_sealevel(pcols,pver)

  real(r8),pointer :: tpert(:), pblh(:)

  integer :: idx, kk, lchnk, ncol, msg

  ! variables returned by buoyan_dilute but not needed here

  real(r8) ::   ztp(pcols,pver) ! parcel temperatures.
  real(r8) :: zqstp(pcols,pver) ! grid slice of parcel temp. saturation mixing ratio.
  real(r8) ::   ztl(pcols)      ! parcel temperature at lcl.
  integer  ::  zlcl(pcols)      ! base level index of deep cumulus convection.
  integer  ::  zlel(pcols)      ! index of highest theoretical convective plume.
  integer  ::  zlon(pcols)      ! index of onset level for deep convection.
  integer  :: zmaxi(pcols)      ! index of level with largest moist static energy.

  !----------------------------------------------------------------------- 
  ncol  = state%ncol
  lchnk = state%lchnk

  msg = limcnv - 1  ! limcnv is the top interface level limit for convection

  idx = pbuf_get_index('tpert') ; call pbuf_get_field( pbuf, idx, tpert )

  pmid_in_hPa(1:ncol,:) = state%pmid(1:ncol,:) * 0.01_r8
  pint_in_hPa(1:ncol,:) = state%pint(1:ncol,:) * 0.01_r8

  qv   => state%q(:,:,1)
  temp => state%t

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

  !----------------
  ! Calculate CAPE
  !----------------
  call buoyan_dilute(lchnk ,ncol, qv, temp,            &! in
                     pmid_in_hPa, zmid_above_sealevel, &! in
                     pint_in_hPa,                      &! in
                     ztp, zqstp, ztl,                  &! out
                     latvap, cape, pblt,               &! in, out, in
                     zlcl, zlel, zlon, zmaxi,          &! out
                     rair, gravit, cpair, msg, tpert   )! in

 end subroutine compute_cape
!---------------------------

subroutine compute_dcape_env( state, pbuf, pcols, pver, dcape_env )
!-------------------------------------------------------------------------------------------
! Purpose: compute dCAPE_env, the change in convecitve available potential energy
!          caused by environment change (i.e. assuming fixed parcel property)
!          using subroutine buoyan_dilute from the ZM deep convection
!          parameterization (module file zm_conv.F90)
! History: first version by Hui Wan and Xiaoliang Song, 2021-08
!-------------------------------------------------------------------------------------------

  use physics_types,  only: physics_state
  use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
  use physconst,      only: cpair, gravit, rair, latvap
  use zm_conv,        only: buoyan_dilute, limcnv

  type(physics_state),intent(in),target:: state
  type(physics_buffer_desc),pointer    :: pbuf(:)
  integer,                  intent(in) :: pver
  integer,                  intent(in) :: pcols
  real(r8),                 intent(out) :: dcape_env(pcols)

  ! local variables used for providing the same input to two calls of subroutine buoyan_dilute

  real(r8) :: pmid_in_hPa(pcols,pver)
  real(r8) :: pint_in_hPa(pcols,pver+1)

  real(r8) ::   zs(pcols)
  real(r8) :: pblt(pcols)
  real(r8) :: zmid_above_sealevel(pcols,pver)

  real(r8),pointer :: tpert(:), pblh(:)

  integer :: idx, kk, lchnk, ncol, msg

  ! variables that distinguish the two calls of buoyan_dilute 

  real(r8),pointer ::    qv_new(:,:)  ! new qv   from current state
  real(r8),pointer ::  temp_new(:,:)  ! new temp from current state

  real(r8),pointer ::    qv_old(:,:)  ! old qv   from pbuf
  real(r8),pointer ::  temp_old(:,:)  ! old temp from pbuf

  logical :: l_find_lnch_lvl    ! whether or not to let buoyan_dilute find new launching level

  integer  :: zmaxi_new(pcols)  ! index of launching level in new environment 
  real(r8) :: zq_mx_new(pcols)  ! new specific humidity at new launching level
  real(r8) :: zt_mx_new(pcols)  ! new temperature       at new launching level

  ! variables returned by buoyan_dilute but not needed here

  real(r8) ::   ztp(pcols,pver) ! parcel temperatures.
  real(r8) :: zqstp(pcols,pver) ! grid slice of parcel temp. saturation mixing ratio.
  real(r8) ::   ztl(pcols)      ! parcel temperature at lcl.
  integer  ::  zlcl(pcols)      ! base level index of deep cumulus convection.
  integer  ::  zlel(pcols)      ! index of highest theoretical convective plume.
  integer  ::  zlon(pcols)      ! index of onset level for deep convection.

  real(r8) ::  cape_new(pcols)        ! cape in new environment (assuming new launching level and parcel properties) 
  real(r8) ::  cape_in_old_env(pcols) ! cape in old environment (assuming new launching level and parcel properties) 

  !----------------------------------------------------------------------- 
  ncol  = state%ncol
  lchnk = state%lchnk

  !--------------------------------------------------------------
  ! Temperature and specific humidity in new and old environment
  !--------------------------------------------------------------
  qv_new   => state%q(:,:,1)
  temp_new => state%t

  idx = pbuf_get_index('Q_old_4CAPE')  ; call pbuf_get_field( pbuf, idx,   qv_old )
  idx = pbuf_get_index('T_old_4CAPE')  ; call pbuf_get_field( pbuf, idx, temp_old )

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
                     cape_new,                         &! out !!
                     pblt,                             &! in
                     zlcl, zlel, zlon,                 &! out
                     zmaxi_new,                        &! out !!
                     rair, gravit, cpair, msg, tpert,  &! in
                     l_find_lnch_lvl,                  &! in  !!
                     zq_mx_new, zt_mx_new              )! out !!

  !---------------------------------------------------------------------
  ! Calculate CAPE using 
  !  - an old state (T, qv profiles)
  !  - newly diagnosed launching level and parcel T, qv
  !---------------------------------------------------------------------
  l_find_lnch_lvl = .false.
  call buoyan_dilute(lchnk ,ncol,                      &! in
                     qv_old, temp_old,                 &! in  !!!
                     pmid_in_hPa, zmid_above_sealevel, &! in
                     pint_in_hPa,                      &! in
                     ztp, zqstp, ztl,                  &! out
                     latvap,                           &! in
                     cape_in_old_env,                  &! out !!!
                     pblt,                             &! in
                     zlcl, zlel, zlon,                 &! out
                     zmaxi_new,                        &! in  !!!
                     rair, gravit, cpair, msg, tpert,  &! in
                     l_find_lnch_lvl,                  &! in  !!!
                     zq_mx_new, zt_mx_new              )! in  !!!

  !---------------------------------------------------------------------
  ! Calculate CAPE difference caused by environment change
  !---------------------------------------------------------------------
  dcape_env(:ncol) = cape_new(:ncol) - cape_in_old_env(:ncol)

  !---------------------------------------------------------------------
  ! Update the "old" temperature and specific humidity values in pbuf
  ! for next call
  !---------------------------------------------------------------------
  temp_old(:ncol) = temp_new(:ncol)
    qv_old(:ncol) =   qv_new(:ncol)

 end subroutine compute_dcape_env
!---------------------------

subroutine compute_cape_pcl( state, pbuf, pcols, pver, cape_pcl )
!-------------------------------------------------------------------------------------------
! Purpose: compute CAPE_pcl, the convecitve available potential energy
!          assuming fixed environment but evolving parcel property.
! History: first version by Hui Wan and Xiaoliang Song, 2021-08
!-------------------------------------------------------------------------------------------

  use physics_types,  only: physics_state
  use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
  use physconst,      only: cpair, gravit, rair, latvap
  use zm_conv,        only: buoyan_dilute, limcnv

  type(physics_state),intent(in),target:: state
  type(physics_buffer_desc),pointer    :: pbuf(:)
  integer,                  intent(in) :: pver
  integer,                  intent(in) :: pcols
  real(r8),                 intent(out) :: cape_pcl(pcols)

  ! local variables used for providing the same input to two calls of subroutine buoyan_dilute

  real(r8) :: pmid_in_hPa(pcols,pver)
  real(r8) :: pint_in_hPa(pcols,pver+1)

  real(r8) ::   zs(pcols)
  real(r8) :: pblt(pcols)
  real(r8) :: zmid_above_sealevel(pcols,pver)

  real(r8),pointer :: tpert(:), pblh(:)

  integer :: idx, kk, lchnk, ncol, msg

  ! variables that distinguish the two calls of buoyan_dilute 

  real(r8),pointer ::    qv_new(:,:)  ! new qv   from current state
  real(r8),pointer ::  temp_new(:,:)  ! new temp from current state

  real(r8),pointer ::    qv_old(:,:)  ! old qv   from pbuf
  real(r8),pointer ::  temp_old(:,:)  ! old temp from pbuf

  logical :: l_find_lnch_lvl    ! whether or not to let buoyan_dilute find new launching level

  integer  :: zmaxi_new(pcols)  ! index of launching level in new environment 
  real(r8) :: zq_mx_new(pcols)  ! new specific humidity at new launching level
  real(r8) :: zt_mx_new(pcols)  ! new temperature       at new launching level

  ! variables returned by buoyan_dilute but not needed here

  real(r8) ::   ztp(pcols,pver) ! parcel temperatures.
  real(r8) :: zqstp(pcols,pver) ! grid slice of parcel temp. saturation mixing ratio.
  real(r8) ::   ztl(pcols)      ! parcel temperature at lcl.
  integer  ::  zlcl(pcols)      ! base level index of deep cumulus convection.
  integer  ::  zlel(pcols)      ! index of highest theoretical convective plume.
  integer  ::  zlon(pcols)      ! index of onset level for deep convection.

  real(r8) ::  cape_new(pcols)  ! cape in new environment (assuming new launching level and parcel properties) 

  !----------------------------------------------------------------------- 
  ncol  = state%ncol
  lchnk = state%lchnk

  !--------------------------------------------------------------
  ! Temperature and specific humidity in new and old environment
  !--------------------------------------------------------------
  qv_new   => state%q(:,:,1)
  temp_new => state%t

  idx = pbuf_get_index('Q_fixed_4CAPE')  ; call pbuf_get_field( pbuf, idx,   qv_old )
  idx = pbuf_get_index('T_fixed_4CAPE')  ; call pbuf_get_field( pbuf, idx, temp_old )

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
                     cape_new,                         &! out !!
                     pblt,                             &! in
                     zlcl, zlel, zlon,                 &! out
                     zmaxi_new,                        &! out !!
                     rair, gravit, cpair, msg, tpert,  &! in
                     l_find_lnch_lvl,                  &! in  !!
                     zq_mx_new, zt_mx_new              )! out !!

  !---------------------------------------------------------------------
  ! Calculate CAPE using 
  !  - an old state (T, qv profiles)
  !  - newly diagnosed launching level and parcel T, qv
  !---------------------------------------------------------------------
  l_find_lnch_lvl = .false.
  call buoyan_dilute(lchnk ,ncol,                      &! in
                     qv_old, temp_old,                 &! in  !!!
                     pmid_in_hPa, zmid_above_sealevel, &! in
                     pint_in_hPa,                      &! in
                     ztp, zqstp, ztl,                  &! out
                     latvap,                           &! in
                     cape_pcl,                         &! out !!!
                     pblt,                             &! in
                     zlcl, zlel, zlon,                 &! out
                     zmaxi_new,                        &! in  !!!
                     rair, gravit, cpair, msg, tpert,  &! in
                     l_find_lnch_lvl,                  &! in  !!!
                     zq_mx_new, zt_mx_new              )! in  !!!

 end subroutine compute_cape_pcl
!--------------------------------


end module misc_diagnostics
