module cloud_rad_props

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

use shr_kind_mod,     only: r8 => shr_kind_r8
use radconstants,     only: nswbands, nlwbands
use rad_constituents, only: iceopticsfile, liqopticsfile
use interpolate_data, only: interp_type, lininterp_init, lininterp, &
                            extrap_method_bndry, lininterp_finish
use cam_logfile,      only: iulog
use cam_abortutils,   only: endrun

implicit none
private
save

public :: &
   cloud_rad_props_init,     &
   mitchell_ice_optics_sw,   & ! return Mitchell SW ice radiative properties
   mitchell_ice_optics_lw,   & ! Mitchell LW ice rad props
   conley_liquid_optics_sw,  & ! return Conley SW rad props
   conley_liquid_optics_lw     ! return Conley LW rad props

integer :: nmu, nlambda
real(r8), allocatable :: g_mu(:)           ! mu samples on grid
real(r8), allocatable :: g_lambda(:,:)     ! lambda scale samples on grid
real(r8), allocatable :: ext_sw_liq(:,:,:)
real(r8), allocatable :: ssa_sw_liq(:,:,:)
real(r8), allocatable :: asm_sw_liq(:,:,:)
real(r8), allocatable :: abs_lw_liq(:,:,:)

integer :: n_g_d
real(r8), allocatable :: g_d_eff(:)        ! radiative effective diameter samples on grid
real(r8), allocatable :: ext_sw_ice(:,:)
real(r8), allocatable :: ssa_sw_ice(:,:)
real(r8), allocatable :: asm_sw_ice(:,:)
real(r8), allocatable :: abs_lw_ice(:,:)

!==============================================================================
contains
!==============================================================================

subroutine cloud_rad_props_init()

   use netcdf
   use spmd_utils,     only: masterproc
   use ioFileMod,      only: getfil
   use error_messages, only: handle_ncerr
#if ( defined SPMD )
   use mpishorthand
#endif

   character(len=256) :: locfn
   integer :: ncid, dimid, varid, f_nlwbands, f_nswbands, ierr
   integer :: vdimids(NF90_MAX_VAR_DIMS)

   ! read liquid cloud optics
   if(masterproc) then
      call getfil( trim(liqopticsfile), locfn, 0)
      call handle_ncerr( nf90_open(locfn, NF90_NOWRITE, ncid), 'liquid optics file missing')
      write(iulog,*)' reading liquid cloud optics from file ',locfn

      call handle_ncerr(nf90_inq_dimid( ncid, 'lw_band', dimid), 'getting lw_band dim')
      call handle_ncerr(nf90_inquire_dimension( ncid, dimid, len=f_nlwbands), 'getting n lw bands')
      if (f_nlwbands /= nlwbands) call endrun('number of lw bands does not match')

      call handle_ncerr(nf90_inq_dimid( ncid, 'sw_band', dimid), 'getting sw_band_dim')
      call handle_ncerr(nf90_inquire_dimension( ncid, dimid, len=f_nswbands), 'getting n sw bands')
      if (f_nswbands /= nswbands) call endrun('number of sw bands does not match')

      call handle_ncerr(nf90_inq_dimid( ncid, 'mu', dimid), 'getting mu dim')
      call handle_ncerr(nf90_inquire_dimension( ncid, dimid, len=nmu), 'getting n mu samples')

      call handle_ncerr(nf90_inq_dimid( ncid, 'lambda_scale', dimid), 'getting lambda dim')
      call handle_ncerr(nf90_inquire_dimension( ncid, dimid, len=nlambda), 'getting n lambda samples')
   endif ! if (masterproc)

#if ( defined SPMD )
   call mpibcast(nmu, 1, mpiint, 0, mpicom, ierr)
   call mpibcast(nlambda, 1, mpiint, 0, mpicom, ierr)
#endif

   if (.not.allocated(g_mu)) allocate(g_mu(nmu))
   if (.not.allocated(g_lambda)) allocate(g_lambda(nmu,nlambda))
   if (.not.allocated(ext_sw_liq)) allocate(ext_sw_liq(nmu,nlambda,nswbands) )
   if (.not.allocated(ssa_sw_liq)) allocate(ssa_sw_liq(nmu,nlambda,nswbands))
   if (.not.allocated(asm_sw_liq)) allocate(asm_sw_liq(nmu,nlambda,nswbands))
   if (.not.allocated(abs_lw_liq)) allocate(abs_lw_liq(nmu,nlambda,nlwbands))

   if(masterproc) then
   call handle_ncerr( nf90_inq_varid(ncid, 'mu', varid),&
      'cloud optics mu get')
   call handle_ncerr( nf90_get_var(ncid, varid, g_mu),&
      'read cloud optics mu values')

   call handle_ncerr( nf90_inq_varid(ncid, 'lambda', varid),&
      'cloud optics lambda get')
   call handle_ncerr( nf90_get_var(ncid, varid, g_lambda),&
      'read cloud optics lambda values')

   call handle_ncerr( nf90_inq_varid(ncid, 'k_ext_sw', varid),&
      'cloud optics ext_sw_liq get')
   call handle_ncerr( nf90_get_var(ncid, varid, ext_sw_liq),&
      'read cloud optics ext_sw_liq values')

   call handle_ncerr( nf90_inq_varid(ncid, 'ssa_sw', varid),&
      'cloud optics ssa_sw_liq get')
   call handle_ncerr( nf90_get_var(ncid, varid, ssa_sw_liq),&
      'read cloud optics ssa_sw_liq values')

   call handle_ncerr( nf90_inq_varid(ncid, 'asm_sw', varid),&
      'cloud optics asm_sw_liq get')
   call handle_ncerr( nf90_get_var(ncid, varid, asm_sw_liq),&
      'read cloud optics asm_sw_liq values')

   call handle_ncerr( nf90_inq_varid(ncid, 'k_abs_lw', varid),&
      'cloud optics abs_lw_liq get')
   call handle_ncerr( nf90_get_var(ncid, varid, abs_lw_liq),&
      'read cloud optics abs_lw_liq values')

   call handle_ncerr( nf90_close(ncid), 'liquid optics file missing')
   endif ! if masterproc

#if ( defined SPMD )
    call mpibcast(g_mu, nmu, mpir8, 0, mpicom, ierr)
    call mpibcast(g_lambda, nmu*nlambda, mpir8, 0, mpicom, ierr)
    call mpibcast(ext_sw_liq, nmu*nlambda*nswbands, mpir8, 0, mpicom, ierr)
    call mpibcast(ssa_sw_liq, nmu*nlambda*nswbands, mpir8, 0, mpicom, ierr)
    call mpibcast(asm_sw_liq, nmu*nlambda*nswbands, mpir8, 0, mpicom, ierr)
    call mpibcast(abs_lw_liq, nmu*nlambda*nlwbands, mpir8, 0, mpicom, ierr)
#endif
   ! I forgot to convert kext from m^2/Volume to m^2/Kg
   ext_sw_liq = ext_sw_liq / 0.9970449e3_r8 
   abs_lw_liq = abs_lw_liq / 0.9970449e3_r8 

   ! read ice cloud optics
   if(masterproc) then
      call getfil( trim(iceopticsfile), locfn, 0)
      call handle_ncerr( nf90_open(locfn, NF90_NOWRITE, ncid), 'ice optics file missing')
      write(iulog,*)' reading ice cloud optics from file ',locfn

      call handle_ncerr(nf90_inq_dimid( ncid, 'lw_band', dimid), 'getting lw_band dim')
      call handle_ncerr(nf90_inquire_dimension( ncid, dimid, len=f_nlwbands), 'getting n lw bands')
      if (f_nlwbands /= nlwbands) call endrun('number of lw bands does not match')

      call handle_ncerr(nf90_inq_dimid( ncid, 'sw_band', dimid), 'getting sw_band_dim')
      call handle_ncerr(nf90_inquire_dimension( ncid, dimid, len=f_nswbands), 'getting n sw bands')
      if (f_nswbands /= nswbands) call endrun('number of sw bands does not match')

      call handle_ncerr(nf90_inq_dimid( ncid, 'd_eff', dimid), 'getting deff dim')
      call handle_ncerr(nf90_inquire_dimension( ncid, dimid, len=n_g_d), 'getting n deff samples')

   endif ! if (masterproc)

#if ( defined SPMD )
   call mpibcast(n_g_d, 1, mpiint, 0, mpicom, ierr)
#endif

   if (.not.allocated(g_d_eff)) allocate(g_d_eff(n_g_d))
   if (.not.allocated(ext_sw_ice)) allocate(ext_sw_ice(n_g_d,nswbands))
   if (.not.allocated(ssa_sw_ice)) allocate(ssa_sw_ice(n_g_d,nswbands))
   if (.not.allocated(asm_sw_ice)) allocate(asm_sw_ice(n_g_d,nswbands))
   if (.not.allocated(abs_lw_ice)) allocate(abs_lw_ice(n_g_d,nlwbands))

   if(masterproc) then
      call handle_ncerr( nf90_inq_varid(ncid, 'd_eff', varid),&
         'cloud optics deff get')
      call handle_ncerr( nf90_get_var(ncid, varid, g_d_eff),&
         'read cloud optics deff values')

      call handle_ncerr( nf90_inq_varid(ncid, 'sw_ext', varid),&
         'cloud optics ext_sw_ice get')
      call handle_ncerr(nf90_inquire_variable ( ncid, varid, dimids=vdimids),&
          'checking dimensions of ext_sw_ice')
      call handle_ncerr(nf90_inquire_dimension( ncid, vdimids(1)),&
          'getting first dimension sw_ext')
      call handle_ncerr(nf90_inquire_dimension( ncid, vdimids(2)),&
          'getting first dimension sw_ext')
      call handle_ncerr( nf90_get_var(ncid, varid, ext_sw_ice),&
         'read cloud optics ext_sw_ice values')

      call handle_ncerr( nf90_inq_varid(ncid, 'sw_ssa', varid),&
         'cloud optics ssa_sw_ice get')
      call handle_ncerr( nf90_get_var(ncid, varid, ssa_sw_ice),&
         'read cloud optics ssa_sw_ice values')

      call handle_ncerr( nf90_inq_varid(ncid, 'sw_asm', varid),&
         'cloud optics asm_sw_ice get')
      call handle_ncerr( nf90_get_var(ncid, varid, asm_sw_ice),&
         'read cloud optics asm_sw_ice values')

      call handle_ncerr( nf90_inq_varid(ncid, 'lw_abs', varid),&
         'cloud optics abs_lw_ice get')
      call handle_ncerr( nf90_get_var(ncid, varid, abs_lw_ice),&
         'read cloud optics abs_lw_ice values')

      call handle_ncerr( nf90_close(ncid), 'ice optics file missing')

   endif ! if masterproc
#if ( defined SPMD )
   call mpibcast(g_d_eff, n_g_d, mpir8, 0, mpicom, ierr)
   call mpibcast(ext_sw_ice, n_g_d*nswbands, mpir8, 0, mpicom, ierr)
   call mpibcast(ssa_sw_ice, n_g_d*nswbands, mpir8, 0, mpicom, ierr)
   call mpibcast(asm_sw_ice, n_g_d*nswbands, mpir8, 0, mpicom, ierr)
   call mpibcast(abs_lw_ice, n_g_d*nlwbands, mpir8, 0, mpicom, ierr)
#endif

   return

end subroutine cloud_rad_props_init

!==============================================================================

subroutine conley_liquid_optics_sw(ncol, nlev, lamc, pgam, iclwpth, tau, ssa, asm)
   integer, intent(in) :: ncol, nlev
   real(r8), intent(in), dimension(ncol,nlev) :: lamc, pgam, iclwpth
   real(r8), intent(out), dimension(nswbands,ncol,nlev) :: tau  ! extinction optical depth
   real(r8), intent(out), dimension(nswbands,ncol,nlev) :: ssa  ! single scattering albedo * tau
   real(r8), intent(out), dimension(nswbands,ncol,nlev) :: asm  ! asymetry parameter * tau * w
   integer :: i, k, swband
   type(interp_type) :: mu_wgts
   type(interp_type) :: lambda_wgts

   do k = 1,nlev
      do i = 1,ncol
         if (lamc(i,k) > 0._r8 .and. iclwpth(i,k) >= 1.e-80_r8) then
            call get_mu_lambda_weights(lamc(i,k), pgam(i,k), mu_wgts, lambda_wgts)
            do swband = 1,nswbands
               call lininterp(ext_sw_liq(:,:,swband), nmu, nlambda, &
                    tau(swband:swband,i,k), 1, mu_wgts, lambda_wgts)
               call lininterp(ssa_sw_liq(:,:,swband), nmu, nlambda, &
                    ssa(swband:swband,i,k), 1, mu_wgts, lambda_wgts)
               call lininterp(asm_sw_liq(:,:,swband), nmu, nlambda, &
                    asm(swband:swband,i,k), 1, mu_wgts, lambda_wgts)
               ! compute radiative properties
               tau(swband,i,k) = iclwpth(i,k) * tau(swband,i,k)
            enddo
            call lininterp_finish(mu_wgts)
            call lininterp_finish(lambda_wgts)
         else
            tau(1:nswbands,i,k) = 0._r8
            ssa(1:nswbands,i,k) = 0._r8
            asm(1:nswbands,i,k) = 0._r8
         endif
      enddo
   enddo
end subroutine conley_liquid_optics_sw

!==============================================================================

subroutine conley_liquid_optics_lw(ncol, nlev, lamc, pgam, iclwpth, abs_od)
   integer, intent(in) :: ncol, nlev
   real(r8), intent(in), dimension(ncol,nlev) :: lamc, pgam, iclwpth
   real(r8), intent(out), dimension(nlwbands,ncol,nlev) :: abs_od
   integer i, k, lwband
   type(interp_type) :: mu_wgts
   type(interp_type) :: lambda_wgts

   abs_od = 0._r8
   do k = 1,nlev
      do i = 1,ncol
         if(lamc(i,k) > 0._r8 .and. iclwpth(i,k) >= 1.e-80_r8) then
            call get_mu_lambda_weights(lamc(i,k), pgam(i,k), mu_wgts, lambda_wgts)
            do lwband = 1,nlwbands
               call lininterp(abs_lw_liq(:,:,lwband), nmu, nlambda, &
                    abs_od(lwband:lwband,i,k), 1, mu_wgts, lambda_wgts)
               abs_od(lwband,i,k) = iclwpth(i,k) * abs_od(lwband,i,k)
            enddo
            call lininterp_finish(mu_wgts)
            call lininterp_finish(lambda_wgts)
         else
            abs_od(1:nlwbands,i,k) = 0._r8
         endif
      enddo
   enddo
end subroutine conley_liquid_optics_lw

!==============================================================================

subroutine mitchell_ice_optics_sw(ncol, nlev, iciwpth, dei, tau, ssa, asm)

  integer, intent(in) :: ncol, nlev
  real(r8), intent(in) :: iciwpth(ncol,nlev)
  real(r8), intent(in) :: dei(ncol,nlev)
  real(r8),intent(out) :: tau(nswbands,ncol,nlev) ! extinction optical depth
  real(r8),intent(out) :: ssa(nswbands,ncol,nlev) ! single scattering albedo
  real(r8),intent(out) :: asm(nswbands,ncol,nlev) ! assymetry parameter

  type(interp_type) :: dei_wgts
  integer :: i, k, swband

  do k = 1,nlev
     do i = 1,ncol
        if( iciwpth(i,k) < 1.e-80_r8 .or. dei(i,k) == 0._r8) then
           ! if ice water path is too small, OD := 0
           tau(:,i,k) = 0._r8
           ssa(:,i,k) = 0._r8
           asm(:,i,k) = 0._r8
        else
           ! for each cell interpolate to find weights in g_d_eff grid.
           call lininterp_init(g_d_eff, n_g_d, dei(i:i,k), 1, &
                extrap_method_bndry, dei_wgts)
           ! interpolate into grid and extract radiative properties
           do swband = 1, nswbands
              call lininterp(ext_sw_ice(:,swband), n_g_d, &
                   tau(swband:swband,i,k), 1, dei_wgts)
              call lininterp(ssa_sw_ice(:,swband), n_g_d, &
                   ssa(swband:swband,i,k), 1, dei_wgts)
              call lininterp(asm_sw_ice(:,swband), n_g_d, &
                   asm(swband:swband,i,k), 1, dei_wgts)
              ! Convert extinction to optical depth
              tau(swband,i,k) = iciwpth(i,k) * tau(swband,i,k)
           end do
           call lininterp_finish(dei_wgts)
        endif
     enddo
  enddo

end subroutine mitchell_ice_optics_sw

!==============================================================================

subroutine mitchell_ice_optics_lw(ncol, nlev, iciwpth, dei, abs_od)

  integer, intent(in) :: ncol, nlev
  real(r8), intent(in) :: iciwpth(ncol,nlev)
  real(r8), intent(in) :: dei(ncol,nlev)
  real(r8),intent(out) :: abs_od(nlwbands,ncol,nlev)

  type(interp_type) :: dei_wgts
  integer :: i, k, lwband
  real(r8) :: absor(nlwbands)

  do k = 1,nlev
     do i = 1,ncol
        ! if ice water path is too small, OD := 0
        if( iciwpth(i,k) < 1.e-80_r8 .or. dei(i,k) == 0._r8) then
           abs_od (:,i,k) = 0._r8
        else
           ! for each cell interpolate to find weights in g_d_eff grid.
           call lininterp_init(g_d_eff, n_g_d, dei(i:i,k), 1, &
                extrap_method_bndry, dei_wgts)
           ! interpolate into grid and extract radiative properties
           do lwband = 1, nlwbands
              call lininterp(abs_lw_ice(:,lwband), n_g_d, &
                   absor(lwband:lwband), 1, dei_wgts)
           enddo
           abs_od(:,i,k) = iciwpth(i,k) * absor
           call lininterp_finish(dei_wgts)
        endif
     enddo
  enddo

end subroutine mitchell_ice_optics_lw

!==============================================================================

subroutine get_mu_lambda_weights(lamc, pgam, mu_wgts, lambda_wgts)
  real(r8), intent(in) :: lamc   ! prognosed value of lambda for cloud
  real(r8), intent(in) :: pgam   ! prognosed value of mu for cloud
  ! Output interpolation weights. Caller is responsible for freeing these.
  type(interp_type), intent(out) :: mu_wgts
  type(interp_type), intent(out) :: lambda_wgts

  integer :: ilambda
  real(r8) :: g_lambda_interp(nlambda)

  ! Make interpolation weights for mu.
  ! (Put pgam in a temporary array for this purpose.)
  call lininterp_init(g_mu, nmu, [pgam], 1, extrap_method_bndry, mu_wgts)

  ! Use mu weights to interpolate to a row in the lambda table.
  do ilambda = 1, nlambda
     call lininterp(g_lambda(:,ilambda), nmu, &
          g_lambda_interp(ilambda:ilambda), 1, mu_wgts)
  end do

  ! Make interpolation weights for lambda.
  call lininterp_init(g_lambda_interp, nlambda, [lamc], 1, &
       extrap_method_bndry, lambda_wgts)

end subroutine get_mu_lambda_weights

!==============================================================================

end module cloud_rad_props
