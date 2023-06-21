
module mo_usrrxt

  use shr_kind_mod, only : r8 => shr_kind_r8
  use cam_logfile,  only : iulog
  use ppgrid,       only : pver, pcols
  use cam_abortutils,   only : endrun

  implicit none

  private
  public :: usrrxt, usrrxt_inti

  save

  integer :: usr_O_O2_ndx
  integer :: usr_HO2_HO2_ndx 
  integer :: usr_N2O5_M_ndx
  integer :: usr_HNO3_OH_ndx 
  integer :: usr_HO2NO2_M_ndx
  integer :: usr_N2O5_aer_ndx
  integer :: usr_NO3_aer_ndx
  integer :: usr_NO2_aer_ndx
  integer :: usr_CO_OH_b_ndx
  integer :: usr_PAN_M_ndx
  integer :: usr_CH3COCH3_OH_ndx
  integer :: usr_MCO3_NO2_ndx
  integer :: usr_MPAN_M_ndx
  integer :: usr_XOOH_OH_ndx 
  integer :: usr_SO2_OH_ndx
  integer :: usr_DMS_OH_ndx
  integer :: usr_HO2_aer_ndx
  
  integer :: tag_NO2_NO3_ndx
  integer :: tag_NO2_OH_ndx
  integer :: tag_NO2_HO2_ndx
  integer :: tag_C2H4_OH_ndx
  integer :: tag_C3H6_OH_ndx
  integer :: tag_CH3CO3_NO2_ndx

  integer :: h2o_ndx

  real(r8), parameter :: t0     = 300._r8                ! K

contains

  subroutine usrrxt_inti
    !-----------------------------------------------------------------
    !        ... intialize the user reaction constants module
    !-----------------------------------------------------------------

    use mo_chem_utls,   only : get_rxt_ndx, get_spc_ndx
    use spmd_utils,     only : masterproc
    use physics_buffer, only : pbuf_get_index

    implicit none
!
! full tropospheric chemistry
!

    usr_O_O2_ndx         = get_rxt_ndx( 'usr_O_O2' )
    usr_HO2_HO2_ndx      = get_rxt_ndx( 'usr_HO2_HO2' )
    usr_N2O5_M_ndx       = get_rxt_ndx( 'usr_N2O5_M' )
    usr_HNO3_OH_ndx      = get_rxt_ndx( 'usr_HNO3_OH' )
    usr_HO2NO2_M_ndx     = get_rxt_ndx( 'usr_HO2NO2_M' )
    usr_N2O5_aer_ndx     = get_rxt_ndx( 'usr_N2O5_aer' )
    usr_NO3_aer_ndx      = get_rxt_ndx( 'usr_NO3_aer' )
    usr_NO2_aer_ndx      = get_rxt_ndx( 'usr_NO2_aer' )
    usr_CO_OH_b_ndx      = get_rxt_ndx( 'usr_CO_OH_b' )
    usr_PAN_M_ndx        = get_rxt_ndx( 'usr_PAN_M' )
    usr_CH3COCH3_OH_ndx  = get_rxt_ndx( 'usr_CH3COCH3_OH' )
    usr_MCO3_NO2_ndx     = get_rxt_ndx( 'usr_MCO3_NO2' )
    usr_MPAN_M_ndx       = get_rxt_ndx( 'usr_MPAN_M' )
    usr_XOOH_OH_ndx      = get_rxt_ndx( 'usr_XOOH_OH' )
    usr_SO2_OH_ndx       = get_rxt_ndx( 'usr_SO2_OH' )
    usr_DMS_OH_ndx       = get_rxt_ndx( 'usr_DMS_OH' )
    usr_HO2_aer_ndx      = get_rxt_ndx( 'usr_HO2_aer' )
 !
    tag_NO2_NO3_ndx      = get_rxt_ndx( 'tag_NO2_NO3' )
    tag_NO2_OH_ndx       = get_rxt_ndx( 'tag_NO2_OH' )
    tag_NO2_HO2_ndx      = get_rxt_ndx( 'tag_NO2_HO2' )
    tag_C2H4_OH_ndx      = get_rxt_ndx( 'tag_C2H4_OH' )
    tag_C3H6_OH_ndx      = get_rxt_ndx( 'tag_C3H6_OH' )
    tag_CH3CO3_NO2_ndx   = get_rxt_ndx( 'tag_CH3CO3_NO2' )     

    h2o_ndx    = get_spc_ndx( 'H2O' )

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'usrrxt_inti: diagnostics '
       write(iulog,'(10i5)') usr_O_O2_ndx,usr_HO2_HO2_ndx,tag_NO2_NO3_ndx,usr_N2O5_M_ndx,tag_NO2_OH_ndx,usr_HNO3_OH_ndx &
                            ,tag_NO2_HO2_ndx,usr_HO2NO2_M_ndx,usr_N2O5_aer_ndx,usr_NO3_aer_ndx,usr_NO2_aer_ndx &
                            ,usr_CO_OH_b_ndx,tag_C2H4_OH_ndx,tag_C3H6_OH_ndx,tag_CH3CO3_NO2_ndx,usr_PAN_M_ndx,usr_CH3COCH3_OH_ndx &
                            ,usr_MCO3_NO2_ndx,usr_MPAN_M_ndx,usr_XOOH_OH_ndx,usr_SO2_OH_ndx,usr_DMS_OH_ndx,usr_HO2_aer_ndx
    endif

  end subroutine usrrxt_inti

  subroutine usrrxt( rxt, & ! inout
                    temp, invariants, mtot, ncol )  ! in

!-----------------------------------------------------------------
!        ... set the user specified reaction rates
!-----------------------------------------------------------------
    
    use chem_mods,     only : nfs, rxntot
    use mo_setinv,     only : inv_h2o_ndx=>h2o_ndx
    implicit none

!-----------------------------------------------------------------
!        ... dummy arguments
!-----------------------------------------------------------------
    integer, intent(in)     :: ncol                       ! number of columns
    real(r8), intent(in)    :: temp(pcols,pver)           ! neutral temperature [K]
    real(r8), intent(in)    :: mtot(ncol,pver)            ! total atm density [/cm^3]
    real(r8), intent(in)    :: invariants(ncol,pver,nfs)  ! invariants density [/cm^3]
    real(r8), intent(inout) :: rxt(ncol,pver,rxntot)      ! gas phase rates
      
!-----------------------------------------------------------------
!        ... local variables
!-----------------------------------------------------------------
    
    integer  ::  kk                             ! vertical level index
    real(r8) ::  tp(ncol)                       ! 300/t [dimensionless]
    real(r8) ::  tinv(ncol)                     ! 1/t [K^{-1}]
    real(r8) ::  ko(ncol)   
    real(r8) ::  kinf(ncol)   
    real(r8) ::  fc(ncol)   
    real(r8) ::  sqrt_t(ncol)                   ! sqrt( temp )
    real(r8) ::  exp_fac(ncol)                  ! vector exponential




    level_loop : do kk = 1,pver
       tinv(:)           = 1._r8 / temp(:ncol,kk)
!BJG       tp(:)             = 300._r8 * tinv(:)
       tp(:)             = t0 * tinv(:)
       sqrt_t(:)         = sqrt( temp(:ncol,kk) )

!-----------------------------------------------------------------
!	... ho2 + ho2 --> h2o2
!	note: this rate involves the water vapor number density
!-----------------------------------------------------------------
       if( usr_HO2_HO2_ndx > 0 ) then

          call comp_exp( exp_fac, 430._r8*tinv, ncol )
          ko(:)   = 3.5e-13_r8 * exp_fac(:)
          call comp_exp( exp_fac, 1000._r8*tinv, ncol )
          kinf(:) = 1.7e-33_r8 * mtot(:,kk) * exp_fac(:)
          call comp_exp( exp_fac, 2200._r8*tinv, ncol )
          fc(:) = 1._r8 + 1.4e-21_r8 * invariants(:,kk,inv_h2o_ndx) * exp_fac(:)

          rxt(:,kk,usr_HO2_HO2_ndx) = (ko(:) + kinf(:)) * fc(:)

       endif

!-----------------------------------------------------------------
!       ... DMS + OH  --> .5 * SO2
!-----------------------------------------------------------------
       if( usr_DMS_OH_ndx > 0 ) then
          call comp_exp( exp_fac, 7460._r8*tinv, ncol )
          ko(:) = 1._r8 + 5.5e-31_r8 * exp_fac * mtot(:,kk) * 0.21_r8
          call comp_exp( exp_fac, 7810._r8*tinv, ncol )
          rxt(:,kk,usr_DMS_OH_ndx) = 1.7e-42_r8 * exp_fac * mtot(:,kk) * 0.21_r8 / ko(:)
       endif

!-----------------------------------------------------------------
!       ... SO2 + OH  --> SO4  (REFERENCE?? - not Liao)
!-----------------------------------------------------------------
       if( usr_SO2_OH_ndx > 0 ) then
          fc(:) = 3.0e-31_r8 *(300._r8*tinv(:))**3.3_r8
          ko(:) = fc(:)*mtot(:,kk)/(1._r8 + fc(:)*mtot(:,kk)/1.5e-12_r8) 
          rxt(:,kk,usr_SO2_OH_ndx) = ko(:)*.6_r8**(1._r8 + (log10(fc(:)*mtot(:,kk)/1.5e-12_r8))**2._r8)**(-1._r8)
       endif

    enddo level_loop

  end subroutine usrrxt

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
  subroutine comp_exp( x, y, n )

    implicit none

    real(r8), intent(out) :: x(:)
    real(r8), intent(in)  :: y(:)
    integer,  intent(in)  :: n
    
#ifdef IBM
    call vexp( x, y, n )
#else
    x(:n) = exp( y(:n) )
#endif

  end subroutine comp_exp

end module mo_usrrxt
