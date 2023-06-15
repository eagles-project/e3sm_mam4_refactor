
module mo_usrrxt

  use shr_kind_mod, only : r8 => shr_kind_r8
  use cam_logfile,  only : iulog
  use ppgrid,       only : pver, pcols
  use cam_abortutils,   only : endrun

  implicit none

  private
  public :: usrrxt, usrrxt_inti, usrrxt_hrates

  save

  integer :: usr_O_O2_ndx
  integer :: usr_HO2_HO2_ndx
  integer :: usr_N2O5_M_ndx
  integer :: usr_HNO3_OH_ndx
  integer :: usr_HO2NO2_M_ndx
  integer :: usr_N2O5_aer_ndx
  integer :: usr_NO3_aer_ndx
  integer :: usr_NO2_aer_ndx
  integer :: usr_CO_OH_a_ndx
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
  
  integer :: usr_OA_O2_NDX
  integer :: usr_XNO2NO3_M_ndx
  integer :: usr_NO2XNO3_M_ndx
  integer :: usr_XHNO3_OH_ndx
  integer :: usr_XHO2NO2_M_ndx
  integer :: usr_XNO2NO3_aer_ndx
  integer :: usr_NO2XNO3_aer_ndx
  integer :: usr_XNO3_aer_ndx
  integer :: usr_XNO2_aer_ndx
  integer :: usr_XPAN_M_ndx
  integer :: usr_XMPAN_M_ndx
  integer :: usr_MCO3_XNO2_ndx
  
  integer :: usr_C2O3_NO2_ndx
  integer :: usr_C2H4_OH_ndx
  integer :: usr_XO2N_HO2_ndx
  integer :: usr_C2O3_XNO2_ndx
  
  integer :: tag_XO2N_NO_ndx
  integer :: tag_XO2_HO2_ndx
  integer :: tag_XO2_NO_ndx
  
  integer :: usr_O_O_ndx
  integer :: usr_CL2O2_M_ndx
  integer :: usr_SO3_H2O_ndx
  integer :: tag_CLO_CLO_ndx
  
  integer :: ion1_ndx, ion2_ndx, ion3_ndx, ion11_ndx
  integer :: elec1_ndx, elec2_ndx, elec3_ndx
  integer :: het1_ndx

  integer :: usr_oh_co_ndx, het_no2_h2o_ndx, usr_oh_dms_ndx, aq_so2_h2o2_ndx, aq_so2_o3_ndx

  integer :: h2o_ndx, so4_ndx, cb2_ndx, oc2_ndx, soa_ndx, nit_ndx

!lke++
  integer :: usr_COhc_OH_ndx
  integer :: usr_COme_OH_ndx
  integer :: usr_CO01_OH_ndx
  integer :: usr_CO02_OH_ndx
  integer :: usr_CO03_OH_ndx
  integer :: usr_CO04_OH_ndx
  integer :: usr_CO05_OH_ndx
  integer :: usr_CO06_OH_ndx
  integer :: usr_CO07_OH_ndx
  integer :: usr_CO08_OH_ndx
  integer :: usr_CO09_OH_ndx
  integer :: usr_CO10_OH_ndx
  integer :: usr_CO11_OH_ndx
  integer :: usr_CO12_OH_ndx
  integer :: usr_CO13_OH_ndx
  integer :: usr_CO14_OH_ndx
  integer :: usr_CO15_OH_ndx
  integer :: usr_CO16_OH_ndx
  integer :: usr_CO17_OH_ndx
  integer :: usr_CO18_OH_ndx
  integer :: usr_CO19_OH_ndx
  integer :: usr_CO20_OH_ndx
  integer :: usr_CO21_OH_ndx
  integer :: usr_CO22_OH_ndx
  integer :: usr_CO23_OH_ndx
  integer :: usr_CO24_OH_ndx
  integer :: usr_CO25_OH_ndx
  integer :: usr_CO26_OH_ndx
  integer :: usr_CO27_OH_ndx
  integer :: usr_CO28_OH_ndx
  integer :: usr_CO29_OH_ndx
  integer :: usr_CO30_OH_ndx
  integer :: usr_CO31_OH_ndx
  integer :: usr_CO32_OH_ndx
  integer :: usr_CO33_OH_ndx
  integer :: usr_CO34_OH_ndx
  integer :: usr_CO35_OH_ndx
  integer :: usr_CO36_OH_ndx
  integer :: usr_CO37_OH_ndx
  integer :: usr_CO38_OH_ndx
  integer :: usr_CO39_OH_ndx
  integer :: usr_CO40_OH_ndx
  integer :: usr_CO41_OH_ndx
  integer :: usr_CO42_OH_ndx
!lke--

  logical :: has_aerosols

  real(r8), parameter :: t0     = 300._r8                ! K
  real(r8), parameter :: trlim2 = 17._r8/3._r8           ! K
  real(r8), parameter :: trlim3 = 15._r8/3._r8           ! K

  logical :: has_ion_rxts

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
    usr_CO_OH_a_ndx      = get_rxt_ndx( 'usr_CO_OH_a' )
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
 !
 ! additional reactions for O3A/XNO
 !
    usr_OA_O2_ndx        = get_rxt_ndx( 'usr_OA_O2' )
    usr_XNO2NO3_M_ndx    = get_rxt_ndx( 'usr_XNO2NO3_M' )
    usr_NO2XNO3_M_ndx    = get_rxt_ndx( 'usr_NO2XNO3_M' )
    usr_XNO2NO3_aer_ndx  = get_rxt_ndx( 'usr_XNO2NO3_aer' )
    usr_NO2XNO3_aer_ndx  = get_rxt_ndx( 'usr_NO2XNO3_aer' )
    usr_XHNO3_OH_ndx     = get_rxt_ndx( 'usr_XHNO3_OH' )
    usr_XNO3_aer_ndx     = get_rxt_ndx( 'usr_XNO3_aer' )
    usr_XNO2_aer_ndx     = get_rxt_ndx( 'usr_XNO2_aer' )
    usr_MCO3_XNO2_ndx    = get_rxt_ndx( 'usr_MCO3_XNO2' )
    usr_XPAN_M_ndx       = get_rxt_ndx( 'usr_XPAN_M' )
    usr_XMPAN_M_ndx      = get_rxt_ndx( 'usr_XMPAN_M' )
    usr_XHO2NO2_M_ndx    = get_rxt_ndx( 'usr_XHO2NO2_M' )
!
! reduced hydrocarbon chemistry
!
    usr_C2O3_NO2_ndx     = get_rxt_ndx( 'usr_C2O3_NO2' )
    usr_C2H4_OH_ndx      = get_rxt_ndx( 'usr_C2H4_OH' )
    usr_XO2N_HO2_ndx     = get_rxt_ndx( 'usr_XO2N_HO2' )
    usr_C2O3_XNO2_ndx    = get_rxt_ndx( 'usr_C2O3_XNO2' )
!
    tag_XO2N_NO_ndx      = get_rxt_ndx( 'tag_XO2N_NO' )
    tag_XO2_HO2_ndx      = get_rxt_ndx( 'tag_XO2_HO2' )
    tag_XO2_NO_ndx       = get_rxt_ndx( 'tag_XO2_NO' )
!
! stratospheric chemistry
!
    usr_O_O_ndx          = get_rxt_ndx( 'usr_O_O' )
    usr_CL2O2_M_ndx      = get_rxt_ndx( 'usr_CL2O2_M' )
    usr_SO3_H2O_ndx      = get_rxt_ndx( 'usr_SO3_H2O' )
!    
    tag_CLO_CLO_ndx      = get_rxt_ndx( 'tag_CLO_CLO' )
!
! stratospheric aerosol chemistry
!
    het1_ndx             = get_rxt_ndx( 'het1' )
!
! ion chemistry
!   
    ion1_ndx  = get_rxt_ndx( 'ion_Op_O2' )
    ion2_ndx  = get_rxt_ndx( 'ion_Op_N2' )
    ion3_ndx  = get_rxt_ndx( 'ion_N2p_Oa' )
    ion11_ndx = get_rxt_ndx( 'ion_N2p_Ob' )

    elec1_ndx  = get_rxt_ndx( 'elec1' )
    elec2_ndx  = get_rxt_ndx( 'elec2' )
    elec3_ndx  = get_rxt_ndx( 'elec3' )

    has_ion_rxts = ion1_ndx>0 .and. ion2_ndx>0 .and. ion3_ndx>0 .and. elec1_ndx>0 &
                 .and. elec2_ndx>0 .and. elec3_ndx>0

    so4_ndx    = get_spc_ndx( 'SO4' )
    cb2_ndx    = get_spc_ndx( 'CB2' )
    oc2_ndx    = get_spc_ndx( 'OC2' )
    soa_ndx    = get_spc_ndx( 'SOA' )
    nit_ndx    = get_spc_ndx( 'NH4NO3' )
    h2o_ndx    = get_spc_ndx( 'H2O' )

    !
    ! llnl super fast
    !
    usr_oh_co_ndx  = get_rxt_ndx( 'usr_oh_co' )
    het_no2_h2o_ndx  = get_rxt_ndx( 'het_no2_h2o' )
    usr_oh_dms_ndx  = get_rxt_ndx( 'usr_oh_dms' )
    aq_so2_h2o2_ndx  = get_rxt_ndx( 'aq_so2_h2o2' )
    aq_so2_o3_ndx  = get_rxt_ndx( 'aq_so2_o3' )
    
!lke++
! CO tags
!
    usr_COhc_OH_ndx      = get_rxt_ndx( 'usr_COhc_OH' )
    usr_COme_OH_ndx      = get_rxt_ndx( 'usr_COme_OH' )
    usr_CO01_OH_ndx      = get_rxt_ndx( 'usr_CO01_OH' )
    usr_CO02_OH_ndx      = get_rxt_ndx( 'usr_CO02_OH' )
    usr_CO03_OH_ndx      = get_rxt_ndx( 'usr_CO03_OH' )
    usr_CO04_OH_ndx      = get_rxt_ndx( 'usr_CO04_OH' )
    usr_CO05_OH_ndx      = get_rxt_ndx( 'usr_CO05_OH' )
    usr_CO06_OH_ndx      = get_rxt_ndx( 'usr_CO06_OH' )
    usr_CO07_OH_ndx      = get_rxt_ndx( 'usr_CO07_OH' )
    usr_CO08_OH_ndx      = get_rxt_ndx( 'usr_CO08_OH' )
    usr_CO09_OH_ndx      = get_rxt_ndx( 'usr_CO09_OH' )
    usr_CO10_OH_ndx      = get_rxt_ndx( 'usr_CO10_OH' )
    usr_CO11_OH_ndx      = get_rxt_ndx( 'usr_CO11_OH' )
    usr_CO12_OH_ndx      = get_rxt_ndx( 'usr_CO12_OH' )
    usr_CO13_OH_ndx      = get_rxt_ndx( 'usr_CO13_OH' )
    usr_CO14_OH_ndx      = get_rxt_ndx( 'usr_CO14_OH' )
    usr_CO15_OH_ndx      = get_rxt_ndx( 'usr_CO15_OH' )
    usr_CO16_OH_ndx      = get_rxt_ndx( 'usr_CO16_OH' )
    usr_CO17_OH_ndx      = get_rxt_ndx( 'usr_CO17_OH' )
    usr_CO18_OH_ndx      = get_rxt_ndx( 'usr_CO18_OH' )
    usr_CO19_OH_ndx      = get_rxt_ndx( 'usr_CO19_OH' )
    usr_CO20_OH_ndx      = get_rxt_ndx( 'usr_CO20_OH' )
    usr_CO21_OH_ndx      = get_rxt_ndx( 'usr_CO21_OH' )
    usr_CO22_OH_ndx      = get_rxt_ndx( 'usr_CO22_OH' )
    usr_CO23_OH_ndx      = get_rxt_ndx( 'usr_CO23_OH' )
    usr_CO24_OH_ndx      = get_rxt_ndx( 'usr_CO24_OH' )
    usr_CO25_OH_ndx      = get_rxt_ndx( 'usr_CO25_OH' )
    usr_CO26_OH_ndx      = get_rxt_ndx( 'usr_CO26_OH' )
    usr_CO27_OH_ndx      = get_rxt_ndx( 'usr_CO27_OH' )
    usr_CO28_OH_ndx      = get_rxt_ndx( 'usr_CO28_OH' )
    usr_CO29_OH_ndx      = get_rxt_ndx( 'usr_CO29_OH' )
    usr_CO30_OH_ndx      = get_rxt_ndx( 'usr_CO30_OH' )
    usr_CO31_OH_ndx      = get_rxt_ndx( 'usr_CO31_OH' )
    usr_CO32_OH_ndx      = get_rxt_ndx( 'usr_CO32_OH' )
    usr_CO33_OH_ndx      = get_rxt_ndx( 'usr_CO33_OH' )
    usr_CO34_OH_ndx      = get_rxt_ndx( 'usr_CO34_OH' )
    usr_CO35_OH_ndx      = get_rxt_ndx( 'usr_CO35_OH' )
    usr_CO36_OH_ndx      = get_rxt_ndx( 'usr_CO36_OH' )
    usr_CO37_OH_ndx      = get_rxt_ndx( 'usr_CO37_OH' )
    usr_CO38_OH_ndx      = get_rxt_ndx( 'usr_CO38_OH' )
    usr_CO39_OH_ndx      = get_rxt_ndx( 'usr_CO39_OH' )
    usr_CO40_OH_ndx      = get_rxt_ndx( 'usr_CO40_OH' )
    usr_CO41_OH_ndx      = get_rxt_ndx( 'usr_CO41_OH' )
    usr_CO42_OH_ndx      = get_rxt_ndx( 'usr_CO42_OH' )
!lke--

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'usrrxt_inti: diagnostics '
       write(iulog,'(10i5)') usr_O_O2_ndx,usr_HO2_HO2_ndx,tag_NO2_NO3_ndx,usr_N2O5_M_ndx,tag_NO2_OH_ndx,usr_HNO3_OH_ndx &
                            ,tag_NO2_HO2_ndx,usr_HO2NO2_M_ndx,usr_N2O5_aer_ndx,usr_NO3_aer_ndx,usr_NO2_aer_ndx &
                            ,usr_CO_OH_b_ndx,tag_C2H4_OH_ndx,tag_C3H6_OH_ndx,tag_CH3CO3_NO2_ndx,usr_PAN_M_ndx,usr_CH3COCH3_OH_ndx &
                            ,usr_MCO3_NO2_ndx,usr_MPAN_M_ndx,usr_XOOH_OH_ndx,usr_SO2_OH_ndx,usr_DMS_OH_ndx,usr_HO2_aer_ndx
    end if

  end subroutine usrrxt_inti

! BJG  subroutine usrrxt( rxt, temp, tempi, tempe, invariants, h2ovmr,  ps, &
  subroutine usrrxt( rxt, temp, invariants, &
! BJG                     pmid, m, sulfate, mmr, relhum, strato_sad, &
! BJG                     m, sulfate, relhum, &
                     mtot,  &
! BJG                     ltrop, ncol, sad_total, cwat, mbar, pbuf )
                     ncol )

!-----------------------------------------------------------------
!        ... set the user specified reaction rates
!-----------------------------------------------------------------
    
    use mo_constants,  only : pi, avo => avogadro, boltz_cgs, rgas
    use chem_mods,     only : nfs, rxntot, gas_pcnst, inv_m_ndx=>indexm
!BJG    use mo_chem_utls,  only : get_rxt_ndx, get_spc_ndx
    use mo_setinv,     only : inv_o2_ndx=>o2_ndx, inv_h2o_ndx=>h2o_ndx
!BJG    use physics_buffer,only : physics_buffer_desc
!BJG    use rad_constituents,only : rad_cnst_get_info
    use modal_aero_data,   only: ntot_amode
    implicit none

!-----------------------------------------------------------------
!        ... dummy arguments
!-----------------------------------------------------------------
    integer, intent(in)     :: ncol
!BJG    integer, intent(in)     :: ltrop(pcols)               ! tropopause vertical index
    real(r8), intent(in)    :: temp(pcols,pver)           ! temperature (K); neutral temperature
!BJG    real(r8), intent(in)    :: tempi(pcols,pver)          ! ionic temperature (K); only used if ion chemistry
!BJG    real(r8), intent(in)    :: tempe(pcols,pver)          ! electronic temperature (K); only used if ion chemistry
    real(r8), intent(in)    :: mtot(ncol,pver)               ! total atm density (/cm^3)
!BJG    real(r8), intent(in)    :: sulfate(ncol,pver)         ! sulfate aerosol (mol/mol)
!BJG    real(r8), intent(in)    :: strato_sad(pcols,pver)     ! stratospheric aerosol sad (1/cm)
!BJG    real(r8), intent(in)    :: h2ovmr(ncol,pver)          ! water vapor (mol/mol)
!BJG    real(r8), intent(in)    :: relhum(ncol,pver)          ! relative humidity
!BJG    real(r8), intent(in)    :: pmid(pcols,pver)           ! midpoint pressure (Pa)
!BJG    real(r8), intent(in)    :: ps(pcols)                  ! surface pressure (Pa)
    real(r8), intent(in)    :: invariants(ncol,pver,nfs)  ! invariants density (/cm^3)
!BJG    real(r8), intent(in)    :: mmr(pcols,pver,gas_pcnst)  ! species concentrations (kg/kg)
!BJG    real(r8), intent(in)    :: cwat(ncol,pver) !PJC Condensed Water (liquid+ice) (kg/kg)
!BJG    real(r8), intent(in)    :: mbar(ncol,pver) !PJC Molar mass of air (g/mol)
    real(r8), intent(inout) :: rxt(ncol,pver,rxntot)      ! gas phase rates
!BJG    real(r8), intent(out)   :: sad_total(pcols,pver)      ! total surface area density (cm2/cm3)
!BJG    type(physics_buffer_desc), pointer :: pbuf(:)
      
!-----------------------------------------------------------------
!        ... local variables
!-----------------------------------------------------------------
    
!BJG    real(r8), parameter :: dg = 0.1_r8            ! mole diffusion =0.1 cm2/s (Dentener, 1993)

!-----------------------------------------------------------------
! 	... reaction probabilities for heterogeneous reactions
!-----------------------------------------------------------------
!BJG    real(r8), parameter :: gamma_n2o5 = 0.10_r8         ! from Jacob, Atm Env, 34, 2131, 2000
!BJG    real(r8), parameter :: gamma_ho2  = 0.20_r8         ! 
!BJG    real(r8), parameter :: gamma_no2  = 0.0001_r8       ! 
!BJG    real(r8), parameter :: gamma_no3  = 0.001_r8        ! 

!BJG    integer  ::  icol, kk
    integer  ::  kk
    real(r8) ::  tp(ncol)                       ! 300/t
    real(r8) ::  tinv(ncol)                     ! 1/t
    real(r8) ::  ko(ncol)   
!BJG    real(r8) ::  term1(ncol)
!BJG    real(r8) ::  term2(ncol)
    real(r8) ::  kinf(ncol)   
    real(r8) ::  fc(ncol)   
!BJG    real(r8) ::  xr(ncol)   
!BJG    real(r8) ::  sur(ncol)   
    real(r8) ::  sqrt_t(ncol)                   ! sqrt( temp )
    real(r8) ::  exp_fac(ncol)                  ! vector exponential
!BJG    real(r8) ::  lwc(ncol)   
!BJG    real(r8) ::  ko_m(ncol)   
!BJG    real(r8) ::  k0(ncol)   
!BJG    real(r8) ::  kinf_m(ncol)   
!BJG    real(r8) ::  o2(ncol)
!BJG    real(r8) ::  c_n2o5, c_ho2, c_no2, c_no3
!BJG    real(r8) ::  amas
    !-----------------------------------------------------------------
    !	... density of sulfate aerosol
    !-----------------------------------------------------------------
!BJG    real(r8), parameter :: gam1 = 0.04_r8                 ! N2O5+SUL ->2HNO3
!BJG    real(r8), parameter :: wso4 = 98._r8
!BJG    real(r8), parameter :: den  = 1.15_r8                 ! each molecule of SO4(aer) density g/cm3
    !-------------------------------------------------
    ! 	... volume of sulfate particles
    !           assuming mean rm 
    !           continient 0.05um  0.07um  0.09um
    !           ocean      0.09um  0.25um  0.37um
    !                      0.16um                  Blake JGR,7195, 1995
    !-------------------------------------------------
!BJG    real(r8), parameter :: rm1  = 0.16_r8*1.e-4_r8             ! mean radii in cm
!BJG    real(r8), parameter :: fare = 4._r8*pi*rm1*rm1             ! each mean particle(r=0.1u) area   cm2/cm3

    !-----------------------------------------------------------------------
    !        ... Aqueous phase sulfur quantities for SO2 + H2O2 and SO2 + O3
    !-----------------------------------------------------------------------
!BJG    real(r8), parameter  :: HENRY298_H2O2 =  7.45e+04_r8
!BJG    real(r8), parameter  :: H298_H2O2     = -1.45e+04_r8
!BJG    real(r8), parameter  :: HENRY298_SO2  =  1.23e+00_r8
!BJG    real(r8), parameter  :: H298_SO2      = -6.25e+03_r8
!BJG    real(r8), parameter  :: K298_SO2_HSO3 =  1.3e-02_r8
!BJG    real(r8), parameter  :: H298_SO2_HSO3 = -4.16e+03_r8
!BJG    real(r8), parameter  :: R_CONC        =  82.05e+00_r8 / avo
!BJG    real(r8), parameter  :: R_CAL         =  rgas * 0.239006e+00_r8
!BJG    real(r8), parameter  :: K_AQ          =  7.57e+07_r8
!BJG    real(r8), parameter  :: ER_AQ         =  4.43e+03_r8

!BJG    real(r8), parameter  :: HENRY298_O3   =  1.13e-02_r8
!BJG    real(r8), parameter  :: H298_O3       = -5.04e+03_r8
!BJG    real(r8), parameter  :: K298_HSO3_SO3 =  6.6e-08_r8
!BJG    real(r8), parameter  :: H298_HSO3_SO3 = -2.23e+03_r8
!BJG    real(r8), parameter  :: K0_AQ         =  2.4e+04_r8
!BJG    real(r8), parameter  :: ER0_AQ        =  0.0e+00_r8
!BJG    real(r8), parameter  :: K1_AQ         =  3.7e+05_r8
!BJG    real(r8), parameter  :: ER1_AQ        =  5.53e+03_r8
!BJG    real(r8), parameter  :: K2_AQ         =  1.5e+09_r8
!BJG    real(r8), parameter  :: ER2_AQ        =  5.28e+03_r8

!BJG    real(r8), parameter  :: pH            =  4.5e+00_r8
    
!BJG    real(r8), pointer :: sfc(:), dm_aer(:)
!BJG    integer :: ntot_amode

!BJG    real(r8), pointer :: sfc_array(:,:,:), dm_array(:,:,:)
    
    ! get info about the modal aerosols
    ! get ntot_amode
!BJG    call rad_cnst_get_info(0, nmodes=ntot_amode)

!BJG    if (ntot_amode>0) then
!BJG       allocate(sfc_array(pcols,pver,ntot_amode), dm_array(pcols,pver,ntot_amode) )
!BJG    else
!BJG       allocate(sfc_array(pcols,pver,5), dm_array(pcols,pver,5) )
!BJG    endif

!BJG    sfc_array(:,:,:) = 0._r8
!BJG    dm_array(:,:,:) = 0._r8
! BJG    sad_total(:,:) = 0._r8


    level_loop : do kk = 1,pver
       tinv(:)           = 1._r8 / temp(:ncol,kk)
       tp(:)             = 300._r8 * tinv(:)
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

!BJG          if( h2o_ndx > 0 ) then
!BJG             fc(:) = 1._r8 + 1.4e-21_r8 * m(:,kk) * h2ovmr(:,kk) * exp_fac(:)
!BJG          else
             fc(:) = 1._r8 + 1.4e-21_r8 * invariants(:,kk,inv_h2o_ndx) * exp_fac(:)
!BJG          endif
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

!BJG      deallocate( sfc_array, dm_array )

  end subroutine usrrxt

      subroutine usrrxt_hrates( rxt, tempn, tempi, tempe, invariants, &
				h2ovmr, pmid, m, ncol, kbot )
!-----------------------------------------------------------------
!        ... set the user specified reaction rates for heating
!-----------------------------------------------------------------

      use shr_kind_mod,  only : r8 => shr_kind_r8
      use chem_mods,     only : nfs, rxntot
      use ppgrid,        only : pver, pcols

      implicit none

!-----------------------------------------------------------------
!        ... dummy arguments
!-----------------------------------------------------------------
      integer, intent(in)     :: ncol                         ! number columns in chunk
      integer, intent(in)     :: kbot                         ! heating levels
      real(r8), intent(in)    :: tempn(pcols,pver)            ! neutral temperature (K)
      real(r8), intent(in)    :: tempi(pcols,pver)            ! ion temperature (K)
      real(r8), intent(in)    :: tempe(pcols,pver)            ! electron temperature (K)
      real(r8), intent(in)    :: m(ncol,pver)                 ! total atm density (1/cm^3)
      real(r8), intent(in)    :: h2ovmr(ncol,pver)            ! water vapor (vmr)
      real(r8), intent(in)    :: pmid(pcols,pver)             ! midpoint pressure (Pa)
      real(r8), intent(in)    :: invariants(ncol,pver,nfs)    ! invariants density (1/cm^3)
      real(r8), intent(inout) :: rxt(ncol,pver,rxntot)        ! gas phase rates
      
!-----------------------------------------------------------------
!        ... local variables
!-----------------------------------------------------------------

      integer  ::  k
      real(r8), dimension(ncol) :: &
                   tp, &
                   tinv, &
                   ko, &
                   kinf, &
                   fc, &
                   xr                       ! factor to increase particle radii depending on rel hum

!-----------------------------------------------------------------
!	... o + o2 + m --> o3 + m
!-----------------------------------------------------------------
      do k = 1,kbot
         tinv(:ncol)       = 1._r8 / tempn(:ncol,k)
         tp(:)             = 300._r8 * tinv(:)
         rxt(:,k,usr_O_O2_ndx) = 6.e-34_r8 * tp(:)**2.4_r8

!-----------------------------------------------------------------
!	... o + o + m -> o2 + m
!-----------------------------------------------------------------
         rxt(:,k,usr_O_O_ndx) = 2.76e-34_r8 * exp( 720.0_r8*tinv(:) )

!-----------------------------------------------------------------
!	... ho2 + ho2 --> h2o2
!	Note: this rate involves the water vapor number density
!-----------------------------------------------------------------
         ko(:)   = 3.5e-13_r8 * exp( 430._r8*tinv(:) )
         kinf(:) = 1.7e-33_r8 * m(:,k) * exp( 1000._r8*tinv(:) )
         fc(:)   = 1._r8 + 1.4e-21_r8 * m(:,k) * h2ovmr(:,k) * exp( 2200._r8*tinv(:) )
         rxt(:,k,usr_HO2_HO2_ndx) = (ko(:) + kinf(:)) * fc(:)

      end do

!-----------------------------------------------------------------
! 	... the ionic rates
!-----------------------------------------------------------------
      if ( has_ion_rxts ) then
         level_loop2 :  do k = 1,kbot
            tp(:ncol)         = (2._r8*tempi(:ncol,k) + tempn(:ncol,k)) / ( 3._r8 * t0 )
            tp(:)             = max( min( tp(:),20._r8 ),1._r8 )
            rxt(:,k,ion1_ndx) = 2.82e-11_r8 + tp(:)*(-7.74e-12_r8 + tp(:)*(1.073e-12_r8  &
                 + tp(:)*(-5.17e-14_r8 + 9.65e-16_r8*tp(:))))
            tp(:ncol)         = (.6363_r8*tempi(:ncol,k) + .3637_r8*tempn(:ncol,k)) / t0
            tp(:)             = max( min( tp(:),trlim2 ),1._r8 )
            rxt(:,k,ion2_ndx) = 1.533e-12_r8 + tp(:)*(-5.92e-13_r8 + tp(:)*8.6e-14_r8)
            tp(:ncol)         = 2._r8 * t0 /(tempi(:ncol,k) + tempn(:ncol,k))
            where( tp(:ncol) < trlim3 )
               rxt(:,k,ion3_ndx)  = 1.4e-10_r8 * tp(:)**.44_r8
            elsewhere
               rxt(:,k,ion3_ndx)  = 5.2e-11_r8 / tp(:)**.2_r8
            endwhere
            tp(:ncol)          = t0 / tempe(:ncol,k)
            rxt(:,k,elec1_ndx) = 4.e-7_r8 * tp(:)**.85_r8
            rxt(:,k,elec3_ndx) = 1.8e-7_r8 * tp(:)**.39_r8
            where( tp(:ncol) < 4._r8 )
               rxt(:,k,elec2_ndx) = 2.7e-7_r8 * tp(:)**.7_r8
            elsewhere
               rxt(:,k,elec2_ndx) = 1.6e-7_r8 * tp(:)**.55_r8
            endwhere
         end do level_loop2
      endif
      end subroutine usrrxt_hrates

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

  !-------------------------------------------------------------------------
  !  Heterogeneous reaction rates for uptake of a gas on an aerosol:
  !-------------------------------------------------------------------------
  function hetrxtrate( sfc, dm_aer, dg_gas, c_gas, gamma_gas ) result(rate)

    real(r8), intent(in) :: sfc(:)
    real(r8), intent(in) :: dm_aer(:)
    real(r8), intent(in) :: dg_gas
    real(r8), intent(in) :: c_gas
    real(r8), intent(in) :: gamma_gas
    real(r8) :: rate

    real(r8),allocatable :: rxt(:)
    integer :: n, i

    n = size(sfc)

    allocate(rxt(n))
    do i=1,n
       rxt(i) = sfc(i) / (0.5_r8*dm_aer(i)/dg_gas + (4._r8/(c_gas*gamma_gas)))
    enddo

    rate = sum(rxt)

    deallocate(rxt)

  endfunction hetrxtrate

end module mo_usrrxt
