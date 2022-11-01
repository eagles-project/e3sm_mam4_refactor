!----------------------------------------------------------------------
!BOP
!
! !MODULE: modal_aero_amicphys_control 
!
   module modal_aero_amicphys_control

! !USES:
  use shr_kind_mod,    only:  r8 => shr_kind_r8
  use chem_mods,       only:  gas_pcnst
  use physconst,       only:  pi
  use modal_aero_data, only:  ntot_aspectype, ntot_amode, nsoa, npoa, nbc
! use ref_pres,        only:  top_lev => clim_modal_aero_top_lev  ! this is for gg02a
  use ref_pres,        only:  top_lev => trop_cloud_top_lev       ! this is for ee02c

! !DESCRIPTION: This module contains the constants, parameters, and control variables 
!  separated from modal_aero_amicphys.
!
! !REVISION HISTORY:
!
!   RCE 07.04.13:  Adapted from MIRAGE2 code
!
!EOP
!----------------------------------------------------------------------

  implicit none

  public

! !PUBLIC DATA MEMBERS:
  type :: misc_vars_aa_type
! using this derived type reduces the number of changes needed to add more mosaic diagnostics to history
     real(r8) :: ncluster_tend_nnuc_1grid
#if ( defined( MOSAIC_SPECIES ) )
     real(r8) :: cnvrg_fail_1grid
     real(r8) :: max_kelvin_iter_1grid
     real(r8), dimension(5,4) :: xnerr_astem_negative_1grid
#endif
  end type misc_vars_aa_type

  logical, public :: mosaic = .true. !BSINGH -  Added logical for mosaic model

  integer, parameter :: pcnstxx = gas_pcnst

  real(r8), public :: n_so4_monolayers_pcage = huge(1.0_r8)
! number of so4(+nh4) monolayers needed to "age" a carbon particle

  real(r8), public :: dr_so4_monolayers_pcage = huge(1.0_r8)
! thickness of the so4 monolayers (m)
! for so4(+nh4), use bi-sulfate mw and 1.77 g/cm3,
!    --> 1 mol so4(+nh4)  = 65 cm^3 --> 1 molecule = (4.76e-10 m)^3
! aging criterion is approximate so do not try to distinguish
!    sulfuric acid, bisulfate, ammonium sulfate

#if ( defined( CAMBOX_ACTIVATE_THIS ) )
  integer, public :: cldy_rh_sameas_clear = 0
! this is only used for some specific box model tests
#endif

  integer, public :: mdo_gaexch_cldy_subarea = 0
! controls if gas condensation is done in cloudy subarea
!    1 = yes ; 0 = no

  integer, public :: gaexch_h2so4_uptake_optaa = 2
! controls treatment of h2so4 condensation in mam_gasaerexch_1subarea
!    1 = sequential   calc. of gas-chem prod then condensation loss
!    2 = simultaneous calc. of gas-chem prod and  condensation loss

  integer, public :: newnuc_h2so4_conc_optaa = 2
! controls treatment of h2so4 concentrationin mam_newnuc_1subarea
!    1 = use average value calculated in standard cam5.2.10 and earlier
!    2 = use average value calculated in mam_gasaerexch_1subarea
!   11 = use average of initial and final values from mam_gasaerexch_1subarea
!   12 = use final value from mam_gasaerexch_1subarea

  integer, public :: rename_method_optaa = 40
! controls renaming parameterization

  integer, public :: update_qaerwat = 0
  integer, public :: update_dgncur_a = 0
  integer, public :: update_dgncur_awet = 0
! controls updating of qaerwat
! controls updating of dgncur_a
! controls updating of dgncur_awet and wetdens_host

  real (r8) :: newnuc_adjust_factor_dnaitdt = 1.0_r8
  real (r8) :: newnuc_adjust_factor_pbl     = 1.0_r8


#if ( defined MODAL_AERO_3MODE )
  integer, parameter :: max_gas = nsoa + 1
  ! the +3 in max_aer are dst, ncl, so4
  integer, parameter :: max_aer = nsoa + npoa + nbc + 3
#elif ( defined MODAL_AERO_4MODE )
  integer, parameter :: max_gas = nsoa + 1
  ! the +3 in max_aer are dst, ncl, so4
  integer, parameter :: max_aer = nsoa + npoa + nbc + 3
#elif ( defined MODAL_AERO_4MODE_MOM )
  integer, parameter :: max_gas = nsoa + 1
  ! the +4 in max_aer are dst, ncl, so4, mom
  integer, parameter :: max_aer = nsoa + npoa + nbc + 4
#elif ( ( defined MODAL_AERO_7MODE ) && ( defined MOSAIC_SPECIES ) )
  integer, parameter :: max_gas = nsoa + 4
  ! the +8 in max_aer are dst, ncl(=na), so4, no3, cl, nh4, ca, co3 
  integer, parameter :: max_aer = nsoa + npoa + nbc + 8
#elif ( defined MODAL_AERO_7MODE )
  integer, parameter :: max_gas = nsoa + 2
  ! the +4 in max_aer are dst, ncl, so4, nh4
  integer, parameter :: max_aer = nsoa + npoa + nbc + 4
#elif ( defined MODAL_AERO_8MODE )
  integer, parameter :: max_gas = nsoa + 2
  ! the +4 in max_aer are dst, ncl, so4, mom ???
  integer, parameter :: max_aer = nsoa + npoa + nbc + 4
#elif ( defined MODAL_AERO_9MODE )
  integer, parameter :: max_gas = nsoa + 2
  ! the +4+5 in max_aer are dst, ncl, so4, nh4 and 5 marine organics
  integer, parameter :: max_aer = nsoa + npoa + nbc + 4 + 5
#endif

#if (( defined MODAL_AERO_8MODE ) || ( defined MODAL_AERO_4MODE ) || ( defined MODAL_AERO_4MODE_MOM ))
  integer, parameter :: ntot_amode_extd = ntot_amode
#else
  integer, parameter :: ntot_amode_extd = ntot_amode + 1
! integer, parameter :: ntot_amode_extd = ntot_amode
#endif

  integer, parameter :: max_mode_fresh = 1

  integer, parameter :: max_mode = ntot_amode_extd + max_mode_fresh
  public max_mode !BSINGH - used in module_mosaic_cam_init.F90


  integer, parameter :: max_agepair = 1
  integer, parameter :: max_coagpair = 3 

  integer, parameter :: maxsubarea = 2

  integer, parameter :: nqtendaa = 4
  integer, parameter :: iqtend_cond = 1
  integer, parameter :: iqtend_rnam = 2
  integer, parameter :: iqtend_nnuc = 3
  integer, parameter :: iqtend_coag = 4
  integer, parameter :: nqqcwtendaa = 1
  integer, parameter :: iqqcwtend_rnam = 1

  integer, parameter :: iqqcwtend_match_iqtend(nqtendaa) = (/ 0, iqqcwtend_rnam, 0, 0 /)

  logical, parameter :: aging_include_seasalt = .false.
                      ! when .true., aging (by coagulation) includes contribution of seasalt
                      ! early versions of mam neglected the seasalt contribution

  ! species indices for various qgas_--- arrays
  integer :: igas_soa, igas_h2so4, igas_nh3, igas_hno3, igas_hcl
  ! species indices for various qaer_--- arrays
  !    when nsoa > 1, igas_soa and iaer_soa are indices of the first soa species
  !    when nbc  > 1, iaer_bc  is index of the first bc  species
  !    when npom > 1, iaer_pom is index of the first pom species
  integer :: iaer_bc, iaer_dst, iaer_ncl, iaer_nh4, iaer_pom, iaer_soa, iaer_so4, &
             iaer_mpoly, iaer_mprot, iaer_mlip, iaer_mhum, iaer_mproc, iaer_mom, &
             iaer_no3, iaer_cl, iaer_ca, iaer_co3
  integer :: i_agepair_pca, i_agepair_macc, i_agepair_mait
  integer :: lmap_gas(max_gas)
  integer :: lmap_aer(max_aer,max_mode), lmapbb_aer(max_aer,max_mode), &
             lmap_aercw(max_aer,max_mode)
  integer :: lmap_num(max_mode), lmap_numcw(max_mode)
  integer :: lmapcc_all(gas_pcnst)
  integer, parameter :: lmapcc_val_gas = 1, lmapcc_val_aer = 2, lmapcc_val_num = 3
  integer :: ngas, naer
  integer :: nacc, nait, npca, nufi, nmacc, nmait

  integer :: n_agepair, n_coagpair
  integer :: modefrm_agepair(max_agepair), modetoo_agepair(max_agepair)
  integer :: mode_aging_optaa(max_mode)
  integer :: modefrm_coagpair(max_coagpair), modetoo_coagpair(max_coagpair), &
             modeend_coagpair(max_coagpair)

  integer :: lun82,   lun97,   lun98,   lun13n,   lun15n
  logical :: ldiag82, ldiag97, ldiag98, ldiag13n, ldiag15n
  logical :: ldiagd1

  real(r8) :: accom_coef_gas(max_gas)
  real(r8) :: alnsg_aer(max_mode)
  real(r8) :: dgnum_aer(max_mode), dgnumhi_aer(max_mode), dgnumlo_aer(max_mode)
  real(r8) :: dens_aer(max_aer)
  real(r8) :: dens_so4a_host
  real(r8) :: fac_m2v_aer(max_aer)        ! converts (mol-aero/mol-air) to (m3-aero/mol-air)
  real(r8) :: fac_eqvso4hyg_aer(max_aer)  ! converts a species volume to a volume of so4
                                          !    (or nh4hso4) having same hygroscopicity
  real(r8) :: fac_m2v_eqvhyg_aer(max_aer) ! = fac_m2v_aer * fac_eqvso4hyg_aer

  real(r8) :: fcvt_gas(max_gas), fcvt_aer(max_aer), fcvt_num, fcvt_wtr
  real(r8) :: fcvt_dgnum_dvolmean(max_mode)
  real(r8) :: hygro_aer(max_aer)
  real(r8) :: mw_gas(max_gas), mw_aer(max_aer)
  real(r8) :: mwhost_gas(max_gas), mwhost_aer(max_aer), mwhost_num
  real(r8) :: mw_nh4a_host, mw_so4a_host
  real(r8) :: mwuse_soa(nsoa), mwuse_poa(npoa)
  real(r8) :: sigmag_aer(max_mode)
  real(r8) :: vol_molar_gas(max_gas)

! following were used in aging calcs but are no longer needed
!    fac_m2v_so4, fac_m2v_nh4, fac_m2v_soa(:)
!    fac_m2v_pcarbon(:)
!    soa_equivso4_factor(:)


  character(len=16) :: name_gas(max_gas), name_aerpfx(max_aer), &
     name_aer(max_aer,max_mode), name_aercw(max_aer,max_mode), &
     name_num(max_mode), name_numcw(max_mode)

  character(len=8) :: suffix_q_coltendaa(nqtendaa) = &
     (/ '_sfgaex1', '_sfgaex2', '_sfnnuc1', '_sfcoag1' /)
  character(len=8) :: suffix_qqcw_coltendaa(nqqcwtendaa) = &
                    '_sfgaex2'

  logical :: do_q_coltendaa(gas_pcnst,nqtendaa) = .false.
  logical :: do_qqcw_coltendaa(gas_pcnst,nqqcwtendaa) = .false.

! *** following 3 variables should eventually be in modal_aero_data
  real(r8) :: specmw2_amode(ntot_aspectype,ntot_amode)
  real(r8) :: specdens2_amode(ntot_aspectype,ntot_amode)
  real(r8) :: spechygro2(ntot_aspectype,ntot_amode)

end module modal_aero_amicphys_control
