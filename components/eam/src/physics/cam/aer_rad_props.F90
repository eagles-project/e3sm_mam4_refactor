module aer_rad_props
#include "../../chemistry/yaml/common_files/common_uses.ymlf90"

  !------------------------------------------------------------------------------------------------
  ! Converts aerosol masses to bulk optical properties for sw and lw radiation
  ! computations.
  !------------------------------------------------------------------------------------------------

  use shr_kind_mod,     only: r8 => shr_kind_r8
  use ppgrid,           only: pcols, pver

  use radconstants,     only: nswbands
  use cam_history,      only: addfld, horiz_only, outfld, add_default
  use mam_support,      only: ptr2d_t

  implicit none
  private
  save

  public :: &
       aer_rad_props_init,        &
       aer_rad_props_sw,          & ! return SW optical props of aerosols
       aer_rad_props_lw             ! return LW optical props of aerosols

  ! Private data
  real(r8), parameter :: km_inv_to_m_inv = 0.001_r8      !1/km to 1/m
  !==============================================================================
contains
  !==============================================================================

  subroutine aer_rad_props_init()
    use phys_control, only: phys_getopts

    !Local variables
    logical                    :: history_amwg         ! output the variables used by the AMWG diag package
    logical                    :: history_aero_optics  ! Output aerosol optics diagnostics

    !----------------------------------------------------------------------------

    call phys_getopts( history_aero_optics_out    = history_aero_optics, &
         history_amwg_out           = history_amwg)

    call addfld ('AEROD_v', horiz_only, 'A', '1', &
         'Total Aerosol Optical Depth in visible band', flag_xyfill=.true.)

    !For testing puposes only, the following addfld call should be removed before merging to master
    call addfld ('extinct_lw_bnd7',(/ 'lev' /),    'A','1/m','EXTINCT LW H2O window band 7 output', flag_xyfill=.true.)
    call addfld ('extinct_lw_inp',(/ 'lev' /),    'A','1/km',&
         'EXTINCT LW H2O window band 7 output directly read from prescribed input file', flag_xyfill=.true.)
    call addfld ('extinct_sw_inp',(/ 'lev' /),    'A','1/km',&
         'Aerosol extinction directly read from prescribed input file', flag_xyfill=.true.)

    ! Contributions to AEROD_v from individual aerosols (climate species).

    ! Determine default fields
    if (history_amwg ) then
       call add_default ('AEROD_v', 1, ' ')
    endif

    if ( history_aero_optics ) then
       call add_default ('AEROD_v', 1, ' ')
    endif

  end subroutine aer_rad_props_init

  !==============================================================================

  subroutine aer_rad_props_sw(dt, lchnk, ncol, zi, pmid, pint, & !in
       temperature, zm, state_q, pdel, pdeldry, cldn, ssa_cmip6_sw, & !in
       af_cmip6_sw, ext_cmip6_sw,  nnite, idxnite, is_cmip6_volc, & !in
       qqcw, tau, tau_w, tau_w_g, tau_w_f) !out

    use modal_aer_opt,    only: modal_aero_sw
    use radconstants,     only: idx_sw_diag

    ! Return bulk layer tau, omega, g, f for all spectral intervals.

    ! Arguments
    integer,          intent(in) :: lchnk            ! number of chunks
    integer,          intent(in) :: ncol             ! number of columns
    real(r8),         intent(in) :: pmid(:,:)        ! midpoint pressure [Pa]
    real(r8),         intent(in) :: pint(:,:)        ! interface pressure [Pa]
    real(r8),         intent(in) :: temperature(:,:) ! temperature [K]
    real(r8),         intent(in) :: zm(:,:)          ! geopotential height above surface at midpoints [m]
    real(r8),         intent(in) :: zi(:,:)          ! geopotential height above surface at interfaces [m]
    real(r8), target, intent(in) :: state_q(:,:,:)
    real(r8),         intent(in) :: pdel(:,:)
    real(r8),         intent(in) :: pdeldry(:,:)
    real(r8),         intent(in) :: cldn(:,:)
    real(r8),         intent(in) :: ext_cmip6_sw(:,:,:)
    real(r8),         intent(in) :: ssa_cmip6_sw(:,:,:)
    real(r8),         intent(in) :: af_cmip6_sw(:,:,:)

    integer,          intent(in) :: nnite                ! number of night columns
    integer,          intent(in) :: idxnite(:)           ! local column indices of night columns
    logical,          intent(in) :: is_cmip6_volc        ! true if cmip6 style volcanic file is read otherwise false
    real(r8),         intent(in) :: dt                   ! time step (s)

    type(ptr2d_t),    intent(inout) :: qqcw(:)               ! Cloud borne aerosols mixing ratios [kg/kg or 1/kg]
    real(r8),         intent(out) :: tau    (pcols,0:pver,nswbands) ! aerosol extinction optical depth
    real(r8),         intent(out) :: tau_w  (pcols,0:pver,nswbands) ! aerosol single scattering albedo * tau
    real(r8),         intent(out) :: tau_w_g(pcols,0:pver,nswbands) ! aerosol assymetry parameter * tau * w
    real(r8),         intent(out) :: tau_w_f(pcols,0:pver,nswbands) ! aerosol forward scattered fraction * tau * w

    ! Local variables
    ! for cmip6 style volcanic file
    integer  :: trop_level(pcols)
    real(r8) :: ext_cmip6_sw_inv_m(pcols,pver,nswbands)! short wave extinction in the units of [1/m]
#include "../../chemistry/yaml/aer_rad_props/f90_yaml/aer_rad_props_sw_beg_yml.f90"

    !-----------------------------------------------------------------------------

    call outfld('extinct_sw_inp',ext_cmip6_sw(:,:,idx_sw_diag), pcols, lchnk)

    !FORTRAN REFACTOR: This is done to fill invalid values in columns where pcols>ncol
    !C++ port can ignore this as C++ model is a single column model
    ! initialize to conditions that would cause failure
    tau     (:,:,:) = -100._r8
    tau_w   (:,:,:) = -100._r8
    tau_w_g (:,:,:) = -100._r8
    tau_w_f (:,:,:) = -100._r8

    ! top layer (ilev = 0) has no aerosol (ie tau = 0)
    ! also initialize rest of layers to accumulate od's
    tau    (1:ncol,:,:) = 0._r8
    tau_w  (1:ncol,:,:) = 0._r8
    tau_w_g(1:ncol,:,:) = 0._r8
    tau_w_f(1:ncol,:,:) = 0._r8

    !Converting it from 1/km to 1/m
    ext_cmip6_sw_inv_m = ext_cmip6_sw * km_inv_to_m_inv

    !Find tropopause (or quit simulation if not found) as extinction should be applied only above tropopause
    trop_level(1:pcols) = tropopause_or_quit(lchnk, ncol, pmid, pint, temperature, zm, zi)

    !Special treatment for CMIP6 volcanic aerosols, where extinction, ssa
    !and af are directly read from the prescribed volcanic aerosol file
    call modal_aero_sw(dt, lchnk, ncol, state_q, zm, temperature, pmid, pdel, & !in
         pdeldry, cldn, nnite, idxnite, .true., & !in
         ext_cmip6_sw_inv_m(:,:,idx_sw_diag), trop_level, & !in
         qqcw, tau, tau_w, tau_w_g, tau_w_f) !out

    !Update tau, tau_w, tau_w_g, and tau_w_f with the read in values of extinction, ssa and asymmetry factors
    call volcanic_cmip_sw(ncol, zi, trop_level, ext_cmip6_sw_inv_m, ssa_cmip6_sw, af_cmip6_sw, & ! in
         tau, tau_w, tau_w_g, tau_w_f)  ! inout

    ! Diagnostic output of total aerosol optical properties
    ! currently implemented for climate list only
    call aer_vis_diag_out(lchnk, ncol, nnite, idxnite, tau(:,:,idx_sw_diag)) !in
#include "../../chemistry/yaml/aer_rad_props/f90_yaml/aer_rad_props_sw_end_yml.f90"

  end subroutine aer_rad_props_sw

  !==============================================================================
  subroutine aer_rad_props_lw(is_cmip6_volc, dt, lchnk, ncol, pmid, pint, & !in
       temperature, zm, zi, state_q, pdel, pdeldry, cldn, ext_cmip6_lw, &    !in
       qqcw, odap_aer) !out

    use modal_aer_opt,    only: modal_aero_lw
    use radconstants,     only: nlwbands, idx_lw_diag

    ! Purpose: Compute aerosol transmissions needed in absorptivity/
    ! emissivity calculations

    !Intent-ins
    logical,          intent(in) :: is_cmip6_volc    ! flag for using cmip6 style volc emissions
    real(r8),         intent(in) :: dt               ! time step[s]
    integer,          intent(in) :: lchnk            ! number of chunks
    integer,          intent(in) :: ncol             ! number of columns
    real(r8),         intent(in) :: pmid(:,:)        ! midpoint pressure [Pa]
    real(r8),         intent(in) :: pint(:,:)        ! interface pressure [Pa]
    real(r8),         intent(in) :: temperature(:,:) ! temperature [K]
    real(r8),         intent(in) :: zm(:,:)          ! geopotential height above surface at midpoints [m]
    real(r8),         intent(in) :: zi(:,:)          ! geopotential height above surface at interfaces [m]
    real(r8), target, intent(in) :: state_q(:,:,:)
    real(r8),         intent(in) :: pdel(:,:)
    real(r8),         intent(in) :: pdeldry(:,:)
    real(r8),         intent(in) :: cldn(:,:)
    real(r8),         intent(in) :: ext_cmip6_lw(:,:,:) !long wave extinction in the units of [1/km] 
    type(ptr2d_t),    intent(inout) :: qqcw(:)   ! Cloud borne aerosols mixing ratios [kg/kg or 1/kg]

    !intent-outs
    real(r8),         intent(out) :: odap_aer(pcols,pver,nlwbands) ! [fraction] absorption optical depth, per layer [unitless]

    ! Local variables
    !For cmip6 volcanic file
    integer  :: trop_level(pcols)
    real(r8) :: ext_cmip6_lw_inv_m(pcols,pver,nlwbands)!long wave extinction in the units of [1/m]
    !-----------------------------------------------------------------------------
#include "../../chemistry/yaml/aer_rad_props/f90_yaml/aer_rad_props_lw_beg_yml.f90"

    !Compute contributions from the modal aerosols.
    call modal_aero_lw(dt, lchnk, ncol, state_q, temperature, pmid, pdel, pdeldry, cldn, &! in
         qqcw, odap_aer) !inout/out

    !write out ext from the volcanic input file
    call outfld('extinct_lw_inp',ext_cmip6_lw(:,:,idx_lw_diag), pcols, lchnk)

    !convert from 1/km to 1/m
    ext_cmip6_lw_inv_m = ext_cmip6_lw * km_inv_to_m_inv

    !Find tropopause or quit simulation if not found
    trop_level(1:pcols) = tropopause_or_quit(lchnk, ncol, pmid, pint, temperature, zm, zi)

    !We are here because tropopause is found, update taus with 50% contributuions from the volcanic input
    !file and 50% from the existing model computed values at the tropopause layer
    call compute_odap_volcanic_at_troplayer_lw(ncol, trop_level, zi, ext_cmip6_lw_inv_m, & !in
         odap_aer) !inout

    !Above the tropopause, the read in values from the file include both the stratospheric
    !and volcanic aerosols. Therefore, we need to zero out odap_aer above the tropopause
    !and populate it exclusively from the read in values.
    call compute_odap_volcanic_above_troplayer_lw(pver, ncol, trop_level, zi, ext_cmip6_lw_inv_m, & !in
         odap_aer) !inout

    call outfld('extinct_lw_bnd7',odap_aer(:,:,idx_lw_diag), pcols, lchnk)

#include "../../chemistry/yaml/aer_rad_props/f90_yaml/aer_rad_props_lw_end_yml.f90"
  end subroutine aer_rad_props_lw

  !=============================================================================================
  subroutine compute_odap_volcanic_at_troplayer_lw(ncol, trop_level, zi, ext_cmip6_lw_inv_m, & !in
       odap_aer) !inout

    !Update odap_aer with a combination read in volcanic aerosol extinction [1/m] (50%)
    !and module computed values (50%).

    !intent-ins
    integer, intent(in) :: ncol
    integer, intent(in) :: trop_level(:)

    real(r8), intent(in) :: zi(:,:) !geopotential height above surface at interfaces [m]
    real(r8), intent(in) :: ext_cmip6_lw_inv_m(:,:,:) !long wave extinction in the units of [1/m]

    !intent-inouts
    real(r8), intent(inout) :: odap_aer(:,:,:)  ! [fraction] absorption optical depth, per layer [unitless]

    !local
    integer :: icol, ilev_tropp
    real(r8) :: lyr_thk !layer thickness [m]
#include "../../chemistry/yaml/aer_rad_props/f90_yaml/compute_odap_volcanic_at_troplayer_lw_beg_yml.f90"

    do icol = 1, ncol
       ilev_tropp = trop_level(icol) !tropopause level
       lyr_thk    = zi(icol,ilev_tropp) - zi(icol,ilev_tropp+1)! compute layer thickness in meters

       !update taus with 50% contributuions from the volcanic input file
       !and 50% from the existing model computed values at the tropopause layer
       odap_aer(icol,ilev_tropp,:) = 0.5_r8*( odap_aer(icol,ilev_tropp,:) + (lyr_thk * ext_cmip6_lw_inv_m(icol,ilev_tropp,:)) )
    enddo
#include "../../chemistry/yaml/aer_rad_props/f90_yaml/compute_odap_volcanic_at_troplayer_lw_end_yml.f90"
  end subroutine compute_odap_volcanic_at_troplayer_lw

  !==============================================================================
  subroutine compute_odap_volcanic_above_troplayer_lw(pver, ncol, trop_level, zi, ext_cmip6_lw_inv_m, & !in
       odap_aer) !inout

    !Above the tropopause, the read in values from the file include both the stratospheric
    !and volcanic aerosols. Therefore, we need to zero out odap_aer above the tropopause
    !and populate it exclusively from the read in values.

    !intent-ins
    integer, intent(in) :: pver, ncol
    integer, intent(in) :: trop_level(:)

    real(r8), intent(in) :: zi(:,:) !geopotential height above surface at interfaces [m]
    real(r8), intent(in) :: ext_cmip6_lw_inv_m(:,:,:) !long wave extinction in the units of [1/m]

    !intent-inouts
    real(r8), intent(inout) :: odap_aer(:,:,:) ! [fraction] absorption optical depth, per layer [unitless]

    !local
    integer :: ipver, icol, ilev_tropp
    real(r8) :: lyr_thk !layer thickness [m]
#include "../../chemistry/yaml/aer_rad_props/f90_yaml/compute_odap_volcanic_above_troplayer_lw_beg_yml.f90"

    !As it will be more efficient for FORTRAN to loop over levels and then columns, the following loops
    !are nested keeping that in mind
    do ipver = 1 , pver
       do icol = 1, ncol
          ilev_tropp = trop_level(icol) !tropopause level
          if (ipver < ilev_tropp) then
             lyr_thk = zi(icol,ipver) - zi(icol,ipver+1) ! compute layer thickness in meters
             odap_aer(icol,ipver,:) = lyr_thk * ext_cmip6_lw_inv_m(icol,ipver,:)
          endif
       enddo
    enddo
#include "../../chemistry/yaml/aer_rad_props/f90_yaml/compute_odap_volcanic_above_troplayer_lw_end_yml.f90"
  end subroutine compute_odap_volcanic_above_troplayer_lw

  !==============================================================================
  function tropopause_or_quit (lchnk, ncol, pmid, pint, temperature, zm, zi) result (trop_level)

    use tropopause,           only: tropopause_find
    use cam_logfile,          only: iulog
    use cam_abortutils,       only: endrun

    !Find tropopause or quit the simulation if not found

    integer, intent(in)  :: lchnk            ! number of chunks
    integer, intent(in)  :: ncol             ! number of columns
    real(r8), intent(in) :: pmid(:,:)        ! midpoint pressure [Pa]
    real(r8), intent(in) :: pint(:,:)        ! interface pressure [Pa]
    real(r8), intent(in) :: temperature(:,:) ! temperature [K]
    real(r8), intent(in) :: zm(:,:)          ! geopotential height above surface at midpoints [m]
    real(r8), intent(in) :: zi(:,:)          ! geopotential height above surface at interfaces [m]

    !return value [out]
    integer  :: trop_level(pcols) !return value

    !Local
    integer :: icol

    !trop_level has a value for tropopause for each column
    call tropopause_find(lchnk, ncol, pmid, pint, temperature, zm, zi, & !in
         trop_level) !out

    !Quit if tropopause is not found
    if (any(trop_level(1:ncol) == -1)) then
       do icol = 1, ncol
          write(iulog,*)'tropopause level,lchnk,column:',trop_level(icol),lchnk,icol
       enddo
       call endrun('aer_rad_props: tropopause not found')
    endif
  end function tropopause_or_quit
  !==============================================================================

  subroutine volcanic_cmip_sw (ncol, zi, trop_level, ext_cmip6_sw_inv_m, ssa_cmip6_sw, af_cmip6_sw, & ! in
       tau, tau_w, tau_w_g, tau_w_f)  ! inout

    !Intent-in
    integer,  intent(in) :: ncol       ! Number of columns
    real(r8), intent(in) :: zi(:,:)    ! Height above surface at interfaces [m]
    integer,  intent(in) :: trop_level(pcols)  ! tropopause level index
    real(r8), intent(in) :: ext_cmip6_sw_inv_m(pcols,pver,nswbands)  ! short wave extinction [m^{-1}]
    real(r8), intent(in) :: ssa_cmip6_sw(:,:,:),af_cmip6_sw(:,:,:)

    !Intent-inout
    real(r8), intent(inout) :: tau    (pcols,0:pver,nswbands) ! aerosol extinction optical depth
    real(r8), intent(inout) :: tau_w  (pcols,0:pver,nswbands) ! aerosol single scattering albedo * tau
    real(r8), intent(inout) :: tau_w_g(pcols,0:pver,nswbands) ! aerosol assymetry parameter * tau * w
    real(r8), intent(inout) :: tau_w_f(pcols,0:pver,nswbands) ! aerosol forward scattered fraction * tau * w

    !Local variables
    integer   :: icol, ipver, ilev_tropp
    real(r8)  :: lyr_thk ! thickness between level interfaces [m]
    real(r8)  :: ext_unitless(nswbands), asym_unitless(nswbands)
    real(r8)  :: ext_ssa(nswbands),ext_ssa_asym(nswbands)
#include "../../chemistry/yaml/aer_rad_props/f90_yaml/volcanic_cmip_sw_beg_yml.f90"

    !Logic:
    !Update taus, tau_w, tau_w_g and tau_w_f with the read in volcanic
    !aerosol extinction (1/km), single scattering albedo and asymmtry factors.

    !Above the tropopause, the read in values from the file include both the stratospheric
    !and volcanic aerosols. Therefore, we need to zero out taus above the tropopause
    !and populate them exclusively from the read in values.

    !If tropopause is found, update taus with 50% contributuions from the volcanic input
    !file and 50% from the existing model computed values

    !First handle the case of tropopause layer itself:
    do icol = 1, ncol
       ilev_tropp = trop_level(icol) !tropopause level

       lyr_thk = zi(icol,ilev_tropp) - zi(icol,ilev_tropp+1)

       ext_unitless(:)  = lyr_thk * ext_cmip6_sw_inv_m(icol,ilev_tropp,:)
       asym_unitless(:) = af_cmip6_sw (icol,ilev_tropp,:)
       ext_ssa(:)       = ext_unitless(:) * ssa_cmip6_sw(icol,ilev_tropp,:)
       ext_ssa_asym(:)  = ext_ssa(:) * asym_unitless(:)

       tau    (icol,ilev_tropp,:) = 0.5_r8 * ( tau    (icol,ilev_tropp,:) + ext_unitless(:) )
       tau_w  (icol,ilev_tropp,:) = 0.5_r8 * ( tau_w  (icol,ilev_tropp,:) + ext_ssa(:))
       tau_w_g(icol,ilev_tropp,:) = 0.5_r8 * ( tau_w_g(icol,ilev_tropp,:) + ext_ssa_asym(:))
       tau_w_f(icol,ilev_tropp,:) = 0.5_r8 * ( tau_w_f(icol,ilev_tropp,:) + ext_ssa_asym(:) * asym_unitless(:))
    enddo

    !As it will be more efficient for FORTRAN to loop over levels and then columns, the following loops
    !are nested keeping that in mind
    do ipver = 1 , pver
       do icol = 1, ncol
          ilev_tropp = trop_level(icol) !tropopause level
          if (ipver < ilev_tropp) then !BALLI: see if this is right!

             lyr_thk = zi(icol,ipver) - zi(icol,ipver+1)

             ext_unitless(:)  = lyr_thk * ext_cmip6_sw_inv_m(icol,ipver,:)
             asym_unitless(:) = af_cmip6_sw(icol,ipver,:)
             ext_ssa(:)       = ext_unitless(:) * ssa_cmip6_sw(icol,ipver,:)
             ext_ssa_asym(:)  = ext_ssa(:) * asym_unitless(:)

             tau    (icol,ipver,:) = ext_unitless(:)
             tau_w  (icol,ipver,:) = ext_ssa(:)
             tau_w_g(icol,ipver,:) = ext_ssa_asym(:)
             tau_w_f(icol,ipver,:) = ext_ssa_asym(:) * asym_unitless(:)
          endif
       enddo
    enddo

#include "../../chemistry/yaml/aer_rad_props/f90_yaml/volcanic_cmip_sw_end_yml.f90"
  end subroutine volcanic_cmip_sw

  !==============================================================================

  subroutine aer_vis_diag_out(lchnk, ncol, nnite, idxnite, tau)

    use cam_history_support, only : fillvalue
    ! output aerosol optical depth for the visible band

    integer,          intent(in) :: lchnk
    integer,          intent(in) :: ncol           ! number of columns
    integer,          intent(in) :: nnite          ! number of night columns
    integer,          intent(in) :: idxnite(:)     ! local column indices of night columns
    real(r8),         intent(in) :: tau(:,:)       ! aerosol optical depth for the visible band [unitless]

    ! Local variables
    integer  :: ii
    real(r8) :: tmp(pcols)
    !-----------------------------------------------------------------------------
#include "../../chemistry/yaml/aer_rad_props/f90_yaml/aer_vis_diag_out_beg_yml.f90"

    ! currently only implemented for climate calc

    ! compute total column aerosol optical depth
    tmp(1:ncol) = sum(tau(1:ncol,:), 2)
    ! use fillvalue to indicate night columns
    do ii = 1, nnite
       tmp(idxnite(ii)) = fillvalue
    enddo

    call outfld('AEROD_v', tmp, pcols, lchnk)
#include "../../chemistry/yaml/aer_rad_props/f90_yaml/aer_vis_diag_out_end_yml.f90"

  end subroutine aer_vis_diag_out

  !==============================================================================

end module aer_rad_props
