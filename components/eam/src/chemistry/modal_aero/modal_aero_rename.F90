module modal_aero_rename
#include "../yaml/common_files/common_uses.ymlf90"
  use shr_kind_mod,    only:  r8 => shr_kind_r8
  use modal_aero_amicphys_control, only: naer, max_mode, max_aer, mass_2_vol
  use modal_aero_data, only: alnsg_amode, dgnumlo_amode, dgnumhi_amode, &
       dgnum_amode
  use physconst,       only:  pi
  use shr_log_mod ,    only: errMsg => shr_log_errMsg

  implicit none

  private

  public:: mam_rename_1subarea


  !NOTE:dryvol_smallest is a very small molar mixing ratio [m3-spc/kmol-air] (where m3-spc
  !is meter cubed volume of a species "spc") used for avoiding overflow.  it corresponds to dp = 1 nm
  !and number = 1e-5 #/mg-air ~= 1e-5 #/cm3-air
  real(r8), parameter :: smallest_dryvol_value = 1.0e-25

contains
  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine mam_rename_1subarea(iscldy, dest_mode_of_mode, nmode, &
       qnum_cur, qaer_cur, qaer_del_grow4rnam, qnumcw_cur,         &
       qaercw_cur,        qaercw_del_grow4rnam                    )

    use modal_aero_data, only: ntot_amode

    implicit none
    !-------------------------------------------------------------------------------------
    ! DESCRIPTION:
    !-------------------------------------------------------------------------------------
    ! Computes TMR (tracer mixing ratio) tendencies for "mode renaming" (i.e. transferring
    ! particles from one mode to another or "renaming" mode of the particles) during a
    ! continuous growth process.
    !
    ! Currently this transfers number and mass (and surface) from the aitken to accumulation
    ! mode after gas condensation or stratiform-cloud aqueous chemistry
    ! (convective cloud aqueous chemistry not yet implemented)
    !-------------------------------------------------------------------------------------

    !input
    logical,  intent(in) :: iscldy                      !true if sub-area is cloudy
    integer,  intent(in) :: dest_mode_of_mode(max_mode) !destination mode of a mode
    integer,  intent(in) :: nmode                       !number of modes
    real(r8), intent(in) :: qaer_del_grow4rnam(1:max_aer, 1:max_mode)            !growth in aerosol molar mixing ratio [kmol/kmol-air]
    real(r8), intent(in), optional :: qaercw_del_grow4rnam(1:max_aer, 1:max_mode)!growth in aerosol molar mixing ratio (cld borne) [kmol/kmol-air]

    !output
    real(r8), intent(inout) :: qnum_cur(1:max_mode)           !aerosol number mixing ratio [#/kmol-air]
    real(r8), intent(inout) :: qaer_cur(1:max_aer, 1:max_mode)!aerosol molar mixing ratio   [kmol/kmol-air]
    real(r8), intent(inout), optional :: qnumcw_cur(1:max_mode)           !aerosol number mixing ratio (cld borne) [#/kmol-air]
    real(r8), intent(inout), optional :: qaercw_cur(1:max_aer, 1:max_mode)!aerosol molar mixing ratio (cld borne) [kmol/kmol-air]

    ! local variables
    integer :: npair !number of pairs in different modes for the transfer

    real(r8) :: deldryvol_a(ntot_amode)   !change in dry volume [m3/kmol-air]
    real(r8) :: deldryvol_c(ntot_amode)   !change in dry volume (cld borne)[m3/kmol-air]
    real(r8) :: diameter_cutoff(max_mode) !cutoff for threshold [m]
    real(r8) :: diameter_threshold(max_mode) !Threshold to decide arosol transfer (99% of cutoff) [m]
    real(r8) :: dryvol_a(ntot_amode)      !dry volume [m3/kmol-air]
    real(r8) :: dryvol_c(ntot_amode)      !dry volume (cld borne)[m3/kmol-air]
    real(r8) :: sz_factor(ntot_amode)     !size factor for each mode [unitless]
    real(r8) :: fmode_dist_tail_fac(ntot_amode) !tail distribution factor for each mode [unitless]
    real(r8) :: lndiameter_cutoff(max_mode) !log of diamter cutoff [m]
    real(r8) :: ln_diameter(max_mode)     !log of diameter [m]
    real(r8) :: v2nhirlx(ntot_amode), v2nlorlx(ntot_amode) !high and low volume to num ratios[m^-3]
#include "../yaml/modal_aero_rename/f90_yaml/mam_rename_1subarea_beg.ymlf90"

    !------------------------------------------------------------------------
    !Find mapping between different modes, so that we can move aerosol
    !particles from one mode to another
    !------------------------------------------------------------------------

    !FIXME: All the arrays in find_renaming_pairs subroutine call should be
    !initialized to HUGE or NaNs as they are partially populated

    !Find (src->destination) pairs of modes (e.g., if only mode #1 and mode #2 can participate in the
    !inter-mode transfer, number of pairs will be 1 and so on) which can participate in
    !inter-mode species transfer

    call find_renaming_pairs (ntot_amode, dest_mode_of_mode, & !input
         npair, sz_factor, fmode_dist_tail_fac, v2nlorlx, &    !output
         v2nhirlx, ln_diameter, diameter_cutoff, &             !output
         lndiameter_cutoff, diameter_threshold)                !output

    if (npair <= 0) return ! if no transfer required, return

    !Interstitial aerosols: Compute initial (before growth) aerosol dry volume and
    !also the growth in dryvolume of the "src" mode
    call compute_dryvol_change_in_src_mode(ntot_amode, naer, dest_mode_of_mode,  & !input
         qaer_cur, qaer_del_grow4rnam, & !input
         dryvol_a, deldryvol_a                 ) !output

    !Cloudborne aerosols: Compute initial (before growth) aerosol dry volume and
    !also the growth in dryvolume of the "src" mode
    if (iscldy) then
       call compute_dryvol_change_in_src_mode(ntot_amode, naer, dest_mode_of_mode,  & !input
            qaercw_cur, qaercw_del_grow4rnam, & !input
            dryvol_c, deldryvol_c                     ) !output
    endif

    !Find fractions (mass and number) to transfer and complete the transfer
    call do_inter_mode_transfer(ntot_amode, naer, dest_mode_of_mode, &                         !input
         iscldy, v2nlorlx, v2nhirlx, dryvol_a, dryvol_c, deldryvol_a, deldryvol_c, &           !input
         sz_factor, fmode_dist_tail_fac, ln_diameter, lndiameter_cutoff, diameter_threshold, & !input
         qaer_cur, qnum_cur, qaercw_cur, qnumcw_cur ) !output
#include "../yaml/modal_aero_rename/f90_yaml/mam_rename_1subarea_end.ymlf90"
  end subroutine mam_rename_1subarea

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------

  subroutine find_renaming_pairs (nmodes, dest_mode_of_mode, &  !input
       num_pairs, sz_factor, fmode_dist_tail_fac, v2n_lo_rlx, & !output
       v2n_hi_rlx, ln_diameter_tail_fac, diameter_cutoff, &     !output
       ln_dia_cutoff, diameter_threshold)    !output

    !------------------------------------------------------------------------
    !Find number of pairs which participates in the inter-mode transfer
    !------------------------------------------------------------------------
    !Number of mode pairs allowed to do inter-mode particle transfer
    !(e.g. if we have a pair "mode_1<-->mode_2", mode_1 and mode_2 can
    !participate in the inter-mode aerosol particle transfer where particles
    !of the same species in mode_1 can be transferred to mode_2 and vice-versa)
    !------------------------------------------------------------------------

    !arguments (intent-ins)
    integer, intent(in) :: nmodes               !total number of modes
    integer, intent(in) :: dest_mode_of_mode(:) !array carry info about the destination mode of a particular mode

    !intent-outs
    integer,  intent(out) :: num_pairs          !total number of pairs to be found
    real(r8), intent(out) :: sz_factor(:), fmode_dist_tail_fac(:) !precomputed factors to be used later [unitless]
    real(r8), intent(out) :: v2n_lo_rlx(:), v2n_hi_rlx(:)         !relaxed volume to num high and low ratio limits [m^-3]
    real(r8), intent(out) :: ln_diameter_tail_fac(:)              !log of diameter factor for distribution tail [unitless]
    real(r8), intent(out) :: diameter_cutoff(:), ln_dia_cutoff(:) !cutoff (threshold) for deciding the  do inter-mode transfer [m]
    real(r8), intent(out) :: diameter_threshold(:)                !threshold diameter(99% of cutoff)[m]

    ! local variables
    integer :: dest_mode, src_mode, imode

    !some parameters
    real(r8), parameter :: sqrt_half = sqrt(0.5)
    real(r8), parameter :: frelax = 27.0_r8 !(3^3): relaxing 3*diameter, which makes it 3^3 for volume
#include "../yaml/modal_aero_rename/f90_yaml/find_renaming_pairs_beg.ymlf90"
    ! Let us assume there are none to start with
    num_pairs = 0

    !if there can be no possible pairs, just return
    if (all(dest_mode_of_mode(:)<=0)) return

    !Go through all the modes to find if we have atleast one or more than one pairs
    do imode = 1, nmodes
       dest_mode   = dest_mode_of_mode(imode) ! "destination" mode for mode "imode"

       !if dest_mode is <=0, transfer is not possible for this mode, cycle the loop for the next mode
       if(dest_mode <= 0)cycle

       src_mode = imode                  ! transfer "src" mode is the current mode (i.e. imode)

       !^^At this point, we know that particles can be transferred from the
       ! "src_mode" to "dest_mode". "src_mode" is the current mode (i.e. imode)

       !update number of pairs found so far
       num_pairs = num_pairs + 1    !increment npair

       !-------------------------------------------------------
       !Now precompute some common factors to be used later
       !-------------------------------------------------------

       call compute_size_factor (src_mode,  sz_factor) !size factor for "src mode"
       call compute_size_factor (dest_mode, sz_factor) !size factor for "dest mode"

       !---------------------------------------------------------------------------------------------------------
       ! We compute few factors below for the "src_mode", which will be used for inter-mode particle transfer
       !---------------------------------------------------------------------------------------------------------

       fmode_dist_tail_fac(src_mode) = sqrt_half/alnsg_amode(src_mode) !factor for computing distribution tails of the  "src mode"

       !compute volume to number high and low limits with relaxation coefficients (watch out for the repeated calculations)
       v2n_lo_rlx(src_mode) = vol_to_num_ratio(src_mode, dgnumlo_amode) * frelax
       v2n_hi_rlx(src_mode) = vol_to_num_ratio(src_mode, dgnumhi_amode) / frelax

       !A factor for computing diameter at the tails of the distribution
       ln_diameter_tail_fac(src_mode) = 3.0 * (alnsg_amode(src_mode)**2)

       !cutoff (based on geometric mean) for making decision to do inter-mode transfers
       !We took geommetric mean of the participating modes (source and destination)
       !to find a cutoff or threshold from moving particles from the source to the
       !destination mode.
       diameter_cutoff(src_mode) = sqrt(   &
            dgnum_amode(src_mode)*exp(1.5*(alnsg_amode(src_mode)**2)) *   &
            dgnum_amode(dest_mode)*exp(1.5*(alnsg_amode(dest_mode)**2)) )

       ln_dia_cutoff(src_mode) = log(diameter_cutoff(src_mode)) !log of cutoff
       diameter_threshold(src_mode) = 0.99*diameter_cutoff(src_mode) !99% of the cutoff

    end do
#include "../yaml/modal_aero_rename/f90_yaml/find_renaming_pairs_end.ymlf90"
  end subroutine find_renaming_pairs
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------

  subroutine compute_size_factor(imode, size_factor)

    !--------------------------------------------------
    ! Compute a common factor to be used in determining
    ! size of a particle
    !--------------------------------------------------

    implicit none

    integer,  intent(in) :: imode     !mode number
    real(r8), intent(inout) :: size_factor(:) !size factor [unitless]

    size_factor(imode) = (pi/6.)*exp(4.5*(alnsg_amode(imode)**2))

  end subroutine compute_size_factor

  !----------------------------------------------------------------------
  !----------------------------------------------------------------------

  pure function vol_to_num_ratio(imode, diameter) result(v2n)

    !Compute volume to number ratio for a mode
    !Defined as the number of particles larger than a specified, critical
    !diameter divided by the  aerosol volume

    implicit none
    integer,  intent(in) :: imode
    real(r8), intent(in) :: diameter(:) !particle diameter [m]

    real(r8) :: v2n !return value [m^-3]

    v2n = ( 1._r8 / ( (pi/6._r8)* &
         (diameter(imode)**3._r8)*exp(4.5_r8*alnsg_amode(imode)**2._r8) ) )

  end function vol_to_num_ratio

  !----------------------------------------------------------------------
  !----------------------------------------------------------------------

  subroutine compute_dryvol_change_in_src_mode(nmode, nspec, dest_mode_of_mode, & !input
       q_mmr, q_del_growth, & !input
       dryvol, deldryvol ) !output

    !--------------------------------------------------------
    !Compute the change in the dryvolume of the source mode
    !--------------------------------------------------------

    !inputs
    integer,  intent(in):: nmode !total number of modes
    integer,  intent(in):: nspec !total number of species in a mode
    integer,  intent(in):: dest_mode_of_mode(:) !destination mode for a mode

    real(r8), intent(in) :: q_mmr(:,:)           !molar mixing ratios (mmr) [kmol/kmol-air]
    real(r8), intent(in) :: q_del_growth(:,:)    !growth in mmr [kmol/kmol-air]

    !intent-outs
    real(r8), intent(out) :: dryvol(:)    !dry volumes (before growth)[m3/kmol-air]
    real(r8), intent(out) :: deldryvol(:) !change in dry volumes [m3/kmol-air]

    !local
    integer :: imode
    integer :: dest_mode
#include "../yaml/modal_aero_rename/f90_yaml/ompute_dryvol_change_in_src_mode_beg_yml.f90"

    !For each mode, compute the initial (before growth) dryvolume and the growth in dryvolume
    do imode = 1, nmode
       !compute dry volume only for modes participating in inter-modal transfer
       dest_mode = dest_mode_of_mode(imode)
       if (dest_mode <= 0) cycle

       !compute dry volumes (before growth) and its change for interstitial aerosols
       call compute_dryvolume_change(imode, nspec, q_mmr, q_del_growth, & !input
            dryvol(imode), deldryvol(imode)) !output
    end do
#include "../yaml/modal_aero_rename/f90_yaml/ompute_dryvol_change_in_src_mode_end_yml.f90"
  end subroutine compute_dryvol_change_in_src_mode

  !----------------------------------------------------------------------
  !----------------------------------------------------------------------

  subroutine compute_dryvolume_change (imode, nspec, q_mmr, q_del_growth, &!input
       dryvol, deldryvol) !output

    !--------------------------------------------------------
    !Compute dry volume and change in dry volume due to growth
    !--------------------------------------------------------

    !intent-ins
    integer,  intent(in) :: imode           !current mode number
    integer,  intent(in) :: nspec           !number of species in the current mode
    real(r8), intent(in) :: q_mmr(:,:)        !molar mixing ratio [kmol/kmol-air]
    real(r8), intent(in) :: q_del_growth(:,:) !change (delta) in molar mixing ratio [kmol/kmol-air]

    !intent-outs
    real(r8), intent(out) :: dryvol, deldryvol !dry volume (before growth) and its grwoth [m3/kmol-air]

    !local variables
    integer  :: ispec, s_spec_ind, e_spec_ind
    real(r8) :: tmp_dryvol, tmp_del_dryvol![m3/kmol-air]

    !For each mode, we compute a dry volume by combining (accumulating) mass/density for each species in that mode.
        !conversion from mass to volume is accomplished by multiplying with precomputed "mass_2_vol" factor

    s_spec_ind = 1     !start species index for this mode [These will be subroutine args]
    e_spec_ind = nspec !end species index for this mode

    !initialize tmp accumulators
    tmp_dryvol     = 0.0_r8 !dry volume accumulator
    tmp_del_dryvol = 0.0_r8 !dry volume growth(change) accumulator

    ! Notes on mass_2_vol factor:Units:[m3/kmol-species]; where kmol-species is the amount of a species "s"
    ! This factor is obtained by  (molecular_weight/density) of a specie. That is,
    ! [ (g/mol-species) / (kg-species/m3) ]; where molecular_weight has units [g/mol-species] and density units are [kg-specie/m3]
    ! which results in the units of m3/kmol-specie

    do ispec = s_spec_ind, e_spec_ind
       !Multiply by mass_2_vol[m3/kmol-species] to convert q_mmr[kmol-specie/kmol-air]) to volume units[m3/kmol-air]
       tmp_dryvol     = tmp_dryvol     + q_mmr(ispec,imode)*mass_2_vol(ispec)        !compute current dryvolume
       !accumulate the "grwoth" in volume units as well
       tmp_del_dryvol = tmp_del_dryvol + q_del_growth(ispec,imode)*mass_2_vol(ispec) !compute dryvolume growth
    end do

    dryvol    = tmp_dryvol-tmp_del_dryvol ! This is dry volume before the growth
    deldryvol = tmp_del_dryvol          ! change in dry volume due to growth

  end subroutine compute_dryvolume_change

  !----------------------------------------------------------------------
  !----------------------------------------------------------------------

  subroutine do_inter_mode_transfer(nmode, nspec, dest_mode_of_mode, &
       iscldy, v2nlorlx, v2nhirlx, dryvol_a, dryvol_c, deldryvol_a, deldryvol_c, &
       sz_factor, fmode_dist_tail_fac, ln_diameter_tail_fac, ln_dia_cutoff, diameter_threshold, &
       qaer_cur, qnum_cur, qaercw_cur, qnumcw_cur)

    !-------------------------------------------------------------------------
    !Compute how much to transfer from one mode to another and do the transfer
    !-------------------------------------------------------------------------

    !intent-ins
    integer,  intent(in) :: nmode, nspec, dest_mode_of_mode(:)
    logical,  intent(in) :: iscldy !true if it is cloudy cell
    real(r8), intent(in) :: v2nlorlx(:), v2nhirlx(:) !volume to number relaxation limits [m^-3]
    real(r8), intent(in) :: dryvol_a(:), dryvol_c(:), deldryvol_a(:), deldryvol_c(:)!dryvolume and the change in dryvolume[m3/kmol-air]
    real(r8), intent(in) :: sz_factor(:), fmode_dist_tail_fac(:) !Some precomputed factors [unitless]
    real(r8), intent(in) :: ln_diameter_tail_fac(:), ln_dia_cutoff(:), diameter_threshold(:)! diameters and diameter thresholds [m]

    !intent-inouts
    real(r8), intent(inout) :: qaer_cur(:,:)!aerosol molar mixing ratio    [kmol/kmol-air]
    real(r8), intent(inout) :: qnum_cur(:)  !aerosol number mixing ratio [#/kmol-air]

    real(r8), intent(inout), optional :: qaercw_cur(:,:) !aerosol molar mixing ratio (cld borne) [kmol/kmol-air]
    real(r8), intent(inout), optional :: qnumcw_cur(:)   !aerosol number mixing ratio (cld borne) [#/kmol-air]

    !local variables
    integer :: src_mode, dest_mode, imode, ispec
    logical :: is_xfer_frac_zero

    !before growth dryvolumes [m3/kmol-air]
    real(r8) :: bef_grwth_dryvol, dryvol_del, bef_grwth_tail_fr_vol, xfer_vol_frac, bef_grwth_dryvolbnd

    !after growth dryvolumes [m3/kmol-air]
    real(r8) :: aft_grwth_tail_fr_vol, aft_grwth_dryvol

    !before growth numbers [#/kmol-air]
    real(r8) :: bef_grwth_num, bef_grwth_numbnd, bef_grwth_tail_fr_num, xfer_num_frac
    real(r8) :: aft_grwth_tail_fr_num

    !diameters [m]
    real(r8) :: aft_grwth_diameter, ln_dia_aftgrwth, bef_grwth_diameter

    !Loop through the modes and do the transfer
    pair_loop:  do imode = 1, nmode

       src_mode = imode                     !source mode
       dest_mode = dest_mode_of_mode(imode) !destination mode

       !if destination mode doesn't exist for the source mode, cycle loop
       if (dest_mode <= 0) cycle pair_loop

       !compute before growth dry volume and number

       call compute_before_growth_dryvol_and_num(iscldy, src_mode, dryvol_a, dryvol_c,   & !input
            qnum_cur, qnumcw_cur, v2nhirlx(src_mode), v2nlorlx(src_mode),        & !input
            bef_grwth_dryvol, bef_grwth_dryvolbnd, bef_grwth_numbnd)               !output

       !change (delta) in dryvol
       dryvol_del = total_inter_cldbrn(iscldy, src_mode, deldryvol_a, deldryvol_c)

       !Total dryvolume after growth (add delta growth)
       aft_grwth_dryvol = bef_grwth_dryvol + dryvol_del

       !Skip inter-mode transfer for this mode if dry after grwoth is ~ 0
       if (aft_grwth_dryvol <= smallest_dryvol_value) cycle pair_loop

       !compute before growth diameter
       bef_grwth_diameter = mode_diameter(bef_grwth_dryvolbnd, bef_grwth_numbnd, sz_factor(src_mode))

       !if the before growth diameter is more than the threshold (diameter_threshold), we restrict diameter
       !to the threshold and change dry volume accorindgly
       if (bef_grwth_diameter > diameter_threshold(src_mode)) then
          ! this revised volume corresponds to bef_grwth_diameter == diameter_threshold, and same number conc
          bef_grwth_dryvol = bef_grwth_dryvol * (diameter_threshold(src_mode)/bef_grwth_diameter)**3
          bef_grwth_diameter = diameter_threshold(src_mode)
       end if

       if ((aft_grwth_dryvol-bef_grwth_dryvol) <= 1.0e-6_r8*bef_grwth_dryvolbnd) cycle pair_loop

       !Compute after growth diameter; if it is less than the "nominal" or "base" diameter for
       !the source mode, skip inter-mode transfer
       aft_grwth_diameter = mode_diameter(aft_grwth_dryvol,bef_grwth_numbnd,sz_factor(src_mode))

       if (aft_grwth_diameter <= dgnum_amode(src_mode)) cycle pair_loop

       !compute before growth number fraction in the tail
       call compute_tail_fraction(bef_grwth_diameter,ln_dia_cutoff(src_mode), fmode_dist_tail_fac(src_mode), & !input
            tail_fraction = bef_grwth_tail_fr_num ) !output

       !compute before growth volume (or mass) fraction in the tail
       call compute_tail_fraction(bef_grwth_diameter,ln_dia_cutoff(src_mode), fmode_dist_tail_fac(src_mode), & !input
            log_dia_tail_fac = ln_diameter_tail_fac(src_mode), & !optional input
            tail_fraction = bef_grwth_tail_fr_vol ) !output

       !compute after growth number fraction in the tail
       call compute_tail_fraction(aft_grwth_diameter,ln_dia_cutoff(src_mode), fmode_dist_tail_fac(src_mode), & !input
            tail_fraction = aft_grwth_tail_fr_num ) !output

       !compute after growth volume (or mass) fraction in the tail
       call compute_tail_fraction(aft_grwth_diameter,ln_dia_cutoff(src_mode), fmode_dist_tail_fac(src_mode), & !input
            log_dia_tail_fac = ln_diameter_tail_fac(src_mode), & !optional input
            tail_fraction = aft_grwth_tail_fr_vol ) !output

       !compute transfer fraction (volume and mass) - if less than zero, cycle loop
       call compute_xfer_fractions(bef_grwth_dryvol, aft_grwth_dryvol, bef_grwth_tail_fr_vol, aft_grwth_tail_fr_vol, & !input
            aft_grwth_tail_fr_num, bef_grwth_tail_fr_num, &
            is_xfer_frac_zero, xfer_vol_frac, xfer_num_frac) !output

       if (is_xfer_frac_zero) cycle pair_loop

       !do the transfer for the interstitial species
       call do_num_and_mass_transfer(nspec, src_mode, dest_mode, xfer_vol_frac, xfer_num_frac, & !input
            qaer_cur, qnum_cur) !output

       !do the traner for the cloud borne species
       if ( iscldy ) then
          call do_num_and_mass_transfer(nspec, src_mode, dest_mode, xfer_vol_frac, xfer_num_frac, & !input
               qaercw_cur, qnumcw_cur) !output
       end if

    enddo pair_loop
  end subroutine do_inter_mode_transfer

  !----------------------------------------------------------------------
  !----------------------------------------------------------------------

  function total_inter_cldbrn(iscldy, imode, interstitial, cldbrn) result (total)

    use cam_abortutils,  only:  endrun
    !compute total (dry volume or number) of interstitial and cloud borne species

    logical,  intent(in) :: iscldy       !TRUE, if a cell has cloud
    integer,  intent(in) :: imode
    real(r8), intent(in) :: interstitial(:)     !interstital part [unit depends on the input]
    real(r8), intent(in), optional :: cldbrn(:) !cloud borne part [unit depends on the input]

    !return value
    real(r8) :: total

    total = interstitial(imode) !if there is no cloud, total is just the interstitial value

    if(iscldy) then ! If there is cloud, add cloud borne part as well
       if(.not.present(cldbrn))then
          call endrun("If a grid cell is cloudy, cloud borne aerosol values must be present:"//errmsg(__FILE__,__LINE__))
       end if
       total = total + cldbrn(imode)
    end if
  end function total_inter_cldbrn

  !----------------------------------------------------------------------
  !----------------------------------------------------------------------

  function mode_diameter(volume, number, size_factor) result(diameter)

    !compute diameter
    real(r8), intent(in) :: volume      ![m3]
    real(r8), intent(in) :: number      ![#/kmol-air]
    real(r8), intent(in) :: size_factor ![unitless]

    !return value
    real(r8) :: diameter ![m]
    !local
    real(r8), parameter :: onethird = 1.0_r8/3.0_r8

    diameter = (volume/(number*size_factor))**onethird

  end function mode_diameter
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------

  subroutine compute_before_growth_dryvol_and_num(iscldy, src_mode, dryvol_a, dryvol_c, & !input
       qnum_cur, qnumcw_cur, v2nhi, v2nlo,                                 & !input
       bef_grwth_dryvol, bef_grwth_dryvolbnd, bef_grwth_numbnd)              !output

    use mam_support, only: min_max_bound

    !Compute before growth dry volume and number
    implicit none

    logical,  intent(in) :: iscldy       !TRUE, if a cell has cloud
    integer,  intent(in) :: src_mode

    real(r8), intent(in) :: dryvol_a(:), dryvol_c(:) ![m3/kmol-air]
    real(r8), intent(in) :: qnum_cur(:)              ![#/kmol-air]

    real(r8), intent(in), optional :: qnumcw_cur(:)  ![#/kmol-air]

    real(r8), intent(in) :: v2nlo, v2nhi             ![m^-3]

    !output
    real(r8), intent(out) :: bef_grwth_dryvol,bef_grwth_dryvolbnd ![m3/kmol-air]
    real(r8), intent(out) :: bef_grwth_numbnd ! [#/kmol-air]

    !local
    real(r8) :: bef_grwth_num ![#/kmol-air]

    !Compute total(i.e. cloud borne and interstitial) of dry volume (before growth)
    ! and delta in dry volume in the source mode [units: (m3 of specie)/(kmol of air)]
    !NOTE: cloudborne input can be optional, so we are sending "src_mode" as a argument
    !as we cannot refecence a member of an optional array if it is not present
    bef_grwth_dryvol = total_inter_cldbrn(iscldy, src_mode, dryvol_a , dryvol_c)

    bef_grwth_dryvolbnd = max( bef_grwth_dryvol, smallest_dryvol_value )

    !Compute total before growth number [units: #/kmol-air]
    bef_grwth_num    = total_inter_cldbrn(iscldy, src_mode, qnum_cur, qnumcw_cur)
    bef_grwth_num = max( 0.0_r8, bef_grwth_num ) !bound to have minimum of 0

    !bound number within min and max of the source mode
    bef_grwth_numbnd = min_max_bound(bef_grwth_dryvolbnd*v2nhi, & !min value
         bef_grwth_dryvolbnd*v2nlo, bef_grwth_num) !max value and input

  end subroutine compute_before_growth_dryvol_and_num


  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  subroutine compute_tail_fraction(diameter,log_dia_cutoff, tail_dist_fac, & !input
       log_dia_tail_fac, & !optional input
       tail_fraction ) !output

    !Compute tail fraction to be used for inter-mode species transfer

    use shr_spfn_mod, only: erfc_shr => shr_spfn_erfc  !E3SM implementation of the erro function

    implicit none

    real(r8), intent(in) :: diameter       ![m]
    real(r8), intent(in) :: log_dia_cutoff ![m]
    real(r8), intent(in) :: tail_dist_fac  ![unitless]

    real(r8), intent(in), optional :: log_dia_tail_fac! [m]

    real(r8), intent(out) :: tail_fraction ![unitless]

    real(r8) :: log_diameter, tail ![m]

    log_diameter  = log(diameter)
    if (present(log_dia_tail_fac)) log_diameter  = log_diameter + log_dia_tail_fac
    tail          = ( log_dia_cutoff - log_diameter ) * tail_dist_fac
    tail_fraction = 0.5_r8*erfc_shr( tail )

  end subroutine compute_tail_fraction

  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  subroutine compute_xfer_fractions(bef_grwth_dryvol, aft_grwth_dryvol, bef_grwth_tail_fr_vol, aft_grwth_tail_fr_vol, & !input
       aft_grwth_tail_fr_num, bef_grwth_tail_fr_num, & !input
       is_xfer_frac_zero, xfer_vol_frac, xfer_num_frac) !output

    implicit none

    !Compute number and volume to transfer for inter-mode species transfer

    !input
    real(r8), intent(in) :: bef_grwth_dryvol      ![m3/kmol-air]
    real(r8), intent(in) :: aft_grwth_dryvol      ![m3/kmol-air]
    real(r8), intent(in) :: bef_grwth_tail_fr_vol ![unitless]
    real(r8), intent(in) :: aft_grwth_tail_fr_vol ![unitless]
    real(r8), intent(in) :: aft_grwth_tail_fr_num ![#/kmol-air]
    real(r8), intent(in) :: bef_grwth_tail_fr_num ![#/kmol-air]

    !output
    logical,  intent(out) :: is_xfer_frac_zero
    real(r8), intent(out) :: xfer_vol_frac ![m3/kmol-air]
    real(r8), intent(out) :: xfer_num_frac ![#/kmol-air]

    !local
    real(r8), parameter :: xferfrac_max = 1.0_r8 - 10.0_r8*epsilon(1.0_r8) !1-eps (this number is little less than 1, e.g. 0.99)
    real(r8) :: volume_fraction

    !assume we have fractions to transfer, so we will not skip the rest of the calculations
    is_xfer_frac_zero = .false.

    !   transfer fraction is difference between new and old tail-fractions
    volume_fraction = aft_grwth_tail_fr_vol*aft_grwth_dryvol - bef_grwth_tail_fr_vol*bef_grwth_dryvol

    if (volume_fraction <= 0.0_r8) then
       is_xfer_frac_zero = .true.
       return
    endif

    xfer_vol_frac = min( volume_fraction, aft_grwth_dryvol )/aft_grwth_dryvol
    xfer_vol_frac = min( xfer_vol_frac, xferfrac_max )
    xfer_num_frac = aft_grwth_tail_fr_num - bef_grwth_tail_fr_num

    !   transfer fraction for number cannot exceed that of mass
    xfer_num_frac = max( 0.0_r8, min( xfer_num_frac, xfer_vol_frac ) )

  end subroutine compute_xfer_fractions

  !----------------------------------------------------------------------
  !----------------------------------------------------------------------

  subroutine do_num_and_mass_transfer(nspec, src_mode, dest_mode, xfer_vol_frac, xfer_num_frac, & !input
       qaer, qnum) !output

    !Transfer species from source to destination model and update the mixing ratios

    !input
    integer,  intent(in) :: nspec, src_mode, dest_mode
    real(r8), intent(in) :: xfer_vol_frac ![m3/kmol-air]
    real(r8), intent(in) :: xfer_num_frac ![#/kmol-air]

    !output
    real(r8), intent(inout) :: qaer(:,:) ![kmol/kmol-air]
    real(r8), intent(inout) :: qnum(:)   ![#/kmol-air]

    !local
    integer :: ispec
    real(r8) :: num_trans, vol_trans

    !compute changes to number and species masses
    num_trans = qnum(src_mode)*xfer_num_frac
    qnum(src_mode) = qnum(src_mode) - num_trans
    qnum(dest_mode) = qnum(dest_mode) + num_trans
    do ispec = 1, nspec
       vol_trans = qaer(ispec,src_mode)*xfer_vol_frac
       qaer(ispec,src_mode) = qaer(ispec,src_mode) - vol_trans
       qaer(ispec,dest_mode) = qaer(ispec,dest_mode) + vol_trans
    enddo

  end subroutine do_num_and_mass_transfer
end module modal_aero_rename
