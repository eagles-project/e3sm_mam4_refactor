!----------------------------------------------------------------------------------
! low level utility module for cloud aerosols
!
! Created by Francis Vitt
!----------------------------------------------------------------------------------
module cldaero_mod
#include "../yaml/common_files/common_uses.ymlf90"

  use shr_kind_mod, only : r8 => shr_kind_r8
  use ppgrid,       only : pcols, pver

  implicit none
  private

  public :: cldaero_uptakerate
  public :: cldaero_conc_t
  public :: cldaero_allocate
  public :: cldaero_deallocate

  type cldaero_conc_t
     real(r8), pointer :: so4c(:,:)
     real(r8), pointer :: nh4c(:,:)
     real(r8), pointer :: no3c(:,:)
     real(r8), pointer :: xlwc(:,:)
     real(r8) :: so4_fact
  end type cldaero_conc_t

contains

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
  function cldaero_allocate( ) result( cldconc )
    type(cldaero_conc_t), pointer:: cldconc

    allocate( cldconc )
    allocate( cldconc%so4c(pcols,pver) )
    allocate( cldconc%nh4c(pcols,pver) )
    allocate( cldconc%no3c(pcols,pver) )
    allocate( cldconc%xlwc(pcols,pver) )

    cldconc%so4c(:,:) = 0._r8
    cldconc%nh4c(:,:) = 0._r8
    cldconc%no3c(:,:) = 0._r8
    cldconc%xlwc(:,:) = 0._r8
    cldconc%so4_fact  = 2._r8

  end function cldaero_allocate

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
  subroutine cldaero_deallocate( cldconc )
    type(cldaero_conc_t), pointer :: cldconc

    if ( associated(cldconc%so4c) ) then
       deallocate(cldconc%so4c)
       nullify(cldconc%so4c)
    endif

    if ( associated(cldconc%nh4c) ) then
       deallocate(cldconc%nh4c)
       nullify(cldconc%nh4c)
    endif
    
    if ( associated(cldconc%no3c) ) then
       deallocate(cldconc%no3c)
       nullify(cldconc%no3c)
    endif

    if ( associated(cldconc%xlwc) ) then
       deallocate(cldconc%xlwc)
       nullify(cldconc%xlwc)
    endif

    deallocate( cldconc )
    nullify( cldconc )

  end subroutine cldaero_deallocate

!----------------------------------------------------------------------------------
! utility function for cloud-borne aerosols
!----------------------------------------------------------------------------------

!==================================================================================
  function cldaero_uptakerate( xl, cldnum, cfact, cldfrc, tfld,  press ) result( uptkrate )
!-----------------------------------------------------------------------
! compute uptake of h2so4 and msa to cloud water
!
! first-order uptake rate is
! 4*pi*(drop radius)*(drop number conc)
! *(gas diffusivity)*(fuchs sutugin correction)
!-----------------------------------------------------------------------

    use mo_constants, only : pi

    ! input arguments
    real(r8), intent(in) :: xl          ! liquid water volume [cm^3/cm^3]
    real(r8), intent(in) :: cldnum      ! droplet number concentration [#/kg]
    real(r8), intent(in) :: cfact       ! total atms density [kg/L]
    real(r8), intent(in) :: cldfrc      ! cloud fraction [fraction]
    real(r8), intent(in) :: tfld        ! temperature [K]
    real(r8), intent(in) :: press       ! pressure [Pa]
    ! output arguments
    real(r8) :: uptkrate        ! uptake rate [cm/cm/s]
    ! local variables
    real(r8) :: rad_cd          ! droplet radius [cm]
    real(r8) :: radxnum_cd      ! radius * number conc [cm/cm^3]
    real(r8) :: num_cd          ! droplet number conc [1/cm^3]
    real(r8) :: gasdiffus       ! H2SO4 gas diffusivity [cm^2/s]
    real(r8) :: gasspeed        ! H2SO4 gas mean molecular speed [cm/s]
    real(r8) :: knudsen         ! knudsen number [unitless]
    real(r8) :: fuchs_sutugin   ! another dimensionless number
    real(r8) :: volx34pi_cd     ! droplet volume * 3/4*pi [cm^3/cm^3]

    real(r8),parameter :: one_third = 0.3333333_r8      ! 1/3
    real(r8),parameter :: three_forth = 0.75_r8         ! 3/4
    real(r8),parameter :: pi4 = 12.56637_r8             ! pi*4
    real(r8),parameter :: cm3_to_L = 1.0e-3_r8          ! conversion factor from cm^3 to L (or from 1/L to 1/cm^3)
    ! artificial thresholds that assumes (radxnum_cd/volx34pi_cd < min)
    ! and (radxnum_cd/volx34pi_cd > max) as unphysical
    real(r8),parameter :: min_factor_volx34pi_radxnum = 4.0e4_r8
    real(r8),parameter :: max_factor_volx34pi_radxnum = 4.0e8_r8
#include "../yaml/cldaero_mod/f90_yaml/cldaero_uptakerate_beg_yml.f90"


! change drop number conc from #/kg to #/cm^3
        num_cd = cm3_to_L*cldnum*cfact/cldfrc
        num_cd = max( num_cd, 0.0_r8 )

! volx34pi_cd = (3/4*pi) * (liquid water volume in cm^3/cm^3)
        volx34pi_cd = xl*three_forth/pi
! radxnum_cd = (drop radius)*(drop number conc)
! following holds because volx34pi_cd = num_cd*(rad_cd**3)
        radxnum_cd = (volx34pi_cd*num_cd*num_cd)**one_third

! rad_cd = (drop radius in cm), computed from liquid water and drop number,
! then bounded by 0.5 and 50.0 micrometers to avoid the occasional unphysical value
        if (radxnum_cd <= volx34pi_cd*min_factor_volx34pi_radxnum) then
            radxnum_cd = volx34pi_cd*min_factor_volx34pi_radxnum
            rad_cd = 50.0e-4_r8
        elseif (radxnum_cd >= volx34pi_cd*max_factor_volx34pi_radxnum) then
            radxnum_cd = volx34pi_cd*max_factor_volx34pi_radxnum
            rad_cd = 0.5e-4_r8
        else
            rad_cd = radxnum_cd/num_cd
        endif

! gasdiffus = h2so4 gas diffusivity from mosaic code (cm^2/s)
! (pmid must be Pa)
        gasdiffus = 0.557_r8 * (tfld**1.75_r8) / press

! gasspeed = h2so4 gas mean molecular speed from mosaic code (cm/s)
        gasspeed = 1.455e4_r8 * sqrt(tfld/98.0_r8)

! knudsen number
        knudsen = 3.0_r8*gasdiffus/(gasspeed*rad_cd)

! following assumes accomodation coefficient accom=0.65
! (Adams & Seinfeld, 2002, JGR, and references therein)
! fuchs_sutugin = (0.75*accom*(1. + knudsen)) /
! (knudsen*(1.0 + knudsen + 0.283*accom) + 0.75*accom)
        fuchs_sutugin = (0.4875_r8*(1._r8 + knudsen)) / &
                        (knudsen*(1.184_r8 + knudsen) + 0.4875_r8)

! instantaneous uptake rate
        uptkrate = pi4*radxnum_cd*gasdiffus*fuchs_sutugin

#include "../yaml/cldaero_mod/f90_yaml/cldaero_uptakerate_end_yml.f90"
  end function cldaero_uptakerate
!==================================================================================

end module cldaero_mod
