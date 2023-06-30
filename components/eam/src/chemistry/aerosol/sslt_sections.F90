!===============================================================================
! used to compute sea salt surface emissions for modal and sectional aerosol models
!===============================================================================
module sslt_sections
  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none

  private

  public :: sslt_sections_init
  public :: fluxes
  public :: nsections
  public :: Dg
  public :: rdry

  integer, parameter :: nsections = 31
  integer, parameter :: sec1_beg = 1
  integer, parameter :: sec1_end = 9
  integer, parameter :: sec2_beg = 10
  integer, parameter :: sec2_end = 13
  integer, parameter :: sec3_beg = 14
  integer, parameter :: sec3_end = 21
  integer, parameter :: sec4_beg = 22
  integer, parameter :: sec4_end = nsections

  ! only use up to ~20um
  real(r8),parameter :: Dg(nsections) = (/  &
       0.0020e-5_r8, 0.0025e-5_r8, 0.0032e-5_r8,  &
       0.0040e-5_r8, 0.0051e-5_r8, 0.0065e-5_r8,  &
       0.0082e-5_r8, 0.0104e-5_r8, 0.0132e-5_r8,  &
       0.0167e-5_r8, 0.0211e-5_r8, 0.0267e-5_r8,  &
       0.0338e-5_r8, 0.0428e-5_r8, 0.0541e-5_r8,  &
       0.0685e-5_r8, 0.0867e-5_r8, 0.1098e-5_r8,  &
       0.1389e-5_r8, 0.1759e-5_r8, 0.2226e-5_r8,  &
       0.2818e-5_r8, 0.3571e-5_r8, 0.4526e-5_r8,  &
       0.5735e-5_r8, 0.7267e-5_r8, 0.9208e-5_r8,  &
       1.1668e-5_r8, 1.4786e-5_r8, 1.8736e-5_r8,  &
       2.3742e-5_r8 /)

  ! use Ekmam's ss
  real(r8),parameter :: rdry(1:nsections)=Dg(1:nsections)/2._r8          ! radius [m]

  ! multiply rm with 1.814 because it should be RH=80% and not dry particles
  ! for the parameterization
  real(r8),parameter :: rm(1:nsections)=1.814_r8*rdry(1:nsections)*1.e6_r8   ! um
  real(r8),parameter :: bm(1:nsections)=(0.380_r8-log10(rm(1:nsections)))/0.65_r8  ! use in Manahan

  real(r8), dimension(4,nsections) :: consta, constb  !constants for calculating emission polynomial

contains
  
  !===========================================================================
  !===========================================================================
  subroutine sslt_sections_init()

    integer :: isec

    ! use Ekman's ss
   
    ! multiply rm with 1.814 because it should be RH=80% and not dry particles
    ! for the parameterization
   

    ! calculate constants form emission polynomials
    do isec = sec1_beg, sec1_end
       consta(1,isec) = (-2.576_r8)*10._r8**35 * Dg(isec)**4 + &
                           5.932_r8*10._r8**28 * Dg(isec)**3 + &
                        (-2.867_r8)*10._r8**21 * Dg(isec)**2 + &
                        (-3.003_r8)*10._r8**13 * Dg(isec) + &
                        (-2.881_r8)*10._r8**6

       constb(1,isec) = 7.188_r8*10._r8**37 * Dg(isec)**4 + &
                     (-1.616_r8)*10._r8**31 * Dg(isec)**3 + &
                        6.791_r8*10._r8**23 * Dg(isec)**2 + &
                        1.829_r8*10._r8**16 * Dg(isec) + &
                        7.609_r8*10._r8**8
    enddo

    do isec = sec2_beg, sec2_end
       consta(2,isec) = (-2.452_r8)*10._r8**33 * Dg(isec)**4 + &
                           2.404_r8*10._r8**27 * Dg(isec)**3 + &
                        (-8.148_r8)*10._r8**20 * Dg(isec)**2 + &
                         (1.183_r8)*10._r8**14 * Dg(isec) + &
                        (-6.743_r8)*10._r8**6

       constb(2,isec) = 7.368_r8*10._r8**35 * Dg(isec)**4 + &
                     (-7.310_r8)*10._r8**29 * Dg(isec)**3 + &
                        2.528_r8*10._r8**23 * Dg(isec)**2 + &
                     (-3.787_r8)*10._r8**16 * Dg(isec)+ &
                        2.279_r8*10._r8**9
    enddo

    do isec = sec3_beg, sec3_end
       consta(3,isec) = (1.085_r8)*10._r8**29 * Dg(isec)**4 + &
                       (-9.841_r8)*10._r8**23 * Dg(isec)**3 + &
                        (3.132_r8)*10._r8**18 * Dg(isec)**2 + &
                       (-4.165_r8)*10._r8**12 * Dg(isec) + &
                        (2.181_r8)*10._r8**6

       constb(3,isec) = (-2.859_r8)*10._r8**31 * Dg(isec)**4 + &
                         (2.601_r8)*10._r8**26 * Dg(isec)**3 + &
                        (-8.297_r8)*10._r8**20 * Dg(isec)**2 + &
                         (1.105_r8)*10._r8**15 * Dg(isec) + &
                        (-5.800_r8)*10._r8**8
    enddo

    do isec = sec4_beg, sec4_end
       ! use monahan
          consta(4,isec) = (1.373_r8*rm(isec)**(-3)*(1+0.057_r8*rm(isec)**1.05_r8)  &
                          * 10**(1.19_r8*exp(-bm(isec)**2)))  &
                          * (rm(isec)-rm(isec-1))
    enddo

  end subroutine sslt_sections_init

  !===========================================================================
  !===========================================================================
  function fluxes ( sst, u10cubed, ncol ) result(fi)

    real (r8),intent(in) :: sst(:)
    real (r8),intent(in) :: u10cubed(:)
    integer  ,intent(in) :: ncol

    real (r8) :: fi(ncol,nsections)

    integer :: isec
    real (r8) :: warea(ncol)

    ! Calculations of source strength and size distribution
    ! NB the 0.1 is the dlogDp we have to multiplie with to get the flux, but the value dependence
    ! of course on what dlogDp you have. You will also have to change the sections of Dg if you use
    ! a different number of size bins with different intervals.

    warea(:ncol)=3.84e-6_r8*u10cubed(:ncol)*0.1_r8 ! whitecap area

    ! calculate number flux fi (#/m2/s)
    fi(:,:)=0._r8
    do isec = sec1_beg, sec1_end
       fi(:ncol,isec)=warea(:ncol)*((sst(:ncol))*consta(1,isec)+constb(1,isec))
    enddo

    do isec = sec2_beg, sec2_end
       fi(:ncol,isec)=warea(:ncol)*((sst(:ncol))*consta(2,isec)+constb(2,isec))
    enddo

    do isec = sec3_beg, sec3_end
       fi(:ncol,isec)=warea(:ncol)*((sst(:ncol))*consta(3,isec)+constb(3,isec))
    enddo

    do isec = sec4_beg, sec4_end
       fi(:ncol,isec)=consta(4,isec)*u10cubed(:ncol)
    enddo

  end function fluxes

end module sslt_sections
