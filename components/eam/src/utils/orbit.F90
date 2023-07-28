module orbit
#include "../chemistry/yaml/common_files/common_uses.ymlf90"

contains

subroutine zenith(calday, clat, clon, & ! in
                  coszrs,             & ! out
                  ncol,               & ! in
                  dt_avg ) ! optional in
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute cosine of solar zenith angle for albedo and radiation
!   computations.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: J. Kiehl
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use shr_orb_mod,  only: shr_orb_decl, shr_orb_cosz
   use cam_control_mod, only: lambm0, obliqr, eccen, mvelpp
   implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer,  intent(in) :: ncol                ! number of positions
   real(r8), intent(in) :: calday              ! Calendar day, including fraction
   real(r8), intent(in) :: clat(ncol)          ! Current centered latitude (radians)
   real(r8), intent(in) :: clon(ncol)          ! Centered longitude (radians)
   real(r8), intent(in), optional :: dt_avg    ! if present, time step to use for the shr_orb_cosz calculation
!
! Output arguments
!
   real(r8), intent(out) :: coszrs(ncol)       ! Cosine solar zenith angle
!
!---------------------------Local variables-----------------------------
!
   integer  icol     ! Position loop index
   real(r8) delta    ! Solar declination angle  in radians
   real(r8) eccf     ! Earth orbit eccentricity factor
#include "../chemistry/yaml/orbit/f90_yaml/zenith_beg_yml.f90"
!
!-----------------------------------------------------------------------
!
   call shr_orb_decl (calday  ,eccen     ,mvelpp  ,lambm0  ,obliqr  , & ! in
                      delta   ,eccf      )      ! out
!
! Compute local cosine solar zenith angle,
!
   do icol=1,ncol
     coszrs(icol) = shr_orb_cosz( calday, clat(icol), clon(icol), delta, dt_avg )
   enddo
#include "../chemistry/yaml/orbit/f90_yaml/zenith_end_yml.f90"

end subroutine zenith
end module orbit
