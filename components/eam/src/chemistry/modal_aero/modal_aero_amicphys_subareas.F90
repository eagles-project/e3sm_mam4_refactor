module modal_aero_amicphys_subareas

  use shr_kind_mod, only: wp => shr_kind_r8

  implicit none

  public

contains

  subroutine setup_subareas( cld, maxsubarea,                         &! in
                             nsubarea, ncldy_subarea, jclea, jcldy,   &! out
                             iscldy_subarea, afracsub, fclea, fcldy )! out
  !--------------------------------------------------------------------------------------
  ! Determine the number of sub-areas and their fractional areas.
  ! Assign values to some bookkeeping variables.
  !--------------------------------------------------------------------------------------

      implicit none

      real(wp), intent(in)  :: cld
      integer,  intent(in)  :: maxsubarea 

      integer,  intent(out) :: nsubarea, ncldy_subarea
      integer,  intent(out) :: jclea, jcldy
      logical,  intent(out) :: iscldy_subarea(maxsubarea)
      real(wp), intent(out) :: afracsub(maxsubarea)
      real(wp), intent(out) :: fclea, fcldy 

      ! Cloud chemistry is only on when cld(i,k) >= 1.0e-5_wp
      ! It may be that the macrophysics has a higher threshold that this
      real(wp), parameter :: fcld_locutoff = 1.0e-5_wp

      real(wp), parameter :: fcld_hicutoff = 0.999_wp


      ! if cloud fraction ~= 0, the grid-cell has a single clear  sub-area      (nsubarea = 1)
      ! if cloud fraction ~= 1, the grid-cell has a single cloudy sub-area      (nsubarea = 1)
      ! otherwise,              the grid-cell has a clear and a cloudy sub-area (nsubarea = 2)

      if (cld < fcld_locutoff) then
         fcldy = 0.0_wp
         nsubarea = 1 ; ncldy_subarea = 0
         jclea = 1 ; jcldy = 0
      else if (cld > fcld_hicutoff) then
         fcldy = 1.0_wp
         nsubarea = 1 ; ncldy_subarea = 1
         jclea = 0 ; jcldy = 1
      else
         fcldy = cld
         nsubarea = 2 ; ncldy_subarea = 1
         jclea = 1 ; jcldy = 2
      end if
      fclea = 1.0_wp - fcldy

      ! Set up a logical array to indicate whether the subareas are clear or cloudy

      iscldy_subarea(:) = .false.
      if (jcldy > 0) iscldy_subarea(jcldy) = .true.

      ! Save the area fractions to an array

      afracsub(:) = 0.0_wp
      if (jclea > 0) afracsub(jclea) = fclea
      if (jcldy > 0) afracsub(jcldy) = fcldy

  end subroutine setup_subareas

  subroutine set_subarea_relhum( maxsubarea,ncldy_subarea,jclea,jcldy, &! in
                                 afracsub,relhumgcm,                   &! in
                                 relhumsub                             )! out

      integer,  intent(in) :: maxsubarea
      integer,  intent(in) :: ncldy_subarea
      integer,  intent(in) :: jclea, jcldy
      real(wp), intent(in) :: afracsub(maxsubarea)
      real(wp), intent(in) :: relhumgcm

      real(wp), intent(out) :: relhumsub(maxsubarea)

      real(wp) :: relhum_tmp

      if (ncldy_subarea <= 0) then
      ! Entire grid cell is cloud-free. There is only one (i.e., clear) subarea; RH = grid cell mean.
         relhumsub(:) = relhumgcm

#if ( defined( MAM_STANDALONE ) )
      else if (cldy_rh_sameas_clear > 0) then
         relhumsub(:) = relhumgcm
#endif

      else
         ! Grid cell has a cloudy subarea. Set RH in that part to 1.0.
         relhumsub(jcldy) = 1.0_wp

         ! If the grid cell also has a clear portion, back out the RH from the
         ! grid-cell mean RH and the cloud fraction.

         if (jclea > 0) then
            relhum_tmp = (relhumgcm - afracsub(jcldy))/afracsub(jclea)
            relhumsub(jclea) = max( 0.0_wp, min( 1.0_wp, relhum_tmp ) )
         end if
      end if

  end subroutine set_subarea_relhum

end module modal_aero_amicphys_subareas
