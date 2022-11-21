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

  subroutine copy_cnst( q_in, q_copy, ncnst, lcopy )

      integer, intent(in)    :: ncnst         ! # of constituents that might need to be copied
      logical, intent(in)    :: lcopy(ncnst)  ! whether individual constituents should be copied
      real(wp),intent(in)    :: q_in(ncnst)   ! values to be copied
      real(wp),intent(inout) :: q_copy(ncnst) ! copy of input values 

      where( lcopy )
        q_copy = q_in
      end where

  end subroutine copy_cnst

  subroutine set_subarea_q_numb_for_cldbrn_aerosols( loffset, maxsubarea, gas_pcnst, jclea, jcldy, fcldy, qqcwgcm, qqcwsub )

     use modal_aero_data, only: ntot_amode, nspec_amode, numptrcw_amode

     integer, intent(in)    :: loffset 
     integer, intent(in)    :: maxsubarea, gas_pcnst
     integer, intent(in)    :: jclea, jcldy
     real(wp),intent(in)    :: fcldy
     real(wp),intent(in)    :: qqcwgcm(gas_pcnst)
     real(wp),intent(inout) :: qqcwsub(gas_pcnst,maxsubarea)

     integer :: imode    ! mode index
     integer :: icnst    ! consitituent index 

     !----------------------------------------------------------------
     do imode = 1, ntot_amode  ! loop thru log-normal modes

        icnst = numptrcw_amode(imode) - loffset

        qqcwsub(icnst,jclea) = 0.0_wp
        qqcwsub(icnst,jcldy) = qqcwgcm(icnst)/fcldy

     end do
     !----------------------------------------------------------------

  end subroutine set_subarea_q_numb_for_cldbrn_aerosols

  subroutine set_subarea_q_mass_for_cldbrn_aerosols( loffset,maxsubarea, gas_pcnst, jclea, jcldy, fcldy, qqcwgcm, qqcwsub )

     use modal_aero_data, only: ntot_amode, nspec_amode, lmassptrcw_amode

     integer, intent(in)    :: loffset 
     integer, intent(in)    :: maxsubarea, gas_pcnst
     integer, intent(in)    :: jclea, jcldy
     real(wp),intent(in)    :: fcldy
     real(wp),intent(in)    :: qqcwgcm(gas_pcnst)
     real(wp),intent(inout) :: qqcwsub(gas_pcnst,maxsubarea)

     integer :: imode    ! mode index
     integer :: ispec    ! aerosol species index
     integer :: icnst    ! consitituent index 

     !----------------------------------------------------------------
     do imode = 1, ntot_amode          ! loop thru log-normal modes
      do ispec = 1, nspec_amode(imode) ! mass of individual species in a mode

           icnst = lmassptrcw_amode(ispec,imode) - loffset

           qqcwsub(icnst,jclea) = 0.0_wp
           qqcwsub(icnst,jcldy) = qqcwgcm(icnst)/fcldy

      end do ! ispec - species loop
     end do ! imode - mode loop
     !----------------------------------------------------------------

  end subroutine set_subarea_q_mass_for_cldbrn_aerosols

!--------------------------------------------------------------------------------
  subroutine get_partition_factors(  q_intrst_gcm, q_cldbrn_gcm, fcldy, fclea, &! in
                                     part_fac_q_intrst_clea, part_fac_q_intrst_cldy  )! out

      real(wp), intent(in)  ::  q_intrst_gcm  ! grid cell mean interstitial aerosol mixing ratio
      real(wp), intent(in)  ::  q_cldbrn_gcm  ! grid cell mean cloud-borne aerosol mixing ratio

      real(wp), intent(in)  ::  fcldy           ! cloudy fraction of the grid cell
      real(wp), intent(in)  ::  fclea           ! clear  fraction of the grid cell

      real(wp), intent(out) ::  part_fac_q_intrst_clea
      real(wp), intent(out) ::  part_fac_q_intrst_cldy

      real(wp) :: tmp_q_intrst_clea, tmp_q_intrst_cldy
      real(wp) :: tmp_q_cldbrn_cldy
      real(wp) :: tmp_aa

      ! Calculate mixing ratios of each subarea

      tmp_q_cldbrn_cldy = q_cldbrn_gcm/fcldy ! cloud-borne,  cloudy subarea
      tmp_q_intrst_cldy = max( 0.0_wp, ((q_intrst_gcm+q_cldbrn_gcm) - tmp_q_cldbrn_cldy) ) ! interstitial, cloudy subarea

      tmp_q_intrst_clea = (q_intrst_gcm - fcldy*tmp_q_intrst_cldy)/fclea ! interstitial, clear  subarea

      ! Calculate the corresponding paritioning factors for interstitial aerosols
      ! using the above-derived subarea mixing ratios plus the constraint that
      ! the cloud fraction weighted average of subarea mean need to match grid box mean.

      tmp_aa = max( 1.e-35_wp, tmp_q_intrst_clea*fclea ) / max( 1.e-35_wp, q_intrst_gcm )
      tmp_aa = max( 0.0_wp, min( 1.0_wp, tmp_aa ) )

      part_fac_q_intrst_clea = tmp_aa/fclea
      part_fac_q_intrst_cldy = (1.0_wp-tmp_aa)/fcldy

  end subroutine get_partition_factors


end module modal_aero_amicphys_subareas
