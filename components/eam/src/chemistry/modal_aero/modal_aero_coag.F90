!--------------------------------------------------------------------------------------
! Purpose:
!
! Modal aerosol coagulation parameterization.
!
! Revision history:
!   Richard C. Easter, 07.04.13:  Adapted from MIRAGE2 code. Some of the 
!                                 code came from the CMAQ model v4.6.
!   Hui Wan 2022: Refactored code: memoved non-4-mode codes; restructured subroutines.
!--------------------------------------------------------------------------------------
module modal_aero_coag

  use shr_kind_mod,    only:  wp => shr_kind_r8
  use cam_logfile,     only:  iulog

  use modal_aero_amicphys_control

  implicit none
  private

  public :: set_coagulation_pairs
  public :: mam_coag_1subarea

contains

subroutine set_coagulation_pairs( big_neg_int )
!----------------------------------------------------------------------
! Purpose: Set up coagulation pairs during model initialization.
!
! The coagulation pairs are
! mam version   modes involved in coagulation          # of coag pairs
! -----------   -----------------------------          ---------------
! 4 mode        accum, aitken, pcarbon/pca             3
!----------------------------------------------------------------------

   integer, intent(in) :: big_neg_int
   integer :: ip

   modefrm_coagpair(:) = big_neg_int
   modetoo_coagpair(:) = big_neg_int
   modeend_coagpair(:) = big_neg_int

   n_coagpair = 3
   ip=1; modefrm_coagpair(ip)=nait; modetoo_coagpair(ip)=nacc; modeend_coagpair(ip)=nacc
   ip=2; modefrm_coagpair(ip)=npca; modetoo_coagpair(ip)=nacc; modeend_coagpair(ip)=nacc
   ip=3; modefrm_coagpair(ip)=nait; modetoo_coagpair(ip)=npca; modeend_coagpair(ip)=nacc

end subroutine set_coagulation_pairs

subroutine mam_coag_1subarea(                                   &
     deltat,                                                    &
     temp,              pmid,             aircon,               &
     dgn_a,             dgn_awet,         wetdens,              &
     qnum_cur,                                                  &
     qaer_cur,          qaer_del_coag_out                       )
!-----------------------------------------------------------------------------------
! Purpose: 
! Considers the coagulation between aitken, pcarbon, and accum modes and
! updates aerosol mass and number mixing ratios.
! 
! This subroutine is called by MAM4's microphysics driver for clear-air conditions.
!-----------------------------------------------------------------------------------

      implicit none

      ! Arguments

      real(wp), intent(in) :: deltat                ! model timestep (s)
      real(wp), intent(in) :: temp                  ! temperature at model levels (K)
      real(wp), intent(in) :: pmid                  ! pressure at layer center (Pa)
      real(wp), intent(in) :: aircon                ! air molar concentration (kmol/m3)
      real(wp), intent(in) :: dgn_a(max_mode)
      real(wp), intent(in) :: dgn_awet(max_mode) ! dry & wet geo. mean dia. (m) of number distrib.
      real(wp), intent(in) :: wetdens(max_mode)  ! interstitial aerosol wet density (kg/m3)
                                                 ! dry & wet geo. mean dia. (m) of number distrib.

      real(wp), intent(inout), dimension(1:max_mode)           :: qnum_cur  ! current number mixing ratios (#/kmol-air)
      real(wp), intent(inout), dimension(1:max_aer,1:max_mode) :: qaer_cur  ! current mass mixing ratios (kmol/kmol-air)

      ! Mass mixing ratio changes (kmol-aer/kmol-air) to be passed to the aging parameterization
      real(wp), intent(out),   dimension(1:max_aer,1:max_agepair) :: qaer_del_coag_out  

      ! Local variables

      integer :: ip                 ! coagulation pair index
      integer :: modefrm, modetoo   ! mode indices corresponding to the origin and desination of mass transfer

      ! Coagulation rate coefficients corresponding the 0th and 3rd moments of the size distribution.

      real(wp) :: ybetaij0(max_coagpair), ybetaij3(max_coagpair), &
                  ybetaii0(max_coagpair), ybetajj0(max_coagpair)

      real(wp), dimension(1:max_mode) :: qnum_bgn  ! number mixing ratios before coagulation
      real(wp), dimension(1:max_mode) :: qnum_tavg ! average number mixing ratios during time step, calculated as
                                                     ! as the arithmetic mean of the begin- and end-values

      real(wp), dimension(1:max_aer,1:max_mode) :: qaer_bgn ! mass mixing ratios before coagulation

      !----------------------------------------------------
      ! Preparation
      !----------------------------------------------------
      ! Clip negative values

      qnum_cur = max( 0.0_wp, qnum_cur )
      qaer_cur = max( 0.0_wp, qaer_cur )

      ! Set initial values before coag

      qnum_bgn = qnum_cur
      qaer_bgn = qaer_cur

      !--------------------------------------------------------------
      ! Compute coagulation rates using the CMAQ models "fast" method
      ! (based on E. Whitby's approximation approach)
      ! Here subr. arguments are all in mks unit.
      !--------------------------------------------------------------
      do ip = 1, n_coagpair

         modefrm = modefrm_coagpair(ip)
         modetoo = modetoo_coagpair(ip)

         call getcoags_wrapper_f(                      &
            temp,                pmid,                 &! in
            dgn_awet(modefrm),   dgn_awet(modetoo),    &! in
            sigmag_aer(modefrm), sigmag_aer(modetoo),  &! in
            alnsg_aer(modefrm),  alnsg_aer(modetoo),   &! in
            wetdens(modefrm),    wetdens(modetoo),     &! in
            ybetaij0(ip),        ybetaij3(ip),         &! out
            ybetaii0(ip),        ybetajj0(ip)          )! out

      end do

      ! Convert coag coefficients from (m3/s) to (kmol-air/s)

      ybetaij0(1:n_coagpair) = ybetaij0(1:n_coagpair)*aircon
      ybetaij3(1:n_coagpair) = ybetaij3(1:n_coagpair)*aircon
      ybetaii0(1:n_coagpair) = ybetaii0(1:n_coagpair)*aircon
      ybetajj0(1:n_coagpair) = ybetajj0(1:n_coagpair)*aircon

      !---------------------------------------------------------------------------------------
      ! Advance solutions in time assuming the coag coefficients are fixed within one timestep
      !---------------------------------------------------------------------------------------
      ! First update number mixing ratios 

      call mam_coag_num_update( ybetaij0, ybetaii0, ybetajj0, deltat, qnum_bgn, &! in
                                qnum_cur, qnum_tavg                             )! inout, out 

      ! Then calculate mass transfers between modes and update mass mixing ratios

      call mam_coag_aer_update( ybetaij3, deltat, qnum_tavg, qaer_bgn, &! in
                                qaer_cur, qaer_del_coag_out            )! inout, out

end subroutine mam_coag_1subarea

subroutine mam_coag_num_update( ybetaij0, ybetaii0, ybetajj0, deltat, qnum_bgn, qnum_end, qnum_tavg )
!----------------------------------------------------------------------------------------------------
! Purpose: update aerosol number mixing ratios by taking into account self-coagulation (i.e., 
!          intra-modal coagulation) and inter-modal coagulation.
!
! Numerical treatment:
!  - Note that the updates for different modes are calculated in a sequential manner using a specific 
!    ordering because
!    - accum   number loss depends on accum number
!    - pcarbon number loss depends on pcarbon and accum number
!
!  - The average number mixing ratio over current time step
!    of other modes are used to calculate the number loss of a mode.
!----------------------------------------------------------------------------------------------------

      real(wp), intent(in) :: ybetaij0(max_coagpair) ! coag rate coefficient
      real(wp), intent(in) :: ybetaii0(max_coagpair) ! coag rate coefficient
      real(wp), intent(in) :: ybetajj0(max_coagpair) ! coag rate coefficient
      real(wp), intent(in) :: deltat

      real(wp), intent(in)    :: qnum_bgn(1:max_mode)  ! beginning values of number mixing ratios
      real(wp), intent(inout) :: qnum_end(1:max_mode)  ! end values of number mixing ratios
      real(wp), intent(out)   :: qnum_tavg(1:max_mode) ! time average defined as 0.5*(bgn+end)

      ! local variables
      real(wp) :: tmpa, rateij
      real(wp) ::       rateii

      !--------------------------------------------
      ! accum mode number loss due to self-coag
      !--------------------------------------------
      call qnum_update_selfcoag( ybetajj0(1), deltat, qnum_bgn(nacc), qnum_end(nacc) )

      qnum_tavg(nacc) = (qnum_bgn(nacc) + qnum_end(nacc))*0.5_wp

      !----------------------------------------------------------------------------
      ! pcarbon mode number loss - approximate analytical solution 
      ! using average number conc. for accum mode
      !----------------------------------------------------------------------------
      rateij = max( 0.0_wp, deltat*ybetaij0(2)*qnum_tavg(nacc) )
      rateii = max( 0.0_wp, deltat*ybetaii0(2) )

      call qnum_update_self_and_intermodal_coag( rateij, rateii, qnum_bgn(npca), qnum_end(npca) )

      qnum_tavg(npca) = (qnum_bgn(npca) + qnum_end(npca))*0.5_wp

      !-----------------------------------------------------------------------------------------
      ! aitken mode number loss - approximate analytical solution
      ! using average number conc. for accum, pcarbon, and marine-org accum modes
      !-----------------------------------------------------------------------------------------
      tmpa = ybetaij0(1)*qnum_tavg(nacc)
      tmpa = tmpa + ybetaij0(3)*qnum_tavg(npca)
      rateij = max( 0.0_wp, deltat*tmpa )
      rateii = max( 0.0_wp, deltat*ybetaii0(1) )

      call qnum_update_self_and_intermodal_coag( rateij, rateii, qnum_bgn(nait), qnum_end(nait) )

      qnum_tavg(nait) = (qnum_bgn(nait) + qnum_end(nait))*0.5_wp


end subroutine mam_coag_num_update


subroutine qnum_update_selfcoag( ybetajj0, deltat, qnum_bgn, qnum_end )
!-----------------------------------------------------------------------------------------
! Purpose: update the number mixing ratio of a single mode by considering self-coagulation
!-----------------------------------------------------------------------------------------

    real(wp),intent(in)  :: qnum_bgn   ! qnum (number mixing ratio) start value
    real(wp),intent(in)  :: ybetajj0   ! self-coagulation coefficient
    real(wp),intent(in)  :: deltat     ! timestep length
    real(wp),intent(out) :: qnum_end   ! qnum (number mixing ratio) end value

    ! Analytical solution
    qnum_end = qnum_bgn / ( 1.0_wp + ybetajj0*deltat*qnum_bgn )

end subroutine qnum_update_selfcoag

subroutine qnum_update_self_and_intermodal_coag( rateij, rateii, qnum_bgn, qnum_end )
!-----------------------------------------------------------------------------------------
! Purpose: update the number mixing ratio of a single mode considering self-coagulation 
!          and inter-modal coagulation
!-----------------------------------------------------------------------------------------

    real(wp),intent(in)  :: rateij, rateii  ! inter-modal and self-coagulation rates
    real(wp),intent(in)  :: qnum_bgn        ! qnum (number mixing ratio) start value
    real(wp),intent(out) :: qnum_end        ! qnum (number mixing ratio) end value

    real(wp) :: tmpc

    if (rateij < 1.0e-5_wp) then
       qnum_end = qnum_bgn / ( 1.0_wp + (rateij + rateii*qnum_bgn)*(1.0_wp + 0.5_wp*rateij) )
    else
       tmpc = exp(-rateij)
       qnum_end = qnum_bgn*tmpc / ( 1.0_wp + (rateii*qnum_bgn/rateij)*(1.0_wp-tmpc) )
    end if

end subroutine qnum_update_self_and_intermodal_coag


subroutine mam_coag_aer_update( ybetaij3, deltat, qnum_tavg, qaer_bgn, qaer_end, qaer_del_coag_out)
!=====================================================================================================
! Purpose: update aerosol mass mixing ratios by taking into account coagulation-induced inter-modal
!          mass transfer.
! 
! Assumed possible mass transfer pathways: 
! - coag pair 1: aitken + accumulation -> accumulation
! - coag pair 2: pca    + accumulation -> accumulation
! - coag pair 3: aitken + pca          -> pca
!
! Sorted by mass source:
!  - From aitken mode:
!      - coag pair 1: aitken + accumulation -> accumulation
!      - coag pair 3: aitken + pca          -> pca
!  - From pca mode:
!      - coag pair 2: pca    + accumulation -> accumulation
!--------------------------------------------------------------------------------
! Numerical treatment:
!
! The transfer amounts are calculated using as an exponential decay of
! the initial mass mixing ratios,
! where the decay rate is calculated using the average (over one timestep)
! number mixing ratios for each mode
!
! The mass transfer calculations are first-order accurate in time,
! because the mass transferred out of a mode does not
! include any mass transferred in during the time step.
! With this approach, the ordering is not important, but the mass transfer
! calculations are done in the reverse order of the number loss calculations
!=====================================================================================================

      ! Arguments

      real(wp), intent(in), dimension( 1:max_coagpair ) :: ybetaij3(max_coagpair)
      real(wp), intent(in) :: deltat

      real(wp), intent(in)    :: qnum_tavg(1:max_mode)
      real(wp), intent(in)    :: qaer_bgn(1:max_aer, 1:max_mode) 
      real(wp), intent(inout) :: qaer_end(1:max_aer, 1:max_mode)
      real(wp), intent(out)   :: qaer_del_coag_out(1:max_aer, 1:max_agepair)

      ! Local variables

      real(wp) :: rate1, rate2  ! coag rates = coag coeff * qnum
      real(wp) :: ratesum       ! coag rate summed over all modes
      real(wp) :: prtn1, prtn2  ! portions of mass going into different desitination modes 
      real(wp) :: tmpc
      real(wp) :: tmp_dq, tmp_xf

      integer :: iaer   ! aerosol species index

      real(wp), parameter :: epsilonx2 = epsilon( 1.0_wp )*2.0_wp

      !--------------------------------------------------------------------
      ! Initialize the array that will be passed onto aging
      !--------------------------------------------------------------------
      qaer_del_coag_out = 0._wp

      !--------------------------------------------------------------------
      ! Mass transfer out of aitken mode. Two coag pairs are involved:
      ! - coag pair 1: aitken + accumulation -> accumulation
      ! - coag pair 3: aitken + pca          -> pca
      !--------------------------------------------------------------------
      ! Calculate the rate of mass transfer into different destination modes and the sum over all modes

      rate1 = max( 0.0_wp, ybetaij3(1)*qnum_tavg(nacc) )
      rate2 = max( 0.0_wp, ybetaij3(3)*qnum_tavg(npca) )
      ratesum = rate1 + rate2

      tmpc = deltat*ratesum ! calculate coag-induced changes only when this number is not ~= zero
      if (tmpc > epsilonx2) then

         ! Portions of mass going into different modes
         prtn2 = rate2/ratesum
         prtn1 = 1.0_wp - prtn2

         tmp_xf = 1.0_wp - exp(-tmpc)    ! total fraction lost from aitken mode
         do iaer = 1, naer

            tmp_dq = tmp_xf*qaer_bgn(iaer,nait)  ! total amount lost from aitken mode

            qaer_end(iaer,nait) = qaer_end(iaer,nait) - tmp_dq          ! subtract from aitken mode
            qaer_end(iaer,nacc) = qaer_end(iaer,nacc) + tmp_dq*prtn1    ! add a portion to accumulation mode
            qaer_end(iaer,npca) = qaer_end(iaer,npca) + tmp_dq*prtn2    ! add a portion to pca mode

            ! prtn2 (pair 3) corresponds to mass transfer to pca mode, which will lead to aging.
            ! Add this amount to the total mass gained by pca mode, to be used in the aging parameterization.

              qaer_del_coag_out(iaer,i_agepair_pca) &
            = qaer_del_coag_out(iaer,i_agepair_pca) + tmp_dq*prtn2
         end do

      end if

      !--------------------------------------------------------------------
      ! Mass transfer out of pcarbon mode. Only one coag pair is involved:
      ! - coag pair 2: pca + accumulation -> accumulation
      !--------------------------------------------------------------------
      ratesum = max( 0.0_wp, ybetaij3(2)*qnum_tavg(nacc) )  ! there is only 1 destination

      tmpc = deltat*ratesum ! calculate coag-induced changes only when this number is not ~= zero
      if (tmpc > epsilonx2) then

         tmp_xf = 1.0_wp - exp(-tmpc) ! total fraction lost from pca mode
         do iaer = 1, naer
            tmp_dq = tmp_xf*qaer_bgn(iaer,npca)                   ! total amount lost from pca mode
            qaer_end(iaer,npca) = qaer_end(iaer,npca) - tmp_dq    ! subtract from pca mode
            qaer_end(iaer,nacc) = qaer_end(iaer,nacc) + tmp_dq    ! add to accumulaiton mode
         end do
      end if

end subroutine mam_coag_aer_update

subroutine getcoags_wrapper_f(&
             airtemp, airprs,   dgatk,   dgacc,    sgatk, sgacc, &
             xxlsgat, xxlsgac,  pdensat, pdensac,                &
             betaij0, betaij3,  betaii0, betajj0                 )
!-------------------------------------------------------------------------------
! Purpose: interface to subr. getcoags adapted from subr. aeroproc of CMAQ v4.6,
! with some of the parameter values from module aero_info_ae4
!-------------------------------------------------------------------------------

    use physconst, only: p0 => pstd, tmelt, boltz

    implicit none

    ! Arguments 

    real(wp), intent(in) :: airtemp  ! air temperature [ k ]
    real(wp), intent(in) :: airprs   ! air pressure in [ pa ]

    real(wp), intent(in) :: dgatk    ! aitken mode geometric mean diameter [m]
    real(wp), intent(in) :: dgacc    ! accumulation mode geometric mean diam [m]

    real(wp), intent(in) :: sgatk    ! aitken mode geometric standard deviation
    real(wp), intent(in) :: sgacc    ! accumulation mode geometric standard deviation

    real(wp), intent(in) :: xxlsgat  ! natural log of geometric standard
    real(wp), intent(in) :: xxlsgac  !  deviations

    real(wp), intent(in) :: pdensat  ! aitken mode particle density [ kg / m**3 ]
    real(wp), intent(in) :: pdensac  ! accumulation mode density [ kg / m**3 ]

    real(wp), intent(out) :: betaij0, betaii0, betajj0, betaij3

    ! Local variables and parameters 

    real(wp) :: t0  ! standard surface temperature (15 deg C) [ k ]
    real(wp), parameter :: two3 = 2.0_wp/3.0_wp

    real(wp) amu            ! atmospheric dynamic viscosity [ kg/m s ]
    real(wp) sqrt_temp      ! square root of ambient temperature
    real(wp) lamda          ! mean free path [ m ]


    ! Near-continnuum regime (independent of mode)

    real(wp) knc         ! knc = two3 * boltz *  airtemp / amu

    ! Free-molecular regime (depends upon modal density)

    real(wp) kfmat       ! kfmat = sqrt(3.0*boltz*airtemp/pdensat)
    real(wp) kfmac       ! kfmac = sqrt(3.0*boltz*airtemp/pdensac)
    real(wp) kfmatac     ! kfmatac = sqrt( 6.0 * boltz * airtemp /( pdensat + pdensac ) )

    ! Output from getcoags

    real(wp)  qn11  ! aitken mode intramodal coagulation rate [ m**3/s ] for 0th moment 
    real(wp)  qn22  ! accum. mode intramodal coagulation rate [ m**3/s ] for 0th moment 

    real(wp)  qn12  ! aitken to accumulation intermodal coagulation rate [ m**3/s ] for 0th moment
    real(wp)  qv12  ! aitken to accumulation intermodal coagulation rate [ m**3/s ] for 3rd moment

    ! For unit conversion
    real(wp)  dumatk3

    !-----------------------------------------------
    ! Prepare input to subr. getcoags 
    !-----------------------------------------------
    t0 = tmelt + 15._wp
    sqrt_temp = sqrt( airtemp)

    ! Calculate mean free path [ m ]:
    ! 6.6328e-8 is the sea level value given in table i.2.8
    ! on page 10 of u.s. standard atmosphere 1962

    lamda = 6.6328e-8_wp * p0 * airtemp  / ( t0 * airprs )

    ! Calculate dynamic viscosity [ kg m**-1 s**-1 ]:
    ! u.s. standard atmosphere 1962 page 14 expression
    ! for dynamic viscosity is:
    ! dynamic viscosity =  beta * t * sqrt(t) / ( t + s)
    ! where beta = 1.458e-6 [ kg sec^-1 k**-0.5 ], s = 110.4 [ k ].

    amu = 1.458e-6_wp * airtemp * sqrt_temp / ( airtemp + 110.4_wp )

    ! Term used in equation a6 of binkowski & shankar (1995)

    knc = two3 * boltz *  airtemp / amu

    ! Terms used in equation a5 of binkowski & shankar (1995)

    kfmat    = sqrt( 3.0_wp * boltz * airtemp / pdensat )
    kfmac    = sqrt( 3.0_wp * boltz * airtemp / pdensac )
    kfmatac  = sqrt( 6.0_wp * boltz * airtemp / ( pdensat + pdensac ) )

    !--------------------------------------------------------------------
    ! Call subr. getcoags from CMAQ to 
    ! calculate intermodal and intramodal coagulation coefficients
    ! for zeroth and intermodal coagulation coefficient for third moment
    !--------------------------------------------------------------------
    call getcoags( lamda, kfmatac, kfmat, kfmac, knc,            &! in
                   dgatk, dgacc, sgatk, sgacc, xxlsgat, xxlsgac, &! in
                   qn11, qn22, qn12, qv12                        )! out

    !--------------------------------------------------------------------
    ! Adjustments to the output from subr. getcoags
    !--------------------------------------------------------------------
    ! Clip negative values

    betaii0  = max( 0.0_wp, qn11 )
    betajj0  = max( 0.0_wp, qn22 )
    betaij0  = max( 0.0_wp, qn12 )

    ! For the mass transfer, convert from the "cmaq" coag rate parameters 
    ! to the "mirage2" parameters

    dumatk3 = ( (dgatk**3) * exp( 4.5_wp*xxlsgat*xxlsgat ) )
    betaij3  = max( 0.0_wp, qv12 / dumatk3 )

end subroutine getcoags_wrapper_f


subroutine getcoags( lamda, kfmatac, kfmat, kfmac, knc,           &
                     dgatk, dgacc, sgatk, sgacc, xxlsgat,xxlsgac, &
                     qn11, qn22, qn12, qv12 )

!//////////////////////////////////////////////////////////////////
!  subroutine getcoags calculates the coagulation rates using a new
!     approximate algorithm for the 2nd moment.  the 0th and 3rd moments
!     are done by analytic expressions from whitby et al. (1991).  the
!     correction factors are also similar to those from whitby et al.
!     (1991), but are derived from the gauss-hermite numerical
!     quadratures used by binkowski and roselle (2003).
!
!  Revision history:
!   fsb 08/25/03 coded by dr. francis s. binkowksi
!   fsb 08/25/04 added in-line documentation
!   rce 04/15/2007
!    - code taken from cmaq v4.6 code; converted to f90;
!    - added "intent" to subr arguments;
!    - renamed "r4" & "dp" variables to "rx4" & "rx8";
!    - changed "real*N" declarations to "real(rN)" (N = 4 or 8)
!   Hui Wan, 2022 : removed unused calculations
!
!  References:
!   1. whitby, e. r., p. h. mcmurry, u. shankar, and f. s. binkowski,
!   modal aerosol dynamics modeling, rep. 600/3-91/020, atmospheric
!   research and exposure assessment laboratory,
!   u.s. environmental protection agency, research triangle park, n.c.,
!   (ntis pb91-161729/as), 1991
!
!   2. binkowski, f.s. an u. shankar, the regional particulate matter
!   model 1. model decsription and preliminary results, journal of
!   geophysical research, 100, d12, pp 26,191-26,209,
!   december 20, 1995.
!
!   3. binkowski, f.s. and s.j. roselle, models-3 community
!   multiscale air quality (cmaq) model aerosol component 1:
!   model description.  j. geophys. res., vol 108, no d6, 4183
!   doi:10.1029/2001jd001409, 2003.
!//////////////////////////////////////////////////////////////////

    implicit none

    real(wp), intent(in) ::  lamda     ! mean free path [ m ]

    ! Coefficients for free molecular regime

    real(wp), intent(in) ::  kfmat     ! aitken mode
    real(wp), intent(in) ::  kfmac     ! accumulation mode
    real(wp), intent(in) ::  kfmatac   ! aitken to accumulation mode

    real(wp), intent(in) ::  knc   ! coefficient for near continnuum regime

    ! Modal geometric mean diameters: [ m ]

    real(wp), intent(in) :: dgatk          ! aitken mode
    real(wp), intent(in) :: dgacc          ! accumulation mode

    ! Modal geometric standard deviation

    real(wp), intent(in) :: sgatk          ! atken mode
    real(wp), intent(in) :: sgacc          ! accumulation mode

    ! Natural log of modal geometric standard deviation

    real(wp), intent(in) :: xxlsgat         ! aitken mode
    real(wp), intent(in) :: xxlsgac         ! accumulation mode

    ! Coagulation coefficients
    real(wp), intent(out) :: qn11, qn22, qn12, qv12 

    integer ibeta, n1, n2a, n2n ! indices for correction factors

    real(wp) :: kngat, kngac
    real(wp),parameter :: one = 1.0_wp, two = 2.0_wp, half = 0.5_wp
    real(wp),parameter :: a = 1.246_wp
    real(wp),parameter :: two3rds = 2._wp / 3._wp

    real(wp) sqrttwo  !  sqrt(two)
    real(wp) dlgsqt2  !  1/ln( sqrt( 2 ) )

    real(wp) esat01         ! aitken mode exp( log^2( sigmag )/8 )
    real(wp) esac01         ! accumulation mode exp( log^2( sigmag )/8 )

    real(wp) esat04,esat05,esat08,esat09,esat16,esat20,esat24,esat25,esat36,esat49,esat64,esat100
    real(wp) esac04,esac05,esac08,esac09,esac16,esac20,esac24,esac25,esac36,       esac64

    real(wp) dgat2, dgac2, dgat3, dgac3
    real(wp) sqdgat, sqdgac
    real(wp) sqdgat5, sqdgac5
    real(wp) sqdgat7
    real(wp) r, r2, r3, rx4
    real(wp) ri1, ri2, ri3
    real(wp) rat
    real(wp) coagfm0, coagnc0
    real(wp) coagfm3, coagnc3
    real(wp) coagfm_at, coagfm_ac
    real(wp) coagnc_at, coagnc_ac

    ! Correction factors for coagulation rates
    real(wp), save :: bm0( 10 )           ! m0 intramodal fm - rpm values
    real(wp), save :: bm0ij( 10, 10, 10 ) ! m0 intermodal fm
    real(wp), save :: bm3i( 10, 10, 10 )  ! m3 intermodal fm- rpm values

    ! Populate the arrays for the correction factors *************************************

    ! rpm 0th moment correction factors for unimodal fm coagulation  rates
    data      bm0  /   &
      0.707106785165097_wp, 0.726148960080488_wp, 0.766430744110958_wp,   &
      0.814106389441342_wp, 0.861679526483207_wp, 0.903600509090092_wp,   &
      0.936578814219156_wp, 0.960098926735545_wp, 0.975646823342881_wp,   &
      0.985397173215326_wp   /

    ! fsb new fm correction factors for m0 intermodal coagulation

    data (bm0ij (  1,  1,ibeta), ibeta = 1,10) /   &
      0.628539_wp,  0.639610_wp,  0.664514_wp,  0.696278_wp,  0.731558_wp,   &
      0.768211_wp,  0.804480_wp,  0.838830_wp,  0.870024_wp,  0.897248_wp/
    data (bm0ij (  1,  2,ibeta), ibeta = 1,10) /   &
      0.639178_wp,  0.649966_wp,  0.674432_wp,  0.705794_wp,  0.740642_wp,   &
      0.776751_wp,  0.812323_wp,  0.845827_wp,  0.876076_wp,  0.902324_wp/
    data (bm0ij (  1,  3,ibeta), ibeta = 1,10) /   &
      0.663109_wp,  0.673464_wp,  0.697147_wp,  0.727637_wp,  0.761425_wp,   &
      0.796155_wp,  0.829978_wp,  0.861419_wp,  0.889424_wp,  0.913417_wp/
    data (bm0ij (  1,  4,ibeta), ibeta = 1,10) /   &
      0.693693_wp,  0.703654_wp,  0.726478_wp,  0.755786_wp,  0.787980_wp,   &
      0.820626_wp,  0.851898_wp,  0.880459_wp,  0.905465_wp,  0.926552_wp/
    data (bm0ij (  1,  5,ibeta), ibeta = 1,10) /   &
      0.727803_wp,  0.737349_wp,  0.759140_wp,  0.786870_wp,  0.816901_wp,   &
      0.846813_wp,  0.874906_wp,  0.900060_wp,  0.921679_wp,  0.939614_wp/
    data (bm0ij (  1,  6,ibeta), ibeta = 1,10) /   &
      0.763461_wp,  0.772483_wp,  0.792930_wp,  0.818599_wp,  0.845905_wp,   &
      0.872550_wp,  0.897051_wp,  0.918552_wp,  0.936701_wp,  0.951528_wp/
    data (bm0ij (  1,  7,ibeta), ibeta = 1,10) /   &
      0.799021_wp,  0.807365_wp,  0.826094_wp,  0.849230_wp,  0.873358_wp,   &
      0.896406_wp,  0.917161_wp,  0.935031_wp,  0.949868_wp,  0.961828_wp/
    data (bm0ij (  1,  8,ibeta), ibeta = 1,10) /   &
      0.833004_wp,  0.840514_wp,  0.857192_wp,  0.877446_wp,  0.898147_wp,   &
      0.917518_wp,  0.934627_wp,  0.949106_wp,  0.960958_wp,  0.970403_wp/
    data (bm0ij (  1,  9,ibeta), ibeta = 1,10) /   &
      0.864172_wp,  0.870734_wp,  0.885153_wp,  0.902373_wp,  0.919640_wp,   &
      0.935494_wp,  0.949257_wp,  0.960733_wp,  0.970016_wp,  0.977346_wp/
    data (bm0ij (  1, 10,ibeta), ibeta = 1,10) /   &
      0.891658_wp,  0.897227_wp,  0.909343_wp,  0.923588_wp,  0.937629_wp,   &
      0.950307_wp,  0.961151_wp,  0.970082_wp,  0.977236_wp,  0.982844_wp/
    data (bm0ij (  2,  1,ibeta), ibeta = 1,10) /   &
      0.658724_wp,  0.670587_wp,  0.697539_wp,  0.731890_wp,  0.769467_wp,   &
      0.807391_wp,  0.843410_wp,  0.875847_wp,  0.903700_wp,  0.926645_wp/
    data (bm0ij (  2,  2,ibeta), ibeta = 1,10) /   &
      0.667070_wp,  0.678820_wp,  0.705538_wp,  0.739591_wp,  0.776758_wp,   &
      0.814118_wp,  0.849415_wp,  0.881020_wp,  0.908006_wp,  0.930121_wp/
    data (bm0ij (  2,  3,ibeta), ibeta = 1,10) /   &
      0.686356_wp,  0.697839_wp,  0.723997_wp,  0.757285_wp,  0.793389_wp,   &
      0.829313_wp,  0.862835_wp,  0.892459_wp,  0.917432_wp,  0.937663_wp/
    data (bm0ij (  2,  4,ibeta), ibeta = 1,10) /   &
      0.711425_wp,  0.722572_wp,  0.747941_wp,  0.780055_wp,  0.814518_wp,   &
      0.848315_wp,  0.879335_wp,  0.906290_wp,  0.928658_wp,  0.946526_wp/
    data (bm0ij (  2,  5,ibeta), ibeta = 1,10) /   &
      0.739575_wp,  0.750307_wp,  0.774633_wp,  0.805138_wp,  0.837408_wp,   &
      0.868504_wp,  0.896517_wp,  0.920421_wp,  0.939932_wp,  0.955299_wp/
    data (bm0ij (  2,  6,ibeta), ibeta = 1,10) /   &
      0.769143_wp,  0.779346_wp,  0.802314_wp,  0.830752_wp,  0.860333_wp,   &
      0.888300_wp,  0.913014_wp,  0.933727_wp,  0.950370_wp,  0.963306_wp/
    data (bm0ij (  2,  7,ibeta), ibeta = 1,10) /   &
      0.798900_wp,  0.808431_wp,  0.829700_wp,  0.855653_wp,  0.882163_wp,   &
      0.906749_wp,  0.928075_wp,  0.945654_wp,  0.959579_wp,  0.970280_wp/
    data (bm0ij (  2,  8,ibeta), ibeta = 1,10) /   &
      0.827826_wp,  0.836542_wp,  0.855808_wp,  0.878954_wp,  0.902174_wp,   &
      0.923316_wp,  0.941345_wp,  0.955989_wp,  0.967450_wp,  0.976174_wp/
    data (bm0ij (  2,  9,ibeta), ibeta = 1,10) /   &
      0.855068_wp,  0.862856_wp,  0.879900_wp,  0.900068_wp,  0.919956_wp,   &
      0.937764_wp,  0.952725_wp,  0.964726_wp,  0.974027_wp,  0.981053_wp/
    data (bm0ij (  2, 10,ibeta), ibeta = 1,10) /   &
      0.879961_wp,  0.886755_wp,  0.901484_wp,  0.918665_wp,  0.935346_wp,   &
      0.950065_wp,  0.962277_wp,  0.971974_wp,  0.979432_wp,  0.985033_wp/
    data (bm0ij (  3,  1,ibeta), ibeta = 1,10) /   &
      0.724166_wp,  0.735474_wp,  0.761359_wp,  0.794045_wp,  0.828702_wp,   &
      0.862061_wp,  0.891995_wp,  0.917385_wp,  0.937959_wp,  0.954036_wp/
    data (bm0ij (  3,  2,ibeta), ibeta = 1,10) /   &
      0.730416_wp,  0.741780_wp,  0.767647_wp,  0.800116_wp,  0.834344_wp,   &
      0.867093_wp,  0.896302_wp,  0.920934_wp,  0.940790_wp,  0.956237_wp/
    data (bm0ij (  3,  3,ibeta), ibeta = 1,10) /   &
      0.745327_wp,  0.756664_wp,  0.782255_wp,  0.814026_wp,  0.847107_wp,   &
      0.878339_wp,  0.905820_wp,  0.928699_wp,  0.946931_wp,  0.960977_wp/
    data (bm0ij (  3,  4,ibeta), ibeta = 1,10) /   &
      0.765195_wp,  0.776312_wp,  0.801216_wp,  0.831758_wp,  0.863079_wp,   &
      0.892159_wp,  0.917319_wp,  0.937939_wp,  0.954145_wp,  0.966486_wp/
    data (bm0ij (  3,  5,ibeta), ibeta = 1,10) /   &
      0.787632_wp,  0.798347_wp,  0.822165_wp,  0.850985_wp,  0.880049_wp,   &
      0.906544_wp,  0.929062_wp,  0.947218_wp,  0.961288_wp,  0.971878_wp/
    data (bm0ij (  3,  6,ibeta), ibeta = 1,10) /   &
      0.811024_wp,  0.821179_wp,  0.843557_wp,  0.870247_wp,  0.896694_wp,   &
      0.920365_wp,  0.940131_wp,  0.955821_wp,  0.967820_wp,  0.976753_wp/
    data (bm0ij (  3,  7,ibeta), ibeta = 1,10) /   &
      0.834254_wp,  0.843709_wp,  0.864356_wp,  0.888619_wp,  0.912245_wp,   &
      0.933019_wp,  0.950084_wp,  0.963438_wp,  0.973530_wp,  0.980973_wp/
    data (bm0ij (  3,  8,ibeta), ibeta = 1,10) /   &
      0.856531_wp,  0.865176_wp,  0.883881_wp,  0.905544_wp,  0.926290_wp,   &
      0.944236_wp,  0.958762_wp,  0.969988_wp,  0.978386_wp,  0.984530_wp/
    data (bm0ij (  3,  9,ibeta), ibeta = 1,10) /   &
      0.877307_wp,  0.885070_wp,  0.901716_wp,  0.920729_wp,  0.938663_wp,   &
      0.953951_wp,  0.966169_wp,  0.975512_wp,  0.982442_wp,  0.987477_wp/
    data (bm0ij (  3, 10,ibeta), ibeta = 1,10) /   &
      0.896234_wp,  0.903082_wp,  0.917645_wp,  0.934069_wp,  0.949354_wp,   &
      0.962222_wp,  0.972396_wp,  0.980107_wp,  0.985788_wp,  0.989894_wp/
    data (bm0ij (  4,  1,ibeta), ibeta = 1,10) /   &
      0.799294_wp,  0.809144_wp,  0.831293_wp,  0.858395_wp,  0.885897_wp,   &
      0.911031_wp,  0.932406_wp,  0.949642_wp,  0.963001_wp,  0.973062_wp/
    data (bm0ij (  4,  2,ibeta), ibeta = 1,10) /   &
      0.804239_wp,  0.814102_wp,  0.836169_wp,  0.862984_wp,  0.890003_wp,   &
      0.914535_wp,  0.935274_wp,  0.951910_wp,  0.964748_wp,  0.974381_wp/
    data (bm0ij (  4,  3,ibeta), ibeta = 1,10) /   &
      0.815910_wp,  0.825708_wp,  0.847403_wp,  0.873389_wp,  0.899185_wp,   &
      0.922275_wp,  0.941543_wp,  0.956826_wp,  0.968507_wp,  0.977204_wp/
    data (bm0ij (  4,  4,ibeta), ibeta = 1,10) /   &
      0.831348_wp,  0.840892_wp,  0.861793_wp,  0.886428_wp,  0.910463_wp,   &
      0.931614_wp,  0.948993_wp,  0.962593_wp,  0.972872_wp,  0.980456_wp/
    data (bm0ij (  4,  5,ibeta), ibeta = 1,10) /   &
      0.848597_wp,  0.857693_wp,  0.877402_wp,  0.900265_wp,  0.922180_wp,   &
      0.941134_wp,  0.956464_wp,  0.968298_wp,  0.977143_wp,  0.983611_wp/
    data (bm0ij (  4,  6,ibeta), ibeta = 1,10) /   &
      0.866271_wp,  0.874764_wp,  0.892984_wp,  0.913796_wp,  0.933407_wp,   &
      0.950088_wp,  0.963380_wp,  0.973512_wp,  0.981006_wp,  0.986440_wp/
    data (bm0ij (  4,  7,ibeta), ibeta = 1,10) /   &
      0.883430_wp,  0.891216_wp,  0.907762_wp,  0.926388_wp,  0.943660_wp,   &
      0.958127_wp,  0.969499_wp,  0.978070_wp,  0.984351_wp,  0.988872_wp/
    data (bm0ij (  4,  8,ibeta), ibeta = 1,10) /   &
      0.899483_wp,  0.906505_wp,  0.921294_wp,  0.937719_wp,  0.952729_wp,   &
      0.965131_wp,  0.974762_wp,  0.981950_wp,  0.987175_wp,  0.990912_wp/
    data (bm0ij (  4,  9,ibeta), ibeta = 1,10) /   &
      0.914096_wp,  0.920337_wp,  0.933373_wp,  0.947677_wp,  0.960579_wp,   &
      0.971111_wp,  0.979206_wp,  0.985196_wp,  0.989520_wp,  0.992597_wp/
    data (bm0ij (  4, 10,ibeta), ibeta = 1,10) /   &
      0.927122_wp,  0.932597_wp,  0.943952_wp,  0.956277_wp,  0.967268_wp,   &
      0.976147_wp,  0.982912_wp,  0.987882_wp,  0.991450_wp,  0.993976_wp/
    data (bm0ij (  5,  1,ibeta), ibeta = 1,10) /   &
      0.865049_wp,  0.872851_wp,  0.889900_wp,  0.909907_wp,  0.929290_wp,   &
      0.946205_wp,  0.959991_wp,  0.970706_wp,  0.978764_wp,  0.984692_wp/
    data (bm0ij (  5,  2,ibeta), ibeta = 1,10) /   &
      0.868989_wp,  0.876713_wp,  0.893538_wp,  0.913173_wp,  0.932080_wp,   &
      0.948484_wp,  0.961785_wp,  0.972080_wp,  0.979796_wp,  0.985457_wp/
    data (bm0ij (  5,  3,ibeta), ibeta = 1,10) /   &
      0.878010_wp,  0.885524_wp,  0.901756_wp,  0.920464_wp,  0.938235_wp,   &
      0.953461_wp,  0.965672_wp,  0.975037_wp,  0.982005_wp,  0.987085_wp/
    data (bm0ij (  5,  4,ibeta), ibeta = 1,10) /   &
      0.889534_wp,  0.896698_wp,  0.912012_wp,  0.929395_wp,  0.945647_wp,   &
      0.959366_wp,  0.970227_wp,  0.978469_wp,  0.984547_wp,  0.988950_wp/
    data (bm0ij (  5,  5,ibeta), ibeta = 1,10) /   &
      0.902033_wp,  0.908713_wp,  0.922848_wp,  0.938648_wp,  0.953186_wp,   &
      0.965278_wp,  0.974729_wp,  0.981824_wp,  0.987013_wp,  0.990746_wp/
    data (bm0ij (  5,  6,ibeta), ibeta = 1,10) /   &
      0.914496_wp,  0.920599_wp,  0.933389_wp,  0.947485_wp,  0.960262_wp,   &
      0.970743_wp,  0.978839_wp,  0.984858_wp,  0.989225_wp,  0.992348_wp/
    data (bm0ij (  5,  7,ibeta), ibeta = 1,10) /   &
      0.926281_wp,  0.931761_wp,  0.943142_wp,  0.955526_wp,  0.966600_wp,   &
      0.975573_wp,  0.982431_wp,  0.987485_wp,  0.991128_wp,  0.993718_wp/
    data (bm0ij (  5,  8,ibeta), ibeta = 1,10) /   &
      0.937029_wp,  0.941877_wp,  0.951868_wp,  0.962615_wp,  0.972112_wp,   &
      0.979723_wp,  0.985488_wp,  0.989705_wp,  0.992725_wp,  0.994863_wp/
    data (bm0ij (  5,  9,ibeta), ibeta = 1,10) /   &
      0.946580_wp,  0.950819_wp,  0.959494_wp,  0.968732_wp,  0.976811_wp,   &
      0.983226_wp,  0.988047_wp,  0.991550_wp,  0.994047_wp,  0.995806_wp/
    data (bm0ij (  5, 10,ibeta), ibeta = 1,10) /   &
      0.954909_wp,  0.958581_wp,  0.966049_wp,  0.973933_wp,  0.980766_wp,   &
      0.986149_wp,  0.990166_wp,  0.993070_wp,  0.995130_wp,  0.996577_wp/
    data (bm0ij (  6,  1,ibeta), ibeta = 1,10) /   &
      0.914182_wp,  0.919824_wp,  0.931832_wp,  0.945387_wp,  0.957999_wp,   &
      0.968606_wp,  0.976982_wp,  0.983331_wp,  0.988013_wp,  0.991407_wp/
    data (bm0ij (  6,  2,ibeta), ibeta = 1,10) /   &
      0.917139_wp,  0.922665_wp,  0.934395_wp,  0.947580_wp,  0.959792_wp,   &
      0.970017_wp,  0.978062_wp,  0.984138_wp,  0.988609_wp,  0.991843_wp/
    data (bm0ij (  6,  3,ibeta), ibeta = 1,10) /   &
      0.923742_wp,  0.928990_wp,  0.940064_wp,  0.952396_wp,  0.963699_wp,   &
      0.973070_wp,  0.980381_wp,  0.985866_wp,  0.989878_wp,  0.992768_wp/
    data (bm0ij (  6,  4,ibeta), ibeta = 1,10) /   &
      0.931870_wp,  0.936743_wp,  0.946941_wp,  0.958162_wp,  0.968318_wp,   &
      0.976640_wp,  0.983069_wp,  0.987853_wp,  0.991330_wp,  0.993822_wp/
    data (bm0ij (  6,  5,ibeta), ibeta = 1,10) /   &
      0.940376_wp,  0.944807_wp,  0.954004_wp,  0.963999_wp,  0.972928_wp,   &
      0.980162_wp,  0.985695_wp,  0.989779_wp,  0.992729_wp,  0.994833_wp/
    data (bm0ij (  6,  6,ibeta), ibeta = 1,10) /   &
      0.948597_wp,  0.952555_wp,  0.960703_wp,  0.969454_wp,  0.977181_wp,   &
      0.983373_wp,  0.988067_wp,  0.991507_wp,  0.993977_wp,  0.995730_wp/
    data (bm0ij (  6,  7,ibeta), ibeta = 1,10) /   &
      0.956167_wp,  0.959648_wp,  0.966763_wp,  0.974326_wp,  0.980933_wp,   &
      0.986177_wp,  0.990121_wp,  0.992993_wp,  0.995045_wp,  0.996495_wp/
    data (bm0ij (  6,  8,ibeta), ibeta = 1,10) /   &
      0.962913_wp,  0.965937_wp,  0.972080_wp,  0.978552_wp,  0.984153_wp,   &
      0.988563_wp,  0.991857_wp,  0.994242_wp,  0.995938_wp,  0.997133_wp/
    data (bm0ij (  6,  9,ibeta), ibeta = 1,10) /   &
      0.968787_wp,  0.971391_wp,  0.976651_wp,  0.982148_wp,  0.986869_wp,   &
      0.990560_wp,  0.993301_wp,  0.995275_wp,  0.996675_wp,  0.997657_wp/
    data (bm0ij (  6, 10,ibeta), ibeta = 1,10) /   &
      0.973822_wp,  0.976047_wp,  0.980523_wp,  0.985170_wp,  0.989134_wp,   &
      0.992215_wp,  0.994491_wp,  0.996124_wp,  0.997277_wp,  0.998085_wp/
    data (bm0ij (  7,  1,ibeta), ibeta = 1,10) /   &
      0.947410_wp,  0.951207_wp,  0.959119_wp,  0.967781_wp,  0.975592_wp,   &
      0.981981_wp,  0.986915_wp,  0.990590_wp,  0.993266_wp,  0.995187_wp/
    data (bm0ij (  7,  2,ibeta), ibeta = 1,10) /   &
      0.949477_wp,  0.953161_wp,  0.960824_wp,  0.969187_wp,  0.976702_wp,   &
      0.982831_wp,  0.987550_wp,  0.991057_wp,  0.993606_wp,  0.995434_wp/
    data (bm0ij (  7,  3,ibeta), ibeta = 1,10) /   &
      0.954008_wp,  0.957438_wp,  0.964537_wp,  0.972232_wp,  0.979095_wp,   &
      0.984653_wp,  0.988907_wp,  0.992053_wp,  0.994330_wp,  0.995958_wp/
    data (bm0ij (  7,  4,ibeta), ibeta = 1,10) /   &
      0.959431_wp,  0.962539_wp,  0.968935_wp,  0.975808_wp,  0.981882_wp,   &
      0.986759_wp,  0.990466_wp,  0.993190_wp,  0.995153_wp,  0.996552_wp/
    data (bm0ij (  7,  5,ibeta), ibeta = 1,10) /   &
      0.964932_wp,  0.967693_wp,  0.973342_wp,  0.979355_wp,  0.984620_wp,   &
      0.988812_wp,  0.991974_wp,  0.994285_wp,  0.995943_wp,  0.997119_wp/
    data (bm0ij (  7,  6,ibeta), ibeta = 1,10) /   &
      0.970101_wp,  0.972517_wp,  0.977428_wp,  0.982612_wp,  0.987110_wp,   &
      0.990663_wp,  0.993326_wp,  0.995261_wp,  0.996644_wp,  0.997621_wp/
    data (bm0ij (  7,  7,ibeta), ibeta = 1,10) /   &
      0.974746_wp,  0.976834_wp,  0.981055_wp,  0.985475_wp,  0.989280_wp,   &
      0.992265_wp,  0.994488_wp,  0.996097_wp,  0.997241_wp,  0.998048_wp/
    data (bm0ij (  7,  8,ibeta), ibeta = 1,10) /   &
      0.978804_wp,  0.980591_wp,  0.984187_wp,  0.987927_wp,  0.991124_wp,   &
      0.993617_wp,  0.995464_wp,  0.996795_wp,  0.997739_wp,  0.998403_wp/
    data (bm0ij (  7,  9,ibeta), ibeta = 1,10) /   &
      0.982280_wp,  0.983799_wp,  0.986844_wp,  0.989991_wp,  0.992667_wp,   &
      0.994742_wp,  0.996273_wp,  0.997372_wp,  0.998149_wp,  0.998695_wp/
    data (bm0ij (  7, 10,ibeta), ibeta = 1,10) /   &
      0.985218_wp,  0.986503_wp,  0.989071_wp,  0.991711_wp,  0.993945_wp,   &
      0.995669_wp,  0.996937_wp,  0.997844_wp,  0.998484_wp,  0.998932_wp/
    data (bm0ij (  8,  1,ibeta), ibeta = 1,10) /   &
      0.968507_wp,  0.970935_wp,  0.975916_wp,  0.981248_wp,  0.985947_wp,   &
      0.989716_wp,  0.992580_wp,  0.994689_wp,  0.996210_wp,  0.997297_wp/
    data (bm0ij (  8,  2,ibeta), ibeta = 1,10) /   &
      0.969870_wp,  0.972210_wp,  0.977002_wp,  0.982119_wp,  0.986619_wp,   &
      0.990219_wp,  0.992951_wp,  0.994958_wp,  0.996405_wp,  0.997437_wp/
    data (bm0ij (  8,  3,ibeta), ibeta = 1,10) /   &
      0.972820_wp,  0.974963_wp,  0.979339_wp,  0.983988_wp,  0.988054_wp,   &
      0.991292_wp,  0.993738_wp,  0.995529_wp,  0.996817_wp,  0.997734_wp/
    data (bm0ij (  8,  4,ibeta), ibeta = 1,10) /   &
      0.976280_wp,  0.978186_wp,  0.982060_wp,  0.986151_wp,  0.989706_wp,   &
      0.992520_wp,  0.994636_wp,  0.996179_wp,  0.997284_wp,  0.998069_wp/
    data (bm0ij (  8,  5,ibeta), ibeta = 1,10) /   &
      0.979711_wp,  0.981372_wp,  0.984735_wp,  0.988263_wp,  0.991309_wp,   &
      0.993706_wp,  0.995499_wp,  0.996801_wp,  0.997730_wp,  0.998389_wp/
    data (bm0ij (  8,  6,ibeta), ibeta = 1,10) /   &
      0.982863_wp,  0.984292_wp,  0.987172_wp,  0.990174_wp,  0.992750_wp,   &
      0.994766_wp,  0.996266_wp,  0.997352_wp,  0.998125_wp,  0.998670_wp/
    data (bm0ij (  8,  7,ibeta), ibeta = 1,10) /   &
      0.985642_wp,  0.986858_wp,  0.989301_wp,  0.991834_wp,  0.993994_wp,   &
      0.995676_wp,  0.996923_wp,  0.997822_wp,  0.998460_wp,  0.998910_wp/
    data (bm0ij (  8,  8,ibeta), ibeta = 1,10) /   &
      0.988029_wp,  0.989058_wp,  0.991116_wp,  0.993240_wp,  0.995043_wp,   &
      0.996440_wp,  0.997472_wp,  0.998214_wp,  0.998739_wp,  0.999108_wp/
    data (bm0ij (  8,  9,ibeta), ibeta = 1,10) /   &
      0.990046_wp,  0.990912_wp,  0.992640_wp,  0.994415_wp,  0.995914_wp,   &
      0.997073_wp,  0.997925_wp,  0.998536_wp,  0.998968_wp,  0.999271_wp/
    data (bm0ij (  8, 10,ibeta), ibeta = 1,10) /   &
      0.991732_wp,  0.992459_wp,  0.993906_wp,  0.995386_wp,  0.996633_wp,   &
      0.997592_wp,  0.998296_wp,  0.998799_wp,  0.999154_wp,  0.999403_wp/
    data (bm0ij (  9,  1,ibeta), ibeta = 1,10) /   &
      0.981392_wp,  0.982893_wp,  0.985938_wp,  0.989146_wp,  0.991928_wp,   &
      0.994129_wp,  0.995783_wp,  0.996991_wp,  0.997857_wp,  0.998473_wp/
    data (bm0ij (  9,  2,ibeta), ibeta = 1,10) /   &
      0.982254_wp,  0.983693_wp,  0.986608_wp,  0.989673_wp,  0.992328_wp,   &
      0.994424_wp,  0.995998_wp,  0.997146_wp,  0.997969_wp,  0.998553_wp/
    data (bm0ij (  9,  3,ibeta), ibeta = 1,10) /   &
      0.984104_wp,  0.985407_wp,  0.988040_wp,  0.990798_wp,  0.993178_wp,   &
      0.995052_wp,  0.996454_wp,  0.997474_wp,  0.998204_wp,  0.998722_wp/
    data (bm0ij (  9,  4,ibeta), ibeta = 1,10) /   &
      0.986243_wp,  0.987386_wp,  0.989687_wp,  0.992087_wp,  0.994149_wp,   &
      0.995765_wp,  0.996971_wp,  0.997846_wp,  0.998470_wp,  0.998913_wp/
    data (bm0ij (  9,  5,ibeta), ibeta = 1,10) /   &
      0.988332_wp,  0.989313_wp,  0.991284_wp,  0.993332_wp,  0.995082_wp,   &
      0.996449_wp,  0.997465_wp,  0.998200_wp,  0.998723_wp,  0.999093_wp/
    data (bm0ij (  9,  6,ibeta), ibeta = 1,10) /   &
      0.990220_wp,  0.991053_wp,  0.992721_wp,  0.994445_wp,  0.995914_wp,   &
      0.997056_wp,  0.997902_wp,  0.998513_wp,  0.998947_wp,  0.999253_wp/
    data (bm0ij (  9,  7,ibeta), ibeta = 1,10) /   &
      0.991859_wp,  0.992561_wp,  0.993961_wp,  0.995403_wp,  0.996626_wp,   &
      0.997574_wp,  0.998274_wp,  0.998778_wp,  0.999136_wp,  0.999387_wp/
    data (bm0ij (  9,  8,ibeta), ibeta = 1,10) /   &
      0.993250_wp,  0.993837_wp,  0.995007_wp,  0.996208_wp,  0.997223_wp,   &
      0.998007_wp,  0.998584_wp,  0.998999_wp,  0.999293_wp,  0.999499_wp/
    data (bm0ij (  9,  9,ibeta), ibeta = 1,10) /   &
      0.994413_wp,  0.994903_wp,  0.995878_wp,  0.996876_wp,  0.997716_wp,   &
      0.998363_wp,  0.998839_wp,  0.999180_wp,  0.999421_wp,  0.999591_wp/
    data (bm0ij (  9, 10,ibeta), ibeta = 1,10) /   &
      0.995376_wp,  0.995785_wp,  0.996597_wp,  0.997425_wp,  0.998121_wp,   &
      0.998655_wp,  0.999048_wp,  0.999328_wp,  0.999526_wp,  0.999665_wp/
    data (bm0ij ( 10,  1,ibeta), ibeta = 1,10) /   &
      0.989082_wp,  0.989991_wp,  0.991819_wp,  0.993723_wp,  0.995357_wp,   &
      0.996637_wp,  0.997592_wp,  0.998286_wp,  0.998781_wp,  0.999132_wp/
    data (bm0ij ( 10,  2,ibeta), ibeta = 1,10) /   &
      0.989613_wp,  0.990480_wp,  0.992224_wp,  0.994039_wp,  0.995594_wp,   &
      0.996810_wp,  0.997717_wp,  0.998375_wp,  0.998845_wp,  0.999178_wp/
    data (bm0ij ( 10,  3,ibeta), ibeta = 1,10) /   &
      0.990744_wp,  0.991523_wp,  0.993086_wp,  0.994708_wp,  0.996094_wp,   &
      0.997176_wp,  0.997981_wp,  0.998564_wp,  0.998980_wp,  0.999274_wp/
    data (bm0ij ( 10,  4,ibeta), ibeta = 1,10) /   &
      0.992041_wp,  0.992716_wp,  0.994070_wp,  0.995470_wp,  0.996662_wp,   &
      0.997591_wp,  0.998280_wp,  0.998778_wp,  0.999133_wp,  0.999383_wp/
    data (bm0ij ( 10,  5,ibeta), ibeta = 1,10) /   &
      0.993292_wp,  0.993867_wp,  0.995015_wp,  0.996199_wp,  0.997205_wp,   &
      0.997985_wp,  0.998564_wp,  0.998981_wp,  0.999277_wp,  0.999487_wp/
    data (bm0ij ( 10,  6,ibeta), ibeta = 1,10) /   &
      0.994411_wp,  0.994894_wp,  0.995857_wp,  0.996847_wp,  0.997685_wp,   &
      0.998334_wp,  0.998814_wp,  0.999159_wp,  0.999404_wp,  0.999577_wp/
    data (bm0ij ( 10,  7,ibeta), ibeta = 1,10) /   &
      0.995373_wp,  0.995776_wp,  0.996577_wp,  0.997400_wp,  0.998094_wp,   &
      0.998630_wp,  0.999026_wp,  0.999310_wp,  0.999512_wp,  0.999654_wp/
    data (bm0ij ( 10,  8,ibeta), ibeta = 1,10) /   &
      0.996181_wp,  0.996516_wp,  0.997181_wp,  0.997861_wp,  0.998435_wp,   &
      0.998877_wp,  0.999202_wp,  0.999435_wp,  0.999601_wp,  0.999717_wp/
    data (bm0ij ( 10,  9,ibeta), ibeta = 1,10) /   &
      0.996851_wp,  0.997128_wp,  0.997680_wp,  0.998242_wp,  0.998715_wp,   &
      0.999079_wp,  0.999346_wp,  0.999538_wp,  0.999673_wp,  0.999769_wp/
    data (bm0ij ( 10, 10,ibeta), ibeta = 1,10) /   &
      0.997402_wp,  0.997632_wp,  0.998089_wp,  0.998554_wp,  0.998945_wp,   &
      0.999244_wp,  0.999464_wp,  0.999622_wp,  0.999733_wp,  0.999811_wp/


    ! rpm....   3rd moment nuclei mode corr. fac. for bimodal fm coag rate

    data (bm3i( 1, 1,ibeta ), ibeta=1,10)/   &
      0.70708_wp,0.71681_wp,0.73821_wp,0.76477_wp,0.79350_wp,0.82265_wp,0.85090_wp,0.87717_wp,   &
      0.90069_wp,0.92097_wp/
    data (bm3i( 1, 2,ibeta ), ibeta=1,10)/   &
      0.72172_wp,0.73022_wp,0.74927_wp,0.77324_wp,0.79936_wp,0.82601_wp,0.85199_wp,0.87637_wp,   &
      0.89843_wp,0.91774_wp/
    data (bm3i( 1, 3,ibeta ), ibeta=1,10)/   &
      0.78291_wp,0.78896_wp,0.80286_wp,0.82070_wp,0.84022_wp,0.85997_wp,0.87901_wp,0.89669_wp,   &
      0.91258_wp,0.92647_wp/
    data (bm3i( 1, 4,ibeta ), ibeta=1,10)/   &
      0.87760_wp,0.88147_wp,0.89025_wp,0.90127_wp,0.91291_wp,0.92420_wp,0.93452_wp,0.94355_wp,   &
      0.95113_wp,0.95726_wp/
    data (bm3i( 1, 5,ibeta ), ibeta=1,10)/   &
      0.94988_wp,0.95184_wp,0.95612_wp,0.96122_wp,0.96628_wp,0.97085_wp,0.97467_wp,0.97763_wp,   &
      0.97971_wp,0.98089_wp/
    data (bm3i( 1, 6,ibeta ), ibeta=1,10)/   &
      0.98318_wp,0.98393_wp,0.98551_wp,0.98728_wp,0.98889_wp,0.99014_wp,0.99095_wp,0.99124_wp,   &
      0.99100_wp,0.99020_wp/
    data (bm3i( 1, 7,ibeta ), ibeta=1,10)/   &
      0.99480_wp,0.99504_wp,0.99551_wp,0.99598_wp,0.99629_wp,0.99635_wp,0.99611_wp,0.99550_wp,   &
      0.99450_wp,0.99306_wp/
    data (bm3i( 1, 8,ibeta ), ibeta=1,10)/   &
      0.99842_wp,0.99848_wp,0.99858_wp,0.99861_wp,0.99850_wp,0.99819_wp,0.99762_wp,0.99674_wp,   &
      0.99550_wp,0.99388_wp/
    data (bm3i( 1, 9,ibeta ), ibeta=1,10)/   &
      0.99951_wp,0.99951_wp,0.99949_wp,0.99939_wp,0.99915_wp,0.99872_wp,0.99805_wp,0.99709_wp,   &
      0.99579_wp,0.99411_wp/
    data (bm3i( 1,10,ibeta ), ibeta=1,10)/   &
      0.99984_wp,0.99982_wp,0.99976_wp,0.99962_wp,0.99934_wp,0.99888_wp,0.99818_wp,0.99719_wp,   &
      0.99587_wp,0.99417_wp/
    data (bm3i( 2, 1,ibeta ), ibeta=1,10)/   &
      0.72957_wp,0.73993_wp,0.76303_wp,0.79178_wp,0.82245_wp,0.85270_wp,0.88085_wp,0.90578_wp,   &
      0.92691_wp,0.94415_wp/
    data (bm3i( 2, 2,ibeta ), ibeta=1,10)/   &
      0.72319_wp,0.73320_wp,0.75547_wp,0.78323_wp,0.81307_wp,0.84287_wp,0.87107_wp,0.89651_wp,   &
      0.91852_wp,0.93683_wp/
    data (bm3i( 2, 3,ibeta ), ibeta=1,10)/   &
      0.74413_wp,0.75205_wp,0.76998_wp,0.79269_wp,0.81746_wp,0.84258_wp,0.86685_wp,0.88938_wp,   &
      0.90953_wp,0.92695_wp/
    data (bm3i( 2, 4,ibeta ), ibeta=1,10)/   &
      0.82588_wp,0.83113_wp,0.84309_wp,0.85825_wp,0.87456_wp,0.89072_wp,0.90594_wp,0.91972_wp,   &
      0.93178_wp,0.94203_wp/
    data (bm3i( 2, 5,ibeta ), ibeta=1,10)/   &
      0.91886_wp,0.92179_wp,0.92831_wp,0.93624_wp,0.94434_wp,0.95192_wp,0.95856_wp,0.96409_wp,   &
      0.96845_wp,0.97164_wp/
    data (bm3i( 2, 6,ibeta ), ibeta=1,10)/   &
      0.97129_wp,0.97252_wp,0.97515_wp,0.97818_wp,0.98108_wp,0.98354_wp,0.98542_wp,0.98665_wp,   &
      0.98721_wp,0.98709_wp/
    data (bm3i( 2, 7,ibeta ), ibeta=1,10)/   &
      0.99104_wp,0.99145_wp,0.99230_wp,0.99320_wp,0.99394_wp,0.99439_wp,0.99448_wp,0.99416_wp,   &
      0.99340_wp,0.99217_wp/
    data (bm3i( 2, 8,ibeta ), ibeta=1,10)/   &
      0.99730_wp,0.99741_wp,0.99763_wp,0.99779_wp,0.99782_wp,0.99762_wp,0.99715_wp,0.99636_wp,   &
      0.99519_wp,0.99363_wp/
    data (bm3i( 2, 9,ibeta ), ibeta=1,10)/   &
      0.99917_wp,0.99919_wp,0.99921_wp,0.99915_wp,0.99895_wp,0.99856_wp,0.99792_wp,0.99698_wp,   &
      0.99570_wp,0.99404_wp/
    data (bm3i( 2,10,ibeta ), ibeta=1,10)/   &
      0.99973_wp,0.99973_wp,0.99968_wp,0.99955_wp,0.99928_wp,0.99883_wp,0.99814_wp,0.99716_wp,   &
      0.99584_wp,0.99415_wp/
    data (bm3i( 3, 1,ibeta ), ibeta=1,10)/   &
      0.78358_wp,0.79304_wp,0.81445_wp,0.84105_wp,0.86873_wp,0.89491_wp,0.91805_wp,0.93743_wp,   &
      0.95300_wp,0.96510_wp/
    data (bm3i( 3, 2,ibeta ), ibeta=1,10)/   &
      0.76412_wp,0.77404_wp,0.79635_wp,0.82404_wp,0.85312_wp,0.88101_wp,0.90610_wp,0.92751_wp,   &
      0.94500_wp,0.95879_wp/
    data (bm3i( 3, 3,ibeta ), ibeta=1,10)/   &
      0.74239_wp,0.75182_wp,0.77301_wp,0.79956_wp,0.82809_wp,0.85639_wp,0.88291_wp,0.90658_wp,   &
      0.92683_wp,0.94350_wp/
    data (bm3i( 3, 4,ibeta ), ibeta=1,10)/   &
      0.78072_wp,0.78758_wp,0.80317_wp,0.82293_wp,0.84437_wp,0.86589_wp,0.88643_wp,0.90526_wp,   &
      0.92194_wp,0.93625_wp/
    data (bm3i( 3, 5,ibeta ), ibeta=1,10)/   &
      0.87627_wp,0.88044_wp,0.88981_wp,0.90142_wp,0.91357_wp,0.92524_wp,0.93585_wp,0.94510_wp,   &
      0.95285_wp,0.95911_wp/
    data (bm3i( 3, 6,ibeta ), ibeta=1,10)/   &
      0.95176_wp,0.95371_wp,0.95796_wp,0.96297_wp,0.96792_wp,0.97233_wp,0.97599_wp,0.97880_wp,   &
      0.98072_wp,0.98178_wp/
    data (bm3i( 3, 7,ibeta ), ibeta=1,10)/   &
      0.98453_wp,0.98523_wp,0.98670_wp,0.98833_wp,0.98980_wp,0.99092_wp,0.99160_wp,0.99179_wp,   &
      0.99145_wp,0.99058_wp/
    data (bm3i( 3, 8,ibeta ), ibeta=1,10)/   &
      0.99534_wp,0.99555_wp,0.99597_wp,0.99637_wp,0.99662_wp,0.99663_wp,0.99633_wp,0.99569_wp,   &
      0.99465_wp,0.99318_wp/
    data (bm3i( 3, 9,ibeta ), ibeta=1,10)/   &
      0.99859_wp,0.99864_wp,0.99872_wp,0.99873_wp,0.99860_wp,0.99827_wp,0.99768_wp,0.99679_wp,   &
      0.99555_wp,0.99391_wp/
    data (bm3i( 3,10,ibeta ), ibeta=1,10)/   &
      0.99956_wp,0.99956_wp,0.99953_wp,0.99942_wp,0.99918_wp,0.99875_wp,0.99807_wp,0.99711_wp,   &
      0.99580_wp,0.99412_wp/
    data (bm3i( 4, 1,ibeta ), ibeta=1,10)/   &
      0.84432_wp,0.85223_wp,0.86990_wp,0.89131_wp,0.91280_wp,0.93223_wp,0.94861_wp,0.96172_wp,   &
      0.97185_wp,0.97945_wp/
    data (bm3i( 4, 2,ibeta ), ibeta=1,10)/   &
      0.82299_wp,0.83164_wp,0.85101_wp,0.87463_wp,0.89857_wp,0.92050_wp,0.93923_wp,0.95443_wp,   &
      0.96629_wp,0.97529_wp/
    data (bm3i( 4, 3,ibeta ), ibeta=1,10)/   &
      0.77870_wp,0.78840_wp,0.81011_wp,0.83690_wp,0.86477_wp,0.89124_wp,0.91476_wp,0.93460_wp,   &
      0.95063_wp,0.96316_wp/
    data (bm3i( 4, 4,ibeta ), ibeta=1,10)/   &
      0.76386_wp,0.77233_wp,0.79147_wp,0.81557_wp,0.84149_wp,0.86719_wp,0.89126_wp,0.91275_wp,   &
      0.93116_wp,0.94637_wp/
    data (bm3i( 4, 5,ibeta ), ibeta=1,10)/   &
      0.82927_wp,0.83488_wp,0.84756_wp,0.86346_wp,0.88040_wp,0.89704_wp,0.91257_wp,0.92649_wp,   &
      0.93857_wp,0.94874_wp/
    data (bm3i( 4, 6,ibeta ), ibeta=1,10)/   &
      0.92184_wp,0.92481_wp,0.93136_wp,0.93925_wp,0.94724_wp,0.95462_wp,0.96104_wp,0.96634_wp,   &
      0.97048_wp,0.97348_wp/
    data (bm3i( 4, 7,ibeta ), ibeta=1,10)/   &
      0.97341_wp,0.97457_wp,0.97706_wp,0.97991_wp,0.98260_wp,0.98485_wp,0.98654_wp,0.98760_wp,   &
      0.98801_wp,0.98777_wp/
    data (bm3i( 4, 8,ibeta ), ibeta=1,10)/   &
      0.99192_wp,0.99229_wp,0.99305_wp,0.99385_wp,0.99449_wp,0.99486_wp,0.99487_wp,0.99449_wp,   &
      0.99367_wp,0.99239_wp/
    data (bm3i( 4, 9,ibeta ), ibeta=1,10)/   &
      0.99758_wp,0.99768_wp,0.99787_wp,0.99800_wp,0.99799_wp,0.99777_wp,0.99727_wp,0.99645_wp,   &
      0.99527_wp,0.99369_wp/
    data (bm3i( 4,10,ibeta ), ibeta=1,10)/   &
      0.99926_wp,0.99928_wp,0.99928_wp,0.99921_wp,0.99900_wp,0.99860_wp,0.99795_wp,0.99701_wp,   &
      0.99572_wp,0.99405_wp/
    data (bm3i( 5, 1,ibeta ), ibeta=1,10)/   &
      0.89577_wp,0.90190_wp,0.91522_wp,0.93076_wp,0.94575_wp,0.95876_wp,0.96932_wp,0.97751_wp,   &
      0.98367_wp,0.98820_wp/
    data (bm3i( 5, 2,ibeta ), ibeta=1,10)/   &
      0.87860_wp,0.88547_wp,0.90052_wp,0.91828_wp,0.93557_wp,0.95075_wp,0.96319_wp,0.97292_wp,   &
      0.98028_wp,0.98572_wp/
    data (bm3i( 5, 3,ibeta ), ibeta=1,10)/   &
      0.83381_wp,0.84240_wp,0.86141_wp,0.88425_wp,0.90707_wp,0.92770_wp,0.94510_wp,0.95906_wp,   &
      0.96986_wp,0.97798_wp/
    data (bm3i( 5, 4,ibeta ), ibeta=1,10)/   &
      0.78530_wp,0.79463_wp,0.81550_wp,0.84127_wp,0.86813_wp,0.89367_wp,0.91642_wp,0.93566_wp,   &
      0.95125_wp,0.96347_wp/
    data (bm3i( 5, 5,ibeta ), ibeta=1,10)/   &
      0.79614_wp,0.80332_wp,0.81957_wp,0.84001_wp,0.86190_wp,0.88351_wp,0.90368_wp,0.92169_wp,   &
      0.93718_wp,0.95006_wp/
    data (bm3i( 5, 6,ibeta ), ibeta=1,10)/   &
      0.88192_wp,0.88617_wp,0.89565_wp,0.90728_wp,0.91931_wp,0.93076_wp,0.94107_wp,0.94997_wp,   &
      0.95739_wp,0.96333_wp/
    data (bm3i( 5, 7,ibeta ), ibeta=1,10)/   &
      0.95509_wp,0.95698_wp,0.96105_wp,0.96583_wp,0.97048_wp,0.97460_wp,0.97796_wp,0.98050_wp,   &
      0.98218_wp,0.98304_wp/
    data (bm3i( 5, 8,ibeta ), ibeta=1,10)/   &
      0.98596_wp,0.98660_wp,0.98794_wp,0.98943_wp,0.99074_wp,0.99172_wp,0.99227_wp,0.99235_wp,   &
      0.99192_wp,0.99096_wp/
    data (bm3i( 5, 9,ibeta ), ibeta=1,10)/   &
      0.99581_wp,0.99600_wp,0.99637_wp,0.99672_wp,0.99691_wp,0.99687_wp,0.99653_wp,0.99585_wp,   &
      0.99478_wp,0.99329_wp/
    data (bm3i( 5,10,ibeta ), ibeta=1,10)/   &
      0.99873_wp,0.99878_wp,0.99884_wp,0.99883_wp,0.99869_wp,0.99834_wp,0.99774_wp,0.99684_wp,   &
      0.99558_wp,0.99394_wp/
    data (bm3i( 6, 1,ibeta ), ibeta=1,10)/   &
      0.93335_wp,0.93777_wp,0.94711_wp,0.95764_wp,0.96741_wp,0.97562_wp,0.98210_wp,0.98701_wp,   &
      0.99064_wp,0.99327_wp/
    data (bm3i( 6, 2,ibeta ), ibeta=1,10)/   &
      0.92142_wp,0.92646_wp,0.93723_wp,0.94947_wp,0.96096_wp,0.97069_wp,0.97842_wp,0.98431_wp,   &
      0.98868_wp,0.99186_wp/
    data (bm3i( 6, 3,ibeta ), ibeta=1,10)/   &
      0.88678_wp,0.89351_wp,0.90810_wp,0.92508_wp,0.94138_wp,0.95549_wp,0.96693_wp,0.97578_wp,   &
      0.98243_wp,0.98731_wp/
    data (bm3i( 6, 4,ibeta ), ibeta=1,10)/   &
      0.83249_wp,0.84124_wp,0.86051_wp,0.88357_wp,0.90655_wp,0.92728_wp,0.94477_wp,0.95880_wp,   &
      0.96964_wp,0.97779_wp/
    data (bm3i( 6, 5,ibeta ), ibeta=1,10)/   &
      0.79593_wp,0.80444_wp,0.82355_wp,0.84725_wp,0.87211_wp,0.89593_wp,0.91735_wp,0.93566_wp,   &
      0.95066_wp,0.96255_wp/
    data (bm3i( 6, 6,ibeta ), ibeta=1,10)/   &
      0.84124_wp,0.84695_wp,0.85980_wp,0.87575_wp,0.89256_wp,0.90885_wp,0.92383_wp,0.93704_wp,   &
      0.94830_wp,0.95761_wp/
    data (bm3i( 6, 7,ibeta ), ibeta=1,10)/   &
      0.92721_wp,0.93011_wp,0.93647_wp,0.94406_wp,0.95166_wp,0.95862_wp,0.96460_wp,0.96949_wp,   &
      0.97326_wp,0.97595_wp/
    data (bm3i( 6, 8,ibeta ), ibeta=1,10)/   &
      0.97573_wp,0.97681_wp,0.97913_wp,0.98175_wp,0.98421_wp,0.98624_wp,0.98772_wp,0.98860_wp,   &
      0.98885_wp,0.98847_wp/
    data (bm3i( 6, 9,ibeta ), ibeta=1,10)/   &
      0.99271_wp,0.99304_wp,0.99373_wp,0.99444_wp,0.99499_wp,0.99528_wp,0.99522_wp,0.99477_wp,   &
      0.99390_wp,0.99258_wp/
    data (bm3i( 6,10,ibeta ), ibeta=1,10)/   &
      0.99782_wp,0.99791_wp,0.99807_wp,0.99817_wp,0.99813_wp,0.99788_wp,0.99737_wp,0.99653_wp,   &
      0.99533_wp,0.99374_wp/
    data (bm3i( 7, 1,ibeta ), ibeta=1,10)/   &
      0.95858_wp,0.96158_wp,0.96780_wp,0.97460_wp,0.98073_wp,0.98575_wp,0.98963_wp,0.99252_wp,   &
      0.99463_wp,0.99615_wp/
    data (bm3i( 7, 2,ibeta ), ibeta=1,10)/   &
      0.95091_wp,0.95438_wp,0.96163_wp,0.96962_wp,0.97688_wp,0.98286_wp,0.98751_wp,0.99099_wp,   &
      0.99353_wp,0.99536_wp/
    data (bm3i( 7, 3,ibeta ), ibeta=1,10)/   &
      0.92751_wp,0.93233_wp,0.94255_wp,0.95406_wp,0.96473_wp,0.97366_wp,0.98070_wp,0.98602_wp,   &
      0.98994_wp,0.99278_wp/
    data (bm3i( 7, 4,ibeta ), ibeta=1,10)/   &
      0.88371_wp,0.89075_wp,0.90595_wp,0.92351_wp,0.94028_wp,0.95474_wp,0.96642_wp,0.97544_wp,   &
      0.98220_wp,0.98715_wp/
    data (bm3i( 7, 5,ibeta ), ibeta=1,10)/   &
      0.82880_wp,0.83750_wp,0.85671_wp,0.87980_wp,0.90297_wp,0.92404_wp,0.94195_wp,0.95644_wp,   &
      0.96772_wp,0.97625_wp/
    data (bm3i( 7, 6,ibeta ), ibeta=1,10)/   &
      0.81933_wp,0.82655_wp,0.84279_wp,0.86295_wp,0.88412_wp,0.90449_wp,0.92295_wp,0.93890_wp,   &
      0.95215_wp,0.96281_wp/
    data (bm3i( 7, 7,ibeta ), ibeta=1,10)/   &
      0.89099_wp,0.89519_wp,0.90448_wp,0.91577_wp,0.92732_wp,0.93820_wp,0.94789_wp,0.95616_wp,   &
      0.96297_wp,0.96838_wp/
    data (bm3i( 7, 8,ibeta ), ibeta=1,10)/   &
      0.95886_wp,0.96064_wp,0.96448_wp,0.96894_wp,0.97324_wp,0.97701_wp,0.98004_wp,0.98228_wp,   &
      0.98371_wp,0.98435_wp/
    data (bm3i( 7, 9,ibeta ), ibeta=1,10)/   &
      0.98727_wp,0.98786_wp,0.98908_wp,0.99043_wp,0.99160_wp,0.99245_wp,0.99288_wp,0.99285_wp,   &
      0.99234_wp,0.99131_wp/
    data (bm3i( 7,10,ibeta ), ibeta=1,10)/   &
      0.99621_wp,0.99638_wp,0.99671_wp,0.99700_wp,0.99715_wp,0.99707_wp,0.99670_wp,0.99599_wp,   &
      0.99489_wp,0.99338_wp/
    data (bm3i( 8, 1,ibeta ), ibeta=1,10)/   &
      0.97470_wp,0.97666_wp,0.98064_wp,0.98491_wp,0.98867_wp,0.99169_wp,0.99399_wp,0.99569_wp,   &
      0.99691_wp,0.99779_wp/
    data (bm3i( 8, 2,ibeta ), ibeta=1,10)/   &
      0.96996_wp,0.97225_wp,0.97693_wp,0.98196_wp,0.98643_wp,0.99003_wp,0.99279_wp,0.99482_wp,   &
      0.99630_wp,0.99735_wp/
    data (bm3i( 8, 3,ibeta ), ibeta=1,10)/   &
      0.95523_wp,0.95848_wp,0.96522_wp,0.97260_wp,0.97925_wp,0.98468_wp,0.98888_wp,0.99200_wp,   &
      0.99427_wp,0.99590_wp/
    data (bm3i( 8, 4,ibeta ), ibeta=1,10)/   &
      0.92524_wp,0.93030_wp,0.94098_wp,0.95294_wp,0.96397_wp,0.97317_wp,0.98038_wp,0.98582_wp,   &
      0.98981_wp,0.99270_wp/
    data (bm3i( 8, 5,ibeta ), ibeta=1,10)/   &
      0.87576_wp,0.88323_wp,0.89935_wp,0.91799_wp,0.93583_wp,0.95126_wp,0.96377_wp,0.97345_wp,   &
      0.98072_wp,0.98606_wp/
    data (bm3i( 8, 6,ibeta ), ibeta=1,10)/   &
      0.83078_wp,0.83894_wp,0.85705_wp,0.87899_wp,0.90126_wp,0.92179_wp,0.93950_wp,0.95404_wp,   &
      0.96551_wp,0.97430_wp/
    data (bm3i( 8, 7,ibeta ), ibeta=1,10)/   &
      0.85727_wp,0.86294_wp,0.87558_wp,0.89111_wp,0.90723_wp,0.92260_wp,0.93645_wp,0.94841_wp,   &
      0.95838_wp,0.96643_wp/
    data (bm3i( 8, 8,ibeta ), ibeta=1,10)/   &
      0.93337_wp,0.93615_wp,0.94220_wp,0.94937_wp,0.95647_wp,0.96292_wp,0.96840_wp,0.97283_wp,   &
      0.97619_wp,0.97854_wp/
    data (bm3i( 8, 9,ibeta ), ibeta=1,10)/   &
      0.97790_wp,0.97891_wp,0.98105_wp,0.98346_wp,0.98569_wp,0.98751_wp,0.98879_wp,0.98950_wp,   &
      0.98961_wp,0.98912_wp/
    data (bm3i( 8,10,ibeta ), ibeta=1,10)/   &
      0.99337_wp,0.99367_wp,0.99430_wp,0.99493_wp,0.99541_wp,0.99562_wp,0.99551_wp,0.99501_wp,   &
      0.99410_wp,0.99274_wp/
    data (bm3i( 9, 1,ibeta ), ibeta=1,10)/   &
      0.98470_wp,0.98594_wp,0.98844_wp,0.99106_wp,0.99334_wp,0.99514_wp,0.99650_wp,0.99749_wp,   &
      0.99821_wp,0.99872_wp/
    data (bm3i( 9, 2,ibeta ), ibeta=1,10)/   &
      0.98184_wp,0.98330_wp,0.98624_wp,0.98934_wp,0.99205_wp,0.99420_wp,0.99582_wp,0.99701_wp,   &
      0.99787_wp,0.99848_wp/
    data (bm3i( 9, 3,ibeta ), ibeta=1,10)/   &
      0.97288_wp,0.97498_wp,0.97927_wp,0.98385_wp,0.98789_wp,0.99113_wp,0.99360_wp,0.99541_wp,   &
      0.99673_wp,0.99766_wp/
    data (bm3i( 9, 4,ibeta ), ibeta=1,10)/   &
      0.95403_wp,0.95741_wp,0.96440_wp,0.97202_wp,0.97887_wp,0.98444_wp,0.98872_wp,0.99190_wp,   &
      0.99421_wp,0.99586_wp/
    data (bm3i( 9, 5,ibeta ), ibeta=1,10)/   &
      0.91845_wp,0.92399_wp,0.93567_wp,0.94873_wp,0.96076_wp,0.97079_wp,0.97865_wp,0.98457_wp,   &
      0.98892_wp,0.99206_wp/
    data (bm3i( 9, 6,ibeta ), ibeta=1,10)/   &
      0.86762_wp,0.87533_wp,0.89202_wp,0.91148_wp,0.93027_wp,0.94669_wp,0.96013_wp,0.97062_wp,   &
      0.97855_wp,0.98441_wp/
    data (bm3i( 9, 7,ibeta ), ibeta=1,10)/   &
      0.84550_wp,0.85253_wp,0.86816_wp,0.88721_wp,0.90671_wp,0.92490_wp,0.94083_wp,0.95413_wp,   &
      0.96481_wp,0.97314_wp/
    data (bm3i( 9, 8,ibeta ), ibeta=1,10)/   &
      0.90138_wp,0.90544_wp,0.91437_wp,0.92513_wp,0.93602_wp,0.94615_wp,0.95506_wp,0.96258_wp,   &
      0.96868_wp,0.97347_wp/
    data (bm3i( 9, 9,ibeta ), ibeta=1,10)/   &
      0.96248_wp,0.96415_wp,0.96773_wp,0.97187_wp,0.97583_wp,0.97925_wp,0.98198_wp,0.98394_wp,   &
      0.98514_wp,0.98559_wp/
    data (bm3i( 9,10,ibeta ), ibeta=1,10)/   &
      0.98837_wp,0.98892_wp,0.99005_wp,0.99127_wp,0.99232_wp,0.99306_wp,0.99339_wp,0.99328_wp,   &
      0.99269_wp,0.99161_wp/
    data (bm3i(10, 1,ibeta ), ibeta=1,10)/   &
      0.99080_wp,0.99158_wp,0.99311_wp,0.99471_wp,0.99607_wp,0.99715_wp,0.99795_wp,0.99853_wp,   &
      0.99895_wp,0.99925_wp/
    data (bm3i(10, 2,ibeta ), ibeta=1,10)/   &
      0.98910_wp,0.99001_wp,0.99182_wp,0.99371_wp,0.99533_wp,0.99661_wp,0.99757_wp,0.99826_wp,   &
      0.99876_wp,0.99912_wp/
    data (bm3i(10, 3,ibeta ), ibeta=1,10)/   &
      0.98374_wp,0.98506_wp,0.98772_wp,0.99051_wp,0.99294_wp,0.99486_wp,0.99630_wp,0.99736_wp,   &
      0.99812_wp,0.99866_wp/
    data (bm3i(10, 4,ibeta ), ibeta=1,10)/   &
      0.97238_wp,0.97453_wp,0.97892_wp,0.98361_wp,0.98773_wp,0.99104_wp,0.99354_wp,0.99538_wp,   &
      0.99671_wp,0.99765_wp/
    data (bm3i(10, 5,ibeta ), ibeta=1,10)/   &
      0.94961_wp,0.95333_wp,0.96103_wp,0.96941_wp,0.97693_wp,0.98303_wp,0.98772_wp,0.99119_wp,   &
      0.99371_wp,0.99551_wp/
    data (bm3i(10, 6,ibeta ), ibeta=1,10)/   &
      0.90943_wp,0.91550_wp,0.92834_wp,0.94275_wp,0.95608_wp,0.96723_wp,0.97600_wp,0.98263_wp,   &
      0.98751_wp,0.99103_wp/
    data (bm3i(10, 7,ibeta ), ibeta=1,10)/   &
      0.86454_wp,0.87200_wp,0.88829_wp,0.90749_wp,0.92630_wp,0.94300_wp,0.95687_wp,0.96785_wp,   &
      0.97626_wp,0.98254_wp/
    data (bm3i(10, 8,ibeta ), ibeta=1,10)/   &
      0.87498_wp,0.88048_wp,0.89264_wp,0.90737_wp,0.92240_wp,0.93642_wp,0.94877_wp,0.95917_wp,   &
      0.96762_wp,0.97429_wp/
    data (bm3i(10, 9,ibeta ), ibeta=1,10)/   &
      0.93946_wp,0.94209_wp,0.94781_wp,0.95452_wp,0.96111_wp,0.96704_wp,0.97203_wp,0.97602_wp,   &
      0.97900_wp,0.98106_wp/
    data (bm3i(10,10,ibeta ), ibeta=1,10)/   &
      0.97977_wp,0.98071_wp,0.98270_wp,0.98492_wp,0.98695_wp,0.98858_wp,0.98970_wp,0.99027_wp,   &
      0.99026_wp,0.98968_wp/

    ! *** end of data statements *************************************

    !----------------------------------------------------------
    ! Start calculations
    !----------------------------------------------------------
    ! Constants and parameters

    sqrttwo = sqrt(two)
    dlgsqt2 = one / log( sqrttwo )

    esat01   = exp( 0.125_wp * xxlsgat * xxlsgat )
    esac01   = exp( 0.125_wp * xxlsgac * xxlsgac )

    esat04  = esat01 ** 4
    esac04  = esac01 ** 4

    esat05  = esat04 * esat01
    esac05  = esac04 * esac01

    esat08  = esat04 * esat04
    esac08  = esac04 * esac04

    esat09  = esat08 * esat01
    esac09  = esac08 * esac01

    esat16  = esat08 * esat08
    esac16  = esac08 * esac08

    esat20  = esat16 * esat04
    esac20  = esac16 * esac04

    esat24  = esat20 * esat04
    esac24  = esac20 * esac04

    esat25  = esat20 * esat05
    esac25  = esac20 * esac05

    esat36  = esat20 * esat16
    esac36  = esac20 * esac16

    esat49  = esat24 * esat25

    esat64  = esat20 * esat20 * esat24
    esac64  = esac20 * esac20 * esac24

    esat100 = esat64 * esat36

    dgat2   = dgatk * dgatk
    dgat3   = dgatk * dgatk * dgatk
    dgac2   = dgacc * dgacc
    dgac3   = dgacc * dgacc * dgacc

    sqdgat  = sqrt( dgatk )
    sqdgac  = sqrt( dgacc )
    sqdgat5 = dgat2 * sqdgat
    sqdgac5 = dgac2 * sqdgac
    sqdgat7 = dgat3 * sqdgat

    !------------------------------------------------------------------
    ! For the free molecular regime:  page h.3 of whitby et al. (1991)
    !------------------------------------------------------------------
    r       = sqdgac / sqdgat
    r2      = r * r
    r3      = r2 * r
    rx4     = r2 * r2
    ri1     = one / r
    ri2     = one / r2
    ri3     = one / r3
    kngat   = two * lamda / dgatk
    kngac   = two * lamda / dgacc

    ! Calculate ratio of geometric mean diameters

    rat = dgacc / dgatk

    ! Trap subscripts for bm0 and bm0i, between 1 and 10.
    ! See page h.5 of whitby et al. (1991)

    n2n = max( 1, min( 10,   nint( 4.0_wp * ( sgatk - 0.75_wp ) ) ) )
    n2a = max( 1, min( 10,   nint( 4.0_wp * ( sgacc - 0.75_wp ) ) ) )
    n1  = max( 1, min( 10,   1 + nint( dlgsqt2 * log( rat ) ) ) )

    !----------------------------------------------------
    ! Intermodal coagulation, set up for zeroeth moment
    !----------------------------------------------------
    ! Near-continuum form:  equation h.10a of whitby et al. (1991)

    coagnc0 = knc * (   &
      two + a * ( kngat * ( esat04 + r2 * esat16 * esac04 )   &
      + kngac * ( esac04 + ri2 * esac16 * esat04 ) )   &
      + ( r2 + ri2 ) * esat04 * esac04  )

    ! Free-molecular form:  equation h.7a of whitby et al. (1991)

    coagfm0 = kfmatac * sqdgat * bm0ij(n1,n2n,n2a) * (   &
      esat01 + r * esac01 + two * r2 * esat01 * esac04   &
      + rx4 * esat09 * esac16 + ri3 * esat16 * esac09   &
      + two * ri1 * esat04 + esac01  )

    ! Loss to accumulation mode, harmonic mean

    qn12 = coagnc0 * coagfm0 / ( coagnc0 + coagfm0 )

    !-----------------------------------------------------------------
    ! Intermodal coagulation, zeroeth moment, set up for third moment
    !-----------------------------------------------------------------
    ! Near-continuum form: equation h.10b of whitby et al. (1991)

    coagnc3 = knc * dgat3 * (   &
      two * esat36   &
      + a * kngat * ( esat16 + r2 * esat04 * esac04 )   &
      + a * kngac * ( esat36 * esac04 + ri2 * esat64 * esac16 )   &
      + r2 * esat16 * esac04 + ri2 * esat64 * esac04 )

    ! Free-molecular form: equation h.7b of whitby et al. (1991)

    coagfm3 = kfmatac * sqdgat7 * bm3i( n1, n2n, n2a ) * (   &
      esat49   &
      +  r * esat36  * esac01   &
      + two * r2 * esat25  * esac04   &
      + rx4 * esat09  * esac16   &
      + ri3 * esat100 * esac09   &
      + two * ri1 * esat64  * esac01 )

    ! Gain by accumulation mode = loss from aitken mode. Harmonic mean

    qv12 = coagnc3 * coagfm3 / ( coagnc3 + coagfm3 )

    !-----------------------------------------------------
    ! Intramodal coagulation
    !-----------------------------------------------------
    ! aitken mode
    !--------------
    ! Near-continuum form: equation h.12a of whitby et al. (1991)
    coagnc_at = knc * (one + esat08 + a * kngat * (esat20 + esat04))

    ! Free-molecular form: equation h.11a of whitby et al. (1991)
    coagfm_at = kfmat * sqdgat * bm0(n2n) *  ( esat01 + esat25 + two * esat05 )

    ! Harmonic mean
    qn11 = coagfm_at * coagnc_at / ( coagfm_at + coagnc_at )

    !-------------------
    ! accumulation mode
    !-------------------
    ! Near-continuum form: equation h.12a of whitby et al. (1991)
    coagnc_ac = knc * (one + esac08 + a * kngac * (esac20 + esac04))

    ! Free-molecular form: equation h.11a of whitby et al. (1991)
    coagfm_ac = kfmac * sqdgac * bm0(n2a) * ( esac01 + esac25 + two * esac05 )

    ! Harmonic mean
    qn22 = coagfm_ac * coagnc_ac / ( coagfm_ac + coagnc_ac )

end  subroutine getcoags

end module modal_aero_coag
