
module mo_setsox

  use shr_kind_mod, only : r8 => shr_kind_r8
  use cam_logfile,  only : iulog
    use spmd_utils,   only : masterproc

  private
  public :: sox_inti, setsox
  public :: has_sox

  save
  logical            ::  inv_o3
  integer            ::  id_msa

  integer :: id_so2, id_nh3, id_hno3, id_h2o2, id_o3, id_ho2
  integer :: id_so4, id_h2so4

  logical :: has_sox = .true.
  logical :: inv_so2, inv_nh3, inv_hno3, inv_h2o2, inv_ox, inv_nh4no3, inv_ho2

  logical :: cloud_borne = .false.
  logical :: modal_aerosols = .false.

contains

!-----------------------------------------------------------------------      
!-----------------------------------------------------------------------      
  subroutine sox_inti
    !-----------------------------------------------------------------------      
    !	... initialize the hetero sox routine
    !-----------------------------------------------------------------------      

    use mo_chem_utls, only : get_spc_ndx, get_inv_ndx
    use spmd_utils,   only : masterproc
    use cam_history,  only : addfld
    use cam_history,  only : add_default
    use ppgrid,       only : pver
    use phys_control, only : phys_getopts
    use sox_cldaero_mod, only : sox_cldaero_init

    implicit none

    logical :: history_aerosol   ! Output aerosol diagnostics
    logical :: history_verbose   ! produce verbose history output

    call phys_getopts( &
         history_aerosol_out = history_aerosol, &
         history_verbose_out = history_verbose, &
         prog_modal_aero_out=modal_aerosols )

    cloud_borne = modal_aerosols

    !-----------------------------------------------------------------
    !       ... get species indicies
    !-----------------------------------------------------------------
    
    if (cloud_borne) then
       id_h2so4 = get_spc_ndx( 'H2SO4' )
    else
       id_so4 = get_spc_ndx( 'SO4' )
    endif
    id_msa = get_spc_ndx( 'MSA' )

    inv_so2 = .false.
    id_so2 = get_inv_ndx( 'SO2' )
    inv_so2 = id_so2 > 0
    if ( .not. inv_so2 ) then
       id_so2 = get_spc_ndx( 'SO2' )
    endif

    inv_NH3 = .false.
    id_NH3 = get_inv_ndx( 'NH3' )
    inv_NH3 = id_NH3 > 0
    if ( .not. inv_NH3 ) then
       id_NH3 = get_spc_ndx( 'NH3' )
    endif

    inv_HNO3 = .false.
    id_HNO3 = get_inv_ndx( 'HNO3' )
    inv_HNO3 = id_hno3 > 0
    if ( .not. inv_HNO3 ) then
       id_HNO3 = get_spc_ndx( 'HNO3' )
    endif

    inv_H2O2 = .false.
    id_H2O2 = get_inv_ndx( 'H2O2' )
    inv_H2O2 = id_H2O2 > 0
    if ( .not. inv_H2O2 ) then
       id_H2O2 = get_spc_ndx( 'H2O2' )
    endif

    inv_HO2 = .false.
    id_HO2 = get_inv_ndx( 'HO2' )
    inv_HO2 = id_HO2 > 0
    if ( .not. inv_HO2 ) then
       id_HO2 = get_spc_ndx( 'HO2' )
    endif

    inv_o3 = get_inv_ndx( 'O3' ) > 0
    if (inv_o3) then
       id_o3 = get_inv_ndx( 'O3' )
    else
       id_o3 = get_spc_ndx( 'O3' )
    endif
    inv_ho2 = get_inv_ndx( 'HO2' ) > 0
    if (inv_ho2) then
       id_ho2 = get_inv_ndx( 'HO2' )
    else
       id_ho2 = get_spc_ndx( 'HO2' )
    endif

    has_sox = (id_so2>0) .and. (id_h2o2>0) .and. (id_o3>0) .and. (id_ho2>0)
    if (cloud_borne) then
       has_sox = has_sox .and. (id_h2so4>0)
    else
       has_sox = has_sox .and. (id_so4>0) .and. (id_nh3>0)
    endif

    if (masterproc) then
       write(iulog,*) 'sox_inti: has_sox = ',has_sox
    endif

    if( has_sox ) then
       if (masterproc) then
          write(iulog,*) '-----------------------------------------'
          write(iulog,*) 'mozart will do sox aerosols'
          write(iulog,*) '-----------------------------------------'
       endif
    else 
       return
    end if

    call addfld( 'XPH_LWC',(/ 'lev' /), 'A','kg/kg', 'pH value multiplied by lwc')
    if ( history_aerosol .and. history_verbose ) then    
       call add_default ('XPH_LWC', 1, ' ') 
    endif

    call sox_cldaero_init()

  end subroutine sox_inti
  
!-----------------------------------------------------------------------      
!-----------------------------------------------------------------------      
  subroutine setsox(   &
       ncol,   lchnk,  loffset,   dtime,  & ! in
       press,  pdel,   tfld,      mbar,   & ! in
       lwc,    cldfrc, cldnum,            & ! in
       xhnm,   invariants,                & ! in
       qcw,    qin                        ) ! inout

    !-----------------------------------------------------------------------      
    !          ... Compute heterogeneous reactions of SOX
    !
    !       (0) using initial PH to calculate PH
    !           (a) HENRYs law constants
    !           (b) PARTIONING
    !           (c) PH values
    !
    !       (1) using new PH to repeat
    !           (a) HENRYs law constants
    !           (b) PARTIONING
    !           (c) REACTION rates
    !           (d) PREDICTION
    !-----------------------------------------------------------------------      
    !
    use ppgrid,    only : pcols, pver
    use chem_mods, only : gas_pcnst, nfs
    use chem_mods,    only : adv_mass
    use physconst,    only : mwdry, gravit
    use mo_constants, only : pi
    use cam_history,  only : outfld
    use sox_cldaero_mod, only : sox_cldaero_update, sox_cldaero_create_obj, sox_cldaero_destroy_obj
    use cldaero_mod,     only : cldaero_conc_t
    use phys_control, only : phys_getopts
    use cam_abortutils,   only: endrun

    !
    implicit none
    !
    !-----------------------------------------------------------------------      
    !      ... Dummy arguments
    !-----------------------------------------------------------------------      
    integer,          intent(in)    :: ncol              ! num of columns in chunk
    integer,          intent(in)    :: lchnk             ! chunk id
    integer,          intent(in)    :: loffset           ! offset of chem tracers in the advected tracers array
    real(r8),         intent(in)    :: dtime             ! time step [sec]
    real(r8),         intent(in)    :: press(:,:)        ! midpoint pressure [Pa]
    real(r8),         intent(in)    :: pdel(:,:)         ! pressure thickness of levels [Pa]
    real(r8),         intent(in)    :: tfld(:,:)         ! temperature [K]
    real(r8),         intent(in)    :: mbar(:,:)         ! mean wet atmospheric mass [amu or g/mol]
    real(r8), target, intent(in)    :: lwc(:,:)          ! cloud liquid water content [kg/kg]
    real(r8), target, intent(in)    :: cldfrc(:,:)       ! cloud fraction [fraction]
    real(r8),         intent(in)    :: cldnum(:,:)       ! droplet number concentration [#/kg]
    real(r8),         intent(in)    :: xhnm(:,:)         ! total atms density [#/cm**3]
    real(r8),         intent(in)    :: invariants(:,:,:) ! invariant density [molecules/cm**3]
    real(r8), target, intent(inout) :: qcw(:,:,:)        ! cloud-borne aerosol [vmr]
    real(r8),         intent(inout) :: qin(:,:,:)        ! transported species [vmr]

    !-----------------------------------------------------------------------      
    !      ... Local variables
    !
    !           xhno3 ... in mixing ratio
    !-----------------------------------------------------------------------      
    integer,  parameter :: itermax = 20  ! maximum number of iterations
    real(r8), parameter :: ph0 = 5.0_r8  ! Initial PH values
    real(r8), parameter :: const0 = 1.e3_r8/6.023e23_r8
    real(r8), parameter :: xa0 = 11._r8
    real(r8), parameter :: xb0 = -.1_r8
    real(r8), parameter :: xa1 = 1.053_r8
    real(r8), parameter :: xb1 = -4.368_r8
    real(r8), parameter :: xa2 = 1.016_r8
    real(r8), parameter :: xb2 = -2.54_r8
    real(r8), parameter :: xa3 = .816e-32_r8
    real(r8), parameter :: xb3 = .259_r8

    real(r8), parameter :: kh0 = 9.e3_r8            ! HO2(g)          -> Ho2(a)
    real(r8), parameter :: kh1 = 2.05e-5_r8         ! HO2(a)          -> H+ + O2-
    real(r8), parameter :: kh2 = 8.6e5_r8           ! HO2(a) + ho2(a) -> h2o2(a) + o2
    real(r8), parameter :: kh3 = 1.e8_r8            ! HO2(a) + o2-    -> h2o2(a) + o2
    real(r8), parameter :: Ra = 8314._r8/101325._r8 ! universal constant   (atm)/(M-K)
    real(r8), parameter :: xkw = 1.e-14_r8          ! water acidity

    !
    real(r8) :: xdelso4hp(ncol,pver)
    real(r8) :: xdelso4hp_ik
    real(r8) :: xphlwc(ncol,pver)

    integer  :: k, i, iter, file
    real(r8) :: wrk, delta
    real(r8) :: xph0, aden, xk, xe, x2
    real(r8) :: tz, xl, px, qz, pz, es, qs, patm
    real(r8) :: Eso2, Eso4, Ehno3, Eco2, Eh2o, Enh3
    real(r8) :: so2g, h2o2g, co2g, o3g
    real(r8) :: hno3a, nh3a, so2a, h2o2a, co2a, o3a
    real(r8) :: rah2o2, rao3, pso4, ccc
    real(r8) :: cnh3, chno3, com, com1, com2, xra

    real(r8) :: hno3g(ncol,pver), nh3g(ncol,pver)
    !
    !-----------------------------------------------------------------------      
    !            for Ho2(g) -> H2o2(a) formation 
    !            schwartz JGR, 1984, 11589
    !-----------------------------------------------------------------------      
    real(r8) :: kh4    ! kh2+kh3
    real(r8) :: xam    ! air density /cm3
    real(r8) :: ho2s   ! ho2s = ho2(a)+o2-
    real(r8) :: r1h2o2 ! prod(h2o2) by ho2 in mole/L(w)/s
    real(r8) :: r2h2o2 ! prod(h2o2) by ho2 in mix/s

    real(r8), dimension(ncol,pver)  ::             &
         xhno3, xh2o2, xso2, xso4, xno3, &
         xnh3, xnh4, xo3,         & ! x***: atom concentration [kg/L]
         cfact, &
         xph, xho2,         &
         xh2so4, xmsa, xso4_init, &
         hehno3, &            ! henry law const for hno3
         heh2o2, &            ! henry law const for h2o2
         heso2,  &            ! henry law const for so2
         henh3,  &            ! henry law const for nh3
         heo3              !!,   &            ! henry law const for o3

    real(r8) :: patm_x

    real(r8), dimension(ncol)  :: work1
    logical :: converged

    real(r8), pointer :: xso4c(:,:)
    real(r8), pointer :: xnh4c(:,:)
    real(r8), pointer :: xno3c(:,:)
    type(cldaero_conc_t), pointer :: cldconc

    real(r8) :: fact1_hno3, fact2_hno3, fact3_hno3
    real(r8) :: fact1_so2, fact2_so2, fact3_so2, fact4_so2
    real(r8) :: fact1_nh3, fact2_nh3, fact3_nh3
    real(r8) :: tmp_hp, tmp_hso3, tmp_hco3, tmp_nh4, tmp_no3
    real(r8) :: tmp_oh, tmp_so3, tmp_so4
    real(r8) :: tmp_neg, tmp_pos
    real(r8) :: yph, yph_lo, yph_hi
    real(r8) :: ynetpos, ynetpos_lo, ynetpos_hi


    !-----------------------------------------------------------------
    !       ... NOTE: The press array is in pascals and must be
    !                 mutiplied by 10 to yield dynes/cm**2.
    !-----------------------------------------------------------------
    !==================================================================
    !       ... First set the PH
    !==================================================================
    !      ... Initial values
    !           The values of so2, so4 are after (1) SLT, and CHEM
    !-----------------------------------------------------------------
    ! initial PH value, in H+ concentration
    xph0 = 10._r8**(-ph0) 

    ! calculate total atms density [kg/L]
    do k = 1,pver
       cfact(:,k) = xhnm(:,k)     &          ! /cm3(a)  
            * 1.e6_r8             &          ! /m3(a)
            * 1.38e-23_r8/287._r8 &          ! Kg(a)/m3(a)
            * 1.e-3_r8                       ! Kg(a)/L(a)
    enddo

    if ( inv_so2 .or. id_hno3>0 .or. inv_h2o2 .or. id_nh3>0 .or. inv_o3 &
                 .or. (.not. inv_ho2) .or. (.not. cloud_borne) .or. id_msa>0) then
        call endrun('FORTRAN refactoring: Only keep the code for default MAM4. &
             The following options are removed:  id_nh3>0  id_hno3>0  id_msa>0 &
             inv_h2o2=.T. inv_so2=.T.  inv_o3=.T. inv_ho2=.F. cloud_borne=.F. ')
    endif

    ! initialize species concentrations
    cldconc => sox_cldaero_create_obj( cldfrc,qcw,lwc, cfact, ncol, loffset )
    xso4c => cldconc%so4c
    xnh4c => cldconc%nh4c
    xno3c => cldconc%no3c

    xso4(:,:) = 0._r8
    xno3(:,:) = 0._r8
    xnh4(:,:) = 0._r8
    xnh3(:,:) = 0._r8
    xhno3(:,:)= 0._r8
    do k = 1,pver
       xph(:,k) = xph0                                ! initial PH value
       xso2 (:,k) = qin(:,k,id_so2)                 
       xh2o2 (:,k) = qin(:,k,id_h2o2)               
       xo3  (:,k) = qin(:,k,id_o3)                  
       xho2 (:,k) = invariants(:,k,id_ho2)/xhnm(:,k)
       xh2so4(:,k) = qin(:,k,id_h2so4)
    enddo
    
    ! NO3, are not incorporated in MAM4, remove the related code
    ! assign 0 to input variables needed for subroutine sox_cldaero_update
    ! this assignment can be removed when incorporating refaction of
    ! sox_cldaero_update
    hno3g(:,:) = 0._r8
    xmsa(:,:)  = 0._r8
    

    !-----------------------------------------------------------------
    !       ... Temperature dependent Henry constants
    !-----------------------------------------------------------------
    ver_loop0: do k = 1,pver                               !! pver loop for STEP 0
       col_loop0: do i = 1,ncol
          
          if (cloud_borne .and. cldfrc(i,k)>0._r8) then
             xso4(i,k) = xso4c(i,k) / cldfrc(i,k)
          endif

          ! cloud liquid water content
          xl = cldconc%xlwc(i,k)
          if( xl >= 1.e-8_r8 ) then
 
             call calc_ph_values(               &
                tfld(i,k), press(i,k), xl,      & ! in
                xso2(i,k), xso4(i,k), xhnm(i,k),  cldconc%so4_fact, & ! in
                Ra, xkw, const0,                & ! in
                converged, xph(i,k)             ) ! out
             if( .not. converged ) then
                write(iulog,*) 'setsox: pH failed to converge @ (',i,',',k,').'
             endif

          else
             xph(i,k) =  1.e-7_r8
          endif
       enddo col_loop0
    enddo ver_loop0 ! end pver loop for STEP 0

    !==============================================================
    !          ... Now use the actual PH
    !==============================================================
    ver_loop1: do k = 1,pver
       col_loop1: do i = 1,ncol
          work1(i) = 1._r8 / tfld(i,k) - 1._r8 / 298._r8
          tz = tfld(i,k)

          xl = cldconc%xlwc(i,k)

          patm = press(i,k)/101300._r8        ! press is in pascal
          xam  = press(i,k)/(1.38e-23_r8*tz)  ! air density /M3

          !-----------------------------------------------------------------
          !        ... h2o2
          !-----------------------------------------------------------------
          xk = 7.4e4_r8   *EXP( 6621._r8*work1(i) )
          xe = 2.2e-12_r8 *EXP(-3730._r8*work1(i) )
          heh2o2(i,k)  = xk*(1._r8 + xe/xph(i,k))

          !-----------------------------------------------------------------
          !         ... so2
          !-----------------------------------------------------------------
          xk = 1.23_r8  *EXP( 3120._r8*work1(i) )
          xe = 1.7e-2_r8*EXP( 2090._r8*work1(i) )
          x2 = 6.0e-8_r8*EXP( 1120._r8*work1(i) )

          wrk = xe/xph(i,k)
          heso2(i,k)  = xk*(1._r8 + wrk*(1._r8 + x2/xph(i,k)))

          !-----------------------------------------------------------------
          !        ... o3
          !-----------------------------------------------------------------
          xk = 1.15e-2_r8 *EXP( 2560._r8*work1(i) )
          heo3(i,k) = xk

          !------------------------------------------------------------------------
          !       ... for Ho2(g) -> H2o2(a) formation 
          !           schwartz JGR, 1984, 11589
          !------------------------------------------------------------------------
          kh4 = (kh2 + kh3*kh1/xph(i,k)) / ((1._r8 + kh1/xph(i,k))**2)
          ho2s = kh0*xho2(i,k)*patm*(1._r8 + kh1/xph(i,k))  ! ho2s = ho2(a)+o2-
          r1h2o2 = kh4*ho2s*ho2s                         ! prod(h2o2) in mole/L(w)/s

          r2h2o2 = r1h2o2*xl        &    ! mole/L(w)/s   * L(w)/fm3(a) = mole/fm3(a)/s
                 / const0*1.e+6_r8  &    ! FIXME: correct a bug here ????
                 / xam                   ! /cm3(a)/s    / air-den     = mix-ratio/s

          !-----------------------------------------------
          !       ... Partioning 
          !-----------------------------------------------

          !------------------------------------------------------------------------
          !        ... h2o2
          !------------------------------------------------------------------------
          px = heh2o2(i,k) * Ra * tz * xl
          h2o2g =  xh2o2(i,k)/(1._r8+ px)

          !------------------------------------------------------------------------
          !         ... so2
          !------------------------------------------------------------------------
          px = heso2(i,k) * Ra * tz * xl
          so2g =  xso2(i,k)/(1._r8+ px)

          !------------------------------------------------------------------------
          !         ... o3
          !------------------------------------------------------------------------
          px = heo3(i,k) * Ra * tz * xl
          o3g =  xo3(i,k)/(1._r8+ px)

          !-----------------------------------------------
          !       ... Aqueous phase reaction rates
          !           SO2 + H2O2 -> SO4
          !           SO2 + O3   -> SO4
          !-----------------------------------------------

          !------------------------------------------------------------------------
          !       ... S(IV) (HSO3) + H2O2
          !------------------------------------------------------------------------
          rah2o2 = 8.e4_r8 * EXP( -3650._r8*work1(i) )  &
               / (.1_r8 + xph(i,k))

          !------------------------------------------------------------------------
          !        ... S(IV)+ O3
          !------------------------------------------------------------------------
          rao3   = 4.39e11_r8 * EXP(-4131._r8/tz)  &
               + 2.56e3_r8  * EXP(-996._r8 /tz) /xph(i,k)

          !-----------------------------------------------------------------
          !       ... Prediction after aqueous phase
          !       so4
          !       When Cloud is present 
          !   
          !       S(IV) + H2O2 = S(VI)
          !       S(IV) + O3   = S(VI)
          !
          !       reference:
          !           (1) Seinfeld
          !           (2) Benkovitz
          !-----------------------------------------------------------------
          
          !............................
          !       S(IV) + H2O2 = S(VI)
          !............................
          
          if (xl >= 1.e-8_r8) then    !! WHEN CLOUD IS PRESENTED          

             call calc_sox_aqueous( modal_aerosols,       &
                rah2o2, h2o2g, so2g, o3g,      rao3,   &
                patm, dtime, work1(i), xl, const0, &
                xhnm(i,k), heo3(i,k), heso2(i,k),      &
                xso2(i,k), xso4(i,k), xso4_init(i,k), xh2o2(i,k), &
                xdelso4hp_ik)
             xdelso4hp(i,k) = xdelso4hp_ik

          endif !! WHEN CLOUD IS PRESENTED

       end do col_loop1
    end do ver_loop1

    call sox_cldaero_update( &
         ncol, lchnk, loffset, dtime, mbar, pdel, press, tfld, cldnum, cldfrc, cfact, cldconc%xlwc, &
         xdelso4hp, xh2so4, xso4, xso4_init, nh3g, hno3g, xnh3, xhno3, xnh4c,  xno3c, xmsa, xso2, xh2o2, qcw, qin )
    
    xphlwc(:,:) = 0._r8
    do k = 1, pver
       do i = 1, ncol
          if (cldfrc(i,k)>=1.e-5_r8 .and. lwc(i,k)>=1.e-8_r8) then
             xphlwc(i,k) = -1._r8*log10(xph(i,k)) * lwc(i,k)
          endif
       end do
    end do
    call outfld( 'XPH_LWC', xphlwc(:ncol,:), ncol , lchnk )

    call sox_cldaero_destroy_obj(cldconc)

  end subroutine setsox 

!===========================================================================
  subroutine calc_ph_values(                    &
                temperature, pressure, xlwc,    & ! in
                xso2, xso4, xhnm,  so4_fact,    & ! in
                Ra,   xkw,  const0,             & ! in
                converged, xph                  ) ! out
!---------------------------------------------------------------------------
! calculate PH value and H+ concentration
!
! 21-mar-2011 changes by rce
! now uses bisection method to solve the electro-neutrality equation
! 3-mode aerosols (where so4 is assumed to be nh4hso4)
!       old code set xnh4c = so4c
!       new code sets xnh4c = 0, then uses a -1 charge (instead of -2)
 !      for so4 when solving the electro-neutrality equation
!---------------------------------------------------------------------------
    implicit none

    real(r8),  intent(in) :: temperature        ! temperature [K]
    real(r8),  intent(in) :: pressure           ! pressure [Pa]
    real(r8),  intent(in) :: xso2               ! SO2 [kg/L]
    real(r8),  intent(in) :: xso4               ! SO4 [kg/L]
    real(r8),  intent(in) :: xhnm               ! [kg/L]
    real(r8),  intent(in) :: xlwc               ! cloud LWC [kg/L]
    real(r8),  intent(in) :: so4_fact           ! factor for SO4
    real(r8),  intent(in) :: Ra                 ! constant parameter
    real(r8),  intent(in) :: xkw                ! constant parameter
    real(r8),  intent(in) :: const0             ! constant parameter

    logical,  intent(out) :: converged          ! if the method converge
    real(r8), intent(out) :: xph                ! H+ ions concentration [kg/L]


    ! local variables
    integer   :: iter  ! iteration number
    real(r8)  :: yph_lo, yph_hi, yph    ! pH values, lower and upper bounds
    real(r8)  :: ynetpos_lo, ynetpos_hi ! lower and upper bounds of ynetpos
    real(r8)  :: t_factor       ! working variables to convert temperature
    real(r8)  :: patm           ! pressure in atm
    real(r8)  :: xk, xe, x2     ! working variables
    real(r8)  :: fact1_so2, fact2_so2, fact3_so2, fact4_so2  ! SO2 factors
    real(r8)  :: Eh2o, Eco2, Eso4 ! effects of species [1/cm3]
    real(r8)  :: ynetpos        ! net positive ions

    integer,  parameter :: itermax = 20  ! maximum number of iterations
    real(r8), parameter :: co2g = 330.e-6_r8    !330 ppm = 330.e-6 atm


    !----------------------------------------
    ! calculate variables before iterating
    !----------------------------------------

    !  working variable to convert to 25 degC (1/T - 1/[298K])
    t_factor = 1._r8 / temperature - 1._r8 / 298._r8
    ! pressure in atm
    patm = .01_r8*pressure/1013._r8

    !----------------------------------------
    ! effect of chemical species
    !----------------------------------------

    ! -------------- hno3 -------------------
    ! FORTRAN refactoring: not incorporated in MAM4

    ! -------------- nh3 -------------------
    ! FORTRAN refactoring: not incorporated in MAM4

    ! -------------- so2 -------------------
    ! previous code
    !    heso2(i,k)  = xk*(1._r8 + wrk*(1._r8 + x2/xph(i,k)))
    !    px = heso2(i,k) * Ra * tz * xl
    !    so2g =  xso2(i,k)/(1._r8+ px)
    !    Eso2 = xk*xe*so2g *patm
    ! equivalent new code
    !    heso2 = xk + xk*xe/hplus * xk*xe*x2/hplus**2
    !    so2g = xso2/(1 + px)
    !         = xso2/(1 + heso2*ra*tz*xl)
    !         = xso2/(1 + xk*ra*tz*xl*(1 + (xe/hplus)*(1 + x2/hplus))
    !    eso2 = so2g*xk*xe*patm
    !          = xk*xe*patm*xso2/(1 + xk*ra*tz*xl*(1 + (xe/hplus)*(1 + x2/hplus))
    !          = ( fact1_so2    )/(1 + fact2_so2 *(1 + (fact3_so2/hplus)*(1 + fact4_so2/hplus)
    !    [hso3-] + 2*[so3--] = (eso2/hplus)*(1 + 2*x2/hplus)
    xk = 1.23_r8  *EXP( 3120._r8*t_factor )
    xe = 1.7e-2_r8*EXP( 2090._r8*t_factor )
    x2 = 6.0e-8_r8*EXP( 1120._r8*t_factor )
    fact1_so2 = xk*xe*patm*xso2
    fact2_so2 = xk*Ra*temperature*xlwc
    fact3_so2 = xe
    fact4_so2 = x2

    ! -------------- h2o effects -------------------
    Eh2o = xkw

    ! -------------- co2 effects -------------------
    xk = 3.1e-2_r8*EXP( 2423._r8*t_factor )
    xe = 4.3e-7_r8*EXP(-913._r8 *t_factor )
    Eco2 = xk*xe*co2g  *patm

    ! -------------- so4 effects -------------------
    Eso4 = xso4*xhnm   &         ! /cm3(a)
               *const0/xlwc

    !-----------------------------------------------------------------
    ! now use bisection method to solve electro-neutrality equation
    ! to calculate PH value and H+ concentration
    !
    ! during the iteration loop,
    !    yph_lo = lower ph value that brackets the root (i.e., correct ph)
    !    yph_hi = upper ph value that brackets the root (i.e., correct ph)
    !    yph    = current ph value
    !    yposnet_lo and yposnet_hi = net positive ions for
    !       yph_lo and yph_hi
    !-----------------------------------------------------------------

    converged = .false.
    ! ---------  1st iteration: set lower bound ph value ----------
    yph_lo = 2.0_r8
    yph_hi = yph_lo
    yph = yph_lo
    call calc_ynetpos (          yph,                   & ! in
         fact1_so2,   fact2_so2, fact3_so2, fact4_so2,  & ! in
         Eco2,  Eh2o, Eso4,      so4_fact,              & ! in
         xph,   ynetpos                                 ) ! out
    if (ynetpos <= 0.0_r8) then
    ! the lower and upper bound ph values (2.0 and 7.0) do not bracket
    !    the correct ph, so use the lower bound
          converged = .true.
          return
    endif
    ynetpos_lo = ynetpos

    ! ---------  2nd iteration: set upper bound ph value ----------
    yph_hi = 7.0_r8
    yph = yph_hi
    call calc_ynetpos (          yph,                   & ! in
         fact1_so2,   fact2_so2, fact3_so2, fact4_so2,  & ! in
         Eco2,  Eh2o, Eso4,      so4_fact,              & ! in
         xph,   ynetpos                                 ) ! out
    if (ynetpos >= 0.0_r8) then
    ! the lower and upper bound ph values (2.0 and 7.0) do not bracket
    !    the correct ph, so use the lower bound
          converged = .true.
          return
    endif
    ynetpos_hi = ynetpos

    ! --------- 3rd iteration and more ------------
    do iter = 3, itermax
        yph = 0.5_r8*(yph_lo + yph_hi)
        call calc_ynetpos (          yph,                   & ! in
             fact1_so2,   fact2_so2, fact3_so2, fact4_so2,  & ! in
             Eco2,  Eh2o, Eso4,      so4_fact,              & ! in
             xph,   ynetpos                                 ) ! out
        if (ynetpos >= 0.0_r8) then
           ! net positive ions are >= 0 for both yph and yph_lo
           !    so replace yph_lo with yph
           yph_lo = yph
           ynetpos_lo = ynetpos
        else
           ! net positive ions are <= 0 for both yph and yph_hi
           !    so replace yph_hi with yph
           yph_hi = yph
           ynetpos_hi = ynetpos
        endif

        if (abs(yph_hi - yph_lo) .le. 0.005_r8) then
           ! |yph_hi - yph_lo| <= convergence criterion, so set
           !    final ph to their midpoint and exit
           ! (.005 absolute error in pH gives .01 relative error in H+)
           yph = 0.5_r8*(yph_hi + yph_lo)
           xph = 10.0_r8**(-yph)
           converged = .true.
           return
        endif
    enddo

  end subroutine calc_ph_values


!===========================================================================
  subroutine calc_ynetpos(   yph,                               & ! in
                fact1_so2, fact2_so2, fact3_so2, fact4_so2,     & ! in
                Eco2,      Eh2o,      Eso4,      so4_fact,      & ! in
                xph,       ynetpos                              ) ! out
    !-----------------------------------------------------------------
    ! calculate net positive ions (ynetpos) for iterations in calc_ph_values
    ! also calculate H+ concentration (xph) from ph value
    !-----------------------------------------------------------------
    implicit none

    real(r8), intent(in) :: yph         ! pH value
    real(r8), intent(in) :: fact1_so2, fact2_so2, fact3_so2, fact4_so2 ! factors for SO2
    real(r8), intent(in) :: Eco2, Eh2o, Eso4, so4_fact          ! effects from species [1/cm3]

    real(r8), intent(out) :: xph        ! H+ concentration from pH value [kg/L]
    real(r8), intent(out) :: ynetpos    ! net positive ions

    ! local variables
    real(r8) :: Eso2    ! effect of so2, which is related to pH value
    real(r8) :: tmp_hso3, tmp_so3, tmp_hco3, tmp_oh, tmp_so4  ! temporary variables
    real(r8) :: tmp_pos, tmp_neg        ! positive and negative values to calculate ynetpos

    ! calc current [H+] from ph
    xph = 10.0_r8**(-yph)

    !-----------------------------------------------------------------
    !          ... so2
    !-----------------------------------------------------------------
    Eso2 = fact1_so2/(1.0_r8 + fact2_so2*(1.0_r8 +(fact3_so2/xph) &
                    *(1.0_r8 + fact4_so2/xph)))

    tmp_hso3 = Eso2 / xph
    tmp_so3  = tmp_hso3 * 2.0_r8*fact4_so2/xph
    tmp_hco3 = Eco2 / xph
    tmp_oh   = Eh2o / xph
    tmp_so4 = so4_fact*Eso4

    ! positive ions are H+ only
    tmp_pos = xph
    ! all negative ions
    tmp_neg = tmp_oh + tmp_hco3 + tmp_hso3 + tmp_so3 + tmp_so4

    ynetpos = tmp_pos - tmp_neg

  end subroutine calc_ynetpos

!===========================================================================
  subroutine calc_sox_aqueous( modal_aerosols,         &
                rah2o2, h2o2g, so2g, o3g,      rao3,   &
                patm, dtime, t_factor, xlwc, const0, &
                xhnm, heo3, heso2,      &
                xso2, xso4, xso4_init, xh2o2, &
                xdelso4hp_ik)
    !-----------------------------------------------------------------
    !       ... Prediction after aqueous phase
    !       so4
    !       When Cloud is present
    !
    !       S(IV) + H2O2 = S(VI)
    !       S(IV) + O3   = S(VI)
    !
    !       reference:
    !           (1) Seinfeld
    !           (2) Benkovitz
    !-----------------------------------------------------------------
    implicit none

    logical,  intent(in) :: modal_aerosols      ! if using MAM
    real(r8), intent(in) :: rah2o2      ! reaction rate with h2o2
    real(r8), intent(in) :: rao3        ! reaction rate with o3
    real(r8), intent(in) :: h2o2g, so2g, o3g    
    real(r8), intent(in) :: patm        ! pressure [atm]
    real(r8), intent(in) :: dtime       ! time step [s]
    real(r8), intent(in) :: t_factor    ! working variables to convert temperature 
    real(r8), intent(in) :: xlwc
    real(r8), intent(in) :: const0
    real(r8), intent(in) :: xhnm
    real(r8), intent(in) :: heo3, heso2 ! henry law constant
    real(r8), intent(inout) :: xso2, xso4, xso4_init, xh2o2 ! mixing ratios
    real(r8), intent(out) :: xdelso4hp_ik ! change of so4 in (i,k)

    ! local variables
    real(r8) :: pso4    ! production rate of so4
    real(r8) :: delta_s ! so4 production in the time step

          !............................
          !       S(IV) + H2O2 = S(VI)
          !............................

    pso4 = rah2o2 * 7.4e4_r8*exp(6621._r8*t_factor) * h2o2g * patm &
                  * 1.23_r8 *exp(3120._r8*t_factor) * so2g * patm

    pso4 = pso4   & ! [M/s] = [mole/L(w)/s]
         * xlwc   & ! [mole/L(a)/s]
         / const0 & ! [/L(a)/s]
         / xhnm


    delta_s = max(pso4*dtime, 1.e-30_r8)

    xso4_init=xso4

    if (delta_s<=xso2 .and. delta_s<=xh2o2) then
        xso4  = xso4  + delta_s
        xh2o2 = xh2o2 - delta_s
        xso2  = xso2  - delta_s
    elseif (xh2o2 > xso2) then
        xso4=xso4+xso2
        xh2o2=xh2o2-xso2
        xso2=1.e-20_r8
    else
        xso4=xso4+xh2o2
        xso2=xso2-xh2o2
        xh2o2=1.e-20_r8
    endif

    if (modal_aerosols) then
       xdelso4hp_ik  =  xso4 - xso4_init
    endif
             !...........................
             !       S(IV) + O3 = S(VI)
             !...........................

    pso4 = rao3 * heo3*o3g*patm * heso2*so2g*patm  ! [M/s]

    pso4 = pso4        &                                ! [M/s] =[mole/L(w)/s]
         * xlwc        &                                ! [mole/L(a)/s]
         / const0      &                                ! [/L(a)/s]
         / xhnm                                    ! [mixing ratio/s]

    delta_s = max(pso4*dtime, 1.e-30_r8)

    xso4_init=xso4

    if (delta_s > xso2) then
       xso4 = xso4 + xso2
       xso2 = 1.e-20_r8
    else
       xso4 = xso4 + delta_s
       xso2 = xso2 - delta_s
    endif


  end subroutine calc_sox_aqueous

!===========================================================================

end module mo_setsox
