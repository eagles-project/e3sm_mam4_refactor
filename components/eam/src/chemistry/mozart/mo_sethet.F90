
module mo_sethet

!
! LKE (10/11/2010): added HCN, CH3CN, HCOOH  to cesm1_0_beta07_offline version
!                   HCN, CH3CN have new Henry's Law coefficients, HCOOH is set to CH3COOH
! LKE (10/18/2010): SO2 washout corrected based on recommendation of R.Easter, PNNL
!
  use shr_kind_mod,    only: r8 => shr_kind_r8
  use cam_logfile,     only: iulog
  use gas_wetdep_opts, only: gas_wetdep_cnt, gas_wetdep_method, gas_wetdep_list
  use phys_control,    only: phys_getopts

  private
  public :: sethet_inti, sethet

  save

  integer :: h2o2_ndx, hno3_ndx, ch2o_ndx, ch3ooh_ndx, ch3coooh_ndx, &
       ho2no2_ndx, ch3cocho_ndx, xooh_ndx, onitr_ndx, glyald_ndx, &
       ch3cho_ndx, mvk_ndx, macr_ndx, pooh_ndx, c2h5ooh_ndx, &
       c3h7ooh_ndx, rooh_ndx, isopno3_ndx, onit_ndx, Pb_ndx, &
       macrooh_ndx, isopooh_ndx, ch3oh_ndx, c2h5oh_ndx, hyac_ndx, hydrald_ndx
  integer :: spc_h2o2_ndx, spc_hno3_ndx
  integer :: spc_so2_ndx
  integer :: spc_sogm_ndx, spc_sogi_ndx, spc_sogt_ndx, spc_sogb_ndx, spc_sogx_ndx

  integer :: alkooh_ndx, mekooh_ndx, tolooh_ndx, terpooh_ndx, ch3cooh_ndx
  integer :: so2_ndx, soa_ndx, so4_ndx, cb2_ndx, oc2_ndx, nh3_ndx, nh4no3_ndx, &
             sa1_ndx, sa2_ndx, sa3_ndx, sa4_ndx, nh4_ndx, h2so4_ndx
  integer :: xisopno3_ndx,xho2no2_ndx,xonitr_ndx,xhno3_ndx,xonit_ndx
  integer :: clono2_ndx, brono2_ndx, hcl_ndx, n2o5_ndx, hocl_ndx, hobr_ndx, hbr_ndx 
  integer :: ch3cn_ndx, hcn_ndx, hcooh_ndx
  integer, allocatable :: wetdep_map(:)
  integer :: sogm_ndx, sogi_ndx, sogt_ndx, sogb_ndx, sogx_ndx
  logical :: do_wetdep

  ! prognostic modal aerosols
  logical :: prog_modal_aero

contains

!=================================================================================
  subroutine sethet_inti
    !-----------------------------------------------------------------------      
    !       ... intialize the wet removal rate constants routine
    !-----------------------------------------------------------------------      

    use mo_chem_utls, only : get_het_ndx, get_spc_ndx
    use spmd_utils,   only : masterproc
    use cam_abortutils,   only : endrun

    integer :: k, m
    
    do_wetdep = gas_wetdep_cnt>0 .and. gas_wetdep_method=='MOZ'
    if ( .not. do_wetdep) return

    call phys_getopts( prog_modal_aero_out = prog_modal_aero )

    allocate( wetdep_map(gas_wetdep_cnt))

    do k=1,gas_wetdep_cnt
       m = get_het_ndx( trim(gas_wetdep_list(k))) 
       if (m>0) then
          wetdep_map(k) = m
       else
          call endrun('sethet_inti: cannot map '//trim(gas_wetdep_list(k)))
       endif
    enddo

    xisopno3_ndx = get_het_ndx( 'XISOPNO3' )
    xho2no2_ndx  = get_het_ndx( 'XHO2NO2' )
    xonitr_ndx   = get_het_ndx( 'XONITR' )
    xhno3_ndx    = get_het_ndx( 'XHNO3' )
    xonit_ndx    = get_het_ndx( 'XONIT' )

    spc_h2o2_ndx = get_spc_ndx( 'H2O2' )
    spc_hno3_ndx = get_spc_ndx( 'HNO3' )
    spc_so2_ndx  = get_spc_ndx( 'SO2' )

    clono2_ndx = get_het_ndx( 'CLONO2' )
    brono2_ndx = get_het_ndx( 'BRONO2' )
    hcl_ndx    = get_het_ndx( 'HCL' )
    n2o5_ndx   = get_het_ndx( 'N2O5' )
    hocl_ndx   = get_het_ndx( 'HOCL' )
    hobr_ndx   = get_het_ndx( 'HOBR' )
    hbr_ndx    = get_het_ndx( 'HBR' )

    h2o2_ndx   = get_het_ndx( 'H2O2' )
    hno3_ndx   = get_het_ndx( 'HNO3' )
    ch2o_ndx   = get_het_ndx( 'CH2O' )
    ch3ooh_ndx = get_het_ndx( 'CH3OOH' )
    ch3coooh_ndx = get_het_ndx( 'CH3COOOH' )
    ho2no2_ndx  = get_het_ndx( 'HO2NO2' )
    ch3cocho_ndx = get_het_ndx( 'CH3COCHO' )
    xooh_ndx    = get_het_ndx( 'XOOH' )
    onitr_ndx   = get_het_ndx( 'ONITR' )
    glyald_ndx  = get_het_ndx( 'GLYALD' )
    ch3cho_ndx  = get_het_ndx( 'CH3CHO' )
    mvk_ndx     = get_het_ndx( 'MVK' )
    macr_ndx    = get_het_ndx( 'MACR' )
    pooh_ndx    = get_het_ndx( 'POOH' )
    c2h5ooh_ndx = get_het_ndx( 'C2H5OOH' )
    c3h7ooh_ndx = get_het_ndx( 'C3H7OOH' )
    rooh_ndx    = get_het_ndx( 'ROOH' )
    isopno3_ndx = get_het_ndx( 'ISOPNO3' )
    onit_ndx    = get_het_ndx( 'ONIT' )
    Pb_ndx      = get_het_ndx( 'Pb' )
    macrooh_ndx = get_het_ndx( 'MACROOH' )
    isopooh_ndx = get_het_ndx( 'ISOPOOH' )
    ch3oh_ndx   = get_het_ndx( 'CH3OH' )
    c2h5oh_ndx  = get_het_ndx( 'C2H5OH' )
    hyac_ndx    = get_het_ndx( 'HYAC' )
    hydrald_ndx = get_het_ndx( 'HYDRALD' )
    alkooh_ndx  = get_het_ndx( 'ALKOOH' )
    mekooh_ndx  = get_het_ndx( 'MEKOOH' )
    tolooh_ndx  = get_het_ndx( 'TOLOOH' )
    terpooh_ndx = get_het_ndx( 'TERPOOH' )
    ch3cooh_ndx = get_het_ndx( 'CH3COOH' )
    so2_ndx     = get_het_ndx( 'SO2' )
    soa_ndx     = get_het_ndx( 'SOA' )
    sogb_ndx    = get_het_ndx( 'SOGB' )
    sogi_ndx    = get_het_ndx( 'SOGI' )
    sogm_ndx    = get_het_ndx( 'SOGM' )
    sogt_ndx    = get_het_ndx( 'SOGT' )
    sogx_ndx    = get_het_ndx( 'SOGX' )
    so4_ndx     = get_het_ndx( 'SO4' )
    cb2_ndx     = get_het_ndx( 'CB2' )
    oc2_ndx     = get_het_ndx( 'OC2' )
    nh3_ndx     = get_het_ndx( 'NH3' )
    nh4no3_ndx  = get_het_ndx( 'NH4NO3' )
    nh4_ndx     = get_het_ndx( 'NH4' )
    h2so4_ndx   = get_het_ndx( 'H2SO4' )
    sa1_ndx     = get_het_ndx( 'SA1' )
    sa2_ndx     = get_het_ndx( 'SA2' )
    sa3_ndx     = get_het_ndx( 'SA3' )
    sa4_ndx     = get_het_ndx( 'SA4' )
    ch3cn_ndx   = get_het_ndx( 'CH3CN' )
    hcn_ndx     = get_het_ndx( 'HCN' )
    hcooh_ndx   = get_het_ndx( 'HCOOH' )

    if (masterproc) then
       write(iulog,*) 'sethet_inti: new ndx ',so2_ndx,soa_ndx,so4_ndx,cb2_ndx,oc2_ndx, &
            nh3_ndx,nh4no3_ndx,sa1_ndx,sa2_ndx,sa3_ndx,sa4_ndx
       write(iulog,*) ' '
       write(iulog,*) 'sethet_inti: diagnotics '
       write(iulog,'(10i5)') h2o2_ndx, hno3_ndx, ch2o_ndx, ch3ooh_ndx, ch3coooh_ndx, &
            ho2no2_ndx, ch3cocho_ndx, xooh_ndx, onitr_ndx, glyald_ndx, &
            ch3cho_ndx, mvk_ndx, macr_ndx, pooh_ndx, c2h5ooh_ndx, &
            c3h7ooh_ndx, rooh_ndx, isopno3_ndx, onit_ndx, Pb_ndx, &
            macrooh_ndx, isopooh_ndx, ch3oh_ndx, c2h5oh_ndx, hyac_ndx, hydrald_ndx
    endif

  end subroutine sethet_inti

!=================================================================================
  subroutine sethet( het_rates, press, zmid,  phis, tfld, &
                     cmfdqr, nrain, nevapr, delt, xhnm, &
                     qin, ncol, lchnk )
    !-----------------------------------------------------------------------      
    !       ... compute rainout loss rates (1/s)
    !-----------------------------------------------------------------------      

    use physconst,    only : rga,pi
    use chem_mods,    only : gas_pcnst
    use ppgrid,       only : pver, pcols
    use phys_grid,    only : get_rlat_all_p
    use cam_abortutils,   only : endrun
    use mo_constants, only : avo => avogadro, boltz_cgs

    implicit none
    !-----------------------------------------------------------------------      
    !       ... dummy arguments
    !-----------------------------------------------------------------------      
    integer, intent(in)   ::    ncol                        ! columns in chunk
    integer, intent(in)   ::    lchnk                       ! chunk index
    real(r8), intent(in)  ::    delt                        ! time step ( s )
    real(r8), intent(in)  ::    press(pcols,pver)           ! pressure in pascals
    real(r8), intent(in)  ::    cmfdqr(ncol,pver)           ! dq/dt for convection
    real(r8), intent(in)  ::    nrain(ncol,pver)            ! stratoform precip
    real(r8), intent(in)  ::    nevapr(ncol,pver)           ! evaporation
    real(r8), intent(in)  ::    qin(ncol,pver,gas_pcnst)    ! xported species ( vmr )
    real(r8), intent(in)  ::    zmid(ncol,pver)             ! midpoint geopot (km)
    real(r8), intent(in)  ::    phis(pcols)                 ! surf geopot
    real(r8), intent(in)  ::    tfld(pcols,pver)            ! temperature (k)
    real(r8), intent(in)  ::    xhnm(ncol,pver)             ! total atms density ( /cm^3)
    real(r8), intent(out) ::    het_rates(ncol,pver,gas_pcnst) ! rainout loss rates

    !-----------------------------------------------------------------------      
    !       ... local variables
    !-----------------------------------------------------------------------      
    real(r8), parameter ::  xrm   = .189_r8             ! mean diameter of rain drop (cm)
    real(r8), parameter ::  xum   = 748._r8             ! mean rain drop terminal velocity (cm/s)
    real(r8), parameter ::  xvv   = 6.18e-2_r8          ! kinetic viscosity (cm^2/s)
    real(r8), parameter ::  xdg   = .112_r8             ! mass transport coefficient (cm/s)
    real(r8), parameter ::  t0    = 298._r8             ! reference temperature (k)
    real(r8), parameter ::  xph0  = 1.e-5_r8            ! cloud [h+]
    real(r8), parameter ::  satf_hno3  = .016_r8        ! saturation factor for hno3 in clouds 
    real(r8), parameter ::  satf_h2o2  = .016_r8        ! saturation factor for h2o2 in clouds 
    real(r8), parameter ::  satf_so2   = .016_r8        ! saturation factor for so2 in clouds 
    real(r8), parameter ::  satf_ch2o  = .1_r8          ! saturation factor for ch2o in clouds 
    real(r8), parameter ::  satf_sog  =  .016_r8        ! saturation factor for sog in clouds
    real(r8), parameter ::  const0   = boltz_cgs * 1.e-6_r8 ! (atmospheres/deg k/cm^3)
    real(r8), parameter ::  hno3_diss = 15.4_r8         ! hno3 dissociation constant
    real(r8), parameter ::  geo_fac  = 6._r8            ! geometry factor (surf area/volume = geo_fac/diameter)
    real(r8), parameter ::  mass_air = 29._r8           ! mass of background atmosphere (amu)
    real(r8), parameter ::  mass_h2o = 18._r8           ! mass of water vapor (amu)
    real(r8), parameter ::  h2o_mol  = 1.e3_r8/mass_h2o ! (gm/mol water)
    real(r8), parameter ::  km2cm    = 1.e5_r8          ! convert km to cm
    real(r8), parameter ::  m2km     = 1.e-3_r8         ! convert m to km
    real(r8), parameter ::  cm3_2_m3 = 1.e-6_r8         ! convert cm^3 to m^3
    real(r8), parameter ::  m3_2_cm3 = 1.e6_r8          ! convert m^3 to cm^3
    real(r8), parameter ::  liter_per_gram = 1.e-3_r8
    real(r8), parameter ::  avo2  = avo * liter_per_gram * cm3_2_m3 ! (liter/gm/mol*(m/cm)^3)

    integer  ::      i, m, k, kk                 ! indicies
    real(r8) ::      xkgm                        ! mass flux on rain drop
    real(r8) ::      all1, all2                  ! work variables
    real(r8) ::      stay                        ! fraction of layer traversed by falling drop in timestep delt
    real(r8) ::      xeqca1, xeqca2, xca1, xca2, xdtm
    real(r8) ::      xxx1, xxx2, yhno3, yh2o2
    real(r8) ::      all3, xeqca3, xca3, xxx3, yso2, so2_diss(ncol)
    real(r8) ::      all4, xeqca4, xca4, xxx4
    real(r8) ::      all5, xeqca5, xca5, xxx5
    real(r8) ::      all6, xeqca6, xca6, xxx6
    real(r8) ::      all7, xeqca7, xca7, xxx7
    real(r8) ::      all8, xeqca8, xca8, xxx8
    real(r8) ::      ysogm,ysogi,ysogt,ysogb,ysogx

    real(r8) :: t_factor(ncol)   ! temperature factor to calculate henry's law parameters
    real(r8), dimension(ncol)  :: &
         xk0, work1, work2, work3, zsurf
    real(r8), dimension(pver)  :: &
         xgas1, xgas2
    real(r8), dimension(pver)  :: xgas3, xgas4, xgas5, xgas6, xgas7, xgas8
    real(r8), dimension(ncol)  :: &
         tmp0_rates, tmp1_rates
    real(r8), dimension(ncol,pver)  :: &
         delz, &              ! layer depth about interfaces (cm)
         xhno3, &             ! hno3 concentration (molecules/cm^3)
         xh2o2, &             ! h2o2 concentration (molecules/cm^3)
         xso2, &              ! so2 concentration (molecules/cm^3)
         xsogm, &             ! sogm concentration (molecules/cm^3)
         xsogi, &             ! sogi concentration (molecules/cm^3)
         xsogt, &             ! sogt concentration (molecules/cm^3)
         xsogb, &             ! sogb concentration (molecules/cm^3)
         xsogx, &             ! sogx concentration (molecules/cm^3)
         xliq, &              ! liquid rain water content in a grid cell (gm/m^3)
         rain                 ! conversion rate of water vapor into rain water (molecules/cm^3/s)
    real(r8), dimension(ncol,pver)  :: &
         xhen_hno3, xhen_h2o2, xhen_ch2o, xhen_ch3ooh, xhen_ch3co3h, &
         xhen_ch3cocho, xhen_xooh, xhen_onitr, xhen_ho2no2, xhen_glyald, &
         xhen_ch3cho, xhen_mvk, xhen_macr,xhen_sog
    real(r8), dimension(ncol,pver)  :: &
         xhen_nh3, xhen_ch3cooh
    real(r8), dimension(ncol,pver,8) :: tmp_hetrates
    real(r8), dimension(ncol,pver)  :: precip
    real(r8), dimension(ncol,pver)  :: xhen_hcn, xhen_ch3cn, xhen_so2

    integer    ::      ktop_all       
    integer    ::      ktop(ncol)                  ! 100 mb level

    real(r8) :: rlat(pcols)                       ! latitude in radians for columns
    real(r8) :: p_limit
    real(r8), parameter :: d2r = pi/180._r8
!
! jfl : new variables for rescaling sum of positive values to actual amount
!
    real(r8) :: total_rain,total_pos
    character(len=3) :: hetratestrg
    real(r8), parameter :: MISSING = -999999._r8
    integer ::  mm

!
    !-----------------------------------------------------------------
    !        note: the press array is in pascals and must be
    !              mutiplied by 10 to yield dynes/cm**2.
    !-----------------------------------------------------------------
    !       ... set wet deposition for
    !           1. h2o2         2. hno3
    !           3. ch2o         4. ch3ooh
    !           5. pooh         6. ch3coooh
    !           7. ho2no2       8. onit
    !           9. mvk         10. macr
    !          11. c2h5ooh     12. c3h7ooh
    !          13. rooh        14. ch3cocho
    !          15. pb          16. macrooh
    !          17. xooh        18. onitr
    !          19. isopooh     20. ch3oh
    !          21. c2h5oh      22. glyald
    !          23. hyac        24. hydrald
    !          25. ch3cho      26. isopno3
    !-----------------------------------------------------------------
    ! FORTRAN refactor note: current MAM4 only have three species in default:
    ! 'H2O2','H2SO4','SO2'.  Options for other species are then removed
    !-----------------------------------------------------------------


    het_rates(:,:,:) = 0._r8

    if ( .not. do_wetdep) return

    call get_rlat_all_p(lchnk, ncol, rlat)

    do mm = 1,gas_wetdep_cnt
       m = wetdep_map(mm)
       if ( m>0 ) then
          het_rates(:,:,m) = MISSING
       endif
    end do

    !-----------------------------------------------------------------
    !	... the 2 and .6 multipliers are from a formula by frossling (1938)
    !-----------------------------------------------------------------
    xkgm = xdg/xrm * 2._r8 + xdg/xrm * .6_r8 * sqrt( xrm*xum/xvv ) * (xvv/xdg)**(1._r8/3._r8) 

    !-----------------------------------------------------------------
    !	... Find the level index that only calculate het_rates below
    !-----------------------------------------------------------------
    call find_ktop( ncol,  rlat,  press,  & ! in
                    ktop                  ) ! out
    ktop_all = minval( ktop(:) )

    ! this is added to rescale the variable precip (which can only be positive)
    ! to the actual vertical integral of positive and negative values.  This
    ! removes point storms
    call calc_precip_rescale( ncol, cmfdqr, nrain, nevapr,  & ! in
                              precip                        ) ! out

    do k = 1,pver
       rain(:ncol,k)   = mass_air*precip(:ncol,k)*xhnm(:ncol,k) / mass_h2o
       xliq(:ncol,k)   = precip(:ncol,k) * delt * xhnm(:ncol,k) / avo*mass_air * m3_2_cm3
       xh2o2(:ncol,k)  = qin(:ncol,k,spc_h2o2_ndx) * xhnm(:ncol,k)
       xso2(:ncol,k)  = qin(:ncol,k,spc_so2_ndx) * xhnm(:ncol,k)
    end do

    zsurf(:ncol) = m2km * phis(:ncol) * rga
    do k = ktop_all,pver-1
       delz(:ncol,k) = abs( (zmid(:ncol,k) - zmid(:ncol,k+1))*km2cm ) 
    end do
    delz(:ncol,pver) = abs( (zmid(:ncol,pver) - zsurf(:ncol) )*km2cm ) 

    !-----------------------------------------------------------------
    !       ... part 0b,  for temperature dependent of henrys
    !                     xxhe1 = henry con for hno3
    !                     xxhe2 = henry con for h2o2
    !lwh 10/00 -- take henry''s law constants from brasseur et al. [1999],
    !             appendix j. for hno3, also consider dissociation to
    !             get effective henry''s law constant; equilibrium
    !             constant for dissociation from brasseur et al. [1999],
    !             appendix k. assume ph=5 (set as xph0 above).
    !             heff = h*k/[h+] for hno3 (complete dissociation)
    !             heff = h for h2o2 (no dissociation)
    !             heff = h * (1 + k/[h+]) (in general)
    !-----------------------------------------------------------------
    do k = ktop_all,pver
       !-----------------------------------------------------------------
       ! 	... effective henry''s law constants:
       !	hno3, h2o2  (brasseur et al., 1999)
       !-----------------------------------------------------------------
       ! temperature factor
       t_factor(:ncol) = (t0 - tfld(:ncol,k))/(t0*tfld(:ncol,k))
       xhen_h2o2(:,k)     = 7.45e4_r8 * exp( 6620._r8 * t_factor(:) )
       ! HNO3, for calculation of H2SO4 het rate use
       xk0(:)             = 2.1e5_r8 *exp( 8700._r8*t_factor(:) )
       xhen_hno3(:,k)     = xk0(:) * ( 1._r8 + hno3_diss / xph0 )
       ! SO2
       xk0(:)             = 1.23_r8 * exp( 3120._r8 * t_factor(:) )
       so2_diss(:)        = 1.23e-2_r8 * exp( 1960._r8 * t_factor(:) )
       xhen_so2(:,k)   = xk0(:) * ( 1._r8 + so2_diss(:) / xph0 )

       ! initiate temporary array
       tmp_hetrates(:,k,:) = 0._r8
    enddo

    !-----------------------------------------------------------------
    !       ... part 1, solve for high henry constant ( hno3, h2o2)
    !-----------------------------------------------------------------
    col_loop :  do i = 1,ncol
       xgas2(:) = xh2o2(i,:)                     ! different levels wash 
       xgas3(:) = xso2 (i,:)
       level_loop1  : do kk = ktop(i),pver
          stay = 1._r8
          if( rain(i,kk) /= 0._r8 ) then            ! finding rain cloud           
             stay = ((zmid(i,kk) - zsurf(i))*km2cm)/(xum*delt)
             stay = min( stay,1._r8 )
             ! calculate gas washout by cloud
             call gas_washout( kk,  xkgm,   xliq(i,kk),       & ! in
                        xhen_h2o2(i,:), tfld(i,:), delz(i,:), & ! in
                        xgas2                                 ) ! inout
             call gas_washout( kk,  xkgm,   xliq(i,kk),       & ! in
                         xhen_so2(i,:), tfld(i,:), delz(i,:), & ! in
                        xgas3                                 ) ! inout
          endif
          !-----------------------------------------------------------------
          !       ... calculate the lifetime of washout (second)
          !             after all layers washout 
          !             the concentration of hno3 is reduced 
          !             then the lifetime xtt is calculated by
          !
          !                  xtt = (xhno3(ini) - xgas1(new))/(dt*xhno3(ini))
          !                  where dt = passing time (s) in vertical
          !                             path below the cloud
          !                        dt = dz(cm)/um(cm/s)
          !-----------------------------------------------------------------
          xdtm = delz(i,kk) / xum                     ! the traveling time in each dz
          xxx2 = (xh2o2(i,kk) - xgas2(kk))
          if( xxx2 /= 0._r8 ) then                       ! if no washout lifetime = 1.e29
             yh2o2  = xh2o2(i,kk)/xxx2 * xdtm     
          else
             yh2o2  = 1.e29_r8
          end if
          tmp_hetrates(i,kk,1) = max( 1._r8 / yh2o2,0._r8 ) * stay
          xxx3 = (xso2( i,kk) - xgas3(kk))
          if( xxx3 /= 0._r8 ) then                       ! if no washout lifetime = 1.e29
             yso2  = xso2( i,kk)/xxx3 * xdtm     
          else
             yso2  = 1.e29_r8
          end if
          tmp_hetrates(i,kk,3) = max( 1._r8 / yso2, 0._r8 ) * stay
       end do level_loop1
    end do col_loop

    !-----------------------------------------------------------------
    !       ... part 2, in-cloud solve for low henry constant
    !                   hno3 and h2o2 have both in and under cloud
    !-----------------------------------------------------------------
    level_loop2 : do k = ktop_all,pver
       Column_loop2 : do i=1,ncol
          if ( rain(i,k) <= 0._r8 ) then
             het_rates(i,k,:) =  0._r8 
             cycle
          endif

          work1(i) = avo2 * xliq(i,k)
          work2(i) = const0 * tfld(i,k)

          if( h2o2_ndx > 0 ) then
             work3(i) = satf_h2o2 * max( rain(i,k) / (h2o_mol*(work1(i) + 1._r8/(xhen_h2o2(i,k)*work2(i)))),0._r8 )    
             het_rates(i,k,h2o2_ndx) =  work3(i) + tmp_hetrates(i,k,1)
          end if
          if ( prog_modal_aero .and. so2_ndx>0 .and. h2o2_ndx>0 ) then
             het_rates(i,k,so2_ndx) = het_rates(i,k,h2o2_ndx)
          elseif( so2_ndx > 0 ) then
             work3(i) = satf_so2 * max( rain(i,k) / (h2o_mol*(work1(i) + 1._r8/(xhen_so2( i,k)*work2(i)))),0._r8 )    
             het_rates(i,k,so2_ndx ) =  work3(i) + tmp_hetrates(i,k,3)
          endif
!
!
          work3(i) = tmp_hetrates(i,k,2) + satf_hno3 * &
               max( rain(i,k) / (h2o_mol*(work1(i) + 1._r8/(xhen_hno3(i,k)*work2(i)))),0._r8 )    
          tmp0_rates(i)   = work3(i)

          if( h2so4_ndx > 0 ) then
             het_rates(i,k,h2so4_ndx) = tmp0_rates(i) 
          end if
       end do Column_loop2
    end do level_loop2

    !-----------------------------------------------------------------
    !	... Set rates above tropopause = 0.
    !-----------------------------------------------------------------
    do mm = 1,gas_wetdep_cnt
       m = wetdep_map(mm)
       do i = 1,ncol
          do k = 1,ktop(i)
             het_rates(i,k,m) = 0._r8
          end do
       end do
       if ( any( het_rates(:ncol,:,m) == MISSING) ) then
          write(hetratestrg,'(I3)') m
          call endrun('sethet: het_rates (wet dep) not set for het reaction number : '//hetratestrg)
       endif
    end do

  end subroutine sethet

!=================================================================================
  subroutine find_ktop( ncol,  rlat,  press,  & ! in
                        ktop                  ) ! out 
  !---------------------------------------------------------------------------
  ! -------- find the top level that het_rates are set as 0 above it ---------
  !--------------------------------------------------------------------------- 

    use ppgrid,       only : pver, pcols
    use physconst,    only : pi

    implicit none
    integer,  intent(in) :: ncol
    real(r8), intent(in) :: rlat(pcols)          ! latitude in radians for columns
    real(r8), intent(in) :: press(pcols,pver)    ! pressure [Pa]
    integer, intent(out) :: ktop(ncol)           ! index that only calculate het_rates above this level

    integer  :: icol, kk
    real(r8) :: p_limit     ! pressure limit [Pa]
    real(r8), parameter :: d2r = pi/180._r8   ! degree to radian


    do icol = 1,ncol

       if ( abs(rlat(icol)) > 60._r8*d2r ) then
          p_limit = 300.e2_r8   ! 300hPa for high latitudes
       else
          p_limit = 100.e2_r8   ! 100hPa for low latitudes
       endif

       k_loop: do kk = pver,1,-1
          if( press(icol,kk) < p_limit ) then
             ktop(icol) = kk
             exit k_loop
          endif
       enddo k_loop

    enddo

  end subroutine find_ktop

!=================================================================================
  subroutine calc_precip_rescale( ncol, cmfdqr, nrain, nevapr,  & ! in
                                  precip                        ) ! out
  ! -----------------------------------------------------------------------
  ! calculate precipitation rate at each grid
  ! this is added to rescale the variable precip (which can only be positive)
  ! to the actual vertical integral of positive and negative values. 
  ! This removes point storms
  ! -----------------------------------------------------------------------
    use ppgrid,       only : pver, pcols
    implicit none
    integer,  intent(in) :: ncol
    real(r8), intent(in) :: cmfdqr(ncol,pver)           ! dq/dt for convection [kg/kg/s]
    real(r8), intent(in) :: nrain(ncol,pver)            ! stratoform precip [kg/kg/s]
    real(r8), intent(in) :: nevapr(ncol,pver)           ! evaporation [kg/kg/s]
    real(r8),intent(out) :: precip(ncol,pver)           ! precipitation [kg/kg/s]

    integer  :: icol, kk
    real(r8) :: total_rain      ! total rain rate (both pos and neg) in the column
    real(r8) :: total_pos       ! total positive rain rate in the column

    do icol = 1,ncol

       total_rain = 0._r8
       total_pos  = 0._r8
       do kk = 1,pver
          precip(icol,kk) = cmfdqr(icol,kk) + nrain(icol,kk) - nevapr(icol,kk)
          total_rain = total_rain + precip(icol,kk)
          if ( precip(icol,kk) < 0._r8 ) then
                 precip(icol,kk) = 0._r8
          endif
          total_pos  = total_pos  + precip(icol,kk)
       enddo

       if ( total_rain <= 0._r8 ) then
          precip(icol,:) = 0._r8        ! set all levels to zero
       else
          do kk = 1,pver
             precip(icol,kk) = precip(icol,kk) * total_rain/total_pos
          enddo
       endif
    enddo

  end subroutine calc_precip_rescale

!=================================================================================
  subroutine gas_washout ( plev,  xkgm,   xliq_ik,      & ! in
                           xhen_i, tfld_i, delz_i,      & ! in
                           xgas                         ) ! inout
   !------------------------------------------------------------------------
   ! calculate gas washout by cloud if not saturated
   !------------------------------------------------------------------------
    use ppgrid,       only : pver
    use mo_constants, only : avo => avogadro, boltz_cgs

    implicit none
    integer,  intent(in) :: plev   ! calculate from this level below
    real(r8), intent(in) :: xliq_ik ! liquid rain water content [gm/m^3]
    real(r8), intent(in) :: xhen_i(pver) ! henry's law constant 
    real(r8), intent(in) :: tfld_i(pver) ! temperature [K]
    real(r8), intent(in) :: delz_i(pver) ! layer depth about interfaces [cm]
    real(r8), intent(in) :: xkgm
    real(r8), intent(inout) :: xgas(pver)   ! gas concentration

    integer  :: kk
    real(r8) :: allca   ! total of ca between level plev and kk [#/cm3]
    real(r8) :: xca, xeqca
    real(r8), parameter ::  const0   = boltz_cgs * 1.e-6_r8 ! [atmospheres/deg k/cm^3]
    real(r8), parameter ::  geo_fac  = 6._r8            ! geometry factor (surf area/volume = geo_fac/diameter)
    real(r8), parameter ::  xrm   = .189_r8             ! mean diameter of rain drop [cm]
    real(r8), parameter ::  xum   = 748._r8             ! mean rain drop terminal velocity [cm/s]
    real(r8), parameter ::  cm3_2_m3 = 1.e-6_r8         ! convert cm^3 to m^3
    real(r8), parameter ::  liter_per_gram = 1.e-3_r8
    real(r8), parameter ::  avo2  = avo * liter_per_gram * cm3_2_m3 ! [L/gm/mol*(m/cm)^3]

     allca = 0._r8
     !-----------------------------------------------------------------
     !       ... calculate the saturation concentration eqca
     !-----------------------------------------------------------------
     do kk = plev,pver                      ! cal washout below cloud
        xeqca =  xgas(kk) &
               / (xliq_ik*avo2 + 1._r8/(xhen_i(kk)*const0*tfld_i(kk))) &
               *  xliq_ik*avo2

        !-----------------------------------------------------------------
        !       ... calculate ca; inside cloud concentration in  #/cm3(air)
        !-----------------------------------------------------------------
        xca = geo_fac*xkgm*xgas(kk)/(xrm*xum)*delz_i(kk) * xliq_ik * cm3_2_m3

        !-----------------------------------------------------------------
        !       ... if is not saturated (take hno3 as an example)
        !               hno3(gas)_new = hno3(gas)_old - hno3(h2o)
        !           otherwise
        !               hno3(gas)_new = hno3(gas)_old
        !-----------------------------------------------------------------
        allca = allca + xca
        if( allca < xeqca ) then
           xgas(kk) = max( xgas(kk) - xca, 0._r8 )
        endif
     enddo

  end subroutine gas_washout

!=================================================================================
end module mo_sethet
