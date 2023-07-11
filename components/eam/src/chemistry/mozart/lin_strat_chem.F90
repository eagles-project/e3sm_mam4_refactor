
!--------------------------------------------------------------------
! linearized ozone chemistry LINOZ
! from Hsu and Prather, JGR, 2008
!
! written by Jean-Francois Lamarque (September 2008)
! modified by
!     24 Oct 2008 -- Francis Vitt
!      9 Dec 2008 -- Philip Cameron-Smith, LLNL, -- added ltrop
!      4 Jul 2019 -- Qi Tang (LLNL), Juno Hsu (UCI), -- added sfcsink
!--------------------------------------------------------------------
module lin_strat_chem

  use shr_kind_mod, only : r8 => shr_kind_r8
  use physics_types,only : physics_state
  use cam_logfile,  only : iulog
  use cam_abortutils,   only : endrun
  use spmd_utils,   only : masterproc
  use physconst,     only : pi
  use chlorine_loading_data, only: chlorine_loading
  !
  implicit none
  !
  private  ! all unless made public

  save
  !
  ! define public components of module
  !
  public :: lin_strat_chem_inti, lin_strat_chem_solve, lin_strat_sfcsink
  public :: do_lin_strat_chem
  public :: linoz_readnl   ! read linoz_nl namelist

  integer :: index_o3

  logical :: do_lin_strat_chem

  real(r8), parameter :: chlorine_loading_bgnd    = 0.0000_r8     ! EESC value [ppbv] for background conditions
  real(r8), parameter :: radians_to_degrees       = 180._r8/pi

  real(r8), parameter :: unset_r8   = huge(1.0_r8)
  integer , parameter :: unset_int  = huge(1)
  integer  :: linoz_lbl = unset_int ! number of layers with ozone decay from the surface
  real(r8) :: linoz_sfc = unset_r8  ! boundary layer concentration (ppb) to which Linoz ozone e-fold
  real(r8) :: linoz_tau = unset_r8  ! Linoz e-fold time scale (in seconds) in the boundary layer
  real(r8) :: linoz_psc_T = unset_r8  ! PSC ozone loss T (K) threshold

  integer  :: o3_lbl ! set from namelist input linoz_lbl
  real(r8) :: o3_sfc ! set from namelist input linoz_sfc
  real(r8) :: o3_tau ! set from namelist input linoz_tau
  real(r8) :: psc_T  ! set from namelist input linoz_psc_T

contains

subroutine linoz_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'linoz_readnl'

   namelist /linoz_nl/ linoz_lbl, linoz_sfc, linoz_tau, linoz_psc_T
   !-----------------------------------------------------------------------------

   ! Set default values
   linoz_lbl = 4
   linoz_sfc = 30.0e-9_r8
   linoz_tau = 172800.0_r8
   linoz_psc_T = 193.0_r8

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'linoz_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, linoz_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)

      ! set local variables
      o3_lbl = linoz_lbl
      o3_sfc = linoz_sfc
      o3_tau = linoz_tau
      psc_T  = linoz_psc_T

      ! check
      write(iulog,*) subname // ', linoz_lbl:', o3_lbl
      write(iulog,*) subname // ', linoz_sfc:', o3_sfc
      write(iulog,*) subname // ', linoz_tau:', o3_tau
      write(iulog,*) subname // ', linoz_psc_T:', psc_T

   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(o3_lbl,            1, mpiint, 0, mpicom)
   call mpibcast(o3_sfc,            1, mpir8,  0, mpicom)
   call mpibcast(o3_tau,            1, mpir8,  0, mpicom)
   call mpibcast(psc_T,             1, mpir8,  0, mpicom)
#endif

end subroutine linoz_readnl

!--------------------------------------------------------------------
!--------------------------------------------------------------------
  subroutine lin_strat_chem_inti()
    !
    ! initialize linearized stratospheric chemistry by reading
    ! input parameters from netcdf file and interpolate to
    ! present model grid
    !
    use linoz_data,   only : linoz_data_init, has_linoz_data
    use ppgrid,       only : pver
    use mo_chem_utls, only : get_spc_ndx
    use cam_history,  only : addfld, horiz_only, add_default

    implicit none

    if (.not.has_linoz_data) return

    !
    ! find index of ozone
    !
    index_o3  =  get_spc_ndx('O3')

    do_lin_strat_chem = has_linoz_data
     if ( index_o3 <= 0 ) then !added index_o3l
       write(iulog,*) ' No ozone in the chemical mechanism, skipping lin_strat_chem'
       do_lin_strat_chem = .false.
       return
    end if

    ! check for synoz

    if( get_spc_ndx( 'SYNOZ' ) > 0 .and. has_linoz_data) then
       call endrun('lin_strat_chem_inti: cannot have both synoz and linoz')
    endif

    ! initialize the linoz data

    call linoz_data_init()

    ! define additional output

    call addfld( 'LINOZ_DO3'    , (/ 'lev' /), 'A', '1/s'     , 'ozone vmr tendency by linearized ozone chemistry'   )
    call addfld( 'LINOZ_DO3_PSC', (/ 'lev' /), 'A', '1/s'     , 'ozone vmr loss by PSCs using Carille et al. (1990)' )
    call addfld( 'LINOZ_SSO3'   , (/ 'lev' /), 'A', 'kg'     , 'steady state ozone in LINOZ'                        )
    call addfld( 'LINOZ_O3COL'  , (/ 'lev' /), 'A', 'DU'     , 'ozone column above'                                 )
    call addfld( 'LINOZ_O3CLIM' , (/ 'lev' /), 'A', 'mol/mol', 'climatology of ozone in LINOZ'                      )
    call addfld( 'LINOZ_SZA'    ,    horiz_only, 'A', 'degrees', 'solar zenith angle in LINOZ'                      )
    call addfld( 'LINOZ_SFCSINK', horiz_only, 'A', 'Tg/yr/m2'   , 'surface o3 sink in LINOZ with an e-fold to a fixed concentration  ' )

    call add_default( 'LINOZ_DO3'    , 1, ' ' )
    call add_default( 'LINOZ_DO3_PSC', 1, ' ' )
    call add_default( 'LINOZ_SSO3'   , 1, ' ' )
    call add_default( 'LINOZ_O3COL'  , 1, ' ' )
    call add_default( 'LINOZ_O3CLIM' , 1, ' ' )
    call add_default( 'LINOZ_SZA'    , 1, ' ' )
    call add_default( 'LINOZ_SFCSINK', 1, ' ' )
    return
  end subroutine lin_strat_chem_inti


  !--------------------------------------------------------------------
  !--------------------------------------------------------------------

  subroutine lin_strat_chem_solve( ncol, lchnk, o3col, temp, sza, pmid, delta_t, rlats, ltrop, & !in
       linoz_o3_clim, linoz_t_clim, linoz_o3col_clim, linoz_PmL_clim, linoz_dPmL_dO3, linoz_dPmL_dT, & !in
       linoz_dPmL_dO3col, linoz_cariolle_psc, & !in
       o3_vmr) !in-out
    !
    ! this subroutine updates the ozone mixing ratio in the stratosphere
    ! using linearized chemistry
    !
    use ppgrid,        only : pcols, pver
    use cam_history,   only : outfld

    !!inten ins
    integer,  intent(in)                           :: ncol                ! number of columns in chunk
    integer,  intent(in)                           :: lchnk               ! chunk index
    real(r8), intent(in)   , dimension(ncol ,pver) :: o3col               ! ozone column above box [mol/cm^2]
    real(r8), intent(in)   , dimension(pcols,pver) :: temp                ! temperature [K]
    real(r8), intent(in)   , dimension(ncol )      :: sza                 ! local solar zenith angle
    real(r8), intent(in)   , dimension(pcols,pver) :: pmid                ! midpoint pressure [Pa]
    real(r8), intent(in)                           :: delta_t             ! timestep size [secs]
    real(r8), intent(in)                           :: rlats(ncol)         ! column latitudes (radians)
    integer,  intent(in)   , dimension(pcols)      :: ltrop               ! chunk index

    real(r8),  intent(in), dimension(:,:) :: linoz_o3_clim      ! ozone (climatology) [vmr]
    real(r8),  intent(in), dimension(:,:) :: linoz_t_clim       ! temperature (climatology) [K]
    real(r8),  intent(in), dimension(:,:) :: linoz_o3col_clim   ! Column O3 above box (climatology) [Dobson Units or DU]
    real(r8),  intent(in), dimension(:,:) :: linoz_PmL_clim     ! P minus L (climatology) [vmr/s]
    real(r8),  intent(in), dimension(:,:) :: linoz_dPmL_dO3     ! Sensitivity of P minus L to O3 [1/s]
    real(r8),  intent(in), dimension(:,:) :: linoz_dPmL_dT      ! Sensitivity of P minus L to T [K]
    real(r8),  intent(in), dimension(:,:) :: linoz_dPmL_dO3col  ! Sensitivity of P minus L to overhead O3 column [vmr/DU]
    real(r8),  intent(in), dimension(:,:) :: linoz_cariolle_psc ! Cariolle parameter for PSC loss of ozone [1/s]

    !inten in-outs
    real(r8), intent(inout), dimension(ncol ,pver) :: o3_vmr              ! ozone volume mixing ratio [vmr]

    !Local
    integer :: icol,kk !indices
    real(r8) :: o3col_du,delta_temp,delta_o3col,o3_old,o3_new,delta_o3
    real(r8) :: max_sza, psc_loss
    real(r8) :: o3_clim
    real(r8), dimension(ncol) :: lats
    real(r8), dimension(ncol,pver) :: do3_linoz,do3_linoz_psc,ss_o3,o3col_du_diag,o3clim_linoz_diag
    logical :: excess_chlorine
    !
    ! parameters
    !
    real(r8), parameter :: convert_to_du            = 1._r8/(2.687e16_r8)      ! convert ozone column from [mol/cm^2] to [DU]

    ! skip if no ozone field available
    if ( .not. do_lin_strat_chem ) return

    ! initialize output arrays
    do3_linoz         = 0._r8
    do3_linoz_psc     = 0._r8
    o3col_du_diag     = 0._r8
    o3clim_linoz_diag = 0._r8
    ss_o3             = 0._r8

    lats = rlats * radians_to_degrees ! convert lats from radians to degrees

    !is there more than the background chlorine?
    excess_chlorine = (chlorine_loading-chlorine_loading_bgnd) > 0._r8

    LOOP_COL: do icol = 1, ncol
       LOOP_LEV: do kk = 1, ltrop(icol)

          o3_clim = linoz_o3_clim(icol,kk)     ! climatological ozone
          o3clim_linoz_diag(icol,kk) = o3_clim ! diagnostic for output
          o3_old = o3_vmr(icol,kk)             ! old ozone mixing ratio
          o3col_du = o3col(icol,kk) * convert_to_du ! convert o3col from mol/cm2
          o3col_du_diag(icol,kk) = o3col_du    !update diagnostic output

          ! compute differences from climatology
          delta_temp  = temp(icol,kk) - linoz_t_clim    (icol,kk)
          delta_o3col = o3col_du  - linoz_o3col_clim(icol,kk)

          ! steady state ozone
          ss_o3(icol,kk) = o3_clim - (               linoz_PmL_clim   (icol,kk)   &
               + delta_o3col * linoz_dPmL_dO3col(icol,kk)   &
               + delta_temp  * linoz_dPmL_dT    (icol,kk)   &
               ) / linoz_dPmL_dO3   (icol,kk)

          delta_o3 = (ss_o3(icol,kk)-o3_old) * (1._r8 - exp(linoz_dPmL_dO3(icol,kk)*delta_t)) ! ozone change

          o3_new = o3_old + delta_o3 ! define new ozone mixing ratio

          do3_linoz(icol,kk) = delta_o3/delta_t ! output diagnostic

          ! PSC activation (follows Cariolle et al 1990.)
          call psc_activation( lats(icol), temp(icol,kk), pmid(icol,kk), sza(icol), linoz_cariolle_psc(icol,kk), delta_t, & !in
               excess_chlorine, o3_old, &  !in
               o3_new, do3_linoz_psc(icol,kk)) !out

          o3_vmr(icol,kk) = o3_new ! update ozone vmr

       enddo LOOP_LEV
    enddo LOOP_COL

    ! output
    call outfld( 'LINOZ_DO3'    , do3_linoz              , ncol, lchnk )
    call outfld( 'LINOZ_DO3_PSC', do3_linoz_psc          , ncol, lchnk )
    call outfld( 'LINOZ_SSO3'   , ss_o3                  , ncol, lchnk )
    call outfld( 'LINOZ_O3COL'  , o3col_du_diag          , ncol, lchnk )
    call outfld( 'LINOZ_O3CLIM' , o3clim_linoz_diag      , ncol, lchnk )
    call outfld( 'LINOZ_SZA'    ,(sza*radians_to_degrees), ncol, lchnk )

    return
  end subroutine lin_strat_chem_solve

!------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------

  subroutine psc_activation( lats, temp, pmid, sza, linoz_cariolle_psc, delta_t, & !in
       excess_chlorine, o3_old, &     !in
       o3_new, do3_linoz_psc) !out

    ! PSC activation (follows Cariolle et al 1990.)
    implicit none

    !intent-ins
    real(r8), intent(in) :: lats !lattitude [degree]
    real(r8), intent(in) :: temp !Temperature [K]
    real(r8), intent(in) :: pmid !Midpoint Pressure [Pa]
    real(r8), intent(in) :: sza  !local solar zenith angle
    real(r8), intent(in) :: linoz_cariolle_psc ! Cariolle parameter for PSC loss of ozone [1/s]
    real(r8), intent(in) :: delta_t ! timestep size [secs]
    logical,  intent(in) :: excess_chlorine !.TRUE. if chlorine is more than the backgroud chlorine
    real(r8), intent(in) :: o3_old ! ozone volume mixing ratio [vmr]

    !intent-outs
    real(r8), intent(out) :: o3_new ! ozone volume mixing ratio [vmr]
    real(r8), intent(out) :: do3_linoz_psc ![vmr/s]

    !local
    real(r8), parameter :: lats_threshold        = 40.0_r8
    real(r8), parameter :: chlorine_loading_1987 = 2.5977_r8     ! EESC value [ppbv]

    real(r8) :: max_sza, psc_loss

    ! use only if abs(latitude) > lats_threshold
    if ( abs(lats) > lats_threshold ) then
       if ( excess_chlorine ) then
          if ( temp <= psc_T ) then
             ! define maximum SZA for PSC loss (= tangent height at sunset)
             max_sza = (90._r8 + sqrt( max( 16._r8*log10(100000._r8/pmid),0._r8)))

             if ( (sza*radians_to_degrees) <= max_sza ) then

                psc_loss = exp(-linoz_cariolle_psc &
                     * (chlorine_loading/chlorine_loading_1987)**2 &
                     * delta_t )

                o3_new = o3_old * psc_loss

                do3_linoz_psc = (o3_new-o3_old)/delta_t ! output diagnostic
             endif
          endif
       endif
    endif

  end subroutine psc_activation

!------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------

  subroutine lin_strat_sfcsink( ncol, lchnk, delta_t, pdel, & !in
       o3l_vmr) !in-out

    use ppgrid,        only : pcols, pver
    use physconst,     only : mw_air => mwdry, &
         mw_o3  => mwo3   !in grams

    use mo_constants, only : rgrav
    use cam_history,   only : outfld

    implicit none

    !intent-ins
    integer,  intent(in)   :: ncol                ! number of columns in chunk
    integer,  intent(in)   :: lchnk               ! chunk index
    real(r8), intent(in)   :: delta_t             ! timestep size [secs]
    real(r8), intent(in)   :: pdel(ncol,pver)     ! pressure delta about midpoints [Pa]

    !inten in-outs
    real(r8), intent(inout):: o3l_vmr(ncol ,pver)             ! ozone volume mixing ratio [vmr]

    !local
    !Two parameters are applied to Linoz O3 for surface sink, O3l is not coupled to real O3
    real(r8), parameter :: KgtoTg                   = 1.0e-9_r8
    real(r8), parameter :: peryear                  = 86400._r8* 365.0_r8 ! to multiply to convert per second to per year

    real(r8) :: mass(ncol,pver)
    real(r8) :: o3l_old, o3l_new, efactor, do3
    real(r8), dimension(ncol)  :: do3mass, o3l_sfcsink
    integer icol, kk

    if ( .not. do_lin_strat_chem ) return

    !initializing array
    o3l_sfcsink = 0._r8

    do kk = 1, pver
       mass(:ncol,kk) = pdel(:ncol,kk) * rgrav  ! air mass in kg/m2
    enddo

    efactor  = (1.d0 - exp(-delta_t/o3_tau)) !compute time scale factor

    do icol = 1, ncol

       do3mass(icol) =0._r8

       do kk = pver, pver-o3_lbl+1, -1

          o3l_old = o3l_vmr(icol,kk)  !vmr

          do3 =  (o3_sfc - o3l_old)* efactor !vmr

          o3l_new  = o3l_old + do3

          do3mass(icol) = do3mass(icol) + do3* mass(icol,kk) * mw_o3/mw_air  ! loss in kg/m2 summed over boundary layers within one time step

          o3l_vmr(icol,kk) = o3l_new

       enddo
    enddo

    o3l_sfcsink(:ncol) = do3mass(:ncol)/delta_t * KgtoTg * peryear ! saved in Tg/yr/m2 unit

    call outfld('LINOZ_SFCSINK', o3l_sfcsink, ncol, lchnk)

    return

  end subroutine lin_strat_sfcsink

end module lin_strat_chem
