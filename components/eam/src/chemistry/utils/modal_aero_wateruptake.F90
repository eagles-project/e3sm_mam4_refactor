module modal_aero_wateruptake

!   RCE 07.04.13:  Adapted from MIRAGE2 code

use shr_kind_mod,     only: r8 => shr_kind_r8
use physconst,        only: pi, rhoh2o
use ppgrid,           only: pcols, pver
use physics_types,    only: physics_state
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field
use modal_aero_data,  only: nmodes => ntot_amode
use modal_aero_data,  only: nspec_amode, specdens_amode, spechygro, sigmag_amode,  &
                            rhcrystal_amode, rhdeliques_amode, lspectype_amode
use cam_history,      only: addfld, add_default, outfld
use ref_pres,         only: top_lev => clim_modal_aero_top_lev
use phys_control,     only: phys_getopts
use cam_abortutils,   only: endrun
use mam_support,      only: min_max_bound, assign_la_lc


implicit none
private
save

public :: &
   modal_aero_wateruptake_init, &
   modal_aero_wateruptake_dr, &
   modal_aero_kohler, &
   modal_aero_wateruptake_reg


! Physics buffer indices
integer :: cld_idx        = 0
integer :: dgnum_idx      = 0
integer :: dgnumwet_idx   = 0
integer :: wetdens_ap_idx = 0
integer :: qaerwat_idx    = 0

logical :: pergro_mods         = .false.

real(r8), parameter   :: third = 1._r8/3._r8
real(r8), parameter   :: pi43  = pi*4.0_r8/3.0_r8

real(r8), allocatable :: maer(:,:,:)      ! aerosol wet mass MR (including water) (kg/kg-air)
real(r8), allocatable :: hygro(:,:,:)     ! volume-weighted mean hygroscopicity (--)
real(r8), allocatable :: naer(:,:,:)      ! aerosol number MR (bounded!) (#/kg-air)
real(r8), allocatable :: dryvol(:,:,:)    ! single-particle-mean dry volume (m3)
real(r8), allocatable :: drymass(:,:,:)   ! single-particle-mean dry mass  (kg)
real(r8), allocatable :: dryrad(:,:,:)    ! dry volume mean radius of aerosol (m)

real(r8), allocatable :: wetrad(:,:,:)    ! wet radius of aerosol (m)
real(r8), allocatable :: wetvol(:,:,:)    ! single-particle-mean wet volume (m3)
real(r8), allocatable :: wtrvol(:,:,:)    ! single-particle-mean water volume in wet aerosol (m3)

real(r8), allocatable :: rhcrystal(:)
real(r8), allocatable :: rhdeliques(:)
real(r8), allocatable :: specdens_1(:)
!$OMP THREADPRIVATE(maer, hygro, naer, dryvol, drymass, dryrad, wetrad, wetvol, wtrvol, rhcrystal, rhdeliques, specdens_1)

!===============================================================================
contains
!===============================================================================

subroutine modal_aero_wateruptake_reg()

  use physics_buffer,   only: pbuf_add_field, dtype_r8

   
   call pbuf_add_field('DGNUMWET',   'global',  dtype_r8, (/pcols, pver, nmodes/), dgnumwet_idx)
   call pbuf_add_field('WETDENS_AP', 'physpkg', dtype_r8, (/pcols, pver, nmodes/), wetdens_ap_idx)

   ! 1st order rate for direct conversion of strat. cloud water to precip (1/s)
   call pbuf_add_field('QAERWAT',    'physpkg', dtype_r8, (/pcols, pver, nmodes/), qaerwat_idx)  

end subroutine modal_aero_wateruptake_reg

!===============================================================================
!===============================================================================

subroutine modal_aero_wateruptake_init(pbuf2d)
   use time_manager,   only: is_first_step
   use physics_buffer, only: pbuf_set_field
   use shr_log_mod ,   only: errmsg => shr_log_errmsg

   type(physics_buffer_desc), pointer :: pbuf2d(:,:)

   integer :: imode, istat
   logical :: history_aerosol      ! Output the MAM aerosol variables and tendencies
   logical :: history_verbose      ! produce verbose history output

   character(len=3) :: trnum       ! used to hold mode number (as characters)
   !----------------------------------------------------------------------------

   cld_idx        = pbuf_get_index('CLD')    
   dgnum_idx      = pbuf_get_index('DGNUM')    


   !$OMP PARALLEL
   allocate(maer(pcols,pver,nmodes),stat=istat)
   if (istat .ne. 0) call endrun("Unable to allocate maer:       "//errmsg(__FILE__,__LINE__) )
   allocate(hygro(pcols,pver,nmodes),   stat=istat)
   if (istat .ne. 0) call endrun("Unable to allocate hygro:      "//errmsg(__FILE__,__LINE__) )
   allocate(naer(pcols,pver,nmodes),    stat=istat)
   if (istat .ne. 0) call endrun("Unable to allocate naer:       "//errmsg(__FILE__,__LINE__) )
   allocate(dryvol(pcols,pver,nmodes),  stat=istat)
   if (istat .ne. 0) call endrun("Unable to allocate dryvol:     "//errmsg(__FILE__,__LINE__) )
   allocate(drymass(pcols,pver,nmodes), stat=istat)
   if (istat .ne. 0) call endrun("Unable to allocate drymass:    "//errmsg(__FILE__,__LINE__) )
   allocate(dryrad(pcols,pver,nmodes),  stat=istat)
   if (istat .ne. 0) call endrun("Unable to allocate dryradr:    "//errmsg(__FILE__,__LINE__) )
   allocate(wetrad(pcols,pver,nmodes),  stat=istat)
   if (istat .ne. 0) call endrun("Unable to allocate wetrad:     "//errmsg(__FILE__,__LINE__) )
   allocate(wetvol(pcols,pver,nmodes),  stat=istat)
   if (istat .ne. 0) call endrun("Unable to allocate wetvol:     "//errmsg(__FILE__,__LINE__) )
   allocate(wtrvol(pcols,pver,nmodes),  stat=istat)
   if (istat .ne. 0) call endrun("Unable to allocate wtrvol:     "//errmsg(__FILE__,__LINE__) )
   allocate(rhcrystal(nmodes),          stat=istat)
   if (istat .ne. 0) call endrun("Unable to allocate rhcrystal:  "//errmsg(__FILE__,__LINE__) )
   allocate(rhdeliques(nmodes),         stat=istat)
   if (istat .ne. 0) call endrun("Unable to allocate rhdeliques: "//errmsg(__FILE__,__LINE__) )
   allocate(specdens_1(nmodes),         stat=istat)
   if (istat .ne. 0) call endrun("Unable to allocate specdens_1: "//errmsg(__FILE__,__LINE__) )
   !$OMP END PARALLEL


   ! determine default variables
   call phys_getopts(history_aerosol_out = history_aerosol, &
                     history_verbose_out = history_verbose, &
                     pergro_mods_out = pergro_mods)

   do imode = 1, nmodes
      write(trnum, '(i3.3)') imode
      call addfld('dgnd_a'//trnum(2:3), (/ 'lev' /), 'A', 'm', &
         'dry dgnum, interstitial, mode '//trnum(2:3))
      call addfld('dgnw_a'//trnum(2:3), (/ 'lev' /), 'A', 'm', &
         'wet dgnum, interstitial, mode '//trnum(2:3))
      call addfld('wat_a'//trnum(3:3), (/ 'lev' /), 'A', 'm', &
         'aerosol water, interstitial, mode '//trnum(2:3))
      
      if (history_aerosol) then  
         if (history_verbose) then
            call add_default('dgnd_a'//trnum(2:3), 1, ' ')
            call add_default('dgnw_a'//trnum(2:3), 1, ' ')
            call add_default('wat_a'//trnum(3:3),  1, ' ')
         endif
      endif

   enddo

   ! Add total aerosol water
   if (history_aerosol .and. .not. history_verbose) then
      call addfld('aero_water', (/ 'lev' /), 'A', 'm', &
         'sum of aerosol water of interstitial modes wat_a1+wat_a2+wat_a3+wat_a4' )
      call add_default( 'aero_water',  1, ' ')
   endif
   
   if (is_first_step()) then
      ! initialize fields in physics buffer
      call pbuf_set_field(pbuf2d, dgnumwet_idx, 0.0_r8)
   endif

end subroutine modal_aero_wateruptake_init

!===============================================================================
subroutine modal_aero_wateruptake_dr(state, pbuf, list_idx_in, dgnumdry_m, & ! in
                                        dgnumwet_m, qaerwat_m  ) ! inout
   !----------------------------------------------------------------------------
   !
   ! CAM specific driver for modal aerosol water uptake code.
   !
   ! *** N.B. *** The calculation has been enabled for diagnostic mode lists
   !              via optional arguments.  If the list_idx arg is present then
   !              all the optional args must be present.
   !
   !----------------------------------------------------------------------------

   ! Arguments
   type(physics_state), target, intent(in)    :: state          ! Physics state variables
   type(physics_buffer_desc),   pointer       :: pbuf(:)        ! physics buffer
   ! Optional inputs for diagnostic mode
   integer,  optional,          intent(in)    :: list_idx_in
   real(r8), optional, allocatable, target, intent(in)    :: dgnumdry_m(:,:,:)
   real(r8), optional, allocatable, target, intent(inout)   :: dgnumwet_m(:,:,:)
   real(r8), optional, allocatable, target, intent(inout)   :: qaerwat_m(:,:,:)

   ! local variables
   integer :: lchnk              ! chunk index
   integer :: ncol               ! number of columns
   integer :: list_idx           ! radiative constituents list index
   integer :: imode              ! mode index
   integer :: itim_old           ! index

   real(r8), pointer :: h2ommr(:,:)      ! specific humidity [kg/kg]
   real(r8), pointer :: temperature(:,:) ! temperatures [K]
   real(r8), pointer :: pmid(:,:)        ! layer pressure [Pa]
   real(r8), pointer :: cldn(:,:)        ! layer cloud fraction [fraction]
   real(r8), pointer :: dgncur_a(:,:,:)  ! aerosol particle diameter [m]
   real(r8), pointer :: dgncur_awet(:,:,:) ! wet aerosol diameter [m]
   real(r8), pointer :: wetdens(:,:,:)   ! wet aerosol density [kg/m3]
   real(r8), pointer :: qaerwat(:,:,:)   ! aerosol water [kg/kg]

   real(r8) :: rh(pcols,pver)        ! relative humidity [fraction]
   logical  :: compute_wetdens

   ! variables for writing history output
   real(r8) :: aerosol_water(pcols,pver) !sum of aerosol water (wat_a1 + wat_a2 + wat_a3 + wat_a4) [kg/kg]
   logical  :: history_aerosol      ! Output the MAM aerosol variables and tendencies
   logical  :: history_verbose      ! produce verbose history output
   character(len=3)  :: trnum       ! used to hold mode number (as characters)

   !----------------------------------------------------------------------------
   lchnk = state%lchnk
   ncol  = state%ncol
   h2ommr       => state%q(:,:,1)
   temperature  => state%t
   pmid         => state%pmid

   list_idx = 0
   if (present(list_idx_in))     list_idx = list_idx_in

   !initialize to an invalid value
   naer(:,:,:)    = huge(1.0_r8)
   dryvol(:,:,:)  = huge(1.0_r8)
   drymass(:,:,:) = huge(1.0_r8)
   dryrad(:,:,:)  = huge(1.0_r8)
   wetrad(:,:,:)  = huge(1.0_r8)
   wetvol(:,:,:)  = huge(1.0_r8)
   wtrvol(:,:,:)  = huge(1.0_r8)
   rhcrystal(:)   = huge(1.0_r8)
   rhdeliques(:)  = huge(1.0_r8)
   specdens_1(:)  = huge(1.0_r8)
   maer(:,:,:)    = 0._r8
   hygro(:,:,:)   = 0._r8

   ! determine default variables
   call phys_getopts(history_aerosol_out = history_aerosol, &
                     history_verbose_out = history_verbose)

   !by default set compute_wetdens to be true
   compute_wetdens = .true.
   if (.not. present(list_idx_in)) then
      call pbuf_get_field(pbuf, dgnum_idx,      dgncur_a )
      call pbuf_get_field(pbuf, dgnumwet_idx,   dgncur_awet )
      call pbuf_get_field(pbuf, wetdens_ap_idx, wetdens)
      call pbuf_get_field(pbuf, qaerwat_idx,    qaerwat)
   else
      dgncur_a    => dgnumdry_m
      dgncur_awet => dgnumwet_m
      qaerwat     => qaerwat_m
      !set compute_wetdens to flase if wetdens is not present
      compute_wetdens = .false.
   endif

   itim_old    =  pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, cld_idx, cldn, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

   !----------------------------------------------------------------------------
   ! retreive aerosol properties
   call modal_aero_wateruptake_dryaer( ncol, list_idx,            & ! in
                                      state, pbuf,   dgncur_a     ) ! in

   !----------------------------------------------------------------------------
   ! estimate clear air relative humidity using cloud fraction

   call modal_aero_wateruptake_rh_clearair( ncol,         & ! in
                        temperature, pmid,  h2ommr, cldn, & ! in
                        rh                                ) ! out

   !----------------------------------------------------------------------------
   ! compute wet aerosol properties

   ! compute aerosol wet radius, volume, diameter and aerosol water
   call modal_aero_wateruptake_wetaer( ncol,                    & ! in
                rhcrystal,      rhdeliques,     dgncur_a,       & ! in
                dryrad, hygro,  rh,     naer,   dryvol,         & ! in
                wetrad, wetvol, wtrvol, dgncur_awet,  qaerwat   ) ! out

   ! compute wet aerosol density
   if (compute_wetdens) then
        call modal_aero_wateruptake_wetdens( ncol,           & ! in
             wetvol, wtrvol, drymass,      specdens_1,       & ! in
             wetdens                                         ) ! inout
   endif

   !----------------------------------------------------------------------------
   ! write history output if not in diagnostic mode

   if (.not.present(list_idx_in)) then

      aerosol_water(:ncol,:) = 0._r8
      do imode = 1, nmodes
         ! output to history
         write( trnum, '(i3.3)' ) imode
         call outfld( 'wat_a'//trnum(3:3),  qaerwat(:,:,imode),     pcols, lchnk)
         call outfld( 'dgnd_a'//trnum(2:3), dgncur_a(:,:,imode),    pcols, lchnk)
         call outfld( 'dgnw_a'//trnum(2:3), dgncur_awet(:,:,imode), pcols, lchnk)
         if (history_aerosol .and. .not. history_verbose) &
         aerosol_water(:ncol,:) = aerosol_water(:ncol,:) + qaerwat(:ncol,:,imode)
      enddo ! nmodes

      if (history_aerosol .and. .not. history_verbose) &
         call outfld( 'aero_water',  aerosol_water(:ncol,:),    ncol, lchnk)

   endif

end subroutine modal_aero_wateruptake_dr

!===============================================================================
subroutine modal_aero_wateruptake_dryaer( ncol, list_idx,            & ! in
                                     state,     pbuf,   dgncur_a     ) ! in
!-----------------------------------------------------------------------
! retreive dry aerosol properties
! output variables are declared in the module so no output variables in this subroutine
!-----------------------------------------------------------------------

   integer,  intent(in)  :: ncol                ! number of columns
   integer,  intent(in)  :: list_idx            ! radiative constituents list index
   real(r8), intent(in)  :: dgncur_a(:,:,:)     ! dry aerosol diameter [m]
   type(physics_state), target, intent(in)  :: state          ! Physics state variables
   type(physics_buffer_desc), pointer       :: pbuf(:)        ! physics buffer

   ! local variables
   integer      :: imode, ispec, kk, icol,la,lc
   integer      :: nspec
   real(r8)     :: sigmag
   real(r8)     :: alnsg                        ! log(sigmag)
   real(r8)     :: specdens                     ! aerosol density [kg/m3]
   real(r8)     :: spechygro_1, spechygro_i     ! aerosol hygroscopicity [unitless]
   real(r8)     :: v2ncur_a                     ! 1 / mean particle volume [1/m3]
   real(r8)     :: drydens                      ! dry particle density [kg/m3]
   real(r8)     :: dryvolmr(pcols,pver)         ! volume MR for aerosol mode [m3/kg]
   real(r8)     :: vol_tmp                      ! temporary aerosol volume [m3/kg]
   real(r8)     :: raer(pcols,pver)             ! aerosol species mixing ratio [kg/kg]
   real(r8), parameter  :: small_value = 1.0e-30_r8
   real(r8), parameter  :: small_value_31 = 1.0e-31_r8


   do imode = 1, nmodes
      ! initiate variable
      dryvolmr(:,:) = 0._r8

      ! get mode properties
      sigmag = sigmag_amode(imode)
      rhcrystal = rhcrystal_amode(imode)
      rhdeliques = rhdeliques_amode(imode)
      nspec = nspec_amode(imode)
      ! specdens_1 is defined in module and will be used in later subroutine
      specdens_1 = specdens_amode(lspectype_amode(1,imode)) 
      spechygro_1 = spechygro(lspectype_amode(1,imode))

      alnsg = log(sigmag)

      do ispec = 1, nspec
         ! get species interstitial mixing ratio ('a')
         spechygro_i = spechygro(lspectype_amode(ispec,imode))
         specdens = specdens_amode(lspectype_amode(ispec,imode))
         call assign_la_lc( imode,      ispec,          & ! in
                            la,         lc              ) ! out
         raer = state%q(:,:,la)

         do kk = top_lev, pver
            do icol = 1, ncol
               vol_tmp             = raer(icol,kk)/specdens
               maer(icol,kk,imode) = maer(icol,kk,imode) + raer(icol,kk)
               dryvolmr(icol,kk)   = dryvolmr(icol,kk) + vol_tmp
               ! hygro currently is sum(hygro * volume) of each species,
               ! need to divided by sum(volume) later to get mean hygro for all species.
               hygro(icol,kk,imode)  = hygro(icol,kk,imode) + vol_tmp*spechygro_i
            enddo ! ncol
         enddo ! kk
      enddo ! nspec


      do kk = top_lev, pver
         do icol = 1, ncol

            if (dryvolmr(icol,kk) > small_value) then
               ! divided by sum(volume) to get mean hygro
               hygro(icol,kk,imode) = hygro(icol,kk,imode)/dryvolmr(icol,kk)
            else
               hygro(icol,kk,imode) = spechygro_1
            endif

            v2ncur_a = 1._r8 / ((pi/6._r8)*(dgncur_a(icol,kk,imode)**3._r8)*exp(4.5_r8*alnsg**2._r8) )
            ! naer = aerosol number (#/kg)
            naer(icol,kk,imode) = dryvolmr(icol,kk)*v2ncur_a

            ! compute mean (1 particle) dry volume and mass for each mode
            ! old coding is replaced because the new (1/v2ncur_a) is equal to
            ! the mean particle volume
            ! also moletomass forces maer >= 1.0e-30, so (maer/dryvolmr)
            ! should never cause problems (but check for maer < 1.0e-31 anyway)
            if (maer(icol,kk,imode) .gt. small_value_31) then
               drydens = maer(icol,kk,imode)/dryvolmr(icol,kk)
            else
               drydens = 1.0_r8
            endif

            ! C++ porting note: these are output but defined in the module
            ! thus not in the subroutine output
            dryvol(icol,kk,imode)   = 1.0_r8/v2ncur_a
            drymass(icol,kk,imode)  = drydens*dryvol(icol,kk,imode)
            dryrad(icol,kk,imode)   = (dryvol(icol,kk,imode)/pi43)**third

         enddo !  ncol
      enddo ! kk

   enddo    ! modes

end subroutine modal_aero_wateruptake_dryaer

!===============================================================================
subroutine modal_aero_wateruptake_rh_clearair( ncol,         & ! in
                     temperature,    pmid,   h2ommr,  cldn,  & ! in
                     rh                                      ) ! out
!-----------------------------------------------------------------------
! estimate clear air relative humidity using cloud fraction
!-----------------------------------------------------------------------
use wv_saturation,    only: qsat_water

   integer,  intent(in) :: ncol                 ! number of columns
   real(r8), intent(in) :: temperature(:,:)     ! temperature [K]
   real(r8), intent(in) :: pmid(:,:)            ! layer pressure [Pa]
   real(r8), intent(in) :: h2ommr(:,:)          ! water mass mixing ratio [kg/kg]
   real(r8), intent(in) :: cldn(:,:)            ! layer cloud fraction [fraction]
   real(r8), intent(out) :: rh(:,:)           ! relative humidity [Fraction]

   ! local variables
   integer   :: icol, kk
   real(r8)  :: rh_max = 0.98_r8        ! maximum relative humidity
   real(r8)  :: es(pcols)               ! saturation vapor pressure [Pa]
   real(r8)  :: qs(pcols)               ! saturation specific humidity [kg/kg]
   real(r8), parameter  :: cldn_thresh = 1.0_r8      ! threshold cloud fraction

   do kk = top_lev, pver
      call qsat_water(temperature(:ncol,kk), pmid(:ncol,kk), es(:ncol), qs(:ncol))
      do icol = 1, ncol
         if (qs(icol) > h2ommr(icol,kk)) then
            rh(icol,kk) = h2ommr(icol,kk)/qs(icol)
         else
            rh(icol,kk) = rh_max
         endif
         rh(icol,kk) = min_max_bound(0.0_r8, rh_max, rh(icol,kk))

         if (cldn(icol,kk) .lt. cldn_thresh) then
            rh(icol,kk) = (rh(icol,kk) - cldn(icol,kk)) / &
                          (1.0_r8 - cldn(icol,kk))  ! RH of clear portion
         endif
         rh(icol,kk) = max(rh(icol,kk), 0.0_r8)
      enddo ! i = 1, ncol
   enddo ! k = top_lev, pver

end subroutine modal_aero_wateruptake_rh_clearair

!===============================================================================
subroutine modal_aero_wateruptake_wetaer(   ncol,         & ! in
             rhcrystal,      rhdeliques,    dgncur_a,     & ! in
             dryrad, hygro,  rh,     naer,  dryvol,       & ! in
             wetrad, wetvol, wtrvol, dgncur_awet, qaerwat ) ! out
!-----------------------------------------------------------------------
!
! Purpose: Compute aerosol wet radius and other properties
!
! Method:  Kohler theory
!
! Author:  S. Ghan
!
!-----------------------------------------------------------------------

   ! Arguments
   integer, intent(in)  :: ncol                    ! number of columns

   real(r8), intent(in) :: rhcrystal(:)          ! crystal RH [fraction]
   real(r8), intent(in) :: rhdeliques(:)         ! deliques RH [Fraction]
   real(r8), intent(in) :: dgncur_a(:,:,:)       ! dry aerosol diameter [m]
   real(r8), intent(in) :: dryrad(:,:,:)         ! dry volume mean radius of aerosol [m]
   real(r8), intent(in) :: hygro(:,:,:)          ! volume-weighted mean hygroscopicity [unitless]
   real(r8), intent(in) :: rh(:,:)               ! relative humidity [fraction]
   real(r8), intent(in) :: naer(:,:,:)           ! number of aerosols [#/kg-air]
   real(r8), intent(in) :: dryvol(:,:,:)         ! dry aerosol volume [m^3]

   real(r8), intent(out) :: wetrad(:,:,:)        ! wet radius of aerosol [m]
   real(r8), intent(out) :: wetvol(:,:,:)        ! single-particle-mean wet volume [m^3]
   real(r8), intent(out) :: wtrvol(:,:,:)        ! single-particle-mean water volume in wet aerosol [m^3]
   real(r8), intent(out) :: dgncur_awet(:,:,:)   ! wet aerosol diameter [m]
   real(r8), intent(out) :: qaerwat(:,:,:)       ! aerosol water [kg/kg]
   ! local variables
   integer :: icol, kk, imode
   real(r8) :: hystfac                ! working variable for hysteresis

   !-----------------------------------------------------------------------
   ! loop over all aerosol modes
   do imode = 1, nmodes

      hystfac = 1.0_r8 / max(1.0e-5_r8, (rhdeliques(imode) - rhcrystal(imode)))

      do kk = top_lev, pver
         do icol = 1, ncol

            ! compute wet radius for each mode
            call modal_aero_kohler(dryrad(icol,kk,imode),               & ! in
                                hygro(icol,kk,imode), rh(icol,kk),      & ! in
                                wetrad(icol,kk,imode)                   ) ! out

            wetrad(icol,kk,imode) = max(wetrad(icol,kk,imode), dryrad(icol,kk,imode))
            wetvol(icol,kk,imode) = pi43*wetrad(icol,kk,imode)**3
            wetvol(icol,kk,imode) = max(wetvol(icol,kk,imode), dryvol(icol,kk,imode))
            wtrvol(icol,kk,imode) = wetvol(icol,kk,imode) - dryvol(icol,kk,imode)
            wtrvol(icol,kk,imode) = max(wtrvol(icol,kk,imode), 0.0_r8)

            ! apply simple treatment of deliquesence/crystallization hysteresis
            ! for rhcrystal < rh < rhdeliques, aerosol water is a fraction of
            ! the "upper curve" value, and the fraction is a linear function of rh
            if (rh(icol,kk) < rhcrystal(imode)) then
               wetrad(icol,kk,imode) = dryrad(icol,kk,imode)
               wetvol(icol,kk,imode) = dryvol(icol,kk,imode)
               wtrvol(icol,kk,imode) = 0.0_r8
            elseif (rh(icol,kk) < rhdeliques(imode)) then
               wtrvol(icol,kk,imode) = wtrvol(icol,kk,imode)*hystfac*(rh(icol,kk) - rhcrystal(imode))
               wtrvol(icol,kk,imode) = max(wtrvol(icol,kk,imode), 0.0_r8)
               wetvol(icol,kk,imode) = dryvol(icol,kk,imode) + wtrvol(icol,kk,imode)
               wetrad(icol,kk,imode) = (wetvol(icol,kk,imode)/pi43)**third
            endif

            ! calculate wet aerosol diameter and aerosol water
            dgncur_awet(icol,kk,imode) = dgncur_a(icol,kk,imode) * &
                                        (wetrad(icol,kk,imode)/dryrad(icol,kk,imode))
            qaerwat(icol,kk,imode) = rhoh2o*naer(icol,kk,imode)*wtrvol(icol,kk,imode)

         enddo  ! columns
      enddo     ! levels

   enddo ! modes

end subroutine modal_aero_wateruptake_wetaer

!===============================================================================
subroutine modal_aero_kohler( rdry_in, hygro, rh,    & ! in
                              rwet_out               ) ! out
!-----------------------------------------------------------------------
! calculates equlibrium radius r of haze droplets as function of
! dry particle mass and relative humidity s using kohler solution
! given in pruppacher and klett (eqn 6-35)

! for multiple aerosol types, assumes an internal mixture of aerosols
!-----------------------------------------------------------------------
   implicit none

! arguments
   real(r8), intent(in) :: rdry_in    ! aerosol dry radius [m]
   real(r8), intent(in) :: hygro      ! aerosol volume-mean hygroscopicity [unitless]
   real(r8), intent(in) :: rh          ! relative humidity (1 = saturated) [fraction]
   real(r8), intent(out) :: rwet_out   ! aerosol wet radius [m]

! local variables
   integer :: nn, nsol

   real(r8) :: aa, bb, pp   ! parameters
   real(r8) :: p40,p41,p42,p43 ! coefficients of polynomial
   real(r8) :: rwet        ! wet radius [microns]
   real(r8) :: rdry        ! radius of dry particle [microns]
   real(r8) :: ss          ! relative humidity (1 = saturated) [fraction]
   real(r8) :: slog        ! log relative humidity
   real(r8) :: vol         ! total volume of particle [microns**3]

   complex(r8) :: cx4(4) ! output of polynomials

   real(r8), parameter :: eps = 1.e-4_r8
   real(r8), parameter :: mw = 18._r8
   real(r8), parameter :: rhow = 1._r8
   real(r8), parameter :: surften = 76._r8
   real(r8), parameter :: tair = 273._r8    
   real(r8), parameter :: third = 1._r8/3._r8
   real(r8), parameter :: ugascon = 8.3e7_r8
   real(r8), parameter :: factor_um2m = 1.e-6_r8  ! convert micron to m
   real(r8), parameter :: factor_m2um = 1.e6_r8   ! convert m to micron
   real(r8), parameter :: small_rh = 1.e-10_r8    ! small value to avoid zero
   real(r8), parameter :: rmax = 30._r8           ! upper bound of radius [microns]

!  effect of organics on surface tension is neglected
   aa=2.e4_r8*mw*surften/(ugascon*tair*rhow)
 
   rdry = rdry_in*factor_m2um   ! convert (m) to (microns)
   vol = rdry**3          ! vol is r**3, not volume
   bb = vol*hygro

!  quartic
   ss = min_max_bound(small_rh, 1._r8-eps, rh)
   slog = log(ss)
   p43 = -aa/slog
   p42 = 0._r8
   p41 = bb/slog-vol
   p40 = aa*vol/slog

   pp = abs(-bb/aa)/(rdry*rdry)
   if(pp.lt.eps)then
!      approximate solution for small particles
       rwet = rdry * (1._r8 + pp*third/(1._r8-slog*rdry/aa))
   else
       call makoh_quartic(cx4,p43,p42,p41,p40)
!      find smallest real(r8) solution
       call find_real_solution( rdry, cx4, rwet,nsol)
   endif


! bound and convert from microns to m
   rwet = min(rwet, rmax) ! upper bound based on 1 day lifetime
   rwet_out = rwet*factor_um2m

   return
end subroutine modal_aero_kohler

!===============================================================================
subroutine makoh_quartic( cx, p3, p2, p1, p0)
!-----------------------------------------------------------------------
!     solves x**4 + p3 x**3 + p2 x**2 + p1 x + p0 = 0
!     where p0, p1, p2, p3 are real
!-----------------------------------------------------------------------
   real(r8), intent(in)     :: p0, p1, p2, p3 
   complex(r8), intent(out) :: cx(4)

   real(r8) :: third, qq, rr    ! temporary variables
   complex(r8) :: cb, cb0, cb1, crad, cy, czero ! temporary variables

   ! set complex zeros and 1/3 values
   czero=cmplx(0.0_r8,0.0_r8,r8)
   third=1._r8/3._r8


   qq=-p2*p2/36._r8+(p3*p1-4._r8*p0)/12._r8
   rr=-(p2/6._r8)**3+p2*(p3*p1-4._r8*p0)/48._r8 + (4._r8*p0*p2-p0*p3*p3-p1*p1)/16._r8

   crad=rr*rr + qq*qq*qq
   crad=sqrt(crad)

   cb=rr-crad

   if(cb.eq.czero)then
!        insoluble particle
         cx(1)=(-p1)**third
         cx(2)=cx(1)
         cx(3)=cx(1)
         cx(4)=cx(1)
   else
         cb = cb**third
         cy = -cb + qq/cb + p2/6._r8

         cb0 = sqrt(cy*cy - p0)
         cb1 = (p3*cy - p1)/(2._r8*cb0)

         cb = p3/2._r8 + cb1
         crad = cb*cb - 4._r8*(cy+cb0)
         crad = sqrt(crad)
         cx(1) = (-cb+crad)/2._r8
         cx(2) = (-cb-crad)/2._r8

         cb = p3/2._r8-cb1
         crad = cb*cb - 4._r8*(cy-cb0)
         crad = sqrt(crad)
         cx(3) = (-cb+crad)/2._r8
         cx(4) = (-cb-crad)/2._r8
   endif

   return
end subroutine makoh_quartic

!===============================================================================
subroutine find_real_solution(               &
             rdry,   cx,                     & ! in
             rwet,   nsol                    ) ! out
!----------------------------------------------------------------------
!  find the smallest real solution from the polynomial solver
!----------------------------------------------------------------------
!
   real(r8),    intent(in) :: rdry      ! dry radius
   complex(r8), intent(in) :: cx(:)     ! polynomial output
   real(r8),    intent(out) :: rwet     ! wet radius
   integer,     intent(out) :: nsol     ! index of the solution

   ! local variables
   integer  :: n_cx     ! length of cx
   integer  :: nn       ! index of cx
   real(r8) :: xr, xi   ! real and image part of cx
   real(r8), parameter :: eps = 1.e-4_r8

   n_cx = size(cx, dim=1)

   rwet = 1000._r8*rdry
   nsol = 0
   do nn=1,n_cx
       xr = real(cx(nn))
       xi = aimag(cx(nn))
       if(abs(xi).gt.abs(xr)*eps) cycle
       if(xr.gt.rwet) cycle
       if(xr.lt.rdry*(1._r8-eps)) cycle
       if(xr.ne.xr) cycle
       rwet = xr
       nsol = nn
   enddo

end subroutine find_real_solution

!===============================================================================
subroutine modal_aero_wateruptake_wetdens( ncol,             & ! in
             wetvol, wtrvol, drymass,      specdens_1,       & ! in
             wetdens                                         ) ! inout
!-----------------------------------------------------------------------
! compute aerosol wet density
!-----------------------------------------------------------------------
   integer, intent(in)  :: ncol                 ! number of columns
   real(r8),intent(in)  :: wetvol(:,:,:)        ! single-particle-mean wet aerosol volume [m^3]
   real(r8),intent(in)  :: wtrvol(:,:,:)        ! single-particle-mean water volume in wet aerosol [m^3]
   real(r8),intent(in)  :: drymass(:,:,:)       ! single-particle-mean dry mass [kg]
   real(r8),intent(in)  :: specdens_1(:)        ! specified dry aerosol density [kg/m3]
   real(r8),intent(inout) :: wetdens(:,:,:)     ! wet aerosol density [kg/m3]
   ! local variables
   integer :: icol, kk, imode
   real(r8), parameter  :: small_value = 1.0e-30_r8

   ! compute aerosol wet density (kg/m3)
   do imode = 1, nmodes
      do kk = top_lev, pver
         do icol = 1, ncol
            if (wetvol(icol,kk,imode) > small_value) then
               ! wet density
               wetdens(icol,kk,imode) = (drymass(icol,kk,imode) + rhoh2o*wtrvol(icol,kk,imode)) &
                                        / wetvol(icol,kk,imode)
            else
               ! dry density
               wetdens(icol,kk,imode) = specdens_1(imode)
            endif

         enddo ! i = 1, ncol
      enddo ! k = top_lev, pver
   enddo ! m = 1, nmodes

end subroutine modal_aero_wateruptake_wetdens


!===============================================================================
end module modal_aero_wateruptake


