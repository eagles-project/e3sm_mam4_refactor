module modal_aero_wateruptake

!   RCE 07.04.13:  Adapted from MIRAGE2 code

use shr_kind_mod,     only: r8 => shr_kind_r8
use physconst,        only: pi, rhoh2o
use ppgrid,           only: pcols, pver
use physics_types,    only: physics_state
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field

use wv_saturation,    only: qsat_water
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_aer_mmr, rad_cnst_get_aer_props, &
                            rad_cnst_get_mode_props, rad_cnst_get_mode_num
use cam_history,      only: addfld, add_default, outfld
use cam_logfile,      only: iulog
use ref_pres,         only: top_lev => clim_modal_aero_top_lev
use phys_control,     only: phys_getopts
use cam_abortutils,   only: endrun
use mam_support,      only: min_max_bound


implicit none
private
save

public :: &
   modal_aero_wateruptake_init, &
   modal_aero_wateruptake_dr, &
   modal_aero_kohler


public :: modal_aero_wateruptake_reg

real(r8), parameter :: third = 1._r8/3._r8
real(r8), parameter :: pi43  = pi*4.0_r8/3.0_r8


! Physics buffer indices
integer :: cld_idx        = 0
integer :: dgnum_idx      = 0
integer :: dgnumwet_idx   = 0
integer :: wetdens_ap_idx = 0
integer :: qaerwat_idx    = 0

logical :: pergro_mods         = .false.

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
  use rad_constituents, only: rad_cnst_get_info

   integer :: nmodes
   
   call rad_cnst_get_info(0, nmodes=nmodes)
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

   integer :: m, nmodes, istat
   logical :: history_aerosol      ! Output the MAM aerosol variables and tendencies
   logical :: history_verbose      ! produce verbose history output

   character(len=3) :: trnum       ! used to hold mode number (as characters)
   !----------------------------------------------------------------------------

   cld_idx        = pbuf_get_index('CLD')    
   dgnum_idx      = pbuf_get_index('DGNUM')    

   ! assume for now that will compute wateruptake for climate list modes only

   call rad_cnst_get_info(0, nmodes=nmodes)

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

   do m = 1, nmodes
      write(trnum, '(i3.3)') m
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

   end do

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

subroutine modal_aero_wateruptake_dr(state, pbuf, list_idx_in, dgnumdry_m, dgnumwet_m, &
                                     qaerwat_m, wetdens_m, clear_rh_in)

  use shr_log_mod ,   only: errmsg => shr_log_errmsg
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
   real(r8), optional, allocatable, target, intent(inout)   :: wetdens_m(:,:,:)
   ! optional input relative humidty (overrides clearsky RH estimate below)
   real(r8), optional,          intent(in)    :: clear_rh_in(pcols,pver)

   !----------------------------------------------------------------------------
   ! local variables

   integer  :: lchnk              ! chunk index
   integer  :: ncol               ! number of columns
   integer  :: list_idx           ! radiative constituents list index
   integer  :: stat

   integer :: i, k, l, m
   integer :: itim_old
   integer :: nmodes
   integer :: nspec

   real(r8), pointer :: h2ommr(:,:) ! specific humidity
   real(r8), pointer :: t(:,:)      ! temperatures (K)
   real(r8), pointer :: pmid(:,:)   ! layer pressure (Pa)
   real(r8), pointer :: raer(:,:)   ! aerosol species MRs (kg/kg and #/kg)

   real(r8), pointer :: cldn(:,:)      ! layer cloud fraction (0-1)
   real(r8), pointer :: dgncur_a(:,:,:)
   real(r8), pointer :: dgncur_awet(:,:,:)
   real(r8), pointer :: wetdens(:,:,:)
   real(r8), pointer :: qaerwat(:,:,:)

   real(r8) :: dryvolmr(pcols,pver)          ! volume MR for aerosol mode (m3/kg)
   real(r8) :: specdens
   real(r8) :: spechygro, spechygro_1
   real(r8) :: duma, dumb
   real(r8) :: sigmag
   real(r8) :: alnsg
   real(r8) :: v2ncur_a
   real(r8) :: drydens               ! dry particle density  (kg/m^3)
   real(r8) :: rh(pcols,pver)        ! relative humidity (0-1)

   real(r8) :: es(pcols)             ! saturation vapor pressure
   real(r8) :: qs(pcols)             ! saturation specific humidity
   real(r8) :: cldn_thresh
   real(r8) :: aerosol_water(pcols,pver) !sum of aerosol water (wat_a1 + wat_a2 + wat_a3 + wat_a4)
   logical :: history_aerosol      ! Output the MAM aerosol variables and tendencies
   logical :: history_verbose      ! produce verbose history output
   logical :: compute_wetdens

   character(len=3) :: trnum       ! used to hold mode number (as characters)
   !----------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol = state%ncol

   ! determine default variables
   call phys_getopts(history_aerosol_out = history_aerosol, &
                     history_verbose_out = history_verbose)

   list_idx = 0
   if (present(list_idx_in)) then
      list_idx = list_idx_in

      ! check that all optional args for diagnostic mode are present
      if (.not. present(dgnumdry_m) .or. .not. present(dgnumwet_m) .or. &
          .not. present(qaerwat_m)) then
         call endrun('modal_aero_wateruptake_dr called '// &
              'with list_idx_in but required args not present '//errmsg(__FILE__,__LINE__))
      end if

      ! arrays for diagnostic calculations must be allocated
      if (.not. allocated(dgnumdry_m) .or. .not. allocated(dgnumwet_m) .or. &
          .not. allocated(qaerwat_m)) then
         call endrun('modal_aero_wateruptake_dr called '// &
              'with list_idx_in but required args not allocated '//errmsg(__FILE__,__LINE__))
      end if
   end if ! if present(list_idx_in)

   ! loop over all aerosol modes
   call rad_cnst_get_info(list_idx, nmodes=nmodes)

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

   maer(:,:,:)     = 0._r8
   hygro(:,:,:)    = 0._r8

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
      if(present(wetdens_m)) then
         if (.not. allocated(wetdens_m)) then
            call endrun('modal_aero_wateruptake_dr called '// &
                 'with list_idx_in but wetdens_m is not allocated '//errmsg(__FILE__,__LINE__))
         endif
         wetdens     => wetdens_m
      else
         !set compute_wetdens to flase if wetdens is not present
         compute_wetdens = .false.
      endif
   end if

   !----------------------------------------------------------------------------
   ! retreive aerosol properties

   do m = 1, nmodes

      dryvolmr(:,:) = 0._r8

      ! get mode properties
      call rad_cnst_get_mode_props(list_idx, m, sigmag=sigmag,  &
         rhcrystal=rhcrystal(m), rhdeliques=rhdeliques(m))

      ! get mode info
      call rad_cnst_get_info(list_idx, m, nspec=nspec)

      do l = 1, nspec

         ! get species interstitial mixing ratio ('a')
         call rad_cnst_get_aer_mmr(list_idx, m, l, 'a', state, pbuf, raer)
         call rad_cnst_get_aer_props(list_idx, m, l, density_aer=specdens, hygro_aer=spechygro)

         if (l == 1) then
            ! save off these values to be used as defaults
            specdens_1(m)  = specdens
            spechygro_1    = spechygro
         end if

         do k = top_lev, pver
            do i = 1, ncol
               duma          = raer(i,k)
               maer(i,k,m)   = maer(i,k,m) + duma
               dumb          = duma/specdens
               dryvolmr(i,k) = dryvolmr(i,k) + dumb
               hygro(i,k,m)  = hygro(i,k,m) + dumb*spechygro
            end do ! i = 1, ncol
         end do ! k = top_lev, pver
      end do ! l = 1, nspec

      alnsg = log(sigmag)

      do k = top_lev, pver
         do i = 1, ncol

            if (dryvolmr(i,k) > 1.0e-30_r8) then
               hygro(i,k,m) = hygro(i,k,m)/dryvolmr(i,k)
            else
               hygro(i,k,m) = spechygro_1
            end if

            ! dry aerosol properties

            v2ncur_a = 1._r8 / ( (pi/6._r8)*(dgncur_a(i,k,m)**3._r8)*exp(4.5_r8*alnsg**2._r8) )
            ! naer = aerosol number (#/kg)
            naer(i,k,m) = dryvolmr(i,k)*v2ncur_a

            ! compute mean (1 particle) dry volume and mass for each mode
            ! old coding is replaced because the new (1/v2ncur_a) is equal to
            ! the mean particle volume
            ! also moletomass forces maer >= 1.0e-30, so (maer/dryvolmr)
            ! should never cause problems (but check for maer < 1.0e-31 anyway)
            if (maer(i,k,m) .gt. 1.0e-31_r8) then
               drydens = maer(i,k,m)/dryvolmr(i,k)
            else
               drydens = 1.0_r8
            end if
            dryvol(i,k,m)   = 1.0_r8/v2ncur_a
            drymass(i,k,m)  = drydens*dryvol(i,k,m)
            dryrad(i,k,m)   = (dryvol(i,k,m)/pi43)**third

         end do ! i = 1, ncol
      end do ! k = top_lev, pver

   end do    ! modes

   !----------------------------------------------------------------------------
   ! specify clear air relative humidity

   if (present(clear_rh_in)) then

      ! use input relative humidity
      rh(1:ncol,1:pver) = clear_rh_in(1:ncol,1:pver)

      ! check that values are reasonable and apply upper limit
      do k = top_lev, pver
         do i = 1, ncol
            if ( rh(i,k)<0 ) then
               write(iulog,*) 'modal_aero_wateruptake_dr: clear_rh_in is negative - rh:',rh(i,k),' k=',k
               call endrun('modal_aero_wateruptake_dr: clear_rh_in cannot be negative')
            end if
            ! limit RH to 98% to be consistent with behavior when clear_rh_in is not provided
            rh(i,k) = min(rh(i,k), 0.98_r8)
         end do ! i
      end do ! k

   else

      ! estimate clear air relative humidity using cloud fraction
      h2ommr => state%q(:,:,1)
      t      => state%t
      pmid   => state%pmid

      itim_old    =  pbuf_old_tim_idx()
      call pbuf_get_field(pbuf, cld_idx, cldn, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

      do k = top_lev, pver
         call qsat_water(t(:ncol,k), pmid(:ncol,k), es(:ncol), qs(:ncol))
         do i = 1, ncol
            if (qs(i) > h2ommr(i,k)) then
               rh(i,k) = h2ommr(i,k)/qs(i)
            else
               rh(i,k) = 0.98_r8
            endif
            rh(i,k) = max(rh(i,k), 0.0_r8)
            rh(i,k) = min(rh(i,k), 0.98_r8)
            if(pergro_mods) then
               cldn_thresh = 0.9998_r8
            else
               cldn_thresh = 1.0_r8 !original code
            endif
            if (cldn(i,k) .lt. cldn_thresh) then
               rh(i,k) = (rh(i,k) - cldn(i,k)) / (1.0_r8 - cldn(i,k))  ! RH of clear portion
            end if
            rh(i,k) = max(rh(i,k), 0.0_r8)
         end do ! i = 1, ncol
      end do ! k = top_lev, pver

   end if ! if present(clear_rh_in)

   !----------------------------------------------------------------------------
   ! compute aerosol wet radius and aerosol water

   call modal_aero_wateruptake_sub( &
      ncol, nmodes, rhcrystal, rhdeliques, dryrad, &
      hygro, rh, dryvol, wetrad, wetvol,           &
      wtrvol)

   do m = 1, nmodes

      do k = top_lev, pver
         do i = 1, ncol

            dgncur_awet(i,k,m) = dgncur_a(i,k,m) * (wetrad(i,k,m)/dryrad(i,k,m))
            qaerwat(i,k,m)     = rhoh2o*naer(i,k,m)*wtrvol(i,k,m)

            ! compute aerosol wet density (kg/m3)
            if(compute_wetdens) then
               if (wetvol(i,k,m) > 1.0e-30_r8) then
                  wetdens(i,k,m) = (drymass(i,k,m) + rhoh2o*wtrvol(i,k,m))/wetvol(i,k,m)
               else
                  wetdens(i,k,m) = specdens_1(m)
               end if
            endif
         end do ! i = 1, ncol
      end do ! k = top_lev, pver

   end do ! m = 1, nmodes

   !----------------------------------------------------------------------------
   ! write history output if not in diagnostic mode

   if (.not.present(list_idx_in)) then

      aerosol_water(:ncol,:) = 0._r8
      do m = 1, nmodes
         ! output to history
         write( trnum, '(i3.3)' ) m
         call outfld( 'wat_a'//trnum(3:3),  qaerwat(:,:,m),     pcols, lchnk)
         call outfld( 'dgnd_a'//trnum(2:3), dgncur_a(:,:,m),    pcols, lchnk)
         call outfld( 'dgnw_a'//trnum(2:3), dgncur_awet(:,:,m), pcols, lchnk)
         if (history_aerosol .and. .not. history_verbose) &
         aerosol_water(:ncol,:) = aerosol_water(:ncol,:) + qaerwat(:ncol,:,m)
      end do ! m = 1, nmodes

      if (history_aerosol .and. .not. history_verbose) &
         call outfld( 'aero_water',  aerosol_water(:ncol,:),    ncol, lchnk)

   end if

end subroutine modal_aero_wateruptake_dr

!===============================================================================

   subroutine modal_aero_wateruptake_sub(               &
           ncol,    nmodes, rhcrystal, rhdeliques,      & ! in
           dryrad,  hygro,  rh,        dryvol,          & ! in
           wetrad,  wetvol, wtrvol                      ) ! out

!-----------------------------------------------------------------------
!
! Purpose: Compute aerosol wet radius
!
! Method:  Kohler theory
!
! Author:  S. Ghan
!
!-----------------------------------------------------------------------

   ! Arguments
   integer, intent(in)  :: ncol                    ! number of columns
   integer, intent(in)  :: nmodes

   real(r8), intent(in) :: rhcrystal(:)          ! crystal RH [fraction]
   real(r8), intent(in) :: rhdeliques(:)         ! deliques RH [Fraction]
   real(r8), intent(in) :: dryrad(:,:,:)         ! dry volume mean radius of aerosol [m]
   real(r8), intent(in) :: hygro(:,:,:)          ! volume-weighted mean hygroscopicity [unitless]
   real(r8), intent(in) :: rh(:,:)               ! relative humidity [fraction]
   real(r8), intent(in) :: dryvol(:,:,:)

   real(r8), intent(out) :: wetrad(:,:,:)        ! wet radius of aerosol [m]
   real(r8), intent(out) :: wetvol(:,:,:)        ! single-particle-mean wet volume [m^3]
   real(r8), intent(out) :: wtrvol(:,:,:)        ! single-particle-mean water volume in wet aerosol [m^3]

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

         enddo  ! columns
      enddo     ! levels

   enddo ! modes

   end subroutine modal_aero_wateruptake_sub

!-----------------------------------------------------------------------
   subroutine modal_aero_kohler( rdry_in, hygro, rh,    & ! in
                                 rwet_out               ) ! out

! calculates equlibrium radius r of haze droplets as function of
! dry particle mass and relative humidity s using kohler solution
! given in pruppacher and klett (eqn 6-35)

! for multiple aerosol types, assumes an internal mixture of aerosols

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
   real(r8) :: p30,p31,p32 ! coefficients of polynomial
   real(r8) :: r3, r4      ! tempararily store wet radius from polynomials [microns]
   real(r8) :: rwet        ! wet radius [microns]
   real(r8) :: rdry        ! radius of dry particle [microns]
   real(r8) :: ss          ! relative humidity (1 = saturated) [fraction]
   real(r8) :: slog        ! log relative humidity
   real(r8) :: vol         ! total volume of particle [microns**3]

   complex(r8) :: cx4(4),cx3(3) ! output of polynomials

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
   real(r8), parameter :: small_vol = 1.e-12_r8   ! small value to avoid zero
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
!  cubic for rh=1
   p32 = 0._r8
   p31 = -bb/aa
   p30 = -vol

   if(vol .le. small_vol) then
         rwet=rdry
         rwet_out = rwet*factor_um2m
         return
   endif

   pp = abs(p31)/(rdry*rdry)
   if(pp.lt.eps)then
!      approximate solution for small particles
       rwet = rdry * (1._r8 + pp*third/(1._r8-slog*rdry/aa))
   else
       call makoh_quartic(cx4,p43,p42,p41,p40)
!      find smallest real(r8) solution
       call find_real_solution( rdry, cx4, rwet,nsol)
       if(nsol.eq.0)then
              write(iulog,*)   &
               'ccm kohlerc - no real(r8) solution found (quartic)'
              write(iulog,*)'roots =', (cx4(nn),nn=1,4)
              write(iulog,*)'p0-p3 =', p40, p41, p42, p43
              write(iulog,*)'rh=',rh
              write(iulog,*)'setting radius to dry radius=',rdry
              rwet = rdry
       endif
   endif

   if(rh.gt.1._r8-eps)then
!      save quartic solution at s=1-eps
       r4=rwet
!      cubic for rh=1
       pp = abs(p31)/(rdry*rdry)
       if(pp.lt.eps)then
            rwet = rdry*(1._r8+pp*third)
       else
            call makoh_cubic(cx3,p32,p31,p30)
!           find smallest real(r8) solution
            call find_real_solution( rdry, cx3, rwet,nsol) 
            if(nsol.eq.0)then
                 write(iulog,*)   &
                  'ccm kohlerc - no real(r8) solution found (cubic)'
                 write(iulog,*)'roots =', (cx3(nn),nn=1,3)
                 write(iulog,*)'p0-p2 =', p30, p31, p32
                 write(iulog,*)'rh=',rh
                 write(iulog,*)'setting radius to dry radius=',rdry
                 rwet=rdry
            endif
       endif
       r3=rwet
!      now interpolate between quartic, cubic solutions
       rwet = (r4*(1._r8-rh) + r3*(rh-1._r8+eps))/eps
   endif

! bound and convert from microns to m
   rwet = min(rwet, rmax) ! upper bound based on 1 day lifetime
   rwet_out = rwet*factor_um2m

   return
   end subroutine modal_aero_kohler


!-----------------------------------------------------------------------
   subroutine makoh_cubic( cx, p2, p1, p0)
!
!     solves  x**3 + p2 x**2 + p1 x + p0 = 0
!     where p0, p1, p2 are real
!
   real(r8) :: p0, p1, p2   ! input parameters
   complex(r8) :: cx(3)     ! output

   real(r8) :: qq, rr, sqrt3, third  ! temparary variables
   complex(r8) :: ci, cq, crad, cw, cwsq, cy, cz  ! temparary variables
   real(r8), parameter :: eps=1.e-20_r8


   third=1._r8/3._r8
   ci=cmplx(0._r8,1._r8,r8)
   sqrt3=sqrt(3._r8)
   cw=0.5_r8*(-1._r8 + ci*sqrt3)
   cwsq=0.5_r8*(-1._r8 - ci*sqrt3)

   if(p1.eq.0._r8)then
!        completely insoluble particle
         cx(1) = (-p0)**third
         cx(2) = cx(1)
         cx(3) = cx(1)
   else
         qq = p1/3._r8
         rr = p0/2._r8
         crad = rr*rr + qq*qq*qq
         crad = sqrt(crad)

         cy = rr-crad
         if (abs(cy).gt.eps) cy=cy**third
         cq = qq
         cz = -cq/cy

         cx(1) = -cy-cz
         cx(2) = -cw*cy - cwsq*cz
         cx(3) = -cwsq*cy - cw*cz
   endif

   return
   end subroutine makoh_cubic


!-----------------------------------------------------------------------
   subroutine makoh_quartic( cx, p3, p2, p1, p0)

!     solves x**4 + p3 x**3 + p2 x**2 + p1 x + p0 = 0
!     where p0, p1, p2, p3 are real
!
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

!----------------------------------------------------------------------
   subroutine find_real_solution(               &
                rdry,   cx,                     & ! in
                rwet,   nsol                    ) ! out
!
!  find the smallest real solution from the polynomial solver
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

!----------------------------------------------------------------------

   end module modal_aero_wateruptake


