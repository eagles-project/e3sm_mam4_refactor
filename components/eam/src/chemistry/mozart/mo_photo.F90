module mo_photo
  !----------------------------------------------------------------------
  !	... photolysis interp table and related arrays
  !----------------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
  use ppgrid,       only : pcols, pver, begchunk, endchunk
  use cam_abortutils,   only : endrun
  use pio
  use cam_pio_utils,only : cam_pio_openfile
  use spmd_utils,   only : masterproc
  use cam_logfile,  only : iulog

  implicit none

  private

  public :: photo_inti, table_photo
  public :: set_ub_col
  public :: setcol 
  public :: photo_timestep_init

  save

  integer, parameter  :: pverm = pver - 1
!  Need a larger maximum zenith angle for WACCM-X extended to high altitudes
  real(r8), parameter :: max_zen_angle = 88.85_r8  !  degrees

  integer ::  jonitr_ndx
  integer ::  jho2no2_ndx
  integer ::  o3_ndx, o3_inv_ndx
  integer ::  o2_ndx
  integer, allocatable :: lng_indexer(:)
  integer, allocatable :: sht_indexer(:)
  integer, allocatable :: euv_indexer(:)

  integer              :: ki
  integer              :: last
  integer              :: next
  integer              :: n_exo_levs
  real(r8)                 :: delp
  real(r8)                 :: dels
  real(r8), allocatable    :: days(:)
  real(r8), allocatable    :: levs(:)
  real(r8), allocatable    :: o2_exo_coldens(:,:,:,:)
  real(r8), allocatable    :: o3_exo_coldens(:,:,:,:)
  logical              :: o2_is_inv
  logical              :: o3_is_inv
  logical              :: has_o2_col
  logical              :: has_o3_col
  logical              :: has_fixed_press


contains

!====================================================================================
  subroutine photo_inti( xs_coef_file, xs_short_file, xs_long_file, rsf_file, &
       euvacdat_file, photon_file, electron_file, &
       exo_coldens_file, tuv_xsect_file, o2_xsect_file, xactive_prates, is_waccm )
    !----------------------------------------------------------------------
    !	... initialize photolysis module
    !----------------------------------------------------------------------

    use mo_photoin,    only : photoin_inti
    use mo_tuv_inti,   only : tuv_inti
    use mo_tuv_inti,   only : nlng
    use mo_seto2,      only : o2_xsect_inti      
    use interpolate_data, only: lininterp_init, lininterp, lininterp_finish, interp_type
    use chem_mods,     only : phtcnt
    use chem_mods,     only : nabscol
    use chem_mods,     only : rxt_tag_lst, pht_alias_lst, pht_alias_mult
    use time_manager,  only : get_calday
    use ioFileMod,     only : getfil
    use mo_chem_utls,  only : get_spc_ndx, get_rxt_ndx, get_inv_ndx
    use mo_jlong,      only : jlong_init
    use mo_constants, only : d2r
    use ref_pres,     only : num_pr_lev, ptop_ref
    use seasalt_model, only : sslt_names=>seasalt_names, sslt_ncnst=>seasalt_nbin
    use mo_jeuv,       only : jeuv_init
    use dyn_grid,      only : get_dyn_grid_parm
    use phys_grid,     only : get_ncols_p, get_rlat_all_p    
    use solar_data,    only : has_spectrum

    implicit none

    !----------------------------------------------------------------------
    !	... dummy arguments
    !----------------------------------------------------------------------
    character(len=*), intent(in) :: xs_long_file, rsf_file
    character(len=*), intent(in) :: exo_coldens_file
    character(len=*), intent(in) :: tuv_xsect_file
    character(len=*), intent(in) :: o2_xsect_file
    logical, intent(in)          :: xactive_prates
    ! waccm 
    character(len=*), intent(in) :: xs_coef_file
    character(len=*), intent(in) :: xs_short_file
    character(len=*), intent(in) :: euvacdat_file
    character(len=*), intent(in) :: photon_file
    character(len=*), intent(in) :: electron_file

    logical, optional, intent(in) :: is_waccm

    !----------------------------------------------------------------------
    !	... local variables
    !----------------------------------------------------------------------
    real(r8), parameter   :: hPa2Pa = 100._r8
    integer           :: k, n
    type(file_desc_t) :: ncid
    type(var_desc_t)  :: vidO2, vidO3, vid
    type(interp_type) :: lat_wgts
    integer           :: dimid
    integer           :: nlat
    integer           :: ntimes
    integer           :: astat
    integer           :: gndx
    integer           :: ndx
    integer           :: spc_ndx
    integer           :: ierr
    integer           :: c, ncols
    integer           :: plev, plevp
    integer, allocatable :: dates(:)
    real(r8)              :: pinterp
    real(r8), allocatable :: lats(:)
    real(r8), allocatable :: coldens(:,:,:)
    character(len=256)    :: locfn
    character(len=256)    :: filespec
    real(r8), parameter :: trop_thrshld = 1._r8 ! Pa
    real(r8) :: to_lats(pcols)


    if( phtcnt < 1 ) then
       return
    end if

    if (present(is_waccm) .and. is_waccm==.true. ) then
       call endrun('FORTRAN refactoring: code related to waccm is removed in MAMxx. is_waccm=true is currently not supported') 
    endif
    if (xactive_prates) call endrun('FORTRAN refactoring: xactive_prates=true is removed in MAMxx')


    ! jeuv_1,,, jeuv_25 --> need euv calculations  
    ! how to determine if shrt calc is needed ?? -- use top level pressure 

    if ( .not. has_spectrum ) then
       write(iulog,*) 'photo_inti: solar_data file needs to contain irradiance spectrum'
       call endrun('photo_inti: ERROR -- solar irradiance spectrum is missing')
    endif
    
    plev = get_dyn_grid_parm('plev')
    plevp = get_dyn_grid_parm('plevp')

    !----------------------------------------------------------------------
    !	... allocate indexers
    !----------------------------------------------------------------------
    allocate( lng_indexer(phtcnt),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'photo_inti: lng_indexer allocation error = ',astat
       call endrun
    end if
    lng_indexer(:) = 0
    allocate( sht_indexer(phtcnt),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'photo_inti: Failed to allocate sht_indexer; error = ',astat
       call endrun
    end if
    sht_indexer(:) = 0
    allocate( euv_indexer(phtcnt),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'photo_inti: Failed to allocate euv_indexer; error = ',astat
       call endrun
    end if
    euv_indexer(:) = 0


    o3_ndx     = get_spc_ndx( 'O3' )
    o3_inv_ndx = get_inv_ndx( 'O3' )

    o2_ndx     = get_inv_ndx( 'O2' )
    o2_is_inv  = o2_ndx > 0
    if( .not. o2_is_inv ) then
       o2_ndx = get_spc_ndx( 'O2' )
    end if
    o3_is_inv  = o3_ndx < 1


    !----------------------------------------------------------------------
    !	... call module initializers
    !----------------------------------------------------------------------
    call jlong_init( xs_long_file, rsf_file, lng_indexer )
    jho2no2_ndx = get_rxt_ndx( 'jho2no2_b' )

    !----------------------------------------------------------------------
    !        ... check that each photorate is in short or long datasets
    !----------------------------------------------------------------------
    if( any( (abs(sht_indexer(:)) + abs(lng_indexer(:)) + abs(euv_indexer(:))) == 0 ) ) then
       write(iulog,*) ' '
       write(iulog,*) 'photo_inti: the following photorate(s) are not in'
       write(iulog,*) '            either the short or long datasets'
       write(iulog,*) ' '
       do ndx = 1,phtcnt
          if( abs(sht_indexer(ndx)) + abs(lng_indexer(ndx)) == 0 ) then
             write(iulog,*) '           ',trim( rxt_tag_lst(ndx) )
          end if
       end do
       call endrun
    end if

    if( masterproc ) then
       write(iulog,*) ' '
       write(iulog,*) '*********************************************'
       write(iulog,*) 'photo_inti: euv_indexer'
       write(iulog,'(10i6)') euv_indexer(:)
       write(iulog,*) 'photo_inti: sht_indexer'
       write(iulog,'(10i6)') sht_indexer(:)
       write(iulog,*) 'photo_inti: lng_indexer'
       write(iulog,'(10i6)') lng_indexer(:)
       write(iulog,*) '*********************************************'
       write(iulog,*) ' '
    endif

    !----------------------------------------------------------------------
    !	... check for o2, o3 absorber columns
    !----------------------------------------------------------------------
    if( nabscol > 0 ) then
       spc_ndx = o3_ndx
       if( spc_ndx > 0 ) then
          has_o3_col = .true.
       else
          has_o3_col = .false.
       end if
       if( nabscol > 1 ) then
          if( o2_ndx > 1 ) then
             has_o2_col = .true.
          else
             has_o2_col = .false.
          end if
       else
          has_o2_col = .false.
       end if
    else
       has_o2_col = .false.
       has_o3_col = .false.
    end if

    if ( len_trim(exo_coldens_file) == 0 ) then
       has_o2_col = .false.
       has_o3_col = .false.
    endif

    has_abs_columns : if( has_o2_col .or. has_o3_col ) then
       !-----------------------------------------------------------------------
       !	... open exo coldens file
       !-----------------------------------------------------------------------
       filespec = trim( exo_coldens_file )
       call getfil( filespec, locfn, 0 )
       call cam_pio_openfile( ncid, trim(locfn), PIO_NOWRITE )

       !-----------------------------------------------------------------------
       !       ... get grid dimensions from file
       !-----------------------------------------------------------------------
       !       ... timing
       !-----------------------------------------------------------------------
       ierr = pio_inq_dimid( ncid, 'month', dimid )
       ierr = pio_inq_dimlen( ncid, dimid, ntimes )

       if( ntimes /= 12 ) then
          call endrun('photo_inti: exo coldens is not annual period')
       end if
       allocate( dates(ntimes),days(ntimes),stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'photo_inti: dates,days allocation error = ',astat
          call endrun
       end if
       dates(:) = (/ 116, 214, 316, 415,  516,  615, &
            716, 816, 915, 1016, 1115, 1216 /)
       !-----------------------------------------------------------------------
       !	... initialize the monthly day of year times
       !-----------------------------------------------------------------------
       do n = 1,ntimes
          days(n) = get_calday( dates(n), 0 )
       end do
       deallocate( dates )
       !-----------------------------------------------------------------------
       !       ... latitudes
       !-----------------------------------------------------------------------
       ierr = pio_inq_dimid( ncid, 'lat', dimid )
       ierr = pio_inq_dimlen( ncid, dimid, nlat )
       allocate( lats(nlat), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'photo_inti: lats allocation error = ',astat
          call endrun
       end if
       ierr = pio_inq_varid( ncid, 'lat', vid )
       ierr = pio_get_var( ncid, vid, lats )
       lats(:nlat) = lats(:nlat) * d2r
       !-----------------------------------------------------------------------
       !       ... levels
       !-----------------------------------------------------------------------
       ierr = pio_inq_dimid( ncid, 'lev', dimid )
       ierr = pio_inq_dimlen( ncid, dimid, n_exo_levs )
       allocate( levs(n_exo_levs), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'photo_inti: levs allocation error = ',astat
          call endrun
       end if
       ierr = pio_inq_varid( ncid, 'lev', vid )
       ierr = pio_get_var( ncid, vid, levs )
       levs(:n_exo_levs) = levs(:n_exo_levs) * hPa2Pa
       !-----------------------------------------------------------------------
       !       ... set up regridding
       !-----------------------------------------------------------------------

       allocate( coldens(nlat,n_exo_levs,ntimes),stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'photo_inti: coldens allocation error = ',astat
          call endrun
       end if
       if( has_o2_col ) then
          allocate( o2_exo_coldens(n_exo_levs,pcols,begchunk:endchunk,ntimes),stat=astat )
          if( astat /= 0 ) then
             write(iulog,*) 'photo_inti: o2_exo_coldens allocation error = ',astat
             call endrun
          end if
          ierr = pio_inq_varid( ncid, 'O2_column_density', vid )
          ierr = pio_get_var( ncid, vid,coldens )

          do c=begchunk,endchunk
             ncols = get_ncols_p(c)
             call get_rlat_all_p(c, pcols, to_lats)
             call lininterp_init(lats, nlat, to_lats, ncols, 1, lat_wgts)
             do n=1,ntimes
                do k = 1,n_exo_levs
                   call lininterp(coldens(:,k,n), nlat, o2_exo_coldens(k,:,c,n), ncols, lat_wgts)
                end do
             end do
             call lininterp_finish(lat_wgts)
          enddo


       end if
       if( has_o3_col ) then
          allocate( o3_exo_coldens(n_exo_levs,pcols,begchunk:endchunk,ntimes),stat=astat )
          if( astat /= 0 ) then
             write(iulog,*) 'photo_inti: o3_exo_coldens allocation error = ',astat
             call endrun
          end if
          ierr = pio_inq_varid( ncid, 'O3_column_density', vid )
          ierr = pio_get_var( ncid, vid,coldens )

          do c=begchunk,endchunk
             ncols = get_ncols_p(c)
             call get_rlat_all_p(c, pcols, to_lats)
             call lininterp_init(lats, nlat, to_lats, ncols, 1, lat_wgts)
             do n=1,ntimes
                do k = 1,n_exo_levs
                   call lininterp(coldens(:,k,n), nlat, o3_exo_coldens(k,:,c,n), ncols, lat_wgts)
                end do
             end do
             call lininterp_finish(lat_wgts)
          enddo
       end if
       call pio_closefile (ncid)
       deallocate( coldens,stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'photo_inti: failed to deallocate coldens; error = ',astat
          call endrun
       end if
       has_fixed_press = (num_pr_lev .ne. 0)
       !-----------------------------------------------------------------------
       !	... setup the pressure interpolation
       !-----------------------------------------------------------------------
       if( has_fixed_press ) then
          pinterp =  ptop_ref
          if( pinterp <= levs(1) ) then
             ki   = 1
             delp = 0._r8
          else
             do ki = 2,n_exo_levs
                if( pinterp <= levs(ki) ) then
                   delp = log( pinterp/levs(ki-1) )/log( levs(ki)/levs(ki-1) )
                   exit
                end if
             end do
          end if
#ifdef DEBUG
          write(iulog,*) '-----------------------------------'
          write(iulog,*) 'photo_inti: diagnostics'
          write(iulog,*) 'ki, delp = ',ki,delp
          write(iulog,*) 'pinterp,levs(ki-1:ki) = ',pinterp,levs(ki-1:ki)
          write(iulog,*) '-----------------------------------'
#endif
       end if
    end if has_abs_columns

  end subroutine photo_inti

!====================================================================================
  subroutine table_photo( photos, pmid, pdel, temper, zmid, zint, &
                          col_dens, zen_angle, srf_alb, lwc, clouds, &
                          esfact, vmr, invariants, ncol, lchnk, pbuf )
!-----------------------------------------------------------------
!   	... table photorates for wavelengths > 200nm
!-----------------------------------------------------------------

    use chem_mods,   only : nabscol, phtcnt, gas_pcnst, nfs
    use chem_mods,   only : pht_alias_mult, indexm
    use mo_jlong,    only : nlng => numj, jlong
    use mo_jeuv,     only : neuv, jeuv, nIonRates
    use mo_constants, only : r2d
    use physics_buffer, only : physics_buffer_desc, pbuf_set_field

    implicit none

!-----------------------------------------------------------------
!   	... dummy arguments
!-----------------------------------------------------------------
    integer,  intent(in)    :: lchnk
    integer,  intent(in)    :: ncol
    real(r8), intent(in)    :: esfact                       ! earth sun distance factor
    real(r8), intent(in)    :: vmr(ncol,pver,max(1,gas_pcnst)) ! vmr
    real(r8), intent(in)    :: col_dens(ncol,pver,nabscol) ! column densities (molecules/cm^2)
    real(r8), intent(in)    :: zen_angle(ncol)              ! solar zenith angle (radians)
    real(r8), intent(in)    :: srf_alb(pcols)               ! surface albedo
    real(r8), intent(in)    :: lwc(ncol,pver)               ! liquid water content (kg/kg)
    real(r8), intent(in)    :: clouds(ncol,pver)            ! cloud fraction
    real(r8), intent(in)    :: pmid(pcols,pver)             ! midpoint pressure (Pa)
    real(r8), intent(in)    :: pdel(pcols,pver)             ! pressure delta about midpoint (Pa)
    real(r8), intent(in)    :: temper(pcols,pver)           ! midpoint temperature (K)
    real(r8), intent(in)    :: zmid(ncol,pver)              ! midpoint height (km)
    real(r8), intent(in)    :: zint(ncol,pver)              ! interface height (km)
    real(r8), intent(in)    :: invariants(ncol,pver,max(1,nfs)) ! invariant densities (molecules/cm^3)
    real(r8), intent(inout) :: photos(ncol,pver,phtcnt)     ! photodissociation rates (1/s)
    type(physics_buffer_desc),pointer :: pbuf(:)

!-----------------------------------------------------------------
!    	... local variables
!-----------------------------------------------------------------
    real(r8), parameter :: Pa2mb         = 1.e-2_r8       ! pascals to mb

    integer ::  icol, kk, mm                 ! indicies
    integer ::  astat
    real(r8) ::  sza
    real(r8) ::  colo3(pver)               ! vertical o3 column density
    real(r8) ::  parg(pver)                ! vertical pressure array (hPa)

    real(r8) ::  eff_alb(pver)             ! effective albedo from cloud modifications
    real(r8) ::  cld_mult(pver)            ! clould multiplier
    real(r8), allocatable ::  lng_prates(:,:) ! photorates matrix (1/s)
    real(r8), allocatable :: tline(:)               ! vertical temperature array


    if( phtcnt < 1 ) then
       return
    endif

    allocate( tline(pver) )

!-----------------------------------------------------------------
!	... allocate long temp arrays
!-----------------------------------------------------------------
    if (nlng>0) then
       allocate( lng_prates(nlng,pver),stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'photo: Failed to allocate lng_prates; error = ',astat
          call endrun
       endif
    endif

!-----------------------------------------------------------------
!	... zero all photorates
!-----------------------------------------------------------------
    do mm = 1,max(1,phtcnt)
       do kk = 1,pver
          photos(:,kk,mm) = 0._r8
       enddo
    enddo

    col_loop : do icol = 1,ncol
       sza = zen_angle(icol)*r2d
       daylight : if( sza >= 0._r8 .and. sza < max_zen_angle ) then
          !-----------------------------------------------------------------
          !     ... compute eff_alb and cld_mult -- needs to be before jlong
          !-----------------------------------------------------------------
          call cloud_mod( zen_angle(icol), clouds(icol,:), lwc(icol,:), & ! in
                          pdel(icol,:), srf_alb(icol), & ! in
                          eff_alb, cld_mult ) ! out
          cld_mult(:) = esfact * cld_mult(:)

          !-----------------------------------------------------------------
          !	... long wave length component
          !-----------------------------------------------------------------
          parg(:)     = Pa2mb*pmid(icol,:)
          colo3(:)    = col_dens(icol,:,1)
          tline(1:pver) = temper(icol,:pver)
          call jlong( pver, sza, eff_alb, parg, tline, colo3, & ! in
                      lng_prates ) ! out

          do mm = 1,phtcnt
             if( lng_indexer(mm) > 0 ) then
                photos(icol,:,mm) = cld_mult(:) * (photos(icol,:,mm) + &
                               pht_alias_mult(mm,2)*lng_prates(lng_indexer(mm),:) ) 
             endif
          enddo

       endif daylight
    enddo col_loop

    if ( allocated(lng_prates) ) deallocate( lng_prates )
    if ( allocated(tline) )   deallocate( tline )

  end subroutine table_photo

!====================================================================================
  subroutine cloud_mod( zen_angle, clouds, lwc, delp, srf_alb, & ! in
                        eff_alb, cld_mult                      ) ! out
    !-----------------------------------------------------------------------
    ! 	... cloud alteration factors for photorates and albedo
    !-----------------------------------------------------------------------

    implicit none

    !-----------------------------------------------------------------------
    ! 	... dummy arguments
    !-----------------------------------------------------------------------
    real(r8), intent(in)    ::  zen_angle         ! zenith angle [deg]
    real(r8), intent(in)    ::  srf_alb           ! surface albedo [fraction]
    real(r8), intent(in)    ::  clouds(pver)       ! cloud fraction [fraction]
    real(r8), intent(in)    ::  lwc(pver)          ! liquid water content [kg/kg]
    real(r8), intent(in)    ::  delp(pver)         ! del press about midpoint [Pa]
    real(r8), intent(out)   ::  eff_alb(pver)      ! effective albedo [fraction]
    real(r8), intent(out)   ::  cld_mult(pver)     ! photolysis mult factor

    !-----------------------------------------------------------------------
    ! 	... local variables
    !-----------------------------------------------------------------------
    integer :: kk
    real(r8)    :: coschi  ! cos (solar zenith angle)
    real(r8)    :: del_lwp(pver)        ! liquid water path in each layer [g/m2]
    real(r8)    :: del_tau(pver)        ! cloud optical depth in each layer
    real(r8)    :: above_tau(pver)      ! cloud optical depth above this layer
    real(r8)    :: below_tau(pver)      ! cloud optical depth below this layer
    real(r8)    :: above_cld(pver)      ! cloud cover above this layer
    real(r8)    :: below_cld(pver)      ! cloud cover below this layer
    real(r8)    :: above_tra(pver)      ! transmission factor above this layer
    real(r8)    :: below_tra(pver)      ! transmission factor below this layer
    real(r8)    :: fac1(pver)           ! factor to calculate cld_mult
    real(r8)    :: fac2(pver)           ! factor to calculate cld_mult

    real(r8), parameter :: rgrav = 1._r8/9.80616_r8  ! 1/g [s^2/m]
    real(r8), parameter :: f_lwp2tau = .155_r8  ! factor converting LWP to tau [unknown source and unit]
    real(r8), parameter :: tau_min = 5._r8      ! tau threshold below which assign cloud as zero

    !---------------------------------------------------------
    !	... modify lwc for cloud fraction and form
    !	    liquid water path and tau for each layer
    !---------------------------------------------------------
    where( clouds(:) /= 0._r8 )
       del_lwp(:) = rgrav * lwc(:) * delp(:) * 1.e3_r8 / clouds(:) ! the unit is (likely) g/m^2
       del_tau(:) = del_lwp(:) * f_lwp2tau * clouds(:)**1.5_r8
    elsewhere
       del_lwp(:) = 0._r8
       del_tau(:) = 0._r8
    endwhere
    !---------------------------------------------------------
    !    	... form integrated tau and cloud cover from top down
    !---------------------------------------------------------
    above_tau(1) = 0._r8
    above_cld(1) = 0._r8
    do kk = 1,pverm
       above_tau(kk+1) = del_tau(kk) + above_tau(kk)
       above_cld(kk+1) = clouds(kk) * del_tau(kk) + above_cld(kk)
    enddo
    do kk = 2,pver
       if( above_tau(kk) /= 0._r8 ) then
          above_cld(kk) = above_cld(kk) / above_tau(kk)
       else
          above_cld(kk) = above_cld(kk-1)
       endif
    enddo
    !---------------------------------------------------------
    !    	... form integrated tau and cloud cover from bottom up
    !---------------------------------------------------------
    below_tau(pver) = 0._r8
    below_cld(pver) = 0._r8
    do kk = pverm,1,-1
       below_tau(kk) = del_tau(kk+1) + below_tau(kk+1)
       below_cld(kk) = clouds(kk+1) * del_tau(kk+1) + below_cld(kk+1)
    enddo
    do kk = pverm,1,-1
       if( below_tau(kk) /= 0._r8 ) then
          below_cld(kk) = below_cld(kk) / below_tau(kk)
       else
          below_cld(kk) = below_cld(kk+1)
       endif
    enddo
    !---------------------------------------------------------
    !	... modify above_tau and below_tau via jfm
    !---------------------------------------------------------
    where( above_cld(2:pver) /= 0._r8 )
       above_tau(2:pver) = above_tau(2:pver) / above_cld(2:pver)
    endwhere
    where( below_cld(:pverm) /= 0._r8 )
       below_tau(:pverm) = below_tau(:pverm) / below_cld(:pverm)
    endwhere
    where( above_tau(2:pver) < tau_min )
       above_cld(2:pver) = 0._r8
    endwhere
    where( below_tau(:pverm) < tau_min )
       below_cld(:pverm) = 0._r8
    endwhere
    !---------------------------------------------------------
    !	... form transmission factors
    !---------------------------------------------------------
    above_tra(:) = 11.905_r8 / (9.524_r8 + above_tau(:))
    below_tra(:) = 11.905_r8 / (9.524_r8 + below_tau(:))
    !---------------------------------------------------------
    !	... form effective albedo
    !---------------------------------------------------------
    where( below_cld(:) /= 0._r8 )
       eff_alb(:) = srf_alb + below_cld(:) * (1._r8 - below_tra(:)) * (1._r8 - srf_alb)
    elsewhere
       eff_alb(:) = srf_alb
    endwhere

    coschi = max( cos( zen_angle ), .5_r8 )
    where( del_lwp(:)*f_lwp2tau < tau_min )
       fac1(:) = 0._r8
    elsewhere
       fac1(:) = 1.4_r8 * coschi - 1._r8
    endwhere

    fac2(:)     = min( 0._r8,  (1.6_r8*coschi*above_tra(:))-1._r8 )
    cld_mult(:) = 1._r8 + fac1(:) * clouds(:) + fac2(:) * above_cld(:)
    cld_mult(:) = max( .05_r8, cld_mult(:) )

  end subroutine cloud_mod

!====================================================================================
  subroutine set_ub_col( col_delta,             & ! out
             vmr, invariants, pdel, ncol, lchnk ) ! in
    !---------------------------------------------------------------
    !        ... set the column densities at the upper boundary
    !---------------------------------------------------------------
    use chem_mods, only : nabscol, nfs, gas_pcnst

    implicit none

    !---------------------------------------------------------------
    !        ... dummy args
    !---------------------------------------------------------------
    integer,  intent(in)    ::  ncol                                   ! number of columns in current chunk
    integer,  intent(in)    ::  lchnk                                  ! latitude indicies in chunk
    real(r8), intent(in)    ::  vmr(ncol,pver,gas_pcnst)               ! xported species vmr [mol/mol]
    real(r8), intent(in)    ::  pdel(pcols,pver)                       ! pressure delta about midpoints [Pa]
    real(r8), intent(in)    ::  invariants(ncol,pver,nfs)              ! invariant densities [molecules/cm^3]
    real(r8), intent(out)   ::  col_delta(ncol,0:pver,max(1,nabscol))  ! o2,o3 col dens above model [1/cm^2]

    !---------------------------------------------------------------
    !        ... local variables
    !---------------------------------------------------------------
    real(r8)    :: o2_exo_col(ncol)
    real(r8)    :: o3_exo_col(ncol)
 
    !---------------------------------------------------------------
    !        ... assign column density at the upper boundary
    !            the first column is o3 and the second is o2.
    !---------------------------------------------------------------

    !	... set exo absorber columns
    o2_exo_col(:) = 0._r8
    o3_exo_col(:) = 0._r8
    has_abs_cols : if( has_o2_col .and. has_o3_col ) then
       if( has_fixed_press ) then

          if( has_o2_col ) then
              call calc_exo_col(ncol, o2_exo_coldens(:,:,lchnk,:), & ! in
                                o2_exo_col) ! out
          endif

          if( has_o3_col ) then
              call calc_exo_col(ncol, o3_exo_coldens(:,:,lchnk,:), & ! in
                                o3_exo_col) ! out
          endif

       endif   ! has_fixed_press
    endif has_abs_cols

    !---------------------------------------------------------------
    col_delta(:,:,:) = 0._r8

    if (o3_ndx > 0 .or. o3_inv_ndx > 0) then
        call calc_col_delta(   col_delta(:,0:,1),  & ! out
                (o3_inv_ndx>0),o3_ndx, o3_inv_ndx, & ! in
                o3_exo_col,    vmr,    invariants, & ! in
                pdel, ncol                         ) ! in
    endif

    if( nabscol > 1 .and. o2_ndx > 1 ) then
        call calc_col_delta( col_delta(:,0:,2), & ! out
             o2_is_inv,    o2_ndx, o2_ndx,      & ! in   o2_ndx is inv if o2_is_inv=.true.
             o2_exo_col,   vmr,    invariants,  & ! in
             pdel, ncol                         ) ! in
    endif

  end subroutine set_ub_col

!====================================================================================
  subroutine calc_exo_col ( ncol,  exo_coldens, & ! in
                            spc_exo_col         ) ! out
    !--------------------------------------------------------------------------------
    ! calculate exo absorber columns for o2 or o3
    !--------------------------------------------------------------------------------

    implicit none

    integer,  intent(in)  :: ncol         ! number of columns in current chunk
    real(r8), intent(in)  :: exo_coldens(:,:,:)     ! [molecules/cm^2]
    real(r8), intent(out) :: spc_exo_col(ncol)      ! exo absorber columns [molecules/cm^2]

    integer     :: icol
    integer     :: kl           ! ki - 1
    real(r8)    :: tint_vals(2) ! [molecules/cm^2]

    ! get ki-1
    kl = ki-1

    ! note that ki, last, next, delp, dels, are all from module variables
    do icol = 1,ncol
        if ( kl > 0 ) then
            tint_vals(1) = exo_coldens(kl,icol,last) &
                 + delp * (exo_coldens(ki,icol,last) &
                 - exo_coldens(kl,icol,last))
            tint_vals(2) = exo_coldens(kl,icol,next) &
                 + delp * (exo_coldens(ki,icol,next) &
                 - exo_coldens(kl,icol,next))
        else
            tint_vals(1) = exo_coldens( 1,icol,last)
            tint_vals(2) = exo_coldens( 1,icol,next)
        endif
        spc_exo_col(icol) = tint_vals(1) + dels * (tint_vals(2) - tint_vals(1))

      enddo

  end subroutine calc_exo_col

!====================================================================================
  subroutine calc_col_delta(  col_delta_s,          & ! out
                spc_is_inv,   spc_ndx, spc_inv_ndx, & ! in
                spc_exo_col,  vmr,     invariants,  & ! in
                pdel,         ncol                  ) ! in
    !--------------------------------------------------------------------------------
    ! calculate o2,o3 (o) col density above model layer 
    !--------------------------------------------------------------------------------
    use chem_mods, only : gas_pcnst, indexm, nfs

    implicit none
    logical, intent(in) :: spc_is_inv  ! if the index is inversed
    integer, intent(in) :: ncol
    integer, intent(in) :: spc_ndx, spc_inv_ndx  ! index or inverse index of the species
    real(r8),intent(in) :: spc_exo_col(ncol)     ! exo absorber columns [molecules/cm^2]
    real(r8),intent(in) :: vmr(ncol,pver,gas_pcnst)       ! xported species vmr [mol/mol]
    real(r8),intent(in) :: pdel(pcols,pver)               ! pressure delta about midpoints [Pa]
    real(r8),intent(in) :: invariants(ncol,pver,nfs)      ! invariant densities [molecules/cm^3]
    real(r8),intent(out):: col_delta_s(ncol,0:pver)       ! layer column densities [molecules/cm^2]

    !        ... local variables
    integer  :: kk
    !---------------------------------------------------------------
    !        note: xfactor = 10.*r/(k*g) in cgs units.
    !              the factor 10. is to convert pdel
    !              from pascals to dyne/cm**2.
    !---------------------------------------------------------------
    real(r8), parameter :: xfactor = 2.8704e21_r8/(9.80616_r8*1.38044_r8)


     col_delta_s(:,0) = spc_exo_col(:)
     if( spc_is_inv ) then
        do kk = 1,pver
           col_delta_s(:ncol,kk) = xfactor * pdel(:ncol,kk) * invariants(:ncol,kk,spc_inv_ndx)/invariants(:ncol,kk,indexm)
        enddo
     else
        do kk = 1,pver
           col_delta_s(:ncol,kk) = xfactor * pdel(:ncol,kk) * vmr(:ncol,kk,spc_ndx)
        enddo
     endif

  end subroutine calc_col_delta

!====================================================================================
  subroutine setcol( col_delta, & ! in
                     col_dens ) ! out
    !---------------------------------------------------------------
    !     	... set the column densities
    !---------------------------------------------------------------

    use chem_mods, only : nabscol

    implicit none

    !---------------------------------------------------------------
    !     	... dummy arguments
    !---------------------------------------------------------------
    real(r8), intent(in)    :: col_delta(:,0:,:)                 ! layer column densities [molecules/cm^2]
    real(r8), intent(out)   :: col_dens(:,:,:)                   ! column densities [ 1/cm**2 ]

    !---------------------------------------------------------------
    !        the local variables
    !---------------------------------------------------------------
    integer  ::   kk, km1, mm      ! indicies

    !---------------------------------------------------------------
    !   	... compute column densities down to the
    !           current eta index in the calling routine.
    !           the first column is o3 and the second is o2.
    !---------------------------------------------------------------
    do mm = 1,nabscol
       col_dens(:,1,mm) = col_delta(:,0,mm) + .5_r8 * col_delta(:,1,mm)  ! kk=1
       do kk = 2,pver
          km1 = kk - 1
          col_dens(:,kk,mm) = col_dens(:,km1,mm) + .5_r8 * (col_delta(:,km1,mm) + col_delta(:,kk,mm))
       enddo
    enddo

  end subroutine setcol

!====================================================================================
  subroutine photo_timestep_init( calday )
    use time_manager,   only : is_end_curr_day
    use euvac,          only : euvac_set_etf
    use mo_solar_parms, only : solar_parms_get
    use mo_jlong,       only : jlong_timestep_init

    !-----------------------------------------------------------------------------
    !	... setup the time interpolation
    !-----------------------------------------------------------------------------

    implicit none

    !-----------------------------------------------------------------------------
    !	... dummy arguments
    !-----------------------------------------------------------------------------
    real(r8), intent(in)    ::  calday                                   ! day of year at end of present time step

    !-----------------------------------------------------------------------------
    !	... local variables
    !-----------------------------------------------------------------------------
    integer :: m
    real(r8) :: f107
    real(r8) :: f107a

    if( has_o2_col .or. has_o3_col ) then
       if( calday < days(1) ) then
          next = 1
          last = 12
          dels = (365._r8 + calday - days(12)) / (365._r8 + days(1) - days(12))
       else if( calday >= days(12) ) then
          next = 1
          last = 12
          dels = (calday - days(12)) / (365._r8 + days(1) - days(12))
       else
          do m = 11,1,-1
             if( calday >= days(m) ) then
                exit
             end if
          end do
          last = m
          next = m + 1
          dels = (calday - days(m)) / (days(m+1) - days(m))
       end if
#ifdef DEBUG
       write(iulog,*) '-----------------------------------'
       write(iulog,*) 'photo_timestep_init: diagnostics'
       write(iulog,*) 'calday, last, next, dels = ',calday,last,next,dels
       write(iulog,*) '-----------------------------------'
#endif
    end if

    !-----------------------------------------------------------------------
    ! Set jlong etf
    !-----------------------------------------------------------------------
    call jlong_timestep_init
    
  end subroutine photo_timestep_init

end module mo_photo
