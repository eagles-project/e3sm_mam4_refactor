module mo_drydep

  !---------------------------------------------------------------------
  !       ... Dry deposition velocity input data and code for netcdf input
  !---------------------------------------------------------------------

!LKE (10/11/2010): added HCN, CH3CN, HCOOH

  use shr_kind_mod, only : r8 => shr_kind_r8, shr_kind_cl
  use chem_mods,    only : gas_pcnst
  use pmgrid,       only : plev, plevp
  use spmd_utils,   only : masterproc, iam
  use ppgrid,       only : pcols, begchunk, endchunk
  use mo_tracname,  only : solsym
  use cam_abortutils,   only : endrun
  use ioFileMod,    only : getfil
#ifdef SPMD
  use mpishorthand, only : mpicom, mpir8, mpiint, mpilog
#endif
  use pio
  use cam_pio_utils,only : cam_pio_openfile
  use cam_logfile,  only : iulog
  use dyn_grid,     only : get_dyn_grid_parm, get_horiz_grid_d
  use scamMod,      only:  single_column

  use seq_drydep_mod, only : nddvels =>  n_drydep, drydep_list, mapping
  use physconst,    only : karman

  implicit none

  save

  private
  public :: drydep_inti_xactive, drydep_xactive, set_soilw, chk_soilw, has_drydep
  public :: n_land_type, fraction_landuse, drydep_srf_file

  real(r8)              :: dels
  real(r8), allocatable :: days(:)          ! day of year for soilw
  real(r8), allocatable :: dvel(:,:,:,:)    ! depvel array interpolated to model grid
  real(r8), allocatable :: dvel_interp(:,:,:) ! depvel array interpolated to grid and time
  integer :: last, next                     ! day indicies
  integer :: ndays                          ! # of days in soilw file
  integer :: map(gas_pcnst)                 ! indices for drydep species
  integer :: nspecies                       ! number of depvel species in input file

  integer :: no2_ndx, hno3_ndx, o3_ndx, &
             h2o2_ndx, onit_ndx, onitr_ndx, ch4_ndx, ch2o_ndx, &
             ch3ooh_ndx, pooh_ndx, ch3coooh_ndx, c2h5ooh_ndx, eooh_ndx, &
             c3h7ooh_ndx, rooh_ndx, ch3cocho_ndx, co_ndx, ch3coch3_ndx, &
             no_ndx, ho2no2_ndx, glyald_ndx, hyac_ndx, ch3oh_ndx, c2h5oh_ndx, &
             hydrald_ndx, h2_ndx, Pb_ndx, o3s_ndx, o3inert_ndx, macrooh_ndx, &
             xooh_ndx
  integer :: soa_ndx, so4_ndx
  logical :: soa_dd, so4_dd

  logical :: no2_dd, hno3_dd, o3_dd, isopooh_dd, ch4_dd,&
             h2o2_dd, onit_dd, onitr_dd, ch2o_dd, macrooh_dd, xooh_dd, &
             c3h7ooh_dd, rooh_dd, ch3cocho_dd, co_dd, ch3coch3_dd, &
             glyald_dd, hyac_dd, ch3oh_dd, hydrald_dd, h2_dd, Pb_dd, o3s_dd, o3inert_dd

  integer :: so2_ndx
  integer :: hcooh_ndx
  logical :: hcn_dd, hcooh_dd

  integer :: o3a_ndx,xno2_ndx,xhno3_ndx,xonit_ndx,xonitr_ndx,xno_ndx
  logical :: o3a_dd,xno2_dd, xhno3_dd, xonit_dd, xonitr_dd, xno_dd



  integer :: &
       o3_tab_ndx = -1, &
       h2o2_tab_ndx = -1, &
       ch3ooh_tab_ndx = -1, &
       co_tab_ndx = -1
  logical :: &
       o3_in_tab = .false., &
       h2o2_in_tab = .false., &
       ch3ooh_in_tab = .false., &
       co_in_tab = .false.

  real(r8), parameter    :: small_value = 1.e-36_r8
  real(r8), parameter    :: large_value = 1.e36_r8
  real(r8), parameter    :: diffm       = 1.789e-5_r8
  real(r8), parameter    :: diffk       = 1.461e-5_r8
  real(r8), parameter    :: difft       = 2.060e-5_r8
  real(r8), parameter    :: vonkar      = karman
  real(r8), parameter    :: ric         = 0.2_r8
  real(r8), parameter    :: r           = 287.04_r8
  real(r8), parameter    :: cp          = 1004._r8
  real(r8), parameter    :: grav        = 9.81_r8
  real(r8), parameter    :: p00         = 100000._r8
  real(r8), parameter    :: wh2o        = 18.0153_r8
  real(r8), parameter    :: ph          = 1.e-5_r8
  real(r8), parameter    :: ph_inv      = 1._r8/ph
  real(r8), parameter    :: rovcp = r/cp

  integer, pointer :: index_season_lai(:,:)

  logical, public :: has_dvel(gas_pcnst) = .false.
  integer         :: map_dvel(gas_pcnst) = 0
  real(r8) , allocatable            :: soilw_3d(:,:,:)

  logical, parameter :: dyn_soilw = .false.

  real(r8), allocatable  :: fraction_landuse(:,:,:)
  real(r8), allocatable, dimension(:,:,:) :: dep_ra ! [s/m] aerodynamic resistance
  real(r8), allocatable, dimension(:,:,:) :: dep_rb ! [s/m] resistance across sublayer
  integer, parameter :: n_land_type = 11

  integer, allocatable :: spc_ndx(:) ! nddvels
  real(r8), public :: crb 

  type lnd_dvel_type
     real(r8), pointer :: dvel(:,:)   ! deposition velocity over land (cm/s)
  end type lnd_dvel_type

  type(lnd_dvel_type), allocatable :: lnd(:)
  character(len=SHR_KIND_CL) :: drydep_srf_file

contains

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  subroutine dvel_inti_fromlnd
    use mo_chem_utls,         only : get_spc_ndx
    use cam_abortutils,           only : endrun
    use chem_mods,            only : adv_mass
    use seq_drydep_mod,       only : dfoxd

    implicit none

    integer :: ispc, l

    allocate(spc_ndx(nddvels))
    allocate( lnd(begchunk:endchunk) )

    do ispc = 1,nddvels

       spc_ndx(ispc) = get_spc_ndx(drydep_list(ispc))
       if (spc_ndx(ispc) < 1) then
          write(*,*) 'drydep_inti: '//trim(drydep_list(ispc))//' is not included in species set'
          call endrun('drydep_init: invalid dry deposition species')
       endif

    enddo

    crb = (difft/diffm)**(2._r8/3._r8) !.666666_r8

  endsubroutine dvel_inti_fromlnd

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  subroutine drydep_inti_xactive( depvel_lnd_file, clim_soilw_file, season_wes_file )
    !-------------------------------------------------------------------------------------
    ! 	... intialize interactive drydep
    !-------------------------------------------------------------------------------------
    use dycore,        only : dycore_is
    use mo_constants,  only : r2d
    use chem_mods,     only : adv_mass
    use mo_chem_utls,  only : get_spc_ndx
    use seq_drydep_mod,only : drydep_method, DD_XATM, DD_XLND
    use phys_control,  only : phys_getopts

    implicit none

    !-------------------------------------------------------------------------------------
    ! 	... dummy arguments
    !-------------------------------------------------------------------------------------
    character(len=*), intent(in) :: depvel_lnd_file, clim_soilw_file, season_wes_file 

    !-------------------------------------------------------------------------------------
    ! 	... local variables
    !-------------------------------------------------------------------------------------
    integer :: i, j, ii, jj, jl, ju
    integer :: nlon_veg, nlat_veg, npft_veg
    integer :: nlat_lai, npft_lai, pos_min, imin
    integer :: dimid
    integer :: m, n, l, id
    integer :: length1, astat
    integer, allocatable :: wk_lai(:,:,:)
    integer, allocatable :: index_season_lai_j(:,:)
    integer :: k, num_max, k_max
    integer :: num_seas(5)
    integer :: plon, plat
    integer :: ierr

    real(r8)              :: spc_mass
    real(r8)              :: diff_min, target_lat
    real(r8), allocatable :: vegetation_map(:,:,:)
    real(r8), pointer     :: soilw_map(:,:,:)
    real(r8), allocatable :: work(:,:)
    real(r8), allocatable :: landmask(:,:)
    real(r8), allocatable :: urban(:,:)
    real(r8), allocatable :: lake(:,:)
    real(r8), allocatable :: wetland(:,:)
    real(r8), allocatable :: lon_veg(:)
    real(r8), allocatable :: lon_veg_edge(:)
    real(r8), allocatable :: lat_veg(:)
    real(r8), allocatable :: lat_veg_edge(:)
    real(r8), allocatable :: lat_lai(:)
    real(r8), allocatable :: clat(:)
    character(len=32) :: test_name
    type(file_desc_t) :: piofile
    type(var_desc_t) :: vid
    logical :: do_soilw

    character(len=shr_kind_cl) :: locfn
    logical :: prog_modal_aero

    ! determine if modal aerosols are active so that fraction_landuse array is initialized for modal aerosal dry dep
    call phys_getopts(prog_modal_aero_out=prog_modal_aero)

    call dvel_inti_fromlnd()

    if( masterproc ) then
       write(iulog,*) 'drydep_inti: following species have dry deposition'
       do i=1,nddvels
          if( len_trim(drydep_list(i)) > 0 ) then
             write(iulog,*) 'drydep_inti: '//trim(drydep_list(i))//' is requested to have dry dep'
          endif
       enddo
       write(iulog,*) 'drydep_inti:'
    endif

    !-------------------------------------------------------------------------------------
    ! 	... get species indices
    !-------------------------------------------------------------------------------------
    ch4_ndx      = get_spc_ndx( 'CH4' )
    h2_ndx       = get_spc_ndx( 'H2' )
    co_ndx       = get_spc_ndx( 'CO' )
    if( o3_ndx < 0 ) then
       o3_ndx  = get_spc_ndx( 'O3' )
    end if
    so2_ndx     = get_spc_ndx( 'SO2' )
    soa_ndx     = get_spc_ndx( 'SOA' )
    so4_ndx     = get_spc_ndx( 'SO4' )

    soa_dd     = has_drydep( 'SOA' )
    so4_dd     = has_drydep( 'SO4' )


    do i=1,nddvels
       if ( mapping(i) > 0 ) then
          test_name = drydep_list(i)
          m = get_spc_ndx( test_name )
          has_dvel(m) = .true.
          map_dvel(m) = i
       endif
    enddo

    if( all( .not. has_dvel(:) ) ) then
       return
    end if

    !---------------------------------------------------------------------------
    ! 	... allocate module variables
    !---------------------------------------------------------------------------
    allocate( dep_ra(pcols,n_land_type,begchunk:endchunk),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'dvel_inti: failed to allocate dep_ra; error = ',astat
       call endrun
    end if
    allocate( dep_rb(pcols,n_land_type,begchunk:endchunk),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'dvel_inti: failed to allocate dep_rb; error = ',astat
       call endrun
    end if

    if (drydep_method == DD_XLND .and. (.not.prog_modal_aero)) then
       return
    endif

    do_soilw = .not. dyn_soilw .and. (has_drydep( 'H2' ) .or. has_drydep( 'CO' ))
    allocate( fraction_landuse(pcols,n_land_type, begchunk:endchunk),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'dvel_inti: failed to allocate fraction_landuse; error = ',astat
       call endrun
    end if
    if(do_soilw) then
       allocate(soilw_3d(pcols,12,begchunk:endchunk),stat=astat)
       if( astat /= 0 ) then
          write(iulog,*) 'dvel_inti: failed to allocate soilw_3d error = ',astat
          call endrun
       end if
    end if

    plon = get_dyn_grid_parm('plon')
    plat = get_dyn_grid_parm('plat')
    allocate( index_season_lai_j(n_land_type,12),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'dvel_inti: failed to allocate index_season_lai_j; error = ',astat
       call endrun
    end if
    if(dycore_is('UNSTRUCTURED') ) then
       call get_landuse_and_soilw_from_file(do_soilw)
       allocate( index_season_lai(plon,12),stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'dvel_inti: failed to allocate index_season_lai; error = ',astat
          call endrun
       end if
    else
       allocate( index_season_lai(plat,12),stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'dvel_inti: failed to allocate index_season_lai; error = ',astat
          call endrun
       end if
       !---------------------------------------------------------------------------
       ! 	... read landuse map
       !---------------------------------------------------------------------------
       call getfil (depvel_lnd_file, locfn, 0)
       call cam_pio_openfile (piofile, trim(locfn), PIO_NOWRITE)
       !---------------------------------------------------------------------------
       ! 	... get the dimensions
       !---------------------------------------------------------------------------
       ierr = pio_inq_dimid( piofile, 'lon', dimid )
       ierr = pio_inq_dimlen( piofile, dimid, nlon_veg )
       ierr = pio_inq_dimid( piofile, 'lat', dimid )
       ierr = pio_inq_dimlen( piofile, dimid, nlat_veg )
       ierr = pio_inq_dimid( piofile, 'pft', dimid )
       ierr = pio_inq_dimlen( piofile, dimid, npft_veg )
       !---------------------------------------------------------------------------
       ! 	... allocate arrays
       !---------------------------------------------------------------------------
       allocate( vegetation_map(nlon_veg,nlat_veg,npft_veg), work(nlon_veg,nlat_veg), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'dvel_inti: failed to allocate vegation_map; error = ',astat
          call endrun
       end if
       allocate( urban(nlon_veg,nlat_veg), lake(nlon_veg,nlat_veg), &
            landmask(nlon_veg,nlat_veg), wetland(nlon_veg,nlat_veg), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'dvel_inti: failed to allocate vegation_map; error = ',astat
          call endrun
       end if
       allocate( lon_veg(nlon_veg), lat_veg(nlat_veg), &
            lon_veg_edge(nlon_veg+1), lat_veg_edge(nlat_veg+1), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'dvel_inti: failed to allocate vegation lon, lat arrays; error = ',astat
          call endrun
       end if
       !---------------------------------------------------------------------------
       ! 	... read the vegetation map and landmask
       !---------------------------------------------------------------------------
       ierr = pio_inq_varid( piofile, 'PCT_PFT', vid )
       ierr = pio_get_var( piofile, vid, vegetation_map )

       ierr = pio_inq_varid( piofile, 'LANDMASK', vid )
       ierr = pio_get_var( piofile, vid, landmask )

       ierr = pio_inq_varid( piofile, 'PCT_URBAN', vid )
       ierr = pio_get_var( piofile, vid, urban )

       ierr = pio_inq_varid( piofile, 'PCT_LAKE', vid )
       ierr = pio_get_var( piofile, vid, lake )

       ierr = pio_inq_varid( piofile, 'PCT_WETLAND', vid )
       ierr = pio_get_var( piofile, vid, wetland )

       call pio_closefile( piofile )

       !---------------------------------------------------------------------------
       ! scale vegetation, urban, lake, and wetland to fraction
       !---------------------------------------------------------------------------
       vegetation_map(:,:,:) = .01_r8 * vegetation_map(:,:,:)
       wetland(:,:)          = .01_r8 * wetland(:,:)
       lake(:,:)             = .01_r8 * lake(:,:)
       urban(:,:)            = .01_r8 * urban(:,:)
#ifdef DEBUG
       if(masterproc) then
          write(iulog,*) 'minmax vegetation_map ',minval(vegetation_map),maxval(vegetation_map)
          write(iulog,*) 'minmax wetland        ',minval(wetland),maxval(wetland)
          write(iulog,*) 'minmax landmask       ',minval(landmask),maxval(landmask)
       end if
#endif
       !---------------------------------------------------------------------------
       ! 	... define lat-lon of vegetation map (1x1)
       !---------------------------------------------------------------------------
       lat_veg(:)      = (/ (-89.5_r8 + (i-1),i=1,nlat_veg  ) /)
       lon_veg(:)      = (/ (  0.5_r8 + (i-1),i=1,nlon_veg  ) /)
       lat_veg_edge(:) = (/ (-90.0_r8 + (i-1),i=1,nlat_veg+1) /)
       lon_veg_edge(:) = (/ (  0.0_r8 + (i-1),i=1,nlon_veg+1) /)
       !---------------------------------------------------------------------------
       ! 	... read soilw table if necessary
       !---------------------------------------------------------------------------

       if( do_soilw ) then
          call soilw_inti( clim_soilw_file, nlon_veg, nlat_veg, soilw_map )
       end if

       !---------------------------------------------------------------------------
       ! 	... regrid to model grid
       !---------------------------------------------------------------------------

       call interp_map( plon, plat, nlon_veg, nlat_veg, npft_veg, lat_veg, lat_veg_edge, &
            lon_veg, lon_veg_edge, landmask, urban, lake, &
            wetland, vegetation_map, soilw_map, do_soilw )

       deallocate( vegetation_map, work, stat=astat )
       deallocate( lon_veg, lat_veg, lon_veg_edge, lat_veg_edge, stat=astat )
       deallocate( landmask, urban, lake, wetland, stat=astat )
       if( do_soilw ) then
          deallocate( soilw_map, stat=astat )
       end if
    endif  ! Unstructured grid

    if (drydep_method == DD_XLND) then
       return
    endif

    !---------------------------------------------------------------------------
    ! 	... read LAI based season indeces
    !---------------------------------------------------------------------------
    call getfil (season_wes_file, locfn, 0)
    call cam_pio_openfile (piofile, trim(locfn), PIO_NOWRITE)
    !---------------------------------------------------------------------------
    ! 	... get the dimensions
    !---------------------------------------------------------------------------
    ierr = pio_inq_dimid( piofile, 'lat', dimid )
    ierr = pio_inq_dimlen( piofile, dimid, nlat_lai )
    ierr = pio_inq_dimid( piofile, 'pft', dimid )
    ierr = pio_inq_dimlen( piofile, dimid, npft_lai )
    !---------------------------------------------------------------------------
    ! 	... allocate arrays
    !---------------------------------------------------------------------------
    allocate( lat_lai(nlat_lai), wk_lai(nlat_lai,npft_lai,12), stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'dvel_inti: failed to allocate vegation_map; error = ',astat
       call endrun
    end if
    !---------------------------------------------------------------------------
    ! 	... read the latitude and the season indicies
    !---------------------------------------------------------------------------
    ierr = pio_inq_varid( piofile, 'lat', vid )
    ierr = pio_get_var( piofile, vid, lat_lai )

    ierr = pio_inq_varid( piofile, 'season_wes', vid )
    ierr = pio_get_var( piofile, vid, wk_lai )

    call pio_closefile( piofile )


    if(dycore_is('UNSTRUCTURED') ) then
       ! For unstructured grids plon is the 1d horizontal grid size and plat=1
       ! So this code averages at the latitude of each grid point - not an ideal solution
       allocate(clat(plon))
       call get_horiz_grid_d(plon, clat_d_out=clat)
       jl = 1
       ju = plon
    else
       allocate(clat(plat))
       call get_horiz_grid_d(plat, clat_d_out=clat)
       jl = 1
       ju = plat
    end if
    imin = 1
    do j = 1,ju
       diff_min = 10._r8
       pos_min  = -99
       target_lat = clat(j)*r2d
       do i = imin,nlat_lai
          if( abs(lat_lai(i) - target_lat) < diff_min ) then
             diff_min = abs(lat_lai(i) - target_lat)
             pos_min  = i
          end if
       end do
       if( pos_min < 0 ) then
          write(iulog,*) 'dvel_inti: cannot find ',target_lat,' at j,pos_min,diff_min = ',j,pos_min,diff_min
          write(iulog,*) 'dvel_inti: imin,nlat_lai = ',imin,nlat_lai
          write(iulog,*) 'dvel_inti: lat_lai'
          write(iulog,'(1p,10g12.5)') lat_lai(:)
          call endrun
       end if
       if(dycore_is('UNSTRUCTURED') ) then
          imin=1
       else
          imin = pos_min
       end if
       index_season_lai_j(:,:) = wk_lai(pos_min,:,:)

       !---------------------------------------------------------------------------
       ! specify the season as the most frequent in the 11 vegetation classes
       ! this was done to remove a banding problem in dvel (JFL Oct 04)
       !---------------------------------------------------------------------------
       do m = 1,12
          num_seas = 0
          do l = 1,11
             do k = 1,5
                if( index_season_lai_j(l,m) == k ) then
                   num_seas(k) = num_seas(k) + 1
                   exit
                end if
             end do
          end do

          num_max = -1
          do k = 1,5
             if( num_seas(k) > num_max ) then
                num_max = num_seas(k)
                k_max = k
             endif
          end do

          index_season_lai(j,m) = k_max
       end do
    end do

    deallocate( lat_lai, wk_lai, clat, index_season_lai_j)

  end subroutine drydep_inti_xactive

  !-------------------------------------------------------------------------------------
  subroutine get_landuse_and_soilw_from_file(do_soilw)
    use cam_pio_utils, only : cam_pio_openfile
    use ncdio_atm, only : infld
    use cam_control_mod, only: aqua_planet
    logical, intent(in) :: do_soilw
    logical :: readvar
    
    type(file_desc_t) :: piofile
    character(len=shr_kind_cl) :: locfn
    logical :: lexist
    
    if (aqua_planet) then
      fraction_landuse = 0.
    else

      call getfil (drydep_srf_file, locfn, 1, lexist)
      if(lexist) then
         call cam_pio_openfile(piofile, locfn, PIO_NOWRITE)

         call infld('fraction_landuse', piofile, 'ncol','class',' ',1,pcols,1,n_land_type, begchunk,endchunk, &
              fraction_landuse, readvar, gridname='physgrid')

         if(do_soilw) then
            call infld('soilw', piofile, 'ncol','month',' ',1,pcols,1,12, begchunk,endchunk, &
                 soilw_3d, readvar, gridname='physgrid')
         end if

         call pio_closefile(piofile)
      else
         call endrun('Unstructured grids require drydep_srf_file ')
      end if
    
    end if ! aqua_planet
  end subroutine get_landuse_and_soilw_from_file

  !-------------------------------------------------------------------------------------
  subroutine interp_map( plon, plat, nlon_veg, nlat_veg, npft_veg, lat_veg, lat_veg_edge, &
                         lon_veg, lon_veg_edge, landmask, urban, lake, &
                         wetland, vegetation_map, soilw_map, do_soilw )

    use mo_constants, only : r2d
    use scamMod, only : latiop,loniop,scmlat,scmlon,use_replay
    use shr_scam_mod  , only: shr_scam_getCloseLatLon  ! Standardized system subroutines
    use filenames, only: ncdata
    use dycore, only : dycore_is
    use phys_grid, only : scatter_field_to_chunk
    implicit none

    !-------------------------------------------------------------------------------------
    ! 	... dummy arguments
    !-------------------------------------------------------------------------------------
    integer, intent(in)      ::  plon, plat, nlon_veg, nlat_veg, npft_veg
    real(r8), pointer            :: soilw_map(:,:,:)
    real(r8), intent(in)         :: landmask(nlon_veg,nlat_veg)
    real(r8), intent(in)         :: urban(nlon_veg,nlat_veg)
    real(r8), intent(in)         :: lake(nlon_veg,nlat_veg)
    real(r8), intent(in)         :: wetland(nlon_veg,nlat_veg)
    real(r8), intent(in)         :: vegetation_map(nlon_veg,nlat_veg,npft_veg)
    real(r8), intent(in)         :: lon_veg(nlon_veg)
    real(r8), intent(in)         :: lon_veg_edge(nlon_veg+1)
    real(r8), intent(in)         :: lat_veg(nlat_veg)
    real(r8), intent(in)         :: lat_veg_edge(nlat_veg+1)
    logical,  intent(in)         :: do_soilw

    !-------------------------------------------------------------------------------------
    ! 	... local variables
    !-------------------------------------------------------------------------------------
    real(r8) :: closelat,closelon
    integer :: latidx,lonidx

    integer, parameter           :: veg_ext = 20
    type(file_desc_t)            :: piofile
    integer                      :: i, j, ii, jj, jl, ju, i_ndx, n
    integer, dimension(plon+1)   :: ind_lon
    integer, dimension(plat+1)  :: ind_lat
    real(r8)                         :: total_land
    real(r8), dimension(plon+1)      :: lon_edge
    real(r8), dimension(plat+1)     :: lat_edge
    real(r8)                         :: lat1, lat2, lon1, lon2
    real(r8)                         :: x1, x2, y1, y2, dx, dy
    real(r8)                         :: area, total_area
    real(r8), dimension(npft_veg+3)  :: fraction
    real(r8)                         :: total_soilw_area
    real(r8)                         :: fraction_soilw
    real(r8)                         :: total_soilw(12)
    
    real(r8),    dimension(-veg_ext:nlon_veg+veg_ext) :: lon_veg_edge_ext
    integer, dimension(-veg_ext:nlon_veg+veg_ext) :: mapping_ext

    real(r8), allocatable :: lam(:), phi(:), garea(:)

    logical, parameter :: has_npole = .true.
    integer :: ploniop,platiop
    character(len=shr_kind_cl) :: ncdata_loc
    real(r8) :: tmp_frac_lu(plon,n_land_type,plat), tmp_soilw_3d(plon,12,plat)

    allocate(lam(plon), phi(plat))
    call get_horiz_grid_d(plon, clon_d_out=lam)
    call get_horiz_grid_d(plat, clat_d_out=phi)



    jl = 1
    ju = plon

    if (single_column) then
       if (use_replay) then
          call getfil (ncdata, ncdata_loc)
          call cam_pio_openfile (piofile, trim(ncdata_loc), PIO_NOWRITE)
          call shr_scam_getCloseLatLon(piofile,scmlat,scmlon,closelat,closelon,latidx,lonidx)
          call pio_closefile ( piofile)
          ploniop=size(loniop)
          platiop=size(latiop)
       else 
          latidx=1
          lonidx=1
          ploniop=1
          platiop=1
       end if
      
       lon_edge(1) = loniop(lonidx) * r2d - .5_r8*(loniop(2) - loniop(1)) * r2d

       if (lonidx.lt.ploniop) then
          lon_edge(2) = loniop(lonidx+1) * r2d - .5_r8*(loniop(2) - loniop(1)) * r2d
       else
          lon_edge(2) = lon_edge(1) + (loniop(2) - loniop(1)) * r2d
       end if

       lat_edge(1) = latiop(latidx) * r2d - .5_r8*(latiop(2) - latiop(1)) * r2d

       if (latidx.lt.platiop) then
          lat_edge(2) = latiop(latidx+1) * r2d - .5_r8*(latiop(2) - latiop(1)) * r2d
       else
          lat_edge(2) = lat_edge(1) + (latiop(2) - latiop(1)) * r2d
       end if       
    else
       do i = 1,plon
          lon_edge(i) = lam(i) * r2d - .5_r8*(lam(2) - lam(1)) * r2d
       end do
       lon_edge(plon+1) = lon_edge(plon) + (lam(2) - lam(1)) * r2d
       if( .not. has_npole ) then
          do j = 1,plat+1
             lat_edge(j) = phi(j) * r2d - .5_r8*(phi(2) - phi(1)) * r2d
          end do
       else
          do j = 1,plat
             lat_edge(j) = phi(j) * r2d - .5_r8*(phi(2) - phi(1)) * r2d
          end do
          lat_edge(plat+1) = lat_edge(plat) + (phi(2) - phi(1)) * r2d
       end if
    end if
    do j = 1,plat+1
       lat_edge(j) = min( lat_edge(j), 90._r8 )
       lat_edge(j) = max( lat_edge(j),-90._r8 )
    end do

    !-------------------------------------------------------------------------------------
    ! wrap around the longitudes
    !-------------------------------------------------------------------------------------
    do i = -veg_ext,0
       lon_veg_edge_ext(i) = lon_veg_edge(nlon_veg+i) - 360._r8
       mapping_ext     (i) =              nlon_veg+i
    end do
    do i = 1,nlon_veg
       lon_veg_edge_ext(i) = lon_veg_edge(i)
       mapping_ext     (i) =              i
    end do
    do i = nlon_veg+1,nlon_veg+veg_ext
       lon_veg_edge_ext(i) = lon_veg_edge(i-nlon_veg) + 360._r8
       mapping_ext     (i) =              i-nlon_veg
    end do
#ifdef DEBUG
    write(iulog,*) 'interp_map : lon_edge ',lon_edge
    write(iulog,*) 'interp_map : lat_edge ',lat_edge
    write(iulog,*) 'interp_map : mapping_ext ',mapping_ext
#endif
    do j = 1,plon+1
       lon1 = lon_edge(j) 
       do i = -veg_ext,nlon_veg+veg_ext
          dx = lon_veg_edge_ext(i  ) - lon1
          dy = lon_veg_edge_ext(i+1) - lon1
          if( dx*dy <= 0._r8 ) then
             ind_lon(j) = i
             exit
          end if
       end do
    end do

    do j = 1,plat+1
       lat1 = lat_edge(j)
       do i = 1,nlat_veg
          dx = lat_veg_edge(i  ) - lat1
          dy = lat_veg_edge(i+1) - lat1
          if( dx*dy <= 0._r8 ) then
             ind_lat(j) = i
             exit
          end if
       end do
    end do
#ifdef DEBUG
    write(iulog,*) 'interp_map : ind_lon ',ind_lon
    write(iulog,*) 'interp_map : ind_lat ',ind_lat
#endif
    lat_loop : do j = 1,plat
       lon_loop : do i = 1,plon
          total_area       = 0._r8
          fraction         = 0._r8
          total_soilw(:)   = 0._r8
          total_soilw_area = 0._r8
          do jj = ind_lat(j),ind_lat(j+1)
             y1 = max( lat_edge(j),lat_veg_edge(jj) )
             y2 = min( lat_edge(j+1),lat_veg_edge(jj+1) ) 
             dy = (y2 - y1)/(lat_veg_edge(jj+1) - lat_veg_edge(jj))
             do ii =ind_lon(i),ind_lon(i+1)
                i_ndx = mapping_ext(ii)
                x1 = max( lon_edge(i),lon_veg_edge_ext(ii) )
                x2 = min( lon_edge(i+1),lon_veg_edge_ext(ii+1) ) 
                dx = (x2 - x1)/(lon_veg_edge_ext(ii+1) - lon_veg_edge_ext(ii))
                area = dx * dy
                total_area = total_area + area
                !-----------------------------------------------------------------
                ! 	... special case for ocean grid point 
                !-----------------------------------------------------------------
                if( nint(landmask(i_ndx,jj)) == 0 ) then
                   fraction(npft_veg+1) = fraction(npft_veg+1) + area
                else
                   do n = 1,npft_veg
                      fraction(n) = fraction(n) + vegetation_map(i_ndx,jj,n) * area
                   end do
                   fraction(npft_veg+1) = fraction(npft_veg+1) + area * lake   (i_ndx,jj)
                   fraction(npft_veg+2) = fraction(npft_veg+2) + area * wetland(i_ndx,jj)
                   fraction(npft_veg+3) = fraction(npft_veg+3) + area * urban  (i_ndx,jj)
                   !-----------------------------------------------------------------
                   ! 	... check if land accounts for the whole area.
                   !           If not, the remaining area is in the ocean
                   !-----------------------------------------------------------------
                   total_land = sum(vegetation_map(i_ndx,jj,:)) &
                              + urban  (i_ndx,jj) &
                              + lake   (i_ndx,jj) &
                              + wetland(i_ndx,jj)
                   if( total_land < 1._r8 ) then
                      fraction(npft_veg+1) = fraction(npft_veg+1) + (1._r8 - total_land) * area
                   end if
                   !-------------------------------------------------------------------------------------
                   ! 	... compute weighted average of soilw over grid (non-water only)
                   !-------------------------------------------------------------------------------------
                   if( do_soilw ) then
                      fraction_soilw = total_land  - (lake(i_ndx,jj) + wetland(i_ndx,jj))
                      total_soilw_area = total_soilw_area + fraction_soilw * area
                      total_soilw(:)   = total_soilw(:) + fraction_soilw * area * soilw_map(i_ndx,jj,:)
                   end if
                end if
             end do
          end do
          !-------------------------------------------------------------------------------------
          ! 	... divide by total area of grid box
          !-------------------------------------------------------------------------------------
          fraction(:) = fraction(:)/total_area
          !-------------------------------------------------------------------------------------
          ! 	... make sure we don't have too much or too little
          !-------------------------------------------------------------------------------------
          if( abs( sum(fraction) - 1._r8) > .001_r8 ) then
             fraction(:) = fraction(:)/sum(fraction)
          end if
          !-------------------------------------------------------------------------------------
          ! 	... map to Wesely land classification
          !-------------------------------------------------------------------------------------

          


          tmp_frac_lu(i, 1, j) =     fraction(20)
          tmp_frac_lu(i, 2, j) = sum(fraction(16:17))
          tmp_frac_lu(i, 3, j) = sum(fraction(13:15))
          tmp_frac_lu(i, 4, j) = sum(fraction( 5: 9))
          tmp_frac_lu(i, 5, j) = sum(fraction( 2: 4))
          tmp_frac_lu(i, 6, j) =     fraction(19)
          tmp_frac_lu(i, 7, j) =     fraction(18)
          tmp_frac_lu(i, 8, j) =     fraction( 1)
          tmp_frac_lu(i, 9, j) = 0._r8
          tmp_frac_lu(i,10, j) = 0._r8
          tmp_frac_lu(i,11, j) = sum(fraction(10:12))
          if( do_soilw ) then
             if( total_soilw_area > 0._r8 ) then
                tmp_soilw_3d(i,:,j) = total_soilw(:)/total_soilw_area
             else
                tmp_soilw_3d(i,:,j) = -99._r8
             end if
          end if
       end do lon_loop
    end do lat_loop
    !-------------------------------------------------------------------------------------
    ! 	... reshape according to lat-lon blocks
    !-------------------------------------------------------------------------------------
    call scatter_field_to_chunk(1,n_land_type,1,plon,tmp_frac_lu,fraction_landuse)
    if(do_soilw) call scatter_field_to_chunk(1,12,1,plon,tmp_soilw_3d,soilw_3d)
    !-------------------------------------------------------------------------------------
    ! 	... make sure there are no out of range values
    !-------------------------------------------------------------------------------------
    where (fraction_landuse < 0._r8) fraction_landuse = 0._r8
    where (fraction_landuse > 1._r8) fraction_landuse = 1._r8

  end subroutine interp_map
  
  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  subroutine drydep_xactive( ncdate, sfc_temp, pressure_sfc,  &
                             wind_speed, spec_hum, air_temp, pressure_10m, rain, &
                             snow, solar_flux, dvel, dflx, mmr, &
                             tv, soilw, rh, ncol, lonndx, latndx, lchnk, &
                             ocnfrc, icefrc, beglandtype, endlandtype )
    !-------------------------------------------------------------------------------------
    !   code based on wesely (atmospheric environment, 1989, vol 23, p. 1293-1304) for
    !   calculation of r_c, and on walcek et. al. (atmospheric enviroment, 1986,
    !   vol. 20, p. 949-964) for calculation of r_a and r_b
    !
    !   as suggested in walcek (u_i)(u*_i) = (u_a)(u*_a)
    !   is kept constant where i represents a subgrid environment and a the
    !   grid average environment. thus the calculation proceeds as follows:
    !   va the grid averaged wind is calculated on dots
    !   z0(i) the grid averaged roughness coefficient is calculated
    !   ri(i) the grid averaged richardson number is calculated
    !   --> the grid averaged (u_a)(u*_a) is calculated
    !   --> subgrid scale u*_i is calculated assuming (u_i) given as above
    !   --> final deposotion velocity is weighted average of subgrid scale velocities
    !
    ! code written by P. Hess, rewritten in fortran 90 by JFL (August 2000)
    ! modified by JFL to be used in MOZART-2 (October 2002)
    !-------------------------------------------------------------------------------------

    use seq_drydep_mod, only: z0, rgso, rgss, h2_a, h2_b, h2_c, ri, rclo, rcls, rlu, rac
    use seq_drydep_mod, only: seq_drydep_setHCoeff, foxd, drat
    use physconst,      only: tmelt
    use seq_drydep_mod, only: drydep_method,  DD_XLND

    implicit none

    !-------------------------------------------------------------------------------------
    ! 	... dummy arguments
    !-------------------------------------------------------------------------------------
    integer, intent(in)   :: ncol
    integer, intent(in)   :: ncdate                   ! present date (yyyymmdd)
    real(r8), intent(in)      :: sfc_temp(pcols)          ! surface temperature (K)
    real(r8), intent(in)      :: pressure_sfc(pcols)      ! surface pressure (Pa)
    real(r8), intent(in)      :: wind_speed(pcols)        ! 10 meter wind speed (m/s)
    real(r8), intent(in)      :: spec_hum(pcols)          ! specific humidity (kg/kg)
    real(r8), intent(in)      :: rh(ncol,1)               ! relative humidity
    real(r8), intent(in)      :: air_temp(pcols)          ! surface air temperature (K)
    real(r8), intent(in)      :: pressure_10m(pcols)      ! 10 meter pressure (Pa)
    real(r8), intent(in)      :: rain(pcols)              
    real(r8), intent(in)      :: snow(pcols)              ! snow height (m)
    real(r8), intent(in)      :: soilw(pcols)             ! soil moisture fraction
    real(r8), intent(in)      :: solar_flux(pcols)        ! direct shortwave radiation at surface (W/m^2)
    real(r8), intent(in)      :: tv(pcols)                ! potential temperature
    real(r8), intent(in)      :: mmr(pcols,plev,gas_pcnst)    ! constituent concentration (kg/kg)
    real(r8), intent(out)     :: dvel(ncol,gas_pcnst)        ! deposition velocity (cm/s)
    real(r8), intent(inout)   :: dflx(pcols,gas_pcnst)        ! deposition flux (/cm^2/s)

    integer, intent(in)     ::   latndx(pcols)           ! chunk latitude indicies
    integer, intent(in)     ::   lonndx(pcols)           ! chunk longitude indicies
    integer, intent(in)     ::   lchnk                   ! chunk number

    integer, intent(in), optional     ::  beglandtype
    integer, intent(in), optional     ::  endlandtype

    real(r8), intent(in), optional      :: ocnfrc(pcols) 
    real(r8), intent(in), optional      :: icefrc(pcols) 

    !-------------------------------------------------------------------------------------
    ! 	... local variables
    !-------------------------------------------------------------------------------------
    real(r8), parameter :: scaling_to_cm_per_s = 100._r8
    real(r8), parameter :: rain_threshold      = 1.e-7_r8  ! of the order of 1cm/day expressed in m/s

    integer :: i, ispec, lt, m
    integer :: sndx
    integer :: month

    real(r8) :: slope = 0._r8
    real(r8) :: z0water ! revised z0 over water
    real(r8) :: p       ! pressure at midpoint first layer
    real(r8) :: pg      ! surface pressure
    real(r8) :: es      ! saturation vapor pressure
    real(r8) :: ws      ! saturation mixing ratio
    real(r8) :: hvar    ! constant to compute xmol
    real(r8) :: h       ! constant to compute xmol
    real(r8) :: psih    ! stability correction factor
    real(r8) :: rs      ! constant for calculating rsmx
    real(r8) :: rmx     ! resistance by vegetation
    real(r8) :: zovl    ! ratio of z to  m-o length
    real(r8) :: cvarb   ! cvar averaged over landtypes
    real(r8) :: bb      ! b averaged over landtypes
    real(r8) :: ustarb  ! ustar averaged over landtypes
    real(r8) :: tc(ncol)  ! temperature in celsius
    real(r8) :: cts(ncol) ! correction to rlu rcl and rgs for frost

    !-------------------------------------------------------------------------------------
    ! local arrays: dependent on location and species
    !-------------------------------------------------------------------------------------
    real(r8), dimension(ncol,nddvels) :: heff

    !-------------------------------------------------------------------------------------
    ! local arrays: dependent on location only
    !-------------------------------------------------------------------------------------
    integer                :: index_season(ncol,n_land_type)
    real(r8), dimension(ncol) :: tha     ! atmospheric virtual potential temperature
    real(r8), dimension(ncol) :: thg     ! ground virtual potential temperature
    real(r8), dimension(ncol) :: z       ! height of lowest level
    real(r8), dimension(ncol) :: va      ! magnitude of v on cross points
    real(r8), dimension(ncol) :: ribn    ! richardson number
    real(r8), dimension(ncol) :: qs      ! saturation specific humidity
    real(r8), dimension(ncol) :: crs     ! multiplier to calculate crs
    real(r8), dimension(ncol) :: rdc     ! part of lower canopy resistance
    real(r8), dimension(ncol) :: uustar  ! u*ustar (assumed constant over grid)
    real(r8), dimension(ncol) :: z0b     ! average roughness length over grid
    real(r8), dimension(ncol) :: wrk     ! work array
    real(r8), dimension(ncol) :: term    ! work array
    real(r8), dimension(ncol) :: resc    ! work array
    real(r8), dimension(ncol) :: lnd_frc ! work array
    logical,  dimension(ncol) :: unstable
    logical,  dimension(ncol) :: has_rain
    logical,  dimension(ncol) :: has_dew

    !-------------------------------------------------------------------------------------
    ! local arrays: dependent on location and landtype
    !-------------------------------------------------------------------------------------
    real(r8), dimension(ncol,n_land_type) :: b     ! buoyancy parameter for unstable conditions
    real(r8), dimension(ncol,n_land_type) :: cvar  ! height parameter
    real(r8), dimension(ncol,n_land_type) :: ustar ! friction velocity
    real(r8), dimension(ncol,n_land_type) :: xmol  ! monin-obukhov length

    !-------------------------------------------------------------------------------------
    ! local arrays: dependent on location, landtype and species
    !-------------------------------------------------------------------------------------
    real(r8), dimension(ncol,n_land_type,gas_pcnst) :: rsmx  ! vegetative resistance (plant mesophyll)
    real(r8), dimension(ncol,n_land_type,gas_pcnst) :: rclx  ! lower canopy resistance
    real(r8), dimension(ncol,n_land_type,gas_pcnst) :: rlux  ! vegetative resistance (upper canopy)
    real(r8), dimension(ncol,n_land_type) :: rlux_o3  ! vegetative resistance (upper canopy)
    real(r8), dimension(ncol,n_land_type,gas_pcnst) :: rgsx  ! ground resistance
    real(r8) :: pmid(ncol,1)                             ! for seasalt aerosols
    real(r8) :: tfld(ncol,1)                             ! for seasalt aerosols
    real(r8) :: fact
    real(r8) :: rc                                    ! combined surface resistance
    real(r8) :: var_soilw, dv_soil_h2, fact_h2        ! h2 dvel wrking variables
    logical  :: fr_lnduse(ncol,n_land_type)           ! wrking array
    real(r8) :: dewm                                  ! multiplier for rs when dew occurs

    real(r8) :: lcl_frc_landuse(ncol,n_land_type) 

    integer :: beglt, endlt


    if (present( beglandtype)) then
      beglt = beglandtype 
    else
      beglt = 1
    endif
    if (present( endlandtype)) then
      endlt = endlandtype 
    else
      endlt = n_land_type
    endif
  
    !-------------------------------------------------------------------------------------
    ! initialize
    !-------------------------------------------------------------------------------------
    do m = 1,gas_pcnst
       dvel(:,m) = 0._r8
    end do

    if( all( .not. has_dvel(:) ) ) then
       return
    end if

    !-------------------------------------------------------------------------------------
    ! define species-dependent parameters (temperature dependent)
    !-------------------------------------------------------------------------------------
    call seq_drydep_setHCoeff( ncol, sfc_temp, heff )

    do lt = 1,n_land_type
       dep_ra (:,lt,lchnk)   = 0._r8
       dep_rb (:,lt,lchnk)   = 0._r8
    end do

    !-------------------------------------------------------------------------------------
    ! 	... set month
    !-------------------------------------------------------------------------------------
    month = mod( ncdate,10000 )/100

    !-------------------------------------------------------------------------------------
    ! define which season (relative to Northern hemisphere climate)
    !-------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------
    ! define season index based on fixed LAI
    !-------------------------------------------------------------------------------------
    do i = 1,ncol
       index_season(i,:) = index_season_lai(latndx(i),month)
    enddo
    !-------------------------------------------------------------------------------------
    ! special case for snow covered terrain
    !-------------------------------------------------------------------------------------
    do i = 1,ncol
       if( snow(i) > .01_r8 ) then
          index_season(i,:) = 4
       endif
    enddo
    !-------------------------------------------------------------------------------------
    ! scale rain and define logical arrays
    !-------------------------------------------------------------------------------------
    has_rain(:ncol) = rain(:ncol) > rain_threshold

    !-------------------------------------------------------------------------------------
    ! loop over longitude points
    !-------------------------------------------------------------------------------------
    col_loop :  do i = 1,ncol
       p   = pressure_10m(i)
       pg  = pressure_sfc(i)
       !-------------------------------------------------------------------------------------
       ! potential temperature
       !-------------------------------------------------------------------------------------
       tha(i) = air_temp(i) * (p00/p )**rovcp * (1._r8 + .61_r8*spec_hum(i))
       thg(i) = sfc_temp(i) * (p00/pg)**rovcp * (1._r8 + .61_r8*spec_hum(i))
       !-------------------------------------------------------------------------------------
       ! height of 1st level
       !-------------------------------------------------------------------------------------
       z(i) = - r/grav * air_temp(i) * (1._r8 + .61_r8*spec_hum(i)) * log(p/pg)
       !-------------------------------------------------------------------------------------
       ! wind speed
       !-------------------------------------------------------------------------------------
       va(i) = max( .01_r8,wind_speed(i) )
       !-------------------------------------------------------------------------------------
       ! Richardson number
       !-------------------------------------------------------------------------------------
       ribn(i) = z(i) * grav * (tha(i) - thg(i))/thg(i) / (va(i)*va(i))
       ribn(i) = min( ribn(i),ric )
       unstable(i) = ribn(i) < 0._r8
       !-------------------------------------------------------------------------------------
       ! saturation vapor pressure (Pascals)
       ! saturation mixing ratio
       ! saturation specific humidity
       !-------------------------------------------------------------------------------------
       es    = 611._r8*exp( 5414.77_r8*(sfc_temp(i) - tmelt)/(tmelt*sfc_temp(i)) )
       ws    = .622_r8*es/(pg - es)
       qs(i) = ws/(1._r8 + ws)
       has_dew(i) = .false.
       if( qs(i) <= spec_hum(i) ) then
          has_dew(i) = .true.
       endif
       if( sfc_temp(i) < tmelt ) then
          has_dew(i) = .false.
       endif
       !-------------------------------------------------------------------------------------
       ! constant in determining rs
       !-------------------------------------------------------------------------------------
       tc(i) = sfc_temp(i) - tmelt
       if( sfc_temp(i) > tmelt .and. sfc_temp(i) < 313.15_r8 ) then
          crs(i) = (1._r8 + (200._r8/(solar_flux(i) + .1_r8))**2) * (400._r8/(tc(i)*(40._r8 - tc(i))))
       else
          crs(i) = large_value
       endif
       !-------------------------------------------------------------------------------------
       ! rdc (lower canopy res)
       !-------------------------------------------------------------------------------------
       rdc(i) = 100._r8*(1._r8 + 1000._r8/(solar_flux(i) + 10._r8))/(1._r8 + 1000._r8*slope)
    enddo col_loop

    !-------------------------------------------------------------------------------------
    ! 	... form working arrays
    !-------------------------------------------------------------------------------------
    do lt = 1,n_land_type
       do i=1,ncol
          lcl_frc_landuse(i,lt) = fraction_landuse(i,lt,lchnk)
       enddo
    enddo
    if ( present(ocnfrc) .and. present(icefrc) ) then
       do i=1,ncol
          ! land type 7 is used for ocean
          ! land type 8 is used for sea ice
          lcl_frc_landuse(i,7) = ocnfrc(i)
          lcl_frc_landuse(i,8) = icefrc(i)
       enddo
    endif
    do lt = 1,n_land_type
       do i=1,ncol
          fr_lnduse(i,lt) = lcl_frc_landuse(i,lt) > 0._r8
       enddo
    enddo

    !-------------------------------------------------------------------------------------
    ! find grid averaged z0: z0bar (the roughness length) z_o=exp[S(f_i*ln(z_oi))]
    ! this is calculated so as to find u_i, assuming u*u=u_i*u_i
    !-------------------------------------------------------------------------------------
    z0b(:) = 0._r8
    do lt = 1,n_land_type
       do i = 1,ncol
          if( fr_lnduse(i,lt) ) then
             z0b(i) = z0b(i) + lcl_frc_landuse(i,lt) * log( z0(index_season(i,lt),lt) )
          endif
       enddo
    enddo

    !-------------------------------------------------------------------------------------
    ! find the constant velocity uu*=(u_i)(u*_i)
    !-------------------------------------------------------------------------------------
    do i = 1,ncol
       z0b(i) = exp( z0b(i) )
       cvarb  = vonkar/log( z(i)/z0b(i) )
       !-------------------------------------------------------------------------------------
       ! unstable and stable cases
       !-------------------------------------------------------------------------------------
       if( unstable(i) ) then
          bb = 9.4_r8*(cvarb**2)*sqrt( abs(ribn(i))*z(i)/z0b(i) )
          ustarb = cvarb * va(i) * sqrt( 1._r8 - (9.4_r8*ribn(i)/(1._r8 + 7.4_r8*bb)) )
       else
          ustarb = cvarb * va(i)/(1._r8 + 4.7_r8*ribn(i))
       endif
       uustar(i) = va(i)*ustarb
    enddo

    !-------------------------------------------------------------------------------------
    ! calculate the friction velocity for each land type u_i=uustar/u*_i
    !-------------------------------------------------------------------------------------
    do lt = beglt,endlt
       do i = 1,ncol
          if( fr_lnduse(i,lt) ) then
             if( unstable(i) ) then
                cvar(i,lt)  = vonkar/log( z(i)/z0(index_season(i,lt),lt) )
                b(i,lt)     = 9.4_r8*(cvar(i,lt)**2)* sqrt( abs(ribn(i))*z(i)/z0(index_season(i,lt),lt) )
                ustar(i,lt) = sqrt( cvar(i,lt)*uustar(i)*sqrt( 1._r8 - (9.4_r8*ribn(i)/(1._r8 + 7.4_r8*b(i,lt))) ) )
             else
                cvar(i,lt)  = vonkar/log( z(i)/z0(index_season(i,lt),lt) )
                ustar(i,lt) = sqrt( cvar(i,lt)*uustar(i)/(1._r8 + 4.7_r8*ribn(i)) )
             endif
          endif
       enddo
    enddo

    !-------------------------------------------------------------------------------------
    ! revise calculation of friction velocity and z0 over water
    !-------------------------------------------------------------------------------------
    lt = 7    
    do i = 1,ncol
       if( fr_lnduse(i,lt) ) then
          if( unstable(i) ) then
             z0water     = (.016_r8*(ustar(i,lt)**2)/grav) + diffk/(9.1_r8*ustar(i,lt))
             cvar(i,lt)  = vonkar/(log( z(i)/z0water ))
             b(i,lt)     = 9.4_r8*(cvar(i,lt)**2)*sqrt( abs(ribn(i))*z(i)/z0water )
             ustar(i,lt) = sqrt( cvar(i,lt)*uustar(i)* sqrt( 1._r8 - (9.4_r8*ribn(i)/(1._r8+ 7.4_r8*b(i,lt))) ) )
          else
             z0water     = (.016_r8*(ustar(i,lt)**2)/grav) + diffk/(9.1_r8*ustar(i,lt))
             cvar(i,lt)  = vonkar/(log(z(i)/z0water))
             ustar(i,lt) = sqrt( cvar(i,lt)*uustar(i)/(1._r8 + 4.7_r8*ribn(i)) )
          endif
       endif
    enddo

    !-------------------------------------------------------------------------------------
    ! compute monin-obukhov length for unstable and stable conditions/ sublayer resistance
    !-------------------------------------------------------------------------------------
    do lt = beglt,endlt
       do i = 1,ncol
          if( fr_lnduse(i,lt) ) then
             hvar = (va(i)/0.74_r8) * (tha(i) - thg(i)) * (cvar(i,lt)**2)
             if( unstable(i) ) then                      ! unstable
                h = hvar*(1._r8 - (9.4_r8*ribn(i)/(1._r8 + 5.3_r8*b(i,lt))))
             else
                h = hvar/((1._r8+4.7_r8*ribn(i))**2)
             endif
             xmol(i,lt) = thg(i) * ustar(i,lt) * ustar(i,lt) / (vonkar * grav * h)
          endif
       enddo
    enddo

    !-------------------------------------------------------------------------------------
    ! psih
    !-------------------------------------------------------------------------------------
    do lt = beglt,endlt
       do i = 1,ncol
          if( fr_lnduse(i,lt) ) then
             if( xmol(i,lt) < 0._r8 ) then
                zovl = z(i)/xmol(i,lt)
                zovl = max( -1._r8,zovl )
                psih = exp( .598_r8 + .39_r8*log( -zovl ) - .09_r8*(log( -zovl ))**2 )
             else
                zovl = z(i)/xmol(i,lt)
                zovl = min( 1._r8,zovl )
                psih = -5._r8 * zovl
             endif
             dep_ra (i,lt,lchnk) = (vonkar - psih*cvar(i,lt))/(ustar(i,lt)*vonkar*cvar(i,lt))
             dep_rb (i,lt,lchnk) = (2._r8/(vonkar*ustar(i,lt))) * crb
          endif
       enddo
    enddo

    !-------------------------------------------------------------------------------------
    ! surface resistance : depends on both land type and species
    ! land types are computed seperately, then resistance is computed as average of values
    ! following wesely rc=(1/(rs+rm) + 1/rlu +1/(rdc+rcl) + 1/(rac+rgs))**-1
    !
    ! compute rsmx = 1/(rs+rm) : multiply by 3 if surface is wet
    !-------------------------------------------------------------------------------------
    rlux = 0._r8
    species_loop1 :  do ispec = 1,gas_pcnst
       if( has_dvel(ispec) ) then
          m = map_dvel(ispec)
          do lt = beglt,endlt
             do i = 1,ncol
                if( fr_lnduse(i,lt) ) then
                   sndx = index_season(i,lt)
                   if( ispec == o3_ndx .or. ispec == o3a_ndx .or. ispec == so2_ndx ) then
                      rmx = 0._r8
                   else
                      rmx = 1._r8/(heff(i,m)/3000._r8 + 100._r8*foxd(m))
                   end if
                   cts(i) = 1000._r8*exp( - tc(i) - 4._r8 )                 ! correction for frost
                   rgsx(i,lt,ispec) = cts(i) + 1._r8/((heff(i,m)/(1.e5_r8*rgss(sndx,lt))) + (foxd(m)/rgso(sndx,lt)))
                   !-------------------------------------------------------------------------------------
                   ! special case for H2 and CO;; CH4 is set ot a fraction of dv(H2)
                   !-------------------------------------------------------------------------------------
                   if( ispec == h2_ndx .or. ispec == co_ndx .or. ispec == ch4_ndx ) then
                      if( ispec == co_ndx ) then
                         fact_h2 = 1.0_r8
                      elseif ( ispec == h2_ndx ) then
                         fact_h2 = 0.5_r8
                      elseif ( ispec == ch4_ndx ) then
                         fact_h2 = 50.0_r8
                      end if
                      !-------------------------------------------------------------------------------------
                      ! no deposition on snow, ice, desert, and water
                      !-------------------------------------------------------------------------------------
                      if( lt == 1 .or. lt == 7 .or. lt == 8 .or. sndx == 4 ) then
                         rgsx(i,lt,ispec) = large_value
                      else
                         var_soilw = max( .1_r8,min( soilw(i),.3_r8 ) )
                         if( lt == 3 ) then
                            var_soilw = log( var_soilw )
                         end if
                         dv_soil_h2 = h2_c(lt) + var_soilw*(h2_b(lt) + var_soilw*h2_a(lt))
                         if( dv_soil_h2 > 0._r8 ) then
                            rgsx(i,lt,ispec) = fact_h2/(dv_soil_h2*1.e-4_r8)
                         end if
                      end if
                   end if
                   if( lt == 7 ) then
                      rclx(i,lt,ispec) = large_value
                      rsmx(i,lt,ispec) = large_value
                      rlux(i,lt,ispec) = large_value
                   else
                      rs = ri(sndx,lt)*crs(i)
                      if ( has_dew(i) .or. has_rain(i) ) then
                         dewm = 3._r8
                      else
                         dewm = 1._r8
                      end if
                      rsmx(i,lt,ispec) = (dewm*rs*drat(m) + rmx)
                      rclx(i,lt,ispec) = cts(i) + 1._r8/((heff(i,m)/(1.e5_r8*rcls(sndx,lt))) + (foxd(m)/rclo(sndx,lt)))
                      rlux(i,lt,ispec) = cts(i) + rlu(sndx,lt)/(1.e-5_r8*heff(i,m) + foxd(m))
                   end if
                end if
             end do
          end do
       end if
    end do species_loop1

    do lt = beglt,endlt
       if( lt /= 7 ) then
          do i = 1,ncol
             if( fr_lnduse(i,lt) ) then
                sndx = index_season(i,lt)
                !-------------------------------------------------------------------------------------
                ! 	... no effect if sfc_temp < O C
                !-------------------------------------------------------------------------------------
                if( sfc_temp(i) > tmelt ) then
                   if( has_dew(i) ) then
                      rlux_o3(i,lt)     = 3000._r8*rlu(sndx,lt)/(1000._r8 + rlu(sndx,lt))
                      if( o3_ndx > 0 ) then
                         rlux(i,lt,o3_ndx) = rlux_o3(i,lt)
                      endif
                   endif
                   if( has_rain(i) ) then
                      ! rlux(i,lt,o3_ndx) = 1./(1.e-3 + (1./(3.*rlu(sndx,lt))))
                      rlux_o3(i,lt)     = 3000._r8*rlu(sndx,lt)/(1000._r8 + 3._r8*rlu(sndx,lt))
                      if( o3_ndx > 0 ) then
                         rlux(i,lt,o3_ndx) = rlux_o3(i,lt)
                      endif
                      if( o3a_ndx > 0 ) then
                         rlux(i,lt,o3a_ndx) = rlux_o3(i,lt)
                      endif
                   endif
                endif

                if ( o3_ndx > 0 ) then
                   rclx(i,lt,o3_ndx) = cts(i) + rclo(index_season(i,lt),lt)
                   rlux(i,lt,o3_ndx) = cts(i) + rlux(i,lt,o3_ndx)
                endif
                
             endif
          enddo
       endif
    enddo

    species_loop2 : do ispec = 1,gas_pcnst
       m = map_dvel(ispec)
       if( has_dvel(ispec) ) then
          if( ispec /= o3_ndx .and. ispec /= o3a_ndx .and. ispec /= so2_ndx ) then
             do lt = beglt,endlt
                if( lt /= 7 ) then
                   do i = 1,ncol
                      if( fr_lnduse(i,lt) ) then
                         !-------------------------------------------------------------------------------------
                         ! no effect if sfc_temp < O C
                         !-------------------------------------------------------------------------------------
                         if( sfc_temp(i) > tmelt ) then
                            if( has_dew(i) ) then
                               rlux(i,lt,ispec) = 1._r8/((1._r8/(3._r8*rlux(i,lt,ispec))) &
                                    + 1.e-7_r8*heff(i,m) + foxd(m)/rlux_o3(i,lt))
                            endif
                         endif
                      endif
                   enddo
                endif
             enddo
          else if( ispec == so2_ndx ) then
             do lt = beglt,endlt
                if( lt /= 7 ) then
                   do i = 1,ncol
                      if( fr_lnduse(i,lt) ) then
                         !-------------------------------------------------------------------------------------
                         ! no effect if sfc_temp < O C
                         !-------------------------------------------------------------------------------------
                         if( sfc_temp(i) > tmelt ) then
                            if( qs(i) <= spec_hum(i) ) then
                               rlux(i,lt,ispec) = 100._r8
                            end if
                            if( has_rain(i) ) then
                               !                               rlux(i,lt,ispec) = 1./(2.e-4 + (1./(3.*rlu(index_season(i,lt),lt))))
                               rlux(i,lt,ispec) = 15._r8*rlu(index_season(i,lt),lt)/(5._r8 + 3.e-3_r8*rlu(index_season(i,lt),lt))
                            end if
                         end if
                         rclx(i,lt,ispec) = cts(i) + rcls(index_season(i,lt),lt)
                         rlux(i,lt,ispec) = cts(i) + rlux(i,lt,ispec)

                      endif
                   enddo
                endif
             enddo
             do i = 1,ncol
                if( fr_lnduse(i,1) .and. (has_dew(i) .or. has_rain(i)) ) then
                   rlux(i,1,ispec) = 50._r8
                endif
             enddo
          endif
       endif
    enddo species_loop2

    !-------------------------------------------------------------------------------------
    ! compute rc
    !-------------------------------------------------------------------------------------
    term(:ncol) = 1.e-2_r8 * pressure_10m(:ncol) / (r*tv(:ncol))
    species_loop3 : do ispec = 1,gas_pcnst
       if( has_dvel(ispec) ) then
          wrk(:) = 0._r8
          lt_loop: do lt = beglt,endlt
             do i = 1,ncol
                if (fr_lnduse(i,lt)) then
                   resc(i) = 1._r8/( 1._r8/rsmx(i,lt,ispec) + 1._r8/rlux(i,lt,ispec) &
                                   + 1._r8/(rdc(i) + rclx(i,lt,ispec)) &
                                   + 1._r8/(rac(index_season(i,lt),lt) + rgsx(i,lt,ispec)))

                   resc(i) = max( 10._r8,resc(i) )

                   lnd_frc(i) = lcl_frc_landuse(i,lt)
                endif
             enddo
             !-------------------------------------------------------------------------------------
             ! 	... compute average deposition velocity
             !-------------------------------------------------------------------------------------
             select case( solsym(ispec) )
             case( 'SO2' )
                if( lt == 7 ) then
                   where( fr_lnduse(:ncol,lt) )
                      ! assume no surface resistance for SO2 over water`
                      wrk(:) = wrk(:) + lnd_frc(:)/(dep_ra(:ncol,lt,lchnk) + dep_rb(:ncol,lt,lchnk)) 
                   endwhere
                else
                   where( fr_lnduse(:ncol,lt) )
                      wrk(:) = wrk(:) + lnd_frc(:)/(dep_ra(:ncol,lt,lchnk) + dep_rb(:ncol,lt,lchnk) + resc(:))
                   endwhere
                end if
             case default
                where( fr_lnduse(:ncol,lt) )
                   wrk(:ncol) = wrk(:ncol) + lnd_frc(:ncol)/(dep_ra(:ncol,lt,lchnk) + dep_rb(:ncol,lt,lchnk) + resc(:ncol))
                endwhere
             end select
          end do lt_loop
          dvel(:ncol,ispec) = wrk(:ncol) * scaling_to_cm_per_s
          dflx(:ncol,ispec) = term(:ncol) * dvel(:ncol,ispec) * mmr(:ncol,plev,ispec)
       endif

    enddo species_loop3

    if ( beglt > 1 ) return

  end subroutine drydep_xactive


  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  subroutine soilw_inti( ncfile, nlon_veg, nlat_veg, soilw_map )
    !------------------------------------------------------------------
    !	... read primary soil moisture table
    !------------------------------------------------------------------

    use time_manager,  only : get_calday

    implicit none

    !------------------------------------------------------------------
    !	... dummy args
    !------------------------------------------------------------------
    integer, intent(in) :: &
         nlon_veg, &
         nlat_veg
    real(r8), pointer :: soilw_map(:,:,:)
    character(len=*), intent(in) :: ncfile ! file name of netcdf file containing data

    !------------------------------------------------------------------
    !	... local variables
    !------------------------------------------------------------------
    integer :: gndx = 0
    integer :: nlat, &             ! # of lats in soilw file
               nlon                ! # of lons in soilw file
    integer :: i, ip, k, m
    integer :: j, jl, ju
    integer :: lev, day, ierr
    type(file_desc_t) :: piofile
    type(var_desc_t) :: vid
    
    integer :: dimid_lat, dimid_lon, dimid_time
    integer :: dates(12) = (/ 116, 214, 316, 415,  516,  615, &
                              716, 816, 915, 1016, 1115, 1216 /)

    character(len=shr_kind_cl) :: locfn

    !-----------------------------------------------------------------------
    !       ... open netcdf file
    !-----------------------------------------------------------------------
    call getfil (ncfile, locfn, 0)
    call cam_pio_openfile (piofile, trim(locfn), PIO_NOWRITE)

    !-----------------------------------------------------------------------
    !       ... get longitudes
    !-----------------------------------------------------------------------
    ierr = pio_inq_dimid( piofile, 'lon', dimid_lon )
    ierr = pio_inq_dimlen( piofile, dimid_lon, nlon )
    if( nlon /= nlon_veg ) then
       write(iulog,*) 'soilw_inti: soil and vegetation lons differ; ',nlon, nlon_veg
       call endrun
    end if
    !-----------------------------------------------------------------------
    !       ... get latitudes
    !-----------------------------------------------------------------------
    ierr = pio_inq_dimid( piofile, 'lat', dimid_lat )
    ierr = pio_inq_dimlen( piofile, dimid_lat, nlat )
    if( nlat /= nlat_veg ) then
       write(iulog,*) 'soilw_inti: soil and vegetation lats differ; ',nlat, nlat_veg
       call endrun
    end if
    !-----------------------------------------------------------------------
    !       ... set times (days of year)
    !-----------------------------------------------------------------------
    ierr = pio_inq_dimid( piofile, 'time', dimid_time )
    ierr = pio_inq_dimlen( piofile, dimid_time, ndays )
    if( ndays /= 12 ) then
       write(iulog,*) 'soilw_inti: dataset not a cyclical year'
       call endrun
    end if
    allocate( days(ndays),stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'soilw_inti: days allocation error = ',ierr
       call endrun
    end if
    do m = 1,min(12,ndays)
       days(m) = get_calday( dates(m), 0 )
    end do

    !------------------------------------------------------------------
    !	... allocate arrays
    !------------------------------------------------------------------
    allocate( soilw_map(nlon,nlat,ndays), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'soilw_inti: soilw_map allocation error = ',ierr
       call endrun
    end if

    !------------------------------------------------------------------
    !	... read in the soil moisture
    !------------------------------------------------------------------
    ierr = pio_inq_varid( piofile, 'SOILW', vid )
    ierr = pio_get_var( piofile, vid, soilw_map )
    !------------------------------------------------------------------
    !	... close file
    !------------------------------------------------------------------
    call pio_closefile( piofile )

  end subroutine soilw_inti
  
  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  subroutine chk_soilw( calday )
    !--------------------------------------------------------------------
    !	... check timing for ub values
    !--------------------------------------------------------------------

    use mo_constants, only : dayspy

    implicit none

    !--------------------------------------------------------------------
    !	... dummy args
    !--------------------------------------------------------------------
    real(r8), intent(in)    :: calday

    !--------------------------------------------------------------------
    !	... local variables
    !--------------------------------------------------------------------
    integer  ::  m, upper
    real(r8)     ::  numer, denom

    !--------------------------------------------------------
    !	... setup the time interpolation
    !--------------------------------------------------------
    if( calday < days(1) ) then
       next = 1
       last = ndays
    else
       if( days(ndays) < dayspy ) then
          upper = ndays
       else
          upper = ndays - 1
       end if
       do m = upper,1,-1
          if( calday >= days(m) ) then
             exit
          end if
       end do
       last = m
       next = mod( m,ndays ) + 1
    end if
    numer = calday - days(last)
    denom = days(next) - days(last)
    if( numer < 0._r8 ) then
       numer = dayspy + numer
    end if
    if( denom < 0._r8 ) then
       denom = dayspy + denom
    end if
    dels = max( min( 1._r8,numer/denom ),0._r8 )

  end subroutine chk_soilw

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  subroutine set_soilw( soilw, lchnk, calday )
    !--------------------------------------------------------------------
    !	... set the soil moisture
    !--------------------------------------------------------------------

    implicit none

    !--------------------------------------------------------------------
    !	... dummy args
    !--------------------------------------------------------------------
    real(r8), intent(inout) :: soilw(pcols)
    integer,  intent(in)    :: lchnk           ! chunk indice
    real(r8), intent(in)    :: calday


    integer :: i, ilon,ilat

    call chk_soilw( calday )

    soilw(:) = soilw_3d(:,last,lchnk) + dels *( soilw_3d(:,next,lchnk) - soilw_3d(:,last,lchnk))

  end subroutine set_soilw

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  function has_drydep( name )

    implicit none

    character(len=*), intent(in) :: name

    logical :: has_drydep
    integer :: i

    has_drydep = .false.

    do i=1,nddvels
       if ( trim(name) == trim(drydep_list(i)) ) then
         has_drydep = .true.
         exit
       endif
    enddo

  endfunction has_drydep

end module mo_drydep
