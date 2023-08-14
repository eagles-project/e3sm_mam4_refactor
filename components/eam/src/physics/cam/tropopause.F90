
! This module is used to diagnose the location of the tropopause. Multiple
! algorithms are provided, some of which may not be able to identify a
! tropopause in all situations. To handle these cases, an analytic
! definition and a climatology are provided that can be used to fill in
! when the original algorithm fails. The tropopause temperature and
! pressure are determined and can be output to the history file.
!
! These routines are based upon code in the WACCM chemistry module
! including mo_tropoause.F90 and llnl_set_chem_trop.F90. The code
! for the Reichler et al. [2003] algorithm is from:
!
!   http://www.gfdl.noaa.gov/~tjr/TROPO/tropocode.htm
!
! Author: Charles Bardeen
! Created: April, 2009

module tropopause
  !---------------------------------------------------------------
  ! ... variables for the tropopause module
  !---------------------------------------------------------------

  use shr_kind_mod,  only : r8 => shr_kind_r8
  use shr_const_mod, only : pi => shr_const_pi
  use ppgrid,        only : pcols, pver, begchunk, endchunk
  use cam_abortutils,    only : endrun
  use cam_logfile,   only : iulog
  use cam_history_support,   only : fillvalue
  use physics_types, only : physics_state
  use physconst,     only : cappa, rair, gravit
  use spmd_utils,    only : masterproc

  implicit none

  private
  
  public  :: tropopause_readnl, tropopause_init, tropopause_find, tropopause_output
  public  :: TROP_ALG_NONE, TROP_ALG_ANALYTIC, TROP_ALG_CLIMATE
  public  :: TROP_ALG_STOBIE, TROP_ALG_HYBSTOB, TROP_ALG_TWMO, TROP_ALG_WMO
  public  :: TROP_ALG_CPP
  public  :: NOTFOUND

  save

  ! These parameters define and enumeration to be used to define the primary
  ! and backup algorithms to be used with the tropopause_find() method. The
  ! backup algorithm is meant to provide a solution when the primary algorithm
  ! fail. The algorithms that can't fail are: TROP_ALG_ANALYTIC, TROP_ALG_CLIMATE
  ! and TROP_ALG_STOBIE.
  ! NOTE:  TROP_ALG_ANALYTIC, TROP_ALG_STOBIE, and TROP_ALG_WMO are not currently
  ! invoked anywhere in the code, so these algorithms have been removed from the source code.
  
  integer, parameter    :: TROP_ALG_NONE      = 1    ! Don't evaluate
  integer, parameter    :: TROP_ALG_ANALYTIC  = 2    ! Analytic Expression
  integer, parameter    :: TROP_ALG_CLIMATE   = 3    ! Climatology
  integer, parameter    :: TROP_ALG_STOBIE    = 4    ! Stobie Algorithm
  integer, parameter    :: TROP_ALG_TWMO      = 5    ! WMO Definition, Reichler et al. [2003]
  integer, parameter    :: TROP_ALG_WMO       = 6    ! WMO Definition
  integer, parameter    :: TROP_ALG_HYBSTOB   = 7    ! Hybrid Stobie Algorithm
  integer, parameter    :: TROP_ALG_CPP       = 8    ! Cold Point Parabolic
  
  integer, parameter    :: TROP_NALG          = 8    ! Number of Algorithms  
  character,parameter   :: TROP_LETTER(TROP_NALG) = (/ ' ', 'A', 'C', 'S', 'T', 'W', 'H', 'F' /)
                                                     ! unique identifier for output, don't use P

  ! These variables should probably be controlled by namelist entries.
  logical ,parameter    :: output_all         = .False.              ! output tropopause info from all algorithms
  integer ,parameter    :: default_primary    = TROP_ALG_TWMO        ! default primary algorithm
  integer ,parameter    :: default_backup     = TROP_ALG_CLIMATE     ! default backup algorithm

  ! Namelist variables
  character(len=256)    :: tropopause_climo_file = 'trop_climo'      ! absolute filepath of climatology file

  ! These variables are used to store the climatology data.
  integer, parameter :: days_in_year_real = 365._r8
  integer, parameter :: months_in_year = 12
  integer, parameter :: first_month = 1
  real(r8)              :: days(months_in_year)                      ! days in the climatology
  real(r8), pointer     :: tropp_p_loc(:,:,:)                        ! climatological tropopause pressures

  integer, parameter :: NOTFOUND = -1

  real(r8),parameter :: ALPHA  = 0.03_r8
    
  ! physical constants
  ! These constants are set in module variables rather than as parameters 
  ! to support the aquaplanet mode in which the constants have values determined
  ! by the experiment protocol
  real(r8) :: cnst_kap     ! = cappa
  real(r8) :: cnst_faktor  ! = -gravit/rair
  real(r8) :: cnst_ka1     ! = cnst_kap - 1._r8

!================================================================================================
contains
!================================================================================================

   ! Read namelist variables.
   subroutine tropopause_readnl(nlfile)

      use namelist_utils,  only: find_group_name
      use units,           only: getunit, freeunit
      use mpishorthand

      character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

      ! Local variables
      integer :: unitn, ierr
      character(len=*), parameter :: subname = 'tropopause_readnl'

      namelist /tropopause_nl/ tropopause_climo_file
      !-----------------------------------------------------------------------------

      if (masterproc) then
         unitn = getunit()
         open( unitn, file=trim(nlfile), status='old' )
         call find_group_name(unitn, 'tropopause_nl', status=ierr)
         if (ierr == 0) then
            read(unitn, tropopause_nl, iostat=ierr)
            if (ierr /= 0) then
               call endrun(subname // ':: ERROR reading namelist')
            endif
         endif
         close(unitn)
         call freeunit(unitn)
      endif

#ifdef SPMD
      ! Broadcast namelist variables
      call mpibcast(tropopause_climo_file, len(tropopause_climo_file), mpichar, 0, mpicom)
#endif

   end subroutine tropopause_readnl


  ! This routine is called during intialization and must be called before the
  ! other methods in this module can be used. Its main tasks are to read in the
  ! climatology from a file and to define the output fields. Much of this code
  ! is taken from mo_tropopause.
  subroutine tropopause_init()
  

    use ppgrid,        only: pver
    use cam_pio_utils, only: 
    use cam_history,   only: addfld, horiz_only, add_default


    implicit none

    ! define physical constants
    cnst_kap    = cappa
    cnst_faktor = -gravit/rair
    cnst_ka1    = cnst_kap - 1._r8

    ! Define the output fields.
    call addfld('TROP_P',          horiz_only,    'A',  'Pa', 'Tropopause Pressure', flag_xyfill=.True.)
    call addfld('TROP_T',           horiz_only,    'A',  'K', 'Tropopause Temperature', flag_xyfill=.True.)
    call addfld('TROP_Z',           horiz_only,    'A',  'm', 'Tropopause Height', flag_xyfill=.True.)
    call addfld('TROP_DZ',           (/ 'lev' /), 'A', 'm', 'Relative Tropopause Height')
    call addfld('TROP_PD', (/ 'lev' /), 'A', 'probability', 'Tropopause Probabilty')
    call addfld('TROP_FD', horiz_only,    'A', 'probability', 'Tropopause Found')
    
    call addfld('TROPP_P',          horiz_only,    'A',  'Pa', 'Tropopause Pressure (primary)', flag_xyfill=.True.)
    call addfld('TROPP_T',           horiz_only,    'A',  'K', 'Tropopause Temperature (primary)', flag_xyfill=.True.)
    call addfld('TROPP_Z',           horiz_only,    'A',  'm', 'Tropopause Height (primary)', flag_xyfill=.True.)
    call addfld('TROPP_DZ',         (/ 'lev' /),   'A', 'm', 'Relalive Tropopause Height (primary)')
    call addfld('TROPP_PD', (/ 'lev' /), 'A', 'probability', 'Tropopause Distribution (primary)')
    call addfld('TROPP_FD', horiz_only,    'A', 'probability', 'Tropopause Found (primary)')

    call addfld('TROPF_P',         horiz_only,  'A',  'Pa',         'Tropopause Pressure (cold point)',    flag_xyfill=.True.)
    call addfld('TROPF_T',         horiz_only,  'A',  'K',          'Tropopause Temperature (cold point)', flag_xyfill=.True.)
    call addfld('TROPF_Z',         horiz_only,  'A',  'm',          'Tropopause Height (cold point)',      flag_xyfill=.True.)
    call addfld('TROPF_DZ',        (/ 'lev' /),  'A', 'm',          'Relative Tropopause Height (cold point)', flag_xyfill=.True.)
    call addfld('TROPF_PD',        (/ 'lev' /), 'A', 'probability', 'Tropopause Distribution (cold point)')
    call addfld('TROPF_FD',        horiz_only,  'A', 'probability', 'Tropopause Found (cold point)')
    
    call addfld( 'hstobie_trop', (/ 'lev' /), 'I',   'fraction of model time', 'Lowest level with stratospheric chemsitry' )
    call addfld( 'hstobie_linoz', (/ 'lev' /), 'I',  'fraction of model time', 'Lowest possible Linoz level' )
    call addfld( 'hstobie_tropop', (/ 'lev' /), 'I', 'fraction of model time', &
         'Troposphere boundary calculated in chemistry' )

    call add_default('TROP_P', 1, ' ')
    call add_default('TROP_T', 1, ' ')
    call add_default('hstobie_linoz', 1, ' ')

    call tropopause_read_file()


  end subroutine tropopause_init
  

  subroutine tropopause_read_file
    !------------------------------------------------------------------
    ! ... initialize upper boundary values
    !------------------------------------------------------------------
    use interpolate_data,  only : lininterp_init, lininterp, interp_type, lininterp_finish
    use dyn_grid,     only : get_dyn_grid_parm
    use phys_grid,    only : get_ncols_p, get_rlat_all_p, get_rlon_all_p	
    use ioFileMod,    only : getfil
    use time_manager, only : get_calday
    use physconst,    only : pi
    use cam_pio_utils, only: cam_pio_openfile
    use pio,          only : file_desc_t, var_desc_t, pio_inq_dimid, pio_inq_dimlen, &
         pio_inq_varid, pio_get_var, pio_closefile, pio_nowrite

    !------------------------------------------------------------------
    ! ... local variables
    !------------------------------------------------------------------
    integer :: i, j, n
    integer :: ierr
    type(file_desc_t) :: pio_id
    integer :: dimid
    type(var_desc_t) :: vid
    integer :: nlon, nlat, ntimes
    integer :: start(3)
    integer :: count(3)
    integer, parameter :: dates(months_in_year) = (/ 116, 214, 316, 415,  516,  615, &
         716, 816, 915, 1016, 1115, 1216 /)
    integer :: plon, plat
    type(interp_type) :: lon_wgts, lat_wgts
    real(r8), allocatable :: tropp_p_in(:,:,:)
    real(r8), allocatable :: lat(:)
    real(r8), allocatable :: lon(:)
    real(r8) :: to_lats(pcols), to_lons(pcols)
    real(r8), parameter :: d2r=pi/180._r8, zero=0._r8, twopi=pi*2._r8
    character(len=256) :: locfn
    integer  :: c, ncols


    plon = get_dyn_grid_parm('plon')
    plat = get_dyn_grid_parm('plat')


    !-----------------------------------------------------------------------
    !       ... open netcdf file
    !-----------------------------------------------------------------------
    call getfil (tropopause_climo_file, locfn, 0)
    call cam_pio_openfile(pio_id, trim(locfn), PIO_NOWRITE)

    !-----------------------------------------------------------------------
    !       ... get time dimension
    !-----------------------------------------------------------------------
    ierr = pio_inq_dimid( pio_id, 'time', dimid )
    ierr = pio_inq_dimlen( pio_id, dimid, ntimes )
    if( ntimes /= months_in_year )then
       write(iulog,*) 'tropopause_init: number of months = ',ntimes,'; expecting 12'
       call endrun
    endif
    !-----------------------------------------------------------------------
    !       ... get latitudes
    !-----------------------------------------------------------------------
    ierr = pio_inq_dimid( pio_id, 'lat', dimid )
    ierr = pio_inq_dimlen( pio_id, dimid, nlat )
    allocate( lat(nlat), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'tropopause_init: lat allocation error = ',ierr
       call endrun
    endif
    ierr = pio_inq_varid( pio_id, 'lat', vid )
    ierr = pio_get_var( pio_id, vid, lat )
    lat(:nlat) = lat(:nlat) * d2r
    !-----------------------------------------------------------------------
    !       ... get longitudes
    !-----------------------------------------------------------------------
    ierr = pio_inq_dimid( pio_id, 'lon', dimid )
    ierr = pio_inq_dimlen( pio_id, dimid, nlon )
    allocate( lon(nlon), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'tropopause_init: lon allocation error = ',ierr
       call endrun
    endif
    ierr = pio_inq_varid( pio_id, 'lon', vid )
    ierr = pio_get_var( pio_id, vid, lon )
    lon(:nlon) = lon(:nlon) * d2r

    !------------------------------------------------------------------
    !  ... allocate arrays
    !------------------------------------------------------------------
    allocate( tropp_p_in(nlon,nlat,ntimes), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'tropopause_init: tropp_p_in allocation error = ',ierr
       call endrun
    endif
    !------------------------------------------------------------------
    !  ... read in the tropopause pressure
    !------------------------------------------------------------------
    ierr = pio_inq_varid( pio_id, 'trop_p', vid )
    start = (/ 1, 1, 1 /)
    count = (/ nlon, nlat, ntimes /)
    ierr = pio_get_var( pio_id, vid, start, count, tropp_p_in )

    !------------------------------------------------------------------
    !  ... close the netcdf file
    !------------------------------------------------------------------
    call pio_closefile( pio_id )

    !--------------------------------------------------------------------
    !  ... regrid
    !--------------------------------------------------------------------

    allocate( tropp_p_loc(pcols,begchunk:endchunk,ntimes), stat=ierr )

    if( ierr /= 0 ) then
      write(iulog,*) 'tropopause_init: tropp_p_loc allocation error = ',ierr
      call endrun
    endif

    do c=begchunk,endchunk
       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, pcols, to_lats)
       call get_rlon_all_p(c, pcols, to_lons)
       call lininterp_init(lon, nlon, to_lons, ncols, 2, lon_wgts, zero, twopi)
       call lininterp_init(lat, nlat, to_lats, ncols, 1, lat_wgts)
       do n=1,ntimes
          call lininterp(tropp_p_in(:,:,n), nlon, nlat, tropp_p_loc(1:ncols,c,n), ncols, lon_wgts, lat_wgts)    
       enddo
       call lininterp_finish(lon_wgts)
       call lininterp_finish(lat_wgts)
    enddo
    deallocate(lon)
    deallocate(lat)
    deallocate(tropp_p_in)

    !--------------------------------------------------------
    ! ... initialize the monthly day of year times
    !--------------------------------------------------------

    do n = 1,months_in_year
       days(n) = get_calday( dates(n), 0 )
    enddo
    if (masterproc) then
       write(iulog,*) 'tropopause_init : days'
       write(iulog,'(1p,5g15.8)') days(:)
    endif

  end subroutine tropopause_read_file
  
  ! Read the tropopause pressure in from a file containging a climatology. The
  ! data is interpolated to the current dat of year and latitude.
  !
  ! NOTE: The data is read in during tropopause_init and stored in the module
  ! variable trop
  subroutine tropopause_climate(lchnk,ncol,pmid,pint,temp,zm,zi,    &  ! in
             tropLev,tropP,tropT,tropZ)   ! inout
    use time_manager, only : get_curr_calday
    use mam_support,  only : min_max_bound

    implicit none

    integer,           intent(in)       :: lchnk  ! chunk identifier
    integer,           intent(in)       :: ncol ! number of columns in the chunk
    real(r8),          intent(in)       :: pmid(:,:)  ! midpoint pressure [Pa]
    real(r8),          intent(in)       :: pint(:,:)  ! top interface pressure [Pa]
    real(r8),          intent(in)       :: temp(:,:)  ! temperature [K]
    real(r8),          intent(in)       :: zm(:,:)    ! midpoint geopotential height above the surface [m]
    real(r8),          intent(in)       :: zi(:,:)    ! interface geopotential height above sfc [m]
    integer,            intent(inout)  :: tropLev(pcols)            ! tropopause level index   
    real(r8), optional, intent(inout)  :: tropP(pcols)              ! tropopause pressure [Pa]   
    real(r8), optional, intent(inout)  :: tropT(pcols)              ! tropopause temperature [K]
    real(r8), optional, intent(inout)  :: tropZ(pcols)              ! tropopause height [m]
 
    ! Local Variables
    integer       :: icol                     ! index of column
    integer       :: kk                       ! index of vertical level
    integer       :: mn                       ! index of month
    real(r8)      :: tP                       ! tropopause pressure (Pa)
    real(r8)      :: calday                   ! day of year including fraction
    real(r8)      :: dels
    integer       :: last
    integer       :: next
    integer       :: tropLevVal               ! tropLevVal at a particular column
    integer       :: tropLevValp1             ! tropLevVal at a particular column at level + 1
    integer       :: tropLevValm1             ! tropLevVal at a particular column at level - 1
 


    ! If any columns remain to be indentified, the nget the current
    ! day from the calendar.
    
    if (any(tropLev == NOTFOUND)) then
    
      ! Determine the calendar day.
      calday = get_curr_calday()
      
      !--------------------------------------------------------
      ! ... setup the time interpolation
      !--------------------------------------------------------
      if( calday < days(first_month) ) then
        next = first_month
        last = months_in_year
        dels = (days_in_year_real + calday - days(months_in_year)) / (days_in_year_real + days(1) - days(months_in_year))
      else if( calday >= days(months_in_year) ) then
        next = first_month
        last = months_in_year
        dels = (calday - days(months_in_year)) / (days_in_year_real + days(first_month) - days(months_in_year))
      else
        do mn = months_in_year-1,first_month,-1
           if( calday >= days(mn) ) then
              exit
           endif
        enddo
        last = mn
        next = mn + 1
        dels = (calday - days(mn)) / (days(mn+1) - days(mn))
      endif
      
      dels = min_max_bound(0._r8,1._r8,dels)  

      ! Iterate over all of the columns.
      do icol = 1, ncol
       
        ! Skip column in which the tropopause has already been found.
        if (tropLev(icol) == NOTFOUND) then
        
        !--------------------------------------------------------
        ! ... get tropopause level from climatology
        !--------------------------------------------------------
          ! Interpolate the tropopause pressure.
          !
          ! NOTE FROM FORTRAN REFACTORING:  tropp_p_loc has location-based data as a function of icol and lchnk.
          ! C++ ported code needs to read in the data and use the data based on the column it is operating on.
          ! END NOTE
          tP = tropp_p_loc(icol,lchnk,last) &
            + dels * (tropp_p_loc(icol,lchnk,next) - tropp_p_loc(icol,lchnk,last))
                
          ! Find the associated level.
          do kk = pver, 2, -1
            if (tP >= pint(icol, kk)) then
              tropLevVal = kk
              tropLevValp1 = kk + 1
              tropLevValm1 = kk - 1
              exit
            endif
          enddo
          tropLev(icol) = tropLevVal

          ! Return the optional outputs
          if (present(tropP)) tropP(icol) = tP
          
          if (present(tropT)) then
            tropT(icol) = tropopause_interpolateT(pmid(icol,tropLevVal),pmid(icol,tropLevValp1),  &
                           pmid(icol,tropLevValm1),temp(icol,tropLevVal),temp(icol,tropLevValp1), &
                           temp(icol,tropLevValm1),tropLevVal,tP)
          endif

          if (present(tropZ)) then
            tropZ(icol) = tropopause_interpolateZ(pmid(icol,tropLevVal),pint(icol,tropLevVal),    & 
                           pint(icol,tropLevValp1),zm(icol,tropLevVal),zi(icol,tropLevVal),       &
                           zi(icol,tropLevValp1),tP)
          endif
        endif
      enddo
    endif        

    return    
  end subroutine tropopause_climate
  
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine tropopause_hybridstobie(lchnk,ncol,pmid,temp,zm,    &  ! in
             tropLev,tropP,tropT,tropZ)    ! inout
    use cam_history,  only : outfld
  
    !-----------------------------------------------------------------------
    ! Originally written by Philip Cameron-Smith, LLNL
    !
    !   Stobie-Linoz hybrid: the highest altitude of 
    !          a) Stobie algorithm, or 
    !          b) minimum Linoz pressure.
    !
    ! NOTE: the ltrop(i) gridbox itself is assumed to be a STRATOSPHERIC gridbox.
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    !        ... Local variables
    !-----------------------------------------------------------------------
    
    implicit none

    integer,           intent(in)       :: lchnk  ! chunk identifier
    integer,           intent(in)       :: ncol ! number of columns in the chunk
    real(r8),          intent(in)       :: pmid(:,:)  ! midpoint pressure [Pa]
    real(r8),          intent(in)       :: temp(:,:)  ! temperature [K]
    real(r8),          intent(in)       :: zm(:,:)    ! midpoint geopotential height above the surface [m]
    integer,            intent(inout)   :: tropLev(pcols)             ! tropopause level index   
    real(r8), optional, intent(inout)   :: tropP(pcols)               ! tropopause pressure [Pa]  
    real(r8), optional, intent(inout)   :: tropT(pcols)               ! tropopause temperature [K]
    real(r8), optional, intent(inout)   :: tropZ(pcols)               ! tropopause height [m]
    
    real(r8),parameter  ::  min_Stobie_Pressure= 40.E2_r8 !For case 2 & 4.  [Pa]
    real(r8),parameter  ::  max_Linoz_Pressure =208.E2_r8 !For case     4.  [Pa]

    integer      :: icol, kk  ! column index, vertical level index
    real(r8)     :: stobie_min, shybrid_temp      !temporary variable for case 2 & 3.
    integer      :: ltrop_linoz(pcols)            !Lowest possible Linoz vertical level
    integer      :: ltrop_trop(pcols)             !Tropopause level for hybrid case.
    logical      :: ltrop_linoz_set               !Flag that lowest linoz level already found.
    real(r8)     :: trop_output(pcols,pver)        !For output purposes only.
    real(r8)     :: trop_linoz_output(pcols,pver)  !For output purposes only.
    real(r8)     :: trop_trop_output(pcols,pver)   !For output purposes only.

    ltrop_linoz(:) = 1  ! Initialize to default value.
    ltrop_trop(:) = 1   ! Initialize to default value.

    LOOP_COL4: do icol=1,ncol

       ! Skip column in which the tropopause has already been found.
       not_found: if (tropLev(icol) == NOTFOUND) then

          stobie_min = 1.e10_r8    ! An impossibly large number
          ! NOTE FROM FORTRAN REFACTORING:  above should be set to huge(stobie_min),
          ! but changing might effect BFB convergence of our tests.
          ! END NOTE
          ltrop_linoz_set = .FALSE.
          LOOP_LEV: do kk=pver,1,-1
             IF (pmid(icol,kk) < min_stobie_pressure) cycle
             shybrid_temp = ALPHA * temp(icol,kk) - Log10(pmid(icol,kk))
             !PJC_NOTE: the units of pmid won't matter, because it is just an additive offset.
             IF (shybrid_temp<stobie_min) then 
                ltrop_trop(icol)=kk     
                stobie_min = shybrid_temp
             ENDIF
             IF (pmid(icol,kk) < max_Linoz_pressure .AND. .NOT. ltrop_linoz_set) THEN
                ltrop_linoz(icol) = kk
                ltrop_linoz_set = .TRUE.
             ENDIF
          enddo LOOP_LEV

          tropLev(icol) = MIN(ltrop_trop(icol),ltrop_linoz(icol))

          if (present(tropP)) then
             tropP(icol) = pmid(icol,tropLev(icol))
          endif
          if (present(tropT)) then
             tropT(icol) = temp(icol,tropLev(icol))
          endif
          if (present(tropZ)) then
             tropZ(icol) = zm(icol,tropLev(icol))
          endif

       endif not_found

    enddo LOOP_COL4

    trop_output(:,:)=0._r8
    trop_linoz_output(:,:)=0._r8
    trop_trop_output(:,:)=0._r8
    do icol=1,ncol
       trop_output(icol,tropLev(icol))=1._r8
       trop_linoz_output(icol,ltrop_linoz(icol))=1._r8
       trop_trop_output(icol,ltrop_trop(icol))=1._r8
    enddo

    call outfld( 'hstobie_trop',   trop_output(:ncol,:),       ncol, lchnk )
    call outfld( 'hstobie_linoz',  trop_linoz_output(:ncol,:), ncol, lchnk )
    call outfld( 'hstobie_tropop', trop_trop_output(:ncol,:),  ncol, lchnk )

  endsubroutine tropopause_hybridstobie
  
  ! This routine is an implementation of Reichler et al. [2003] done by
  ! Reichler and downloaded from his web site. Minimal modifications were
  ! made to have the routine work within the CAM framework (i.e. using
  ! CAM constants and types).
  !
  ! NOTE: I am not a big fan of the goto's and multiple returns in this
  ! code, but for the moment I have left them to preserve as much of the
  ! original and presumably well tested code as possible.
  ! UPDATE: The most "obvious" substitutions have been made to replace
  ! goto/return statements with cycle/exit. The structure is still
  ! somewhat tangled.
  ! UPDATE 2: "gamma" renamed to "gam" in order to avoid confusion
  ! with the Fortran 2008 intrinsic. "level" argument removed because
  ! a physics column is not contiguous, so using explicit dimensions
  ! will cause the data to be needlessly copied.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! determination of tropopause height from gridded temperature data
  !
  ! reference: Reichler, T., M. Dameris, and R. Sausen (2003)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine twmo(temp1d, pmid1d, plimu, pliml, gam, &  ! in
             trp)   ! out

    real(r8), intent(in), dimension(:)      :: temp1d   !  temperature in column [K]
    real(r8), intent(in), dimension(:)      :: pmid1d   !  midpoint pressure in column [Pa]
    real(r8), intent(in)                    :: plimu    ! upper limit of tropopause pressure [Pa]
    real(r8), intent(in)                    :: pliml    ! lower limit of tropopause pressure [Pa]
    real(r8), intent(in)                    :: gam      ! lapse rate to indicate tropopause [K/m]
    real(r8), intent(out)                   :: trp      ! tropopause pressure [Pa]
    
    real(r8), parameter                     :: deltaz = 2000.0_r8   ! [m]

    real(r8)                                :: pmk, pmk2
    real(r8)                                :: pm, pm2      ! mean pressure [Pa]
    real(r8)                                :: tm, tm2      ! mean temperature [K]
    real(r8)                                :: dtdz     ! temperature lapse rate vs. height [K/m]
    real(r8)                                :: ag, bg
    real(r8)                                :: ptph     ! pressure at tropopause height [Pa]
    real(r8)                                :: pm0, pmk0
    real(r8)                                :: dtdz0, dtdz2, asum, aquer    ! [K/m]
    real(r8)                                :: p2km     ! pressure at 2 km above tropopause height [Pa]
    integer                                 :: level   ! size of temp1d array
    integer                                 :: icount, jj
    integer                                 :: kk      ! vertical level index
!! BJG below
    real(r8)                                :: a1, b1
    real(r8)                                :: dtdp     ! temperature lapse rate vs pressure [K/Pa]
! BJG end
    trp=-99.0_r8                           ! negative means not valid
    
    ! initialize start level
    ! dt/dz
    level = size(temp1d)
    pmk= .5_r8 * (pmid1d(level-1)**cnst_kap+pmid1d(level)**cnst_kap)
    pm = pmk**(1/cnst_kap)               
!! BJG    call get_dtdz(pm,pmk,pmid1d(level-1),pmid1d(level),temp1d(level-1),temp1d(level),  & ! in
!! BJG         dtdz,tm)  ! out
      a1 = (temp1d(level-1)-temp1d(level))/(pmid1d(level-1)**cnst_kap-pmid1d(level)**cnst_kap)
      b1 = temp1d(level)-(a1*pmid1d(level)**cnst_kap)
      tm = a1 * pmk + b1               
      dtdp = a1 * cnst_kap * (pm**cnst_ka1)
      dtdz = cnst_faktor*dtdp*pm/tm
!!BJG end

    main_loop: do kk=level-1,2,-1
      pm0 = pm
      pmk0 = pmk
      dtdz0  = dtdz
    
      ! dt/dz
      pmk= .5_r8 * (pmid1d(kk-1)**cnst_kap+pmid1d(kk)**cnst_kap)
      pm = pmk**(1/cnst_kap)               
!!BJG      call get_dtdz(pm,pmk,pmid1d(kk-1),pmid1d(kk),temp1d(kk-1),temp1d(kk),   & ! in 
!!BJG           dtdz,tm)  ! out
      a1 = (temp1d(kk-1)-temp1d(kk))/(pmid1d(kk-1)**cnst_kap-pmid1d(kk)**cnst_kap)
      b1 = temp1d(kk)-(a1*pmid1d(kk)**cnst_kap)
      tm = a1 * pmk + b1               
      dtdp = a1 * cnst_kap * (pm**cnst_ka1)
      dtdz = cnst_faktor*dtdp*pm/tm
!!BJG end
      ! dt/dz valid?
      if (dtdz<=gam)   cycle main_loop    ! no, dt/dz < -2 K/km
      if (pm>plimu)   cycle main_loop    ! no, too low
  
      ! dtdz is valid, calculate tropopause pressure
      if (dtdz0<gam) then
        ag = (dtdz-dtdz0) / (pmk-pmk0)     
        bg = dtdz0 - (ag * pmk0)          
        ptph = exp(log((gam-bg)/ag)/cnst_kap)
      else
        ptph = pm
      endif
  
      if (ptph<pliml) cycle main_loop
      if (ptph>plimu) cycle main_loop
  
      ! 2nd test: dtdz above 2 km must not exceed gam
      p2km = ptph + deltaz*(pm/tm)*cnst_faktor     ! p at ptph + 2km
      asum = 0.0_r8                                ! dtdz above
      icount = 0                                   ! number of levels above
  
      ! test until apm < p2km
      in_loop: do jj=kk,2,-1
    
        pmk2 = .5_r8 * (pmid1d(jj-1)**cnst_kap+pmid1d(jj)**cnst_kap) ! p mean ^kappa
        pm2 = pmk2**(1/cnst_kap)                           ! p mean
        if(pm2>ptph) cycle in_loop            ! doesn't happen
        if(pm2<p2km) exit in_loop             ! ptropo is valid

!!BJG        call get_dtdz(pm2,pmk2,pmid1d(jj-1),pmid1d(jj),temp1d(jj-1),temp1d(jj),  & ! in
!!BJG             dtdz2,tm2) ! out
        a1 = (temp1d(jj-1)-temp1d(jj))                     ! a
        a1 = a1/(pmid1d(jj-1)**cnst_kap-pmid1d(jj)**cnst_kap)
        b1 = temp1d(jj)-(a1*pmid1d(jj)**cnst_kap)          ! b
        tm2 = a1 * pmk2 + b1                     ! T mean
        dtdp2 = a1 * cnst_kap * (pm2**(cnst_kap-1))  ! dt/dp
        dtdz2 = cnst_faktor*dtdp2*pm2/tm2
!!BJG end
        asum = asum+dtdz2
        icount = icount+1
        aquer = asum/float(icount)               ! dt/dz mean
   
        ! discard ptropo ?
        if (aquer.le.gam) cycle main_loop      ! dt/dz above < gam
    
      enddo in_loop  ! test next level
    
      trp = ptph
      exit main_loop
    enddo main_loop
    
  end subroutine twmo
  
  subroutine get_dtdz(pm,pmk,pmid1d_up,pmid1d_down,temp1d_up,temp1d_down,    & ! in
             dtdz,tm)   ! out

    implicit none

    real(r8), intent(in)                    :: pm      ! mean pressure [Pa]
    real(r8), intent(in)                    :: pmk
    real(r8), intent(in)                    :: pmid1d_up     !  midpoint pressure in column at upper level [Pa]
    real(r8), intent(in)                    :: pmid1d_down   !  midpoint pressure in column at lower level [Pa]
    real(r8), intent(in)                    :: temp1d_up     !  temperature in column at upper level [K]
    real(r8), intent(in)                    :: temp1d_down   !  temperature in column at lower level [K]
    real(r8), intent(out)                   :: dtdz     ! temperature lapse rate vs. height [K/m]
    real(r8), intent(out)                   :: tm       ! mean temperature [K] -- needed to find pressure at trop + 2 km

    real(r8)                                :: a1, b1
    real(r8)                                :: dtdp     ! temperature lapse rate vs pressure [K/Pa]


    a1 = (temp1d_up-temp1d_down)/(pmid1d_up**cnst_kap-pmid1d_down**cnst_kap)
    b1 = temp1d_down-(a1*pmid1d_down**cnst_kap)
    tm = a1 * pmk + b1
    dtdp = a1 * cnst_kap * (pm**cnst_ka1)
    dtdz = cnst_faktor*dtdp*pm/tm

  end subroutine get_dtdz

  ! This routine uses an implementation of Reichler et al. [2003] done by
  ! Reichler and downloaded from his web site. This is similar to the WMO
  ! routines, but is designed for GCMs with a coarse vertical grid.
  subroutine tropopause_twmo(ncol,pmid,pint,temp,zm,zi,   &  ! in
             tropLev,tropP,tropT,tropZ)  ! inout
    implicit none

    integer,           intent(in)       :: ncol ! number of columns in the chunk
    real(r8),          intent(in)       :: pmid(:,:)  ! midpoint pressure [Pa]
    real(r8),          intent(in)       :: pint(:,:)  ! top interface pressure [Pa]
    real(r8),          intent(in)       :: temp(:,:)  ! temperature [K]
    real(r8),          intent(in)       :: zm(:,:)    ! midpoint geopotential height above the surface [m]
    real(r8),          intent(in)       :: zi(:,:)    ! interface geopotential height above sfc [m]
    integer,            intent(inout)  :: tropLev(pcols)            ! tropopause level index   
    real(r8), optional, intent(inout)  :: tropP(pcols)              ! tropopause pressure [Pa]   
    real(r8), optional, intent(inout)  :: tropT(pcols)              ! tropopause temperature [K]
    real(r8), optional, intent(inout)  :: tropZ(pcols)              ! tropopause height [m]

    ! Local Variables 
    real(r8), parameter     :: gam    = -0.002_r8         ! lapse rate to indicate tropopause [K/m]
    real(r8), parameter     :: plimu    = 45000._r8       ! upper limit of tropopause pressure [Pa]
    real(r8), parameter     :: pliml    = 7500._r8        ! lower limit of tropopause pressure [Pa]
     
    integer                 :: icol                     ! column index     
    integer                 :: kk                       ! vertical level index
    real(r8)                :: tP                       ! tropopause pressure [Pa]
    integer       :: tropLevVal               ! tropLevVal at a particular column
    integer       :: tropLevValp1             ! tropLevVal at a particular column at level + 1
    integer       :: tropLevValm1             ! tropLevVal at a particular column at level - 1


    ! Iterate over all of the columns.
    do icol = 1, ncol
     
      ! Skip column in which the tropopause has already been found.
      if (tropLev(icol) == NOTFOUND) then

        ! Use the routine from Reichler.
        call twmo(temp(icol, :), pmid(icol, :), plimu, pliml, gam, &  ! in
             tP)  ! out
     
        ! if successful, store of the results and find the level and temperature.
        if (tP > 0) then
        
          ! Find the associated level.
          do kk = pver, 2, -1
            if (tP >= pint(icol, kk)) then
              tropLevVal = kk
              tropLevValp1 = kk + 1
              tropLevValm1 = kk - 1
              exit
            endif
          enddo
          tropLev(icol) = tropLevVal
          
          ! Return the optional outputs
          if (present(tropP)) tropP(icol) = tP
          
          if (present(tropT)) then
            tropT(icol) = tropopause_interpolateT(pmid(icol,tropLevVal),pmid(icol,tropLevValp1),   &
                            pmid(icol,tropLevValm1),temp(icol,tropLevVal),temp(icol,tropLevValp1), &
                            temp(icol,tropLevValm1),tropLevVal,tP)
          endif

          if (present(tropZ)) then
            tropZ(icol) = tropopause_interpolateZ(pmid(icol,tropLevVal),pint(icol,tropLevVal),     & 
                           pint(icol,tropLevValp1),zm(icol,tropLevVal),zi(icol,tropLevVal),        &
                           zi(icol,tropLevValp1),tP)
          endif
        endif
      endif
    enddo
    
    return
  end subroutine tropopause_twmo
  
  ! This routine searches for the cold point tropopause, and uses a parabolic
  ! fit of the coldest point and two adjacent points to interpolate the cold point
  ! between model levels.
  subroutine tropopause_cpp(ncol,pmid,temp,zm,    &  ! in
             tropLev,tropP,tropT,tropZ)   ! inout
    implicit none

    integer,           intent(in)       :: ncol ! number of columns in the chunk
    real(r8),          intent(in)       :: pmid(:,:)  ! midpoint pressure [Pa]
    real(r8),          intent(in)       :: temp(:,:)  ! temperature [K]
    real(r8),          intent(in)       :: zm(:,:)    ! midpoint geopotential height above the surface [m]
    integer,            intent(inout)  :: tropLev(pcols)            ! tropopause level index   
    real(r8), optional, intent(inout)  :: tropP(pcols)              ! tropopause pressure [Pa]   
    real(r8), optional, intent(inout)  :: tropT(pcols)              ! tropopause temperature [K]
    real(r8), optional, intent(inout)  :: tropZ(pcols)              ! tropopause height [m]

    ! Local Variables 
    real(r8), parameter    :: ztrop_low   = 5000._r8        ! lowest tropopause level allowed [m]
    real(r8), parameter    :: ztrop_high  = 25000._r8       ! highest tropopause level allowed [m]

    integer                 :: icol
    integer                 :: kk, firstk, lastk
    integer                 :: k2
    real(r8)                :: tZ                           ! tropopause height [m]
    real(r8)                :: tmin
    real(r8)                :: f0, f1, f2
    real(r8)                :: x0, x1, x2
    real(r8)                :: c0, c1, c2
    real(r8)                :: aterm, bterm, cterm
    integer       :: tropLevVal               ! tropLevVal at a particular column
    integer       :: tropLevValp1             ! tropLevVal at a particular column at level + 1
    integer       :: tropLevValm1             ! tropLevVal at a particular column at level - 1

    ! Iterate over all of the columns.
    do icol = 1, ncol
     
      firstk = 0
      lastk  = pver+1

      ! Skip column in which the tropopause has already been found.
      if (tropLev(icol) == NOTFOUND) then
        tmin = 1e6_r8
     ! NOTE FROM FORTRAN REFACTORING:  above should be set to huge(tmin), 
     ! but changing might effect BFB convergence
     ! END NOTE
        kloop: do kk = pver-1, 2, -1
         
          ! Skip levels below the minimum and stop if nothing is found
          ! before the maximum.
          if (zm(icol, kk) < ztrop_low) then
            firstk = kk
            cycle kloop
          elseif (zm(icol, kk) > ztrop_high) then
            lastk = kk
            exit kloop
          endif
          
          ! Find the coldest point
          if (temp(icol, kk) < tmin) then
            tropLevVal   = kk
            tropLevValp1 = kk + 1
            tropLevValm1 = kk - 1
            tmin = temp(icol,kk)
          endif
        enddo kloop
        tropLev(icol) = tropLevVal

        ! If the minimum is at the edge of the search range, then don't
        ! consider this to be a minima
        if ((tropLev(icol) >= (firstk-1)) .or. (tropLev(icol) <= (lastk+1))) then
          tropLev(icol) = NOTFOUND
        else

          ! If returning P, Z, or T, then do a parabolic fit using the
          ! cold point and it its 2 surrounding points to interpolate
          ! between model levels.
          if (present(tropP) .or. present(tropZ) .or. present(tropT)) then
            f0 = temp(icol, tropLevValm1)
            f1 = temp(icol, tropLevVal)
            f2 = temp(icol, tropLevValp1)

            x0 = zm(icol, tropLevValm1)
            x1 = zm(icol, tropLevVal)
            x2 = zm(icol, tropLevValp1)

            c0 = (x0-x1)*(x0-x2)
            c1 = (x1-x0)*(x1-x2)
            c2 = (x2-x0)*(x2-x1)

            ! Determine the quadratic coefficients of:
            !   T = a * z^2 - b*z + c
            aterm = (f0/c0 + f1/c1 + f2/c2)
            bterm = (f0/c0*(x1+x2) + f1/c1*(x0+x2) + f2/c2*(x0+x1))
            cterm = f0/c0*x1*x2 + f1/c1*x0*x2 + f2/c2*x0*x1

            ! Find the altitude of the minimum temperature
            tZ = 0.5_r8 * bterm / aterm
            
            ! The fit should be between the upper and lower points,
            ! so skip the point if the fit fails.
            if ((tZ >= x0) .or. (tZ <= x2)) then
              tropLev(icol) = NOTFOUND
            else
              ! Return the optional outputs
              if (present(tropP)) then

                tropP(icol) = tropopause_interpolateP(pmid(icol,tropLevVal),pmid(icol,tropLevValp1),  &
                                   pmid(icol,tropLevValm1),zm(icol,tropLevVal),zm(icol,tropLevValp1), &
                                   zm(icol,tropLevValm1),tropLev(icol),tZ)
              endif

              if (present(tropT)) then
                tropT(icol) = aterm * tZ*tZ - bterm*tZ + cterm
              endif

              if (present(tropZ)) then
                tropZ(icol) = tZ
              endif

            endif  !   if ((tZ >= x0) .or. (tZ <= x2))

          endif  !  if (present(tropP) .or. present(tropZ) .or. present(tropT)) 

        endif   ! if ((tropLev(icol) >= (firstk-1)) .or. (tropLev(icol) <= (lastk+1))) 

      endif    ! if (tropLev(icol) == NOTFOUND)

    enddo
    
    return
  end subroutine tropopause_cpp
  
  ! Searches all the columns in the chunk and attempts to identify the tropopause.
  ! Two routines can be specifed, a primary routine which is tried first and a
  ! backup routine which will be tried only if the first routine fails. If the
  ! tropopause can not be identified by either routine, then a NOTFOUND is returned
  ! for the tropopause level, temperature and pressure.
    subroutine tropopause_find(lchnk,ncol,pmid,pint,temp,zm,zi,   &   ! in
               tropLev,tropP,tropT,tropZ,    &   ! out / optional out
               primary,backup)   !  optional in

    implicit none

    integer,           intent(in)       :: lchnk
    integer,           intent(in)       :: ncol
    real(r8),          intent(in)       :: pmid(:,:)  ! midpoint pressure [Pa]
    real(r8),          intent(in)       :: pint(:,:)  ! top interface pressure [Pa]
    real(r8),          intent(in)       :: temp(:,:)  ! temperature [K]
    real(r8),          intent(in)       :: zm(:,:)    ! midpoint geopotential height above the surface [m]
    real(r8),          intent(in)       :: zi(:,:)    ! interface geopotential height above sfc [m]

    integer, optional, intent(in)       :: primary                   ! primary detection algorithm
    integer, optional, intent(in)       :: backup                    ! backup detection algorithm
    integer,            intent(out)     :: tropLev(pcols)            ! tropopause level index   
    real(r8), optional, intent(out)     :: tropP(pcols)              ! tropopause pressure [Pa]  
    real(r8), optional, intent(out)     :: tropT(pcols)              ! tropopause temperature [K]
    real(r8), optional, intent(out)     :: tropZ(pcols)              ! tropopause height [m]
    
    ! Local Variable
    integer       :: primAlg            ! Primary algorithm  
    integer       :: backAlg            ! Backup algorithm  
  
    ! Initialize the results to a missing value, so that the algorithms will
    ! attempt to find the tropopause for all of them.
    tropLev(:) = NOTFOUND
    if (present(tropP)) tropP(:) = fillvalue
    if (present(tropT)) tropT(:) = fillvalue
    if (present(tropZ)) tropZ(:) = fillvalue
    
    ! Set the algorithms to be used, either the ones provided or the defaults.
    if (present(primary)) then
      primAlg = primary
    else
      primAlg = default_primary
    endif
    
    if (present(backup)) then
      backAlg = backup
    else
      backAlg = default_backup
    endif
    
    ! Try to find the tropopause using the primary algorithm.
    if (primAlg /= TROP_ALG_NONE) then
      call tropopause_findUsing(lchnk,ncol,pmid,pint,temp,zm,zi,primAlg, & ! in
           tropLev, tropP, tropT, tropZ)  ! inout
    endif
 
    if ((backAlg /= TROP_ALG_NONE) .and. any(tropLev(:) == NOTFOUND)) then
      call tropopause_findUsing(lchnk,ncol,pmid,pint,temp,zm,zi,backAlg, & ! in
           tropLev, tropP, tropT, tropZ)  ! inout
    endif
    
    return
  end subroutine tropopause_find
  
  ! Call the appropriate tropopause detection routine based upon the algorithm
  ! specifed.
  !
  ! NOTE: It is assumed that the output fields have been initialized by the
  ! caller, and only output values set to fillvalue will be detected.
  subroutine tropopause_findUsing(lchnk,ncol,pmid,pint,temp,zm,zi,algorithm,    &  ! in
             tropLev,tropP,tropT,tropZ)   ! inout

    implicit none

    integer,           intent(in)       :: lchnk
    integer,           intent(in)       :: ncol
    real(r8),          intent(in)       :: pmid(:,:)  ! midpoint pressure [Pa]
    real(r8),          intent(in)       :: pint(:,:)  ! top interface pressure [Pa]
    real(r8),          intent(in)       :: temp(:,:)  ! temperature [K]
    real(r8),          intent(in)       :: zm(:,:)    ! midpoint geopotential height above the surface [m]
    real(r8),          intent(in)       :: zi(:,:)    ! interface geopotential height above sfc [m]
    integer,            intent(in)      :: algorithm                 ! detection algorithm
    integer,            intent(inout)   :: tropLev(pcols)            ! tropopause level index   
    real(r8), optional, intent(inout)   :: tropP(pcols)              ! tropopause pressure [Pa]  
    real(r8), optional, intent(inout)   :: tropT(pcols)              ! tropopause temperature [K]
    real(r8), optional, intent(inout)   :: tropZ(pcols)              ! tropopause height [m]

    ! Dispatch the request to the appropriate routine.
    select case(algorithm)

      case(TROP_ALG_CLIMATE)
        call tropopause_climate(lchnk,ncol,pmid,pint,temp,zm,zi, & ! in
             tropLev,tropP,tropT,tropZ)  ! inout
      case(TROP_ALG_HYBSTOB)
        call tropopause_hybridstobie(lchnk,ncol,pmid,temp,zm, & ! in
             tropLev,tropP,tropT,tropZ)  ! inout
      case(TROP_ALG_TWMO)
        call tropopause_twmo(ncol,pmid,pint,temp,zm,zi,  & ! in
             tropLev,tropP,tropT,tropZ)  ! inout
      case(TROP_ALG_CPP)
        call tropopause_cpp(ncol,pmid,temp,zm,  & ! in
             tropLev,tropP,tropT,tropZ)  ! inout
      case default
        write(iulog, *) 'tropopause: Invalid detection algorithm (',  algorithm, ') specified.'
        call endrun
    end select
    
    return
  end subroutine tropopause_findUsing

  ! This routine interpolates the pressures in the physics state to
  ! find the pressure at the specified tropopause altitude.
  function tropopause_interpolateP(pmid,pmid_p1,pmid_m1,zm,zm_p1,zm_m1,tropLev,tropZ)
 
    implicit none

    real(r8),          intent(in)       :: pmid     ! midpoint pressure [Pa]
    real(r8),          intent(in)       :: pmid_p1  ! midpoint pressure at level + 1 [Pa]
    real(r8),          intent(in)       :: pmid_m1  ! midpoint pressure at level - 1 [Pa]
    real(r8),          intent(in)       :: zm       ! midpoint geopotential height above the surface [m]
    real(r8),          intent(in)       :: zm_p1    ! midpoint geopotential height above the surface at level + 1 [m]
    real(r8),          intent(in)       :: zm_m1    ! midpoint geopotential height above the surface at level - 1 [m]
    integer, intent(in)                 :: tropLev  ! tropopause level index   
    real(r8), optional, intent(in)      :: tropZ    ! tropopause height [m]
    real(r8)                            :: tropopause_interpolateP   ! function output tropopause pressure [Pa]
    
    ! Local Variables
    real(r8)   :: tropP              ! tropopause pressure [Pa]
    real(r8)   :: dlogPdZ            ! dlog(p)/dZ


    ! Interpolate the temperature linearly against log(P)
    
    ! Is the tropopause at the midpoint?
    if (tropZ == zm) then
      tropP = pmid
    
    elseif (tropZ > zm) then
    
      ! It is above the midpoint? Make sure we aren't at the top.
      if (tropLev > 1) then
        dlogPdZ = (log(pmid) - log(pmid_m1)) / &
          (zm - zm_m1) 
        tropP = pmid + exp((tropZ - zm) * dlogPdZ)
      endif
    else
      
      ! It is below the midpoint. Make sure we aren't at the bottom.
      if (tropLev < pver) then
        dlogPdZ =  (log(pmid_p1) - log(pmid)) / &
          (zm_p1 - zm)
        tropP = pmid + exp((tropZ - zm) * dlogPdZ)
      endif
    endif
    
    tropopause_interpolateP = tropP
  end function tropopause_interpolateP

  ! This routine interpolates the temperatures in the physics state to
  ! find the temperature at the specified tropopause pressure.
  function tropopause_interpolateT(pmid,pmid_p1,pmid_m1,temp,temp_p1,temp_m1,tropLev,tropP)
 
    implicit none

    real(r8),          intent(in)       :: pmid               ! midpoint pressure [Pa]
    real(r8),          intent(in)       :: pmid_p1            ! midpoint pressure at level + 1 [Pa]
    real(r8),          intent(in)       :: pmid_m1            ! midpoint pressure at level - 1 [Pa]
    real(r8),          intent(in)       :: temp               ! temperature [K]
    real(r8),          intent(in)       :: temp_p1            ! temperature at level + 1 [K]
    real(r8),          intent(in)       :: temp_m1            ! temperature at level - 1 [K]
    integer, intent(in)                 :: tropLev            ! tropopause level index   
    real(r8), optional, intent(in)      :: tropP              ! tropopause pressure [Pa]
    real(r8)                            :: tropopause_interpolateT  ! function output tropopause temperature [K]
    
    ! Local Variables
    real(r8)   :: tropT              ! tropopause temperature [K]
    real(r8)   :: dTdlogP            ! dT/dlog(P)
    
    ! Intrepolate the temperature linearly against log(P)

    ! Is the tropopause at the midpoint?
    if (tropP == pmid) then
      tropT = temp
    
    elseif (tropP < pmid) then
    
      ! It is above the midpoint? Make sure we aren't at the top.
      if (tropLev > 1) then
        dTdlogP = (temp - temp_m1) / & 
          (log(pmid) - log(pmid_m1))
        tropT = temp + (log(tropP) - log(pmid)) * dTdlogP
      endif
    else
      
      ! It is below the midpoint. Make sure we aren't at the bottom.
      if (tropLev < pver) then
        dTdlogP = (temp_p1 - temp) / &
          (log(pmid_p1) - log(pmid))
        tropT = temp + (log(tropP) - log(pmid)) * dTdlogP
      endif
    endif
    
    tropopause_interpolateT = tropT
  end function tropopause_interpolateT

  
  ! This routine interpolates the geopotential height in the physics state to
  ! find the geopotential height at the specified tropopause pressure.
  function tropopause_interpolateZ(pmid,pint,pint_p1,zm,zi,zi_p1,tropP)
 
    implicit none

    real(r8),          intent(in)       :: pmid     ! midpoint pressure [Pa]
    real(r8),          intent(in)       :: pint     ! top interface pressure [Pa]
    real(r8),          intent(in)       :: pint_p1  ! top interface pressure at level + 1 [Pa]
    real(r8),          intent(in)       :: zm       ! midpoint geopotential height above the surface [m]
    real(r8),          intent(in)       :: zi       ! interface geopotential height above sfc [m]
    real(r8),          intent(in)       :: zi_p1    ! interface geopotential height above sfc at level + 1 [m]
    real(r8), optional, intent(in)      :: tropP              ! tropopause pressure [Pa]
    real(r8)                            :: tropopause_interpolateZ  ! function output tropopause geopotential height [m]
    
    ! Local Variables
    real(r8)   :: tropZ              ! tropopause geopotential height [m]
    real(r8)   :: dZdlogP            ! dZ/dlog(P)
    
    ! Intrepolate the geopotential height linearly against log(P)

    ! Is the tropoause at the midpoint?
    if (tropP == pmid) then
      tropZ = zm
    
    elseif (tropP < pmid) then
    
      ! It is above the midpoint? Make sure we aren't at the top.
      dZdlogP = (zm - zi) / &
        (log(pmid) - log(pint))
      tropZ = zm + (log(tropP) - log(pmid)) * dZdlogP
    else
      
      ! It is below the midpoint. Make sure we aren't at the bottom.
      dZdlogP = (zm - zi_p1) / &
        (log(pmid) - log(pint_p1))
      tropZ = zm + (log(tropP) - log(pmid)) * dZdlogP
    endif
    
    tropopause_interpolateZ = tropZ
  end function tropopause_interpolateZ

  
  ! Output the tropopause pressure and temperature to the history files. Two sets
  ! of output will be generated, one for the default algorithm and another one
  ! using the default routine, but backed by a climatology when the default
  ! algorithm fails.
  subroutine tropopause_output(pstate)
    use cam_history,  only : outfld
    
    implicit none

    type(physics_state), intent(in)     :: pstate
  
    ! Local Variables
    integer       :: icol
    integer       :: alg
    integer       :: ncol                     ! number of columns in the chunk
    integer       :: lchnk                    ! chunk identifier
    integer       :: tropLev(pcols)           ! tropopause level index   
    real(r8)      :: tropP(pcols)             ! tropopause pressure [Pa]  
    real(r8)      :: tropT(pcols)             ! tropopause temperature [K] 
    real(r8)      :: tropZ(pcols)             ! tropopause height [m] 
    real(r8)      :: tropFound(pcols)         ! tropopause found  
    real(r8)      :: tropDZ(pcols, pver)      ! relative tropopause height [m] 
    real(r8)      :: tropPdf(pcols, pver)     ! tropopause probability distribution  

    ! Information about the chunk.  
    lchnk = pstate%lchnk
    ncol  = pstate%ncol

    ! Find the tropopause using the default algorithm backed by the climatology.
    call tropopause_find(lchnk,ncol,pstate%pint,pstate%pmid,pstate%t,pstate%zm,pstate%zi,  & ! in
         tropLev,tropP=tropP,tropT=tropT,tropZ=tropZ)   ! out / optional out
    
    tropPdf(:,:) = 0._r8
    tropFound(:) = 0._r8
    tropDZ(:,:) = fillvalue 
    do icol = 1, ncol
      if (tropLev(icol) /= NOTFOUND) then
        tropPdf(icol, tropLev(icol)) = 1._r8
        tropFound(icol) = 1._r8
        tropDZ(icol,:) = pstate%zm(icol,:) - tropZ(icol) 
      endif
    enddo

    call outfld('TROP_P',   tropP(:ncol),      ncol, lchnk)
    call outfld('TROP_T',   tropT(:ncol),      ncol, lchnk)
    call outfld('TROP_Z',   tropZ(:ncol),      ncol, lchnk)
    call outfld('TROP_DZ',  tropDZ(:ncol, :), ncol, lchnk)
    call outfld('TROP_PD',  tropPdf(:ncol, :), ncol, lchnk)
    call outfld('TROP_FD',  tropFound(:ncol),  ncol, lchnk)
    
    
    ! Find the tropopause using just the primary algorithm.
    call tropopause_find(lchnk,ncol,pstate%pint,pstate%pmid,pstate%t,pstate%zm,pstate%zi,  & ! in
         tropLev,tropP=tropP,tropT=tropT,tropZ=tropZ,   &  ! out / optional out
         backup=TROP_ALG_NONE)  ! optional in

    tropPdf(:,:) = 0._r8
    tropFound(:) = 0._r8
    tropDZ(:,:) = fillvalue 
    
    do icol = 1, ncol
      if (tropLev(icol) /= NOTFOUND) then
        tropPdf(icol, tropLev(icol)) = 1._r8
        tropFound(icol) = 1._r8
        tropDZ(icol,:) = pstate%zm(icol,:) - tropZ(icol) 
      endif
    enddo

    call outfld('TROPP_P',   tropP(:ncol),      ncol, lchnk)
    call outfld('TROPP_T',   tropT(:ncol),      ncol, lchnk)
    call outfld('TROPP_Z',   tropZ(:ncol),      ncol, lchnk)
    call outfld('TROPP_DZ',  tropDZ(:ncol, :), ncol, lchnk)
    call outfld('TROPP_PD',  tropPdf(:ncol, :), ncol, lchnk)
    call outfld('TROPP_FD',  tropFound(:ncol),  ncol, lchnk)
    
    ! Find the tropopause using just the cold point algorithm.
    call tropopause_find(lchnk,ncol,pstate%pint,pstate%pmid,pstate%t,pstate%zm,pstate%zi,  & ! in
         tropLev,tropP=tropP,tropT=tropT,tropZ=tropZ,  &  ! out / optional out
         primary=TROP_ALG_CPP,backup=TROP_ALG_NONE) ! optional in

    tropPdf(:,:) = 0._r8
    tropFound(:) = 0._r8
    tropDZ(:,:) = fillvalue

    do icol = 1, ncol
      if (tropLev(icol) /= NOTFOUND) then
        tropPdf(icol, tropLev(icol)) = 1._r8
        tropFound(icol) = 1._r8
        tropDZ(icol,:) = pstate%zm(icol,:) - tropZ(icol)
      endif
    enddo

    call outfld('TROPF_P',   tropP(:ncol),      ncol, lchnk)
    call outfld('TROPF_T',   tropT(:ncol),      ncol, lchnk)
    call outfld('TROPF_Z',   tropZ(:ncol),      ncol, lchnk)
    call outfld('TROPF_DZ',  tropDZ(:ncol, :), ncol, lchnk)
    call outfld('TROPF_PD',  tropPdf(:ncol, :), ncol, lchnk)
    call outfld('TROPF_FD',  tropFound(:ncol),  ncol, lchnk)
    
    
    return
  end subroutine tropopause_output
end module tropopause
