!===============================================================================
! Seasalt for Modal Aerosol Model
!===============================================================================
module seasalt_model
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use ppgrid,         only: pcols
  use cam_abortutils, only: endrun
  use modal_aero_data,only: ntot_amode
  use tracer_data,    only: trfld, trfile
  use cam_logfile,    only: iulog

  implicit none
  private

  public :: seasalt_nbin
  public :: seasalt_nnum
  public :: seasalt_names
  public :: seasalt_indices
  public :: seasalt_init
  public :: seasalt_emis
  public :: seasalt_active
  public :: marine_organic_emis

  public :: n_ocean_data
  public :: nslt_om
  public :: F_eff_out                 ! Output effective enrichment ratio?
                                      !  (logical, currently set to FALSE) 
  public :: has_mam_mom               ! run with marine organics?
                                      !  (logical, set to TRUE if user supplies file)
  public :: advance_ocean_data        ! advance ocean data in time
  public :: init_ocean_data           ! initialize ocean data variables
  public :: ocean_data_readnl         ! read ocean data namelist

  integer, parameter :: nslt = max(3,ntot_amode-3)
  integer, parameter :: nnum = nslt

  integer, parameter :: nslt_om = 3
  integer, parameter :: nnum_om = 1
  integer, parameter :: om_num_modes = 3
  character(len=6),parameter :: seasalt_names(nslt+nslt_om+nnum+nnum_om) = &
       (/ 'ncl_a1', 'ncl_a2', 'ncl_a3', &
       'mom_a1', 'mom_a2', 'mom_a4', &
       'num_a1', 'num_a2', 'num_a3', 'num_a4'/)
  integer, dimension(om_num_modes), parameter :: om_num_ind =  (/ 1, 2, 4 /)

  integer, parameter :: seasalt_nbin = nslt+nslt_om
  integer, parameter :: seasalt_nnum = nnum+nnum_om

  integer, parameter :: & ! number of ocean data fields
       n_ocean_data = 4

!  logical, parameter   :: F_eff_out = .true.
  logical, parameter   :: F_eff_out = .false.
  type(trfld), pointer :: fields(:)
  type(trfile)         :: file

  ! Settings for marine organics code

  real(r8), parameter :: small_oceanorg = 1.0e-30 ! smallest ocean organic concentration allowed

  integer :: seasalt_indices(seasalt_nbin+seasalt_nnum)

  logical :: seasalt_active = .false.

! Parameters for organic sea salt emissions
    real(r8), parameter :: mw_carbon = 12.0107_r8              ! molecular weight for carbon element [g/mol]
    real(r8), parameter :: dens_vol_NaCl_in_seawat = 35875._r8 ! approx volume density of salt in seawater [g/m^3]
    real(r8) :: dens_srf_NaCl_in_bubsrf                        ! salt per area bubble surface [g/m^2]
    integer, parameter  :: n_org_max = 3                       ! max number of organic compound classes
    integer, parameter  :: n_org  = 3                          ! actual number of organic compound classes (scheme dependent)

! Marine organics namelist variables

! Namelist variables related to dataset specification
   character(len=32)   :: specifier(n_ocean_data) = ''
   character(len=256)  :: filename = ' '
   character(len=256)  :: filelist = ' '
   character(len=256)  :: datapath = ' '
   character(len=32)   :: datatype = 'CYCLICAL'
   integer             :: data_cycle_yr = 0
   logical             :: rmv_file = .false.
   integer             :: fixed_ymd = 0
   integer             :: fixed_tod = 0

! Namelist variables for parameterization specification
  real(r8) :: l_bub = 0.1e-6_r8 ! Bubble film thickness [m]

  ! Determine mixing state for MOM emissions.
  ! Currently implemented options:
  ! mixing_state = 0 : total external mixture, replace mass
  !                1 : total external mixture, add to mass
  !                2 : total internal mixture, replace mass
  !                3 : total internal mixture, add to mass
   integer             :: mixing_state = 1

  ! Selection of alternate parameterizations
  ! Set fmoa=1 for Burrows et al., 2014 parameterization
  !     fmoa=2 for Gantt et al., 2011 parameterization
  !     fmoa=3 for simple parameterization based on Quinn et al., 2014
  !     fmoa=4 for Rinaldi et al. (JGR, 2013)
   integer             :: fmoa = 1

! TODO SMB: Implement better mechanism for setting this switch.
   logical :: has_mam_mom = .true.

! Order: mpoly, mprot, mlip
    real(r8), dimension(n_org), parameter :: & ! OM:OC mass ratios for input fields (mpoly, mprot, mlip)
         OM_to_OC_in = (/ 2.3_r8, 2.2_r8, 1.3_r8 /)
    real(r8), dimension(n_org), parameter :: & ! Langmuir parameters (inverse C_1/2)  [m3 mol-1]
         alpha_org = (/ 90.58_r8, 25175._r8, 18205._r8 /)
! Molecular weights needed for output of optional diagnostic variable F_eff
    real(r8), dimension(n_org), parameter :: & ! Molecular weights [g mol-1]
         mw_org   = (/ 250000._r8, 66463._r8, 284._r8 /)
    real(r8), dimension(n_org), parameter :: & ! mass per sq. m at saturation
         dens_srf_org = (/ 0.1376_r8, 0.00219_r8, 0.002593_r8 /) ! Mw_org / a_org

    real(r8), parameter :: sst_sz_range_lo (nslt+nslt_om) = &
         (/ 0.08e-6_r8,  0.02e-6_r8,  1.0e-6_r8, &  ! accu, aitken, coarse
            0.08e-6_r8,  0.02e-6_r8,  0.08e-6_r8 /) ! accu, aitken, POM accu
    real(r8), parameter :: sst_sz_range_hi (nslt+nslt_om) = &
         (/ 1.0e-6_r8,   0.08e-6_r8, 10.0e-6_r8, &  ! accu, aitken, coarse
            1.0e-6_r8,   0.08e-6_r8,  1.0e-6_r8 /)  ! accu, aitken, POM accu

contains
  
  !=============================================================================
  !=============================================================================
  subroutine seasalt_init
    use sslt_sections, only: sslt_sections_init
    use constituents,  only: cnst_get_ind

    integer :: m

    do m = 1, seasalt_nbin
       call cnst_get_ind(seasalt_names(m), seasalt_indices(m),abrtf=.false.)
    enddo
    do m = 1, seasalt_nnum
       call cnst_get_ind(seasalt_names(seasalt_nbin+m), seasalt_indices(seasalt_nbin+m),abrtf=.false.)
    enddo

    seasalt_active = any(seasalt_indices(:) > 0)

    if (.not.seasalt_active) return

    call sslt_sections_init()

  end subroutine seasalt_init

  !=============================================================================

!-------------------------------------------------------------------
! Advance ocean data fields to the current time step
!
! Adapted from prescribed_aero_adv
!
! Author: Susannah M. Burrows
! Date: 13 Jan 2015
!-------------------------------------------------------------------
subroutine advance_ocean_data(state, pbuf2d)
    use physics_types,  only : physics_state
    use tracer_data,    only : advance_trcdata, get_fld_data, put_fld_data
    use ppgrid,         only : begchunk, endchunk
    use physics_buffer, only : physics_buffer_desc, pbuf_get_chunk
    use cam_history,    only : outfld

    implicit none

    type(physics_state), intent(in)    :: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    type(physics_buffer_desc), pointer :: pbuf_chnk(:)

    integer :: i,c,ncol
!    real(r8),pointer :: outdata(:,:)
!    real(r8) :: outdata(pcols,begchunk:endchunk)
    real(r8) :: outdata(pcols,1)
    integer lchnk

!    write(iulog,*) 'Advancing ocean data ...' ! for debugging
    call advance_trcdata( fields, file, state, pbuf2d )
!    write(iulog,*) 'Done advancing ocean data ...' ! for debugging

! Add new values to history files
    fldloop:do i = 1,n_ocean_data

       chnkloop: do c = begchunk,endchunk
          ncol = state(c)%ncol
          pbuf_chnk => pbuf_get_chunk(pbuf2d, c)
          lchnk = state(c)%lchnk

          call get_fld_data( fields, fields(i)%fldnam, outdata(:ncol,:), ncol, lchnk, pbuf_chnk)

          ! work-around for interpolation errors that introduce negative values
          ! near coasts: reset negative values to zero.
          where (outdata(:ncol,1) < small_oceanorg)
             outdata(:ncol,1) = 0.0_r8
          end where

          call put_fld_data( fields, fields(i)%fldnam, outdata(:ncol,:), ncol, lchnk, pbuf_chnk)

          ! The following line is probably redundant but is included for safety
          call get_fld_data( fields, fields(i)%fldnam, outdata(:ncol,:), ncol, lchnk, pbuf_chnk)

          call outfld( trim(fields(i)%fldnam), outdata(:ncol,1), ncol, lchnk )
       enddo chnkloop

    enddo fldloop

end subroutine advance_ocean_data

subroutine ocean_data_readnl(nlfile)
   use spmd_utils,      only: masterproc
   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'ocean_data_readnl'

   character(len=32)   :: mam_mom_specifier(n_ocean_data)
   character(len=256)  :: mam_mom_filename
   character(len=256)  :: mam_mom_filelist
   character(len=256)  :: mam_mom_datapath
   character(len=32)   :: mam_mom_datatype
   integer             :: mam_mom_cycle_yr
   logical             :: mam_mom_rmv_file
   integer             :: mam_mom_fixed_ymd
   integer             :: mam_mom_fixed_tod

   real(r8)            :: mam_mom_bubble_thickness
   integer             :: mam_mom_mixing_state
   integer             :: mam_mom_parameterization

   namelist /mam_mom_nl/ &
      mam_mom_specifier, &
      mam_mom_filename,  &
      mam_mom_filelist,  &
      mam_mom_datapath,  &
      mam_mom_datatype,  &
      mam_mom_rmv_file,  &
      mam_mom_cycle_yr,  &
      mam_mom_fixed_ymd, &
      mam_mom_fixed_tod, &
      mam_mom_bubble_thickness, &
      mam_mom_mixing_state, &
      mam_mom_parameterization

   !-----------------------------------------------------------------------------

   ! Initialize namelist variables from local module variables.
   mam_mom_specifier= specifier
   mam_mom_filename = filename
   mam_mom_filelist = filelist
   mam_mom_datapath = datapath
   mam_mom_datatype = datatype
   mam_mom_rmv_file = rmv_file
   mam_mom_cycle_yr = data_cycle_yr
   mam_mom_fixed_ymd= fixed_ymd
   mam_mom_fixed_tod= fixed_tod

   mam_mom_bubble_thickness = l_bub
   mam_mom_mixing_state     = mixing_state
   mam_mom_parameterization = fmoa

   ! Read aerosol namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'mam_mom_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, mam_mom_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      call freeunit(unitn)
      close(unitn)
   endif

!      mam_mom_specifier, & ! Names of variables containing aerosol data in the prescribed aerosol datasets.
!      mam_mom_filename,  & ! Filename of dataset for prescribed marine organic matter emissions.
!      mam_mom_filelist,  & ! Filename of file that contains a sequence of filenames for prescribed
!                         & ! aerosols.  The filenames in this file are relative to the directory specied
!                         & ! by mam_mom_datapath.
!      mam_mom_datapath,  & ! Full pathname of the directory that contains the files specified in mam_mom_filelist.
!      mam_mom_datatype,      & ! Type of time interpolation for data in mam_mom files.
!                         & ! Can be set to 'CYCLICAL', 'SERIAL', 'INTERP_MISSING_MONTHS', or 'FIXED'.
!      mam_mom_rmv_file,  & ! Remove the file containing prescribed aerosol deposition fluxes from local disk when no longer needed.
!      mam_mom_cycle_yr,  & ! The  cycle year of the prescribed aerosol flux data
!                         & ! if mam_mom_datatype  is 'CYCLICAL'.
!      mam_mom_fixed_ymd, & ! The date at which the prescribed aerosol flux data is fixed
!                         & ! if mam_mom_datatype is 'FIXED'.
!      mam_mom_fixed_tod, & ! The time of day (seconds) corresponding to mam_mom_fixed_ymd
!                           ! at which the prescribed aerosol flux data is fixed
!                           ! if mam_mom_datatype is 'FIXED'.
! mam_mom_bubble_thickness, & ! Bubble film thickness (in m) for marine organic aerosol emission
!                             ! mechanism.  The physically reasonable range is approximately
!                             ! (0.1 - 1) x 10 ^-6.
! mam_mom_mixing_state, &   ! Switch to select mixing state assumption in marine organic aerosol
!                           ! code. Currently implemented options: 0 : total external mixture, add
!                           ! to mass; 1 : total external mixture, replace mass; 2 : total
!                           ! internal mixture, add to mass; 3 : total internal mixture, replace
!                           ! mass.
! mam_mom_parameterization  ! Selection of alternate parameterizations for marine organic matter
!                           ! emissions.  Set fmoa=1 for Burrows et al., 2014 parameterization;
!                           ! fmoa=2 for Gantt et al. (2011, ACP) parameterization; fmoa=3 for
!                           ! simple parameterization based on Quinn et al., 2014; fmoa=4 for
!                           ! Rinaldi et al. (JGR, 2013).


#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(mam_mom_specifier,len(mam_mom_specifier)*n_ocean_data,   mpichar, 0, mpicom)
   call mpibcast(mam_mom_filename, len(mam_mom_filename),   mpichar, 0, mpicom)
   call mpibcast(mam_mom_filelist, len(mam_mom_filelist),   mpichar, 0, mpicom)
   call mpibcast(mam_mom_datapath, len(mam_mom_datapath),   mpichar, 0, mpicom)
   call mpibcast(mam_mom_datatype, len(mam_mom_datatype),   mpichar, 0, mpicom)
   call mpibcast(mam_mom_rmv_file, 1, mpilog, 0, mpicom)
   call mpibcast(mam_mom_cycle_yr, 1, mpiint, 0, mpicom)
   call mpibcast(mam_mom_fixed_ymd,1, mpiint, 0, mpicom)
   call mpibcast(mam_mom_fixed_tod,1, mpiint, 0, mpicom)

   call mpibcast(mam_mom_bubble_thickness,1, mpir8, 0, mpicom)
   call mpibcast(mam_mom_mixing_state,1, mpiint, 0, mpicom)
   call mpibcast(mam_mom_parameterization,1, mpiint, 0, mpicom)
#endif

   ! Update module variables with user settings.
   specifier     = mam_mom_specifier
   filename      = mam_mom_filename
   filelist      = mam_mom_filelist
   datapath      = mam_mom_datapath
   datatype      = mam_mom_datatype
   rmv_file      = mam_mom_rmv_file
   data_cycle_yr = mam_mom_cycle_yr
   fixed_ymd     = mam_mom_fixed_ymd
   fixed_tod     = mam_mom_fixed_tod

   l_bub         = mam_mom_bubble_thickness
   mixing_state  = mam_mom_mixing_state
   fmoa          = mam_mom_parameterization



end subroutine ocean_data_readnl

  !=============================================================================
subroutine seasalt_emis(lchnk, ncol, u10cubed, srf_temp, ocnfrc, emis_scale, & ! in
                        cflx)                                                  ! inout

    use sslt_sections, only: nsections, fluxes

    ! input
    integer, intent(in)  :: lchnk
    integer, intent(in)  :: ncol
    real(r8), intent(in) :: u10cubed(pcols)   ! 3.41 power of 10m wind 
    real(r8), intent(in) :: srf_temp(pcols)   ! sea surface temperature [K]
    real(r8), intent(in) :: ocnfrc(pcols)     ! ocean fraction [unitless]
    real(r8), intent(in) :: emis_scale        ! sea salt emission tuning factor [unitless]

    ! output
    real(r8), intent(inout) :: cflx(:,:)      ! mass and number emission fluxes for aerosols [kg/m2/s or #/m2/s]

    ! local vars
    integer  :: ispec, ibin
    integer  :: mass_mode_idx, num_mode_idx
    real(r8) :: fi(pcols,nsections)           ! sea salt number fluxes in each size bin [#/m2/s]
    real(r8) :: cflx_tmp1(pcols)              ! temp array for calculating emission fluxes [kg/m2/s or #/m2/s]
   
    integer, parameter :: num_flx_flag  = 0
    integer, parameter :: mass_flx_flag = 1

    fi(:ncol,:nsections) = fluxes( srf_temp, u10cubed, ncol )

    ! calculate seasalt number emission fluxes
    call seasalt_emisflx_calc(ncol, fi, ocnfrc, emis_scale, num_flx_flag, &  ! in
                              cflx)                                          ! inout
 
    ! calculate seasalt mass emission fluxes
    call seasalt_emisflx_calc(ncol, fi, ocnfrc, emis_scale, mass_flx_flag, & ! in
                              cflx)                                          ! inout

end subroutine seasalt_emis


subroutine seasalt_emisflx_calc(ncol, fi, ocnfrc, emis_scale, flx_type, & ! in
                                cflx)                                     ! inout
   
    use sslt_sections, only: nsections, fluxes, Dg, rdry
    use mo_constants,  only: dns_aer_sst=>seasalt_density, pi

    ! input
    integer, intent(in)  :: ncol
    real(r8), intent(in) :: fi(pcols, nsections)        ! sea salt number fluxes in each size bin [#/m2/s]
    real(r8), intent(in) :: ocnfrc(pcols)               ! ocean fraction [unitless]
    real(r8), intent(in) :: emis_scale                  ! sea salt emission tuning factor [unitless]
    integer, intent(in)  :: flx_type                    ! number ==0 or mass ==1 flux to be calculated

    ! output
    real(r8), intent(inout) :: cflx(:,:)      ! mass and number emission fluxes for aerosols [kg/m2/s or #/m2/s]

    ! local vars
    integer  :: ispec, ibin
    integer  :: mode_idx
    real(r8) :: cflx_tmp1(pcols)              ! temp array for calculating emission fluxes [kg/m2/s or #/m2/s]
    integer  :: num_idx_append

    num_idx_append = 0

    if (flx_type == 0) num_idx_append = nslt+nslt_om
   
    do ispec = 1, nslt
       mode_idx  = seasalt_indices(num_idx_append+ispec)  ! Index of number mode
       
       if (mode_idx>0) then 
          ! Note for C++ port, ideally initialization of clfx should be done
          ! for both number and mass modes. We only zero-out number, because
          ! zero-out mass fluxes cause NBFB                 
          if (flx_type == 1) cflx(:ncol,mode_idx) = 0.0_r8 ! for number fluxes, zero-out the fluxes 

          ! Total number flux per mode
          do ibin=1, nsections
             cflx_tmp1(:ncol) = 0.0_r8

             if (Dg(ibin) >= sst_sz_range_lo(ispec) .and. Dg(ibin) < sst_sz_range_hi(ispec)) then
                cflx_tmp1(:ncol) = fi(:ncol,ibin)*ocnfrc(:ncol)*emis_scale 
              
                ! For mass fluxes, multiply by the diameter
                ! Note for C++ port, 4._r8/3._r8*pi*rdry(ibin)**3*dns_aer_sst
                ! can be factored out in a function, but it is not done here as
                ! the results might not stay BFB
                if (flx_type == 1) cflx_tmp1(:ncol) = cflx_tmp1(:ncol) * 4._r8/3._r8*pi*rdry(ibin)**3*dns_aer_sst

                ! Mixing state 3: internal mixture, add OM to mass and number
                cflx(:ncol,mode_idx) = cflx(:ncol,mode_idx) + cflx_tmp1(:ncol)
             endif
          enddo          
       endif       
    enddo

end subroutine seasalt_emisflx_calc


subroutine marine_organic_emis(lchnk, ncol, u10cubed, srf_temp, ocnfrc, emis_scale, &  ! in
                               cflx)                                                   ! inout

    use sslt_sections, only: nsections, fluxes
    use spmd_utils,    only: masterproc

    ! input
    integer, intent(in)  :: lchnk
    integer, intent(in)  :: ncol
    real(r8), intent(in) :: u10cubed(pcols)    ! 3.41 power of 10m wind
    real(r8), intent(in) :: srf_temp(pcols)    ! sea surface temperature [K]
    real(r8), intent(in) :: ocnfrc(pcols)      ! ocean fraction [unitless] 
    real(r8), intent(in) :: emis_scale         ! sea salt emission tuning factor [unitless]
    
    ! output
    real(r8), intent(inout) :: cflx(:,:)       ! mass and number emission fluxes for aerosols [kg/m2/s or #/m2/s]

    ! local vars
    integer  :: ifld
    real(r8) :: fi(pcols,nsections)            ! sea salt number fluxes in each size bin [#/m2/s]

   real(r8), pointer :: chla(:)          ! for Gantt et al. (2011) organic mass fraction
   real(r8), pointer :: mpoly(:)         ! for Burrows et al. (2014) organic mass fraction
   real(r8), pointer :: mprot(:)         ! for Burrows et al. (2014) organic mass fraction
   real(r8), pointer :: mlip(:)          ! for Burrows et al. (2014) organic mass fraction

   logical, parameter :: emit_this_mode(om_num_modes) = (/ .true., .true., .false. /)

   real(r8) :: mass_frac_bub_section(pcols, n_org_max, nsections)
   real(r8) :: om_ssa(pcols, nsections)



   fi(:ncol,:nsections) = fluxes( srf_temp, u10cubed, ncol )

   nullify(chla)
   nullify(mpoly)
   nullify(mprot)
   nullify(mlip)

   fldloop: do ifld = 1, n_ocean_data
      select case (trim(fields(ifld)%fldnam))
         case ("chla")
             chla   => fields(ifld)%data(:ncol,1,lchnk)
         case ("mpoly")
             mpoly  => fields(ifld)%data(:ncol,1,lchnk)
         case ("mprot")
             mprot  => fields(ifld)%data(:ncol,1,lchnk)
         case ("mlip")
             mlip   => fields(ifld)%data(:ncol,1,lchnk)
         case default
             if ( masterproc ) then
                write(iulog,*) 'Unknown field name '//fields%fldnam//' in ocean_data fields ...'
             endif
      end select
   enddo fldloop

   mass_frac_bub_section(:ncol,:,:) = 0.0_r8
   om_ssa(:ncol,:) = 0.0_r8

   ! Calculate marine organic aerosol mass fraction based on Burrows et al., ACP (2013)
   call calc_om_ssa(ncol, mpoly(:ncol), mprot(:ncol), mlip(:ncol), &      ! in
                    mass_frac_bub_section(:ncol, :, :), om_ssa(:ncol, :)) ! inout


! Calculate emission of MOM mass.

! Determine which modes to emit MOM in depending on mixing state assumption

    ! OM modes (m in this loop)
    ! m=1 : accu       (internal w/ SS)
    ! m=2 : Aitken     (internal w/ SS)
    ! m=3 : accu MOM   (external)
    ! m=4 : Aitken MOM (external)

    ! Total external mixture: emit only in modes 4, 5
   
    call marine_organic_numflx_calc(ncol, fi, ocnfrc, emis_scale, &  ! in 
                                    om_ssa, emit_this_mode, &        ! in
                                    cflx)                            ! inout

    call marine_organic_massflx_calc(ncol, fi, ocnfrc, emis_scale, om_ssa, &  ! in
                                     mass_frac_bub_section, emit_this_mode, & ! in
                                     cflx)                                    ! inout

end subroutine marine_organic_emis


subroutine marine_organic_numflx_calc(ncol, fi, ocnfrc, emis_scale, &  ! in
                                      om_ssa, emit_this_mode, &        ! in
                                      cflx)                            ! inout

    use sslt_sections, only: nsections, Dg

    ! input
    integer, intent(in)  :: ncol
    real(r8), intent(in) :: fi(pcols, nsections)          ! sea salt number fluxes in each size bin [#/m2/s]    
    real(r8), intent(in) :: ocnfrc(pcols)                 ! ocean fraction [unitless] 
    real(r8), intent(in) :: emis_scale                    ! sea salt emission tuning factor [unitless]
    real(r8), intent(in) :: om_ssa(pcols, nsections)      ! marine organic aerosol fraction per size bin [unitless]
    logical, intent(in)  :: emit_this_mode(om_num_modes)  ! logical flags turn on/off marine organic emission in aerosol modes 

    ! output
    real(r8), intent(inout) :: cflx(:,:)       ! mass and number emission fluxes for aerosols [kg/m2/s or #/m2/s]

    ! local vars
    integer  :: ispec, ibin
    integer  :: num_mode_idx
    integer  :: om_num_idx
    real(r8) :: cflx_tmp1(pcols)               ! temp array for calculating emission fluxes [kg/m2/s or #/m2/s]

    if (size(om_num_ind) .eq. 1) then
       call endrun( "Error: om_num_ind is a scalar, but attempting to calculate MOM. &
                     Something bad happened!!  We should never get here!")
    endif

    ! Loop over OM modes
    do ispec = 1, om_num_modes ! modes in which to emit OM

       om_num_idx = om_num_ind(ispec)
       num_mode_idx=seasalt_indices(nslt+nslt_om+om_num_idx)

      ! add number tracers for organics-only modes
      if (emit_this_mode(ispec)) then
         do ibin=1, nsections
            cflx_tmp1(:ncol) = 0.0_r8
            if (Dg(ibin) >= sst_sz_range_lo(nslt+ispec) .and. Dg(ibin) < sst_sz_range_hi(nslt+ispec)) then
               cflx_tmp1(:ncol)=fi(:ncol,ibin)*ocnfrc(:ncol)*emis_scale

               ! Mixing state 3: internal mixture, add OM to mass and number
               cflx(:ncol,num_mode_idx) = cflx(:ncol,num_mode_idx) + cflx_tmp1(:ncol) * &
                                          (1._r8 / (1._r8 - om_ssa(:ncol, ibin)) - 1._r8)
             endif
          enddo
       endif
    enddo

end subroutine marine_organic_numflx_calc


subroutine marine_organic_massflx_calc(ncol, fi, ocnfrc, emis_scale, om_ssa, &  ! in
                                       mass_frac_bub_section, emit_this_mode, & ! in
                                       cflx)                                    ! inout

    use sslt_sections, only: nsections, Dg, rdry
    use mo_constants,  only: dns_aer_sst=>seasalt_density, pi

    ! input
    integer, intent(in)  :: ncol
    real(r8), intent(in) :: fi(pcols, nsections)             ! sea salt number fluxes in each size bin [#/m2/s]
    real(r8), intent(in) :: ocnfrc(pcols)                    ! ocean fraction [unitless] 
    real(r8), intent(in) :: emis_scale                       ! sea salt emission tuning factor [unitless]
    real(r8), intent(in) :: om_ssa(pcols, nsections)         ! marine organic aerosol fraction per size bin [unitless]
    real(r8), intent(in) :: mass_frac_bub_section(pcols, n_org_max, nsections) ! marine organic aerosol fraction per organic species per size bin [unitless]
    logical, intent(in)  :: emit_this_mode(om_num_modes)     ! logical flags turn on/off marine organic emission in aerosol modes

    ! output
    real(r8), intent(inout) :: cflx(:,:)       ! mass and number emission fluxes for aerosols [kg/m2/s or #/m2/s]

    ! local vars
    integer  :: ispec, ibin, iorg
    integer  :: mass_mode_idx
    real(r8) :: cflx_tmp1(pcols)               ! temp array for calculating emission fluxes [kg/m2/s or #/m2/s]


    do ispec=1,nslt_om
       mass_mode_idx = seasalt_indices(nslt+ispec)

       cflx(:ncol,mass_mode_idx)=0.0_r8

       if (emit_this_mode(ispec)) then
          do iorg = 1, n_org
             do ibin = 1, nsections
                if (Dg(ibin)>=sst_sz_range_lo(nslt+ispec) .and. Dg(ibin)<sst_sz_range_hi(nslt+ispec)) then
                   cflx_tmp1(:ncol) = fi(:ncol,ibin)*ocnfrc(:ncol)*emis_scale &
                                     *4._r8/3._r8*pi*rdry(ibin)**3*dns_aer_sst ! should use dry size, convert from number to mass flux (kg/m2/s)

                   ! Mixing state 3: internal mixture, add OM to mass and number
                   where (om_ssa(:ncol,ibin) > 0.0_r8) ! avoid division by zero
                      cflx(:ncol,mass_mode_idx) = cflx(:ncol,mass_mode_idx) + cflx_tmp1(:ncol) &
                                                  * mass_frac_bub_section(:ncol, iorg, ibin) / om_ssa(:ncol, ibin) * &
                                                 (1._r8 / (1._r8 - om_ssa(:ncol, ibin)) - 1._r8)
                   elsewhere
                      cflx(:ncol,mass_mode_idx) = cflx(:ncol,mass_mode_idx)
                   endwhere
                endif
             enddo
          enddo
       endif
    enddo

end subroutine marine_organic_massflx_calc


subroutine calc_om_ssa(ncol, mpoly_in, mprot_in, mlip_in, & ! in
                       mass_frac_bub_section, om_ssa)       ! inout

   !----------------------------------------------------------------------- 
   ! Purpose:
   ! Calculate OM fraction for five organic classes, and overall
   ! effective organic enrichment, following Burrows et al., ACP (2013).
   !
   ! Author:
   ! Susannah Burrows, 9 Mar 2015
   !----------------------------------------------------------------------- 
   use sslt_sections, only: nsections
   use cam_history,   only: outfld
   implicit none
   !-----------------------------------------------------------------------
   ! Input variables:
   integer,  intent(in) :: ncol
   real(r8), intent(in) :: mpoly_in(ncol)         ! for Burrows et al. (2014) organic mass fraction
   real(r8), intent(in) :: mprot_in(ncol)         ! for Burrows et al. (2014) organic mass fraction
   real(r8), intent(in) :: mlip_in(ncol)          ! for Burrows et al. (2014) organic mass fraction
   !
   ! Output variables
   real(r8), intent(inout) :: mass_frac_bub_section(:,:,:)          ! mass fraction per organic class per size bin [unitless]
   real(r8), intent(inout) :: om_ssa(:,:)                           ! mass fraction per size bin [unitless]

   !
   ! Local variables
   real(r8) :: alpha_help(ncol)
   real(r8) :: mass_frac_bub_tot(ncol)
   real(r8) :: om_conc(ncol, n_org)             ! mole concentration of organic [mol/m^3]
   real(r8) :: theta(ncol, n_org)               ! fractional surface coverage [unitless]
   real(r8) :: mass_frac_bub(ncol, n_org)
   real(r8) :: theta_help(ncol, n_org)
   real(r8) :: mass_frac_bub_help(ncol, n_org)

   real(r8), parameter :: omfrac_max = 0.78  ! OMF maximum and minimum values -- max from Rinaldi et al. (2013)
   real(r8), parameter :: liter_to_m3 = 1.0e-3_r8  

   integer  :: iorg
   !
   !-----------------------------------------------------------------------

   ! Initialize arrays to zero for safety.
   theta(:,:)              = 0.0_r8
   theta_help(:,:)         = 0.0_r8
   mass_frac_bub_tot(:)    = 0.0_r8
   mass_frac_bub(:,:)      = 0.0_r8
   mass_frac_bub_help(:,:) = 0.0_r8
   om_conc(:,:)            = 0.0_r8
   

   ! Convert input fields from [(mol C) L-1] to [(g OM) m-3] and store in single array
   om_conc(:, 1) = mpoly_in(:)  * liter_to_m3 * OM_to_OC_in(1) * mw_carbon
   om_conc(:, 2) = mprot_in(:)  * liter_to_m3 * OM_to_OC_in(2) * mw_carbon
   om_conc(:, 3) = mlip_in(:)   * liter_to_m3 * OM_to_OC_in(3) * mw_carbon


   ! Calculate the surface coverage by class
   do iorg = 1, n_org
      ! Bulk mass concentration [mol m-3] = [g m-3] / [g mol-1]
      om_conc(:, iorg)     = om_conc(:, iorg) / mw_org(iorg)
      ! use theta_help as work array -- theta_help = alpha(i) * x(i)
      theta_help(:, iorg)  = alpha_org(iorg) * om_conc(:, iorg)
   enddo
   alpha_help(:) = sum(theta, dim=2)
   
   do iorg = 1, n_org
      ! complete calculation -- theta = alpha(i) * x(i) / (1 + sum( alpha(i) * x(i) ))
      theta(:, iorg)         = theta_help(:, iorg) / (1.0_r8 + alpha_help(:))
      ! Calculate the organic mass per area (by class) [g m-2]
      !  (use mass_frac_bub_help as local work array -- organic mass per area in g per m2)
      mass_frac_bub_help(:, iorg) = theta(:, iorg) * dens_srf_org(iorg)
   enddo

   ! Calculate g NaCl per m2
   dens_srf_NaCl_in_bubsrf = dens_vol_NaCl_in_seawat*l_bub ! Redundant, but allows for easier adjustment to l_bub

   ! mass_frac_bub = 2*[g OM m-2] / (2*[g OM m-2] * [g NaCl m-2])
   ! Factor 2 for bubble bilayer (coated on both surfaces of film)
   do iorg = 1, n_org
      mass_frac_bub(:, iorg) = 2.0_r8*mass_frac_bub_help(:, iorg) / &
           (2.0_r8*sum(mass_frac_bub_help(:, :), dim=2) + dens_srf_NaCl_in_bubsrf)
   enddo

   mass_frac_bub_tot(:) = sum(mass_frac_bub, dim=2)

   do iorg = 1, n_org
      where (mass_frac_bub_tot(:) > omfrac_max)
         mass_frac_bub(:, iorg) = mass_frac_bub(:, iorg) / mass_frac_bub_tot(:) * omfrac_max
      endwhere
   enddo

   ! Must exceed threshold value small_oceanorg
   where (mass_frac_bub(:, :) < small_oceanorg)
      mass_frac_bub(:, :) = 0.0_r8
   endwhere

   mass_frac_bub_tot(:) = sum(mass_frac_bub, dim=2)


   ! Distribute mass fraction evenly into Aitken and accumulation modes
   call omfrac_accu_aitk(mass_frac_bub_tot(:), & ! in
                         om_ssa(:,:))            ! inout

   ! mass_frac_bub_section(pcols, n_org_max, nsections) -- org classes in dim 2, size nsections in dim 3
   mass_frac_bub_section(:, :, :)   = 0.0_r8

   do iorg = 1, n_org
      call omfrac_accu_aitk(mass_frac_bub(:, iorg), & ! in
                            mass_frac_bub_section(:, iorg, :)) ! inout
   enddo

end subroutine calc_om_ssa


subroutine omfrac_accu_aitk(om_ssa_in, & ! in
                            om_ssa)      ! inout

   !----------------------------------------------------------------------- 
   ! Purpose:
   ! Put OM fraction directly into aitken and accumulation modes and don't
   ! exceed om_ssa_max.
   use sslt_sections, only: nsections, Dg
   implicit none
   !-----------------------------------------------------------------------
   ! Arguments:
   !
   real(r8), intent(in)    :: om_ssa_in(:)     ! mass fraction [unitless]
   real(r8), intent(inout) :: om_ssa(:,:)      ! mass fraction for each size bin [unitless]
   !
   ! Local variables
   !
   integer  :: ibin
   real(r8), parameter :: om_ssa_max = 1.0_r8
   !
   !-----------------------------------------------------------------------

   ! Initialize array and set to zero for "fine sea salt" and "coarse sea salt" modes
   om_ssa(:, :) = 0.0_r8

   ! distribute OM fraction!
   do ibin = 1, nsections
      ! update only in Aitken and accu. modes
      if ((Dg(ibin) >= sst_sz_range_lo(2)) .and. (Dg(ibin) < sst_sz_range_hi(1))) then
             om_ssa(:, ibin) = om_ssa_in(:)
      endif
   enddo

   ! For safety, force fraction to be within bounds [0, 1]
   where (om_ssa(:, :) >  om_ssa_max)
      om_ssa(:, :) = om_ssa_max
   endwhere
   where (om_ssa(:, :) < small_oceanorg)
      om_ssa(:, :) = 0.0_r8
   endwhere

end subroutine omfrac_accu_aitk


!-------------------------------------------------------------------
!! READ INPUT FILES, CREATE FIELDS, and horizontally interpolate ocean data
!-------------------------------------------------------------------
subroutine init_ocean_data()
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: register advected constituents for all aerosols
    ! 
    ! Method: 
    ! <Describe the algorithm(s) used in the routine.> 
    ! <Also include any applicable external references.> 
    ! 
    ! Author: S. M. Burrows, adapted from dust_initialize
    ! 
    !-----------------------------------------------------------------------

    use tracer_data,      only : trcdata_init
    use cam_history,      only : addfld, horiz_only, add_default
    use spmd_utils,       only : masterproc
    use sslt_sections,    only : nsections

    !-----------------------------------------------------------------------
    ! 	... local variables
    !-----------------------------------------------------------------------

!    type(interp_type)     :: lon_wgts, lat_wgts
!    real(r8), parameter   :: zero=0._r8, twopi=2._r8*pi

    integer :: i, m, m_om
    integer :: number_flds

    if ( masterproc ) then
       write(iulog,*) 'ocean organics are prescribed in :'//trim(filename)
    endif

!    allocate (file%in_pbuf(size(specifier)))
    allocate (file%in_pbuf(n_ocean_data))
    file%in_pbuf(:) = .false.
!    file%in_pbuf(:) = .true.
!
!       fields(i)%pbuf_ndx = pbuf_get_index(fields(i)%fldnam,errcode)

    call trcdata_init( specifier, filename, filelist, datapath, fields, file, &
                       rmv_file, data_cycle_yr, fixed_ymd, fixed_tod, datatype)

    number_flds = 0
    if (associated(fields)) number_flds = size( fields )

    if ( number_flds .eq. n_ocean_data ) then
       if ( masterproc ) then
          write(iulog,"(A21,I3,A20)") 'Successfully read in ',number_flds,' ocean data fields'
       endif
    else if( number_flds < 1 ) then
       if ( masterproc ) then
          write(iulog,*) 'Failed to read in any ocean data'
          write(iulog,*) ' '
       endif
       return
    else if ( number_flds .ne. n_ocean_data ) then
       if ( masterproc ) then
          write(iulog,"(A8,I3,A20)") 'Read in ',number_flds,' ocean data fields'
          write(iulog,"(A9,I3,A20)") 'Expected ',n_ocean_data,' ocean data fields'
          write(iulog,*) ' '
          return
       endif
    end if

    ! Following loop adds fields for output.
    !   Note that the units are given in fields(i)%units, avgflag='A' indicates output mean
    fldloop:do i = 1,n_ocean_data

       if ( masterproc ) then
          write(iulog,*) 'adding field '//fields(i)%fldnam//' ...'
       endif

!!$       ! Set names of variable tendencies and declare them as history variables
!!$       !    addfld(fname,                 unite,              numlev, avgflag, long_name, decomp_type, ...)
       if ( trim(fields(i)%fldnam) == "chla" ) then
          call addfld(trim(fields(i)%fldnam), horiz_only, 'A', 'mg L-1 ', 'ocean input data: '//fields(i)%fldnam ) 
          call add_default (fields(i)%fldnam, 1, ' ')
       else
          call addfld(trim(fields(i)%fldnam), horiz_only, 'A', 'uM C ', 'ocean input data: '//fields(i)%fldnam ) 
          call add_default (fields(i)%fldnam, 1, ' ')
       endif

    enddo fldloop

    if ( masterproc ) then
       write(iulog,*) 'Done initializing marine organics data'
    endif

  end subroutine init_ocean_data

end module seasalt_model

