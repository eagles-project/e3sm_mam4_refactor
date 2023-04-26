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

  logical :: debug_mam_mom = .false.

! Parameters for organic sea salt emissions
    real(r8), parameter :: Aw_carbon = 12.0107_r8       ! Atomic weight oc carbon
    real(r8), parameter :: g_per_m3_NaCl_bulk = 35875._r8 ! approx volume density of salt in seawater
    real(r8) :: g_per_m2_NaCl_bub                       ! g salt per area bubble surface
    integer, parameter  :: n_org_max = 3                ! max number of organic compound classes
    integer  :: n_org                ! actual number of organic compound classes (scheme dependent)
    integer, parameter  :: n_org_burrows = 3        ! actual number of organic compound classes (scheme dependent)
    integer, parameter  :: n_org_gantt   = 1        ! actual number of organic compound classes (scheme dependent)
    integer, parameter  :: n_org_rinaldi = 1        ! actual number of organic compound classes (scheme dependent)
    integer, parameter  :: n_org_quinn   = 1        ! actual number of organic compound classes (scheme dependent)

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
  ! Bubble film thickness
   real(r8)            :: l_bub = 0.1e-6_r8
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
    real(r8), dimension(n_org_burrows), parameter :: & ! OM:OC mass ratios for input fields (mpoly, mprot, mlip)
         OM_to_OC_in = (/ 2.3_r8, 2.2_r8, 1.3_r8 /)
    real(r8), dimension(n_org_burrows), parameter :: & ! Langmuir parameters (inverse C_1/2)  [m3 mol-1]
         alpha_org = (/ 90.58_r8, 25175._r8, 18205._r8 /)
! Molecular weights needed for output of optional diagnostic variable F_eff
    real(r8), dimension(n_org_burrows), parameter :: & ! Molecular weights [g mol-1]
         Mw_org   = (/ 250000._r8, 66463._r8, 284._r8 /)
    real(r8), dimension(n_org_burrows), parameter :: & ! mass per sq. m at saturation
         g_per_m2_org = (/ 0.1376_r8, 0.00219_r8, 0.002593_r8 /) ! Mw_org / a_org

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


!   ! Turn on mam_mom aerosols if user has specified an input dataset.
!   has_mam_mom = len_trim(filename) > 0

end subroutine ocean_data_readnl

  !=============================================================================
subroutine seasalt_emis(lchnk, ncol, u10cubed, srf_temp, ocnfrc, emis_scale, cflx)

    use sslt_sections, only: nsections, fluxes, Dg, rdry
    use mo_constants,  only: dns_aer_sst=>seasalt_density, pi
    use cam_history,   only: outfld
    use spmd_utils,    only: masterproc    

    ! input
    integer, intent(in)  :: lchnk
    integer, intent(in)  :: ncol
    real(r8), intent(in) :: u10cubed(pcols)
    real(r8), intent(in) :: srf_temp(pcols)
    real(r8), intent(in) :: ocnfrc(pcols)
    real(r8), intent(in) :: emis_scale   

    ! output
    real(r8), intent(inout) :: cflx(:,:)

    ! local vars
    integer  :: mn, mm, ibin, i
    real(r8) :: fi(pcols,nsections)
    integer :: m, n
    real(r8):: cflx_help2(pcols)

    fi(:ncol,:nsections) = fluxes( srf_temp, u10cubed, ncol )

   tracer_loop: do ibin = 1,nslt
      ! Index of mass mode
      mm = seasalt_indices(ibin)
      ! Index of number mode
      mn = seasalt_indices(nslt+nslt_om+ibin)

      if (mn>0) then
         ! Total number flux per mode
         section_loop_ssa_num: do i=1, nsections
            cflx_help2(:ncol) = 0.0_r8
            if (Dg(i) >= sst_sz_range_lo(ibin) .and. Dg(i) < sst_sz_range_hi(ibin)) then
               cflx_help2(:ncol)=fi(:ncol,i)*ocnfrc(:ncol)*emis_scale  !++ ag: scale sea-salt

               ! Mixing state 1: external mixture, add OM to mass and number
               ! Mixing state 3: internal mixture, add OM to mass and number
               cflx(:ncol,mn) = cflx(:ncol,mn) + cflx_help2(:ncol)
            endif
         enddo section_loop_ssa_num
      endif

      cflx(:ncol,mm)=0.0_r8
      section_loop_sslt_mass: do i=1, nsections
         if (Dg(i) >= sst_sz_range_lo(ibin) .and. Dg(i) < sst_sz_range_hi(ibin)) then
            cflx_help2(:ncol) = 0.0_r8
            cflx_help2(:ncol)=fi(:ncol,i)*ocnfrc(:ncol)*emis_scale  &   !++ ag: scale sea-salt
                  *4._r8/3._r8*pi*rdry(i)**3*dns_aer_sst  ! should use dry size, convert from number to mass flux (kg/m2/s)

            ! Mixing state 1: external mixture, add OM to mass and number
            ! Mixing state 3: internal mixture, add OM to mass and number
            cflx(:ncol,mm)      = cflx(:ncol,mm)      +cflx_help2(:ncol)
         endif
      enddo section_loop_sslt_mass
   enddo tracer_loop    

end subroutine seasalt_emis


subroutine marine_organic_emis(lchnk, ncol, u10cubed, srf_temp, ocnfrc, emis_scale, cflx, F_eff)

    use sslt_sections, only: nsections, fluxes, Dg, rdry
    use mo_constants,  only: dns_aer_sst=>seasalt_density, pi
    use cam_history,   only: outfld
    use spmd_utils,    only: masterproc

    ! input
    integer, intent(in)  :: lchnk
    integer, intent(in)  :: ncol
    real(r8), intent(in) :: u10cubed(pcols)
    real(r8), intent(in) :: srf_temp(pcols)
    real(r8), intent(in) :: ocnfrc(pcols)
    real(r8), intent(in) :: emis_scale
    
    ! output
    real(r8), intent(inout) :: cflx(:,:)   
    real(r8), intent(out) :: F_eff(pcols)   ! optional diagnostic output

    ! local vars
    integer  :: mn, mm, ibin, i
    real(r8) :: fi(pcols,nsections)
    integer :: m, n
    real(r8):: cflx_help2(pcols)

   real(r8), pointer :: chla(:)          ! for Gantt et al. (2011) organic mass fraction
   real(r8), pointer :: mpoly(:)         ! for Burrows et al. (2014) organic mass fraction
   real(r8), pointer :: mprot(:)         ! for Burrows et al. (2014) organic mass fraction
   real(r8), pointer :: mlip(:)          ! for Burrows et al. (2014) organic mass fraction

   logical :: emit_this_mode(om_num_modes)

   real(r8) :: mass_frac_bub_section(pcols, n_org_max, nsections)
   real(r8) :: om_ssa(pcols, nsections)

   integer  :: m_om ! integer for iteration

   fi(:ncol,:nsections) = fluxes( srf_temp, u10cubed, ncol )

   nullify(chla)
   nullify(mpoly)
   nullify(mprot)
   nullify(mlip)

   fldloop: do i=1,n_ocean_data
      select case (trim(fields(i)%fldnam))
         case ("chla")
             chla   => fields(i)%data(:ncol,1,lchnk)
         case ("mpoly")
             mpoly  => fields(i)%data(:ncol,1,lchnk)
         case ("mprot")
             mprot  => fields(i)%data(:ncol,1,lchnk)
         case ("mlip")
             mlip   => fields(i)%data(:ncol,1,lchnk)
         case default
             if ( masterproc ) then
                write(iulog,*) 'Unknown field name '//fields%fldnam//' in ocean_data fields ...'
             endif
      end select
   enddo fldloop

   mass_frac_bub_section(:ncol,:,:) = 0.0_r8
   om_ssa(:ncol,:) = 0.0_r8
   F_eff(:ncol) = 0.0_r8

   !if (fmoa==1) then ! Burrows et al. organic sea spray
   n_org = n_org_burrows
   call calc_om_ssa_burrows(ncol, mpoly(:ncol), mprot(:ncol), mlip(:ncol), &
                            mass_frac_bub_section(:ncol, :, :), om_ssa(:ncol, :), F_eff(:ncol), lchnk)


! Calculate emission of MOM mass.

! Determine which modes to emit MOM in depending on mixing state assumption

    ! OM modes (m in this loop)
    ! m=1 : accu       (internal w/ SS)
    ! m=2 : Aitken     (internal w/ SS)
    ! m=3 : accu MOM   (external)
    ! m=4 : Aitken MOM (external)

    ! Total external mixture: emit only in modes 4, 5
    ! Mixing state 0: external mixture, replace mass and number
    !                 of mode with mass and number in MOM modes
    ! Mixing state 1: external mixture, add OM to mass and number

   emit_this_mode = (/ .true., .true., .false. /)

   ! Loop over OM modes
   om_num_mode_loop: do m_om=1,om_num_modes ! modes in which to emit OM
      if (size(om_num_ind) .eq. 1) then
         call endrun(&
         "Error: om_num_ind is a scalar, but attempting to calculate MOM.  Something bad happened!!  We should never get here!")
      endif
      
      m = om_num_ind(m_om)
      mn=seasalt_indices(nslt+nslt_om+m)

      ! add number tracers for organics-only modes
      if (emit_this_mode(m_om)) then
         section_loop_OM_num: do i=1, nsections
            cflx_help2(:ncol) = 0.0_r8
            if (Dg(i) >= sst_sz_range_lo(nslt+m_om) .and. Dg(i) < sst_sz_range_hi(nslt+m_om)) then
               cflx_help2(:ncol)=fi(:ncol,i)*ocnfrc(:ncol)*emis_scale

               ! Mixing state 3: internal mixture, add OM to mass and number
               cflx(:ncol,mn) = cflx(:ncol,mn) + cflx_help2(:ncol) * &
                               (1._r8 / (1._r8 - om_ssa(:ncol, i)) - 1._r8)
            endif
         enddo section_loop_OM_num
      endif
   enddo om_num_mode_loop

   om_mode_loop: do m_om=1,nslt_om
      mm = seasalt_indices(nslt+m_om)

      cflx(:ncol,mm)=0.0_r8
      
      if (emit_this_mode(m_om)) then
         ! write(iulog,"(A30,A10,I3)") "Constituent name and number: ", trim(seasalt_names(nslt+m_om)), mm ! for debugging
         ! add mass tracers
         om_type_loop: do n=1,n_org
            section_loop_OM_mass: do i=1, nsections
               if (Dg(i)>=sst_sz_range_lo(nslt+m_om) .and. Dg(i)<sst_sz_range_hi(nslt+m_om)) then
                  cflx_help2(:ncol)=fi(:ncol,i)*ocnfrc(:ncol)*emis_scale &
                                    *4._r8/3._r8*pi*rdry(i)**3*dns_aer_sst  ! should use dry size, convert from number to mass flux (kg/m2/s)
                  !  mass_frac_bub_section(pcols, n_org_max, nsections) -- org classes in dim 2, size nsections in dim 3
             
                  ! Mixing state 1: external mixture, add OM to mass and number
                  ! Mixing state 3: internal mixture, add OM to mass and number
                  where (om_ssa(:ncol, i) > 0.0_r8) ! avoid division by zero
                     cflx(:ncol,mm) = cflx(:ncol,mm) + cflx_help2(:ncol) &
                                    * mass_frac_bub_section(:ncol, n, i) / om_ssa(:ncol, i) * &
                                      (1._r8 / (1._r8 - om_ssa(:ncol, i)) - 1._r8)
                  elsewhere
                     cflx(:ncol,mm) = cflx(:ncol,mm)
                  endwhere
               endif
            enddo section_loop_OM_mass
         enddo om_type_loop
      endif
   enddo om_mode_loop

end subroutine marine_organic_emis


subroutine calc_om_ssa_burrows(ncol, mpoly_in, mprot_in, mlip_in, &
                                mass_frac_bub_section, om_ssa, F_eff, lchnk)

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
   real(r8), intent(inout) :: mass_frac_bub_section(:,:,:)
   real(r8), intent(inout) :: om_ssa(:,:)
   real(r8), intent(inout) :: F_eff(:) ! optional diagnostic output
   !
   ! Local variables
   real(r8) :: g_per_m3(ncol, n_org_burrows), mol_per_m3(ncol, n_org_burrows)
   real(r8) :: theta(ncol, n_org_burrows), alpha_help(ncol)
   real(r8) :: mass_frac_bub(ncol, n_org_burrows), mass_frac_bub_tot(ncol)
   real(r8) :: theta_help(ncol, n_org_burrows), mass_frac_bub_help(ncol, n_org_burrows)
   real(r8), parameter :: particle_size_for_OMF_param = 0.5_r8 ! in microns
   !
   ! OMF maximum and minimum values -- max from Rinaldi et al. (2013)
   real(r8), parameter :: omfrac_max = 0.78
   !
   integer  :: i
   integer  :: lchnk
   !
   !-----------------------------------------------------------------------

   ! Initialize arrays to zero for safety.
   theta(:,:) = 0.0_r8
   theta_help(:,:) = 0.0_r8
   mass_frac_bub_tot(:) = 0.0_r8
   mass_frac_bub(:,:) = 0.0_r8
   mass_frac_bub_help(:,:) = 0.0_r8
   mol_per_m3(:,:) = 0.0_r8
   g_per_m3(:,:) = 0.0_r8
   F_eff(:) = 0.0_r8

   ! Convert input fields from [(mol C) L-1] to [(g OM) m-3]
   ! and store in single array
   g_per_m3(:, 1) = mpoly_in(:)  * 1.0e-3_r8 * OM_to_OC_in(1) * Aw_carbon
   g_per_m3(:, 2) = mprot_in(:)  * 1.0e-3_r8 * OM_to_OC_in(2) * Aw_carbon
   g_per_m3(:, 3) = mlip_in(:)   * 1.0e-3_r8 * OM_to_OC_in(3) * Aw_carbon


   ! Calculate the surface coverage by class
   do i=1,n_org_burrows
      ! Bulk mass concentration [mol m-3] = [g m-3] / [g mol-1]
      mol_per_m3(:, i)     = g_per_m3(:, i) / Mw_org(i)
      ! use theta_help as work array -- theta_help = alpha(i) * x(i)
      theta_help(:, i)         = alpha_org(i)*mol_per_m3(:, i)
   enddo
   alpha_help(:) = sum(theta, dim=2)
   
   do i=1,n_org_burrows
      ! complete calculation -- theta = alpha(i) * x(i) / (1 + sum( alpha(i) * x(i) ))
      theta(:, i)         = theta_help(:, i) / (1.0_r8 + alpha_help(:))
      ! Calculate the organic mass per area (by class) [g m-2]
      !  (use mass_frac_bub_help as local work array -- organic mass per area in g per m2)
      mass_frac_bub_help(:, i) = theta(:, i) * g_per_m2_org(i)
   enddo

   ! Calculate g NaCl per m2
   g_per_m2_NaCl_bub = g_per_m3_NaCl_bulk*l_bub ! Redundant, but allows for easier adjustment to l_bub

   ! mass_frac_bub = 2*[g OM m-2] / (2*[g OM m-2] * [g NaCl m-2])
   ! Factor 2 for bubble bilayer (coated on both surfaces of film)
   do i=1,n_org_burrows
      mass_frac_bub(:, i) = 2.0_r8*mass_frac_bub_help(:, i) / &
           (2.0_r8*sum(mass_frac_bub_help(:, :), dim=2) + g_per_m2_NaCl_bub)
   enddo

   mass_frac_bub_tot(:) = sum(mass_frac_bub, dim=2)

   do i=1,n_org_burrows
      where (mass_frac_bub_tot(:) > omfrac_max)
         mass_frac_bub(:, i) = mass_frac_bub(:, i) / mass_frac_bub_tot(:) * omfrac_max
      endwhere
   enddo

   ! Must exceed threshold value small_oceanorg
   where (mass_frac_bub(:, :) < small_oceanorg)
      mass_frac_bub(:, :) = 0.0_r8
   endwhere

   mass_frac_bub_tot(:) = sum(mass_frac_bub, dim=2)


   ! Effective mass enrichment ratio (for bulk OM) -- diagnostic variable
   !
   ! F_eff_mass = organic mass per area (film) [g m-2] / salt mass per area (film) [g m-2] *
   !                 bulk salt concentration [g m-3] / bulk OM concentration [g m-3]

   if ( F_eff_out ) then
      F_eff(:) = sum(g_per_m3(:, :), dim=2)
      ! avoid division by zero
      where( (mass_frac_bub_tot(:) > small_oceanorg) .and. (F_eff(:) > small_oceanorg) )
         F_eff(:) = 2.0_r8*sum(mass_frac_bub(:, :), dim=2) / g_per_m2_NaCl_bub * &
                 g_per_m3_NaCl_bulk / F_eff(:)
      elsewhere
         F_eff(:) = 0.0_r8
      endwhere
   endif

   ! Distribute mass fraction evenly into Aitken and accumulation modes
   call omfrac_accu_aitk(mass_frac_bub_tot(:), om_ssa(:,:))

   ! mass_frac_bub_section(pcols, n_org_max, nsections) -- org classes in dim 2, size nsections in dim 3
   mass_frac_bub_section(:, :, :)   = 0.0_r8

   do i=1,n_org_burrows
      call omfrac_accu_aitk(mass_frac_bub(:, i), mass_frac_bub_section(:, i, :))
   enddo

end subroutine calc_om_ssa_burrows


subroutine omfrac_accu_aitk(om_ssa_in, om_ssa)

   !----------------------------------------------------------------------- 
   ! Purpose:
   ! Put OM fraction directly into aitken and accumulation modes and don't
   ! exceed om_ssa_max.
   use sslt_sections, only: nsections, Dg
   implicit none
   !-----------------------------------------------------------------------
   ! Arguments:
   !
   real(r8), intent(in)    :: om_ssa_in(:)
   real(r8), intent(inout) :: om_ssa(:,:)
   !
   ! Local variables
   !
   integer  :: m
   real(r8), parameter :: om_ssa_max = 1.0_r8
   !
   !-----------------------------------------------------------------------

   ! distribute OM fraction!
   do m=1,nsections
      ! update only in Aitken and accu. modes
      if ((Dg(m) >= sst_sz_range_lo(2)) .and. (Dg(m) < sst_sz_range_hi(1))) then
             om_ssa(:, m) = om_ssa_in(:)
      else
          om_ssa(:, m) = 0.0_r8 ! Set to zero for "fine sea salt" and "coarse sea salt" modes
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

