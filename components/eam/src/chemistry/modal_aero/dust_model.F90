!===============================================================================
! Dust for Modal Aerosol Model
!===============================================================================
module dust_model 
  use shr_kind_mod, only: r8 => shr_kind_r8, cl => shr_kind_cl
  use spmd_utils,   only: masterproc
  use cam_abortutils,   only: endrun

  implicit none
  private

  public :: dust_names
  public :: dust_nbin
  public :: dust_nnum
  public :: dust_indices
  public :: dust_emis
  public :: dust_readnl
  public :: dust_init
  public :: dust_active

  integer, parameter :: dust_nbin = 2
  integer, parameter :: dust_nnum = 2

  character(len=6), parameter :: dust_names(dust_nbin+dust_nnum) = (/ 'dst_a1', 'dst_a3', 'num_a1', 'num_a3' /)
  real(r8),         parameter :: dust_dmt_grd(dust_nbin+1) = (/ 0.1e-6_r8, 1.0e-6_r8, 10.0e-6_r8/)
! Kok11: fractions of bin (0.1-1) and bin (1-10) in size 0.1-10
  real(r8),         parameter :: dust_emis_sclfctr(dust_nbin) = (/ 0.011_r8,0.989_r8 /)

  integer  :: dust_indices(dust_nbin+dust_nnum)
  real(r8) :: dust_dmt_vwr(dust_nbin)
  real(r8) :: dust_stk_crc(dust_nbin)

  real(r8)          :: dust_emis_fact = -1.e36_r8        ! tuning parameter for dust emissions
  character(len=cl) :: soil_erod_file = 'soil_erod_file' ! full pathname for soil erodibility dataset

  logical :: dust_active = .false.

 contains

  !=============================================================================
  ! reads dust namelist options
  !=============================================================================
  subroutine dust_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'dust_readnl'

    namelist /dust_nl/ dust_emis_fact, soil_erod_file

    !-----------------------------------------------------------------------------

    ! Read namelist
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'dust_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, dust_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast(dust_emis_fact, 1,                   mpir8,   0, mpicom)
    call mpibcast(soil_erod_file, len(soil_erod_file), mpichar, 0, mpicom)
#endif

  end subroutine dust_readnl

  !=============================================================================
  !=============================================================================
  subroutine dust_init()
    use soil_erod_mod, only: soil_erod_init
    use constituents,  only: cnst_get_ind
    use dust_common,   only: dust_set_params

    integer :: n

    do n = 1, dust_nbin
       call cnst_get_ind(dust_names(n), dust_indices(n),abrtf=.false.)
    end do
    do n = 1, dust_nnum
       call cnst_get_ind(dust_names(dust_nbin+n), dust_indices(dust_nbin+n),abrtf=.false.)
    enddo 
    dust_active = any(dust_indices(:) > 0)
    if (.not.dust_active) return
   
    call  soil_erod_init( dust_emis_fact, soil_erod_file )

    call dust_set_params( dust_nbin, dust_dmt_grd, dust_dmt_vwr, dust_stk_crc )

  end subroutine dust_init

  !===============================================================================
  !===============================================================================
  subroutine dust_emis( ncol, lchnk, dust_flux_in, &  ! in
                        cflx, &                       ! inout
                        soil_erod )                   ! out
    use soil_erod_mod, only : soil_erod_fact
    use soil_erod_mod, only : soil_erodibility
    use mo_constants,  only : dust_density
    use physconst,     only : pi

  ! args
    integer,  intent(in)    :: ncol, lchnk
    real(r8), intent(in)    :: dust_flux_in(:,:)
    real(r8), intent(inout) :: cflx(:,:)
    real(r8), intent(out)   :: soil_erod(:)

  ! local vars
    integer :: icol, ibin, idx_dst, inum
    real(r8) :: x_mton
    real(r8), parameter :: soil_erod_threshold = 0.1_r8
    real(r8) :: dst_mass_to_num(dust_nbin)
    
    do ibin = 1, dust_nbin
       dst_mass_to_num(ibin) = 6._r8 / (pi * dust_density * (dust_dmt_vwr(ibin)**3._r8))
    enddo

    ! set dust emissions

    col_loop: do icol =1,ncol

       soil_erod(icol) = soil_erodibility( icol, lchnk )

       if( soil_erod(icol) .lt. soil_erod_threshold ) soil_erod(icol) = 0._r8

       ! rebin and adjust dust emissons..
       do ibin = 1,dust_nbin

          idx_dst = dust_indices(ibin)

       ! Correct the dust input flux calculated by CLM, which uses size distribution in Zender03
       ! to calculate fraction of bin (0.1-10um) in range (0.1-20um) = 0.87
       ! based on Kok11, that fraction is 0.73
          cflx(icol,idx_dst) = sum( -dust_flux_in(icol,:) ) * 0.73_r8/0.87_r8 &
                             * dust_emis_sclfctr(ibin)*soil_erod(icol)/soil_erod_fact*1.15_r8

          inum = dust_indices(ibin+dust_nbin)

          cflx(icol,inum) = cflx(icol,idx_dst)*dst_mass_to_num(ibin)

       enddo

    enddo col_loop

  end subroutine dust_emis

end module dust_model
