module mam_support
  !---------------------------------------------------------------------
  !Purpose:
  !This module contains utlity variables/subroutines/functions which
  ! multiple MAM routinesuse.
  !---------------------------------------------------------------------

  use shr_kind_mod,     only: r8 => shr_kind_r8

  implicit none

  private ! make everything private

  !explicitly declare public functions/variables
  public :: min_max_bound
  public :: assign_la_lc
  public :: ptr2d_t
  public :: get_cldbrn_mmr

  type ptr2d_t
     real(r8), pointer :: fld(:,:)
  end type ptr2d_t

contains
  !===============================================================================
  subroutine assign_la_lc( imode,      ispec,          & ! in
       la,         lc,             & ! out
       is_lc_append_in             ) ! optional in
    !-----------------------------------------------------------------------
    ! get the index of interstital (la) and cloudborne (lc) aerosols
    ! from mode index and species index
    ! is_lc_append_in is true when cloudborne aerosols are appended after 
    ! interstitial aerosol array (default is false)
    !-----------------------------------------------------------------------
    use constituents,    only: pcnst
    use modal_aero_data, only: lmassptr_amode, lmassptrcw_amode, &
         numptr_amode, numptrcw_amode

    integer, intent(in)     :: imode            ! index of MAM4 modes
    integer, intent(in)     :: ispec            ! index of species, in which:
    ! 0 = number concentration
    ! other = mass concentration
    integer, intent(out)    :: la               ! index of interstitial aerosol
    integer, intent(out)    :: lc               ! index of cloudborne aerosol
    logical, optional, intent(in) :: is_lc_append_in  ! if cloudborne aerosol is appended after interstitial aerosols

    logical :: is_lc_append   

    ! the default option is treat cloudborne aerosols in separated array
    is_lc_append = .false.
    if (present(is_lc_append_in)) is_lc_append=is_lc_append_in

    if (ispec == 0) then
       la = numptr_amode(imode)
       lc = numptrcw_amode(imode)
    else
       la = lmassptr_amode(ispec,imode)
       lc = lmassptrcw_amode(ispec,imode)
    endif

    ! if true: cloudborne aerosol is append after interstitial aerosol
    if (is_lc_append) then
       lc = lc + pcnst  
    endif

  end subroutine assign_la_lc

  !===============================================================================
  pure function min_max_bound(minlim, maxlim, input) result(bounded)
    !Bound a quantity between a min limit and a max limit
    real(r8), intent(in) :: minlim, maxlim
    real(r8), intent(in) :: input

    !return value
    real(r8) :: bounded

    bounded = max(min(maxlim, input), minlim)

  end function min_max_bound
  !===============================================================================

  subroutine get_cldbrn_mmr(lchnk, pbuf, qqcw)
    !Get MMR for cloud borne aerosols using qqcw_get_field function
    use constituents,    only: pcnst
    use physics_buffer,  only: physics_buffer_desc
    use modal_aero_data, only: qqcw_get_field

    !intent-ins
    integer, intent(in) :: lchnk
    type(physics_buffer_desc), pointer :: pbuf(:)

    !intent - out
    type(ptr2d_t), intent(out) :: qqcw(pcnst) !cloud-borne aerosols mass and number mixing ratios

    !local
    integer :: icnst

    integer, parameter :: AER_START_IND = 16 ! starting index of aerosols in MAM4

    do icnst = AER_START_IND, pcnst
       qqcw(icnst)%fld => qqcw_get_field(pbuf,icnst,lchnk)
    enddo

  end subroutine get_cldbrn_mmr
  !===============================================================================

end module mam_support
