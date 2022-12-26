module modal_aero_amicphys_diags

  use modal_aero_amicphys_control, only: ncnst=>gas_pcnst

  implicit none

  private

  ! Public interfaces

  public :: amicphys_diags_init
  public :: get_gcm_tend_diags_from_subareas
  public :: accumulate_column_tend_integrals
  public :: outfld_1proc_all_cnst

  ! Public parameters

  integer,parameter,public :: nqtendaa = 4
  integer,parameter,public :: iqtend_cond = 1
  integer,parameter,public :: iqtend_rnam = 2
  integer,parameter,public :: iqtend_nnuc = 3
  integer,parameter,public :: iqtend_coag = 4

  integer,parameter,public :: nqqcwtendaa = 1
  integer,parameter,public :: iqqcwtend_rnam = 1

  character(len=8),parameter,public :: suffix_q_coltendaa(nqtendaa) = (/'_sfgaex1','_sfgaex2','_sfnnuc1','_sfcoag1'/)
  character(len=8),parameter,public :: suffix_qqcw_coltendaa(nqqcwtendaa) = '_sfgaex2'

  ! Public variables

  logical,public :: do_q_coltendaa(ncnst,nqtendaa) = .false.
  logical,public :: do_qqcw_coltendaa(ncnst,nqqcwtendaa) = .false.

contains

!-----------------------------------------------------------------------------------
! Purpose: set switches that turn on or off column-integrated tendency diagnostics.
!-----------------------------------------------------------------------------------
subroutine amicphys_diags_init( do_cond, do_rename, do_newnuc, do_coag )

  logical,intent(in) :: do_cond, do_rename, do_newnuc, do_coag

  if ( (.not. do_cond) .and. (.not. do_rename) ) then
     do_q_coltendaa(:,iqtend_cond) = .false.
     do_q_coltendaa(:,iqtend_rnam) = .false.
     do_qqcw_coltendaa(:,iqqcwtend_rnam) = .false.
  end if
  if (.not.do_newnuc) do_q_coltendaa(:,iqtend_nnuc) = .false.
  if (.not.do_coag)   do_q_coltendaa(:,iqtend_coag) = .false.

end subroutine amicphys_diags_init

!--------------------------------------------------------------------------------
! Purpose: Get grid cell mean tendencies by calculating area-weighted averages
!          of the values in different subareas.
!--------------------------------------------------------------------------------
subroutine get_gcm_tend_diags_from_subareas( nsubarea, ncldy_subarea, afracsub, &! in
                                             qsub_tendaa, qqcwsub_tendaa,       &! in
                                             qgcm_tendaa, qqcwgcm_tendaa        )! out

  use modal_aero_amicphys_control, only: wp=>r8, ncnst=>gas_pcnst, maxsubarea

  integer,  intent(in) :: nsubarea                   ! # of active subareas
  integer,  intent(in) :: ncldy_subarea              ! # of cloudy subareas
  real(wp), intent(in) :: afracsub(maxsubarea)       ! area fraction of subareas [unitless]

  real(wp),intent(in) ::    qsub_tendaa(ncnst,   nqtendaa,maxsubarea)
  real(wp),intent(in) :: qqcwsub_tendaa(ncnst,nqqcwtendaa,maxsubarea)

  real(wp),intent(out) ::    qgcm_tendaa(ncnst,   nqtendaa)
  real(wp),intent(out) :: qqcwgcm_tendaa(ncnst,nqqcwtendaa)

  integer :: jsub

  ! Gases and interstitial aerosols

  qgcm_tendaa(:,:) = 0.0_wp
  do jsub = 1, nsubarea
     qgcm_tendaa(:,:) = qgcm_tendaa(:,:) &
                      + qsub_tendaa(:,:,jsub)*afracsub(jsub)
  end do

  ! Cloud-borne aerosols

  if (ncldy_subarea <= 0) then
     qqcwgcm_tendaa(:,:) = 0.0_wp
  else
     qqcwgcm_tendaa(:,:) = 0.0_wp
     do jsub = 1, nsubarea
        qqcwgcm_tendaa(:,:) = qqcwgcm_tendaa(:,:) &
                            + qqcwsub_tendaa(:,:,jsub)*afracsub(jsub)
     end do
  end if

end subroutine get_gcm_tend_diags_from_subareas

!-------------------------------------------------------------------------------
! Purpose: add mass-weighted tendencies to the corresponding vertical integrals
!-------------------------------------------------------------------------------
subroutine accumulate_column_tend_integrals( pdel, gravit,                &! in
                                             qgcm_tendaa, qqcwgcm_tendaa, &! in
                                             q_coltendaa, qqcw_coltendaa  )! inout

  use modal_aero_amicphys_control, only: wp=>r8, ncnst=>gas_pcnst

  real(wp),intent(in) :: pdel       ! layer thickness in pressure coordinate [Pa]
  real(wp),intent(in) :: gravit     ! gravitational acceleration [m/s^2]

  real(wp),intent(in)    ::    qgcm_tendaa(ncnst,   nqtendaa)
  real(wp),intent(in)    :: qqcwgcm_tendaa(ncnst,nqqcwtendaa)

  real(wp),intent(inout) ::    q_coltendaa(ncnst,   nqtendaa)
  real(wp),intent(inout) :: qqcw_coltendaa(ncnst,nqqcwtendaa)

  real(wp) :: pdel_fac   ! = air density * dz

  integer :: iqtend  ! loop index corresponding to different aerosol processes
  integer :: icnst   ! loop index corresponding to different constituents (tracers)

  pdel_fac = pdel/gravit

  ! Gases and interstitial aerosols

  do iqtend = 1,nqtendaa
  do icnst = 1,ncnst
     if ( do_q_coltendaa(icnst,iqtend) ) then
        q_coltendaa(icnst,iqtend) = q_coltendaa(icnst,iqtend) + qgcm_tendaa(icnst,iqtend)*pdel_fac
     end if
  end do ! icnst
  end do ! iqtend

  ! cloud-borne aerosols

  do iqtend = 1,nqqcwtendaa 
  do icnst = 1,ncnst
     if ( do_qqcw_coltendaa(icnst,iqtend) ) then
        qqcw_coltendaa(icnst,iqtend) = qqcw_coltendaa(icnst,iqtend) + qqcwgcm_tendaa(icnst,iqtend)*pdel_fac
     end if
  end do ! icnst
  end do ! iqtend

end subroutine accumulate_column_tend_integrals

!---------------------------------------------------------------------------------
! Purpose: output vertically integrated tendencies of one process and different
!          constituents (tracers)
!---------------------------------------------------------------------------------
subroutine outfld_1proc_all_cnst( qcoltend, pcols, ntend, itend, lout, &
                                  cnst_name, tendname_suffix, &
                                  mwdry, loffset, ncol, lchnk )

  use cam_history,                 only: outfld, fieldname_len
  use chem_mods,                   only: adv_mass
  use modal_aero_amicphys_control, only: wp=>r8, ncnst=>gas_pcnst

  integer, intent(in)    :: pcols
  integer, intent(in)    :: ntend
  real(wp),intent(inout) :: qcoltend(pcols,ncnst,ntend)

  integer, intent(in)    :: itend
  logical, intent(in)    :: lout(ncnst,ntend)

  character(len=*),intent(in) :: cnst_name(:)
  character(len=*),intent(in) :: tendname_suffix(:)

  real(wp),intent(in)    :: mwdry
  integer, intent(in)    :: loffset
  integer, intent(in)    :: ncol
  integer, intent(in)    :: lchnk

  integer :: icnst
  character(len=fieldname_len+3) :: fieldname

  do icnst = 1,ncnst
    if ( lout(icnst,itend) .and. (itend>0) ) then

       ! Convert unit
       qcoltend(1:ncol,icnst,itend) = qcoltend(1:ncol,icnst,itend)* (adv_mass(icnst)/mwdry)

       ! Send to history output
       fieldname = trim(cnst_name(icnst+loffset))//tendname_suffix(itend)
       call outfld( fieldname, qcoltend(1:ncol,icnst,itend), ncol, lchnk )

    end if
  end do ! icnst

end subroutine outfld_1proc_all_cnst

end module modal_aero_amicphys_diags
