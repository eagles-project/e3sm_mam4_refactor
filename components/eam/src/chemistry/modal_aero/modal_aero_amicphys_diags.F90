module modal_aero_amicphys_diags

  use modal_aero_amicphys_control, only: ncnst=>gas_pcnst

  implicit none

  private

  ! Public interfaces

  public :: amicphys_diags_init
  public :: get_gcm_tend_diags_from_subareas
  public :: accumulate_column_tend_integrals
  public :: outfld_1proc_all_cnst
  public :: m_a_amicphys_init_history

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

  logical,public,protected :: do_q_coltendaa(ncnst,nqtendaa) = .false.
  logical,public,protected :: do_qqcw_coltendaa(ncnst,nqqcwtendaa) = .false.

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


!-----------------------------------------------------------------------
!
! Purpose:
!    set do_adjust and do_aitken flags
!    create history fields for column tendencies associated with
!       modal_aero_calcsize
!
! Author: R. Easter
!
!-----------------------------------------------------------------------
subroutine m_a_amicphys_init_history( loffset )

use cam_history, only  :  addfld, horiz_only, add_default, fieldname_len
use cam_logfile, only  :  iulog
use constituents, only :  pcnst, cnst_get_ind, cnst_name
use spmd_utils, only   :  masterproc
use phys_control,only  :  phys_getopts

use modal_aero_data, only : &
    cnst_name_cw, &
    modeptr_accum, modeptr_aitken, modeptr_pcarbon, modeptr_ufine

use modal_aero_coag, only: n_coagpair, src_mode_coagpair, dest_mode_coagpair, end_mode_coagpair
use modal_aero_amicphys_control

implicit none

!-----------------------------------------------------------------------
! arguments
   integer, intent(in)  :: loffset

!-----------------------------------------------------------------------
! local
   integer  :: iaer, igas, ipair, iok
   integer  :: lmz, lmza, lmzb, lmzc
   integer  :: m
   integer  :: n, na, nb, nc

   real(r8) :: tmp1, tmp2

   character(len=fieldname_len)   :: tmpnamea, tmpnameb
   character(len=fieldname_len+3) :: fieldname
   character(128)                 :: long_name
   character(128)                 :: msg
   character(8)                   :: unit
   character(2)                   :: tmpch2

   logical                        :: history_aerosol      ! Output the MAM aerosol tendencies
   logical                        :: history_verbose      ! produce verbose history output
   logical                        :: history_aerocom    ! Output the aerocom history
   !-----------------------------------------------------------------------


      call phys_getopts( history_aerosol_out = history_aerosol, &
                         history_verbose_out = history_verbose   )

      history_aerocom = .false.
!
! set the do_q_coltendaa
!
      do_q_coltendaa(:,:) = .false.

! gas-->aer condensation and resulting aging
      do igas = 1, ngas
         lmz = lmap_gas(igas)
         if (lmz <= 0) cycle
         do_q_coltendaa(lmz,iqtend_cond) = .true.
         iaer = igas
         do n = 1, ntot_amode
            lmz = lmap_aer(iaer,n)
            if (lmz <= 0) cycle
            do_q_coltendaa(lmz,iqtend_cond) = .true.
         end do ! n
      end do ! igas

      do ipair = 1, n_agepair
         na = modefrm_agepair(ipair)
         nb = modetoo_agepair(ipair)
         if (na < 1 .or. nb < 1) cycle

         lmza = lmap_num(na)
         lmzb = lmap_num(nb)
         do_q_coltendaa(lmza,iqtend_cond) = .true.
         do_q_coltendaa(lmzb,iqtend_cond) = .true.
         do iaer = 1, naer
            lmza = lmap_aer(iaer,na)
            lmzb = lmap_aer(iaer,nb)
            if (lmza > 0) then
               do_q_coltendaa(lmza,iqtend_cond) = .true.
               if (lmzb > 0) do_q_coltendaa(lmzb,iqtend_cond) = .true.
            end if
         end do ! iaer
      end do ! ipair

!  define history fields for gas-->aer condensation and resulting aging
      do lmz = 1, gas_pcnst
         if ( do_q_coltendaa(lmz,iqtend_cond)) then
            tmpnamea = cnst_name(lmz+loffset)
            fieldname = trim(tmpnamea) // '_sfgaex1'
            long_name = trim(tmpnamea) // ' gas-aerosol-exchange primary column tendency'
            unit = 'kg/m2/s'
            call addfld( fieldname, horiz_only, 'A', unit, long_name )
            if ( history_aerosol ) call add_default( fieldname, 1, ' ' )
            if ( masterproc ) write(iulog,'(3(a,3x))') 'gasaerexch addfld', fieldname, unit
         end if
      end do ! lmz

!  define history fields for 3d soa production for aerocom
      do igas = 1, nsoa
         lmz = lmap_gas(igas)
         if (lmz <= 0) cycle
         if ( .not. do_q_coltendaa(lmz,iqtend_cond)) cycle
         if ( .not. history_aerocom ) cycle

         tmpnamea = cnst_name(lmz+loffset)
         fieldname = trim(tmpnamea) // '_sfgaex3d'
         long_name = trim(tmpnamea) // ' gas-aerosol-exchange primary 3d tendency'
         unit = 'kg/m2/s'
         call addfld( fieldname, (/ 'lev' /), 'A', unit, long_name )
         call add_default( fieldname, 1, ' ' )
         if ( masterproc ) write(iulog,'(3(a,3x),2i5)') &
            'gasaerexch addfld', fieldname, unit, igas, lmz+loffset
      end do


! renaming during gas-->aer condensation or cloud chemistry
      na = modeptr_aitken
      nb = modeptr_accum
      if (na > 0 .and. nb > 0) then
         lmza = lmap_num(na)
         lmzb = lmap_num(nb)
         do_q_coltendaa(lmza,iqtend_rnam) = .true.
         do_q_coltendaa(lmzb,iqtend_rnam) = .true.
         lmza = lmap_numcw(na)
         lmzb = lmap_numcw(nb)
         do_qqcw_coltendaa(lmza,iqqcwtend_rnam) = .true.
         do_qqcw_coltendaa(lmzb,iqqcwtend_rnam) = .true.
         do iaer = 1, naer
            lmza = lmap_aer(iaer,na)
            lmzb = lmap_aer(iaer,nb)
            if (lmza > 0) then
               do_q_coltendaa(lmza,iqtend_rnam) = .true.
               if (lmzb > 0) do_q_coltendaa(lmzb,iqtend_rnam) = .true.
            end if
            lmza = lmap_aercw(iaer,na)
            lmzb = lmap_aercw(iaer,nb)
            if (lmza > 0) then
               do_qqcw_coltendaa(lmza,iqqcwtend_rnam) = .true.
               if (lmzb > 0) do_qqcw_coltendaa(lmzb,iqqcwtend_rnam) = .true.
            end if
         end do ! iaer
      end if ! (na > 0 .and. nb > 0)

!  define history fields for renaming during gas-->aer condensation or cloud chemistry
      do lmz = 1, gas_pcnst
         if ( do_q_coltendaa(lmz,iqtend_rnam)) then
            tmpnamea = cnst_name(lmz+loffset)
            fieldname = trim(tmpnamea) // '_sfgaex2'
            long_name = trim(tmpnamea) // ' gas-aerosol-exchange renaming column tendency'
            unit = 'kg/m2/s'
            if (tmpnamea(1:4) == 'num_' .or. tmpnamea(1:4) == 'NUM_') unit = '#/m2/s'
            call addfld( fieldname, horiz_only, 'A', unit, long_name )
            if ( history_aerosol .and. history_verbose ) call add_default( fieldname, 1, ' ' )
            if ( masterproc ) write(iulog,'(3(a,3x))') 'gasaerexch addfld', fieldname, unit
         end if
         if ( do_qqcw_coltendaa(lmz,iqqcwtend_rnam)) then
            tmpnamea = cnst_name_cw(lmz+loffset)
            fieldname = trim(tmpnamea) // '_sfgaex2'
            long_name = trim(tmpnamea) // ' gas-aerosol-exchange renaming column tendency'
            unit = 'kg/m2/s'
            if (tmpnamea(1:4) == 'num_' .or. tmpnamea(1:4) == 'NUM_') unit = '#/m2/s'
            call addfld( fieldname, horiz_only, 'A', unit, long_name )
            if ( history_aerosol .and. history_verbose ) call add_default( fieldname, 1, ' ' )
            if ( masterproc ) write(iulog,'(3(a,3x))') 'gasaerexch addfld', fieldname, unit
         end if
      end do ! lmz



! coagulation
      do ipair = 1, n_coagpair
         na =  src_mode_coagpair(ipair)
         nb = dest_mode_coagpair(ipair)
         nc =  end_mode_coagpair(ipair)
         if (na < 1 .or. nb < 1 .or. nc < 1) cycle

         lmza = lmap_num(na)
         lmzb = lmap_num(nb)
         lmzc = lmap_num(nc)
         do_q_coltendaa(lmza,iqtend_coag) = .true.
         do_q_coltendaa(lmzb,iqtend_coag) = .true.
         do_q_coltendaa(lmzc,iqtend_coag) = .true.
         do iaer = 1, naer
            lmza = lmap_aer(iaer,na)
            lmzb = lmap_aer(iaer,nb)
            lmzc = lmap_aer(iaer,nc)
            if (lmza > 0) then
               do_q_coltendaa(lmza,iqtend_coag) = .true.
               if (lmzc > 0) do_q_coltendaa(lmzc,iqtend_coag) = .true.
            end if
            if (nb == nc) cycle
            if (lmzb > 0) then
               do_q_coltendaa(lmzb,iqtend_coag) = .true.
               if (lmzc > 0) do_q_coltendaa(lmzc,iqtend_coag) = .true.
            end if
         end do ! iaer
      end do ! ipair

!  define history fields for coagulation
      do lmz = 1, gas_pcnst
         if ( do_q_coltendaa(lmz,iqtend_coag)) then
            tmpnamea = cnst_name(lmz+loffset)
            fieldname = trim(tmpnamea) // '_sfcoag1'
            long_name = trim(tmpnamea) // ' modal_aero coagulation column tendency'
            unit = 'kg/m2/s'
            if (tmpnamea(1:4) == 'num_' .or. tmpnamea(1:4) == 'NUM_') unit = '#/m2/s'
            call addfld( fieldname, horiz_only, 'A', unit, long_name )
            if ( history_aerosol .and. history_verbose ) call add_default( fieldname, 1, ' ' )
            if ( masterproc ) write(iulog,'(3(a,3x))') 'modal_aero_coag_init addfld', fieldname, unit
         end if
      end do ! lmz


! nucleation
      n = modeptr_aitken
      do igas = 1, ngas
         iok = 0
         if (igas == igas_h2so4) iok = 1
         if (igas == igas_nh3  ) iok = 1
         if (iok <= 0) cycle
         lmz = lmap_gas(igas)
         if (lmz > 0) then
            do_q_coltendaa(lmz,iqtend_nnuc) = .true.
            iaer = igas
            lmz = lmap_aer(iaer,n)
            if (lmz > 0) do_q_coltendaa(lmz,iqtend_nnuc) = .true.
          end if
      end do ! igas
      lmzc = lmap_num(n)
      do_q_coltendaa(lmzc,iqtend_nnuc) = .true.

!  define history fields for nucleation
      do lmz = 1, gas_pcnst
         if ( do_q_coltendaa(lmz,iqtend_nnuc)) then
            tmpnamea = cnst_name(lmz+loffset)
            fieldname = trim(tmpnamea) // '_sfnnuc1'
            long_name = trim(tmpnamea) // ' modal_aero new particle nucleation column tendency'
            unit = 'kg/m2/s'
            if (tmpnamea(1:4) == 'num_' .or. tmpnamea(1:4) == 'NUM_') unit = '#/m2/s'
            call addfld( fieldname, horiz_only, 'A', unit, long_name )
            if ( history_aerosol .and. history_verbose ) call add_default( fieldname, 1, ' ' )
            if ( masterproc ) write(iulog,'(3(a,3x))') 'modal_aero_newnuc_init addfld', fieldname, unit
         end if
      end do ! lmz

      if ( history_aerocom ) then
            tmpnamea = cnst_name(lmzc+loffset)
            fieldname = trim(tmpnamea) // '_nuc1'
            long_name = trim(tmpnamea) // ' modal_aero new particle nucleation tendency'
            unit = '#/m3/s'
            call addfld( fieldname, (/ 'lev' /), 'A', unit, long_name )
            call add_default( fieldname, 1, ' ' )
            if ( masterproc ) write(iulog,'(3(a,2x))') &
                'modal_aero_newnuc_init addfld', fieldname, unit

            fieldname = trim(tmpnamea) // '_nuc2'
            long_name = trim(tmpnamea) // ' modal_aero cluster nucleation rate'
            unit = '#/m3/s'
            call addfld( fieldname, (/ 'lev' /), 'A', unit, long_name )
            call add_default( fieldname, 1, ' ' )
            if ( masterproc ) write(iulog,'(3(a,2x))') &
                'modal_aero_newnuc_init addfld', fieldname, unit
      endif

      end subroutine m_a_amicphys_init_history

end module modal_aero_amicphys_diags
