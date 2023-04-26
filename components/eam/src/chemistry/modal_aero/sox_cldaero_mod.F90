!----------------------------------------------------------------------------------
! Modal aerosol implementation
!----------------------------------------------------------------------------------
module sox_cldaero_mod
#include "../yaml/common_files/common_uses.ymlf90"

  use shr_kind_mod,    only : r8 => shr_kind_r8
  use cam_abortutils,      only : endrun
  use ppgrid,          only : pcols, pver
  use mo_chem_utls,    only : get_spc_ndx
  use cldaero_mod,     only : cldaero_conc_t
  use modal_aero_data, only : ntot_amode, modeptr_accum, lptr_so4_cw_amode, lptr_msa_cw_amode
  use modal_aero_data, only : numptrcw_amode, lptr_nh4_cw_amode
  use modal_aero_data, only : cnst_name_cw, specmw_so4_amode
  use cam_history,     only : outfld
  use cam_history,     only : addfld, horiz_only, add_default
  use chem_mods,       only : adv_mass
  use physconst,       only : gravit
  use phys_control,    only : phys_getopts
  use cldaero_mod,     only : cldaero_uptakerate
  use chem_mods,       only : gas_pcnst
use spmd_utils,   only : masterproc

  implicit none
  private

  public :: sox_cldaero_init
  public :: sox_cldaero_create_obj
  public :: sox_cldaero_update
  public :: sox_cldaero_destroy_obj

  integer :: id_msa, id_h2so4, id_so2, id_h2o2, id_nh3

contains

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------

  subroutine sox_cldaero_init

    integer :: l, m
    logical :: history_aerosol      ! Output the MAM aerosol tendencies
    logical :: history_verbose      ! produce verbose history output

    id_msa = get_spc_ndx( 'MSA' )
    id_h2so4 = get_spc_ndx( 'H2SO4' )
    id_so2 = get_spc_ndx( 'SO2' )
    id_h2o2 = get_spc_ndx( 'H2O2' )
    id_nh3 = get_spc_ndx( 'NH3' )

    if (id_h2so4<1 .or. id_so2<1 .or. id_h2o2<1) then
      call endrun('sox_cldaero_init:MAM mech does not include necessary species' &
                  //' -- should not invoke sox_cldaero_mod ')
    endif

    call phys_getopts( history_aerosol_out        = history_aerosol, &
                       history_verbose_out        = history_verbose  )
    !
    !   add to history
    !
    do m = 1, ntot_amode

       l = lptr_so4_cw_amode(m)
       if (l > 0) then
          call addfld (&
               trim(cnst_name_cw(l))//'AQSO4',horiz_only,  'A','kg/m2/s', &
               trim(cnst_name_cw(l))//' aqueous phase chemistry')
          call addfld (&
               trim(cnst_name_cw(l))//'AQH2SO4',horiz_only,  'A','kg/m2/s', &
               trim(cnst_name_cw(l))//' aqueous phase chemistry')
          if ( history_aerosol .and. history_verbose ) then 
             call add_default (trim(cnst_name_cw(l))//'AQSO4', 1, ' ')
             call add_default (trim(cnst_name_cw(l))//'AQH2SO4', 1, ' ')
          endif
       end if

    end do

    call addfld ('AQSO4_H2O2',horiz_only,  'A','kg/m2/s', &
         'SO4 aqueous phase chemistry due to H2O2')
    call addfld ('AQSO4_O3',horiz_only,  'A','kg/m2/s', &
         'SO4 aqueous phase chemistry due to O3')

    if ( history_aerosol .and. history_verbose) then    
       call add_default ('AQSO4_H2O2', 1, ' ')
       call add_default ('AQSO4_O3', 1, ' ')    
    endif
  
  end subroutine sox_cldaero_init

!===================================================================================
  function sox_cldaero_create_obj(cldfrc, qcw, lwc, cfact, ncol, loffset) result( conc_obj )
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
    use cldaero_mod, only : cldaero_allocate

    ! input variables    
    real(r8), intent(in) :: cldfrc(:,:) ! cloud fraction [fraction]
    real(r8), intent(in) :: qcw(:,:,:)  ! cloud-borne aerosol [vmr]
    real(r8), intent(in) :: lwc(:,:)    ! cloud liquid water content [kg/kg]
    real(r8), intent(in) :: cfact(:,:)  ! total atms density [kg/L]
    integer,  intent(in) :: ncol
    integer,  intent(in) :: loffset     ! # of tracers in the host model that are not part of MAM
    ! output variables
    type(cldaero_conc_t), pointer :: conc_obj

    ! local variable indexes
    integer :: id_so4_1a, id_so4_2a, id_so4_3a
    integer :: icol,kk
#include "../yaml/sox_cldaero_mod/f90_yaml/sox_cldaero_create_obj_beg_yml.f90"

    conc_obj => cldaero_allocate()

    do kk = 1,pver
       do icol = 1,ncol
          ! conc_obj%xlwc is in-cloud LWC with the unit of [kg/L]
          if( cldfrc(icol,kk) >0._r8) then
             ! cloud water L(water)/L(air)
             conc_obj%xlwc(icol,kk) = lwc(icol,kk) *cfact(icol,kk)
             ! liquid water in the cloudy fraction of cell
             conc_obj%xlwc(icol,kk) = conc_obj%xlwc(icol,kk) / cldfrc(icol,kk)
          else
             conc_obj%xlwc(icol,kk) = 0._r8
          endif
       enddo
    enddo

    conc_obj%no3c(:,:) = 0._r8

    ! FORTRAN refactor: remove code of MAM7 and MAM9
    id_so4_1a = lptr_so4_cw_amode(1) - loffset
    id_so4_2a = lptr_so4_cw_amode(2) - loffset
    id_so4_3a = lptr_so4_cw_amode(3) - loffset
    conc_obj%so4c(:ncol,:) = qcw(:,:,id_so4_1a) + qcw(:,:,id_so4_2a) + qcw(:,:,id_so4_3a)

    ! for 3-mode, so4 is assumed to be nh4hso4
    ! the partial neutralization of so4 is handled by using a 
    !    -1 charge (instead of -2) in the electro-neutrality equation
    conc_obj%nh4c(:ncol,:) = 0._r8

    ! with 3-mode, assume so4 is nh4hso4, and so half-neutralized
    conc_obj%so4_fact = 1._r8

#include "../yaml/sox_cldaero_mod/f90_yaml/sox_cldaero_create_obj_end_yml.f90"

  end function sox_cldaero_create_obj

!=================================================================================
  subroutine sox_cldaero_update( ncol, lchnk, loffset,  & ! in
                dtime, mbar, pdel, press, tfld, cldnum, cldfrc, cfact, xlwc, & ! in
                delso4_hprxn, xh2so4, xso4, xso4_init, & ! in
                qcw, qin ) ! inout
!----------------------------------------------------------------------------------
! Update the mixing ratios
!----------------------------------------------------------------------------------
   
    ! args
    integer,  intent(in) :: ncol
    integer,  intent(in) :: lchnk       ! chunk id
    integer,  intent(in) :: loffset     ! # of tracers in the host model that are not part of MAM
    real(r8), intent(in) :: dtime       ! time step [sec]

    real(r8), intent(in) :: mbar(:,:)   ! mean wet atmospheric mass [amu or g/mol]
    real(r8), intent(in) :: pdel(:,:)   ! pressure interval [Pa] 
    real(r8), intent(in) :: press(:,:)  ! pressure [Pa]
    real(r8), intent(in) :: tfld(:,:)   ! temperature [K]

    real(r8), intent(in) :: cldnum(:,:) ! droplet number concentration [#/kg]
    real(r8), intent(in) :: cldfrc(:,:) ! cloud fraction [fraction]
    real(r8), intent(in) :: cfact(:,:)  ! total atms density [kg/L]
    real(r8), intent(in) :: xlwc(:,:)   ! liquid water volume [cm^3/cm^3]

    real(r8), intent(in) :: delso4_hprxn(:,:)   ! change of so4 due to H2O2 chemistry [mol/mol]
    real(r8), intent(in) :: xh2so4(:,:)         ! H2SO4 mass mixing ratio [mol/mol]
    real(r8), intent(in) :: xso4(:,:)           ! final SO4 mass mixing ratio [mol/mol]
    real(r8), intent(in) :: xso4_init(:,:)      ! initial SO4 mass mixing ratio [mol/mol]

    real(r8), intent(inout) :: qcw(:,:,:) ! cloud-borne aerosol [vmr]
    real(r8), intent(inout) :: qin(:,:,:) ! xported species [vmr]

    ! local vars ...

    ! FORTRAN refactor note: aqueous chemistry (aq) here reprent two processes:
    !       S(IV) + H2O2 = S(VI)
    !       S(IV) + O3   = S(VI)
    ! see the parent subroutine in mo_setsox for reference
    real(r8) :: delso4_o3rxn    ! change of so4 due to O3 chemistry [mol/mol]
    real(r8) :: dso4dt_aqrxn    ! so4_c tendency from aqueous chemistry [mol/mol/s]
    real(r8) :: dso4dt_hprxn    ! so4_c tendency from H2O2 chemistry [mol/mol/s]
    real(r8) :: dso4dt_gasuptk  ! so4_c tendency from h2so4 gas uptake [mol/mol/s]
    real(r8) :: dqdt_aq         ! dqdt due to aqueous chemistry [mol/mol/s]
    real(r8) :: dqdt_wr         ! dqdt due to wet removal, currently set as zero [mol/mol/s]
    real(r8) :: dqdt_aqso4(ncol,pver,gas_pcnst)    ! dqdt due to so4 aqueous chemistry [mol/mol/s]
    real(r8) :: dqdt_aqh2so4(ncol,pver,gas_pcnst)  ! dqdt due to h2so4 uptake [mol/mol/s]
    real(r8) :: dqdt_aqhprxn(ncol,pver)         ! dqdt due to H2O2 chemistry [mol/mol/s]
    real(r8) :: dqdt_aqo3rxn(ncol,pver)         ! dqdt due to O3 chemistry [mol/mol/s]
    real(r8) :: sflx(ncol)      ! integrated surface fluxes [kg/m2/s]
    real(r8) :: faqgain_so4(ntot_amode) ! factor of TMR among modes [fraction]
    real(r8) :: uptkrate        ! uptake rate [1/s]

    integer :: ll, mm, imode    ! aerosol mode index
    integer :: icol,kk          ! column and level index

    real(r8), parameter :: small_value_8  = 1.e-8_r8
    real(r8), parameter :: small_value_5  = 1.e-5_r8
#include "../yaml/sox_cldaero_mod/f90_yaml/sox_cldaero_update_beg_yml.f90"

    ! make sure dqdt is zero initially, for budgets
    dqdt_aqso4(:,:,:)   = 0.0_r8
    dqdt_aqh2so4(:,:,:) = 0.0_r8
    dqdt_aqhprxn(:,:)   = 0.0_r8
    dqdt_aqo3rxn(:,:)   = 0.0_r8

    lev_loop: do kk = 1,pver
       col_loop: do icol = 1,ncol
          cloud: if ((cldfrc(icol,kk) >= small_value_5) .and. (xlwc(icol,kk) >= small_value_8)) then

            !-------------------------------------------------------------------------
            ! compute factors for partitioning aerosol mass gains among modes
            ! the factors are proportional to the activated particle MR for each
            ! mode, which is the MR of cloud drops "associated with" the mode
            ! thus we are assuming the cloud drop size is independent of the
            ! associated aerosol mode properties (i.e., drops associated with
            ! Aitken and coarse sea-salt particles are same size)
            call  compute_aer_factor(qcw(icol,kk,:), loffset,   & ! in
                                     faqgain_so4                ) ! out

            uptkrate = cldaero_uptakerate(xlwc(icol,kk), cldnum(icol,kk), &
                        cfact(icol,kk), cldfrc(icol,kk), tfld(icol,kk),  press(icol,kk))
            ! average uptake rate over dtime
            uptkrate = (1.0_r8 - exp(-min(100._r8,dtime*uptkrate))) / dtime
            ! dso4dt_gasuptk = so4_c tendency from h2so4 gas uptake (mol/mol/s)
            dso4dt_gasuptk = xh2so4(icol,kk) * uptkrate

            delso4_o3rxn = xso4(icol,kk) - xso4_init(icol,kk)
            dso4dt_aqrxn = (delso4_o3rxn + delso4_hprxn(icol,kk)) / dtime
            dso4dt_hprxn = delso4_hprxn(icol,kk) / dtime

            !-----------------------------------------------------------------------
            ! now compute TMR tendencies
            ! this includes the above aqueous so2 chemistry AND
            ! the uptake of highly soluble aerosol precursor gases (h2so4, ...)
            ! The wetremoval of dissolved, unreacted so2 and h2o2 are assumed as zero

            ! compute TMR tendencies for so4 aerosol-in-cloud-water
            do imode = 1, ntot_amode
               ll = lptr_so4_cw_amode(imode) - loffset
               if (ll > 0) then
                  dqdt_aqso4(icol,kk,ll) = faqgain_so4(imode)*dso4dt_aqrxn*cldfrc(icol,kk)
                  dqdt_aqh2so4(icol,kk,ll) = faqgain_so4(imode)*dso4dt_gasuptk*cldfrc(icol,kk)
                  dqdt_aq = dqdt_aqso4(icol,kk,ll) + dqdt_aqh2so4(icol,kk,ll)
                  dqdt_wr =  0.0_r8 ! don't have wet removal here
                  call update_tmr ( qcw(icol,kk,ll), dqdt_aq + dqdt_wr, dtime )
               endif
            enddo

            ! For gas species, tendency includes reactive uptake to cloud water
            ! that essentially transforms the gas to a different species.
            ! Need to multiply both these parts by cldfrc
            ! Currently it assumes no wet removal here

            ! h2so4 (g)         
            qin(icol,kk,id_h2so4) = qin(icol,kk,id_h2so4) - dso4dt_gasuptk * dtime * cldfrc(icol,kk)
! FORTRAN refactor: The order of multiplying cldfrc makes the following call
! failing BFB test, so this calculation is not refactored with new subroutine

            ! so2 -- the first order loss rate for so2 is frso2_c*clwlrat(i,k)
            dqdt_wr =  0.0_r8 ! don't have wet removal here
            dqdt_aq = -dso4dt_aqrxn*cldfrc(icol,kk)
            call update_tmr ( qin(icol,kk,id_so2), dqdt_aq + dqdt_wr, dtime )

            ! h2o2 -- the first order loss rate for h2o2 is frh2o2_c*clwlrat(i,k)
            dqdt_wr =  0.0_r8 ! don't have wet removal here
            dqdt_aq = -dso4dt_hprxn*cldfrc(icol,kk)
            call update_tmr ( qin(icol,kk,id_h2o2), dqdt_aq + dqdt_wr, dtime )

            ! for SO4 from H2O2/O3 budgets
            dqdt_aqhprxn(icol,kk) = dso4dt_hprxn*cldfrc(icol,kk)
            dqdt_aqo3rxn(icol,kk) = (dso4dt_aqrxn - dso4dt_hprxn)*cldfrc(icol,kk)

          endif cloud
       enddo col_loop
    enddo lev_loop

    !==============================================================
    ! ... Update the mixing ratios
    !==============================================================
    do imode = 1, ntot_amode
       call update_tmr_nonzero ( qcw, (lptr_so4_cw_amode(imode) - loffset) )
       call update_tmr_nonzero ( qcw, (lptr_nh4_cw_amode(imode) - loffset) )
    enddo
    call update_tmr_nonzero ( qin, id_so2 )


    ! diagnostics

    do imode = 1, ntot_amode
       mm = lptr_so4_cw_amode(imode)
       ll = mm - loffset
       if (ll > 0) then
          call calc_sfc_flux( dqdt_aqso4(:,:,ll)*adv_mass(ll)/mbar, pdel, sflx)
          call outfld( trim(cnst_name_cw(mm))//'AQSO4', sflx(:ncol), ncol, lchnk)

          call calc_sfc_flux( dqdt_aqh2so4(:,:,ll)*adv_mass(ll)/mbar, pdel, sflx)
          call outfld( trim(cnst_name_cw(mm))//'AQH2SO4', sflx(:ncol), ncol, lchnk)
       endif
    enddo

    call calc_sfc_flux( dqdt_aqhprxn*specmw_so4_amode/mbar, pdel, sflx)
    call outfld( 'AQSO4_H2O2', sflx(:ncol), ncol, lchnk)

    call calc_sfc_flux( dqdt_aqo3rxn*specmw_so4_amode/mbar, pdel, sflx)
    call outfld( 'AQSO4_O3', sflx(:ncol), ncol, lchnk)

#include "../yaml/sox_cldaero_mod/f90_yaml/sox_cldaero_update_end_yml.f90"
  end subroutine sox_cldaero_update

  !=============================================================================
  subroutine compute_aer_factor(tmr, loffset,   & ! in
                                faqgain_so4     ) ! out
    !-------------------------------------------------------------------------
    ! compute factors for partitioning aerosol mass gains among modes
    ! the factors are proportional to the activated particle MR for each
    ! mode, which is the MR of cloud drops "associated with" the mode
    ! thus we are assuming the cloud drop size is independent of the
    ! associated aerosol mode properties (i.e., drops associated with
    ! Aitken and coarse sea-salt particles are same size)
    ! qnum_c(n) = activated particle number MR for mode n (these are just
    ! used for partitioning among modes, so don't need to divide by cldfrc)
    !-------------------------------------------------------------------------
    real(r8), intent(in) :: tmr(:)   ! tracer mixing ratio [vmr]
    integer,  intent(in) :: loffset  ! # of tracers in the host model that are not part of MAM
    real(r8), intent(out) :: faqgain_so4(ntot_amode)   ! factor of TMR among modes [fraction]

    ! local variables
    integer  :: imode, ll               ! index
    real(r8) :: sumf                    ! total TMR for all modes [vmr]
    real(r8) :: qnum_c(ntot_amode)      ! tracer mixing ratio [vmr]
    real(r8), parameter :: small_value_10 = 1.e-10_r8
#include "../yaml/sox_cldaero_mod/f90_yaml/compute_aer_factor_beg_yml.f90"

    !-------------------------------------------------------------------------
    sumf  = 0.0_r8
    qnum_c(:)  = 0.0_r8
    faqgain_so4(:) = 0.0_r8

    do imode = 1, ntot_amode

       ll = numptrcw_amode(imode) - loffset
       if (ll > 0) qnum_c(imode) = max( 0.0_r8, tmr(ll) )

       ! force qnum_c(n) to be positive for n=modeptr_accum or n=1
       if (imode == modeptr_accum) qnum_c(imode) = max(small_value_10, qnum_c(imode))

       ! faqgain_so4(n) = fraction of total so4_c gain going to mode n
       ! these are proportional to the activated particle MR for each mode
       if (lptr_so4_cw_amode(imode) > 0) then
          faqgain_so4(imode) = qnum_c(imode)
          sumf = sumf + faqgain_so4(imode)
       endif

    enddo

    ! at this point (sumf <= 0.0) only when all the faqgain_so4 are zero
    if (sumf > 0.0_r8) faqgain_so4(:) = faqgain_so4(:)/sumf
#include "../yaml/sox_cldaero_mod/f90_yaml/compute_aer_factor_end_yml.f90"

  end subroutine compute_aer_factor

  !=============================================================================
  subroutine update_tmr ( tmr, dqdt, dtime )
    !-----------------------------------------------------------------------
    ! update tracer mixing ratio by adding tendencies
    !-----------------------------------------------------------------------
    real(r8), intent(inout) :: tmr   ! tracer mixing ratio [vmr]
    real(r8), intent(in)    :: dqdt  ! tmr tendency [vmr/s]
    real(r8), intent(in)    :: dtime ! time step [s]
#include "../yaml/sox_cldaero_mod/f90_yaml/update_tmr_beg_yml.f90"

    tmr = tmr + dqdt * dtime

#include "../yaml/sox_cldaero_mod/f90_yaml/update_tmr_end_yml.f90"
  end subroutine update_tmr
  !=============================================================================
  subroutine update_tmr_nonzero ( tmr, idx )
    !-----------------------------------------------------------------------
    ! basically it just makes sure the value is greater than zero
    !-----------------------------------------------------------------------
    real(r8), intent(inout) :: tmr(:,:,:) ! tracer mixing ratio [vmr]
    integer,  intent(in)    :: idx        ! index for the third dimension of vmr
    
    real(r8), parameter :: small_value_20 = 1.e-20_r8

    if (idx>0) then
        tmr(:,:,idx) = max(tmr(:,:,idx), small_value_20)
    endif

  end subroutine update_tmr_nonzero

  !=============================================================================
  subroutine calc_sfc_flux(layer_tend, pdel, sflx)
    !-----------------------------------------------------------------------
    ! calculate surface fluxes of wet deposition from vertical integration of tendencies 
    !-----------------------------------------------------------------------
    real(r8), intent(in) :: pdel(:,:)      ! pressure difference between two layers [Pa]
    real(r8), intent(in) :: layer_tend(:,:)! physical tendencies in each layer [kg/kg/s]
    real(r8), intent(out):: sflx(:)        ! integrated surface fluxes [kg/m2/s]

    integer :: kk
#include "../yaml/sox_cldaero_mod/f90_yaml/calc_sfc_flux_beg_yml.f90"

     sflx(:)=0.0_r8
     do kk=1,pver
        sflx(:) = sflx(:) + layer_tend(:,kk)*pdel(:,kk)/gravit
     enddo

#include "../yaml/sox_cldaero_mod/f90_yaml/calc_sfc_flux_end_yml.f90"
  end subroutine calc_sfc_flux

  !=============================================================================
  subroutine sox_cldaero_destroy_obj( conc_obj )
  !----------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------
    use cldaero_mod, only : cldaero_deallocate

    type(cldaero_conc_t), pointer :: conc_obj

    call cldaero_deallocate( conc_obj )

  end subroutine sox_cldaero_destroy_obj

end module sox_cldaero_mod
