module mo_gas_phase_chemdr
#include "../yaml/common_files/common_uses.ymlf90"
  use shr_kind_mod,     only : r8 => shr_kind_r8
  use shr_const_mod,    only : pi => shr_const_pi
  use constituents,     only : pcnst
  use cam_history,      only : fieldname_len
  use chem_mods,        only : phtcnt, rxntot, gas_pcnst
  use chem_mods,        only : rxt_tag_cnt, rxt_tag_lst, rxt_tag_map, extcnt
  use dust_model,       only : dust_names, ndust => dust_nbin
  use ppgrid,           only : pcols, pver
  use phys_control,     only : phys_getopts

  implicit none
  save

  private
  public :: gas_phase_chemdr, gas_phase_chemdr_inti 
  public :: map2chm

  integer :: map2chm(pcnst) = 0           ! index map to/from chemistry/constituents list

  integer :: synoz_ndx
  integer :: o3_ndx
  integer :: ndx_h2so4
  integer :: inv_ndx_cnst_o3

  character(len=fieldname_len),dimension(rxntot-phtcnt) :: rxn_names
  character(len=fieldname_len),dimension(phtcnt)        :: pht_names
  character(len=fieldname_len),dimension(rxt_tag_cnt)   :: tag_names
  character(len=fieldname_len),dimension(extcnt)        :: extfrc_name
  logical :: convproc_do_aer 

contains

  subroutine gas_phase_chemdr_inti(chem_name)

    use mo_chem_utls,      only : get_spc_ndx, get_extfrc_ndx, get_inv_ndx
    use cam_history,       only : addfld, add_default, horiz_only
    use cam_abortutils,    only : endrun
    use mo_chm_diags,      only : chm_diags_inti
    use rate_diags,        only : rate_diags_init

    implicit none

    character(len=*), intent(in) :: chem_name

    character(len=3)  :: string
    integer           :: nn, mm
    logical           :: history_aerosol      ! Output the MAM aerosol tendencies

    !-----------------------------------------------------------------------

    call phys_getopts( history_aerosol_out = history_aerosol, &
         convproc_do_aer_out = convproc_do_aer ) 

    ndx_h2so4 = get_spc_ndx('H2SO4')
    o3_ndx    = get_spc_ndx('O3')
    synoz_ndx = get_extfrc_ndx( 'SYNOZ' )

    do mm = 1,extcnt
       WRITE(UNIT=string, FMT='(I2.2)') mm
       extfrc_name(mm) = 'extfrc_'// trim(string)
       call addfld( extfrc_name(mm), (/ 'lev' /), 'I', ' ', 'ext frcing' )
    enddo

    do nn = 1,rxt_tag_cnt
       tag_names(nn) = trim(rxt_tag_lst(nn))
       if (nn<=phtcnt) then
          call addfld( tag_names(nn), (/ 'lev' /), 'I', '/s', 'photolysis rate' )
       else
          call addfld( tag_names(nn), (/ 'lev' /), 'I', '/cm3/s', 'reaction rate' )
       endif
    enddo

    do nn = 1,phtcnt
       WRITE(UNIT=string, FMT='(I3.3)') nn
       pht_names(nn) = 'J_' // trim(string)
       call addfld( pht_names(nN), (/ 'lev' /), 'I', '/s', 'photolysis rate' )
    enddo

    do nn = 1,rxntot-phtcnt
       WRITE(UNIT=string, FMT='(I3.3)') nn
       rxn_names(nn) = 'R_' // trim(string)
       call addfld( rxn_names(nn), (/ 'lev' /), 'I', '/cm3/s', 'reaction rate' )
    enddo

    call addfld( 'DTCBS',  horiz_only, 'I',   ' ','photolysis diagnostic black carbon OD' )
    call addfld( 'DTOCS',  horiz_only, 'I',   ' ','photolysis diagnostic organic carbon OD' )
    call addfld( 'DTSO4',  horiz_only, 'I',   ' ','photolysis diagnostic SO4 OD' )
    call addfld( 'DTSOA',  horiz_only, 'I',   ' ','photolysis diagnostic SOA OD' )
    call addfld( 'DTANT',  horiz_only, 'I',   ' ','photolysis diagnostic NH4SO4 OD' )
    call addfld( 'DTSAL',  horiz_only, 'I',   ' ','photolysis diagnostic salt OD' )
    call addfld( 'DTDUST',  horiz_only, 'I',  ' ','photolysis diagnostic dust OD' )
    call addfld( 'DTTOTAL',  horiz_only, 'I', ' ','photolysis diagnostic total aerosol OD' )   
    call addfld( 'FRACDAY',  horiz_only, 'I', ' ','photolysis diagnostic fraction of day' )

    call addfld( 'QDSAD', (/ 'lev' /), 'I', '/s', 'water vapor sad delta' )
    call addfld( 'SAD', (/ 'lev' /), 'I', 'cm2/cm3', 'sulfate aerosol SAD' )
    call addfld( 'SAD_SULFC', (/ 'lev' /), 'I', 'cm2/cm3', 'chemical sulfate aerosol SAD' )
    call addfld( 'SAD_SAGE', (/ 'lev' /), 'I', 'cm2/cm3', 'SAGE sulfate aerosol SAD' )
    call addfld( 'SAD_LNAT', (/ 'lev' /), 'I', 'cm2/cm3', 'large-mode NAT aerosol SAD' )
    call addfld( 'SAD_ICE', (/ 'lev' /), 'I', 'cm2/cm3', 'water-ice aerosol SAD' )
    call addfld( 'RAD_SULFC', (/ 'lev' /), 'I', 'cm', 'chemical sad sulfate' )
    call addfld( 'RAD_LNAT', (/ 'lev' /), 'I', 'cm', 'large nat radius' )
    call addfld( 'RAD_ICE', (/ 'lev' /), 'I', 'cm', 'sad ice' )
    call addfld( 'SAD_TROP', (/ 'lev' /), 'I', 'cm2/cm3', 'tropospheric aerosol SAD' )
    call addfld( 'HNO3_STS', (/ 'lev' /), 'I', 'mol/mol', 'STS condensed HNO3' )
    call addfld( 'HNO3_NAT', (/ 'lev' /), 'I', 'mol/mol', 'NAT condensed HNO3' )
    call addfld( 'QDSETT', (/ 'lev' /), 'I', '/s', 'water vapor settling delta' )
    call addfld( 'QDCHEM', (/ 'lev' /), 'I', '/s', 'water vapor chemistry delta' )
    call addfld( 'HNO3_GAS', (/ 'lev' /), 'I', 'mol/mol', 'gas-phase hno3' )
    call addfld( 'H2O_GAS', (/ 'lev' /), 'I', 'mol/mol', 'gas-phase h2o' )
    call addfld( 'SZA', horiz_only, 'I', 'degrees', 'solar zenith angle' )

    call chm_diags_inti()
    call rate_diags_init()

    !-----------------------------------------------------------------------
    ! get fixed oxidant (troposphere) index for Linoz_MAM
    !-----------------------------------------------------------------------

    inv_ndx_cnst_o3 = get_inv_ndx( 'cnst_O3' ) ! prescribed O3 oxidant field

    if ( chem_name == 'linoz_mam3'.or.chem_name == 'linoz_mam4_resus'.or.chem_name == 'linoz_mam4_resus_mom' &
         .or.chem_name == 'linoz_mam4_resus_soag'.or.chem_name == 'linoz_mam4_resus_mom_soag') then
       if ( inv_ndx_cnst_o3 < 1 ) then
          call endrun('ERROR: chem_name = '//trim(chem_name)//&
               ' requies cnst_O3 fixed oxidant field. Use cnst_O3:O3 in namelist tracer_cnst_specifier')
       endif
    endif

  end subroutine gas_phase_chemdr_inti


  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine gas_phase_chemdr(lchnk, ncol, imozart, state_q, &
       phis, zm, zi, calday, &
       tfld, pmid, pdel, pdeldry, pint,  &
       cldw, troplev, &
       ncldwtr, ufld, vfld,  &
       prain, cldfr, cmfdqr, nevapr, &
       delt, ps, linoz_o3_clim, linoz_t_clim, &
       linoz_o3col_clim, linoz_PmL_clim, linoz_dPmL_dO3, &
       linoz_dPmL_dT, linoz_dPmL_dO3col, linoz_cariolle_psc, & !in
       fsds, ts, asdir, &
       precc, precl, snowhland, pblh, &
       drydepflx, cflx, qtend, pbuf, qqcw, &
       dgnum, dgnumwet, wetdens         ) ! inout
    !-----------------------------------------------------------------------
    !     ... Chem_solver advances the volumetric mixing ratio
    !         forward one time step via a combination of explicit,
    !         ebi, hov, fully implicit, and/or rodas algorithms.
    !-----------------------------------------------------------------------

    use chem_mods,         only : nabscol, nfs, indexm
    use physconst,         only : rga
    use mo_photo,          only : set_ub_col, setcol, table_photo
    use mo_imp_sol,        only : imp_sol
    use mo_setrxt,         only : setrxt
    use mo_adjrxt,         only : adjrxt
    use mo_usrrxt,         only : usrrxt
    use mo_setinv,         only : setinv
    use mo_negtrc,         only : negtrc
    use mo_setext,         only : setext
    use mo_sethet,         only : sethet
    use mo_drydep,         only : drydep_xactive
    use phys_grid,         only : get_rlat_all_p, get_rlon_all_p, get_lat_all_p, get_lon_all_p
    use mo_mean_mass,      only : set_mean_mass
    use cam_history,       only : outfld
    use wv_saturation,     only : qsat
    use time_manager,      only : get_ref_date
    use shr_orb_mod,       only : shr_orb_decl
    use cam_control_mod,   only : lambm0, eccen, mvelpp, obliqr
    use mo_chm_diags,      only : chm_diags, het_diags
    use perf_mod,          only : t_startf, t_stopf
    use physics_buffer,    only : physics_buffer_desc
    use infnan,            only : nan, assignment(=)
    use rate_diags,        only : rate_diags_calc
    use mo_mass_xforms,    only : mmr2vmr, vmr2mmr, h2o_to_vmr
    use orbit,             only : zenith
    use mam_support,       only : min_max_bound
    !
    ! LINOZ
    !
    use lin_strat_chem,    only : lin_strat_chem_solve, lin_strat_sfcsink
    !
    ! for aqueous chemistry and aerosol growth
    !
    use aero_model,        only : aero_model_gasaerexch
    use mam_support,       only : ptr2d_t

    implicit none

    !-----------------------------------------------------------------------
    !        ... Dummy arguments
    !-----------------------------------------------------------------------
    integer,        intent(in)    :: lchnk                          ! chunk index
    integer,        intent(in)    :: ncol                           ! number columns in chunk
    integer,        intent(in)    :: imozart                        ! gas phase start index in state_q
    real(r8),       intent(in)    :: delt                           ! timestep (s)
    real(r8),       intent(in)    :: calday                         ! day of year
    real(r8),       intent(in)    :: ps(pcols)                      ! surface pressure
    !pointers to read LINOZ data
    real(r8),       intent(in)    :: linoz_o3_clim(pcols,pver)      ! ozone (climatology) [vmr]              
    real(r8),       intent(in)    :: linoz_t_clim(pcols,pver)       ! temperature (climatology) [K]
    real(r8),       intent(in)    :: linoz_o3col_clim(pcols,pver)   ! Column O3 above box (climatology) [Dobson Units or DU]
    real(r8),       intent(in)    :: linoz_PmL_clim(pcols,pver)     ! P minus L (climatology) [vmr/s]
    real(r8),       intent(in)    :: linoz_dPmL_dO3(pcols,pver)     ! Sensitivity of P minus L to O3 [1/s]
    real(r8),       intent(in)    :: linoz_dPmL_dT(pcols,pver)      ! Sensitivity of P minus L to T [K]
    real(r8),       intent(in)    :: linoz_dPmL_dO3col(pcols,pver)  ! Sensitivity of P minus L to overhead O3 column [vmr/DU]
    real(r8),       intent(in)    :: linoz_cariolle_psc(pcols,pver) ! Cariolle parameter for PSC loss of ozone [1/s]
    real(r8),       intent(in)    :: phis(pcols)                    ! surface geopotential
    real(r8),       intent(in)    :: tfld(pcols,pver)               ! midpoint temperature (K)
    real(r8),       intent(in)    :: pmid(pcols,pver)               ! midpoint pressures (Pa)
    real(r8),       intent(in)    :: pdel(pcols,pver)               ! pressure delta about midpoints (Pa)
    real(r8),       intent(in)    :: pdeldry(pcols,pver)            ! dry pressure delta about midpoints (Pa)
    real(r8),       intent(in)    :: ufld(pcols,pver)               ! zonal velocity (m/s)
    real(r8),       intent(in)    :: vfld(pcols,pver)               ! meridional velocity (m/s)
    real(r8),       intent(in)    :: cldw(pcols,pver)               ! cloud water (kg/kg)
    real(r8),       intent(in)    :: ncldwtr(pcols,pver)            ! droplet number concentration (#/kg)
    real(r8),       intent(in)    :: zm(pcols,pver)                 ! midpoint geopotential height above the surface (m)
    real(r8),       intent(in)    :: zi(pcols,pver+1)               ! interface geopotential height above the surface (m)
    real(r8),       intent(in)    :: pint(pcols,pver+1)             ! interface pressures (Pa)
    real(r8),       intent(in)    :: state_q(pcols,pver,pcnst)      ! species concentrations (kg/kg)
    real(r8),       intent(in)    :: fsds(pcols)                    ! longwave down at sfc
    real(r8),       intent(in)    :: asdir(pcols)                   ! albedo: shortwave, direct
    real(r8),       intent(in)    :: ts(pcols)                      ! sfc temp (merged w/ocean if coupled)
    real(r8),       intent(in)    :: precc(pcols)                   !
    real(r8),       intent(in)    :: precl(pcols)                   !
    real(r8),       intent(in)    :: snowhland(pcols)               !
    real(r8),       intent(in)    :: pblh(:)                        ! pbl height [m]
    real(r8),       intent(in)    :: prain(:,:)
    real(r8),       intent(in)    :: nevapr(:,:)
    real(r8),       intent(in)    :: cmfdqr(:,:)
    real(r8),       intent(in)    :: cldfr(:,:)
    integer,        intent(in)    ::  troplev(pcols)

    real(r8),       intent(inout) :: qtend(pcols,pver,pcnst)        ! species tendencies (kg/kg/s)
    real(r8),       intent(inout) :: cflx(pcols,pcnst)              ! constituent surface flux (kg/m^2/s)
    real(r8),       intent(out)   :: drydepflx(pcols,pcnst)         ! dry deposition flux (kg/m^2/s)
    type(ptr2d_t), intent(inout), target :: qqcw(:)                ! Cloud borne aerosols mixing ratios [kg/kg or 1/kg]
    real(r8),       intent(inout), target    :: dgnum(:,:,:)                   ! aerosol diameter [m]
    real(r8),       intent(inout), target    :: dgnumwet(:,:,:)                ! aerosol wet diameter [m]
    real(r8),       intent(inout), target    :: wetdens(:,:,:)                 ! aerosol wet density [kg/m3] 


    type(physics_buffer_desc), pointer :: pbuf(:)

    !-----------------------------------------------------------------------
    !     	... Local variables
    !-----------------------------------------------------------------------
    real(r8), parameter :: m2km  = 1.e-3_r8
    real(r8), parameter :: rad2deg = 180._r8/pi            ! radians to degrees conversion factor

    integer      ::  ii, kk, mm, nn
    real(r8)     ::  delt_inverse
    real(r8)     ::  esfact
    integer      ::  latndx(pcols)                         ! chunk lat indicies
    integer      ::  lonndx(pcols)                         ! chunk lon indicies
    real(r8)     ::  invariants(ncol,pver,nfs)
    real(r8)     ::  col_dens(ncol,pver,max(1,nabscol))    ! column densities (molecules/cm^2)
    real(r8)     ::  col_delta(ncol,0:pver,max(1,nabscol)) ! layer column densities (molecules/cm^2)
    real(r8)     ::  extfrc(ncol,pver,max(1,extcnt))
    real(r8)     ::  vmr(ncol,pver,gas_pcnst)              ! xported species (vmr)
    real(r8)     ::  reaction_rates(ncol,pver,max(1,rxntot))      ! reaction rates
    real(r8)     ::  depvel(ncol,gas_pcnst)                ! dry deposition velocity (cm/s)
    real(r8)     ::  het_rates(ncol,pver,max(1,gas_pcnst)) ! washout rate (1/s)

    real(r8)     ::  h2ovmr(ncol,pver)                     ! water vapor volume mixing ratio
    real(r8)     ::  mbar(ncol,pver)                       ! mean wet atmospheric mass ( amu )
    real(r8)     ::  zmid(ncol,pver)                       ! midpoint geopotential in km
    real(r8)     ::  cwat(ncol,pver)                       ! cloud water mass mixing ratio (kg/kg)
    real(r8)     ::  zintr(ncol,pver+1)                    ! interface geopotential in km realitive to surf
    real(r8)     ::  zint(ncol,pver+1)                     ! interface geopotential in km
    real(r8)     ::  zen_angle(ncol)                       ! solar zenith angles
    real(r8)     ::  zsurf(ncol)                           ! surface height (m)
    real(r8)     ::  rlats(ncol), rlons(ncol)              ! chunk latitudes and longitudes (radians)
    real(r8)     ::  sza(ncol)                             ! solar zenith angles (degrees)
    real(r8)     ::  relhum(ncol,pver)                     ! relative humidity
    real(r8)     ::  satv(ncol,pver)                       ! wrk array for relative humidity
    real(r8)     ::  satq(ncol,pver)                       ! wrk array for relative humidity

    integer      ::  ltrop_sol(pcols)                 ! tropopause vertical index used in chem solvers
    real(r8)     ::  sad_total(pcols,pver)            ! total trop. SAD (cm2/cm3)

    real(r8) :: tvs(pcols)
    integer  :: ncdate,yr,mon,day,sec
    real(r8) :: wind_speed(pcols)        ! surface wind speed (m/s)
    real(r8) :: soilw(pcols)
    real(r8) :: prect(pcols)
    real(r8) :: sflx(pcols,gas_pcnst)
    real(r8) :: mmr(pcols,pver,gas_pcnst)      ! chem working concentrations (kg/kg)
    real(r8) :: mmr_new(pcols,pver,gas_pcnst)  ! chem working concentrations (kg/kg)
    real(r8) :: mmr_tend(pcols,pver,gas_pcnst) ! chemistry species tendencies (kg/kg/s)
    real(r8) :: qh2o(pcols,pver)               ! specific humidity (kg/kg)
    real(r8) :: delta

    ! for aerosol formation....  
    real(r8) :: del_h2so4_gasprod(ncol,pver)
    real(r8) :: vmr0(ncol,pver,gas_pcnst)
#include "../yaml/mo_gas_phase_chemdr/f90_yaml/gas_phase_chemdr_beg_yml.f90"
    call t_startf('chemdr_init')

    ! initialize to NaN to hopefully catch user defined rxts that go unset
    reaction_rates(:,:,:) = nan

    delt_inverse = 1._r8 / delt ! inverse of time step
    !-----------------------------------------------------------------------      
    !        ... Get chunck latitudes and longitudes
    !-----------------------------------------------------------------------      
    call get_rlat_all_p( lchnk, ncol, rlats )
    call get_rlon_all_p( lchnk, ncol, rlons )

    !-----------------------------------------------------------------------      
    !        ... Calculate cosine of zenith angle
    !            then cast back to angle (radians)
    !-----------------------------------------------------------------------      
    call zenith( calday, rlats, rlons, & !in
         zen_angle, &!out
         ncol ) !in
    zen_angle(:) = acos( zen_angle(:) )

    sza(:) = zen_angle(:) * rad2deg
    call outfld( 'SZA',   sza,    ncol, lchnk )

    !-----------------------------------------------------------------------      
    !        ... Xform geopotential height from m to km 
    !            and pressure from Pa to mb
    !-----------------------------------------------------------------------      
    zsurf(:ncol) = rga * phis(:ncol)
    do kk = 1, pver
       zintr(:ncol,kk) = m2km * zi(:ncol,kk)
       zmid(:ncol,kk) = m2km * (zm(:ncol,kk) + zsurf(:ncol))
       zint(:ncol,kk) = m2km * (zi(:ncol,kk) + zsurf(:ncol))
    enddo
    zint(:ncol,pver+1) = m2km * (zi(:ncol,pver+1) + zsurf(:ncol))
    zintr(:ncol,pver+1)= m2km *  zi(:ncol,pver+1)

    !-----------------------------------------------------------------------      
    !        ... map incoming concentrations to working array
    !-----------------------------------------------------------------------      
    do mm = 1,pcnst
       nn = map2chm(mm)
       if( nn > 0 ) then
          mmr(:ncol,:,nn) = state_q(:ncol,:,mm)
       endif
    enddo

    !-----------------------------------------------------------------------      
    !        ... Set atmosphere mean mass
    !-----------------------------------------------------------------------      
    call set_mean_mass( ncol, & !in 
         mbar ) !out

    !-----------------------------------------------------------------------      
    !        ... Xform from mmr to vmr
    !-----------------------------------------------------------------------      
    call mmr2vmr( mmr, & !in
         vmr, & !in-out
         mbar, ncol ) !in


    qh2o(:ncol,:) = state_q(:ncol,:,1)
    !-----------------------------------------------------------------------      
    !        ... Xform water vapor from mmr to vmr and set upper bndy values
    !-----------------------------------------------------------------------      
    call h2o_to_vmr( state_q(:,:,1), & !in
         h2ovmr, &!in-out
         mbar, ncol ) !in

    !-----------------------------------------------------------------------      
    !        ... Set the "invariants"
    !-----------------------------------------------------------------------  
    call setinv( invariants, &! out
         tfld, h2ovmr, vmr, pmid, ncol, lchnk, pbuf ) !in

    !-----------------------------------------------------------------------      
    !        ... Set the column densities at the upper boundary
    !-----------------------------------------------------------------------      
    call set_ub_col( col_delta, & ! out
         vmr, invariants, pdel, ncol, lchnk) ! in

    !-----------------------------------------------------------------------      
    !       ...  Set rates for "tabular" and user specified reactions
    !-----------------------------------------------------------------------      
    call setrxt( reaction_rates, & ! inout
         tfld, ncol )  ! in

    !-----------------------------------------------------------------
    !	... compute the relative humidity
    !-----------------------------------------------------------------
    call qsat(tfld(:ncol,:), pmid(:ncol,:), satv, satq)

    do kk = 1,pver
       relhum(:,kk) = .622_r8 * h2ovmr(:,kk) / satq(:,kk)
       relhum(:,kk) = min_max_bound(0._r8, 1._r8,relhum(:,kk),size(relhum(:,kk)))
    enddo

    cwat(:ncol,:pver) = cldw(:ncol,:pver)

    call usrrxt( reaction_rates, &  ! inout
         tfld, invariants, invariants(:,:,indexm), ncol ) ! in

    call outfld( 'SAD_TROP', sad_total(:ncol,:), ncol, lchnk )


    do ii = phtcnt+1,rxntot
       call outfld( rxn_names(ii-phtcnt), reaction_rates(:,:,ii), ncol, lchnk )
    enddo

    call adjrxt( reaction_rates, & ! inout
         invariants, invariants(1,1,indexm), ncol )  ! in

    !-----------------------------------------------------------------------
    !        ... Compute the photolysis rates at time = t(n+1)
    !-----------------------------------------------------------------------      
    !-----------------------------------------------------------------------      
    !     	... Set the column densities
    !-----------------------------------------------------------------------      
    call setcol(  col_delta, & ! in
         col_dens ) ! out

    !-----------------------------------------------------------------------      
    !     	... Calculate the photodissociation rates
    !-----------------------------------------------------------------------      

    esfact = 1._r8
    call shr_orb_decl( calday, eccen, mvelpp, lambm0, obliqr  , & !in
         delta, esfact ) !out


    !-----------------------------------------------------------------
    !	... lookup the photolysis rates from table
    !-----------------------------------------------------------------
    ! FORTRAN refactor notes: it looks that reaction_rates is reset in this
    ! subroutine, and does not depend on how it was calculated before
    call table_photo( reaction_rates, & ! out
         pmid, pdel, tfld, & ! in
         col_dens, zen_angle, asdir, cwat, cldfr, & ! in
         esfact,  ncol ) ! in

    do ii = 1,phtcnt
       call outfld( pht_names(ii), reaction_rates(:ncol,:,ii), ncol, lchnk )
       call outfld( tag_names(ii), reaction_rates(:ncol,:,rxt_tag_map(ii)), ncol, lchnk )
    enddo

    !-----------------------------------------------------------------------
    !        ... Compute the extraneous frcing at time = t(n+1)
    !-----------------------------------------------------------------------      
    call setext( extfrc,            & ! out
         lchnk, ncol, zintr ) ! in

    do mm = 1,extcnt
       if( mm /= synoz_ndx ) then
          do kk = 1,pver
             extfrc(:ncol,kk,mm) = extfrc(:ncol,kk,mm) / invariants(:ncol,kk,indexm)
          enddo
       endif
       call outfld( extfrc_name(mm), extfrc(:ncol,:,mm), ncol, lchnk )
    enddo

    !-----------------------------------------------------------------------
    !        ... Form the washout rates
    !-----------------------------------------------------------------------      
    call sethet( het_rates, pmid, zmid, phis, tfld, &
         cmfdqr, prain, nevapr, delt, invariants(:,:,indexm), &
         vmr, ncol, lchnk )

    do ii = phtcnt+1,rxt_tag_cnt
       call outfld( tag_names(ii), reaction_rates(:ncol,:,rxt_tag_map(ii)), ncol, lchnk )
    enddo

    ltrop_sol(:ncol) = 0 ! apply solver to all levels

    ! save h2so4 before gas phase chem (for later new particle nucleation)
    del_h2so4_gasprod(1:ncol,:) = vmr(1:ncol,:,ndx_h2so4)

    vmr0(:ncol,:,:) = vmr(:ncol,:,:) ! mixing ratios before chemistry changes

    call t_stopf('chemdr_init')

    !=======================================================================
    !        ... Call the class solution algorithms
    !=======================================================================

    !-----------------------------------------------------------------------
    !	... Solve for "Implicit" species
    !-----------------------------------------------------------------------
    !
    call t_startf('imp_sol')
    call imp_sol( vmr, reaction_rates, het_rates, extfrc, delt, &
         invariants(1,1,indexm), ncol, lchnk, ltrop_sol(:ncol) )
    call t_stopf('imp_sol')

    if(convproc_do_aer) then 
       call vmr2mmr( vmr, mmr_new, mbar, ncol )  !RCE
       mmr_new(:ncol,:,:) = 0.5_r8*( mmr(:ncol,:,:)+mmr_new(:ncol,:,:) )  !RCE
       !RCE - mmr_new = average of mmr values before and after imp_sol
       call het_diags( het_rates(:ncol,:,:), mmr_new(:ncol,:,:), pdel(:ncol,:), lchnk, ncol )  !RCE
    endif

    ! save h2so4 change by gas phase chem (for later new particle nucleation)
    if (ndx_h2so4 > 0) then
       del_h2so4_gasprod(1:ncol,:) = vmr(1:ncol,:,ndx_h2so4) - del_h2so4_gasprod(1:ncol,:)
    endif

    !
    ! Aerosol processes ...
    !
    !-----------------------------------------------------------------------      
    !        ... Get chunck latitudes and longitudes
    !-----------------------------------------------------------------------      
    call get_lat_all_p( lchnk, ncol, latndx )
    call get_lon_all_p( lchnk, ncol, lonndx )

    call t_startf('aero_model_gasaerexch')
    call aero_model_gasaerexch( imozart-1, ncol, lchnk, delt, latndx, lonndx, & !in
         tfld, pmid, pdel, mbar, zm,  qh2o, cwat,      & !in
         cldfr, ncldwtr, invariants(:,:,indexm), vmr0, & !in
         pblh,                                         & !in
         vmr, qqcw, dgnum, dgnumwet, wetdens         ) ! inout
    call t_stopf('aero_model_gasaerexch')

    !
    ! LINOZ
    !
    call lin_strat_chem_solve( ncol, lchnk, col_dens(:,:,1), tfld, zen_angle, pmid, delt, rlats, troplev, & !in
         linoz_o3_clim, linoz_t_clim, linoz_o3col_clim, linoz_PmL_clim, linoz_dPmL_dO3, linoz_dPmL_dT, & !in
         linoz_dPmL_dO3col, linoz_cariolle_psc, & !in
         vmr(:,:,o3_ndx) ) !in-out

    call   lin_strat_sfcsink (ncol, lchnk, delt, pdel(:ncol,:), & !in
         vmr(:,:,o3_ndx)) !in-outs

    !-----------------------------------------------------------------------      
    !         ... Check for negative values and reset to zero
    !-----------------------------------------------------------------------      
    call negtrc( vmr, ncol )

    !-----------------------------------------------------------------------      
    !         ... Xform from vmr to mmr
    !-----------------------------------------------------------------------      
    call vmr2mmr( vmr, mmr_tend, mbar, ncol )

    !-----------------------------------------------------------------------      
    !         ... Form the tendencies
    !----------------------------------------------------------------------- 
    do mm = 1,gas_pcnst 
       mmr_new(:ncol,:,mm) = mmr_tend(:ncol,:,mm)
       mmr_tend(:ncol,:,mm) = (mmr_tend(:ncol,:,mm) - mmr(:ncol,:,mm))*delt_inverse
    enddo

    do mm = 1,pcnst
       nn = map2chm(mm)
       if( nn > 0 ) then
          qtend(:ncol,:,mm) = qtend(:ncol,:,mm) + mmr_tend(:ncol,:,nn) 
       endif
    enddo

    tvs(:ncol) = tfld(:ncol,pver) * (1._r8 + qh2o(:ncol,pver))

    sflx(:,:) = 0._r8
    call get_ref_date(yr, mon, day, sec)
    ncdate = yr*10000 + mon*100 + day
    wind_speed(:ncol) = sqrt( ufld(:ncol,pver)*ufld(:ncol,pver) + vfld(:ncol,pver)*vfld(:ncol,pver) )
    prect(:ncol) = precc(:ncol) + precl(:ncol)

    call t_startf('drydep')

    !drydep_method == DD_XATM 
    call drydep_xactive( lchnk, ncol, latndx, &                                         ! in
         ncdate, ts, tfld(:,pver), tvs, ps, pmid(:,pver), &                             ! in
         qh2o(:,pver), wind_speed, prect, snowhland, fsds, mmr, &                       ! in
         depvel, &                                                                      ! out
         sflx)                                                                          ! inout
    call t_stopf('drydep')

    drydepflx(:,:) = 0._r8
    do mm = 1,pcnst
       nn = map2chm( mm )
       if ( nn > 0 ) then
          cflx(:ncol,mm)      = cflx(:ncol,mm) - sflx(:ncol,nn)
          drydepflx(:ncol,mm) = sflx(:ncol,nn)
       endif
    enddo

    call t_startf('chemdr_diags')
    call chm_diags( lchnk, ncol, vmr(:ncol,:,:), mmr_new(:ncol,:,:), &       !intent-in
         depvel(:ncol,:),  sflx(:ncol,:), &                       !intent-in
         mmr_tend(:ncol,:,:), pdel(:ncol,:), pdeldry(:ncol,:), &  !intent-in
         qqcw, troplev(:ncol)  )                                  !intent-in

    call rate_diags_calc( reaction_rates(:,:,:), &! inout
         vmr(:,:,:), invariants(:,:,indexm), ncol, lchnk ) !in
    call t_stopf('chemdr_diags')
#include "../yaml/mo_gas_phase_chemdr/f90_yaml/gas_phase_chemdr_end_yml.f90"
  end subroutine gas_phase_chemdr

end module mo_gas_phase_chemdr
