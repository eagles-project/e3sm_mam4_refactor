module mo_chm_diags
#include "../yaml/common_files/common_uses.ymlf90"

  use shr_kind_mod, only : r8 => shr_kind_r8
  use chem_mods,    only : gas_pcnst
  use mo_tracname,  only : solsym
  use chem_mods,    only : gas_pcnst, adv_mass
  use ppgrid,       only : pcols, pver
  use mo_constants, only : pi, rgrav, rearth, avogadro
  use mo_chem_utls, only : get_rxt_ndx, get_spc_ndx
  use cam_history,  only : fieldname_len

  private

  public :: chm_diags_inti
  public :: chm_diags
  public :: het_diags

  integer :: id_o3
  integer :: sox_species(3)
  integer :: aer_species(gas_pcnst)

  character(len=fieldname_len) :: dtchem_name(gas_pcnst)
  character(len=fieldname_len) :: depvel_name(gas_pcnst)
  character(len=fieldname_len) :: depflx_name(gas_pcnst)
  character(len=fieldname_len) :: wetdep_name(gas_pcnst)
  character(len=fieldname_len) :: wtrate_name(gas_pcnst)

  real(r8), parameter :: S_molwgt = 32.066_r8

  ! constants for converting O3 mixing ratio to DU
  real(r8), parameter :: DUfac = 2.687e20_r8   ! 1 DU in molecules per m^2

contains

!========================================================================
  subroutine chm_diags_inti
    !--------------------------------------------------------------------
    !	... initialize utility routine
    !--------------------------------------------------------------------

    use cam_history,  only : addfld, horiz_only, add_default
    use constituents, only : cnst_get_ind, cnst_longname
    use phys_control, only : phys_getopts

    implicit none

    integer :: i, j, k, m, n
    character(len=16) :: spc_name, attr
    character(len=2)  :: unit_basename  ! Units 'kg' or '1' 

    integer :: id_nh4no3
    integer :: id_so2, id_so4, id_h2so4
    integer :: id_dst01, id_dst02, id_dst03, id_dst04, id_sslt01, id_sslt02, id_sslt03, id_sslt04
    integer :: id_soa,  id_oc1, id_oc2, id_cb1, id_cb2
    integer :: id_soam,id_soai,id_soat,id_soab,id_soax

    logical :: history_aerosol      ! Output the MAM aerosol tendencies
    logical :: history_amwg         ! output the variables used by the AMWG diag package
    logical :: history_verbose      ! produce verbose history output
    integer :: bulkaero_species(20)

    !-----------------------------------------------------------------------

    call phys_getopts( history_aerosol_out = history_aerosol, &
                       history_amwg_out    = history_amwg,  &
                       history_verbose_out = history_verbose)

    id_o3      = get_spc_ndx( 'O3' )

    id_so2     = get_spc_ndx( 'SO2' )
    id_so4     = get_spc_ndx( 'SO4' )
    id_h2so4   = get_spc_ndx( 'H2SO4' )

    id_nh4no3  = get_spc_ndx( 'NH4NO3' )
    id_dst01   = get_spc_ndx( 'DST01' )
    id_dst02   = get_spc_ndx( 'DST02' )
    id_dst03   = get_spc_ndx( 'DST03' )
    id_dst04   = get_spc_ndx( 'DST04' )
    id_sslt01  = get_spc_ndx( 'SSLT01' )
    id_sslt02  = get_spc_ndx( 'SSLT02' )
    id_sslt03  = get_spc_ndx( 'SSLT03' )
    id_sslt04  = get_spc_ndx( 'SSLT04' )
    id_soa     = get_spc_ndx( 'SOA' )
    id_oc1     = get_spc_ndx( 'OC1' )
    id_oc2     = get_spc_ndx( 'OC2' )
    id_cb1     = get_spc_ndx( 'CB1' )
    id_cb2     = get_spc_ndx( 'CB2' )


    id_soam = get_spc_ndx( 'SOAM' )
    id_soai = get_spc_ndx( 'SOAI' )
    id_soat = get_spc_ndx( 'SOAT' )
    id_soab = get_spc_ndx( 'SOAB' )
    id_soax = get_spc_ndx( 'SOAX' )

    sox_species = (/ id_so2, id_so4, id_h2so4 /)
    bulkaero_species(:) = -1
    bulkaero_species(1:20) = (/ id_dst01, id_dst02, id_dst03, id_dst04, &
                                id_sslt01, id_sslt02, id_sslt03, id_sslt04, &
                                id_soa, id_so4, id_oc1, id_oc2, id_cb1, id_cb2, id_nh4no3, &
                                id_soam,id_soai,id_soat,id_soab,id_soax /)

    aer_species(:) = -1
    n = 1
    do m = 1,gas_pcnst
       k=0
       if ( any(bulkaero_species(:)==m) ) k=1
       if ( k==0 ) k = index(trim(solsym(m)), '_a')
       if ( k==0 ) k = index(trim(solsym(m)), '_c')
       if ( k>0 ) then ! must be aerosol species
          aer_species(n) = m
          n = n+1
       endif
    enddo

    call addfld( 'NOX', (/ 'lev' /), 'A', 'mol/mol', 'nox volume mixing ratio' )
    call addfld( 'NOY', (/ 'lev' /), 'A', 'mol/mol', 'noy volume mixing ratio' )
    call addfld( 'BROX', (/ 'lev' /), 'A','mol/mol', 'brox volume mixing ratio' )
    call addfld( 'BROY', (/ 'lev' /), 'A','mol/mol', 'total inorganic bromine (Br+BrO+HOBr+BrONO2+HBr+BrCl)' )
    call addfld( 'CLOX', (/ 'lev' /), 'A','mol/mol', 'clox volume mixing ratio' )
    call addfld( 'CLOY', (/ 'lev' /), 'A','mol/mol', 'total inorganic chlorine (Cl+ClO+2Cl2+2Cl2O2+OClO+HOCl+ClONO2+HCl+BrCl)' )
    call addfld( 'TCLY', (/ 'lev' /), 'A','mol/mol', 'total Cl volume mixing ratio' )
    call addfld( 'TOTH', (/ 'lev' /), 'A','mol/mol', 'total H2 volume mixing ratio' )

    call addfld( 'NOY_mmr', (/ 'lev' /), 'A', 'kg/kg', 'NOy mass mixing ratio' )
    call addfld( 'SOX_mmr', (/ 'lev' /), 'A', 'kg/kg', 'SOx mass mixing ratio' )
    call addfld( 'NHX_mmr', (/ 'lev' /), 'A', 'kg/kg', 'NHx mass mixing ratio' )

    do m = 1,gas_pcnst

       spc_name = trim(solsym(m))

       call cnst_get_ind(spc_name, n, abrtf=.false. )
       if ( n > 0 ) then
          attr = cnst_longname(n)
       elseif ( trim(spc_name) == 'H2O' ) then
          attr = 'water vapor'
       else
          attr = spc_name
       endif

       depvel_name(m) = 'DV_'//trim(spc_name)
       depflx_name(m) = 'DF_'//trim(spc_name)
       dtchem_name(m) = 'D'//trim(spc_name)//'CHM'

       call addfld( depvel_name(m),   horiz_only,    'A', 'cm/s', 'deposition velocity ' )
       call addfld( depflx_name(m), horiz_only,    'A', 'kg/m2/s', 'dry deposition flux ' )
       call addfld( dtchem_name(m),   (/ 'lev' /), 'A', 'kg/s', 'net tendency from chem' )

       wetdep_name(m) = 'WD_'//trim(spc_name)
       wtrate_name(m) = 'WDR_'//trim(spc_name)

       call addfld( wetdep_name(m),   horiz_only,    'A', 'kg/s', spc_name//' wet deposition' )
       call addfld( wtrate_name(m),   (/ 'lev' /), 'A',   '/s', spc_name//' wet deposition rate' )
       
       if (spc_name(1:3) == 'num') then
          unit_basename = ' 1'
       else
          unit_basename = 'kg'
       endif

       if ( any( aer_species == m ) ) then
          call addfld( spc_name,   (/ 'lev' /), 'A', unit_basename//'/kg ', trim(attr)//' concentration')
          call addfld( trim(spc_name)//'_SRF', horiz_only, 'A', unit_basename//'/kg', trim(attr)//" in bottom layer")    
       else
          call addfld( spc_name, (/ 'lev' /), 'A', 'mol/mol', trim(attr)//' concentration')
          call addfld( trim(spc_name)//'_SRF', horiz_only, 'A', 'mol/mol', trim(attr)//" in bottom layer")
       endif

       if (history_aerosol) then
          if (history_verbose .or. trim(spc_name) == 'O3' .or. trim(spc_name) == 'SO2' ) then
             call add_default( spc_name, 1, ' ' )
          endif
          call add_default( trim(spc_name)//'_SRF', 1, ' ' )
       endif 
       if (history_amwg) then
          call add_default( trim(spc_name)//'_SRF', 1, ' ' )
       endif

    enddo

    ! Add sum of mass mixing ratios for each aerosol class
    if (history_aerosol .and. .not. history_verbose) then
       call addfld( 'Mass_bc',   (/ 'lev' /), 'A', 'kg/kg ', &
            'sum of bc mass concentration bc_a1+bc_c1+bc_a3+bc_c3+bc_a4+bc_c4')
       call add_default( 'Mass_bc', 1, ' ' )
       call addfld( 'Mass_pom',   (/ 'lev' /), 'A', 'kg/kg ', &
            'sum of pom mass concentration pom_a1+pom_c1+pom_a3+pom_c3+pom_a4+pom_c4')
       call add_default( 'Mass_pom', 1, ' ' )
       call addfld( 'Mass_mom',   (/ 'lev' /), 'A', 'kg/kg ', &
            'sum of mom mass concentration mom_a1+mom_c1+mom_a2+mom_c2+mom_a3+mom_c3+mom_a4+mom_c4')
       call add_default( 'Mass_mom', 1, ' ' )
       call addfld( 'Mass_ncl',   (/ 'lev' /), 'A', 'kg/kg ', &
            'sum of ncl mass concentration ncl_a1+ncl_c1+ncl_a2+ncl_c2+ncl_a3+ncl_c3')
       call add_default( 'Mass_ncl', 1, ' ' )
       call addfld( 'Mass_soa',   (/ 'lev' /), 'A', 'kg/kg ', &
            'sum of soa mass concentration soa_a1+soa_c1+soa_a2+soa_c2+soa_a3+soa_c3')
       call add_default( 'Mass_soa', 1, ' ' )
       call addfld( 'Mass_so4',   (/ 'lev' /), 'A', 'kg/kg ', &
            'sum of so4 mass concentration so4_a1+so4_c1+so4_a2+so4_c2+so4_a3+so4_c3')
       call add_default( 'Mass_so4', 1, ' ' )
       call addfld( 'Mass_dst',   (/ 'lev' /), 'A', 'kg/kg ', &
            'sum of dst mass concentration dst_a1+dst_c1+dst_a3+dst_c3')
       call add_default( 'Mass_dst', 1, ' ' )
    endif

    call addfld( 'MASS', (/ 'lev' /), 'A', 'kg', 'mass of grid box' )
    call addfld( 'DRYMASS', (/ 'lev' /), 'A', 'kg', 'dry air mass of grid box' )
    call addfld( 'AREA', horiz_only,    'A', 'm2', 'area of grid box' )

    call addfld( 'WD_NOY', horiz_only, 'A', 'kg/s', 'NOy wet deposition' )
    call addfld( 'DF_NOY', horiz_only, 'I', 'kg/m2/s', 'NOy dry deposition flux ' )

    call addfld( 'WD_SOX', horiz_only, 'A', 'kg/s', 'SOx wet deposition' )
    call addfld( 'DF_SOX', horiz_only, 'I', 'kg/m2/s', 'SOx dry deposition flux ' )

    call addfld( 'WD_NHX', horiz_only, 'A', 'kg/s', 'NHx wet deposition' )
    call addfld( 'DF_NHX', horiz_only, 'I', 'kg/m2/s', 'NHx dry deposition flux ' )

    call addfld( 'TOZ', horiz_only,    'A', 'DU', 'Total column ozone' )
    call addfld( 'TCO', horiz_only,    'A', 'DU', 'Tropospheric column ozone based on chemistry tropopause' )
    call add_default( 'TCO', 1, ' ' )
    call addfld( 'SCO', horiz_only,    'A', 'DU', 'Stratospheric column ozone based on chemistry tropopause' )
    call add_default( 'SCO', 1, ' ' )

  end subroutine chm_diags_inti

!========================================================================
  subroutine chm_diags( lchnk, ncol, vmr, mmr, depvel, depflx, mmr_tend, pdel, pdeldry, pbuf, ltrop )
    !--------------------------------------------------------------------
    !	... utility routine to output chemistry diagnostic variables
    !--------------------------------------------------------------------
    
    use cam_history,  only : outfld
    use constituents, only : pcnst
    use constituents, only : cnst_get_ind
    use phys_grid,    only : get_area_all_p, pcols
    use physconst,    only : mwdry                   ! molecular weight of dry air
    use physics_buffer, only : physics_buffer_desc
    use modal_aero_data,  only : cnst_name_cw, qqcw_get_field ! for calculate sum of aerosol masses
    use phys_control, only: phys_getopts
    
    implicit none

    !--------------------------------------------------------------------
    !	... dummy arguments
    !--------------------------------------------------------------------
    integer,  intent(in)  :: lchnk
    integer,  intent(in)  :: ncol
    real(r8), intent(in)  :: vmr(ncol,pver,gas_pcnst)
    real(r8), intent(in)  :: mmr(ncol,pver,gas_pcnst)
    real(r8), intent(in)  :: depvel(ncol, gas_pcnst)
    real(r8), intent(in)  :: depflx(ncol, gas_pcnst)
    real(r8), intent(in)  :: mmr_tend(ncol,pver,gas_pcnst)
    real(r8), intent(in)  :: pdel(ncol,pver)
    real(r8), intent(in)  :: pdeldry(ncol,pver)
    integer,  intent(in)  :: ltrop(pcols)  ! index of the lowest stratospheric level
    type(physics_buffer_desc), pointer :: pbuf(:)

    !--------------------------------------------------------------------
    !	... local variables
    !--------------------------------------------------------------------
    integer     :: icol,kk, mm, nn
    real(r8)    :: ozone_layer(ncol,pver)   ! ozone concentration [DU]
    real(r8)    :: ozone_col(ncol)          ! vertical integration of ozone [DU]
    real(r8)    :: ozone_trop(ncol)         ! vertical integration of ozone in troposphere [DU]
    real(r8)    :: ozone_strat(ncol)        ! vertical integration of ozone in stratosphere [DU]
    
    real(r8), dimension(ncol,pver) :: vmr_nox, vmr_noy, vmr_clox, vmr_cloy, vmr_tcly, vmr_brox, vmr_broy, vmr_toth
    real(r8), dimension(ncol,pver) :: mmr_noy, mmr_sox, mmr_nhx, net_chem
    real(r8), dimension(ncol)      :: df_noy, df_sox, df_nhx

    real(r8) :: area(ncol), mass(ncol,pver), drymass(ncol,pver)
    real(r8) :: wgt
    character(len=16) :: spc_name
    real(r8), pointer :: fldcw(:,:)  !working pointer to extract data from pbuf for sum of mass for aerosol classes
    real(r8), dimension(ncol,pver) :: mass_bc, mass_dst, mass_mom, mass_ncl, mass_pom, mass_so4, mass_soa

    logical :: history_aerosol      ! output aerosol variables
    logical :: history_verbose      ! produce verbose history output
#include "../yaml/mo_chm_diags/f90_yaml/chm_diags_beg_yml.f90"

    !-----------------------------------------------------------------------

    call phys_getopts( history_aerosol_out = history_aerosol, &
                       history_verbose_out = history_verbose )
    !--------------------------------------------------------------------
    !	... "diagnostic" groups
    !--------------------------------------------------------------------
    vmr_nox(:ncol,:) = 0._r8
    vmr_noy(:ncol,:) = 0._r8
    vmr_clox(:ncol,:) = 0._r8
    vmr_cloy(:ncol,:) = 0._r8
    vmr_tcly(:ncol,:) = 0._r8
    vmr_brox(:ncol,:) = 0._r8
    vmr_broy(:ncol,:) = 0._r8
    vmr_toth(:ncol,:) = 0._r8
    mmr_noy(:ncol,:) = 0._r8
    mmr_sox(:ncol,:) = 0._r8
    mmr_nhx(:ncol,:) = 0._r8
    df_noy(:ncol) = 0._r8
    df_sox(:ncol) = 0._r8
    df_nhx(:ncol) = 0._r8

    ! Save the sum of mass mixing ratios for each class instea of individual
    ! species to reduce history file size

    ! Mass_bc = bc_a1 + bc_c1 + bc_a3 + bc_c3 + bc_a4 + bc_c4
    ! Mass_dst = dst_a1 + dst_c1 + dst_a3 + dst_c3
    ! Mass_mom = mom_a1 + mom_c1 + mom_a2 + mom_c2 + mom_a3 + mom_c3 + mom_a4 + mom_c4
    ! Mass_ncl = ncl_a1 + ncl_c1 + ncl_a2 + ncl_c2 + ncl_a3 + ncl_c3
    ! Mass_pom = pom_a1 + pom_c1 + pom_a3 + pom_c3 + pom_a4 + pom_c4
    ! Mass_so4 = so4_a1 + so4_c1 + so4_a2 + so4_c2 + so4_a3 + so4_c3
    ! Mass_soa = soa_a1 + soa_c1 + soa_a2 + soa_c2 + soa_a3 + soa_c3

    !initialize the mass arrays
    if (history_aerosol .and. .not. history_verbose) then
       mass_bc(:ncol,:) = 0._r8
       mass_dst(:ncol,:) = 0._r8
       mass_mom(:ncol,:) = 0._r8
       mass_ncl(:ncol,:) = 0._r8
       mass_pom(:ncol,:) = 0._r8
       mass_so4(:ncol,:) = 0._r8
       mass_soa(:ncol,:) = 0._r8
    endif

    call get_area_all_p(lchnk, ncol, area)
    area = area * rearth**2

    do kk = 1,pver
       mass(:ncol,kk) = pdel(:ncol,kk) * area(:ncol) * rgrav
       drymass(:ncol,kk) = pdeldry(:ncol,kk) * area(:ncol) * rgrav
    enddo

    call outfld( 'AREA', area(:ncol),   ncol, lchnk )
    call outfld( 'MASS', mass(:ncol,:), ncol, lchnk )
    call outfld( 'DRYMASS', drymass(:ncol,:), ncol, lchnk )

    ! convert ozone from mol/mol (w.r.t. dry air mass) to DU
    ozone_layer(:ncol,:) = pdeldry(:ncol,:)*vmr(:ncol,:,id_o3)*avogadro*rgrav/mwdry/DUfac*1.e3_r8
    ! total column ozone
    ozone_col(:) = 0._r8
    ozone_trop(:) = 0._r8
    ozone_strat(:) = 0._r8
    do icol = 1,ncol
       do kk = 1,pver
          ozone_col(icol) = ozone_col(icol) + ozone_layer(icol,kk)
          if (kk <= ltrop(icol)) then
             ozone_strat(icol) = ozone_strat(icol) + ozone_layer(icol,kk)
          else
             ozone_trop(icol) = ozone_trop(icol) + ozone_layer(icol,kk)
          endif
       enddo
    enddo
    call outfld( 'TOZ', ozone_col, ncol, lchnk )
    ! stratospheric column ozone
    call outfld( 'SCO', ozone_strat, ncol, lchnk )
    ! tropospheric column ozone
    call outfld( 'TCO', ozone_trop, ncol, lchnk )

    do mm = 1,gas_pcnst

      ! other options of species are not used, only use weight=1
       wgt = 1._r8

       if ( any( sox_species == mm ) ) then
          mmr_sox(:ncol,:) = mmr_sox(:ncol,:) +  wgt * mmr(:ncol,:,mm)
       endif
       
       if ( any( aer_species == mm ) ) then
          call outfld( solsym(mm), mmr(:ncol,:,mm), ncol ,lchnk )
          call outfld( trim(solsym(mm))//'_SRF', mmr(:ncol,pver,mm), ncol ,lchnk )
          if (history_aerosol .and. .not. history_verbose) then
             select case (trim(solsym(mm)))
             case ('bc_a1','bc_a3','bc_a4')
                  mass_bc(:ncol,:) = mass_bc(:ncol,:) + mmr(:ncol,:,mm)
             case ('dst_a1','dst_a3')
                  mass_dst(:ncol,:) = mass_dst(:ncol,:) + mmr(:ncol,:,mm)
             case ('mom_a1','mom_a2','mom_a3','mom_a4')
                  mass_mom(:ncol,:) = mass_mom(:ncol,:) + mmr(:ncol,:,mm)
             case ('ncl_a1','ncl_a2','ncl_a3')
                  mass_ncl(:ncol,:) = mass_ncl(:ncol,:) + mmr(:ncol,:,mm)
             case ('pom_a1','pom_a3','pom_a4')
                  mass_pom(:ncol,:) = mass_pom(:ncol,:) + mmr(:ncol,:,mm)
             case ('so4_a1','so4_a2','so4_a3')
                  mass_so4(:ncol,:) = mass_so4(:ncol,:) + mmr(:ncol,:,mm)
             case ('soa_a1','soa_a2','soa_a3')
                  mass_soa(:ncol,:) = mass_soa(:ncol,:) + mmr(:ncol,:,mm)
             endselect
          endif
       else
          call outfld( solsym(mm), vmr(:ncol,:,mm), ncol ,lchnk )
          call outfld( trim(solsym(mm))//'_SRF', vmr(:ncol,pver,mm), ncol ,lchnk )
       endif

       call outfld( depvel_name(mm), depvel(:ncol,mm), ncol ,lchnk )
       call outfld( depflx_name(mm), depflx(:ncol,mm), ncol ,lchnk )

       if ( any( sox_species == mm ) ) then
          df_sox(:ncol) = df_sox(:ncol) +  wgt * depflx(:ncol,mm)*S_molwgt/adv_mass(mm)
       endif

       do kk=1,pver
          do icol=1,ncol
             net_chem(icol,kk) = mmr_tend(icol,kk,mm) * mass(icol,kk) 
          enddo
       enddo
       call outfld( dtchem_name(mm), net_chem(:ncol,:), ncol, lchnk )

    enddo

    ! diagnostics for cloud-borne aerosols, then add to corresponding mass accumulators
    if (history_aerosol .and. .not. history_verbose) then

       do nn = 1,pcnst
          fldcw => qqcw_get_field(pbuf,nn,lchnk,errorhandle=.true.)
          if(associated(fldcw)) then
             select case (trim(cnst_name_cw(nn)))
                case ('bc_c1','bc_c3','bc_c4')
                     mass_bc(:ncol,:) = mass_bc(:ncol,:) + fldcw(:ncol,:)
                case ('dst_c1','dst_c3')
                     mass_dst(:ncol,:) = mass_dst(:ncol,:) + fldcw(:ncol,:)
                case ('mom_c1','mom_c2','mom_c3','mom_c4')
                     mass_mom(:ncol,:) = mass_mom(:ncol,:) + fldcw(:ncol,:)
                case ('ncl_c1','ncl_c2','ncl_c3')
                     mass_ncl(:ncol,:) = mass_ncl(:ncol,:) + fldcw(:ncol,:)
                case ('pom_c1','pom_c3','pom_c4')
                     mass_pom(:ncol,:) = mass_pom(:ncol,:) + fldcw(:ncol,:)
                case ('so4_c1','so4_c2','so4_c3')
                     mass_so4(:ncol,:) = mass_so4(:ncol,:) + fldcw(:ncol,:)
                case ('soa_c1','soa_c2','soa_c3')
                     mass_soa(:ncol,:) = mass_soa(:ncol,:) + fldcw(:ncol,:)
             endselect
          endif
       enddo
       call outfld( 'Mass_bc', mass_bc(:ncol,:),ncol,lchnk)
       call outfld( 'Mass_dst', mass_dst(:ncol,:),ncol,lchnk)
       call outfld( 'Mass_mom', mass_mom(:ncol,:),ncol,lchnk)
       call outfld( 'Mass_ncl', mass_ncl(:ncol,:),ncol,lchnk)
       call outfld( 'Mass_pom', mass_pom(:ncol,:),ncol,lchnk)
       call outfld( 'Mass_so4', mass_so4(:ncol,:),ncol,lchnk)
       call outfld( 'Mass_soa', mass_soa(:ncol,:),ncol,lchnk)
    endif

    call outfld( 'NOX',  vmr_nox(:ncol,:),  ncol, lchnk )
    call outfld( 'NOY',  vmr_noy(:ncol,:),  ncol, lchnk )
    call outfld( 'CLOX', vmr_clox(:ncol,:), ncol, lchnk )
    call outfld( 'CLOY', vmr_cloy(:ncol,:), ncol, lchnk )
    call outfld( 'BROX', vmr_brox(:ncol,:), ncol, lchnk )
    call outfld( 'BROY', vmr_broy(:ncol,:), ncol, lchnk )
    call outfld( 'TCLY', vmr_tcly(:ncol,:), ncol, lchnk )
    call outfld( 'NOY_mmr', mmr_noy(:ncol,:), ncol ,lchnk )
    call outfld( 'SOX_mmr', mmr_sox(:ncol,:), ncol ,lchnk )
    call outfld( 'NHX_mmr', mmr_nhx(:ncol,:), ncol ,lchnk )
    call outfld( 'DF_NOY', df_noy(:ncol), ncol ,lchnk )
    call outfld( 'DF_SOX', df_sox(:ncol), ncol ,lchnk )
    call outfld( 'DF_NHX', df_nhx(:ncol), ncol ,lchnk )
#include "../yaml/mo_chm_diags/f90_yaml/chm_diags_end_yml.f90"


  end subroutine chm_diags

!========================================================================
  subroutine het_diags( het_rates, mmr, pdel, lchnk, ncol )

    use cam_history,  only : outfld
    use phys_grid,    only : get_wght_all_p
    implicit none

    integer,  intent(in)  :: lchnk
    integer,  intent(in)  :: ncol
    real(r8), intent(in)  :: het_rates(ncol,pver,max(1,gas_pcnst))
    real(r8), intent(in)  :: mmr(ncol,pver,gas_pcnst)
    real(r8), intent(in)  :: pdel(ncol,pver)

    real(r8), dimension(ncol) :: noy_wk, sox_wk, nhx_wk, wrk_wd
    integer  :: mm, kk
    real(r8) :: wght(ncol)
#include "../yaml/mo_chm_diags/f90_yaml/het_diags_beg_yml.f90"
    !
    ! output integrated wet deposition field
    !
    noy_wk(:) = 0._r8
    sox_wk(:) = 0._r8
    nhx_wk(:) = 0._r8

    call get_wght_all_p(lchnk, ncol, wght)

    do mm = 1,gas_pcnst
       !
       ! compute vertical integral
       !
       wrk_wd(:ncol) = 0._r8
       do kk = 1,pver
          wrk_wd(:ncol) = wrk_wd(:ncol) + het_rates(:ncol,kk,mm) * mmr(:ncol,kk,mm) * pdel(:ncol,kk) 
       enddo
       !
       wrk_wd(:ncol) = wrk_wd(:ncol) * rgrav * wght(:ncol) * rearth**2
       !
       ! wrk_wd_mm(:ncol,mm) = wrk_wd  ! only use for testing files
       call outfld( wetdep_name(mm), wrk_wd(:ncol),         ncol, lchnk )
       call outfld( wtrate_name(mm), het_rates(:ncol,:,mm), ncol, lchnk )

       if ( any(sox_species == mm ) ) then
          sox_wk(:ncol) = sox_wk(:ncol) + wrk_wd(:ncol)*S_molwgt/adv_mass(mm)
       endif

    enddo
    
    call outfld( 'WD_NOY', noy_wk(:ncol), ncol, lchnk )
    call outfld( 'WD_SOX', sox_wk(:ncol), ncol, lchnk )
    call outfld( 'WD_NHX', nhx_wk(:ncol), ncol, lchnk )
#include "../yaml/mo_chm_diags/f90_yaml/het_diags_end_yml.f90"

  end subroutine het_diags
!========================================================================
end module mo_chm_diags
