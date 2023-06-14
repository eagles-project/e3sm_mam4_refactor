
module mo_setext

  use cam_logfile, only: iulog

  private
  public :: setext_inti, setext, has_ions

  save

  integer :: co_ndx, no_ndx, synoz_ndx, xno_ndx
  integer :: op_ndx, o2p_ndx, np_ndx, n2p_ndx, n2d_ndx, n_ndx, e_ndx, oh_ndx
  logical :: has_ions = .false.

contains

  subroutine setext_inti
    !--------------------------------------------------------
    !	... Initialize the external forcing module
    !--------------------------------------------------------

    use mo_chem_utls, only : get_extfrc_ndx
    use ppgrid,       only : pver
    use cam_history,  only : addfld
    use spmd_utils,   only : masterproc

    implicit none

    co_ndx    = get_extfrc_ndx( 'CO' )
    no_ndx    = get_extfrc_ndx( 'NO' )
    synoz_ndx = get_extfrc_ndx( 'SYNOZ' )
    xno_ndx   = get_extfrc_ndx( 'XNO' )

    op_ndx   = get_extfrc_ndx( 'Op' )
    o2p_ndx  = get_extfrc_ndx( 'O2p' )
    np_ndx   = get_extfrc_ndx( 'Np' )
    n2p_ndx  = get_extfrc_ndx( 'N2p' )
    n2d_ndx  = get_extfrc_ndx( 'N2D' )
    n_ndx    = get_extfrc_ndx( 'N' )
    e_ndx    = get_extfrc_ndx( 'e' )
    oh_ndx   = get_extfrc_ndx( 'OH' )

    has_ions = op_ndx > 0 .and. o2p_ndx > 0 .and. np_ndx > 0 .and. n2p_ndx > 0 .and. e_ndx > 0

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'setext_inti: diagnostics: co_ndx, no_ndx, synoz_ndx, xno_ndx'
       write(iulog,'(10i5)') co_ndx, no_ndx, synoz_ndx, xno_ndx
    endif

    call addfld( 'NO_Lightning', (/ 'lev' /), 'A','molec/cm3/s', 'lightning NO source' )
    call addfld( 'NO_Aircraft', (/ 'lev' /), 'A', 'molec/cm3/s', 'aircraft NO source' )
    call addfld( 'CO_Aircraft', (/ 'lev' /), 'A', 'molec/cm3/s', 'aircraft CO source' )

    call addfld( 'N4S_SPE', (/ 'lev' /), 'I', 'molec/cm3/s', 'solar proton event N(4S) source' )
    call addfld( 'N2D_SPE', (/ 'lev' /), 'I', 'molec/cm3/s', 'solar proton event N(2S) source' )
    call addfld( 'OH_SPE', (/ 'lev' /), 'I',  'molec/cm3/s', 'solar proton event HOx source' )

    if ( has_ions ) then
       call addfld( 'P_Op', (/ 'lev' /), 'I', '/s', 'production o+' )
       call addfld( 'P_O2p', (/ 'lev' /), 'I', '/s', 'production o2+' )
       call addfld( 'P_N2p', (/ 'lev' /), 'I', '/s', 'production n2+' )
       call addfld( 'P_Np', (/ 'lev' /), 'I', '/s', 'production n+' )
       call addfld( 'P_IONS', (/ 'lev' /), 'I', '/s', 'total ion production' )
    endif

  end subroutine setext_inti

!==============================================================================
  subroutine setext( extfrc, zint_abs, zint_rel, cldtop, &
       zmid, lchnk, tfld, o2mmr, ommr, &
       pmid, mbar, rlats, calday, ncol, rlons, pbuf )
    !--------------------------------------------------------
    !     ... for this latitude slice:
    !         - form the production from datasets
    !         - form the nox (xnox) production from lighing
    !         - form the nox (xnox) production from airplanes
    !         - form the co production from airplanes
    !--------------------------------------------------------

    use cam_history,  only : outfld
    use shr_kind_mod, only : r8 => shr_kind_r8
    use ppgrid,       only : pver, pcols
    use mo_airplane,  only : airpl_set
    use chem_mods,    only : extcnt
    use mo_lightning, only : prod_no

    use mo_extfrc,    only : extfrc_set
    use chem_mods,    only : extcnt
    use tracer_srcs,  only : num_tracer_srcs, tracer_src_flds, get_srcs_data
    use mo_chem_utls, only : get_extfrc_ndx
    use mo_synoz,     only : po3

    use mo_aurora,      only : aurora
    use mo_solarproton, only : spe_prod 
    use physics_buffer, only : physics_buffer_desc

    implicit none

    !--------------------------------------------------------
    !     ... dummy arguments
    !--------------------------------------------------------
    !--------------------------------------------------------
    integer,  intent(in)  ::   lchnk                       ! chunk id
    integer,  intent(in)  ::   ncol                        ! columns in chunk
    real(r8), intent(in)  ::   zint_abs(ncol,pver+1)           ! interface geopot height ( km )
    real(r8), intent(in)  ::   zint_rel(ncol,pver+1)           ! interface geopot height ( km )
    real(r8), intent(in)  ::   cldtop(ncol)                ! cloud top index
    real(r8), intent(out) ::   extfrc(ncol,pver,extcnt)    ! the "extraneous" forcing

    real(r8), intent(in)  ::   calday                      ! calendar day of year
    real(r8), intent(in)  ::   rlats(ncol)                 ! column latitudes (radians)
    real(r8), intent(in)  ::   rlons(ncol)                 ! column longitudes (radians)
    real(r8), intent(in)  ::   zmid(ncol,pver)             ! midpoint geopot height ( km )
    real(r8), intent(in)  ::   pmid(pcols,pver)            ! midpoint pressure (Pa)
    real(r8), intent(in)  ::   tfld(pcols,pver)            ! midpoint temperature (K)
    real(r8), intent(in)  ::   mbar(ncol,pver)             ! mean molecular mass (g/mole)
    real(r8), intent(in)  ::   o2mmr(ncol,pver)            ! o2 concentration (kg/kg)
    real(r8), intent(in)  ::   ommr(ncol,pver)             ! o concentration (kg/kg)

    type(physics_buffer_desc),pointer :: pbuf(:)

    !--------------------------------------------------------
    !     ... local variables
    !--------------------------------------------------------
    integer :: i, k, k1, kk, cldind, ii, jj, nlev
    real(r8) :: srcs_offline( ncol, pver )
    integer :: ndx

    real(r8), dimension(ncol,pver) :: no_lgt, no_air, co_air

    real(r8)    :: spe_nox(ncol,pver)        ! Solar Proton Event NO production
    real(r8)    :: spe_hox(ncol,pver)        ! Solar Proton Event HOx production

    !--------------------------------------------------------
    !     ... set frcing from datasets
    !--------------------------------------------------------
    call extfrc_set( lchnk, zint_rel, ncol, & ! in
                     extfrc ) ! out
    
    !--------------------------------------------------------
    !     ... set nox production from lighting
    !         note: from ground to cloud top production is c shaped
    !
    ! FORTRAN refactor: nox is not included in current MAM4
    ! the related code are removed but outfld is kept for BFB testing
    !--------------------------------------------------------
    no_lgt(:,:) = 0._r8
    call outfld( 'NO_Lightning', no_lgt(:ncol,:), ncol, lchnk )

    ! FORTRAN refactor: in the subroutine airpl_set, has_airpl_src is false,
    ! the subroutine only has two outfld calls that output zero
    ! remove the subroutine call and move out the zero outfld
    no_air(:,:) = 0._r8
    co_air(:,:) = 0._r8
    call outfld( 'NO_Aircraft',  no_air(:ncol,:), ncol, lchnk )
    call outfld( 'CO_Aircraft',  co_air(:ncol,:), ncol, lchnk )

  end subroutine setext

end module mo_setext
