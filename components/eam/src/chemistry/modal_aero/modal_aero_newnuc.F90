! modal_aero_newnuc.F90


!----------------------------------------------------------------------
!BOP
!
! !MODULE: modal_aero_newnuc --- modal aerosol new-particle nucleation
!
! !INTERFACE:
   module modal_aero_newnuc
#if (defined MODAL_AERO)

! !USES:
   use shr_kind_mod,  only:  r8 => shr_kind_r8
   use shr_kind_mod,  only:  r4 => shr_kind_r4
   use cam_logfile,   only:  iulog
   use mo_constants,  only:  pi
   use chem_mods,     only:  gas_pcnst
         use module_perturb
  use time_manager

  implicit none
  private
  save

! !PUBLIC MEMBER FUNCTIONS:
  public :: modal_aero_newnuc_init, &
            mer07_veh02_nuc_mosaic_1box


! !PUBLIC DATA MEMBERS:
  integer, public :: newnuc_h2so4_conc_flag = 1

! min h2so4 vapor for nuc calcs = 4.0e-16 mol/mol-air ~= 1.0e4 molecules/cm3, 
  real(r8), public, parameter :: qh2so4_cutoff = 4.0e-16_r8

! adjustment factors
  real(r8), parameter, public :: adjust_factor_dnaitdt = 1.0_r8       ! applied to final dnait/dt
  real(r8), parameter, public :: adjust_factor_bin_tern_ratenucl = 1.0_r8  !  applied to binary/ternary nucleation rate
! real(r8), parameter, public :: adjust_factor_pbl_ratenucl = 1.0_r8  ! applied to boundary layer nucleation rate
  real(r8),            public :: adjust_factor_pbl_ratenucl = 1.0_r8  ! applied to boundary layer nucleation rate


! !NON-PUBLIC DATA MEMBERS:
  integer, parameter  :: pcnstxx = gas_pcnst
  integer  :: l_h2so4_sv, l_nh3_sv, lnumait_sv, lnh4ait_sv, lso4ait_sv

! max cloud fraction for nuc calcs
  real(r8), parameter :: cld_cutoff = 0.99_r8


! !DESCRIPTION: This module implements ...
!
! !REVISION HISTORY:
!
!   R.Easter 2007.09.14:  Adapted from MIRAGE2 code
!
!EOP
!----------------------------------------------------------------------
!BOC

! list private module data here

!EOC
!----------------------------------------------------------------------


  contains

!----------------------------------------------------------------------
!-----------------------------------------------------------------------
        subroutine mer07_veh02_nuc_mosaic_1box(   print_out, ii,kk,lchnk,&
           newnuc_method_flagaa, dtnuc, temp_in, rh_in, press_in,   &
           zm_in, pblh_in,   &
           qh2so4_cur, qh2so4_avg, qnh3_cur, h2so4_uptkrate,   &
           mw_so4a_host,   &
           nsize, maxd_asize, dplom_sect, dphim_sect,   &
           isize_nuc, qnuma_del, qso4a_del, qnh4a_del,   &
           qh2so4_del, qnh3_del, dens_nh4so4a, ldiagaa,   &
           dnclusterdt )

          use mo_constants, only: rgas, &               ! Gas constant (J/K/kmol)
                                  avogad => avogadro    ! Avogadro's number (1/kmol)
          use physconst,    only: mw_so4a => mwso4, &   ! Molecular weight of sulfate
                                  mw_nh4a => mwnh4      ! Molecular weight of ammonium
!.......................................................................
!
! calculates new particle production from homogeneous nucleation
!    over timestep dtnuc, using nucleation rates from either
!    merikanto et al. (2007) h2so4-nh3-h2o ternary parameterization
!    vehkamaki et al. (2002) h2so4-h2o binary parameterization
!
! the new particles are "grown" to the lower-bound size of the host code's 
!    smallest size bin.  (this "growth" is somewhat ad hoc, and would not be
!    necessary if the host code's size bins extend down to ~1 nm.)
!
!    if the h2so4 and nh3 mass mixing ratios (mixrats) of the grown new 
!    particles exceed the current gas mixrats, the new particle production
!    is reduced so that the new particle mass mixrats match the gas mixrats.
!
!    the correction of kerminen and kulmala (2002) is applied to account
!    for loss of the new particles by coagulation as they are
!    growing to the "host code mininum size"
!
! revision history
!    coded by rc easter, pnnl, xx-apr-2007
!
! key routines called: subr ternary_nuc_napari
!
! references:
!    merikanto, j., i. napari, h. vehkamaki, t. anttila,
!     and m. kulmala, 2007, new parameterization of
!     sulfuric acid-ammonia-water ternary nucleation
!     rates at tropospheric conditions,
!       j. geophys. res., 112, d15207, doi:10.1029/2006jd0027977
!
!    vehkamäki, h., m. kulmala, i. napari, k.e.j. lehtinen,
!       c. timmreck, m. noppel and a. laaksonen, 2002,
!       an improved parameterization for sulfuric acid-water nucleation
!       rates for tropospheric and stratospheric conditions,
!       j. geophys. res., 107, 4622, doi:10.1029/2002jd002184
!
!    kerminen, v., and m. kulmala, 2002,
!	analytical formulae connecting the "real" and the "apparent"
!	nucleation rate and the nuclei number concentration
!	for atmospheric nucleation events
!
!.......................................................................
      implicit none
logical :: print_out
! subr arguments (in)
        real(r8), intent(in) :: dtnuc             ! nucleation time step (s)
        real(r8), intent(in) :: temp_in           ! temperature, in k
        real(r8), intent(in) :: rh_in             ! relative humidity, as fraction
        real(r8), intent(in) :: press_in          ! air pressure (pa)
        real(r8), intent(in) :: zm_in             ! layer midpoint height (m)
        real(r8), intent(in) :: pblh_in           ! pbl height (m)

        real(r8), intent(in) :: qh2so4_cur, qh2so4_avg
                                                  ! gas h2so4 mixing ratios (mol/mol-air)
        real(r8), intent(in) :: qnh3_cur          ! gas nh3 mixing ratios (mol/mol-air)
             ! qxxx_cur = current value (after gas chem and condensation)
             ! qxxx_avg = estimated average value (for simultaneous source/sink calcs)
        real(r8), intent(in) :: h2so4_uptkrate    ! h2so4 uptake rate to aerosol (1/s)
        real(r8), intent(in) :: mw_so4a_host      ! mw of so4 aerosol in host code (g/mol)

        integer, intent(in) :: newnuc_method_flagaa     ! 1=merikanto et al (2007) ternary
                                                        ! 2=vehkamaki et al (2002) binary
        integer, intent(in) :: nsize, ii,kk,lchnk                    ! number of aerosol size bins
        integer, intent(in) :: maxd_asize               ! dimension for dplom_sect, ...
        real(r8), intent(in) :: dplom_sect(maxd_asize)  ! dry diameter at lower bnd of bin (m)
        real(r8), intent(in) :: dphim_sect(maxd_asize)  ! dry diameter at upper bnd of bin (m)
        integer, intent(in) :: ldiagaa

! subr arguments (out)
        integer, intent(out) :: isize_nuc         ! size bin into which new particles go
        real(r8), intent(out) :: qnuma_del        ! change to aerosol number mixing ratio (#/mol-air)
        real(r8), intent(out) :: qso4a_del        ! change to aerosol so4 mixing ratio (mol/mol-air)
        real(r8), intent(out) :: qnh4a_del        ! change to aerosol nh4 mixing ratio (mol/mol-air)
        real(r8), intent(out) :: qh2so4_del       ! change to gas h2so4 mixing ratio (mol/mol-air)
        real(r8), intent(out) :: qnh3_del         ! change to gas nh3 mixing ratio (mol/mol-air)
                                                  ! aerosol changes are > 0; gas changes are < 0
        real(r8), intent(out) :: dens_nh4so4a     ! dry-density of the new nh4-so4 aerosol mass (kg/m3)
        real(r8), intent(out), optional :: &
                                 dnclusterdt      ! cluster nucleation rate (#/m3/s)

! subr arguments (out) passed via common block  
!    these are used to duplicate the outputs of yang zhang's original test driver
!    they are not really needed in wrf-chem
        real(r8) :: ratenuclt        ! j = ternary nucleation rate from napari param. (cm-3 s-1)
        real(r8) :: rateloge         ! ln (j)
        real(r8) :: cnum_h2so4       ! number of h2so4 molecules in the critical nucleus
        real(r8) :: cnum_nh3         ! number of nh3   molecules in the critical nucleus
        real(r8) :: cnum_tot         ! total number of molecules in the critical nucleus
        real(r8) :: radius_cluster   ! the radius of cluster (nm)


! local variables
        integer :: i
        integer :: igrow
        integer, save :: icase = 0, icase_reldiffmax = 0
!       integer, parameter :: ldiagaa = -1
        integer :: lun
        integer :: newnuc_method_flagaa2

        real(r8), parameter :: onethird = 1.0_r8/3.0_r8

        real(r8), parameter :: accom_coef_h2so4 = 0.65_r8   ! accomodation coef for h2so4 conden

! dry densities (kg/m3) molecular weights of aerosol 
! ammsulf, ammbisulf, and sulfacid (from mosaic  dens_electrolyte values)
!       real(r8), parameter :: dens_ammsulf   = 1.769e3
!       real(r8), parameter :: dens_ammbisulf = 1.78e3
!       real(r8), parameter :: dens_sulfacid  = 1.841e3
! use following to match cam3 modal_aero densities
        real(r8), parameter :: dens_ammsulf   = 1.770e3_r8
        real(r8), parameter :: dens_ammbisulf = 1.770e3_r8
        real(r8), parameter :: dens_sulfacid  = 1.770e3_r8

! molecular weights (g/mol) of aerosol ammsulf, ammbisulf, and sulfacid
!    for ammbisulf and sulfacid, use 114 & 96 here rather than 115 & 98
!    because we don't keep track of aerosol hion mass
        real(r8), parameter :: mw_ammsulf   = 132.0_r8
        real(r8), parameter :: mw_ammbisulf = 114.0_r8
        real(r8), parameter :: mw_sulfacid  =  96.0_r8

        real(r8), save :: reldiffmax = 0.0_r8

        real(r8) cair                     ! dry-air molar density (mol/m3)
        real(r8) cs_prime_kk              ! kk2002 "cs_prime" parameter (1/m2)
        real(r8) cs_kk                    ! kk2002 "cs" parameter (1/s)
        real(r8) dens_part                ! "grown" single-particle dry density (kg/m3)
        real(r8) dfin_kk, dnuc_kk         ! kk2002 final/initial new particle wet diameter (nm)
        real(r8) dpdry_clus               ! critical cluster diameter (m)
        real(r8) dpdry_part               ! "grown" single-particle dry diameter (m)
        real(r8) tmpa, tmpb, tmpc, tmpe, tmpq
        real(r8) tmpa1, tmpb1
        real(r8) tmp_m1, tmp_m2, tmp_m3, tmp_n1, tmp_n2, tmp_n3
        real(r8) tmp_spd                  ! h2so4 vapor molecular speed (m/s)
        real(r8) factor_kk
        real(r8) fogas, foso4a, fonh4a, fonuma
        real(r8) freduce                  ! reduction factor applied to nucleation rate
                                          ! due to limited availability of h2so4 & nh3 gases
        real(r8) freducea, freduceb
        real(r8) gamma_kk                 ! kk2002 "gamma" parameter (nm2*m2/h)
        real(r8) gr_kk                    ! kk2002 "gr" parameter (nm/h)
        real(r8) kgaero_per_moleso4a      ! (kg dry aerosol)/(mol aerosol so4)
        real(r8) mass_part                ! "grown" single-particle dry mass (kg)
        real(r8) molenh4a_per_moleso4a    ! (mol aerosol nh4)/(mol aerosol so4)
        real(r8) nh3ppt, nh3ppt_bb        ! actual and bounded nh3 (ppt)
        real(r8) nu_kk                    ! kk2002 "nu" parameter (nm)
        real(r8) qmolnh4a_del_max         ! max production of aerosol nh4 over dtnuc (mol/mol-air)
        real(r8) qmolso4a_del_max         ! max production of aerosol so4 over dtnuc (mol/mol-air)
        real(r8) ratenuclt_bb             ! nucleation rate (#/m3/s)
        real(r8) ratenuclt_kk             ! nucleation rate after kk2002 adjustment (#/m3/s)
        real(r8) rh_bb                    ! bounded value of rh_in
        real(r8) so4vol_in                ! concentration of h2so4 for nucl. calc., molecules cm-3
        real(r8) so4vol_bb                ! bounded value of so4vol_in
        real(r8) temp_bb                  ! bounded value of temp_in
        real(r8) voldry_clus              ! critical-cluster dry volume (m3)
        real(r8) voldry_part              ! "grown" single-particle dry volume (m3)
        real(r8) wetvol_dryvol            ! grown particle (wet-volume)/(dry-volume)
        real(r8) wet_volfrac_so4a         ! grown particle (dry-volume-from-so4)/(wet-volume)

if (print_out .and. kk==kprnt .and. ii==icolprnt(lchnk)) then
        write(106,*)'rgas, avogad:', rgas,avogad 
endif
!
! if h2so4 vapor < qh2so4_cutoff
! exit with new particle formation = 0
!
        isize_nuc = 1
        qnuma_del = 0.0_r8
        qso4a_del = 0.0_r8
        qnh4a_del = 0.0_r8
        qh2so4_del = 0.0_r8
        qnh3_del = 0.0_r8
        if ( present ( dnclusterdt ) ) dnclusterdt = 0.0_r8
!       if (qh2so4_avg .le. qh2so4_cutoff) return   ! this no longer needed
!       if (qh2so4_cur .le. qh2so4_cutoff) return   ! this no longer needed

        if ((newnuc_method_flagaa /=  1) .and. &
            (newnuc_method_flagaa /=  2) .and. &
            (newnuc_method_flagaa /= 11) .and. &
            (newnuc_method_flagaa /= 12)) return


!
! make call to parameterization routine
!

! calc h2so4 in molecules/cm3 and nh3 in ppt
        cair = press_in/(temp_in*rgas)
        so4vol_in  = qh2so4_avg * cair * avogad * 1.0e-6_r8
        if (print_out .and. kk==kprnt .and. ii==icolprnt(lchnk)) then
                write(106,'(A,7(ES24.15e2,","))')"nuc0:",so4vol_in, qh2so4_avg, cair, avogad, press_in, temp_in, rgas
        endif        
        nh3ppt    = qnh3_cur * 1.0e12_r8
        ratenuclt = 1.0e-38_r8
        rateloge = log( ratenuclt )
        if (print_out .and. kk==kprnt .and. ii==icolprnt(lchnk)) then
                write(106,'(A,5(ES24.15e2,","))')"nuc1:", cair, so4vol_in, nh3ppt, ratenuclt, rateloge
        endif
        if ( (newnuc_method_flagaa /=  2) .and. &
             (nh3ppt >= 0.1_r8) ) then
            newnuc_method_flagaa2 = 1

        else
! make call to vehkamaki binary parameterization routine

            if (so4vol_in >= 1.0e4_r8) then
               temp_bb = max( 230.15_r8, min( 305.15_r8, temp_in ) )
               rh_bb = max( 1.0e-4_r8, min( 1.0_r8, rh_in ) )
               so4vol_bb = max( 1.0e4_r8, min( 1.0e11_r8, so4vol_in ) )
               call binary_nuc_vehk2002(   &
                  temp_bb, rh_bb, so4vol_bb,   &
                  ratenuclt, rateloge,   &
                  cnum_h2so4, cnum_tot, radius_cluster )
            end if
            cnum_nh3 = 0.0_r8
            newnuc_method_flagaa2 = 2

        end if
        rateloge  = rateloge &
                  + log( max( 1.0e-38_r8, adjust_factor_bin_tern_ratenucl ) )


! do boundary layer nuc
        if ((newnuc_method_flagaa == 11) .or.   &
            (newnuc_method_flagaa == 12)) then
           if ( zm_in <= max(pblh_in,100.0_r8) ) then
              so4vol_bb = so4vol_in
              call pbl_nuc_wang2008( so4vol_bb,   &
                 newnuc_method_flagaa, newnuc_method_flagaa2,   &
                 ratenuclt, rateloge,   &
                 cnum_tot, cnum_h2so4, cnum_nh3, radius_cluster )
           end if
        end if


! if nucleation rate is less than 1e-6 #/cm3/s ~= 0.1 #/cm3/day,
! exit with new particle formation = 0
        if (rateloge  .le. -13.82_r8) return
!       if (ratenuclt .le. 1.0e-6) return

        ratenuclt = exp( rateloge )
        ratenuclt_bb = ratenuclt*1.0e6_r8  ! ratenuclt_bb is #/m3/s; ratenuclt is #/cm3/s
        if ( present ( dnclusterdt ) ) dnclusterdt = ratenuclt_bb


! wet/dry volume ratio - use simple kohler approx for ammsulf/ammbisulf
        tmpa = max( 0.10_r8, min( 0.95_r8, rh_in ) )
        wetvol_dryvol = 1.0_r8 - 0.56_r8/log(tmpa)


! determine size bin into which the new particles go
! (probably it will always be bin #1, but ...)
        voldry_clus = ( max(cnum_h2so4,1.0_r8)*mw_so4a + cnum_nh3*mw_nh4a ) /   &
                      (1.0e3_r8*dens_sulfacid*avogad)
! correction when host code sulfate is really ammonium bisulfate/sulfate
        voldry_clus = voldry_clus * (mw_so4a_host/mw_so4a)
        dpdry_clus = (voldry_clus*6.0_r8/pi)**onethird

        isize_nuc = 1
        dpdry_part = dplom_sect(1)
        if (dpdry_clus <= dplom_sect(1)) then
           igrow = 1   ! need to clusters to larger size
        else if (dpdry_clus >= dphim_sect(nsize)) then
           igrow = 0
           isize_nuc = nsize
           dpdry_part = dphim_sect(nsize)
        else
           igrow = 0
           do i = 1, nsize
              if (dpdry_clus < dphim_sect(i)) then
                 isize_nuc = i
                 dpdry_part = dpdry_clus
                 dpdry_part = min( dpdry_part, dphim_sect(i) )
                 dpdry_part = max( dpdry_part, dplom_sect(i) )
                 exit
              end if
           end do
        end if
        voldry_part = (pi/6.0_r8)*(dpdry_part**3)


!
! determine composition and density of the "grown particles"
! the grown particles are assumed to be liquid
!    (since critical clusters contain water)
!    so any (nh4/so4) molar ratio between 0 and 2 is allowed
! assume that the grown particles will have 
!    (nh4/so4 molar ratio) = min( 2, (nh3/h2so4 gas molar ratio) )
!
        if (igrow .le. 0) then
! no "growing" so pure sulfuric acid
           tmp_n1 = 0.0_r8
           tmp_n2 = 0.0_r8
           tmp_n3 = 1.0_r8
        else if (qnh3_cur .ge. qh2so4_cur) then
! combination of ammonium sulfate and ammonium bisulfate
! tmp_n1 & tmp_n2 = mole fractions of the ammsulf & ammbisulf
           tmp_n1 = (qnh3_cur/qh2so4_cur) - 1.0_r8
           tmp_n1 = max( 0.0_r8, min( 1.0_r8, tmp_n1 ) )
           tmp_n2 = 1.0_r8 - tmp_n1
           tmp_n3 = 0.0_r8
        else
! combination of ammonium bisulfate and sulfuric acid
! tmp_n2 & tmp_n3 = mole fractions of the ammbisulf & sulfacid
           tmp_n1 = 0.0_r8
           tmp_n2 = (qnh3_cur/qh2so4_cur)
           tmp_n2 = max( 0.0_r8, min( 1.0_r8, tmp_n2 ) )
           tmp_n3 = 1.0_r8 - tmp_n2
	end if

        tmp_m1 = tmp_n1*mw_ammsulf
        tmp_m2 = tmp_n2*mw_ammbisulf
        tmp_m3 = tmp_n3*mw_sulfacid
        dens_part = (tmp_m1 + tmp_m2 + tmp_m3)/   &
           ((tmp_m1/dens_ammsulf) + (tmp_m2/dens_ammbisulf)   &
                                  + (tmp_m3/dens_sulfacid))
        dens_nh4so4a = dens_part
        mass_part  = voldry_part*dens_part 
! (mol aerosol nh4)/(mol aerosol so4)
        molenh4a_per_moleso4a = 2.0_r8*tmp_n1 + tmp_n2  
! (kg dry aerosol)/(mol aerosol so4)
        kgaero_per_moleso4a = 1.0e-3_r8*(tmp_m1 + tmp_m2 + tmp_m3)  
! correction when host code sulfate is really ammonium bisulfate/sulfate
        kgaero_per_moleso4a = kgaero_per_moleso4a * (mw_so4a_host/mw_so4a)

! fraction of wet volume due to so4a
        tmpb = 1.0_r8 + molenh4a_per_moleso4a*17.0_r8/98.0_r8
        wet_volfrac_so4a = 1.0_r8 / ( wetvol_dryvol * tmpb )


!
! calc kerminen & kulmala (2002) correction
!
        if (igrow <=  0) then
            factor_kk = 1.0_r8

        else
! "gr" parameter (nm/h) = condensation growth rate of new particles
! use kk2002 eqn 21 for h2so4 uptake, and correct for nh3 & h2o uptake
            tmp_spd = 14.7_r8*sqrt(temp_in)   ! h2so4 molecular speed (m/s)
            gr_kk = 3.0e-9_r8*tmp_spd*mw_sulfacid*so4vol_in/   &
                    (dens_part*wet_volfrac_so4a)

! "gamma" parameter (nm2/m2/h)
! use kk2002 eqn 22
!
! dfin_kk = wet diam (nm) of grown particle having dry dia = dpdry_part (m)
            dfin_kk = 1.0e9_r8 * dpdry_part * (wetvol_dryvol**onethird)
! dnuc_kk = wet diam (nm) of cluster
            dnuc_kk = 2.0_r8*radius_cluster
            dnuc_kk = max( dnuc_kk, 1.0_r8 )
! neglect (dmean/150)**0.048 factor, 
! which should be very close to 1.0 because of small exponent
            gamma_kk = 0.23_r8 * (dnuc_kk)**0.2_r8   &
                     * (dfin_kk/3.0_r8)**0.075_r8   &
                     * (dens_part*1.0e-3_r8)**(-0.33_r8)   &
                     * (temp_in/293.0_r8)**(-0.75_r8)

! "cs_prime parameter" (1/m2) 
! instead kk2002 eqn 3, use
!     cs_prime ~= tmpa / (4*pi*tmpb * h2so4_accom_coef)
! where
!     tmpa = -d(ln(h2so4))/dt by conden to particles   (1/h units)
!     tmpb = h2so4 vapor diffusivity (m2/h units)
! this approx is generally within a few percent of the cs_prime
!     calculated directly from eqn 2, 
!     which is acceptable, given overall uncertainties
! tmpa = -d(ln(h2so4))/dt by conden to particles   (1/h units)
            tmpa = h2so4_uptkrate * 3600.0_r8
            tmpa1 = tmpa
            tmpa = max( tmpa, 0.0_r8 )
! tmpb = h2so4 gas diffusivity (m2/s, then m2/h)
            tmpb = 6.7037e-6_r8 * (temp_in**0.75_r8) / cair
            tmpb1 = tmpb         ! m2/s
            tmpb = tmpb*3600.0_r8   ! m2/h
            cs_prime_kk = tmpa/(4.0_r8*pi*tmpb*accom_coef_h2so4)
            cs_kk = cs_prime_kk*4.0_r8*pi*tmpb1

! "nu" parameter (nm) -- kk2002 eqn 11
            nu_kk = gamma_kk*cs_prime_kk/gr_kk
! nucleation rate adjustment factor (--) -- kk2002 eqn 13
            factor_kk = exp( (nu_kk/dfin_kk) - (nu_kk/dnuc_kk) )

        end if
        ratenuclt_kk = ratenuclt_bb*factor_kk


! max production of aerosol dry mass (kg-aero/m3-air)
        tmpa = max( 0.0_r8, (ratenuclt_kk*dtnuc*mass_part) )
! max production of aerosol so4 (mol-so4a/mol-air)
        tmpe = tmpa/(kgaero_per_moleso4a*cair)
! max production of aerosol so4 (mol/mol-air)
! based on ratenuclt_kk and mass_part
        qmolso4a_del_max = tmpe

! check if max production exceeds available h2so4 vapor
        freducea = 1.0_r8
        if (qmolso4a_del_max .gt. qh2so4_cur) then
           freducea = qh2so4_cur/qmolso4a_del_max
        end if

! check if max production exceeds available nh3 vapor
        freduceb = 1.0_r8
        if (molenh4a_per_moleso4a .ge. 1.0e-10_r8) then
! max production of aerosol nh4 (ppm) based on ratenuclt_kk and mass_part
           qmolnh4a_del_max = qmolso4a_del_max*molenh4a_per_moleso4a
           if (qmolnh4a_del_max .gt. qnh3_cur) then
              freduceb = qnh3_cur/qmolnh4a_del_max
           end if
        end if
        freduce = min( freducea, freduceb )

! if adjusted nucleation rate is less than 1e-12 #/m3/s ~= 0.1 #/cm3/day,
! exit with new particle formation = 0
        if (freduce*ratenuclt_kk .le. 1.0e-12_r8) return


! note:  suppose that at this point, freduce < 1.0 (no gas-available 
!    constraints) and molenh4a_per_moleso4a < 2.0
! if the gas-available constraints is do to h2so4 availability,
!    then it would be possible to condense "additional" nh3 and have
!    (nh3/h2so4 gas molar ratio) < (nh4/so4 aerosol molar ratio) <= 2 
! one could do some additional calculations of 
!    dens_part & molenh4a_per_moleso4a to realize this
! however, the particle "growing" is a crude approximate way to get
!    the new particles to the host code's minimum particle size,
! are such refinements worth the effort?


! changes to h2so4 & nh3 gas (in mol/mol-air), limited by amounts available
        tmpa = 0.9999_r8
        qh2so4_del = min( tmpa*qh2so4_cur, freduce*qmolso4a_del_max )
        qnh3_del   = min( tmpa*qnh3_cur, qh2so4_del*molenh4a_per_moleso4a )
        qh2so4_del = -qh2so4_del
        qnh3_del   = -qnh3_del

! changes to so4 & nh4 aerosol (in mol/mol-air)
        qso4a_del = -qh2so4_del
        qnh4a_del =   -qnh3_del
! change to aerosol number (in #/mol-air)
        qnuma_del = 1.0e-3_r8*(qso4a_del*mw_so4a + qnh4a_del*mw_nh4a)/mass_part

! do the following (tmpa, tmpb, tmpc) calculations as a check
! max production of aerosol number (#/mol-air)
        tmpa = max( 0.0_r8, (ratenuclt_kk*dtnuc/cair) )
! adjusted production of aerosol number (#/mol-air)
        tmpb = tmpa*freduce
! relative difference from qnuma_del
        tmpc = (tmpb - qnuma_del)/max(tmpb, qnuma_del, 1.0e-35_r8)


!
! diagnostic output to fort.41
! (this should be commented-out or deleted in the wrf-chem version)
!
        if (ldiagaa <= 0) return

        icase = icase + 1
        if (abs(tmpc) .gt. abs(reldiffmax)) then
           reldiffmax = tmpc
           icase_reldiffmax = icase
        end if
!       do lun = 41, 51, 10
        do lun = 6, 6
!          write(lun,'(/)')
           write(lun,'(a,2i9,1p,e10.2)')   &
               'vehkam bin-nuc icase, icase_rdmax =',   &
               icase, icase_reldiffmax, reldiffmax
           if (freduceb .lt. freducea) then
              if (abs(freducea-freduceb) .gt.   &
                   3.0e-7_r8*max(freduceb,freducea)) write(lun,'(a,1p,2e15.7)')   &
                 'freducea, b =', freducea, freduceb
           end if
        end do

! output factors so that output matches that of ternucl03
!       fogas  = 1.0e6                     ! convert mol/mol-air to ppm
!       foso4a = 1.0e9*mw_so4a/mw_air      ! convert mol-so4a/mol-air to ug/kg-air
!       fonh4a = 1.0e9*mw_nh4a/mw_air      ! convert mol-nh4a/mol-air to ug/kg-air
!       fonuma = 1.0e3/mw_air              ! convert #/mol-air to #/kg-air
        fogas  = 1.0_r8
        foso4a = 1.0_r8
        fonh4a = 1.0_r8
        fonuma = 1.0_r8

!       do lun = 41, 51, 10
        do lun = 6, 6

        write(lun,'(a,2i5)') 'newnuc_method_flagaa/aa2',   &
           newnuc_method_flagaa, newnuc_method_flagaa2

        write(lun,9210)
        write(lun,9201) temp_in, rh_in,   &
           ratenuclt, 2.0_r8*radius_cluster*1.0e-7_r8, dpdry_part*1.0e2_r8,   &
           voldry_part*1.0e6_r8, float(igrow)
        write(lun,9215)
        write(lun,9201)   &
           qh2so4_avg*fogas, 0.0_r8,  &
           qh2so4_cur*fogas, qnh3_cur*fogas,  &
           qh2so4_del*fogas, qnh3_del*fogas,  &
           qso4a_del*foso4a, qnh4a_del*fonh4a

        write(lun,9220)
        write(lun,9201)   &
           dtnuc, dens_nh4so4a*1.0e-3_r8,   &
           (qnh3_cur/qh2so4_cur), molenh4a_per_moleso4a,   &
           qnuma_del*fonuma, tmpb*fonuma, tmpc, freduce

        end do

!       lun = 51
        lun = 6
        write(lun,9230)
        write(lun,9201)   &
           press_in, cair*1.0e-6_r8, so4vol_in,   &
           wet_volfrac_so4a, wetvol_dryvol, dens_part*1.0e-3_r8

        if (igrow > 0) then
        write(lun,9240)
        write(lun,9201)   &
           tmp_spd, gr_kk, dnuc_kk, dfin_kk,   &
           gamma_kk, tmpa1, tmpb1, cs_kk

        write(lun,9250)
        write(lun,9201)   &
           cs_prime_kk, nu_kk, factor_kk, ratenuclt,   &
           ratenuclt_kk*1.0e-6_r8
        end if

9201    format ( 1p, 40e10.2  )
9210    format (   &
        '      temp        rh',   &
        '   ratenuc  dia_clus ddry_part',   &
        ' vdry_part     igrow' )
9215    format (   &
        '  h2so4avg  h2so4pre',   &
        '  h2so4cur   nh3_cur',   &
        '  h2so4del   nh3_del',   &
        '  so4a_del  nh4a_del' )
9220    format (    &
        '     dtnuc    dens_a   nh/so g   nh/so a',   &
        '  numa_del  numa_dl2   reldiff   freduce' )
9230    format (   &
        '  press_in      cair so4_volin',   &
        ' wet_volfr wetv_dryv dens_part' )
9240    format (   &
        '   tmp_spd     gr_kk   dnuc_kk   dfin_kk',   &
        '  gamma_kk     tmpa1     tmpb1     cs_kk' )
9250    format (   &
        ' cs_pri_kk     nu_kk factor_kk ratenuclt',   &
        ' ratenu_kk' )


        return
        end subroutine mer07_veh02_nuc_mosaic_1box



!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        subroutine pbl_nuc_wang2008( so4vol,   &
            newnuc_method_flagaa, newnuc_method_flagaa2,   &
            ratenucl, rateloge,   &
            cnum_tot, cnum_h2so4, cnum_nh3, radius_cluster )
!
! calculates boundary nucleation nucleation rate
! using the first or second-order parameterization in  
!     wang, m., and j.e. penner, 2008,
!        aerosol indirect forcing in a global model with particle nucleation,
!        atmos. chem. phys. discuss., 8, 13943-13998
!
        implicit none

! subr arguments (in)
        real(r8), intent(in) :: so4vol            ! concentration of h2so4 (molecules cm-3)
        integer, intent(in)  :: newnuc_method_flagaa  
                                ! [11,12] value selects [first,second]-order parameterization

! subr arguments (inout)
        integer, intent(inout)  :: newnuc_method_flagaa2
        real(r8), intent(inout) :: ratenucl         ! binary nucleation rate, j (# cm-3 s-1)
        real(r8), intent(inout) :: rateloge         ! log( ratenucl )

        real(r8), intent(inout) :: cnum_tot         ! total number of molecules
                                                    ! in the critical nucleus
        real(r8), intent(inout) :: cnum_h2so4       ! number of h2so4 molecules
        real(r8), intent(inout) :: cnum_nh3         ! number of nh3 molecules
        real(r8), intent(inout) :: radius_cluster   ! the radius of cluster (nm)


! local variables
        real(r8) :: tmp_diam, tmp_mass, tmp_volu
        real(r8) :: tmp_rateloge, tmp_ratenucl

! executable


! nucleation rate
        if (newnuc_method_flagaa == 11) then
           tmp_ratenucl = 1.0e-6_r8 * so4vol
        else if (newnuc_method_flagaa == 12) then
           tmp_ratenucl = 1.0e-12_r8 * (so4vol**2)
        else
           return
        end if
        tmp_ratenucl = tmp_ratenucl * adjust_factor_pbl_ratenucl
        tmp_rateloge = log( max( 1.0e-38_r8, tmp_ratenucl ) )

! exit if pbl nuc rate is lower than (incoming) ternary/binary rate
        if (tmp_rateloge <= rateloge) return

        rateloge = tmp_rateloge
        ratenucl = tmp_ratenucl
        newnuc_method_flagaa2 = newnuc_method_flagaa

! following wang 2002, assume fresh nuclei are 1 nm diameter
!    subsequent code will "grow" them to aitken mode size
        radius_cluster = 0.5_r8

! assume fresh nuclei are pure h2so4
!    since aitken size >> initial size, the initial composition 
!    has very little impact on the results
        tmp_diam = radius_cluster * 2.0e-7_r8   ! diameter in cm
        tmp_volu = (tmp_diam**3) * (pi/6.0_r8)  ! volume in cm^3
        tmp_mass = tmp_volu * 1.8_r8            ! mass in g
        cnum_h2so4 = (tmp_mass / 98.0_r8) * 6.023e23_r8   ! no. of h2so4 molec assuming pure h2so4
        cnum_tot = cnum_h2so4
        cnum_nh3 = 0.0_r8


        return
        end subroutine pbl_nuc_wang2008



!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        subroutine binary_nuc_vehk2002( temp, rh, so4vol,   &
            ratenucl, rateloge,   &
            cnum_h2so4, cnum_tot, radius_cluster )
!
! calculates binary nucleation rate and critical cluster size
! using the parameterization in  
!     vehkamäki, h., m. kulmala, i. napari, k.e.j. lehtinen,
!        c. timmreck, m. noppel and a. laaksonen, 2002,
!        an improved parameterization for sulfuric acid-water nucleation
!        rates for tropospheric and stratospheric conditions,
!        j. geophys. res., 107, 4622, doi:10.1029/2002jd002184
!
        implicit none

! subr arguments (in)
        real(r8), intent(in) :: temp              ! temperature (k)  
        real(r8), intent(in) :: rh                ! relative humidity (0-1)
        real(r8), intent(in) :: so4vol            ! concentration of h2so4 (molecules cm-3)

! subr arguments (out)
        real(r8), intent(out) :: ratenucl         ! binary nucleation rate, j (# cm-3 s-1)
        real(r8), intent(out) :: rateloge         ! log( ratenucl )

        real(r8), intent(out) :: cnum_h2so4       ! number of h2so4 molecules
                                                  ! in the critical nucleus
        real(r8), intent(out) :: cnum_tot         ! total number of molecules
                                                  ! in the critical nucleus
        real(r8), intent(out) :: radius_cluster   ! the radius of cluster (nm)


! local variables
        real(r8) :: crit_x
        real(r8) :: acoe, bcoe, ccoe, dcoe, ecoe, fcoe, gcoe, hcoe, icoe, jcoe
        real(r8) :: tmpa, tmpb

! executable


! calc sulfuric acid mole fraction in critical cluster
        crit_x = 0.740997_r8 - 0.00266379_r8 * temp   &
               - 0.00349998_r8 * log (so4vol)   &
               + 0.0000504022_r8 * temp * log (so4vol)   &
               + 0.00201048_r8 * log (rh)   &
               - 0.000183289_r8 * temp * log (rh)   &
               + 0.00157407_r8 * (log (rh)) ** 2.0_r8   &
               - 0.0000179059_r8 * temp * (log (rh)) ** 2.0_r8   &
               + 0.000184403_r8 * (log (rh)) ** 3.0_r8   &
               - 1.50345e-6_r8 * temp * (log (rh)) ** 3.0_r8


! calc nucleation rate
        acoe    = 0.14309_r8+2.21956_r8*temp   &
                - 0.0273911_r8 * temp**2.0_r8   &
                + 0.0000722811_r8 * temp**3.0_r8 + 5.91822_r8/crit_x

        bcoe    = 0.117489_r8 + 0.462532_r8 *temp   &
                - 0.0118059_r8 * temp**2.0_r8   &
                + 0.0000404196_r8 * temp**3.0_r8 + 15.7963_r8/crit_x

        ccoe    = -0.215554_r8-0.0810269_r8 * temp   &
                + 0.00143581_r8 * temp**2.0_r8   &
                - 4.7758e-6_r8 * temp**3.0_r8   &
                - 2.91297_r8/crit_x

        dcoe    = -3.58856_r8+0.049508_r8 * temp   &
                - 0.00021382_r8 * temp**2.0_r8   &
                + 3.10801e-7_r8 * temp**3.0_r8   &
                - 0.0293333_r8/crit_x

        ecoe    = 1.14598_r8 - 0.600796_r8 * temp   &
                + 0.00864245_r8 * temp**2.0_r8   &
                - 0.0000228947_r8 * temp**3.0_r8   &
                - 8.44985_r8/crit_x

        fcoe    = 2.15855_r8 + 0.0808121_r8 * temp   &
                -0.000407382_r8 * temp**2.0_r8   &
                -4.01957e-7_r8 * temp**3.0_r8   &
                + 0.721326_r8/crit_x

        gcoe    = 1.6241_r8 - 0.0160106_r8 * temp   &
                + 0.0000377124_r8 * temp**2.0_r8   &
                + 3.21794e-8_r8 * temp**3.0_r8   &
                - 0.0113255_r8/crit_x

        hcoe    = 9.71682_r8 - 0.115048_r8 * temp   &
                + 0.000157098_r8 * temp**2.0_r8   &
                + 4.00914e-7_r8 * temp**3.0_r8   &
                + 0.71186_r8/crit_x

        icoe    = -1.05611_r8 + 0.00903378_r8 * temp   &
                - 0.0000198417_r8 * temp**2.0_r8   &
                + 2.46048e-8_r8  * temp**3.0_r8   &
                - 0.0579087_r8/crit_x

        jcoe    = -0.148712_r8 + 0.00283508_r8 * temp   &
                - 9.24619e-6_r8  * temp**2.0_r8   &
                + 5.00427e-9_r8 * temp**3.0_r8   &
                - 0.0127081_r8/crit_x

        tmpa     =     (   &
                  acoe   &
                + bcoe * log (rh)   &
                + ccoe * ( log (rh))**2.0_r8   &
                + dcoe * ( log (rh))**3.0_r8   &
                + ecoe * log (so4vol)   &
                + fcoe * (log (rh)) * (log (so4vol))   &
                + gcoe * ((log (rh) ) **2.0_r8)   &
                       * (log (so4vol))   &
                + hcoe * (log (so4vol)) **2.0_r8   &
                + icoe * log (rh)   &
                       * ((log (so4vol)) **2.0_r8)   &
                + jcoe * (log (so4vol)) **3.0_r8   &
                )
        rateloge = tmpa
        tmpa = min( tmpa, log(1.0e38_r8) )
        ratenucl = exp ( tmpa )
!       write(*,*) 'tmpa, ratenucl =', tmpa, ratenucl



! calc number of molecules in critical cluster
        acoe    = -0.00295413_r8 - 0.0976834_r8*temp   &
                + 0.00102485_r8 * temp**2.0_r8   &
                - 2.18646e-6_r8 * temp**3.0_r8 - 0.101717_r8/crit_x

        bcoe    = -0.00205064_r8 - 0.00758504_r8*temp   &
                + 0.000192654_r8 * temp**2.0_r8   &
                - 6.7043e-7_r8 * temp**3.0_r8 - 0.255774_r8/crit_x

        ccoe    = +0.00322308_r8 + 0.000852637_r8 * temp   &
                - 0.0000154757_r8 * temp**2.0_r8   &
                + 5.66661e-8_r8 * temp**3.0_r8   &
                + 0.0338444_r8/crit_x

        dcoe    = +0.0474323_r8 - 0.000625104_r8 * temp   &
                + 2.65066e-6_r8 * temp**2.0_r8   &
                - 3.67471e-9_r8 * temp**3.0_r8   &
                - 0.000267251_r8/crit_x

        ecoe    = -0.0125211_r8 + 0.00580655_r8 * temp   &
                - 0.000101674_r8 * temp**2.0_r8   &
                + 2.88195e-7_r8 * temp**3.0_r8   &
                + 0.0942243_r8/crit_x

        fcoe    = -0.038546_r8 - 0.000672316_r8 * temp   &
                + 2.60288e-6_r8 * temp**2.0_r8   &
                + 1.19416e-8_r8 * temp**3.0_r8   &
                - 0.00851515_r8/crit_x

        gcoe    = -0.0183749_r8 + 0.000172072_r8 * temp   &
                - 3.71766e-7_r8 * temp**2.0_r8   &
                - 5.14875e-10_r8 * temp**3.0_r8   &
                + 0.00026866_r8/crit_x

        hcoe    = -0.0619974_r8 + 0.000906958_r8 * temp   &
                - 9.11728e-7_r8 * temp**2.0_r8   &
                - 5.36796e-9_r8 * temp**3.0_r8   &
                - 0.00774234_r8/crit_x

        icoe    = +0.0121827_r8 - 0.00010665_r8 * temp   &
                + 2.5346e-7_r8 * temp**2.0_r8   &
                - 3.63519e-10_r8 * temp**3.0_r8   &
                + 0.000610065_r8/crit_x

        jcoe    = +0.000320184_r8 - 0.0000174762_r8 * temp   &
                + 6.06504e-8_r8 * temp**2.0_r8   &
                - 1.4177e-11_r8 * temp**3.0_r8   &
                + 0.000135751_r8/crit_x

        cnum_tot = exp (   &
                  acoe   &
                + bcoe * log (rh)   &
                + ccoe * ( log (rh))**2.0_r8   &
                + dcoe * ( log (rh))**3.0_r8   &
                + ecoe * log (so4vol)   &
                + fcoe * (log (rh)) * (log (so4vol))   &
                + gcoe * ((log (rh) ) **2.0_r8)   &
                       * (log (so4vol))   &
                + hcoe * (log (so4vol)) **2.0_r8   &
                + icoe * log (rh)   &
                       * ((log (so4vol)) **2.0_r8)   &
                + jcoe * (log (so4vol)) **3.0_r8   &
                )

        cnum_h2so4 = cnum_tot * crit_x

!   calc radius (nm) of critical cluster
        radius_cluster = exp( -1.6524245_r8 + 0.42316402_r8*crit_x   &
                              + 0.3346648_r8*log(cnum_tot) )
      

      return
      end subroutine binary_nuc_vehk2002



!----------------------------------------------------------------------
!----------------------------------------------------------------------
subroutine modal_aero_newnuc_init( mam_amicphys_optaa )

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

use modal_aero_data

use cam_abortutils,   only:  endrun
use cam_history,  only:  addfld, horiz_only, add_default, fieldname_len
use constituents, only:  pcnst, cnst_get_ind, cnst_name
use spmd_utils,   only:  masterproc
use phys_control, only: phys_getopts


implicit none

!-----------------------------------------------------------------------
! arguments
   integer, intent(in) :: mam_amicphys_optaa

!-----------------------------------------------------------------------
! local
   integer  :: l_h2so4, l_nh3
   integer  :: lnumait, lnh4ait, lso4ait
   integer  :: l
   integer  :: m, mait

   character(len=fieldname_len)   :: tmpname
   character(len=fieldname_len+3) :: fieldname
   character(128)                 :: long_name
   character(8)                   :: unit

   logical                        :: dotend(pcnst)
   logical                        :: history_aerosol      ! Output the MAM aerosol tendencies

   !-----------------------------------------------------------------------     
   

!   set these indices
!   skip if no h2so4 species
!   skip if no aitken mode so4 or num species
	l_h2so4_sv = 0
	l_nh3_sv = 0
	lnumait_sv = 0
	lnh4ait_sv = 0
	lso4ait_sv = 0

	call cnst_get_ind( 'H2SO4', l_h2so4, .false. )
	call cnst_get_ind( 'NH3', l_nh3, .false. )

	mait = modeptr_aitken
	if (mait > 0) then
	    lnumait = numptr_amode(mait)
	    lso4ait = lptr_so4_a_amode(mait)
	    lnh4ait = lptr_nh4_a_amode(mait)
	end if
	if ((l_h2so4  <= 0) .or. (l_h2so4 > pcnst)) then
	    write(iulog,'(/a/)')   &
		'*** modal_aero_newnuc bypass -- l_h2so4 <= 0'
	    return
	else if ((lso4ait <= 0) .or. (lso4ait > pcnst)) then
	    write(iulog,'(/a/)')   &
		'*** modal_aero_newnuc bypass -- lso4ait <= 0'
	    return
	else if ((lnumait <= 0) .or. (lnumait > pcnst)) then
	    write(iulog,'(/a/)')   &
		'*** modal_aero_newnuc bypass -- lnumait <= 0'
	    return
	else if ((mait <= 0) .or. (mait > ntot_amode)) then
	    write(iulog,'(/a/)')   &
		'*** modal_aero_newnuc bypass -- modeptr_aitken <= 0'
	    return
	end if

	l_h2so4_sv = l_h2so4
	l_nh3_sv   = l_nh3
	lnumait_sv = lnumait
	lnh4ait_sv = lnh4ait
	lso4ait_sv = lso4ait

!
!   create history file column-tendency fields
!
        if (mam_amicphys_optaa > 0) return

        call phys_getopts( history_aerosol_out = history_aerosol )

	dotend(:) = .false.
	dotend(lnumait) = .true.
	dotend(lso4ait) = .true.
	dotend(l_h2so4) = .true.
	if ((l_nh3   > 0) .and. (l_nh3   <= pcnst) .and. &
	    (lnh4ait > 0) .and. (lnh4ait <= pcnst)) then
	    dotend(lnh4ait) = .true.
	    dotend(l_nh3) = .true.
	end if

	do l = 1, pcnst
	    if ( .not. dotend(l) ) cycle
	    tmpname = cnst_name(l)
	    unit = 'kg/m2/s'
	    do m = 1, ntot_amode
	        if (l == numptr_amode(m)) unit = '#/m2/s'
	    end do
	    fieldname = trim(tmpname) // '_sfnnuc1'
	    long_name = trim(tmpname) // ' modal_aero new particle nucleation column tendency'
	    call addfld( fieldname, horiz_only, 'A', unit, long_name )
            if ( history_aerosol ) then 
               call add_default( fieldname, 1, ' ' )
            endif
	    if ( masterproc ) write(iulog,'(3(a,2x))') &
		'modal_aero_newnuc_init addfld', fieldname, unit
	end do ! l = ...


      return
      end subroutine modal_aero_newnuc_init

!----------------------------------------------------------------------
#endif
   end module modal_aero_newnuc



