! TODO:
! 1) Set up qneg fields as a namelist variable
! 2) Handle updating qneg state using a subroutine instead of in physics update
!    and wherever qneg3 is called.
! 3) possible handle if qneg write stuff is actually called using if statements,
!    to avoid pointless do loops that aren't needed if no output is used.

module qneg

   use shr_kind_mod, only: r8 => shr_kind_r8
   use phys_grid,    only: get_lat_p, get_lon_p, get_rlat_p, get_rlon_p
   use cam_logfile,  only: iulog
   use constituents,       only: pcnst, cnst_name
   use ppgrid,             only: pcols, pver
!   use cam_history,        only: addfld, outfld

implicit none
private
save

public ::         &
   qneg3,         &
   qneg4,         &
   qneg_register, &
   massborrow,    &
   qqflx_fixer,   &
   qneg_ind

! Public fields
integer, public                   :: qneg3_numflds, qneg4_numflds
integer, public                   :: qneg3_numcnst, qneg4_numcnst
! Private fields
character*40, allocatable :: qneg3_flds(:), qneg3_cnst(:)  ! List of fields/constituents which call qneg3 errors
character*40, allocatable :: qneg4_flds(:), qneg4_cnst(:)  ! List of fields/constituents which call qneg4 errors

contains

!===============================================================================
subroutine qneg_register()
! TODO: change so that qneg3 and qneg4 fields are controlled by a namelist

   use spmd_utils,         only: masterproc
   use cam_logfile,        only: iulog
   use cam_history_support, only: add_hist_coord

   integer :: k
   character*40 :: tmp3_flds(30), tmp4_flds(30)


   !----  Add each QNEG3 field
   k=0

   k=k+1
   tmp3_flds(k) = 'zm_conv_tend'
   k=k+1
   tmp3_flds(k) = 'convect_shallow'
   k=k+1
   tmp3_flds(k) = 'clubb'
   k=k+1
   tmp3_flds(k) = 'cldwat'
   k = k+1
   tmp3_flds(k) = 'TPHYSBCb'
   k = k+1
   tmp3_flds(k) = 'TPHYSBCc'

   qneg3_numflds = k

   !----  Add each QNEG4 field
   k=0
   k=k+1
   tmp4_flds(k) = 'TPHYSAC'

   qneg4_numflds = k

   allocate(qneg3_flds(qneg3_numflds))
   allocate(qneg4_flds(qneg4_numflds))
   do k = 1,qneg3_numflds
      qneg3_flds(k) = tmp3_flds(k)
   end do
   do k = 1,qneg4_numflds
      qneg4_flds(k) = tmp4_flds(k)
   end do

   !----  Add each QNEG3 species
   k=0
   k=k+1
   tmp3_flds(k) = trim('Q') !trim(cnst_name(1))
   k=k+1
   tmp3_flds(k) = trim('CLDLIQ') !trim(cnst_name(2))
   k=k+1
   tmp3_flds(k) = trim('CLDICE') !trim(cnst_name(3))
   
   qneg3_numcnst = k
   !----  Add each QNEG4 species
   k=0
   k=k+1
   tmp4_flds(k) = trim('Q') !trim(cnst_name(1))
   
   qneg4_numcnst = k

   allocate(qneg3_cnst(qneg3_numcnst))
   allocate(qneg4_cnst(qneg4_numcnst))
   do k = 1,qneg3_numcnst
      qneg3_cnst(k) = tmp3_flds(k)
   end do
   do k = 1,qneg4_numcnst
      qneg4_cnst(k) = tmp4_flds(k)
   end do

   ! Register qneg hist coordinates
   call add_hist_coord('qneg3',  qneg3_numflds, 'qneg3 field length')
   call add_hist_coord('qneg4',  qneg4_numflds, 'qneg4 field length')
   call add_hist_coord('qneg3c', qneg3_numcnst, 'qneg3 species length')
   call add_hist_coord('qneg4c', qneg4_numcnst, 'qneg4 species length')

   ! Write key for QNEG error output to ATM logfile
   if (masterproc) then
      write(iulog,*) '===================================================='
      write(iulog,*) 'QNEG3 Errors Output Master List'
      write(iulog,*) 'QNEG3(time,#fields,#const.,nlev,ncols)'
      write(iulog,*) '-------------------------------'
      write(iulog,'(A24,A6)') 'Field', 'Index'
      do k = 1,qneg3_numflds
         write(iulog,'(A24,I6)') qneg3_flds(k), k
      end do ! k
      write(iulog,*) '-------------------------------'
      write(iulog,'(A24,A6)') 'Constituent', 'Index'
      do k = 1,qneg3_numcnst
         write(iulog,'(A24,I6)') qneg3_cnst(k), k
      end do ! k
      write(iulog,*) '-------------------------------'
      write(iulog,*) 'QNEG4 Errors Output Master List'
      write(iulog,*) 'QNEG4(time,#fields,#const,ncols)'
      write(iulog,*) '-------------------------------'
      write(iulog,'(A24,A6)') 'Field', 'Index'
      do k = 1,qneg4_numflds
         write(iulog,'(A24,I6)') qneg4_flds(k), k
      end do ! k
      write(iulog,*) '-------------------------------'
      write(iulog,'(A24,A6)') 'Constituent', 'Index'
      do k = 1,qneg4_numcnst
         write(iulog,'(A24,I6)') qneg4_cnst(k), k
      end do ! k
      write(iulog,*) '===================================================='
   end if ! masterproc

end subroutine qneg_register
!===============================================================================
subroutine qneg_ind(name,m,qtype,ind,indc)

   character (len=*), intent(in) :: name
   integer,intent(in)       :: qtype,m
   integer,intent(out)      :: ind,indc
   integer                  :: k

   ind = 0
   indc = 0

   ! First check for the field
   select case(qtype)
   case(3)
      do k = 1,qneg3_numflds
         if (trim(name).eq.trim(qneg3_flds(k))) then
            ind = k
            goto 999
         end if ! name
      end do ! k
   case(4)
      do k = 1,qneg4_numflds
         if (trim(name).eq.trim(qneg4_flds(k))) then
            ind = k
            goto 999
         end if ! name
      end do ! k
   case default
      ind = 0  ! Need to add logic here for error message
   end select

   ! Then check for the constituent
999 if (ind.eq.0) return ! No need to check, this won't be recorded

   select case(qtype)
      case(3)
         do k = 1,qneg3_numcnst
            if (trim(cnst_name(m)).eq.trim(qneg3_cnst(k))) then
               indc = k
               return
            end if ! name
         end do ! k
      case(4)
         do k = 1,qneg4_numcnst
            if (trim(cnst_name(m)).eq.trim(qneg4_cnst(k))) then
               indc = k
               return
            end if ! name
         end do ! k
      case default
         indc = 0  ! Need to add logic here for error message
   end select

   return

end subroutine qneg_ind
!===============================================================================
subroutine qneg3 (subnam  ,idx     ,ncol    ,ncold   ,lver    ,lconst_beg  , &
                  lconst_end       ,qmin    ,q       ,lfix, qneg3mat )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Check moisture and tracers for minimum value, reset any below
! minimum value to minimum value and return information to allow
! warning message to be printed. The global average is NOT preserved.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: J. Rosinski
! 
! Modifications: 
!
!  2016-08   Kai Zhang (kai.zhang@pnnl.gov) 
!            Added an option to print out the information about negative 
!            values without clipping the tracer concentration. This can
!            be used together with subroutine massborrow.  
!-----------------------------------------------------------------------
   use phys_control, only: use_mass_borrower, print_fixer_message 

   implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   character (len=*), intent(in) :: subnam ! name of calling routine

   integer, intent(in) :: idx          ! chunk/latitude index
   INTEGER, intent(in) :: ncol         ! number of atmospheric columns
   integer, intent(in) :: ncold        ! declared number of atmospheric columns
   integer, intent(in) :: lver         ! number of vertical levels in column
   integer, intent(in) :: lconst_beg   ! beginning constituent
   integer, intent(in) :: lconst_end   ! ending    constituent
   logical, intent(in) :: lfix         ! if true, fix negative tracers

   real(r8), intent(in) :: qmin(lconst_beg:lconst_end)      ! Global minimum constituent concentration

!
! Input/Output arguments
!
   real(r8), intent(inout) :: q(ncold,lver,lconst_beg:lconst_end) ! moisture/tracer field
   real(r8), intent(inout) :: qneg3mat(ncold,lver,lconst_beg:lconst_end) ! moisture/tracer field
!
!---------------------------Local workspace-----------------------------
!

   integer indx(ncol,lver)  ! array of indices of points < qmin
   integer nval(lver)       ! number of points < qmin for 1 level
   integer nvals            ! number of values found < qmin
   integer nn
   integer iwtmp
   integer i,ii,k           ! longitude, level indices
   integer m                ! constituent index
   integer iw,kw            ! i,k indices of worst violator

   logical found            ! true => at least 1 minimum violator found

   real(r8) worst           ! biggest violator

!-----------------------------------------------------------------------
!

   do m=lconst_beg,lconst_end
      nvals = 0
      found = .false.
      worst = 1.e35_r8
      iw = -1
!
! Test all field values for being less than minimum value. Set q = qmin
! for all such points. Trace offenders and identify worst one.
!
!DIR$ preferstream
      do k=1,lver
         nval(k) = 0
!DIR$ prefervector
         nn = 0
         do i=1,ncol
            if (q(i,k,m) < qmin(m)) then
               qneg3mat(i,k,m) = qneg3mat(i,k,m) + 1.0 ! Record qneg error for this constituent and time (AaronDonahue)
               nn = nn + 1
               indx(nn,k) = i
            end if
         end do
         nval(k) = nn
      end do

      do k=1,lver
         if (nval(k) > 0) then
            found = .true.
            nvals = nvals + nval(k)
            iwtmp = -1
!cdir nodep,altcode=loopcnt
            do ii=1,nval(k)
               i = indx(ii,k)
               if (q(i,k,m) < worst) then
                  worst = q(i,k,m)
                  iwtmp = ii
               end if
            end do
            if (iwtmp /= -1 ) kw = k
            if (iwtmp /= -1 ) iw = indx(iwtmp,k)
            if(lfix) then 
!cdir nodep,altcode=loopcnt
               do ii=1,nval(k)
                  i = indx(ii,k)
                  q(i,k,m) = qmin(m)
               end do
            end if
         end if
      end do
      if (lfix) then 
         if (found .and. abs(worst)>max(qmin(m),1.e-8_r8)) then 
            write(iulog,9001)subnam//'/'//trim(cnst_name(m)),m,idx,nvals,qmin(m),worst,get_rlon_p(idx,iw),get_rlat_p(idx,iw),kw
         end if
      else
         if (print_fixer_message .and. found .and. abs(worst)>max(qmin(m),1.e-8_r8)) then 
            write(iulog,8001)subnam//'/'//trim(cnst_name(m)),m,idx,nvals,worst,get_rlon_p(idx,iw),get_rlat_p(idx,iw),kw
         end if
      end if
   end do
!
   return
8001 format(' QNEG3 from ',a,':m=',i3,' lat/lchnk=',i7, &
            ' Min. mixing ratio violated at ',i4,' points. ', &
            ' Worst =',e8.1,' at lon,lat,k=',f8.4,',',f8.4,i3)
9001 format(' QNEG3 from ',a,':m=',i3,' lat/lchnk=',i7, &
            ' Min. mixing ratio violated at ',i4,' points.  Reset to ', &
            1p,e8.1,' Worst =',e8.1,' at lon,lat,k=',f8.4,',',f8.4,i3)
end subroutine qneg3

!====================================================================================
subroutine qneg4 (subnam  ,lchnk   ,ncol    ,ztodt   ,        &
                  qbot    ,srfrpdel,shflx   ,lhflx   ,qflx, qneg4mat    )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Check if moisture flux into the ground is exceeding the total
! moisture content of the lowest model layer (creating negative moisture
! values).  If so, then subtract the excess from the moisture and
! latent heat fluxes and add it to the sensible heat flux.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: J. Olson
! 
!-----------------------------------------------------------------------
   use ppgrid
   use physconst,    only: gravit, latvap
   use constituents, only: qmin, pcnst

   implicit none

!
! Input arguments
!
   character (len=*), intent(in) :: subnam         ! name of calling routine
!
   integer, intent(in) :: lchnk              ! chunk index
   integer, intent(in) :: ncol               ! number of atmospheric columns
!
   real(r8), intent(in) :: ztodt             ! two times model timestep (2 delta-t)
   real(r8), intent(in) :: qbot(pcols,pcnst) ! moisture at lowest model level
   real(r8), intent(in) :: srfrpdel(pcols)   ! 1./(pint(K+1)-pint(K))
!
! Input/Output arguments
!
   real(r8), intent(inout) :: shflx(pcols)   ! Surface sensible heat flux (J/m2/s)
   real(r8), intent(inout) :: lhflx(pcols)   ! Surface latent   heat flux (J/m2/s)
   real(r8), intent(inout) :: qflx (pcols,pcnst)   ! surface water flux (kg/m^2/s)
   real(r8), intent(inout) :: qneg4mat(pcols,qneg4_numcnst,qneg4_numflds)
!
!---------------------------Local workspace-----------------------------
!
   integer :: i,ii              ! longitude indices
   integer :: iw                ! i index of worst violator
   integer :: indxexc(pcols)    ! index array of points with excess flux
   integer :: nptsexc           ! number of points with excess flux
!
   real(r8):: worst             ! biggest violator
   real(r8):: excess(pcols)     ! Excess downward sfc latent heat flux
!
   integer qneg_idx,qneg_idc      ! index in qneg fields matrix (AaronDonahue)
   logical qneg_rec               ! Flag if qneg4 recorded for this field/constituent
!-----------------------------------------------------------------------
   qneg_rec = .false.
   call qneg_ind(subnam,1,4,qneg_idx,qneg_idc) ! Get index for this field (AaronDonahue)
   if (qneg_idx*qneg_idc > 0)  then
      qneg_rec = .true.
   end if
!
! Compute excess downward (negative) q flux compared to a theoretical
! maximum downward q flux.  The theoretical max is based upon the
! given moisture content of lowest level of the model atmosphere.
!
   nptsexc = 0
   do i = 1,ncol
      excess(i) = qflx(i,1) - (qmin(1) - qbot(i,1))/(ztodt*gravit*srfrpdel(i))
!
! If there is an excess downward (negative) q flux, then subtract
! excess from "qflx" and "lhflx" and add to "shflx".
!
      if (excess(i) < 0._r8) then
         if (qneg_rec) qneg4mat(i,qneg_idc,qneg_idx) = 1.0 ! Record qneg error for this constituent and time (AaronDonahue)
         nptsexc = nptsexc + 1
         indxexc(nptsexc) = i
         qflx (i,1) = qflx (i,1) - excess(i)
         lhflx(i) = lhflx(i) - excess(i)*latvap
         shflx(i) = shflx(i) + excess(i)*latvap
      end if
   end do
!
! Write out worst value if excess
!
   if (nptsexc.gt.0) then
      worst = 0._r8
      do ii=1,nptsexc
         i = indxexc(ii)
         if (excess(i) < worst) then
            worst = excess(i)
            iw = i
         end if
      end do
      write(iulog,9001) subnam,nptsexc,worst, lchnk, iw, get_rlat_p(lchnk,iw),get_rlon_p(lchnk,iw)
   end if
!
   return
9001 format(' QNEG4 WARNING from ',a8 &
            ,' Max possible LH flx exceeded at ',i5,' points. ' &
            ,', Worst excess = ',1pe12.4 &
            ,', lchnk = ',i5 &
            ,', i = ',i5 &
            ,', same as lat =', f8.4 &
            ,', lon =', f8.4 &
           )
end subroutine qneg4
! ===============================================================================
subroutine massborrow(subnam,lchnk,ncol,pcols,mbeg,mend,qmin,q,pdel,qneg3mat) 

!!....................................................................... 
!! The mass borrower borrows tracer mass from an adjacent layer. 
!! It conserves the mass and can avoid negative tracers. 
!! 
!! At level k, it will first borrow the mass from the layer k+1 (lower level). 
!! If the mass is not sufficient in layer k+1, it will borrow mass from 
!! layer k+2. The borrower will proceed this process until the bottom layer. 
!! If the tracer mass in the bottom layer goes negative, it will repeat the 
!! process from the bottom to the top. In this way, the borrower works for 
!! any shape of mass profiles.
!! 
!! The code is adapted from the tracer mass borrower implemented in the 
!! global aerosol-climate model ECHAM-HAM (Feichter et al.,1996; 
!! Stier et al., 2005, Zhang et al., 2012).
!! 
!! Author : Kai Zhang (kai.zhang@pnnl.gov) 
!!....................................................................... 

  use shr_kind_mod,    only: r8 => shr_kind_r8
  use ppgrid,          only: pver
  use spmd_utils,      only: masterproc
  use cam_logfile,     only: iulog
  use phys_control,    only: print_fixer_message

  implicit none

!! interface 
!!....................................................................... 

  character (len=*), intent(in) :: subnam                 ! name of calling routine
  integer, intent(in) :: lchnk                        ! chunk identifier
  integer, intent(in) :: ncol                         ! number of atmospheric columns
  integer, intent(in) :: pcols                        ! number of dim members
  integer, intent(in) :: mbeg                         ! first index 
  integer, intent(in) :: mend                         ! last index 
  real(r8), intent(in) :: qmin(mbeg:mend)             ! smallest value
  real(r8), intent(in) :: pdel(pcols,pver)            ! pressure thickness 
  real(r8), intent(inout) :: q(pcols,pver,mbeg:mend)  ! moisture/tracer field
  real(r8), intent(inout) :: qneg3mat(pcols,pver,mbeg:mend) ! moisture/tracer field

!! local 
!!....................................................................... 

  integer :: i, k, m, j 
  integer :: ic(pcols) 
  real(r8):: nmass, zeps
  real(r8):: bmass(pcols)

  integer :: qneg_idx, qneg_idc
  logical :: qneg_rec
  !! init
  !!....................................................................... 

  zeps = epsilon(1.0_r8)

  !! loop over tracers
  !!....................................................................... 

  do m = mbeg, mend
      qneg_rec = .false.
      call qneg_ind(subnam,m,3,qneg_idx,qneg_idc) ! Get index for this field (AaronDonahue)
      if (qneg_idx*qneg_idc > 0) then
          qneg_rec = .true.
      end if

     ic(1:ncol) = 0 

     bmass(1:ncol) = 0.0_r8
     
     !! top to bottom
     !!....................................................................... 

     do k = 1, pver
        do i = 1, ncol

           !! new mass in the current layer
           !!....................................................................... 

           nmass = q(i,k,m) + bmass(i)/pdel(i,k)

           if ( nmass > qmin(m) ) then

              !! if new mass in the current layer is positive, don't borrow mass any more 
              !!....................................................................... 

              q(i,k,m) = nmass
              bmass(i) = 0.0_r8

           else

              !! set mass to qmin in the current layer, and save bmass
              !!....................................................................... 

              bmass(i) = (nmass - qmin(m)) * pdel(i,k)

              q(i,k,m) = qmin(m) 
              if (qneg_rec) qneg3mat(i,k,m) = qneg3mat(i,k,m) + 1.0 ! Record qneg error for this constituent and time (AaronDonahue)

              ic(i) = ic(i) + 1 

           end if !! nmass > 0.0_r8 

        end do !! i 
     end do !! k 

!!     do i = 1, ncol
!!
!!        if(print_fixer_message .and. ic(i).gt.0) then 
!!            write(iulog,*) '### mass borrower T2B ### tracer : ', m, ' column : ', i, ' chunk : ', lchnk  
!!        end if 
!!
!!        if(print_fixer_message .and. bmass(i) < 0._r8 ) then 
!!            write(iulog,*) '### mass borrower B2T ### tracer : ', m, ' column : ', i, ' chunk : ', lchnk  
!!        end if 
!!
!!     end do 

     !!....................................................................... 
     !! bottom to top
     !!....................................................................... 
     
     do k = pver, 1, -1 

        do i = 1, ncol

           !! if the surface layer still needs to borrow mass 
           !!....................................................................... 

           if (bmass(i) < 0._r8 ) then

              !! new mass in the current layer
              !!....................................................................... 

              nmass = q(i,k,m) + bmass(i)/pdel(i,k)

              if ( nmass > qmin(m) ) then

                 !! if new mass in the current layer is positive, don't borrow mass any more 
                 !!....................................................................... 

                 q(i,k,m) = nmass 
                 bmass(i) = 0.0_r8

              else

                 !! if new mass in the current layer is negative, continue to borrow mass
                 !!....................................................................... 

                 bmass(i) = (nmass - qmin(m))*pdel(i,k)
                 q(i,k,m) = qmin(m)
                 if (qneg_idx > 0) qneg3mat(i,k,m) = qneg3mat(i,k,m) + 1.0 ! Record qneg error for this constituent and time (AaronDonahue)

              end if !! nmass > 0.0_r8 

           end if !! bmass(i) < -zeps 

        end do !! i 

     end do !! k 

!!!     do k = 1, pver
!!!     do i = 1, ncol
!!!        if(print_fixer_message .and. q(i,k,m).lt.qmin(m)) then 
!!!            write(iulog,*) '### massborrow ### index : ', m, ' column : ', i, ' level : ', k, ' chunk : ', lchnk, ' tracer conc : ', q(i,k,m) 
!!!        end if 
!!!     end do 
!!!     end do 

  end do !! m

  return 
  end subroutine massborrow
! ===============================================================================

subroutine qqflx_fixer (subnam  ,lchnk   ,ncol    ,ztodt   ,        &
                        q, rpdel,shflx   ,lhflx   ,qflx    )

!!........................................................................
!! Water conservation fixer 
!! 
!! If QFLX is too negative, the condensation or deposition water vapor at 
!! the surface will take all the available moisture in surface layer. This 
!! will cause problems in the vertical diffusion calculation, where QFLX 
!! is applied. In the original CESM model, QNEG4 is called to correct QFLX 
!! so that it won't take out all the available moisture from the surface layer.
!! 
!! The new fixer, named as qqflx fixer, borrows water vapor from the whole 
!! column above the surface layer proportionally and add moisture into the 
!! surface layer, so that it can compensate the downward (negative) QFLX. 
!!
!! The excess downward (negative) q flux is compared to a theoretical
!! maximum downward q flux.  The theoretical max is based upon the
!! given moisture content of lowest level of the model atmosphere.
!!
!! Author: Kai Zhang (kai.zhang@pnnl.gov) and Phil Rasch 
!!........................................................................

   use shr_kind_mod,    only: r8 => shr_kind_r8
   use ppgrid,          only: pcols, pver
   use phys_grid,       only: get_lat_p, get_lon_p, get_wght_all_p 
   use physconst,       only: gravit, latvap
   use constituents,    only: qmin, pcnst
   use cam_logfile,     only: iulog
   use spmd_utils,      only: masterproc
   use cam_abortutils,  only: endrun
   use phys_control,    only: print_fixer_message

   implicit none

   integer, parameter :: i_wv = 1 
   real(r8),parameter :: factor = 1.5_r8 

!! Input arguments
!!........................................................................

   character*8, intent(in) :: subnam         ! name of calling routine

   integer, intent(in) :: lchnk              ! chunk index
   integer, intent(in) :: ncol               ! number of atmospheric columns

   real(r8), intent(in) :: ztodt             ! time step
   real(r8), intent(in) :: rpdel(pcols,pver) ! 1./(pint(k+1)-pint(k))

!! Input/Output arguments
!!........................................................................

   real(r8), intent(in) :: shflx(pcols)        ! surface sensible heat flux (J/m2/s)
   real(r8), intent(in) :: lhflx(pcols)        ! surface latent   heat flux (J/m2/s)
   real(r8), intent(in) :: qflx (pcols,pcnst)  ! surface water flux (kg/m^2/s)
   real(r8), intent(inout) :: q(pcols,pver,pcnst) ! moisture (kg/kg)

!! Local workspace
!!........................................................................


   real(r8):: excess(pcols)     ! excess downward surface latent heat flux 

   integer :: i, k 
   real(r8) :: wst(pver)        ! local
   real(r8) :: wsp              ! local
   real(r8) :: wpnet            ! local 
   real(r8) :: qori             ! local 
   real(r8) :: qbot             ! local 
   real(r8) :: dqbot            ! local 
   real(r8) :: ratio            ! local 
   real(r8) :: dwpbot           ! local 
   real(r8) :: qflx_org         ! local 

!! begin
!!........................................................................


   !! loop over cols
   !!....................................................................... 

   do i = 1, ncol

      !! check if downward water flux is too large
      !!....................................................................... 

      excess(i) = qflx(i,i_wv) - (qmin(i_wv) - q(i,pver,i_wv))/(ztodt*gravit*rpdel(i,pver))

      if (excess(i) < 0.0_r8 ) then

         qori = q(i,pver,i_wv) 

         !! new q at the bottom to balance the negative qflx 
         !!....................................................................... 

         qbot = qmin(i_wv) - qflx(i,i_wv) * ztodt*gravit*rpdel(i,pver) 

         !! due to process splitting, it could happen that the estimated qbot is 
         !! not sufficiently large to compensate the negative qflx in the vertical 
         !! diffusion scheme, so an adjust factor can be applied. 
         !!....................................................................... 

         qbot = qbot * factor

         !! change of q at the bottom
         !!....................................................................... 

         dqbot = qbot - q(i,pver,i_wv) 

         !! change of q * pdel /g  
         !!....................................................................... 

         dwpbot = dqbot / (gravit*rpdel(i,pver))

         if(dwpbot.lt.0._r8) then 
            call endrun('qflx_fixer: dwpbot < 0 ')  
         end if 

         wsp = 0._r8 

         do k = 1, pver-1 
            wst(k) = q(i,k,i_wv) / (gravit*rpdel(i,k))
            wsp    = wsp + wst(k) 
         end do

         wpnet = wsp - dwpbot

         if(wpnet.gt.0._r8 .and. wsp.gt.0._r8) then 
          
            !! if there is sufficient water vapor in the column, scale q at each level
            !! to compensate the potential water vapor sink in the surface layer. 
            !!....................................................................... 
 
            ratio = wpnet / wsp 

            q(i,1:pver-1,i_wv) = q(i,1:pver-1,i_wv) * ratio  
             
            q(i,pver,i_wv) = qbot 

            write(iulog,*) ' ### qflx_fixer ### ', &
                           ' chunk', lchnk, &
                           ' col', i, &
                           ' lat = ', get_lat_p(lchnk,i), &
                           ' lon = ', get_lon_p(lchnk,i), & 
                           ' qflx = ', qflx(i,i_wv), &
                           ' qori = ', qori, &
                           ' qnew = ', qbot, &
                           ' column wp = ', wsp, & 
                           ' dwp by qflx = ', dwpbot, & 
                           ' q scaled to ', 100.*ratio, "%"

         else 

            write(iulog,*) ' ### qqflx_fixer ### ', &
                           ' chunk', lchnk, &
                           ' col', i, &
                           ' lat = ', get_lat_p(lchnk,i), &
                           ' lon = ', get_lon_p(lchnk,i), & 
                           ' original qflx = ', qflx_org, &  
                           ' qsurf = ', qori, & 
                           ' column wp = ', wsp, & 
                           ' dwp by qflx = ', dwpbot
            write(iulog,*) ' column does not have enough water to compensate water vapor sink due to negative qflx !!! '

            call endrun('qflx_fixer: spurious negative qflx ') 


!!            !! if the column above the surface layer doens't have enough water
!!            !! use qneg4 to fix qflx and printout the conservation error message
!!            !!....................................................................... 
!!          
!!            qflx_org   = qflx (i,i_wv) 
!! 
!!            qflx (i,i_wv) = qflx (i,i_wv) - excess(i)
!!            lhflx(i)   = lhflx(i) - excess(i)*latvap
!!            shflx(i)   = shflx(i) + excess(i)*latvap
!!
!!            write(iulog,*) ' ### qqflx_fixer ### ', &
!!                           ' chunk', lchnk, &
!!                           ' col', i, &
!!                           ' lat = ', get_lat_p(lchnk,i), &
!!                           ' lon = ', get_lon_p(lchnk,i), & 
!!                           ' original qflx = ', qflx_org, &  
!!                           ' modified qflx = ', qflx(i,i_wv), & 
!!                           ' qsurf = ', qori, & 
!!                           ' column wp = ', wsp, & 
!!                           ' dwp by qflx = ', dwpbot

         end if !! if(wpnet.gt.0._r8 .and. wsp.gt.0._r8) 
         
     end if !! if (excess(i) < 0.0_r8 ) 
 
   end do !! i 

!! end
!!........................................................................

   return
end subroutine qqflx_fixer




end module qneg
