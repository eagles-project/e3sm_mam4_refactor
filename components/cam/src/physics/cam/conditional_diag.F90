module conditional_diag
!-------------------------------------------------
! Conditional sampling and diagnostics.
! This module contains 
!  - derived types definitions
!  - namelist handling utilities
!
! History:
!  First version by Hui Wan, PNNL, March - Aprill 2021
!-------------------------------------------------
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use spmd_utils,     only: masterproc
  use cam_logfile,    only: iulog
  use cam_abortutils, only: endrun

  implicit none

  private

  ! Derived types

  public cnd_diag_info_t   ! for metadata
  public cnd_diag_t        ! for actual data

  ! Variable(s) of derived type

  public cnd_diag_info  ! metadata

  ! Subroutines

  public conditional_diag_readnl
  public conditional_diag_alloc
 !public conditional_diag_dealloc

  ! module parameters

  integer, parameter :: ncnd_max         = 10 ! max # of conditions allowed in a single simulation
  integer, parameter :: mname_maxlen     = 8  ! string length for metric name

  integer, parameter :: nqoi_max         = 20 ! max # of conditionally sampled QoIs in a single simulation
  integer, parameter :: qoiname_maxlen   = 8  ! string length for QoI name

  integer, parameter :: nchkpt_max       = 99 ! max # of active checkpoints in a single simulation
  integer, parameter :: chkptname_maxlen = 8  ! string length for checkpoint name

  !-------------------------------------------------------------------------------
  ! Derived type for metadata
  !-------------------------------------------------------------------------------
  type cnd_diag_info_t

    ! Do we want to write out the QoIs?
    logical :: l_output_state = .false.

    ! Do we want to write out increments of the QoIs? 
    logical :: l_output_incrm = .false.

    ! Which history tape (1-6) should contain a complete set of the conditional diagnostics?
    ! Note: we allow the user to specify multiple history tapes that will contain all output variables 
    ! associated with the conditional diagnostics. It is envisioned that these multiple history 
    ! tapes might have different output frequencies and averaging flags. In addition,
    ! the user can manually add or remove individual variables using fincl or fexcl, 
    ! just as is the case for all other variables on the master list.

    integer,allocatable :: hist_tape_with_all_output(:)  ! tape indices
    integer             :: ntape                         ! number of tapes

    ! Sampling conditions. 
    ! The current implementation allows the user to define multiple conditions,
    ! each of which will correspond to its own metric, QoIs, and output.
    ! But to keep it simple (at least as a start), we assume that
    ! the QoIs and checkpoints to monitor are the same for different sampling conditions

    integer                      :: ncnd = 0             ! total # of sampling conditions used in this simulation
    character(len=mname_maxlen),&
                     allocatable :: metric_name(:)       ! shape = (ncnd); name of the metric
    integer,allocatable          :: metric_nver(:)       ! shape = (ncnd); # of vertical levels of the metric
    integer,allocatable          :: metric_cmpr_type(:)  ! shape = (ncnd); see module parameters in conditional_diag_main.F90
    real(r8),allocatable         :: metric_threshold(:)  ! shape = (ncnd); threshold value for conditional sampling 
    real(r8),allocatable         :: metric_tolerance(:)  ! shape = (ncnd); tolerance for the "equal to" comparison type
    character(len=chkptname_maxlen),&
                     allocatable :: eval_after(:)        ! shape = (ncnd); checkpoints at which the evaluation of 
                                                         ! the sampling conditions will happen
    character(len=chkptname_maxlen),&
                     allocatable :: sample_after(:)      ! shape = (ncnd); checkpoints at which the sampling (masking) 
                                                         ! of QoIs will happen

    ! QoIs to be monitored. Each QoI can have 1, pver, or pver+1 vertical levels
    integer                                   :: nqoi = 0
    character(len=qoiname_maxlen),allocatable :: qoi_name(:)     ! shape = (nqoi)
    integer,allocatable                       :: qoi_nver(:)     ! shape = (nqoi); # of vertical levels

    ! Active checkpoints at which the QoI will be monitored 
    integer                                      :: nchkpt = 0     ! total # of active checkpoints
    character(len=chkptname_maxlen), allocatable :: chkpt_name(:)  ! checkpoint names 


  end type cnd_diag_info_t

  !-------------------------------------------------------------------------------
  ! Derived types for the sampling conditions and sampled fields
  !-------------------------------------------------------------------------------
  ! Values of a single sampled field at different checkpoints and the 
  ! inrements relative to the previous checkpoint

  type snapshots_and_increments_t

    real(r8), allocatable :: val(:,:,:) ! shape = (pcols,info%qoi_nver(iqoi),info%nchkpt) QoI values at active checkpoints
    real(r8), allocatable :: inc(:,:,:) ! shape = (pcols,info%qoi_nver(iqoi),info%nchkpt) QoI increments between adjacent active checkpoints
    real(r8), allocatable :: old(:,:)   ! shape = (pcols,info%qoi_nver(iqoi)) QoI values at the previous active checkpoint

  end type snapshots_and_increments_t

  !----------------------------------------------------------------------
  ! The collection of all QoIs sampled under the same condition, and
  ! the values of the metric and flag used for sampling

  type metric_and_qois_t

    real(r8),                        allocatable :: metric (:,:)     ! shape = (pcols, info%metric_nver(icnd))
    real(r8),                        allocatable :: flag   (:,:)     ! shape = (pcols, info%metric_nver(icnd))
    type(snapshots_and_increments_t),allocatable :: qoi    (:)       ! shape = (info%nqoi)

  end type metric_and_qois_t

  !----------------------------------------------------------------------
  ! A collection of multiple conditions (including the corresponding
  ! metrics and sampled QoIs)

  type cnd_diag_t

    type(metric_and_qois_t), allocatable :: cnd(:) ! shape = (info%ncnd)

  end type cnd_diag_t
  !-----------------------------------------------------------------------------


!===============================================================================
! Module variables
!===============================================================================
  type(cnd_diag_info_t) :: cnd_diag_info

contains
!===============================================================================
! Procedures
!===============================================================================
subroutine conditional_diag_readnl(nlfile)

   use cam_history_support,only: ptapes
   use infnan,          only: nan, assignment(=), isnan
   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'conditional_diag_readnl'

   ! Local variables for reading namelist

   character(len=mname_maxlen)     :: metric_name     (ncnd_max)
   integer                         :: metric_nver     (ncnd_max)
   integer                         :: metric_cmpr_type(ncnd_max)
   real(r8)                        :: metric_threshold(ncnd_max)
   real(r8)                        :: metric_tolerance(ncnd_max)
   character(len=chkptname_maxlen) :: eval_after      (ncnd_max)
   character(len=chkptname_maxlen) :: sample_after    (ncnd_max)

   character(len=chkptname_maxlen) :: chkpt_name(nchkpt_max)

   character(len=qoiname_maxlen)   :: qoi_name (nqoi_max)
   integer                         :: qoi_nver (nqoi_max)

   logical :: l_output_state, l_output_incrm
   integer :: hist_tape_with_all_output(ptapes)  ! tape indices
   integer :: ntape

   ! other misc local variables
   integer :: ncnd, nchkpt, nqoi, ntape

   integer :: ii

   !-------
   namelist /conditional_diag_nl/  &
            metric_name, metric_nver, metric_cmpr_type, metric_threshold, metric_tolerance, &
            eval_after, sample_after, & 
            chkpt_name, qoi_name, qoi_nver, &
            l_output_state, l_output_incrm, hist_tape_with_all_output

   !----------------------------------------
   !  Default values
   !----------------------------------------
   metric_name      = ' '
   metric_nver      = 0
   metric_cmpr_type = -99
   metric_threshold = nan
   metric_tolerance = 0._r8
   eval_after       = ' '
   sample_after     = ' '

   chkpt_name     = ' '

   qoi_name       = ' '
   qoi_nver       = 0 

   l_output_state = .false.
   l_output_incrm = .false.
   hist_tape_with_all_output = -1

   !----------------------------------------
   ! Read namelist and check validity
   !----------------------------------------
   if (masterproc) then

      ! Read namelist

      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'conditional_diag_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, conditional_diag_nl, iostat=ierr)
         if (ierr /= 0) then
            write(iulog,*) 'read error ',ierr
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)

      !------------------------------------------
      ! Count user-specified sampling conditions

      ii = 0
      do while ( (ii+1) <= ncnd_max .and. metric_name(ii+1) /= ' ')
         ii = ii + 1
      end do
      ncnd = ii

      !----------------------------------------------------------------------
      ! If no condition has been specified, set the other counts to zero
      !----------------------------------------------------------------------
      if (ncnd==0) then

         nqoi   = 0
         nchkpt = 0
         ntape  = 0

      !----------------------------------------------------------------------
      ! If at least one condition has been sepecified, do some sanity check, 
      ! then parse additional namelist settings
      !----------------------------------------------------------------------
      else

         if (any( metric_nver     (1:ncnd) <= 0   )) call endrun(subname//' error: need positive metric_nver for each metric_name')
         if (any( metric_cmpr_type(1:ncnd) == -99 )) call endrun(subname//' error: need valid metric_cmpr_type for each metric_name')
         if (any( isnan(metric_threshold(1:ncnd)) )) call endrun(subname//' error: need valid metric_threshold for each metric_name')
         if (any( eval_after      (1:ncnd) == ' ' )) call endrun(subname//' error: be sure to specify eval_after for each metric_name')
         if (any( sample_after    (1:ncnd) == ' ' )) call endrun(subname//' error: be sure to specify sample_after for each metric_name')

         !-------------------------------------------------------
         ! Count QoIs to be monitored, then do some sanity check

         ii = 0
         do while ( (ii+1) <= nqoi_max .and. qoi_name(ii+1) /= ' ')
            ii = ii + 1
         end do
         nqoi = ii

         if (any(qoi_nver(1:nqoi)<=0)) call endrun(subname//'error: need positive qoi_nver for each qoi_name')

         !--------------------------
         ! Count active checkpoints 

         ii = 0
         do while ( (ii+1) <= nchkpt_max .and. chkpt_name(ii+1) /= ' ')
            ii = ii + 1
         end do
         nchkpt = ii

         if (nqoi==0) nchkpt = 0 ! If user did not specify any QoI, set nchkpt to 0 for consistency
         if (nchkpt==0) nqoi = 0 ! If user did not specify any checkpoint for QoI monitoring, set nqoi to 0 for consistency

         !---------------------------------------------------------------------------------
         ! Count history tapes that will each contain a full suite of the output variables

         ii = 0
         do while ( (ii+1) <= ptapes .and. hist_tape_with_all_output(ii+1) >= 0)
            ii = ii + 1
         end do
         ntape = ii

         ! If the user did not specify any tape, then at least add the output variables to h0

         if (ntape==0) then
            ntape = 1
            hist_tape_with_all_output(ntape) = 1
         end if

      end if ! ncnd = 0 or > 0

   end if ! masterproc
   !--------------------------------------

#ifdef SPMD
   call mpibcast(ncnd,   1, mpiint, 0, mpicom)
   call mpibcast(nqoi,   1, mpiint, 0, mpicom)
   call mpibcast(nchkpt, 1, mpiint, 0, mpicom)
   call mpibcast(ntape,  1, mpiint, 0, mpicom)
#endif

   cnd_diag_info%ncnd = ncnd

   if (ncnd==0) then

      if (masterproc) then
         write(iulog,*)'==========================================================='
         write(iulog,*)'       *** Conditional diagnostics NOT requested ***'
         write(iulog,*)'==========================================================='
      end if

      return
   end if

#ifdef SPMD
   !--------------------------------------
   ! Broadcast namelist variables
   !--------------------------------------
   call mpibcast(metric_name,      ncnd_max*len(metric_name(1)), mpichar, 0, mpicom)
   call mpibcast(metric_nver,      ncnd_max,                     mpiint,  0, mpicom)
   call mpibcast(metric_cmpr_type, ncnd_max,                     mpiint,  0, mpicom)
   call mpibcast(metric_threshold, ncnd_max,                     mpir8,   0, mpicom)
   call mpibcast(metric_tolerance, ncnd_max,                     mpir8,   0, mpicom)
   call mpibcast(eval_after,       ncnd_max*len(eval_after(1)),  mpichar, 0, mpicom)
   call mpibcast(sample_after,     ncnd_max*len(sample_after(1)),mpichar, 0, mpicom)

   call mpibcast(chkpt_name,     nchkpt_max*len(chkpt_name(1)), mpichar, 0, mpicom)

   call mpibcast(qoi_name,  nqoi_max*len(qoi_name(1)),  mpichar, 0, mpicom)
   call mpibcast(qoi_nver,  nqoi_max,                   mpiint,  0, mpicom)

   call mpibcast(l_output_state, 1, mpilog, 0, mpicom)
   call mpibcast(l_output_incrm, 1, mpilog, 0, mpicom)

   call mpibcast(hist_tape_with_all_output, ptapes, mpiint, 0, mpicom)
#endif

   !-------------------------------------------
   ! Pack information into to cnd_diag_info
   !-------------------------------------------
   ! Conditions for conditional sampling

   allocate( cnd_diag_info% metric_name(ncnd), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info% metric_name')
   do ii = 1,ncnd
      cnd_diag_info% metric_name(ii) = trim(adjustl(metric_name(ii)))
   end do

   allocate( cnd_diag_info% metric_nver(ncnd), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info% metric_nver')
   cnd_diag_info% metric_nver(1:ncnd) = metric_nver(1:ncnd)

   allocate( cnd_diag_info% metric_cmpr_type(ncnd), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info% metric_cmpr_type')
   cnd_diag_info% metric_cmpr_type(1:ncnd) = metric_cmpr_type(1:ncnd)

   allocate( cnd_diag_info% metric_threshold(ncnd), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info% metric_threshold')
   cnd_diag_info% metric_threshold(1:ncnd) = metric_threshold(1:ncnd)

   allocate( cnd_diag_info% metric_tolerance(ncnd), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info% metric_tolerance')
   cnd_diag_info% metric_tolerance(1:ncnd) = metric_tolerance(1:ncnd)

   allocate( cnd_diag_info% eval_after(ncnd), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info% eval_after')
   cnd_diag_info% eval_after(1:ncnd) = eval_after(1:ncnd)

   allocate( cnd_diag_info% sample_after(ncnd), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info% sample_after')
   cnd_diag_info% sample_after(1:ncnd) = sample_after(1:ncnd)

   ! QoIs 

   cnd_diag_info%nqoi = nqoi

   allocate( cnd_diag_info% qoi_name(nqoi), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info% qoi_name')
   do ii = 1,nqoi
      cnd_diag_info% qoi_name(ii) = trim(adjustl(qoi_name(ii)))
   end do

   allocate( cnd_diag_info% qoi_nver(nqoi), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info% qoi_nver')
   cnd_diag_info% qoi_nver(1:nqoi) = qoi_nver(1:nqoi)

   ! Active checkpoints at which the QoIs will be monitored

   cnd_diag_info% nchkpt = nchkpt

   allocate( cnd_diag_info% chkpt_name(nchkpt), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info% chkpt_name')
   do ii = 1,nchkpt
      cnd_diag_info% chkpt_name(ii) = trim(adjustl(chkpt_name(ii)))
   end do

   ! output to history tape(s)

   cnd_diag_info%l_output_state = l_output_state
   cnd_diag_info%l_output_incrm = l_output_incrm

   allocate( cnd_diag_info% hist_tape_with_all_output(ntape), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info% hist_tape_with_all_output')
   cnd_diag_info% hist_tape_with_all_output(1:ntape) = hist_tape_with_all_output(1:ntape)

   !-----------------------------------------------
   ! Send information to log file
   !-----------------------------------------------
   if (masterproc) then

      write(iulog,*)'==========================================================='
      write(iulog,*)'       *** Conditional diagnostics requested ***'
      write(iulog,*)'==========================================================='

      write(iulog,*)
      write(iulog,'(4x,2x,a20,a6,a12,4a20)')'metric','nlev','cmpr_type','threshold','tolerance','eval_after','sample_after'
      do ii = 1,cnd_diag_info%ncnd
         write(iulog,'(i4.3,2x,a20,i6,i12,2e20.10,2a20)') ii,                               &
                                                adjustr(cnd_diag_info%metric_name(ii)),     &
                                                        cnd_diag_info%metric_nver(ii),      &
                                                        cnd_diag_info%metric_cmpr_type(ii), &
                                                        cnd_diag_info%metric_threshold(ii), &
                                                        cnd_diag_info%metric_tolerance(ii), &
                                                adjustr(cnd_diag_info%eval_after(ii)),      &
                                                adjustr(cnd_diag_info%sample_after(ii))
      end do

      write(iulog,*)
      write(iulog,*)'--------------------------------------------------'
      write(iulog,'(4x,a20)') 'chkpt_name'
      do ii = 1,cnd_diag_info%nchkpt
         write(iulog,'(i4.3,a20)') ii, adjustr(cnd_diag_info%chkpt_name(ii))
      end do

      write(iulog,*)
      write(iulog,*)'--------------------------------------------------'
      write(iulog,'(4x,a20,a6)')'QoI_name','nlev'
      do ii = 1,cnd_diag_info%nqoi
         write(iulog,'(i4.3,a20,i6)') ii, adjustr(cnd_diag_info%qoi_name(ii)), cnd_diag_info%qoi_nver(ii)
      end do
      write(iulog,*)

      write(iulog,*)'--------------------------------------------------'
      write(iulog,*)' l_output_state = ',l_output_state
      write(iulog,*)' l_output_incrm = ',l_output_incrm
      write(iulog,*)
      write(iulog,'(a,i3)')' hist_tape_with_all_output = ',hist_tape_with_all_output
      write(iulog,*)'==========================================================='
      write(iulog,*)

  end if  ! masterproc

end subroutine conditional_diag_readnl

!===============================================================================
subroutine conditional_diag_alloc( phys_diag, begchunk, endchunk, pcols, cnd_diag_info )

  type(cnd_diag_t), pointer :: phys_diag(:)

  integer, intent(in) :: begchunk, endchunk
  integer, intent(in) :: pcols
  type(cnd_diag_info_t),intent(in) :: cnd_diag_info

  integer :: ierr, lchnk

  character(len=*), parameter :: subname = 'conditional_diag_alloc'
  !-----------------------------------------------------------------

  allocate(phys_diag(begchunk:endchunk), stat=ierr)
  if( ierr /= 0 ) then
     write(iulog,*) subname//': phys_diag allocation error = ',ierr
     call endrun(subname//': failed to allocate phys_diag array')
  end if

  if (cnd_diag_info%ncnd <= 0) return

  do lchnk=begchunk,endchunk
     call single_chunk_cnd_diag_alloc( phys_diag(lchnk), lchnk, pcols, cnd_diag_info )
  end do

  if (masterproc) then
     write(iulog,*)
     write(iulog,*) "====================================================================="
     write(iulog,*) " Finished memory allocation for conditional diagnostics including:"
     write(iulog,*) cnd_diag_info%ncnd,   " conditions"
     write(iulog,*) cnd_diag_info%nqoi,   " quantities of interest (QoIs)"
     write(iulog,*) cnd_diag_info%nchkpt, " active checkpoints to monitor QoIs"
     write(iulog,*) "====================================================================="
     write(iulog,*)
  end if

end subroutine conditional_diag_alloc


!-----------------------------------------------------------------------------
! Allocate memory for conditional diagnostics in a single grid chunk
!-----------------------------------------------------------------------------
subroutine single_chunk_cnd_diag_alloc( diag, lchnk, psetcols, cnd_diag_info )

  type(cnd_diag_t), intent(inout)   :: diag
  integer,intent(in)                :: lchnk

  integer, intent(in)               :: psetcols
  type(cnd_diag_info_t), intent(in) :: cnd_diag_info

  integer :: ierr = 0
  integer :: ncnd, icnd

  character(len=*), parameter :: subname = 'single_chunk_cnd_diag_alloc'

  !----------------------------------------------
  ! Allocate an array for all sampling conditions
  !----------------------------------------------
  ncnd = cnd_diag_info%ncnd

  allocate( diag%cnd(ncnd), stat=ierr)
  if ( ierr /= 0 ) call endrun(subname//': allocation of diag%cnd')

  !----------------------------------------------------------------
  ! Allocate memory for metrics and QoIs of each condition
  !----------------------------------------------------------------
  do icnd = 1,ncnd

     call metrics_and_qois_alloc( diag%cnd(icnd),                  &
                                  cnd_diag_info%metric_nver(icnd), &
                                  cnd_diag_info%nchkpt,            &
                                  cnd_diag_info%nqoi,              &
                                  cnd_diag_info%qoi_nver,          &
                                  cnd_diag_info%l_output_state,    &
                                  cnd_diag_info%l_output_incrm,    &
                                  psetcols                         )
  end do

end subroutine single_chunk_cnd_diag_alloc


!-----------------------------------------------------------------------------
! Allocate memory for metrics and diagnostics for a single sampling condition
!-----------------------------------------------------------------------------
subroutine metrics_and_qois_alloc( cnd, metric_nver, nchkpt, nqoi, qoi_nver, &
                                   l_output_state, l_output_incrm, psetcols )

 !use infnan, only : inf, assignment(=)

  type(metric_and_fields_t), intent(inout) :: cnd

  integer, intent(in) :: metric_nver, nchkpt
  integer, intent(in) :: nqoi
  integer, intent(in) :: qoi_nver(nqoi)
  logical, intent(in) :: l_output_state
  logical, intent(in) :: l_output_incrm
  integer, intent(in) :: psetcols

  integer :: iqoi
  integer :: ierr

  character(len=*), parameter :: subname = 'metrics_and_qois_alloc'

  ! metric and flag, which might have 1, pver, or pver+1 vertical levels

  allocate( cnd% metric(psetcols,metric_nver), stat=ierr)
  if ( ierr /= 0 ) call endrun(subname//': allocation of cnd% metric')

  allocate( cnd% flag  (psetcols,metric_nver), stat=ierr)
  if ( ierr /= 0 ) call endrun(subname//': allocation of cnd% flag')

  ! QoIs 

  if (nqoi > 0) then

     allocate( cnd% qoi( nqoi ), stat=ierr)
     if ( ierr /= 0 ) call endrun(subname//': allocation of cnd% qoi')

     ! snapshots of each QoI 
     if (l_output_state) then
      do iqoi = 1, nqoi

        allocate( cnd%qoi(iqoi)% val(psetcols,qoi_nver(iqoi),nchkpt), stat=ierr)
        if ( ierr /= 0 ) call endrun(subname//': allocation of cnd%qoi% val')

        cnd%qoi(iqoi)% val(:,:,:) = 0._r8

      end do !iqoi
     end if

     ! increments of each QoI
     if (l_output_incrm) then
      do iqoi = 1, nqoi

        allocate( cnd%qoi(iqoi)% old(psetcols,qoi_nver(iqoi)), stat=ierr)
        if ( ierr /= 0 ) call endrun(subname//': allocation of cnd%qoi% old')

        allocate( cnd%qoi(iqoi)% inc(psetcols,qoi_nver(iqoi),nchkpt), stat=ierr)
        if ( ierr /= 0 ) call endrun(subname//': allocation of cnd%qoi% inc')

        cnd%qoi(iqoi)% old(:,:)   = 0._r8
        cnd%qoi(iqoi)% inc(:,:,:) = 0._r8

      end do !iqoi
     end if

   end if

end subroutine metrics_and_qois_alloc

!subroutine conditional_diag_dealloc
!end subroutine conditional_diag_alloc


end module conditional_diag

