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

  public cnd_diag_info_t
  public cnd_diag_t

  ! Variable(s) of derived type

  public cnd_diag_info

  ! Subroutines

  public conditional_diag_readnl
  public conditional_diag_alloc
 !public conditional_diag_dealloc

  ! module parameters

  integer, parameter :: ncnd_max = 10
  integer, parameter :: mname_maxlen = 8

  integer, parameter :: nqoi_max = 20
  integer, parameter :: qoiname_maxlen = 8 

  integer, parameter :: nphysproc_max   = 99
  integer, parameter :: physproc_name_maxlen = 8

  !-------------------------------------------------------------------------------
  ! Derived type for metadata
  !-------------------------------------------------------------------------------
  type cnd_diag_info_t

    ! Do we want to write out the field value after different physical processes?
    logical :: l_output_state = .false.

    ! Do we want to write out increments caused by different physicall processes? 
    logical :: l_output_incrm = .false.

    ! Metrics used for conditional sampling.
    ! The current implementation allows the user to define multiple metrics,
    ! each of which will correspond to its own sample and output.
    ! But to keep it simple (at least as a start), we assume that
    ! the physical processes and physical fields to monitor are the same
    ! for different metrics

    integer                      :: ncnd = 0          ! total # of metrics used in this simulation
    character(len=mname_maxlen),&
                     allocatable :: metric_name(:)       ! shape = (ncnd); name of the metric
    integer,allocatable          :: metric_nver(:)       ! shape = (ncnd); # of vertical levels
    real(r8),allocatable         :: metric_threshold(:)  ! shape = (ncnd); threshold value for conditional sampling 
    integer,allocatable          :: metric_cmpr_type(:)  ! shape = (ncnd); see module parameters
    character(len=physproc_name_maxlen),&
                     allocatable :: sample_after(:)      ! shape = (ncnd); after which atmospheric process
                                                         ! will conditional sampling be applied? The process names
                                                         ! need to match ptend%name)

    ! Physical processes to be monitored
    integer                          :: nphysproc = 0    ! total # of processes
    character(len=physproc_name_maxlen),&
                         allocatable :: physproc_name(:)  ! process labels 

    ! Physical fields to be monitored. Each field can have 1, nlev, or nlev+1 vertical levels

    integer                                   :: nqoi = 0
    character(len=qoiname_maxlen),allocatable :: qoi_name(:)     ! shape = (nqoi)
    integer,allocatable                       :: qoi_nver(:)     ! shape = (nqoi); # of vertical levels

  end type cnd_diag_info_t

  !-------------------------------------------------------------------------------
  ! Derived types for the sampling conditions and sampled fields
  !-------------------------------------------------------------------------------
  ! Values of a single sampled field at different checkpoints and the 
  ! inrements relative to the previous checkpoint

  type snapshots_and_increments_t

    real(r8), allocatable :: val(:,:,:) ! shape = (pcols,info%qoi_nver(iqoi),info%nphysproc) field values after different processes
    real(r8), allocatable :: inc(:,:,:) ! shape = (pcols,info%qoi_nver(iqoi),info%nphysproc) increments caused by different processes
    real(r8), allocatable :: old(:,:)   ! shape = (pcols,info%qoi_nver(iqoi)) old field values

  end type snapshots_and_increments_t

  !----------------------------------------------------------------------
  ! The collection of all fields sampled under the same condition, and
  ! values of the metric used for sampling

  type metric_and_fields_t

    real(r8),                        allocatable :: metric (:,:)     ! shape = (pcols, info%metric_nver(icnd))
    real(r8),                        allocatable :: flag   (:,:)     ! shape = (pcols, info%metric_nver(icnd))
    type(snapshots_and_increments_t),allocatable :: qoi    (:)       ! shape = (info%nqoi)

  end type metric_and_fields_t

  !----------------------------------------------------------------------
  ! A collection of multiple conditions (including the corresponding
  ! metrics and sampled fields

  type cnd_diag_t

    type(metric_and_fields_t), allocatable :: cnd(:) ! shape = (info%ncnd)

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

   use infnan,          only: nan, assignment(=), isnan
   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'conditional_diag_readnl'

   ! Local variables for reading namelist

   character(len=mname_maxlen)         :: metric_name     (ncnd_max)
   integer                             :: metric_nver     (ncnd_max)
   integer                             :: metric_cmpr_type(ncnd_max)
   real(r8)                            :: metric_threshold(ncnd_max)
   character(len=physproc_name_maxlen) :: sample_after    (ncnd_max)

   character(len=physproc_name_maxlen) :: physproc_name(nphysproc_max)

   character(len=qoiname_maxlen)    :: qoi_name (nqoi_max)
   integer                          :: qoi_nver (nqoi_max)

   logical :: l_output_state, l_output_incrm

   ! other misc local variables
   integer :: ncnd
   integer :: nphysproc
   integer :: nqoi

   integer :: ii

   !-------
   namelist /conditional_diag_nl/  &
            metric_name, metric_nver, metric_cmpr_type, metric_threshold, sample_after, &
            physproc_name, qoi_name, qoi_nver, l_output_state, l_output_incrm

   !----------------------------------------
   !  Default values
   !----------------------------------------
   metric_name      = ' '
   metric_nver      = 0
   metric_cmpr_type = 0
   metric_threshold = nan
   sample_after     = ' '

   physproc_name   = ' '

   qoi_name       = ' '
   qoi_nver       = 0 

   l_output_state = .false.
   l_output_incrm = .false.

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

      ! Check validity of namelist variables for user-specified metrics

      ii = 0
      do while ( (ii+1) <= ncnd_max .and. metric_name(ii+1) /= ' ')
         ii = ii + 1
      end do
      ncnd = ii

      if (any( metric_nver     (1:ncnd)<=0     )) call endrun(subname//' error: need positive metric_nver for each metric_name')
      if (any( metric_cmpr_type(1:ncnd)==0     )) call endrun(subname//' error: need valid metric_cmpr_type for each metric_name')
      if (any( isnan(metric_threshold(1:ncnd)) )) call endrun(subname//' error: need valid metric_threshold for each metric_name')
      if (any( sample_after    (1:ncnd)==' '   )) call endrun(subname//' error: be sure to specify sample_after for each metric_name')

      ! Count atmospheric processes to monitor

      ii = 0
      do while ( (ii+1) <= nphysproc_max .and. physproc_name(ii+1) /= ' ')
         ii = ii + 1
      end do
      nphysproc = ii

      ! Check validity of namelist variables for physical fields to monitor

      ii = 0
      do while ( (ii+1) <= nqoi_max .and. qoi_name(ii+1) /= ' ')
         ii = ii + 1
      end do
      nqoi = ii

      if (any(qoi_nver(1:nqoi)<=0)) call endrun(subname//'error: need positive qoi_nver for each qoi_name')

   end if ! masterproc
   !--------------------------------------

#ifdef SPMD
   call mpibcast(ncnd,     1, mpiint, 0, mpicom)
   call mpibcast(nphysproc,1, mpiint, 0, mpicom)
   call mpibcast(nqoi,     1, mpiint, 0, mpicom)
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
   call mpibcast(sample_after,     ncnd_max*len(sample_after(1)),mpichar, 0, mpicom)

   call mpibcast(physproc_name,  nphysproc_max*len(physproc_name(1)), mpichar, 0, mpicom)

   call mpibcast(qoi_name,  nqoi_max*len(qoi_name(1)),  mpichar, 0, mpicom)
   call mpibcast(qoi_nver,  nqoi_max,                   mpiint,  0, mpicom)

   call mpibcast(l_output_state, 1, mpilog, 0, mpicom)
   call mpibcast(l_output_incrm, 1, mpilog, 0, mpicom)
#endif

   !-------------------------------------------
   ! Pack information into to cnd_diag_info
   !-------------------------------------------
   ! Metrics for conditional sampling

   allocate( cnd_diag_info%metric_name(ncnd), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info%metric_name')
   do ii = 1,ncnd
      cnd_diag_info%metric_name(ii) = trim(adjustl(metric_name(ii)))
   end do

   allocate( cnd_diag_info%metric_nver(ncnd), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info%metric_nver')
   cnd_diag_info%metric_nver(1:ncnd) = metric_nver(1:ncnd)

   allocate( cnd_diag_info%metric_cmpr_type(ncnd), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info%metric_cmpr_type')
   cnd_diag_info%metric_cmpr_type(1:ncnd) = metric_cmpr_type(1:ncnd)

   allocate( cnd_diag_info%metric_threshold(ncnd), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info%metric_threshold')
   cnd_diag_info%metric_threshold(1:ncnd) = metric_threshold(1:ncnd)

   allocate( cnd_diag_info%sample_after(ncnd), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info%sample_after')
   cnd_diag_info%sample_after(1:ncnd) = sample_after(1:ncnd)

   ! Atmospheric processes to monitor 

   cnd_diag_info%nphysproc = nphysproc

   allocate( cnd_diag_info%physproc_name(nphysproc), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info%physproc_name')
   do ii = 1,nphysproc
      cnd_diag_info%physproc_name(ii) = trim(adjustl(physproc_name(ii)))
   end do

   ! snapshots and increments of physical fields

   cnd_diag_info%nqoi = nqoi

   allocate( cnd_diag_info%qoi_name(nqoi), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info%qoi_name')
   do ii = 1,nqoi
      cnd_diag_info%qoi_name(ii) = trim(adjustl(qoi_name(ii)))
   end do

   allocate( cnd_diag_info%qoi_nver(nqoi), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info%qoi_nver')
   cnd_diag_info%qoi_nver(1:nqoi) = qoi_nver(1:nqoi)

   ! output to history file(s)

   cnd_diag_info%l_output_state = l_output_state
   cnd_diag_info%l_output_incrm = l_output_incrm

   !-----------------------------------------------
   ! Send information to log file
   !-----------------------------------------------
   if (masterproc) then

      write(iulog,*)'==========================================================='
      write(iulog,*)'       *** Conditional diagnostics requested ***'
      write(iulog,*)'==========================================================='

      write(iulog,*)
      write(iulog,'(4x,2x,a20,a6,a12,a20,a20)')'metric','nlev','cmpr_type','threshold', 'sample_after'
      do ii = 1,cnd_diag_info%ncnd
         write(iulog,'(i4.3,2x,a20,i6,i12,e20.10,a20)') ii,                          &
                                                adjustr(cnd_diag_info%metric_name(ii)),     &
                                                        cnd_diag_info%metric_nver(ii),      &
                                                        cnd_diag_info%metric_cmpr_type(ii), &
                                                        cnd_diag_info%metric_threshold(ii), &
                                                adjustr(cnd_diag_info%sample_after(ii))
      end do

      write(iulog,*)
      write(iulog,*)'--------------------------------------------------'
      write(iulog,'(4x,a20)') 'physproc_name'
      do ii = 1,cnd_diag_info%nphysproc
         write(iulog,'(i4.3,a20)') ii, adjustr(cnd_diag_info%physproc_name(ii))
      end do

      write(iulog,*)
      write(iulog,*)'--------------------------------------------------'
      write(iulog,'(4x,a20,a6)')'field','nlev'
      do ii = 1,cnd_diag_info%nqoi
         write(iulog,'(i4.3,a20,i6)') ii, adjustr(cnd_diag_info%qoi_name(ii)), cnd_diag_info%qoi_nver(ii)
      end do
      write(iulog,*)

      write(iulog,*)'--------------------------------------------------'
      write(iulog,*)' l_output_state = ',l_output_state
      write(iulog,*)' l_output_incrm = ',l_output_incrm
      write(iulog,*)
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
     write(iulog,*) cnd_diag_info%ncnd,      " conditions"
     write(iulog,*) cnd_diag_info%nphysproc, " atmospheric processes"
     write(iulog,*) cnd_diag_info%nqoi,      " quantities of interest to sample and output"
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
  ! Allocate memory for metrics and diagnostics of each condition
  !----------------------------------------------------------------
  do icnd = 1,ncnd

     call metrics_and_fields_alloc( diag%cnd(icnd),                  &
                                    cnd_diag_info%metric_nver(icnd), &
                                    cnd_diag_info%nphysproc, &
                                    cnd_diag_info%nqoi,      &
                                    cnd_diag_info%qoi_nver,  &
                                    cnd_diag_info%l_output_state, &
                                    cnd_diag_info%l_output_incrm, &
                                    psetcols                 )
  end do

end subroutine single_chunk_cnd_diag_alloc


!-----------------------------------------------------------------------------
! Allocate memory for metrics and diagnostics for a single sampling condition
!-----------------------------------------------------------------------------
subroutine metrics_and_fields_alloc( cnd, metric_nver, nphysproc, nqoi, qoi_nver, &
                                     l_output_state, l_output_incrm, psetcols )

 !use infnan, only : inf, assignment(=)

  type(metric_and_fields_t), intent(inout) :: cnd

  integer, intent(in) :: metric_nver, nphysproc
  integer, intent(in) :: nqoi
  integer, intent(in) :: qoi_nver(nqoi)
  logical, intent(in) :: l_output_state
  logical, intent(in) :: l_output_incrm
  integer, intent(in) :: psetcols

  integer :: iqoi
  integer :: ierr

  character(len=*), parameter :: subname = 'metrics_and_fields_alloc'

  ! the metric fields, which might have 1, pver, or pver+1 vertical levels

  allocate( cnd% metric(psetcols,metric_nver), stat=ierr)
  if ( ierr /= 0 ) call endrun(subname//': allocation of cnd%metric')

  allocate( cnd% flag  (psetcols,metric_nver), stat=ierr)
  if ( ierr /= 0 ) call endrun(subname//': allocation of cnd%flag')

  ! diagnostical fields

  if (nqoi > 0) then

     allocate( cnd% qoi( nqoi ), stat=ierr)
     if ( ierr /= 0 ) call endrun(subname//': allocation of cnd%qoi')

     ! snapshots of each field
     if (l_output_state) then
      do iqoi = 1, nqoi

        allocate( cnd%qoi(iqoi)% val(psetcols,qoi_nver(iqoi),nphysproc), stat=ierr)
        if ( ierr /= 0 ) call endrun(subname//': allocation of cnd%qoi%val')

        cnd%qoi(iqoi)% val(:,:,:) = 0._r8

      end do !iqoi
     end if

     ! increments of each field
     if (l_output_incrm) then
      do iqoi = 1, nqoi

        allocate( cnd%qoi(iqoi)% old(psetcols,qoi_nver(iqoi)), stat=ierr)
        if ( ierr /= 0 ) call endrun(subname//': allocation of cnd%qoi%old')

        allocate( cnd%qoi(iqoi)% inc(psetcols,qoi_nver(iqoi),nphysproc), stat=ierr)
        if ( ierr /= 0 ) call endrun(subname//': allocation of cnd%qoi%inc')

        cnd%qoi(iqoi)% old(:,:)   = 0._r8
        cnd%qoi(iqoi)% inc(:,:,:) = 0._r8

      end do !iqoi
     end if

   end if

end subroutine metrics_and_fields_alloc

!subroutine conditional_diag_dealloc
!end subroutine conditional_diag_alloc


end module conditional_diag

