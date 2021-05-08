module conditional_diag_main

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use cam_abortutils, only: endrun

  use conditional_diag, only: cnd_diag_info, cnd_diag_info_t

  implicit none

  private

  public conditional_diag_cal_and_output

  integer, parameter :: GE  =  2
  integer, parameter :: GT  =  1
  integer, parameter :: EQ  =  0
  integer, parameter :: LT  = -1
  integer, parameter :: LE  = -2

  real(r8),parameter :: ON  = 1._r8
  real(r8),parameter :: OFF = 0._r8

  real(r8),parameter :: FILLVALUE = 0._r8

contains

!======================================================
subroutine conditional_diag_cal_and_output( diag, this_chkpt, state, pbuf, cam_in )

 !use time_manager,        only: get_nstep
  use ppgrid,              only: pcols
  use cam_history_support, only: max_fieldname_len
  use cam_history,         only: outfld

  use physics_types,    only: physics_state
  use camsrfexch,       only: cam_in_t
  use physics_buffer,   only: physics_buffer_desc

  use conditional_diag,    only: cnd_diag_t
  use conditional_diag_output_utils, only: get_metric_and_flag_names_for_output, &
                                          get_fld_name_for_output

  type(cnd_diag_t),    intent(inout), target :: diag
  character(len=*),    intent(in)            :: this_chkpt

  type(physics_state),    intent(in) :: state
  type(physics_buffer_desc), pointer :: pbuf(:)
  type(cam_in_t),         intent(in) :: cam_in

  integer :: ncnd, nchkpt, nqoi
  integer :: icnd, ichkpt, ii, iqoi
  integer :: ncol, lchnk
 !integer :: nstep

  real(r8),pointer :: metric(:,:), flag(:,:), inc(:,:), old(:,:)
  real(r8),allocatable :: new(:,:)

  character(len=max_fieldname_len) :: outfldname, flag_name_out, metric_name_out

  !---------------------------------------------------------------------------
  if (cnd_diag_info%ncnd == 0) return  ! no conditional diagnostics requested 

  ncnd   = cnd_diag_info%ncnd
  nchkpt = cnd_diag_info%nchkpt
  nqoi   = cnd_diag_info%nqoi

  lchnk  = state%lchnk
  ncol   = state%ncol

  !=======================================
  ! Obtain QoI values (and/or increments)
  !=======================================
  ! First check if this checkpoint is active for QoI monitoring

  ichkpt = 0  ! 0 = checkpoint inactive; this is the default

  do ii = 1,nchkpt
     if ( trim(cnd_diag_info%chkpt_name(ii)) == trim(this_chkpt) ) then 
        ichkpt = ii
        exit
     end if
  end do

  if (ichkpt>0) then 
  !---------------------------------------------------------------------------
  ! This checkpoint is active for QoI monitoring. Obtain the QoI values 
  ! (and their increments if needed), and save to variable "diag". 
  ! Note that here we only obtain and save the QoIs. 
  ! Conditional sampling won't be applied until the "cnd_end_chkpt" checkpoint
  !---------------------------------------------------------------------------
     do iqoi = 1,nqoi

        !----------------------------------------------------------------
        ! Obtain the most up-to-date values of the QoI 
        !----------------------------------------------------------------
        allocate( new( pcols,cnd_diag_info%qoi_nver(iqoi) )

        call get_values( new, trim(cnd_diag_info%qoi_name(iqoi)), &! inout, in
                         state, pbuf, cam_in )                     ! in

        !----------------------------------------------------------------
        ! The current implementation is such that the same set of 
        ! checkpoints and QoIs are monitored for all different sampling
        ! conditions. Now that we have obtain the QoI values and/or increments 
        ! at checkpoint iqoi, we save the same set of values under all   
        ! conditions. When conditional sampling is applied later,
        ! the different copies will likely be sampled differently.
        ! In the future, if we decide to allow for different QoIs and/or
        ! checkpoints under different sampling conditions, then
        ! the loops in this subroutine will need to be reworked.
        !---------------------------------------------------------------
        if (cnd_diag_info%l_output_state) then
           do icnd = 1,ncnd
              diag%cnd(icnd)%qoi(iqoi)% val(1:ncol,:,ichkpt) = new(1:ncol,:)
           end do
        end if

        !--------------------------------------------------------------
        ! Calculate increments if requested by user
        !--------------------------------------------------------------
        if (cnd_diag_info%l_output_incrm) then

           icnd = 1
           inc => diag%cnd(icnd)%qoi(iqoi)% inc(:,:,ichkpt)
           old => diag%cnd(icnd)%qoi(iqoi)% old

          !if (.not.is_firststep)& 
           inc(1:ncol,:) = new(1:ncol,:) - old(1:ncol,:)

           old(1:ncol,:) = new(1:ncol,:)

           ! Save increments for other sampling conditions; update "old" value

           do icnd = 2,ncnd
             !if (.not.is_firststep)& 
              diag%cnd(icnd)%qoi(iqoi)% inc(1:ncol,:,ichkpt) = inc(1:ncol,:)
              diag%cnd(icnd)%qoi(iqoi)% old(1:ncol,:)        = new(1:ncol,:)
           end do
          
        end if ! l_output_incrm

        ! Calculations done for this QoI. Clean up.
        deallocate( new )

     end do ! iqoi = 1,nqoi
  end if ! ichkpt > 0

  !===============================================================
  ! Evaluate sampling condition if this is the correct checkpoint 
  !===============================================================
  do icnd = 1,ncnd

     ! Check if sampling condition needs to be evaluated at this checkpoint.
     ! The answer could be .t. for multiple icnd values

     if (trim(cnd_diag_info% cnd_eval_chkpt(icnd)) .eq. trim(this_chkpt)) then 

        !---------------------------------
        ! Get metric values and set flags 
        !---------------------------------
        metric => diag%cnd(icnd)%metric
        call get_values( metric, trim(cnd_diag_info%metric_name(icnd)), &! inout, in
                         state, pbuf, cam_in )                           ! in

        flag => diag%cnd(icnd)%flag
        call get_flags( metric, icnd, ncol, cnd_diag_info, flag )

        !--------------------------------------
        ! Apply conditional sampling to metric
        !--------------------------------------
        where(flag.eq.OFF)  metric = FILLVALUE

        !----------------------------------------------------
        ! Send both metric and flag values to history buffer
        !----------------------------------------------------
        call get_metric_and_flag_names_for_output( icnd, cnd_diag_info, metric_name_out, flag_name_out )

        call outfld( trim(metric_name_out), metric, pcols, lchnk )
        call outfld( trim(flag_name_out),     flag, pcols, lchnk )

     end if !right chkpt
  end do    !icnd

  !-------------------------------------------------------------------------------
  ! Apply conditional sampling, then send QoIs to history buffer
  !-------------------------------------------------------------------------------
  do icnd = 1,ncnd

     ! Check if conditional sampling needs to be completed at this checkpoint.
     ! The answer could be .t. for multiple icnd values

     if (trim(cnd_diag_info% cnd_end_chkpt(icnd)).eq.trim(this_chkpt)) then 

        ! Each sampling condition has its own flags that will be applied
        ! to all QoIs and checkpoints

        flag => diag%cnd(icnd)%flag

        !----------------------------------------------------------------
        ! Apply conditional sampling to QoIs (and/or their increments),
        ! then do the outfld calls for output. In subroutine apply_masking,
        ! different actions are taken depending on the vertical 
        ! dimension sizes of the metric and the QoIs.
        !----------------------------------------------------------------
        ! Apply to QoI values

        if (cnd_diag_info%l_output_state) then        
           do iqoi = 1,nqoi

              call apply_masking( flag, diag%cnd(icnd)%qoi(iqoi)%val ) 
  
              do ii = 1,nchkpt
                 call get_fld_name_for_output( '', cnd_diag_info, icnd, iqoi, ii, outfldname)
                 call outfld( trim(outfldname), diag%cnd(icnd)%qoi(iqoi)%val(:,:,ii), pcols, lchnk )
              end do

           end do
        end if

        ! QoI increments

        if (cnd_diag_info%l_output_incrm) then        
           do iqoi = 1,nqoi

              call apply_masking( flag, diag%cnd(icnd)%qoi(iqoi)%inc ) 

              do ii = 1,nchkpt
                 call get_fld_name_for_output( '_inc', cnd_diag_info, icnd, iqoi, ii, outfldname)
                 call outfld( trim(outfldname), diag%cnd(icnd)%qoi(iqoi)%inc(:,:,ii), pcols, lchnk )
              end do

           end do
        end if

     end if  !trim(this_chkpt).eq.trim(cnd_diag_info% cnd_end_chkpt(icnd))
  end do ! icnd = 1,ncnd

end subroutine conditional_diag_cal_and_output

!==================================================================
subroutine apply_masking( flag, array )

    real(r8),intent(in)    ::  flag(:,:)
    real(r8),intent(inout) :: array(:,:,:)

    integer :: flag_nver      ! # of vertical levels the flag array has
    integer :: array_nver     ! # of vertical levels the output array has
    integer :: kk             ! vertical level index for a loop
    integer :: nchkpt, ichkpt ! # of physical processes to process, and the loop index
    integer :: pcols, icol

         pcols = size(flag, 1)
     flag_nver = size(flag, 2) 
    array_nver = size(array,2) 
        nchkpt = size(array,3)

    if (flag_nver == array_nver) then 
    ! same vertical dimension size; simply apply masking

       do ichkpt = 1,nchkpt
          where(flag(:,:).eq.OFF) array(:,:,ichkpt) = FILLVALUE 
       end do

    elseif (flag_nver == 1 .and. array_nver > 1) then 
    ! apply the same masking to all vertical levels

       do icol = 1,pcols
          if (flag(icol,1).eq.OFF) array(icol,:,:) = FILLVALUE
       end do

    elseif (flag_nver > 1 .and. array_nver == 1) then
    ! if any level in a grid column is selected, select that column;
    ! In other words, mask out a column in the output array 
    ! only if flags on all levels in that column are masked out.

       do icol = 1,pcols
          if (all(flag(icol,:).eq.OFF)) then
             array(icol,1,:) = FILLVALUE
          end if
       end do

    else
    ! QoI and flag both have multiple levels but the number of 
    ! vertical levels are different. Do not apply masking.

       continue

    end if

end subroutine apply_masking 


!========================================================
subroutine get_values( arrayout, varname, state, pbuf, cam_in )

  use time_manager,   only: get_nstep
  use physics_types,  only: physics_state
  use camsrfexch,     only: cam_in_t
  use physics_buffer, only: physics_buffer_desc

  real(r8),           intent(inout) :: arrayout(:,:)
  character(len=*),   intent(in)    :: varname
  type(physics_state),intent(in)    :: state
  type(physics_buffer_desc), pointer:: pbuf(:)
  type(cam_in_t),     intent(in)    :: cam_in

  character(len=*),parameter :: subname = 'conditional_diag_main:get_values'

  integer :: ncol

  ncol = state%ncol

  select case (trim(adjustl(varname)))
  case('T')
     arrayout(1:ncol,:) = state%t(1:ncol,:)

  case('U')
     arrayout(1:ncol,:) = state%u(1:ncol,:)

  case('V')
     arrayout(1:ncol,:) = state%v(1:ncol,:)

  case('OMEGA')
     arrayout(1:ncol,:) = state%omega(1:ncol,:)

  case('PMID')
     arrayout(1:ncol,:) = state%pmid(1:ncol,:)

  case('PINT')
     arrayout(1:ncol,:) = state%pint(1:ncol,:)

  case('ZM')
     arrayout(1:ncol,:) = state%zm(1:ncol,:)

  case('ZI')
     arrayout(1:ncol,:) = state%zi(1:ncol,:)

  case('PS')
     arrayout(1:ncol,1) = state%ps(1:ncol)

  case('PHIS')
     arrayout(1:ncol,1) = state%phis(1:ncol)

  case('LANDFRAC')
     arrayout(1:ncol,1) = cam_in%landfrac(1:ncol)

 !elseif (varname.eq.'QSATW') then
 !   call qsatw()

 !elseif (varname.eq.'CAPE') then
 !   call cape()

  case('NSTEP')
     arrayout(1:ncol,:) = get_nstep() 

  case default 
     call endrun(subname//': unknow varname - '//trim(varname))
  end select

end subroutine get_values

subroutine get_flags( metric, icnd, ncol, cnd_diag_info, flag )

  real(r8),              intent(in) :: metric(:,:)
  integer,               intent(in) :: icnd
  integer,               intent(in) :: ncol
  type(cnd_diag_info_t), intent(in) :: cnd_diag_info

  real(r8), intent(inout) :: flag(:,:)

  character(len=*),parameter :: subname = 'get_flags'

  flag(:,:) = OFF

  select case (cnd_diag_info%metric_cmpr_type(icnd))
  case (GT)

    where (metric(1:ncol,:).gt.cnd_diag_info%metric_threshold(icnd))
      flag(1:ncol,:) = ON
    elsewhere
      flag(1:ncol,:) = OFF
    end where

  case (GE)

    where (metric(1:ncol,:).ge.cnd_diag_info%metric_threshold(icnd))
      flag(1:ncol,:) = ON
    elsewhere
      flag(1:ncol,:) = OFF
    end where

  case (LT)

    where (metric(1:ncol,:).lt.cnd_diag_info%metric_threshold(icnd))
      flag(1:ncol,:) = ON
    elsewhere
      flag(1:ncol,:) = OFF
    end where

  case (LE)

    where (metric(1:ncol,:).le.cnd_diag_info%metric_threshold(icnd))
      flag(1:ncol,:) = ON
    elsewhere
      flag(1:ncol,:) = OFF
    end where

  case (EQ)

    where ( abs(metric(1:ncol,:)-cnd_diag_info%metric_threshold(icnd)) &
            .le. cnd_diag_info%metric_tolerance(icnd)                  )
      flag(1:ncol,:) = ON
    elsewhere
      flag(1:ncol,:) = OFF
    end where

  case default
    call endrun(subname//': unknown cnd_diag_info%metric_cmpr_type')
  end select

end subroutine get_flags

end module conditional_diag_main
