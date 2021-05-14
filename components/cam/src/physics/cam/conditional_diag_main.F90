module conditional_diag_main

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use cam_abortutils, only: endrun

  use conditional_diag, only: cnd_diag_info, cnd_diag_info_t

  implicit none

  private

  public cnd_diag_checkpoint

  ! comparison types used in defining sampling condition

  integer, parameter :: GE  =  2
  integer, parameter :: GT  =  1
  integer, parameter :: EQ  =  0
  integer, parameter :: LT  = -1
  integer, parameter :: LE  = -2

  ! flags for grid cell/column masking

  real(r8),parameter :: ON  = 1._r8
  real(r8),parameter :: OFF = 0._r8

  ! fillvalue for cells that are masked out 

  real(r8),parameter :: FILLVALUE = 0._r8

contains

!======================================================
subroutine cnd_diag_checkpoint( diag, this_chkpt, state, pbuf, cam_in, cam_out )

  use time_manager,        only: get_nstep
  use ppgrid,              only: pcols
  use cam_history_support, only: max_fieldname_len
  use cam_history,         only: outfld

  use physics_types,    only: physics_state
  use camsrfexch,       only: cam_in_t, cam_out_t
  use physics_buffer,   only: physics_buffer_desc

  use conditional_diag,    only: cnd_diag_t
  use conditional_diag_output_utils, only: get_metric_and_flag_names_for_output, &
                                           get_fld_name_for_output

  type(cnd_diag_t),    intent(inout), target :: diag
  character(len=*),    intent(in)            :: this_chkpt

  type(physics_state),     intent(in) :: state
  type(physics_buffer_desc),pointer   :: pbuf(:)
  type(cam_in_t),          intent(in) :: cam_in
  type(cam_out_t),optional,intent(in) :: cam_out

  integer :: ncnd, nchkpt, nqoi
  integer :: icnd, ichkpt, ii, iqoi
  integer :: ncol, lchnk
  integer :: nstep
  integer :: mult_by_dp

  real(r8),pointer :: metric(:,:), flag(:,:), inc(:,:), old_x_dp(:,:)
  real(r8),allocatable :: new(:,:), new_x_dp(:,:)

  character(len=max_fieldname_len) :: outfldname, flag_name_out, metric_name_out

  !---------------------------------------------------------------------------
  if (cnd_diag_info%ncnd == 0) return  ! no conditional diagnostics requested 

  ncnd   = cnd_diag_info%ncnd
  nchkpt = cnd_diag_info%nchkpt
  nqoi   = cnd_diag_info%nqoi

  lchnk  = state%lchnk
  ncol   = state%ncol

  nstep  = get_nstep()  ! current time step. Later in this routine, we
                        ! - save QoI values starting from the first step.
                        ! - do increment calculation when nstep > 1.
                        ! - evaluate sampling metrics starting from the first step.
                        ! - do conditional sampling and outfld calls only when 
                        !   nstep > 2, because the sampling time window might 
                        !   involve some checkpints from the previous nstep,
                        !   and valid increments are available only from nstep = 2.

  !=======================================
  ! Obtain QoI values and/or increments
  !=======================================
  ! First determine if this checkpoint is active for QoI monitoring

  ichkpt = 0  ! 0 = checkpoint inactive; this is the default

  do ii = 1,nchkpt
     if ( trim(cnd_diag_info%qoi_chkpt(ii)) == trim(this_chkpt) ) then 
        ichkpt = ii
        exit
     end if
  end do

  if (ichkpt>0) then 
  !---------------------------------------------------------------------------
  ! This checkpoint is active for QoI monitoring. Obtain the QoI values 
  ! and/or their increments if needed, and save to variable "diag". 
  ! Note that here we only obtain and save the QoIs and/or increments.
  ! Conditional sampling won't be applied until the "cnd_end_chkpt" checkpoint
  !---------------------------------------------------------------------------
     do iqoi = 1,nqoi

        ! Allocate memory for tmp arrays. This has to be done inside the iqoi
        ! loop as different QoIs might have different numbers of vertical levels

        allocate( new     ( pcols,cnd_diag_info%qoi_nver(iqoi) ))
        allocate( new_x_dp( pcols,cnd_diag_info%qoi_nver(iqoi) ))

        !------------------------------------------------------------------------
        ! Obtain the most up-to-date values of the QoI (not yet multiplied by dp)
        !------------------------------------------------------------------------
        call get_values( new,                                &! out
                         trim(cnd_diag_info%qoi_name(iqoi)), &! in
                         state, pbuf, cam_in, cam_out )       ! in

        !----------------------------------------------------------------
        ! Save the QoI and its increment to corresponding
        ! components of the diag data structure. Because the QoI and 
        ! increment might need to be multiplied by dp (pressure layer 
        ! thickness) if requested by user, and the request might be 
        ! different for different sampling conditions, we have to 
        ! loop over different conditions
        !----------------------------------------------------------------
        do icnd = 1,ncnd

           x_dp = cnd_diag_info% x_dp(icnd,iqoi,ichkpt)

           !-----------------------------------------------------------------
           ! mutiply the new QoI values by dp and save in tmp array new_x_dp
           !-----------------------------------------------------------------
           call multiply_by_dp( new, state, x_dp, new_x_dp )! in, in, in, out 

           !---------------------------------------------------------------------
           ! save the QoI (possibly multiplied by dp) to the diag data structure
           !---------------------------------------------------------------------
           if (cnd_diag_info%l_output_state) then
              diag%cnd(icnd)%qoi(iqoi)% val(:,:,ichkpt) = new_x_dp(:,:) 
           end if

           !--------------------------------------------------------------
           ! calculate and save the inrements, possibly multiplied by dp 
           !--------------------------------------------------------------
           if (cnd_diag_info%l_output_incrm) then

              old_x_dp => diag%cnd(icnd)%qoi(iqoi)% old

              ! calculate increments

              if (nstep>1) then 
                 inc => diag%cnd(icnd)%qoi(iqoi)% inc(:,:,ichkpt)
                 inc(1:ncol,:) = new_x_dp(1:ncol,:) - old_x_dp(1:ncol,:)
              end if

              ! Save the current value of QoI as "old" for next checkpoint

              old_x_dp(1:ncol,:) = new_x_dp(1:ncol,:)

           end if ! l_output_incrm
           !----------------------

        end do ! icnd

        ! Calculations done for this QoI. Clean up.
        deallocate(new)
        deallocate(new_x_dp)

     end do ! iqoi = 1,nqoi
  end if    ! ichkpt > 0

  !=======================================================
  ! Evaluate sampling condition if this is cnd_eval_chkpt 
  !=======================================================
  do icnd = 1,ncnd

     ! Check if sampling condition needs to be evaluated at this checkpoint.
     ! The answer could be .t. for multiple icnd values

     if (trim(cnd_diag_info% cnd_eval_chkpt(icnd)) .eq. trim(this_chkpt)) then 

        !---------------------------------
        ! Get metric values and set flags 
        !---------------------------------
        metric => diag%cnd(icnd)%metric
        call get_values( metric,                                &! out
                         trim(cnd_diag_info%metric_name(icnd)), &! in
                         state, pbuf, cam_in, cam_out )          ! in

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

  !=======================================================================
  ! Apply conditional sampling, then send QoIs to history buffer.
  ! (Do this only when nstep > 2, because the sampling time window might 
  ! involve some checkpints from the previous nstep,
  ! and valid increments are available only from nstep = 2 onwards)
  !=======================================================================
  if (nstep > 2) then
  do icnd = 1,ncnd

     ! Check if conditional sampling needs to be completed at this checkpoint.
     ! The answer could be .t. for multiple icnd values

     if (trim(cnd_diag_info% cnd_end_chkpt(icnd)).eq.trim(this_chkpt)) then 

        ! Each sampling condition has its own flags that will be applied
        ! to all QoIs and checkpoints

        flag => diag%cnd(icnd)%flag

        !----------------------------------------------------------------
        ! Apply conditional sampling to QoIs and/or their increments,
        ! then do the outfld calls to send the values to history buffer. 
        ! In subroutine apply_masking,
        ! different actions are taken depending on the vertical 
        ! dimension sizes of the metric and the QoIs.
        !----------------------------------------------------------------
        ! Apply to QoI values

        if (cnd_diag_info%l_output_state) then        
           do iqoi = 1,nqoi

              call apply_masking( flag, diag%cnd(icnd)%qoi(iqoi)%val ) 
  
              do ichkpt = 1,nchkpt
                 call get_fld_name_for_output( '', cnd_diag_info, icnd, iqoi, ichkpt, outfldname)
                 call outfld( trim(outfldname), diag%cnd(icnd)%qoi(iqoi)%val(:,:,ichkpt), pcols, lchnk )
              end do

           end do
        end if

        ! QoI increments

        if (cnd_diag_info%l_output_incrm) then        
           do iqoi = 1,nqoi

              call apply_masking( flag, diag%cnd(icnd)%qoi(iqoi)%inc ) 

              do ichkpt = 1,nchkpt
                 call get_fld_name_for_output( '_inc', cnd_diag_info, icnd, iqoi, ichkpt, outfldname)
                 call outfld( trim(outfldname), diag%cnd(icnd)%qoi(iqoi)%inc(:,:,ichkpt), pcols, lchnk )
              end do

           end do
        end if

     end if  !trim(this_chkpt).eq.trim(cnd_diag_info% cnd_end_chkpt(icnd))
  end do ! icnd = 1,ncnd
  end if ! nstep > 2

end subroutine cnd_diag_checkpoint

!==================================================================
subroutine apply_masking( flag, array )

    real(r8),intent(in)    ::  flag(:,:)
    real(r8),intent(inout) :: array(:,:,:)

    integer :: kk             ! vertical level index for a loop
    integer :: flag_nver      ! # of vertical levels the flag array has
    integer :: array_nver     ! # of vertical levels the output array has
    integer :: nchkpt, ichkpt ! # of checkpoints to process, and the loop index
    integer :: pcols, icol    ! # of columns in grid chunk, and the loop index

         pcols = size(flag, 1)
     flag_nver = size(flag, 2) 
    array_nver = size(array,2) 
        nchkpt = size(array,3)

    if (flag_nver == array_nver) then 
    ! same vertical dimension size; simply apply masking - and do this 
    ! for all checkpoints

       do ichkpt = 1,nchkpt
          where(flag(:,:).eq.OFF) array(:,:,ichkpt) = FILLVALUE 
       end do

    elseif (flag_nver == 1 .and. array_nver > 1) then 
    ! apply the same masking to all vertical levels and checkpoints

       do icol = 1,pcols
          if (flag(icol,1).eq.OFF) array(icol,:,:) = FILLVALUE
       end do

    elseif (flag_nver > 1 .and. array_nver == 1) then
    ! if any level in a grid column is selected, select that column;
    ! In other words, mask out a column in the output array 
    ! only if cells on all levels in that column are masked out.

       do icol = 1,pcols
          if (all( flag(icol,:).eq.OFF )) then
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
subroutine get_values( arrayout, varname, state, pbuf, cam_in, cam_out )

  use physics_types,  only: physics_state
  use camsrfexch,     only: cam_in_t, cam_out_t
  use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
  use time_manager,   only: get_nstep
  use constituents,   only: cnst_get_ind

  real(r8),           intent(out)   :: arrayout(:,:)
  character(len=*),   intent(in)    :: varname
  type(physics_state),intent(in)    :: state
  type(physics_buffer_desc), pointer:: pbuf(:)
  type(cam_in_t),     intent(in)    :: cam_in
  type(cam_out_t),    intent(in)    :: cam_out

  real(r8),pointer :: ptr2d(:,:)
  real(r8),pointer :: ptr1d(:)

  character(len=*),parameter :: subname = 'conditional_diag_main:get_values'

  integer :: ncol, idx

  ncol = state%ncol

  !--------------------------------------------------------------------------------
  ! If the requested variable is one of the advected tracers, get it from state%q
  !--------------------------------------------------------------------------------
  ! cnst_get_ind returns the index of a tracer in the host  model's advected 
  ! tracer array. If variable is not found on the tracer list, and index value 
  ! of -1 will be returned

  call cnst_get_ind(trim(adjustl(varname)),idx, abort=.false.)  !in, out, in

  if (idx /= -1) then ! This variable is a tracer field
     arrayout(1:ncol,:) = state%q(1:ncol,:,idx)

  else
  !----------------------
  ! Non-tracer variables
  !----------------------
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

  !------ cam_in -------

  case('LANDFRAC')
     arrayout(1:ncol,1) = cam_in%landfrac(1:ncol)

  !------ cam_out -------

  case('FLDS')
     arrayout(1:ncol,1) = cam_out%flwds(1:ncol)

  !------ pbuf -------

  case('PBLH')
      idx = pbuf_get_index('pblh') ; call pbuf_get_field( pbuf, idx, ptr1d )
      arrayout(:,1) = ptr1d

  case('CLD')
      idx = pbuf_get_index('CLD')  ; call pbuf_get_field( pbuf, idx, ptr2d )
      arrayout(:,:) = ptr2d

 !elseif (varname.eq.'QSATW') then
 !   call qsatw()

 !elseif (varname.eq.'CAPE') then
 !   call cape()

  !-----------------------------------------------------------------------------------
  ! The following were added mostly for testing of the conditional diag functionality
  !-----------------------------------------------------------------------------------
  case('NSTEP')
     arrayout(1:ncol,:) = get_nstep()

  case('LAT')
     arrayout(1:ncol,1) = state%lat(1:ncol)

  case('LON')
     arrayout(1:ncol,1) = state%lon(1:ncol)

  case default 
     call endrun(subname//': unknow varname - '//trim(varname))
  end select

  end if !whether the requested variable is a tracer field

end subroutine get_values

!==============================================================
subroutine multiply_by_dp( arrayin, state, x_dp, arrayout )

  use physics_types,   only: physics_state
  use conditional_diag,only: PDEL,PDELDRY,NODP

  real(r8),           intent(out) :: arrayin (:,:)
  type(physics_state),intent(in)  :: state
  integer,            intent(in)  :: x_dp 
  real(r8),           intent(out) :: arrayout(:,:)

  integer :: ncol
  character(len=*),parameter :: subname = 'conditional_diag_main:multiply_by_dp'

  ncol = state%ncol

  select case (x_dp)
  case (PDEL)
     arrayout(1:ncol,:) = arrayin(1:ncol,:) * state%pdel(1:ncol,:)

  case (PDELDRY)
     arrayout(1:ncol,:) = arrayin(1:ncol,:) * state%pdeldry(1:ncol,:)

  case (NODP)
     arrayout(1:ncol,:) = arrayin(1:ncol,:)

  case default
     call endrun(subname//': unexpected value of x_dp')
  end select

end subroutine multiply_by_dp

!==============================================================
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
