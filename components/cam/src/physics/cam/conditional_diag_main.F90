module conditional_diag_main

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use cam_abortutils, only: endrun

  use physics_types,    only: physics_state
  use camsrfexch,       only: cam_in_t

  use conditional_diag, only: cnd_diag_info, cnd_diag_info_t

  implicit none

  private

  public conditional_diag_cal_and_output

  integer, parameter :: GT  =  1
  integer, parameter :: GE  =  2
  integer, parameter :: LT  = -1
  integer, parameter :: LE  = -2

  real(r8),parameter :: ON  = 1._r8
  real(r8),parameter :: OFF = 0._r8

  real(r8),parameter :: FILLVALUE = 0._r8

contains

!======================================================
subroutine conditional_diag_cal_and_output( diag, proc_name, state, cam_in )

  use ppgrid,              only: pcols
  use cam_history_support, only: max_fieldname_len
  use cam_history,         only: outfld

  use conditional_diag_output_util, only: get_metric_and_flag_names_for_output, &
                                          get_fld_name_for_output

  type(cnd_diag_t),    intent(inout), target :: diag
  character(len=*),    intent(in)            :: proc_name

  type(physics_state), intent(in)          :: state
  type(cam_in_t),      intent(in),optional :: cam_in

  integer :: ncnd, nphysproc, nfld
  integer :: icnd, iphys, ii, ifld
  integer :: ncol, lchnk

  real(r8),pointer :: metric(:,:), flag(:,:), new(:,:), inc(:,:), old(:,:)
  real(r8),pointer :: fld(:,:)

  character(len=max_fieldname_len) :: outfldname, flag_name_out, metric_name_out

  if (cnd_diag_info%ncnd == 0) return  ! no conditional diagnostics requested 

  ncnd      = cnd_diag_info%ncnd
  nphysproc = cnd_diag_info%nphysproc
  nfld      = cnd_diag_info%nfld

  lchnk    = state%lchnk
  ncol     = state%ncol

  !---------------------------------------------------------------------------
  ! Check if the atmospheric process proc_name is being monitored 
  !---------------------------------------------------------------------------
  iphys = 0
  do ii = 1,nphysproc
     if ( trim(cnd_diag_info%physproc_name(ii)) == trim(proc_name) ) then 
        iphys = ii
        exit
     end if
  end do

  if (iphys>0) then 
  !------------------------------------------------------------------------
  ! This atmospheric process is being monitored. Calculate the diagnostics
  ! (and their increments if needed) and save to diag. Note that here 
  ! we only calculate and save the diagnostics and/or increments. 
  ! Conditional sampling won't be applied until the metrics and flags 
  ! have been obtained.
  !------------------------------------------------------------------------
     do ifld = 1,nfld

        ! Note: The current implementation is such that the same set of 
        ! atmospheric processes and diagnostics are monitored for all 
        ! different metrics. Here we calculate the diagnostics and/or
        ! increments just once, for icnd = 1, and then copy the values
        ! for other conditions (icnd = 2,ncnd). When conditional sampling
        ! is applied, these different copies will likely be sampled differently.
        !--------------------------------------------------------------
        ! Obtain the most up-to-date values of the diagnostic variable 

        icnd = 1
        new => diag%cnd(icnd)%fld(ifld)% val(:,:,iphys)
        call get_values( trim(cnd_diag_info%fld_name(ifld)), state, new ) !in, in, out

        do icnd = 2,ncnd
           diag%cnd(icnd)%fld(ifld)% val(1:ncol,:,iphys) = new(1:ncol,:)
        end do

        !--------------------------------------------------------------
        ! Calculate increments if requested by user

        if (cnd_diag_info%l_output_incrm) then

           icnd = 1
           inc => diag%cnd(icnd)%fld(ifld)% inc(:,:,iphys)
           old => diag%cnd(icnd)%fld(ifld)% old

           inc(1:ncol,:) = new(1:ncol,:) - old(1:ncol,:)
           old(1:ncol,:) = new(1:ncol,:)

           ! Save increments for other metrics; update "old" value

           do icnd = 2,ncnd
              diag%cnd(icnd)%fld(ifld)% inc(1:ncol,:,iphys) = inc(1:ncol,:)
              diag%cnd(icnd)%fld(ifld)% old(1:ncol,:)       = new(1:ncol,:)
           end do
          
        end if ! l_output_incrm

     end do ! ifld = 1,nfld
  end if ! iphys > 0

  !-------------------------------------------------------------------------------
  ! Apply conditional sampling, then send values to output
  !-------------------------------------------------------------------------------
  do icnd = 1,ncnd

     ! Check if sampling condition needs to be evaluated at this point 
     ! for condition icnd. The answer could be .t. for multiple icnd values

     if (trim(proc_name).eq.trim(cnd_diag_info%sample_after(icnd))) then 

        !----------------------------------------
        ! Get metric values and set flags 

        metric => diag%cnd(icnd)%metric
        call get_values( trim(cnd_diag_info%metric_name(icnd)), state, metric )

        flag => diag%cnd(icnd)%flag
        call get_flags( metric, icnd, ncol, cnd_diag_info, flag )

        !----------------------------------------
        ! Apply conditional sampling to metric

        where(flag.eq.OFF)  metric = FILLVALUE

        !-------------------------------------------------
        ! Send both the metric values and flags to output

        call get_metric_and_flag_names_for_output( icnd, cnd_diag_info, metric_name_out, flag_name_out )

        call outfld( trim(metric_name_out), metric, pcols, lchnk )
        call outfld( trim(flag_name_out),     flag, pcols, lchnk )

        !----------------------------------------------------------------
        ! Apply conditional sampling to diagnostics and their increments
        ! caused by various atmospheric processes.
        ! Different actions are taken based on the vertical
        ! dimension sizes of the metric and the diagnostic fields.
        !----------------------------------------------------------------
        ! Diagnostic fields

        if (cnd_diag_info%l_output_state) then        
           do ifld = 1,nfld

              call apply_masking( flag, diag%cnd(icnd)%fld(ifld)%val ) 
  
              do ii = 1,nphysproc
                 call get_fld_name_for_output( '_val', cnd_diag_info, icnd, ifld, ii, outfldname)
                 call outfld( trim(outfldname), diag%cnd(icnd)%fld(ifld)%val(:,:,ii), pcols, lchnk )
              end do

           end do
        end if

        ! Increments

        if (cnd_diag_info%l_output_incrm) then        
           do ifld = 1,nfld

              call apply_masking( flag, diag%cnd(icnd)%fld(ifld)%inc ) 

              do ii = 1,nphysproc
                 call get_fld_name_for_output( '_inc', cnd_diag_info, icnd, ifld, ii, outfldname)
                 call outfld( trim(outfldname), diag%cnd(icnd)%fld(ifld)%inc(:,:,ii), pcols, lchnk )
              end do

           end do
        end if

     end if  !trim(proc_name).eq.trim(cnd_diag_info%sample_after(im))
  end do ! icnd = 1,ncnd

end subroutine conditional_diag_cal_and_output

!==================================================================
subroutine apply_masking( flag, array )

    real(r8),intent(in)    ::  flag(:,:)
    real(r8),intent(inout) :: array(:,:,:)

    integer :: flag_nver      ! # of vertical levels the flag array has
    integer :: array_nver     ! # of vertical levels the output array has
    integer :: kk             ! vertical level index for a loop
    integer :: nphysproc, iph ! # of physical processes to process, and the loop index
    
     flag_nver = size(flag, 2) 
    array_nver = size(array,2) 
     nphysproc = size(array,3)

    if (flag_nver == array_nver) then 
    ! same vertical dimension size; simply apply masking

       do iph = 1,nphysproc
          where(flag(:,:).eq.OFF) array(:,:,iph) = FILLVALUE 
       end do

    elseif (flag_nver == 1 .and. array_nver > 1) then 
    ! apply the same masking to all vertical levels

       do iph = 1,nphysproc
       do kk  = 1,array_nver
          where(flag(:,1).eq.OFF) array(:,kk,iph) = FILLVALUE 
       end do
       end do

    else
    ! flag has multiple levels and field has a different
    ! number of levels. Do not apply masking.

       continue

    end if

end subroutine apply_masking 


!========================================================
subroutine get_values( varname, state, arrayout )

  use time_manager, only: get_nstep

  character(len=*),   intent(in)    :: varname
  type(physics_state),intent(inout) :: state
  real(r8),           intent(inout) :: arrayout(:,:)

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

  case('PMID')
     arrayout(1:ncol,:) = state%pmid(1:ncol,:)

  case('PINT')
     arrayout(1:ncol,:) = state%pmid(1:ncol,:)

  case('PS')
     arrayout(1:ncol,1) = state%ps(1:ncol)

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

  case default
    call endrun(subname//': unknown cnd_diag_info%metric_cmpr_type')
  end select

end subroutine get_flags

end module conditional_diag_main
