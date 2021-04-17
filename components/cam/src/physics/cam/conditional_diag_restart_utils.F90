module conditional_diag_restart_utils

  implicit none
  private

  public :: cnd_diag_init_restart

contains
  !==================================================================================================
  subroutine cnd_diag_init_restart( dimids, hdimcnt, pver, pver_id, pverp_id,            &! in
                                    file, cnd_metric_desc,  cnd_flag_desc,               &! inout
                                    cnd_fld_val_desc, cnd_fld_old_desc, cnd_fld_inc_desc )! inout
  !------------------------------------------------------------------------------------------------
  ! Purpose: add variables to the restart file for conditional sampling and diagnostics.
  ! History: First version by Hui Wan, PNNL, 2021-04
  !------------------------------------------------------------------------------------------------

  use cam_abortutils,   only: endrun
  use conditional_diag, only: cnd_diag_info
  use pio,              only: file_desc_t, var_desc_t, pio_def_var, pio_double

  integer, intent(in) :: dimids(4)
  integer, intent(in) :: hdimcnt
  integer, intent(in) :: pver
  integer, intent(in) :: pver_id
  integer, intent(in) :: pverp_id

  type(file_desc_t), intent(inout) :: file

  type(var_desc_t),allocatable,intent(inout) :: cnd_metric_desc(:)
  type(var_desc_t),allocatable,intent(inout) :: cnd_flag_desc(:)
  type(var_desc_t),allocatable,intent(inout) :: cnd_fld_val_desc(:,:,:)
  type(var_desc_t),allocatable,intent(inout) :: cnd_fld_inc_desc(:,:,:)
  type(var_desc_t),allocatable,intent(inout) :: cnd_fld_old_desc(:,:)

  ! Local variables

  integer :: ierr, ndims
  integer :: dimids_local(4)
  character(len=256) :: pname  !variable name in restart file

  integer :: ncnd, nphys, nfld, icnd, iphys, ifld, nver
  character(len=*),parameter :: subname = 'cnd_diag_init_restart'

  if (cnd_diag_info%ncnd <= 0 ) return

  !-------------------------------------------------------------------------------
  ! Copy dimension IDs to a local array.
  !-------------------------------------------------------------------------------
  ! Because the various metrics, flags and fields in the conditional diagnostics 
  ! data structre can have different sizes in the vertical dimension,
  ! the last element of the local variable dimids_local might get different values
  ! during this subroutine. The use of a local variable here ensures that 
  ! the array dimids in the calling subroutine remains intact.
  !-------------------------------------------------------------------------------
  dimids_local(:) = dimids(:)

  !-------------------------------------------------------------------------------
  ! Allocate memeory for the variable description arrays
  !-------------------------------------------------------------------------------
  ! (Question: would it be better to allocate and deallocate at the beginning
  ! and end of each run instead of each time step?)

  ncnd  = cnd_diag_info%ncnd
  nphys = cnd_diag_info%nphysproc
  nfld  = cnd_diag_info%nfld

  allocate( cnd_metric_desc(ncnd) )
  allocate( cnd_flag_desc(ncnd) )

  if (nfld>0) then
     allocate( cnd_fld_old_desc(nfld,ncnd) )
     allocate( cnd_fld_val_desc(nphys,ncnd,nfld) )
     allocate( cnd_fld_inc_desc(nphys,ncnd,nfld) )
  end if

  !-----------------------------------------------------------
  ! Add the metrics and corresponding flags to restart file
  !-----------------------------------------------------------
  do icnd = 1,ncnd

     nver = cnd_diag_info%metric_nver(icnd)

     ! Dimension information

     if ( nver == 1) then
        ndims = hdimcnt

     else if ( nver == pver ) then
        ndims = hdimcnt+1
        dimids_local(ndims) = pver_id

     else if ( nver == pver+1 ) then
        ndims = hdimcnt+1
        dimids_local(ndims) = pverp_id

     else
        call endrun(subname//': check cnd_diag_info%metric_nver')
     end if

     ! Add the metric variable

     write(pname,'(a,i2.2,a)') 'cnd',icnd,'_metric'
     ierr = pio_def_var(File, trim(pname), pio_double, dimids_local(1:ndims), cnd_metric_desc(icnd))

     ! Add the flag variable

     write(pname,'(a,i2.2,a)') 'cnd',icnd,'_flag'
     ierr = pio_def_var(File, trim(pname), pio_double, dimids_local(1:ndims),   cnd_flag_desc(icnd))

  end do

  !-----------------------------------------------------------
  ! Add the conditionally sampled diagnostics to restart file
  !-----------------------------------------------------------
  do ifld = 1,nfld

     ! Dimension information

     nver = cnd_diag_info%fld_nver(ifld)

     if ( nver == 1) then
        ndims = hdimcnt

     else if ( nver == pver ) then
        ndims = hdimcnt+1
        dimids_local(ndims) = pver_id

     else if ( nver == pver+1 ) then
        ndims = hdimcnt+1
        dimids_local(ndims) = pverp_id

     else
        call endrun(subname//': check cnd_diag_info%fld_nver')
     end if

     ! Add to the restart file the variables containing field values after various physics processes

     if (cnd_diag_info%l_output_state) then
        do icnd = 1,ncnd
         do iphys = 1,nphys
            write(pname,'(3(a,i2.2))') 'cnd',icnd, '_fld',ifld, '_val',iphys
            ierr = pio_def_var(File, trim(pname), pio_double, dimids_local(1:ndims), cnd_fld_val_desc(iphys,icnd,ifld))
         end do
        end do
     end if

     ! Add to the restart file the variables containing increments associated with various physics processes

     if (cnd_diag_info%l_output_incrm) then

        do icnd = 1,ncnd
         do iphys = 1,nphys
            write(pname,'(3(a,i2.2))') 'cnd',icnd, '_fld',ifld, '_inc',iphys
            ierr = pio_def_var(File, trim(pname), pio_double, dimids_local(1:ndims), cnd_fld_inc_desc(iphys,icnd,ifld))
         end do
        end do

        ! Add to the restart file the variable containing the "old" value of the field 

        do icnd = 1,ncnd
           write(pname,'(2(a,i2.2),a)') 'cnd',icnd, '_fld',ifld, '_old'
           ierr = pio_def_var(File, trim(pname), pio_double, dimids_local(1:ndims), cnd_fld_old_desc(icnd,ifld))
        end do

     end if

  end do !ifld

  end subroutine cnd_diag_init_restart
  !=========================================================
 
end module conditional_diag_restart_utils
