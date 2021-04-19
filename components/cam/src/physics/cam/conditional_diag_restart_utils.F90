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

  subroutine cnd_diag_write_restart( phys_diag, begchunk, endchunk,         &! in
                                     physgrid, file_hdimsize, file_nhdim,   &! in
                                     pcols, chunk_ncols, fillvalue,         &! in
                                     file, cnd_metric_desc,  cnd_flag_desc, &! inout
                                     cnd_fld_val_desc, cnd_fld_inc_desc,    &! inout
                                     cnd_fld_old_desc                       )! inout
  !------------------------------------------------------------------------------------------------
  ! Purpose: add variables to the restart file for conditional sampling and diagnostics.
  ! History: First version by Hui Wan, PNNL, 2021-04
  !------------------------------------------------------------------------------------------------

  use cam_abortutils,   only: endrun
  use conditional_diag, only: cnd_diag_info
  use pio,              only: file_desc_t, var_desc_t, pio_def_var, pio_double
  use cam_grid_support, only: cam_grid_get_decomp, cam_grid_write_dist_array

  type(cnd_diag_t),intent(in) :: phys_diag(begchunk:endchunk)
  integer, intent(in) :: begchunk, endchunk

  integer, intent(in) :: physgrid
  integer, intent(in) :: file_hdimidsize(*)  ! horizontal dimension sizes in restart file
  integer, intent(in) :: file_nhdims         ! number of horizontal dimensions in restart file (?)

  integer, intent(in) :: pcols
  integer, intent(in) :: chunk_ncols(begchunk,endchunk)
  real(r8),intent(in) :: fillvalue

  type(file_desc_t), intent(inout) :: file

  type(var_desc_t),allocatable,intent(inout) :: cnd_metric_desc(:)
  type(var_desc_t),allocatable,intent(inout) :: cnd_flag_desc(:)
  type(var_desc_t),allocatable,intent(inout) :: cnd_fld_val_desc(:,:,:)
  type(var_desc_t),allocatable,intent(inout) :: cnd_fld_inc_desc(:,:,:)
  type(var_desc_t),allocatable,intent(inout) :: cnd_fld_old_desc(:,:)

  ! Local variables

  integer :: file_dims(3)   ! dimension sizes in restart file, local variable
  integer :: arry_dims(3)   ! dimension sizes of array holding values to be written out, local variable

  integer :: ierr
  integer :: lchnk,ncol
  integer :: ncnd, nphys, nfld, icnd, iphys, ifld, nver
  character(len=*),parameter :: subname = 'cnd_write_init_restart'

  real(r8):: tmpfield_2d_1(pcols, begchunk:endchunk)
  real(r8):: tmpfield_2d_2(pcols, begchunk:endchunk)

  real(r8), allocatable :: tmpfield_3d_1(:,:,:)
  real(r8), allocatable :: tmpfield_3d_2(:,:,:)


  if (cnd_diag_info%ncnd <= 0 ) return

  !---------------------------------------------------------------
  ! Gather dimension info and save in local variables
  !---------------------------------------------------------------
  ncnd  = cnd_diag_info%ncnd
  nphys = cnd_diag_info%nphysproc
  nfld  = cnd_diag_info%nfld

  file_dims(1:file_nhdims) = file_hdimsize(1:file_nhdims)

  !-----------------
  ! metric and flag
  !-----------------
  do icnd = 1,ncnd

     nver = cnd_diag_info%metric_nver(icnd)

     if (nver==1) then

        !--------------------------------------------
        ! get iodesc needed by pio_write_darray calls

        arry_dims(1) = pcols
        arry_dims(2) = endchunk - begchunk + 1

        call cam_grid_get_decomp(physgrid, arry_dims(1:2), file_dims(1:file_nhdims), pio_double, iodesc) ! 4xin, out

        !----------------------------------------------------------------
        ! pack metric and flag values into tmp arrays and write them out

        tmpfield_2d_1 = fillvalue
        tmpfield_2d_2 = fillvalue

        do lchnk = begchunk, endchunk
           ncol = chunk_ncols(lchnk) 
           tmpfield_2d_1(:ncol,lchnk) = phys_diag(lchnk)%cnd(icnd)% metric(:ncol,1)
           tmpfield_2d_2(:ncol,lchnk) = phys_diag(lchnk)%cnd(icnd)%   flag(:ncol,1)
        end do

        call pio_write_darray(File, cnd_metric_desc(icnd), iodesc, tmpfield_2d_1, ierr)
        call pio_write_darray(File,   cnd_flag_desc(icnd), iodesc, tmpfield_2d_2, ierr)

     else ! nver > 1

        !----------------------------------------------------
        ! prepare input for cam_grid_write_dist_array calls

        arry_dims(1) = pcols
        arry_dims(2) = nver
        arry_dims(3) = endchunk - begchunk + 1

        file_dims(file_nhdims+1) = nver

        allocate( tmpfield_3d_1(pcols,nver,begchunk:endchunk) )
        allocate( tmpfield_3d_2(pcols,nver,begchunk:endchunk) )

        tmpfield_3d_1 = fillvalue
        tmpfield_3d_2 = fillvalue
        
        !----------------------------------------------------------------
        ! pack metric and flag values into tmp arrays and write them out

        do lchnk = begchunk, endchunk
           ncol = chunk_ncols(lchnk)
           tmpfield_3d_1(:ncol,:,lchnk) = phys_diag(lchnk)%cnd(icnd)% metric(:ncol,:)
           tmpfield_3d_2(:ncol,:,lchnk) = phys_diag(lchnk)%cnd(icnd)%   flag(:ncol,:)
        end do

        call cam_grid_write_dist_array(File, physgrid, arry_dims(1:3), file_dims(1:grid_nhdims+1), tmpfield_3d_1, cnd_metric_desc(icnd))
        call cam_grid_write_dist_array(File, physgrid, arry_dims(1:3), file_dims(1:grid_nhdims+1), tmpfield_3d_2,   cnd_flag_desc(icnd))

        deallocate(tmpfield_3d_1)
        deallocate(tmpfield_3d_2)

     end if 
  end do

  !-------------------------------------------------------------------------
  ! Conditionally sampled diagnostics (fields): "old", "val", and "inc"
  !-------------------------------------------------------------------------
  do ifld = 1,nfld

     nver = cnd_diag_info%fld_nver(ifld)

     if (nver==1) then

        tmpfield_2d_1 = fillvalue
        tmpfield_2d_2 = fillvalue

        !--------------------------------------------
        ! get iodesc needed by pio_write_darray calls
        !--------------------------------------------
        arry_dims(1) = pcols
        arry_dims(2) = endchunk - begchunk + 1
 
        call cam_grid_get_decomp(physgrid, arry_dims(1:2), file_dims(1:file_nhdims), pio_double, iodesc) ! 4xin, out

        !-------------------------------------------------------
        ! pack field values into tmp arrays and write them out
        !-------------------------------------------------------
        if (cnd_diag_info%l_output_state) then

           do icnd = 1,ncnd
           do iphys = 1,nphys

              do lchnk = begchunk, endchunk
                 ncol = chunk_ncols(lchnk)
                 tmpfield_2d_1(:ncol,lchnk) = phys_diag(lchnk)%cnd(icnd)%fld(ifld)% val(:ncol,1,iphys)
              end do

              call pio_write_darray(File, cnd_fld_val_desc(iphys,icnd,ifld), iodesc, tmpfield_2d_1, ierr)

           end do
           end do

        end if !cnd_diag_info%l_output_state

        !-------------------------------------------------------------------
        ! pack increment and old values into tmp arrays and write them out
        !-------------------------------------------------------------------
        if (cnd_diag_info%l_output_incrm) then

           do icnd = 1,ncnd

              ! increments corresponding to various physical processes

              do iphys = 1,nphys

                 do lchnk = begchunk, endchunk
                    ncol = chunk_ncols(lchnk)
                    tmpfield_2d_1(:ncol,lchnk) = phys_diag(lchnk)%cnd(icnd)%fld(ifld)% inc(:ncol,1,iphys)
                 end do

                 call pio_write_darray(File, cnd_fld_inc_desc(iphys,icnd,ifld), iodesc, tmpfield_2d_1, ierr)

              end do !iphys

              ! the "old" values
 
              do lchnk = begchunk, endchunk
                 ncol = chunk_ncols(lchnk)
                 tmpfield_2d_1(:ncol,lchnk) = phys_diag(lchnk)%cnd(icnd)%fld(ifld)%old(:ncol,1)
              end do

              call pio_write_darray(File, cnd_fld_old_desc(icnd,ifld), iodesc, tmpfield_2d_1, ierr)

           end do !icnd

        end if !cnd_diag_info%l_output_incrm


    else ! nver > 1

        !----------------------------------------------------
        ! prepare input for cam_grid_write_dist_array calls
        !----------------------------------------------------
        arry_dims(1) = pcols
        arry_dims(2) = nver
        arry_dims(3) = endchunk - begchunk + 1

        file_dims(file_nhdims+1) = nver

        allocate( tmpfield_3d_1(pcols,nver,begchunk:endchunk) )
        allocate( tmpfield_3d_2(pcols,nver,begchunk:endchunk) )
        
        tmpfield_3d_1 = fillvalue
        tmpfield_3d_2 = fillvalue
        
        !----------------------------------------------------------------
        ! pack metric and flag values into tmp arrays and write them out
        !----------------------------------------------------------------
        if (cnd_diag_info%l_output_state) then

           do icnd = 1,ncnd
           do iphys = 1,nphys

              do lchnk = begchunk, endchunk
                 ncol = chunk_ncols(lchnk)
                 tmpfield_3d_1(:ncol,:,lchnk) = phys_diag(lchnk)%cnd(icnd)%fld(ifld)% val(:ncol,:,iphys)
              end do

              call cam_grid_write_dist_array( File, physgrid, arry_dims(1:3), file_dims(1:grid_nhdims+1), &
                                              tmpfield_3d_1, cnd_fld_val_desc(iphy,icnd,ifld)             )

           end do
           end do

        end if !cnd_diag_info%l_output_state

        !----------------------------------------------------------------
        ! pack metric and flag values into tmp arrays and write them out
        !----------------------------------------------------------------
        if (cnd_diag_info%l_output_incrm) then

           do icnd = 1,ncnd
           do iphys = 1,nphys

              do lchnk = begchunk, endchunk
                 ncol = chunk_ncols(lchnk)
                 tmpfield_3d_1(:ncol,:,lchnk) = phys_diag(lchnk)%cnd(icnd)%fld(ifld)% inc(:ncol,:,iphys)
              end do

              call cam_grid_write_dist_array( File, physgrid, arry_dims(1:3), file_dims(1:grid_nhdims+1), &
                                              tmpfield_3d_1, cnd_fld_inc_desc(iphy,icnd,ifld)             )

           end do
           end do

           do icnd = 1,ncnd

              do lchnk = begchunk, endchunk
                 ncol = chunk_ncols(lchnk)
                 tmpfield_3d_1(:ncol,:,lchnk) = phys_diag(lchnk)%cnd(icnd)%fld(ifld)% old(:ncol,:)
              end do

              call cam_grid_write_dist_array( File, physgrid, arry_dims(1:3), file_dims(1:grid_nhdims+1), &
                                              tmpfield_3d_1, cnd_fld_old_desc(icnd,ifld)             )

           end do

       end if !cnd_diag_info%l_output_state


        deallocate(tmpfield_3d_1)
        deallocate(tmpfield_3d_2)

    end if 
  end do

  !-------------------------------------------------------------------------
  ! Done writing variables for restart. Dealocate description info arrays.
  ! (Question: would it be better to allocate and deallocate at the beginning
  ! and end of each run instead of each time step?)
  !-------------------------------------------------------------------------
  deallocate( cnd_metric_desc )
  deallocate( cnd_flag_desc )

  if (nfld>0) then
     deallocate( cnd_fld_old_desc )
     deallocate( cnd_fld_val_desc )
     deallocate( cnd_fld_inc_desc )
  end if
 
end module conditional_diag_restart_utils
