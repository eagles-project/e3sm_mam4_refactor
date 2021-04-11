module conditional_diag_output_utils
!-------------------------------------------------
! Utility subroutines used for handling model output 
! for the conditional diagnostics. Currently included
! are a few small subroutines that construct output
! variable names, and a subroutine that makes
! addfld and add_default calls to register output
! variables during model initialization.
! The outfld calls can be found in module 
! conditional_diag_main, subroutine xxx.
!
! History:
!  First version by Hui Wan, PNNL, March-April 2021
!-------------------------------------------------
  use cam_abortutils, only: endrun

  use conditional_diag, only: cnd_diag_info_t

  implicit none
  public

contains

subroutine conditional_diag_output_init(pver, cnd_diag_info)
!----------------------------------------------------------------------- 
! 
! Purpose: Register variables related to conditional diagnostics for output
!
! Method: (1) Add variables to the master field list by doing addfld calls
!         (2) Add variables, by default, to the first history tape (h0 files)
!             by doing add_default calls.
!         These two things are done for each sampling condition. 
!         Registered output variables include, for each condition,
!          - the metric used for sampling,
!          - the flag resulting from conditional sampling,
!          - the various diagnostic variables to which the conditional sampling
!            is applied. Per user's choice, these diagnostics (and their increments
!            if requested) are written out after various physical processes.
!-----------------------------------------------------------------------
  use cam_history,         only: addfld, horiz_only, add_default
  use cam_history_support, only: max_fieldname_len

  integer,intent(in) :: pver
  type(cnd_diag_info_t), intent(in) :: cnd_diag_info

  integer          :: icnd, ifld, iphys, ii
  character(len=4) :: val_inc_suff(2), suff
  logical          :: l_output(2)

  character(len=max_fieldname_len) :: output_fld_name
  character(len=max_fieldname_len) :: output_fld_name2

  character(len=256) :: fld_long_name

  character(len=*),parameter :: subname = 'conditional_diag_output_init'

  if (cnd_diag_info%ncnd==0) return

  ! Loop through all sampling conditions. Each of then will have
  ! its out set of output variables identified by the prefix cndxx_

  do icnd = 1,cnd_diag_info%ncnd

     !-------------------------------------------------------
     ! Register the sampling metric and the flag
     !-------------------------------------------------------
     call get_metric_and_flag_names_for_output( icnd, cnd_diag_info, &!in
                                                output_fld_name,     &!out
                                                output_fld_name2    ) !out

     ! Add the 2 variables to the master list of possible output variables.
     ! The 4 arguments of an addfld call are: variable name, 
     ! vertical dimension name, avg. flag, units, long_name.
     ! Units are set to blank right now; we could add a namelist variable
     ! to let the user provide the info. 

     if (cnd_diag_info%metric_nver(icnd)==1) then

       call addfld(trim(output_fld_name ), horiz_only, 'A',' ','Metric used in conditional sampling')
       call addfld(trim(output_fld_name2), horiz_only, 'A',' ','Flags  used in conditional sampling')

     elseif(cnd_diag_info%metric_nver(icnd)==pver) then

       call addfld(trim(output_fld_name ), (/'lev'/),  'A',' ','Metric used in conditional sampling')
       call addfld(trim(output_fld_name2), (/'lev'/),  'A',' ','Flags  used in conditional sampling')

     elseif(cnd_diag_info%metric_nver(icnd)==pver+1) then

       call addfld(trim(output_fld_name ), (/'ilev'/), 'A',' ','Metric used in conditional sampling')
       call addfld(trim(output_fld_name2), (/'ilev'/), 'A',' ','Flags  used in conditional sampling')

     else 
       call endrun(subname//': invalid number of vertical levels for metric '//trim(cnd_diag_info%metric_name(icnd)))
     end if 

     ! Add the 2 variables to the first history tape (i.e., the h0 file)
     ! so that the user does not need to explicit request them via namelist.
     ! The 3 arguments of an add_default call are (1) variable name, (2) hist. tape index,
     ! and (3) hist. averaging flag

     call add_default(trim(output_fld_name ),1,' ')  
     call add_default(trim(output_fld_name2),1,' ')

     !-------------------------------------------------------------------
     ! Register the diagnostics fields and their increments 
     !-------------------------------------------------------------------
     ! Put the on/off switches for value and increment output into one
     ! array so that we can deal with both using the ii-loop below

     l_output = (/cnd_diag_info%l_output_state, cnd_diag_info%l_output_incrm/)

     ! In terms of variable names in the output files, the values and 
     ! increments of a field are distinguished by different suffixes

     val_inc_suff = (/"_VAL","_INC"/)

     do ii = 1,2 ! field value (state) or increment

        if (.not.l_output(ii)) cycle
        suff = val_inc_suff(ii)

        do ifld  = 1,cnd_diag_info%nfld
        do iphys = 1,cnd_diag_info%nphysproc

           call get_fld_name_for_output( suff, cnd_diag_info, &! in
                                         icnd, ifld, iphys,   &! in
                                         output_fld_name      )! out

           call get_fld_longname_for_output( suff, cnd_diag_info, &! in
                                             icnd, ifld, iphys,   &! in
                                             fld_long_name        )! out

           ! Add the variable to the master list of possible output variables.
           ! The 4 arguments of an addfld call are:
           ! variable name, vertical dimension name, avg. flag, units, long_name.
           ! Units are set to blank right now; we could add a namelist variable
           ! to let the user provide the info. 

           if (cnd_diag_info%fld_nver(ifld)==1) then
              call addfld(trim(output_fld_name), horiz_only, 'A',' ',trim(fld_long_name)) 

           elseif (cnd_diag_info%fld_nver(ifld)==pver) then
              call addfld(trim(output_fld_name), (/'lev'/),  'A',' ',trim(fld_long_name)) 

           elseif (cnd_diag_info%fld_nver(ifld)==pver+1) then
              call addfld(trim(output_fld_name), (/'ilev'/), 'A',' ',trim(fld_long_name)) 
           else
              call endrun(subname//': invalid number of vertical levels for '//cnd_diag_info%fld_name(ifld))
           end if

           ! Add the variable to the first history tape (i.e., the h0 file)
           ! so that the user does not need to explicit request it via namelist.
           ! The 3 arguments of an add_default call are (1) variable name, (2) hist. tape index,
           ! and (3) hist. averaging flag

           call add_default(trim(output_fld_name),1,' ')

        end do ! iphys
        end do ! ifld

     end do ! ii = 1,2, field value (state) or tendency

  end do ! icnd = 1,ncnd

end subroutine conditional_diag_output_init

!======================================================
function get_metric_and_flag_names_for_output( icnd, cnd_diag_info,   &! in
                                               metric_name_in_output, &! out
                                               flag_name_in_output    )! out

   use cam_history_support, only: max_fieldname_len

   integer,               intent(in)  :: icnd
   type(cnd_diag_info_t), intent(in)  :: cnd_diag_info

   character(len=max_fieldname_len) :: metric_name_in_output
   character(len=max_fieldname_len) :: flag_name_in_output

   character(len=2) :: icnd_str ! condition index as a string

   write(icnd_str,'(i2.2)') icnd
   metric_name_in_output = 'cnd'//icnd_str//'_'//trim(cnd_diag_info%metric_name(icnd))
     flag_name_in_output = 'cnd'//icnd_str//'_'//trim(cnd_diag_info%metric_name(icnd))//'_flag'

end function flag_name_in_output

!======================================================
subroutine get_fld_name_for_output( suff, cnd_diag_info,    &!in
                                    icnd, ifld, iphys,      &!in
                                    fld_name_in_output      )!out

   use cam_history_support, only: max_fieldname_len

   character(len=*),      intent(in)  :: suff
   type(cnd_diag_info_t), intent(in)  :: cnd_diag_info
   integer,               intent(in)  :: icnd, ifld, iphys

   character(len=max_fieldname_len),intent(out) :: fld_name_in_output 

   character(len=2) :: icnd_str ! condition index as a string

   write(icnd_str,'(i2.2)') icnd

   fld_name_in_output = 'cnd'//icnd_str//'_'// &
                        trim(cnd_diag_info%fld_name(ifld))//'_'// &
                        trim(cnd_diag_info%physproc_name(iphys))//suff

end subroutine get_fld_name_and_longname_for_output 

!======================================================
subroutine get_fld_longname_for_output( suff, cnd_diag_info,    &!in
                                        icnd, ifld, iphys,      &!in
                                        fld_long_name_in_output )!out

   use cam_history_support, only: max_fieldname_len

   character(len=*),      intent(in)  :: suff
   type(cnd_diag_info_t), intent(in)  :: cnd_diag_info
   integer,               intent(in)  :: icnd, ifld, iphys

   character(len=256),intent(out)     :: fld_long_name_in_output 

   character(len=2) :: icnd_str ! condition index as a string

   write(icnd_str,'(i2.2)') icnd

   fld_long_name_in_output = trim(cnd_diag_info%fld_name(ifld))//suff// &
                             ' at '//trim(cnd_diag_info%physproc_name(iphys))// &
                             ' sampled under condition '//icnd_str// &
                             ' ('//trim(cnd_diag_info%metric_name(icnd))//')' 

end subroutine get_fld_name_and_longname_for_output 

end module conditional_diag_output_utils
