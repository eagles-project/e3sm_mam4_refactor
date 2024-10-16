#ifdef YAML_NDROP
  !<"lchnk" is needed for the following code to work,
  ! temporarily pass it along from upper level subroutines
  ! as y_lchnk and uncomment the following code:
  ! integer, intent(in) :: y_lchnk

  ! **OR** if this subroutine is called in a nested loop of columns and levels,
  ! we also might need column index (y_i or icol) and level
  ! index (y_k or klev) to be passed to this routine and
  ! uncomment the following code:
  ! integer, intent(in) :: y_i, y_k, y_lchnk

  ! Fields are output for all vertical levels, so just need column
  ! y_mm = -1 when explmix is called for droplet number mixing ratio
  ! otherwise, y_mm = the local (compressed) index of the aerosol mode and species combinations in the raercol arrays,
  ! which is given by mm in the dropmixnuc subroutine, and the mam_idx array in the ndrop_init subroutine

  integer, intent(in) :: y_i, y_lchnk, y_mm

  !>

  character(len=200) :: ext_str  ! this is used when multiple sets of yaml output is needed to cover different options (e.g., true and false)
  type(yaml_vars) :: yaml
  integer  :: unit_input, unit_output, y_nstep
  integer,save :: n_calls=0   ! some subroutines are called multiple times in one timestep, record the number of calls
  logical,save :: is_firstsubstep = .true.

! y_mm is set to -1 when explmix is first called in subloop, for qndrop; when called for aerosols it has the value of the local mam index (mm)
! if we are past the first invocation of explmix and y_mm is -1 again, we are no longer in the first pass of the subloop 
  if((y_mm<0).and.(n_calls>0)) then
     is_firstsubstep = .false.  
  endif

  !populate YAML structure
!  yaml%lev_print = <Add hardwired level here> !level (**remove these if generating data for a dependent subroutines**)
!!!  yaml%lev_print = 54 !level (**remove these if generating data for a dependent subroutines**)
!  yaml%nstep_print = <add hardwired time step here> !time step(**remove these if generating data for a dependent subroutines**)
  
  ! Below used with phys_debug_lat = 39.0553303768561, phys_debug_lon = 262.904774388088
  yaml%nstep_print = 1417 !time step(**remove these if generating data for a dependent subroutines**)

  ! set ext_str if there are multiple sets of yaml output to be write out
  ! here gives an example that "flag" in the code can be 0, 1, or 2:
  ! write(ext_str,'(I2)') flag
  if(y_mm<0) then
     ext_str='qndrop'
  else
     write(ext_str,'(I2)') y_mm
     if(is_unact) then
        ext_str = 'mm'//adjustl(ext_str)
     else
        ext_str = 'actmm'//adjustl(ext_str)
     endif
  endif 

  !YAML file input generation code- DO NOT PORT to C++
  !print all inputs one-by-one at column "yaml%col_print"
  yaml%col_print = icolprnt(y_lchnk) !column to write data(**remove these if generating data for a dependent subroutines**)
  y_nstep = get_nstep() !time step (**remove these if generating data for a dependent subroutines**)

  yaml%flag_print = .false. ! to write or not to write data (**remove these if generating data for a dependent subroutines**)
!  if(yaml%col_print == y_i .and. y_nstep==yaml%nstep_print .and. y_k == yaml%lev_print) then ! if this column exists in y_lchnk
  if(yaml%col_print == y_i .and. y_nstep==yaml%nstep_print .and. is_firstsubstep ) then ! if this column exists in y_lchnk
     !<
     !In the case of y_i or y_k are not passed as arguments:
     !if(yaml%col_print >0 .and. y_nstep==yaml%nstep_print) then

     !For generating data for a dependent subroutines where "yaml" derived type is already initialized:
     !if(yaml%flag_print) then
     !>

     ! record number of calls that can output yaml file. you only need to write one set of input/output
     n_calls = n_calls+1
  !   if (n_calls==1) then ! output the first call. change this if writing out other calls
     ! if ((n_calls==1 .and. flag==0) .or. (n_calls==3 .and. flag==1) .or. (n_calls==5 .and. flag==2)) then ! this is an example of writing multiple sets of YAML files

        yaml%lchnk_print = y_lchnk !(**remove these if generating data for a dependent subroutines**)
        yaml%flag_print  = .true.!(**remove these if generating data for a dependent subroutines**)

        !open I/O yaml files (it can have an extra optional argument to pass a unique string to differentiate file names)
        call open_files('explmix', &  !intent-in
             unit_input, unit_output, trim(ext_str)) !intent-out, with the use of ext_str
!             unit_input, unit_output) !intent-out


        !start by adding an input string
        call write_input_output_header(unit_input, unit_output,yaml%lchnk_print,yaml%col_print, &
             'explmix',yaml%nstep_print, yaml%lev_print)

        !< add code for writing data here>
        !call write_var(unit_input, unit_output, fld_name,field)!write a single variable
        call write_var(unit_input, unit_output,'dtmix',dtmix)   !write a single variable
        call write_var(unit_input, unit_output,'is_unact',is_unact)   !write a single variable
!  the following is a module variable:
        call write_var(unit_input, unit_output,'pver',pver)   !write a single variable
        !call write_1d_var(unit_input, unit_output, fld_name,dim,field) ! writes 1D variables of any dimension
        call write_1d_var(unit_input, unit_output, 'qold',pver,qold) ! writes 1D variables of any dimension
        call write_1d_var(unit_input, unit_output, 'src',pver,src) ! writes 1D variables of any dimension
        call write_1d_var(unit_input, unit_output, 'ekkp',pver,ekkp) ! writes 1D variables of any dimension
        call write_1d_var(unit_input, unit_output, 'ekkm',pver,ekkm) ! writes 1D variables of any dimension
        call write_1d_var(unit_input, unit_output, 'overlapp',pver,overlapp) ! writes 1D variables of any dimension
        call write_1d_var(unit_input, unit_output, 'overlapm',pver,overlapm) ! writes 1D variables of any dimension
        !call write_2d_var(unit_input, unit_output, fld_name, dim1, dim2, field) ! writes 2D variables of any dimension: field(dim1,dim2)

! optional variable
        if (is_unact) then
          call write_1d_var(unit_input, unit_output, 'qactold',pver,qactold) ! writes 1D variables of any dimension
        endif
        !writes aerosol mmr from state%q or q vector (cloud borne and interstitial)
        !"aer_num_only" is .ture. if printing aerosol num only
        !call write_aerosol_mmr_from_stateq(unit_input, unit_output, fld_name,field,aer_num_only)

        !close only the input file, not the output file
        close(unit_input)
        call freeunit(unit_input)
!     endif
  endif
#endif
