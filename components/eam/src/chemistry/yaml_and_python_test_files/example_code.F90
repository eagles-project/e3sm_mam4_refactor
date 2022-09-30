  ! This is an example code which shows how to use YAML file generation API
  ! Since all the code is temporary, we will stick with the practice of
  ! specifying minimal code required to accomplish the task.

  !------------------------------------------------------------------------------------------------------
  !make all function available using the "use" statement
  !Include the following statement at the top of the module
  !(No need to specify "only")
  use yaml_input_file_io
  !.....
  !.....
  !..... code
  !.....
  !.....


  !Let's take clddiag subroutine as an example:
  subroutine clddiag(t, pmid, pdel, cmfdqr, evapc, &
       cldt, cldcu, cldst, cme, evapr, &
       prain, cldv, cldvcu, cldvst, rain, &
       ncol, lchnk)

    ! Input arguments:
    real(r8), intent(in) :: t(pcols,pver)        ! temperature (K)
    real(r8), intent(in) :: pmid(pcols,pver)     ! pressure at layer midpoints
    real(r8), intent(in) :: pdel(pcols,pver)     ! pressure difference across layers
    real(r8), intent(in) :: cmfdqr(pcols,pver)   ! dq/dt due to convective rainout
    real(r8), intent(in) :: evapc(pcols,pver)    ! Evaporation rate of convective precipitation ( >= 0 )
    real(r8), intent(in) :: cldt(pcols,pver)    ! total cloud fraction
    real(r8), intent(in) :: cldcu(pcols,pver)    ! Cumulus cloud fraction
    real(r8), intent(in) :: cldst(pcols,pver)    ! Stratus cloud fraction
    real(r8), intent(in) :: cme(pcols,pver)      ! rate of cond-evap within the cloud
    real(r8), intent(in) :: evapr(pcols,pver)    ! rate of evaporation of falling precipitation (kg/kg/s)
    real(r8), intent(in) :: prain(pcols,pver)    ! rate of conversion of condensate to precipitation (kg/kg/s)
    integer, intent(in) :: ncol, lchnk

    ! Output arguments:
    real(r8), intent(out) :: cldv(pcols,pver)     ! fraction occupied by rain or cloud water
    real(r8), intent(out) :: cldvcu(pcols,pver)   ! Convective precipitation volume
    real(r8), intent(out) :: cldvst(pcols,pver)   ! Stratiform precipitation volume
    real(r8), intent(out) :: rain(pcols,pver)     ! mixing ratio of rain (kg/kg)

    !variables need for input/output test file
    character(len=2000):: finp, fout
    integer :: unit_input, unit_output

    !-----------------------------------------------------------------------------------------------------
    ! YAML and Python file generation code starts here
    ! Insert the following code at the start of the subroutine
    ! This example shows input variables of clddiags. Replace them with variables from
    ! respective subroutine
    !-----------------------------------------------------------------------------------------------------

    !To capture input
    if(icolprnt(lchnk) > 0) then ! if this column exists in lchnk

       !print all inputs one-by-one at column "i"
       i = icolprnt(lchnk) !column to write data

       !open I/O yaml and py files (provide a string same as the suborutine name)
       call open_files('subroutine_name', &  !intent-in
            unit_input, unit_output) !intent-out


       !start by adding an input header
       call write_input_header(unit_input, unit_output)

       !now add all the input variables
       call write_var_with_levs(unit_input, unit_output,'temperature',pver,t(i,:))
       !...Rest of the variables go here....

       !close just the input file, leave output file open
       close(unit_input)
       call freeunit(unit_input)
    endif
    !-----------------------------------------------------------------------------------------------------
    ! END - YAML file input generation code only
    !-----------------------------------------------------------------------------------------------------


    ! Subroutine's internal code sits here and does all the computations.
    ! After all the calculations are done, capture the output as shows in the code below
    !...
    !.... subrioutine code
    !....

    !-----------------------------------------------------------------------------------------------------
    ! Python module file output generation code
    !-----------------------------------------------------------------------------------------------------
   if(icolprnt(lchnk) > 0) then ! if this column exists in lchnk

      !write output header
      call write_output_header(unit_output)

      !print all outputs one-by-one at column "i"
      i = icolprnt(lchnk) !get column number

      !add all the output variables here
      call write_output_var_with_levs(unit_output,'cldv',   pver,cldv(i,:))
      !Rest of the output variable go here...

      !close the output file
      close(unit_output)
      call freeunit(unit_output)
   endif
   !------------------------------------------------------------------------------------------------------
   ! END - Python module file output generation code
   !------------------------------------------------------------------------------------------------------
 end subroutine clddiag
