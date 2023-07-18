#ifdef YAML_AERO_DRYDEP
  if(yaml%flag_print) then ! if this column exists in lchnk

     !write output header
     call write_output_header(unit_output)

     !start writing data
     do icnst = 16,pcnst
        write(vname,'(A,I2)')'ptendq_',icnst
        call write_var(unit_output,vname,ptend%q(yaml%col_print,:,icnst))
        write(vname,'(A,I2)')'aerdepdryis_',icnst
        call write_var(unit_output,vname,aerdepdryis(yaml%col_print,icnst))
        write(vname,'(A,I2)')'aerdepdrycw_',icnst
        call write_var(unit_output,vname,aerdepdrycw(yaml%col_print,icnst))
        write(vname,'(A,I2)')'qqcw_',icnst
        call write_var(unit_output,trim(vname),qqcw(icnst)%fld(yaml%col_print,:))
     enddo

     !close the output file
     close(unit_output)
     call freeunit(unit_output)
  endif
#endif
