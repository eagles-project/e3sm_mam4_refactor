subroutine write_pkey_stuff(rank)

  implicit none
  integer, intent(in) :: rank
  integer :: pid, status
  integer, external :: getpid, hostnm
  character*32 :: hostname
  character*32 :: SLURM_NNODES, SLURM_NODEID, SLURM_JOB_ID, SLURM_LOCALID, SLURM_TASK_PID, SLURM_JOB_USER, SLURM_NTASKS
  character*32 :: SLURMD_NODENAME
  character*192 :: SLURM_JOB_NAME
  character*192 :: SLURM_JOB_NODELIST  !-- was looking for the compute nodename

  !#open (unit=94,file="pkey.txt",status='unknown')
  pid=getpid()
  status = hostnm(hostname)
  !write(*,'(a,i10,a,a)') "pkey: getpid=", pid, " hostname=", hostname


  call get_environment_variable("SLURM_NNODES", SLURM_NNODES)
  !write(*,*) "SLURM_NNODES=", trim(SLURM_NNODES)
  call get_environment_variable("SLURM_NTASKS", SLURM_NTASKS)
  !write(*,*) "SLURM_NTASKS=", trim(SLURM_NTASKS)
  call get_environment_variable("SLURM_NODEID", SLURM_NODEID)
  !write(*,*) "SLURM_NODEID=", trim(SLURM_NODEID)
  call get_environment_variable("SLURM_JOB_ID", SLURM_JOB_ID)
  !write(*,*) "SLURM_JOB_ID=", trim(SLURM_JOB_ID)
  call get_environment_variable("SLURM_JOB_NAME", SLURM_JOB_NAME)
  !write(*,*) "SLURM_JOB_NAME=", trim(SLURM_JOB_NAME)
  call get_environment_variable("SLURM_LOCALID", SLURM_LOCALID)
  !write(*,*) "SLURM_LOCALID=", trim(SLURM_LOCALID)
  call get_environment_variable("SLURM_TASK_PID", SLURM_TASK_PID)
  !write(*,*) "SLURM_TASK_PID=", trim(SLURM_TASK_PID)
  call get_environment_variable("SLURM_JOB_USER", SLURM_JOB_USER)
  !write(*,*) "SLURM_JOB_USER=", trim(SLURM_JOB_USER)
  call get_environment_variable("SLURMD_NODENAME", SLURMD_NODENAME)
  !write(*,*) "SLURMD_NODENAME=", trim(SLURMD_NODENAME)

  if (rank==0) then
     write(*,'(7(1x,a))') "pkey rank0:", trim(SLURM_NNODES), trim(SLURM_NTASKS), trim(SLURM_JOB_ID), trim(SLURM_JOB_USER), trim(SLURM_JOB_NAME)
  endif
  write(*,'(a,i10,7(1x,a))') "pkey:", pid, hostname, trim(SLURM_NODEID), trim(SLURM_LOCALID), trim(SLURM_TASK_PID), trim(SLURMD_NODENAME)


end subroutine write_pkey_stuff


program cime_driver

  !-------------------------------------------------------------------------------
  !
  ! Purpose: Main program for a CIME-driven model.  Can have different
  !          land, sea-ice, and ocean models plugged in at compile-time.
  !          These models can be either: stub, dead, data, or active
  !          components or some combination of the above.
  !
  !               stub -------- Do nothing.
  !               dead -------- Send analytic data back.
  !               data -------- Send data back interpolated from input files.
  !               active ------ Prognostically simulate the given component.
  !
  ! Method: Call appropriate initialization, run (time-stepping), and
  !         finalization routines.
  !
  !-------------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! share code & libs
  !----------------------------------------------------------------------------
  use shr_kind_mod,  only : r8 => SHR_KIND_R8
  use shr_kind_mod,  only : i8 => SHR_KIND_I8
  use shr_sys_mod,   only : shr_sys_irtc
  use perf_mod,      only : t_startf, t_adj_detailf, t_stopf, t_startstop_valsf
  use ESMF,          only : ESMF_Initialize, ESMF_Finalize
  use seq_comm_mct,  only : esmf_logfile_kind
  use cime_comp_mod, only : cime_pre_init1
  use cime_comp_mod, only : cime_pre_init2
  use cime_comp_mod, only : cime_init
  use cime_comp_mod, only : cime_run
  use cime_comp_mod, only : cime_final
  use mpi ! for call to  write_pkey_stuff(rank)

  implicit none

  !--------------------------------------------------------------------------
  ! timing variables
  !--------------------------------------------------------------------------
  integer(i8) :: beg_count, end_count, irtc_rate
  real(r8)    :: cime_pre_init1_time, ESMF_Initialize_time, &
       cime_pre_init2_time, cime_init_time_adjustment
  integer :: rank, ierr

  !--------------------------------------------------------------------------
  ! Setup and initialize the communications and logging.
  !--------------------------------------------------------------------------
  beg_count = shr_sys_irtc(irtc_rate)

  call cime_pre_init1()

  end_count = shr_sys_irtc(irtc_rate)
  cime_pre_init1_time = real( (end_count-beg_count), r8)/real(irtc_rate, r8)

  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call write_pkey_stuff(rank) ! ndk
  !print*, "ndk in cime_driver.F90"

  !--------------------------------------------------------------------------
  ! Initialize ESMF.  This is done outside of the ESMF_INTERFACE ifdef
  ! because it is needed for the time manager, even if the ESMF_INTERFACE
  ! is not used.
  !--------------------------------------------------------------------------
  beg_count = shr_sys_irtc(irtc_rate)

  call ESMF_Initialize(logkindflag=esmf_logfile_kind)

  end_count = shr_sys_irtc(irtc_rate)
  ESMF_Initialize_time = real( (end_count-beg_count), r8)/real(irtc_rate, r8)

  !--------------------------------------------------------------------------
  ! Read in the configuration information and initialize the time manager.
  !--------------------------------------------------------------------------
  ! Timer initialization has to be after determination of the maximum number
  ! of threads used across all components, so called inside of
  ! cime_pre_init2, as are t_startf and t_stopf for CPL:INIT and
  ! cime_pre_init2.
  !--------------------------------------------------------------------------
  beg_count = shr_sys_irtc(irtc_rate)

  call cime_pre_init2()

  end_count = shr_sys_irtc(irtc_rate)
  cime_pre_init2_time = real( (end_count-beg_count), r8)/real(irtc_rate, r8)

  !--------------------------------------------------------------------------
  ! Call the initialize, run and finalize routines.
  !--------------------------------------------------------------------------

  call t_startf('CPL:INIT')
  call t_adj_detailf(+1)

  call t_startstop_valsf('CPL:cime_pre_init1',  walltime=cime_pre_init1_time)
  call t_startstop_valsf('CPL:ESMF_Initialize', walltime=ESMF_Initialize_time)
  call t_startstop_valsf('CPL:cime_pre_init2',  walltime=cime_pre_init2_time)

  call cime_init()

  call t_adj_detailf(-1)
  call t_stopf('CPL:INIT')

  cime_init_time_adjustment = cime_pre_init1_time  &
       + ESMF_Initialize_time &
       + cime_pre_init2_time
  call t_startstop_valsf('CPL:INIT',  walltime=cime_init_time_adjustment, &
       callcount=0)

  call cime_run()
  call cime_final()

  !--------------------------------------------------------------------------
  ! Clean-up
  !--------------------------------------------------------------------------
  call ESMF_Finalize( )

end program cime_driver
