
module macrop_driver

  !-------------------------------------------------------------------------------------------------------
  ! Purpose:
  !
  ! Provides the CAM interface to the prognostic cloud macrophysics
  !
  ! Author: Andrew Gettelman, Cheryl Craig October 2010
  ! Origin: modified from stratiform.F90 elements 
  !    (Boville 2002, Coleman 2004, Park 2009, Kay 2010)
  !-------------------------------------------------------------------------------------------------------

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use spmd_utils,    only: masterproc
  use cam_abortutils,    only: endrun

  implicit none
  private
  save

  public :: macrop_driver_readnl
  
  logical, public :: do_cldice             ! .true., park macrophysics is prognosing cldice
  logical, public :: do_cldliq             ! .true., park macrophysics is prognosing cldliq
  logical, public :: do_detrain            ! .true., park macrophysics is detraining ice into stratiform

  ! ------------------------- !
  ! Private Module Parameters !
  ! ------------------------- !

  contains

  ! ===============================================================================
  subroutine macrop_driver_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Namelist variables
   logical  :: macro_park_do_cldice  = .true.   ! do_cldice = .true., park macrophysics is prognosing cldice
   logical  :: macro_park_do_cldliq  = .true.   ! do_cldliq = .true., park macrophysics is prognosing cldliq
   logical  :: macro_park_do_detrain = .true.   ! do_detrain = .true., park macrophysics is detraining ice into stratiform

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'macrop_driver_readnl'

   namelist /macro_park_nl/ macro_park_do_cldice, macro_park_do_cldliq, macro_park_do_detrain
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'macro_park_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, macro_park_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)

      ! set local variables

      do_cldice  = macro_park_do_cldice
      do_cldliq  = macro_park_do_cldliq
      do_detrain = macro_park_do_detrain

   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(do_cldice,             1, mpilog, 0, mpicom)
   call mpibcast(do_cldliq,             1, mpilog, 0, mpicom)
   call mpibcast(do_detrain,            1, mpilog, 0, mpicom)
#endif

end subroutine macrop_driver_readnl
end module macrop_driver
