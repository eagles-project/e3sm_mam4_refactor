module physics_update_mod 
  !============================================================================
  !Author: Balwinder Singh
  !
  !Purpose: An interface for physics update call:
  !For perturbation growth test, we need to output certain variables after every
  !physics_update call (when the state variables are updated). Placing outfld
  !directly in the existing physics_update call resulted in hard to resolve
  !circular dependencies. The alternative was to create this interface which calls 
  !physics_update and outfld subroutines one after another.

  !The circular dependencies, if we add outfld call to physics_types module, were:
  ! [NOTE: '<-' can be read as "depends upon" as the example below can be read as 
  !"physics_type depends upon cam_history" (due to the outfld call in physics_type)]

  !1. physics_type<-cam_history<-subcol_utils<-physics_types
  !2. physics_types<-cam_history<-chem_surfvals<-mo_flbc<-phys_gmean<-physics_type

  !This module helps breaking these circular dependencies but it requires changes 
  !to other parts of the code where we need to replace the "use" statement of
  !physics_update call with physics_update_mod module
  !============================================================================

  use spmd_utils,    only: masterproc
  use cam_abortutils,only: endrun
  use cam_history,   only: outfld, fieldname_len  
  use shr_kind_mod,  only: r8 => shr_kind_r8
  use ppgrid,        only: pcols, pver, begchunk
  use physics_types, only: physics_update_main, physics_ptend, physics_state, physics_tend
  use wv_saturation, only: qsat, qsat_water, qsat_ice, svp_ice

  implicit none
  private
  public  :: physics_update_init, physics_update, get_var_2d, get_var_3d
  
  save

  character(len = 25), parameter :: fname = 'pergro_ptend_names.txt'
  character(len = fieldname_len) :: plist(70)
  logical :: pergro_test_active
  integer :: unitn, pid

  !Following arrays and variables are declared so that we can add all variables in a loop to the history files for pergro test.
  !For adding any new variables, we need to do the following:
  !
  !1. Add variable name to 'hist_var3d' (and increment nvars_prtrb_hist variable accordingly), if a variable is part of the 
  !   constituent array ("q" array), add the _exact_ name as in cnst_add call(e.g.  NUMLIQ, CLDICE etc.)
  !2. If the variable is not present in the constituent array,add a "case" statement for that variable in the "select case" 
  !   construct in get_var_2d and get_var_3d functions in this module

  integer, public, parameter :: nvars_prtrb_hist = 15
 !character(len=6), public, parameter :: hist_var3d(nvars_prtrb_hist) = ['s     ', 't     ', 'Q     ', 'v     ', &
 !     'CLDLIQ', 'NUMLIQ', 'CLDICE', 'NUMICE', 'num_a1','num_a2','num_a3']
  character(len=6), public, parameter :: hist_var3d(nvars_prtrb_hist) = ['t     ', 'Q     ', 'u     ', 'CLDLIQ', 'NUMLIQ', &
        'CLDICE', 'NUMICE', 'RAINQM','NUMRAI','SNOWQM','NUMSNO','QSW   ','QSI   ', 'RHW   ', 'RHI   ']
  character(len=6), public, parameter :: hist_var2d(nvars_prtrb_hist) = ['tvint ', 'QVINT ', 'uvint ', 'CLVINT', 'NLVINT', &
        'CIVINT', 'NIVINT', 'RIVINT','NRVINT','SOVINT','NSVINT','QSWVIN','QSIVIN', 'RHWVIN', 'RHIVIN']

contains 

  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  subroutine physics_update_init()
    !Purpose: Initialize variables for physics update interface and pergro test

    use phys_control, only: phys_getopts
    use units,        only: getunit

    integer :: stat

    call phys_getopts(pergro_test_active_out = pergro_test_active)

    if (pergro_test_active) then
       !open file for writing ptend names
       pid = huge(1)
       if(masterproc) then
          pid = 0   !initialize pid for masterproc only
          unitn = getunit()
          open( unitn,file=fname,status='replace', form='formatted', action='write', position='append' )
          write(unitn,*,iostat=stat)'topphysbc' ! topphysbc record variables values at top of tphysbc
          if( stat > 0 ) then
             call endrun('Error writing topphysbc string  in '//fname//' file')
          endif
       endif
    endif

  end subroutine physics_update_init

  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------

  subroutine physics_update(state, ptend, dt, tend)
    !purpose: This subroutine calls physics_update_main (old physics_update)
    !and also output variables for pergro test

    use time_manager,  only: is_first_step, get_nstep

    
    !Arguments
    type(physics_ptend), intent(inout)  :: ptend   ! Parameterization tendencies
    type(physics_state), intent(inout)  :: state   ! Physics state variables
    real(r8),            intent(in)     :: dt      ! time step
    
    !optional arguments
    type(physics_tend ), intent(inout), optional  :: tend  ! Physics tendencies over timestep
    
    !Local vars
    character(len = fieldname_len)   :: pname, varname, vsuffix
    
    logical                          :: outfld_active, add_pname
    integer                          :: lchnk, stat, ip, ihist, ncol, nstep
    
    real(r8)                         :: timestep(pcols)  ! used for outfld call

    lchnk = state%lchnk
    ncol  = state%ncol
    
    !IMPORTANT:Store ptend%name as it will be modified in physics_update_main call
    pname = ptend%name
    
    !if nothing is to be updated in physics_update_main, DO NOT output (using outfld calls) 
    !PERGRO variables ("t_...", "s_..." etc.) below
    !Note: The following logical flag is required as sometimes "pname" is an empty string
    !      due to some stub routine calls (e.g. "iondrag_calc") where ptend is an 
    !      intent-out. This causes issues as intent-out will cause ptend%name to become undefined

    outfld_active = .true. !decides whether to call outfld calls or not
    if (.not. (any(ptend%lq(:)) .or. ptend%ls .or. ptend%lu .or. ptend%lv)) outfld_active = .false.
    
    !call the old physics update call
    call physics_update_main (state, ptend, dt, tend)
    
    if (pergro_test_active .and. outfld_active) then
       
       !write text file to be used for the post processing
       if(masterproc .and. lchnk == begchunk .and. is_first_step()) then
          !Here we write a text file to use for the post processing. We list
          !all the ptend names in this file. We do not want duplicates as it will
          !confuse the post processing script, therefore we skip ptend names when
          !they already exist in the plist

          !Find if this pname already exist in plist
          add_pname = .true.!decides whether to add pname to plist or not
          do ip = 1 , pid
             if (trim(adjustl(pname)) == trim(adjustl(plist(ip)))) then
                !already exists, do NOT add in plist
                add_pname = .false.
                exit
             endif
          enddo
          if (add_pname) then
             write(unitn,*,iostat=stat) pname
             if( stat > 0 ) then
                call endrun('Error writing '//pname// 'in '//fname//' file')
             endif
             !increment index pid and add pname to the list
             pid = pid + 1
             plist(pid) = trim(adjustl(pname))
          endif
       endif
          
       !call outfld
       do ihist = 1 , nvars_prtrb_hist
          vsuffix  = trim(adjustl(hist_var3d(ihist)))
          varname  = trim(adjustl(vsuffix))//'_'//trim(adjustl(pname)) ! form variable name
          !find the prognostic variable associated with this hist_var3d(ihist) via "get_var_3d" function
          call outfld( trim(adjustl(varname)), get_var_3d(state,vsuffix), pcols, lchnk )
          !!output the vertically integrated values;
          vsuffix  = trim(adjustl(hist_var2d(ihist)))
          varname  = trim(adjustl(vsuffix))//'_'//trim(adjustl(pname)) ! form variable name
          !find the prognostic variable associated with this hist_var2d(ihist) via "get_var_2d" function
          call outfld( trim(adjustl(varname)), get_var_2d(state,vsuffix), pcols, lchnk )
       enddo

    endif
  end subroutine physics_update
  
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------

  function get_var_3d(state,hist_var) result (prg_var)
    !Purpose: Find which state variable to output based on the hist_var string

    use constituents, only: cnst_get_ind
    use physconst,    only: cpair, latvap, gravit, rga
   
    character(len=fieldname_len), intent(in)  :: hist_var
    type(physics_state),          intent(in)  :: state  
    real(r8)                                  :: prg_var(pcols,pver)
    real(r8) esl(pcols,pver)   ! saturation vapor pressures 
    real(r8) esi(pcols,pver)   ! 
    real(r8) ftem(pcols,pver)  ! temporary workspace
    real(r8) gtem(pcols,pver)  ! temporary workspace

    !local vars
    integer :: idx,k

    !see if the hist_var exists in constituent array
    call cnst_get_ind(trim(adjustl(hist_var)), idx, abort=.false.)
    
    if (idx .ne. -1 ) then ! idx == -1  means, variable doesn't exists in the constituent array

       prg_var(1:pcols,1:pver) = state%q(1:pcols,1:pver,idx)*state%pdel(1:pcols,1:pver)*rga

    else !variable doesn't exists in the constituent array

       select case (trim(adjustl(hist_var)))
       case('s')
          prg_var(1:pcols,1:pver) = state%s(1:pcols,1:pver)
       case('t')
          prg_var(1:pcols,1:pver) = cpair*state%t(1:pcols,1:pver) * state%pdel(1:pcols,1:pver)*rga
                                    !state%t(1:pcols,1:pver)
       case('u')
          prg_var(1:pcols,1:pver) = state%u(1:pcols,1:pver) * state%pdel(1:pcols,1:pver)*rga

       case('v')
          prg_var(1:pcols,1:pver) = state%v(1:pcols,1:pver) * state%pdel(1:pcols,1:pver)*rga

       case('QSW')
         ! calculate from CAM q and t using CAM built-in functions
         call qsat_water(state%t(1:pcols,1:pver), state%pmid(1:pcols,1:pver), &
              esl(1:pcols,1:pver), ftem(1:pcols,1:pver))
              prg_var(1:pcols,1:pver) = ftem(1:pcols,1:pver) * state%pdel(1:pcols,1:pver)*rga
       case('RHW')
         ! calculate from CAM q and t using CAM built-in functions
         call qsat_water(state%t(1:pcols,1:pver), state%pmid(1:pcols,1:pver), &
              esl(1:pcols,1:pver), prg_var(1:pcols,1:pver))
              prg_var(1:pcols,1:pver) = state%q(1:pcols,1:pver,1)/prg_var(1:pcols,1:pver) * 100._r8
       case('QSI')
         ! calculate from CAM q and t using CAM built-in functions
         call qsat_ice(state%t(1:pcols,1:pver), state%pmid(1:pcols,1:pver), &
              esi(1:pcols,1:pver), ftem(1:pcols,1:pver))
              prg_var(1:pcols,1:pver) = ftem(1:pcols,1:pver) * state%pdel(1:pcols,1:pver)*rga
       case('RHI')
         ! calculate from CAM q and t using CAM built-in functions
         call qsat_water(state%t(1:pcols,1:pver), state%pmid(1:pcols,1:pver), &
              esl(1:pcols,1:pver), ftem(1:pcols,1:pver))
         gtem(1:pcols,1:pver) = state%q(1:pcols,1:pver,1)/ftem(1:pcols,1:pver) * 100._r8
         ! convert to RHI (ice)
         esi(1:pcols,1:pver)=svp_ice(state%t(1:pcols,1:pver))
         prg_var(1:pcols,1:pver)=gtem(1:pcols,1:pver)*esl(1:pcols,1:pver)/esi(1:pcols,1:pver)
       case('SUPERSAT_FLAG')
         ! calculate from CAM q and t using CAM built-in functions
         call qsat_water(state%t(1:pcols,1:pver), state%pmid(1:pcols,1:pver), &
              esl(1:pcols,1:pver), ftem(1:pcols,1:pver))
         gtem(1:pcols,1:pver) = state%q(1:pcols,1:pver,1)
         prg_var(1:pcols,1:pver) = 1._r8
         where(ftem(1:pcols,1:pver) .gt. gtem(1:pcols,1:pver))
          prg_var = 0._r8
         end where
       case default
         call endrun('physics_update_mod.F90 - func get_var_3d, unrecognized variable: '// trim(adjustl(hist_var)))
       end select
    endif
  end function get_var_3d

  function get_var_2d(state,hist_var) result (prg_var)
    !Purpose: Find which state variable to output based on the hist_var string

    use constituents, only: cnst_get_ind
    use physconst,    only: cpair, latvap, gravit, rga

    character(len=fieldname_len), intent(in)  :: hist_var
    type(physics_state),          intent(in)  :: state
    real(r8)                                  :: prg_var(pcols)
    real(r8) esl(pcols,pver)   ! saturation vapor pressures 
    real(r8) esi(pcols,pver)   ! 
    real(r8) ftem(pcols,pver)  ! temporary workspace
    real(r8) gtem(pcols,pver)  ! temporary workspace

    !local vars
    integer :: idx,k

    !see if the hist_var exists in constituent array
    select case (trim(adjustl(hist_var)))
    case('tvint') 
      ftem(1:pcols,1:pver) = cpair*state%t(1:pcols,1:pver) * state%pdel(1:pcols,1:pver)*rga
      do k=2,pver
         ftem(1:pcols,1) = ftem(1:pcols,1) + ftem(1:pcols,k)
      end do
      prg_var(1:pcols) = ftem(1:pcols,1)
    case('uvint')
      ftem(1:pcols,1:pver) = state%u(1:pcols,1:pver) * state%pdel(1:pcols,1:pver)*rga
      do k=2,pver
         ftem(1:pcols,1) = ftem(1:pcols,1) + ftem(1:pcols,k)
      end do
      prg_var(1:pcols) = ftem(1:pcols,1)
    case('vvint')
      ftem(1:pcols,1:pver) = state%v(1:pcols,1:pver) * state%pdel(1:pcols,1:pver)*rga
      do k=2,pver
         ftem(1:pcols,1) = ftem(1:pcols,1) + ftem(1:pcols,k)
      end do
      prg_var(1:pcols) = ftem(1:pcols,1)
    case('QVINT')
      call cnst_get_ind('Q', idx, abort=.false.)
      ftem(1:pcols,1:pver) = state%q(1:pcols,1:pver,idx) * state%pdel(1:pcols,1:pver)*rga
      do k=2,pver
         ftem(1:pcols,1) = ftem(1:pcols,1) + ftem(1:pcols,k)
      end do
      prg_var(1:pcols) = ftem(1:pcols,1)
    case('CLVINT')
      call cnst_get_ind('CLDLIQ', idx, abort=.false.)
      ftem(1:pcols,1:pver) = state%q(1:pcols,1:pver,idx) * state%pdel(1:pcols,1:pver)*rga
      do k=2,pver
         ftem(1:pcols,1) = ftem(1:pcols,1) + ftem(1:pcols,k)
      end do
      prg_var(1:pcols) = ftem(1:pcols,1)
     case('NLVINT')
      call cnst_get_ind('NUMLIQ', idx, abort=.false.)
      ftem(1:pcols,1:pver) = state%q(1:pcols,1:pver,idx) * state%pdel(1:pcols,1:pver)*rga
      do k=2,pver
         ftem(1:pcols,1) = ftem(1:pcols,1) + ftem(1:pcols,k)
      end do
      prg_var(1:pcols) = ftem(1:pcols,1)
    case('CIVINT')
      call cnst_get_ind('CLDICE', idx, abort=.false.)
      ftem(1:pcols,1:pver) = state%q(1:pcols,1:pver,idx) * state%pdel(1:pcols,1:pver)*rga
      do k=2,pver
         ftem(1:pcols,1) = ftem(1:pcols,1) + ftem(1:pcols,k)
      end do
      prg_var(1:pcols) = ftem(1:pcols,1)
    case('NIVINT')
      call cnst_get_ind('NUMICE', idx, abort=.false.)
      ftem(1:pcols,1:pver) = state%q(1:pcols,1:pver,idx) * state%pdel(1:pcols,1:pver)*rga
      do k=2,pver
         ftem(1:pcols,1) = ftem(1:pcols,1) + ftem(1:pcols,k)
      end do
      prg_var(1:pcols) = ftem(1:pcols,1)
    case('RIVINT')
      call cnst_get_ind('RAINQM', idx, abort=.false.)
      ftem(1:pcols,1:pver) = state%q(1:pcols,1:pver,idx) * state%pdel(1:pcols,1:pver)*rga
      do k=2,pver
         ftem(1:pcols,1) = ftem(1:pcols,1) + ftem(1:pcols,k)
      end do
      prg_var(1:pcols) = ftem(1:pcols,1)
    case('NRVINT')
      call cnst_get_ind('NUMRAI', idx, abort=.false.)
      ftem(1:pcols,1:pver) = state%q(1:pcols,1:pver,idx) * state%pdel(1:pcols,1:pver)*rga
      do k=2,pver
         ftem(1:pcols,1) = ftem(1:pcols,1) + ftem(1:pcols,k)
      end do
      prg_var(1:pcols) = ftem(1:pcols,1)
    case('SOVINT')
      call cnst_get_ind('SNOWQM', idx, abort=.false.)
      ftem(1:pcols,1:pver) = state%q(1:pcols,1:pver,idx) * state%pdel(1:pcols,1:pver)*rga
      do k=2,pver
         ftem(1:pcols,1) = ftem(1:pcols,1) + ftem(1:pcols,k)
      end do
      prg_var(1:pcols) = ftem(1:pcols,1)
    case('NSVINT')
      call cnst_get_ind('NUMSNO', idx, abort=.false.)
      ftem(1:pcols,1:pver) = state%q(1:pcols,1:pver,idx) * state%pdel(1:pcols,1:pver)*rga
      do k=2,pver
         ftem(1:pcols,1) = ftem(1:pcols,1) + ftem(1:pcols,k)
      end do
      prg_var(1:pcols) = ftem(1:pcols,1)
    case('QSWVIN')
      ! calculate from CAM q and t using CAM built-in functions
      call qsat_water(state%t(1:pcols,1:pver), state%pmid(1:pcols,1:pver), &
           esl(1:pcols,1:pver), gtem(1:pcols,1:pver))
      ftem(1:pcols,1:pver) = gtem(1:pcols,1:pver) * state%pdel(1:pcols,1:pver)*rga
      do k=2,pver
         ftem(1:pcols,1) = ftem(1:pcols,1) + ftem(1:pcols,k)
      end do
      prg_var(1:pcols) = ftem(1:pcols,1)
    case('RHWVIN')
      ! calculate from CAM q and t using CAM built-in functions
      call qsat_water(state%t(1:pcols,1:pver), state%pmid(1:pcols,1:pver), &
           esl(1:pcols,1:pver), ftem(1:pcols,1:pver))
      gtem(1:pcols,1:pver) = state%q(1:pcols,1:pver,1)/ftem(1:pcols,1:pver) * 100._r8
      ftem(1:pcols,1:pver) = gtem(1:pcols,1:pver) - 100._r8
      where(gtem(1:pcols,1:pver) .lt. 100._r8)
       ftem = 0._r8
      end where
      do k=2,pver
         ftem(1:pcols,1) = ftem(1:pcols,1) + ftem(1:pcols,k)
      end do
      prg_var(1:pcols) = ftem(1:pcols,1)
    case('QSIVIN')
      ! calculate from CAM q and t using CAM built-in functions
      call qsat_ice(state%t(1:pcols,1:pver), state%pmid(1:pcols,1:pver), &
           esi(1:pcols,1:pver), gtem(1:pcols,1:pver))
      ftem(1:pcols,1:pver) = gtem(1:pcols,1:pver) * state%pdel(1:pcols,1:pver)*rga
      do k=2,pver
         ftem(1:pcols,1) = ftem(1:pcols,1) + ftem(1:pcols,k)
      end do
      prg_var(1:pcols) = ftem(1:pcols,1)
    case('RHIVIN')
      ! calculate from CAM q and t using CAM built-in functions
      call qsat_water(state%t(1:pcols,1:pver), state%pmid(1:pcols,1:pver), &
           esl(1:pcols,1:pver), ftem(1:pcols,1:pver))
      ftem(1:pcols,1:pver) = state%q(1:pcols,1:pver,1)/ftem(1:pcols,1:pver) * 100._r8
      ! convert to RHI (ice)
      esi(1:pcols,1:pver)  = svp_ice(state%t(1:pcols,1:pver))
      gtem(1:pcols,1:pver) = ftem(1:pcols,1:pver)*esl(1:pcols,1:pver)/esi(1:pcols,1:pver) 
      ftem(1:pcols,1:pver) = gtem(1:pcols,1:pver) - 100._r8
      where(gtem(1:pcols,1:pver) .lt. 100._r8)
       ftem = 0._r8
      end where
      do k=2,pver
         ftem(1:pcols,1) = ftem(1:pcols,1) + ftem(1:pcols,k)
      end do
      prg_var(1:pcols) = ftem(1:pcols,1)
    case default
      call endrun('physics_update_mod.F90 - func get_var_2d, unrecognized variable: '// trim(adjustl(hist_var)))
    end select

  end function get_var_2d

end module physics_update_mod
