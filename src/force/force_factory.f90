! ===================================================================
! Factory for all forcing functions
!
! This factory applies chain operation to create multiple forcings
! and forcing update methods for all fields in a pde solver.
!
! Versión de SPECTER del force_factory.f90 de GHOST (rama
! new_paradigm, cb6943d), con los forzados de la velocidad de
! force_velocity.f90 y solo el método constant_fv (forzado constante
! en el tiempo). Sin métodos shuffle no hace falta agrandar el pool,
! así que workspace no se usa.
!
! DATE : 01/18/26 (PDM)
! ===================================================================

module force_factory
  USE forcebase_mod
  USE force_velocity

  IMPLICIT NONE

CONTAINS

  ! ================= Factory function ==============================
  function init_forcing_from_file(infile,workspace) result(new_object)
    USE class_GWorkspace3D
    USE gstate_mod
    USE commtypes
    USE gutils
    type(forceChain),   allocatable :: new_object(:)
    type(GWorkspace), intent(inout) :: workspace
    character(len=*), intent   (in) :: infile
    ! Temporary data:
    real(kind=GP)      :: cort
    integer            :: i,fnmb,unmb
    character(len=512) :: forces
    character(len= 50) :: forcelist(10)
    character( len=25) :: forcename(2)
    ! Required namelist:
    namelist/ forcing / forces,cort

    if ( myrank .eq. 0 ) then
      open(1,file=infile,status='unknown',form="formatted")
      read(1,NML=forcing)
      close(1)
    endif
    call MPI_BCAST(forces,512,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(cort  ,  1,GC_REAL      ,0,MPI_COMM_WORLD,ierr)

    ! Loops through all requested forces and allocates classes.
    call parsestr(forces,';',forcelist ,50,10,fnmb)
    if ( myrank .eq. 0 ) then
       print *,'Found',fnmb,'force methods'
    endif
    allocate( new_object(fnmb) )
    do i = 1,fnmb
      call parsestr(forcelist(i),',',forcename,25,2,unmb)
      if ( myrank .eq. 0 ) print *,'Forcing',i,':', &
           trim(adjustl(forcename(1))),',',trim(adjustl(forcename(2)))
      ! Forcing functions  ------------------------------------------
      select case (trim(adjustl(forcename(1))))
      ! Velocity forcing functions  ------------
      case ('null_fv')
        allocate( forceNull_fv    :: new_object(i)%force )
      case ('uniform_fv')
        allocate( forceUniform_fv :: new_object(i)%force )
      case default
        stop 'Unknown forcing function'
      end select
      ! Update methods ----------------------------------------------
      select case (trim(adjustl(forcename(2))))
      ! Velocity forcing update methods --------
      case ('constant_fv')
        if ( allocated(new_object(i)%update) ) deallocate(new_object(i)%update)
      case default
        stop 'Unknown or undefined forcing update method'
      end select
    end do
  end function init_forcing_from_file

end module force_factory
