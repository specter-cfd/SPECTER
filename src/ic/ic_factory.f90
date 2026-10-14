! ===================================================================
! Factory for all initial conditions (ICs)
!
! This factory applies chain operation to create multiple initial
! conditions for all fields in a pde solver.
!
! Versión de SPECTER del ic_factory.f90 de GHOST (rama new_paradigm,
! afaf9ea), con las condiciones iniciales de la velocidad que ya están
! en ic_velocity.f90.
!
! DATE : 04/08/26 (PDM)
! ===================================================================

module ic_factory
  USE icbase_mod
  USE ic_velocity

  IMPLICIT NONE

CONTAINS

  ! ================= Factory function ==============================
  function init_ic_from_file(infile) result(new_object)
    USE gstate_mod
    USE commtypes
    USE mpivars, only: myrank, ierr
    USE gutils
    type   (icChain), allocatable :: new_object(:)
    character(len=*),  intent(in) :: infile

    ! Temporary data to read from namelist:
    integer            :: i,icnmb
    character(len=256) :: iclist
    character(len= 25) :: icname(10)
    ! Required namelist:
    namelist/ initial_conditions / iclist

    if ( myrank .eq. 0 ) then
      open(1,file=infile,status='unknown',form="formatted")
      read(1,NML=initial_conditions)
      close(1)
    endif
    call MPI_BCAST(iclist,256,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

    ! Loops through all requested ICs and allocates classes
    call parsestr(iclist,';',icname,25,10,icnmb)
    if ( myrank .eq. 0 ) print *,'Found',icnmb,'ICs'
    allocate( new_object(icnmb) )
    do i = 1,icnmb
      if ( myrank .eq. 0 ) print *,'IC',i,':',trim(adjustl(icname(i)))
      select case (trim(adjustl(icname(i))))
      ! Velocity field ICs ---------------------
      case ('read_v')
        allocate( icRead_v      :: new_object(i)%ic )
      case ('null_v')
        allocate( icNull_v      :: new_object(i)%ic )
      case ('cbnoise_v')
        allocate( icCbnoise_v   :: new_object(i)%ic )
      case ('couette_v')
        allocate( icCouette_v   :: new_object(i)%ic )
      case ('poiseuille_v')
        allocate( icPoiseuille_v :: new_object(i)%ic )
      case default
        stop 'Unknown initial conditions'
      end select
    end do
  end function init_ic_from_file

end module ic_factory
