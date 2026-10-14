! ===================================================================
! Factory for all PDEs.
!
! To add a new solver, "USE" the corresponding module, and add a new
! case clause in the init_pdes_from_file function. The clause should
! allocate the solver class, and declare the number of field
! components the solver needs.
!
! Versión de SPECTER del equation_factory.f90 de GHOST (rama
! new_paradigm, afaf9ea), solo con el solver HD. El pool de
! temporarios es el pico de SPECTER (6 complejos y 3 reales, los de
! write_states, medidos con report_peaks), y el pool solo-host no se
! usa.
!
! DATE : 03/29/26 (PDM)
! ===================================================================

module equation_factory
  USE equationbase_mod
  USE hd_mod
! USE userdefinedpde_mod

  IMPLICIT NONE

  ! ================= Global parameters ===============================
  ! Number of complex and real temporary arrays in the workspace pool,
  ! the peak of the nested requests of each solver as measured with
  ! report_peaks on the example inputs.
  integer, public :: NUMTMPCOMP = 0 ! Number of cmplx tmp arrays
  integer, public :: NUMTMPREAL = 0 ! Number of real tmp arrays
  ! Host-only temporaries, for the diagnostics (not used in SPECTER).
  integer, public :: NUMTMPHCOMP = 0 ! Number of host-only cmplx tmp arrays
  integer, public :: NUMTMPHREAL = 0 ! Number of host-only real tmp arrays

CONTAINS

  ! ================= Factory function ==============================
  function init_pdes_from_file(infile) result(new_object)
    USE commtypes
    class(EquationBase), allocatable  :: new_object
    character(len=*)   , intent(in)   :: infile

    ! Temporary data to read from namelist:
    integer           :: nprocs,myrank,ierr
    character(len=64) :: solver
    ! Required namelist:
    namelist/ pdes / solver

    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
    if ( myrank .eq. 0 ) then
      open(1,file=infile,status='unknown',form="formatted")
      read(1,NML=pdes)
      close(1)
    endif
    call MPI_BCAST(solver,64,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

    ! Clauses for each solver class
    select case (trim(adjustl(solver)))
      case ('HD')
        allocate(HDsolver    :: new_object)
        NUMTMPCOMP =  6; NUMTMPREAL = 3
      case default
        stop 'Equation factory :: init_pdes_from_file : Unknown solver name'
    end select
  end function init_pdes_from_file

end module equation_factory
