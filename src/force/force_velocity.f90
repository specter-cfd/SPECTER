! =====================================================================
! NAME       : force_velocity.f90
! DESCRIPTION: Forcing methods for all VelocityBase solver classes.
!              Each force needs an initial set up, and an update.
!              These forcing methods can be used to initialize and
!              update forcings on other solver classes, using
!              composition of operators, as long as no correlation
!              between the forces is needed.
!
! Forces avaliable:
!   null_fv    : Null forcing of the velocity field
!   uniform_fv : Uniform forcing along x (mean pressure gradient)
!
! Update methods available:
!   constant_fv: Constant forcing (no method)
!
! Versión de SPECTER del force_velocity.f90 de GHOST (rama
! new_paradigm, cb6943d). null_fv es el de GHOST y uniform_fv es el
! src/initialfv.f90 de SPECTER. tg_fv, abc_fv y random_fv, y los
! métodos shift_fv y shuffle_fv, no están: suponen un dominio
! periódico en las tres direcciones.
!
! DATE       : 01/18/26 (PDM)
! =====================================================================

module force_velocity
  USE forcebase_mod

  IMPLICIT NONE

  ! ================= Forcing functions supported =====================
  type, extends(forceBase) :: forceNull_fv
    contains
      procedure :: init_GForce => init_nullfv
  end type forceNull_fv
  type, extends(forceBase) :: forceUniform_fv
    contains
      procedure :: init_GForce => init_uniformfv
  end type forceUniform_fv

CONTAINS

  ! ===================================================================
  ! Forcing functions
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Null forcing
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_nullfv(this,solver,state)
    use gstate_mod
    use hd_mod
    use grid
    use mpivars
    implicit none

    class(forceNull_fv), intent   (in) :: this
    class(EquationBase), intent(inout) :: solver
    type   (GStateComp), intent(inout) :: state(:)
    integer                            :: i,j,k

    select type (solver)
    class is (VelocityBase)
!$omp parallel do collapse(2) private (k)
      DO i = ista,iend
         DO j = 1,ny
            DO CONCURRENT (k=1:nz)
              state(solver%VELOCITY  )%ccomp(k,j,i) = 0.0_GP
              state(solver%VELOCITY+1)%ccomp(k,j,i) = 0.0_GP
              state(solver%VELOCITY+2)%ccomp(k,j,i) = 0.0_GP
            END DO
         END DO
      END DO
    class default
      error stop "This solver does not support velocity forcing"
    end select
  end subroutine init_nullfv


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Forzado uniforme en x, f = (f0,0,0): el gradiente de
  !! presión medio que mueve el flujo de Poiseuille. En Fourier
  !! es solo el modo (kz,ky,kx) = (0,0,0), multiplicado por
  !! nx.ny.nz porque la antitransformada no normaliza.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   f0: amplitud del forzado
  subroutine init_uniformfv(this,solver,state)
    use gstate_mod
    use hd_mod
    use grid
    use mpivars
    use commtypes
    implicit none

    class(forceUniform_fv), intent   (in) :: this
    class(EquationBase)   , intent(inout) :: solver
    type   (GStateComp)   , intent(inout) :: state(:)
    real   (kind=GP)                      :: f0
    integer                               :: i,j,k

    namelist/ uniform_fv / f0
    select type (solver)
    class is (VelocityBase)
    ! Read parameters from a namelist in the input file
    if ( myrank .eq. 0 ) then
      open(1,file=solver%infile_,status='unknown',form="formatted")
      read(1,NML=uniform_fv)
      close(1)
    endif
    call mpi_bcast(f0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
!$omp parallel do collapse(2) private (k)
      DO i = ista,iend
         DO j = 1,ny
            DO CONCURRENT (k=1:nz)
              state(solver%VELOCITY  )%ccomp(k,j,i) = 0.0_GP
              state(solver%VELOCITY+1)%ccomp(k,j,i) = 0.0_GP
              state(solver%VELOCITY+2)%ccomp(k,j,i) = 0.0_GP
            END DO
         END DO
      END DO
    ! El modo (0,0,0) lo tiene la tarea con kx = 0 (ista = 1)
    if ( ista .eq. 1 ) then
      state(solver%VELOCITY)%ccomp(1,1,1) = &
           f0*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
    endif
    class default
      error stop "This solver does not support velocity forcing"
    end select
  end subroutine init_uniformfv

end module force_velocity
