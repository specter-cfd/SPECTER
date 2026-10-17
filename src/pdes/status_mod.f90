!=================================================================
! General MODULES for pde methods
!
! Adaptado de status_mod.f90 de GHOST (rama new_paradigm, cb6943d).
! Diferencias: suma tdir (tablas del FC-Gram), no tiene sstep ni los
! contadores de los espectros (SPECTER no los escribe) y status_init
! no usa grid, kes ni ali.
!
! 2026 Pablo D. Mininni.
!      Department of Physics,
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar
!=================================================================

  MODULE status
!
! General status flags for all pde/time integrators
! Reads a "&status" namelist with:
!  idir : default directory for unformatted input
!  odir : default directory for unformatted output
!  tdir : directory for FC-Gram tables
!  stat : = 0 starts a new run
!         OR  gives the number of the file used to continue a run
!  mult : time step multiplier
!  bench: = 0 production run
!         = 1 benchmark run (no I/O)
!         = 2 higher level benchmark run (+time to create plans)
!  outs : = 0 writes velocity
!         = 1 writes vorticity
!  dt   : time step size
!  step : total number of time steps to compute
!  tstep: number of steps between binary output
!  cstep: number of steps between output of global quantities
!  cort : global correlation time for random forcings (when needed)
!  seed : global seed for random number generation (when needed)
      USE fprecision
      REAL(KIND=GP)      :: dt
      INTEGER            :: stat ,ini
      INTEGER            :: step
      INTEGER            :: tstep,cstep
      INTEGER            :: fstep
      INTEGER            :: bench = 0
      INTEGER            :: outs
      INTEGER            :: mult
      INTEGER            :: seed
      INTEGER            :: tind
      INTEGER            :: timet,timec
      INTEGER            :: timef
      CHARACTER(len=128) :: odir,idir,tdir
    CONTAINS
!-----------------------------------------------------------------
      SUBROUTINE status_init(infile_)
!
! Initializes general global status parameters for the pde solvers
      USE commtypes
      USE mpivars
      IMPLICIT NONE

      REAL(KIND=GP)                :: cort
      CHARACTER(len=*), INTENT(IN) :: infile_
      NAMELIST / status / idir,odir,tdir,stat,mult,bench,outs,dt
      NAMELIST / status / step,tstep,cstep,seed,cort

      cort = 1000.0_GP
      IF (myrank.eq.0) THEN
         OPEN(1,file=infile_,status='unknown',form="formatted")
         READ(1,NML=status)
         CLOSE(1)
         dt = dt/real(mult,kind=GP)
         step = step*mult
         tstep = tstep*mult
         cstep = cstep*mult
         fstep = int(cort/dt)
      ENDIF
      CALL MPI_BCAST(idir ,128,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(odir ,128,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(tdir ,128,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(stat ,1  ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mult ,1  ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(bench,1  ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(outs ,1  ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(dt   ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(step ,1  ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(tstep,1  ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cstep,1  ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(fstep,1  ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(seed ,1  ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      END SUBROUTINE status_init

  END MODULE status

