!=================================================================
      PROGRAM SPECTER
!=================================================================
! SPECTER code: Special PEriodic Continuation Turbulence solvER
!
! SPECTER performs numerical integration for a set of PDEs usually
! found in fluid dynamics. A pseudo-spectral method is used to
! compute spatial derivatives, coupled with an adjustable order
! Runge-Kutta time stepping. When non-periodic directions are
! specified, it utilizes FC-Gram continuations to compute Fourier
! transforms with minimal overhead and considerable accuracy.
! To compile, you need the FFTW and MPI libraries installed on
! your system. More info in the attached README file or in 
! https://github.com/mfontanaar/SPECTER.
!
! Notation: index 'i' moves in the 'x' direction (or kx)
!           index 'j' moves in the 'y' direction (or ky)
!           index 'k' moves in the 'z' direction (or kz)
!
! Organizado como el main de GHOST (rama new_paradigm, afaf9ea): el
! solver (&pdes), las condiciones iniciales (&initial_conditions), el
! forzado (&forcing) y el integrador (&stepper) se eligen al correr,
! con sus fábricas. Por ahora el único solver es HD.
!
! 2019 Mauro Fontana
!      Department of Physics,
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mfontana@df.uba.ar
!
! Feb 2020: Initial release.
! Jul 2021: Added MHD support.
!
! References:
! - Fontana M, Bruno OP, Mininni PD, Dmitruk P; Comp. Phys. Comm. 256,
!     107482 (2020) -HD,BOUSS-
! - Fontana M, Mininni PD, Bruno OP, Dmitruk P; arXiv: 2107.07077 (2021) -MHD-
! - Rosenberg DL, Mininni PD, Reddy R, Pouquet A.; Atmosphere 11, 178 (2020)
! - Mininni PD, Rosenberg DL, Reddy R, Pouquet A.; P. Comp. 37, 123 (2011)
! - Bruno OP, Lyon M; Jour. Comp. Phys. 229, 2009 (2010)
!
! Modules
      USE commtypes
      USE filefmt
      USE iovar
      USE fft
      USE fcgram
      USE kes
      USE threads
      USE boxsize
      USE status
      USE gtimer
      USE ic_factory
      USE force_factory
      USE equation_factory
      USE stepper_factory

      IMPLICIT NONE

      TYPE(GStateComp), ALLOCATABLE, TARGET :: field(:), field_nxt(:), force(:)

      TYPE(GWorkspace), TARGET :: workspace

!
! Book keeping variables, options holders and clocks 
      INTEGER :: t
      INTEGER :: num_components
      INTEGER :: ihcpu1,ihcpu2
      INTEGER :: ihomp1,ihomp2
      INTEGER :: ihwtm1,ihwtm2
!$    INTEGER, EXTERNAL     :: omp_get_max_threads

      CLASS(EquationBase), ALLOCATABLE, TARGET :: fluid
      CLASS(GStepperBase), ALLOCATABLE :: stepper
      CLASS(icChain), ALLOCATABLE :: iclist(:)
      CLASS(forceChain), ALLOCATABLE :: forcemethod(:)

      TYPE(IOPLAN), TARGET :: planio

      LOGICAL               :: bbenchexist

! Initializes the MPI and I/O libraries
      CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED,provided,ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

! Initializes the grid and the FC-Gram continuation (&grid, see
! pseudo/pseudospec_mod.f90). This must be done early to have nx, ny, nz.
      CALL grid_init('parameter.inp')

! Create partitions across tasks
     CALL range(1,nx/2+1,nprocs,myrank,ista,iend)
     CALL range(1,nz,nprocs,myrank,ksta,kend)
! Initialize array IO
      CALL io_init(myrank,(/nx-Cx,ny-Cy,nz-Cz/),ksta,kend,planio)

! PDE solver (&pdes, see pdes/equation_factory.f90) and the pool of
! temporary arrays, sized for that solver
      fluid = init_pdes_from_file('parameter.inp')
      CALL workspace%initialize_pool(NUMTMPREAL,NUMTMPCOMP)
      CALL workspace%init_host_entries(NUMTMPHREAL,NUMTMPHCOMP)

! Reads the namelist 'status' on the external file
! 'parameter.inp' (see pdes/status_mod.f90)
      CALL status_init('parameter.inp')

!
! Box size and grids (see pseudo/pseudospec_mod.f90)
      CALL box_init('parameter.inp')

 ! Initializes the FFT library. This must be done at
 ! this stage as it requires the variable "bench" to
 ! be properly initialized.
 ! Use FFTW_ESTIMATE or FFTW_MEASURE in short runs
 ! Use FFTW_PATIENT or FFTW_EXHAUSTIVE in long runs
 ! FFTW 2.x only supports FFTW_ESTIMATE or FFTW_MEASURE

      nth = 1
!$    nth = omp_get_max_threads()
!$    CALL fftp3d_init_threads(ierr)
      IF (bench.eq.2) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL GTStart(ihcpu2,GT_CPUTIME)
         CALL GTStart(ihomp2,GT_OMPTIME)
         CALL GTStart(ihwtm2,GT_WTIME)
      ENDIF

      CALL fcgram_create_plan(planfc, (/nx,ny,nz/), (/Cx,Cy,Cz/), &
                              (/ox,oy,oz/), (/x(2),y(2),z(2)/), &
                              tdir, FFTW_ESTIMATE)
      IF (bench.eq.2) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL GTStop(ihcpu2)
         CALL GTStop(ihomp2)
         CALL GTStop(ihwtm2)
      ENDIF

! Solver construction (reads its namelists and sets up the boundary
! conditions and the pressure; it needs the FC-Gram plan)
      CALL fluid%Solver_ctor('parameter.inp',workspace,planio)

      num_components = fluid%state_size()
      CALL GState_alloc(field    , num_components)
      CALL GState_alloc(field_nxt, num_components)
      CALL GState_alloc(force    , num_components)

! Initial conditions (see ic/ic_factory.f90)
      iclist = init_ic_from_file('parameter.inp')

! Forcing (see force/force_factory.f90)
      forcemethod = init_forcing_from_file('parameter.inp',workspace)

! If stat=0 we start a new run, and if stat.ne.0 a previous run is
! continued: only the counters depend on stat (the initial conditions
! are set below by the list iclist)
 IC : IF (stat.eq.0) THEN
      ini = 1
      tind = 0                          ! index for the binaries
      timet = tstep
      timec = cstep
      timef = 0
      ELSE
! If stat.ne.0 a previous run is continued
      ini = int((stat-1)*tstep) + 1
      tind = int(stat)
      timet = 0
      timec = int(modulo(float(ini-1),float(cstep)))
      timef = int(modulo(float(ini-1),float(fstep)))
      ENDIF IC

! Initial conditions for the fields (iclist in &initial_conditions);
! field_nxt is the state that is advanced and written
      CALL init_allstates(iclist,fluid,field)
      CALL GState_copy(field_nxt,field)

      CALL init_forcing(forcemethod,fluid,force)

! Sets up the time stepper (&stepper, see steppers/stepper_factory.f90)
      CALL build_stepper_from_file('parameter.inp',stepper,workspace,fluid)

!
! Time integration scheme starts here. If
! we are doing a benchmark, we measure
! cputime before starting.

      IF (bench.eq.1) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL GTStart(ihcpu1,GT_CPUTIME)
         CALL GTStart(ihomp1,GT_OMPTIME)
         CALL GTStart(ihwtm1,GT_WTIME)
! Must reset FFT internal timers after initialization:
         comtime = 0.0D0
         ffttime = 0.0D0
         tratime = 0.0D0
      ENDIF


 RK : DO t = ini,step
! Every 'tstep' steps, stores the fields 
! in binary files
 BIN:    IF ((timet.eq.tstep).and.(bench.eq.0)) THEN
            timet = 0
            tind = tind+1
            CALL fluid%write_states(field_nxt,planio)
         ENDIF BIN

! Every 'cstep' steps, generates external files 
! with global quantities.

         IF ((timec.eq.cstep).and.(bench.eq.0)) THEN
            timec = 0
            CALL fluid%global(field_nxt,force,t)
         ENDIF

! Forcing update (nothing to do for a constant forcing)
         CALL update_forcing(forcemethod,fluid,force)

! Time step (GHOST Canuto stepper): field keeps the state at the
! beginning of the step and field_nxt is advanced
         CALL GState_copy(field,field_nxt)
         CALL stepper%gstep((t-1)*dt,field,force,dt,field_nxt)

         timet = timet+1
         timec = timec+1
      END DO RK
!
! End of Runge-Kutta

! Computes the benchmark
      IF (bench.gt.0) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL GTStop(ihcpu1)
         CALL GTStop(ihomp1)
         CALL GTStop(ihwtm1)
         IF (myrank.eq.0) THEN
            INQUIRE( file='benchmark.txt', exist=bbenchexist )
            OPEN(1,file='benchmark.txt',position='append')
            IF ( .NOT. bbenchexist ) THEN
               WRITE(1,*) &
	       '# nx ny nz nsteps nprocs nth TCPU TOMP TWTIME TFFT TTRA TCOM',&
           'TCONT TNEU TROB TTOT'
            ENDIF
            WRITE(1,*) nx,ny,nz,(step-ini+1),nprocs,nth, &
                       GTGetTime(ihcpu1)/(step-ini+1),   &
                       GTGetTime(ihomp1)/(step-ini+1),   &
                       GTGetTime(ihwtm1)/(step-ini+1),   &
                       ffttime/(step-ini+1), tratime/(step-ini+1), &
                       comtime/(step-ini+1), conttime/(step-ini+1),&
                       neutime/(step-ini+1), robtime/(step-ini+1), &
                       tottime/(step-ini+1)
            IF (bench.eq.2) THEN
               WRITE(1,*) 'FFTW: Create_plan = ',      &
                       GTGetTime(ihcpu2)/(step-ini+1), &
                       GTGetTime(ihomp2)/(step-ini+1), &
                       GTGetTime(ihwtm2)/(step-ini+1)
            ENDIF
            CLOSE(1)
         ENDIF
      ENDIF

! Peak usage of the pool of temporary arrays
      CALL workspace%report_peaks()

      CALL MPI_FINALIZE(ierr)

      CALL GTFree(ihcpu1)
      CALL GTFree(ihomp1)
      CALL GTFree(ihwtm1)
      CALL GTFree(ihcpu2)
      CALL GTFree(ihomp2)
      CALL GTFree(ihwtm2)

      CALL fcgram_destroy_plan(planfc)

      CALL GState_dealloc(field)
      CALL GState_dealloc(field_nxt)
      CALL GState_dealloc(force)

      DEALLOCATE( x,y,z )
      DEALLOCATE( kx,ky,kz )
      DEALLOCATE( kk2,khom )

      END PROGRAM SPECTER
