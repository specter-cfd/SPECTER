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
! Solvers (chosen at run time with &pdes):
!           HD         Hydrodynamic solver
!
! 2019 Mauro Fontana
!      Department of Physics,
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mfontana@df.uba.ar
!
! Feb 2020: Initial release.
! Jul 2021: Added MHD support.
! Oct 2026: Ported to the GHOST new paradigm.
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
      USE pstatus
      USE gtimer
      USE ic_factory
      USE force_factory
      USE equation_factory
      USE particle_factory
      USE icp_factory
      USE stepper_factory
      USE gdevice
      USE gmem

      IMPLICIT NONE

      TYPE(GStateComp), ALLOCATABLE :: field(:), field_nxt(:), force(:)
      TYPE(GPStateComp), ALLOCATABLE :: part(:), part_nxt(:)

      TYPE(GWorkspace), TARGET :: workspace

!
! Book keeping variables, options holders and clocks
      INTEGER :: t
      INTEGER :: num_components
      INTEGER :: ihcpu1,ihcpu2
      INTEGER :: ihomp1,ihomp2
      INTEGER :: ihwtm1,ihwtm2

      CLASS(EquationBase), ALLOCATABLE, TARGET :: fluid
      CLASS(ParticleBase), ALLOCATABLE :: particle
      CLASS(GStepperBase), ALLOCATABLE :: stepper
      CLASS(icChain), ALLOCATABLE :: iclist(:)
      CLASS(forceChain), ALLOCATABLE :: forcemethod(:)
      CLASS(icpChain), ALLOCATABLE :: icplist(:)

      TYPE(IOPLAN), TARGET :: planio

      LOGICAL               :: bbenchexist
      LOGICAL               :: fstatic
      REAL                  :: rbal = 0.0

! Initializes the MPI and I/O libraries
      CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED,provided,ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

! Initializes the grid and the FC-Gram continuation (&grid, see
! pseudo/pseudospec_mod.f90). This must be done early to have nx, ny, nz.
      CALL grid_init('parameter.inp')

! Binds this MPI task to a GPU in offload builds (no-op otherwise).
! SPECTER assumes the number of MPI tasks per node is a multiple of the
! number of GPUs per node; the user must ensure this condition.
      CALL device_init(myrank)

! Create partitions across tasks
     CALL range(1,nx/2+1,nprocs,myrank,ista,iend)
     CALL range(1,nz,nprocs,myrank,ksta,kend)
! Initialize array IO
      CALL io_init(myrank,(/nx-Cx,ny-Cy,nz-Cz/),ksta,kend,planio)

! Reads the namelist 'status' on the external file 'parameter.inp'
! (see pdes/status_mod.f90), and the particle status '&pstatus'
! (see particles/pstatus_mod.f90; if absent, dopart is false and
! the run has no particles)
      CALL status_init ('parameter.inp')
      CALL pstatus_init('parameter.inp')

! PDE solver (&pdes, see pdes/equation_factory.f90), particle solver
! (&particles, see particles/particle_factory.f90), and the pool of
! temporary arrays, sized for them
      fluid = init_pdes_from_file('parameter.inp')
      if (dopart) particle = init_particles_from_file('parameter.inp')
      CALL workspace%initialize_pool(NUMTMPREAL,NUMTMPCOMP,NUMTMPPART)
      CALL workspace%init_host_entries(NUMTMPHREAL,NUMTMPHCOMP)

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
! Device copy of the FC-Gram table in z (no-op in CPU builds)
      IF (ALLOCATED(planfc%z%dir)) THEN
         CALL gdev_alloc(C_LOC(planfc%z%dir),INT(SIZE(planfc%z%dir),C_SIZE_T)* &
                         INT(STORAGE_SIZE(planfc%z%dir)/8,C_SIZE_T))
         CALL gdev_update_to(C_LOC(planfc%z%dir),INT(SIZE(planfc%z%dir),C_SIZE_T)* &
                             INT(STORAGE_SIZE(planfc%z%dir)/8,C_SIZE_T))
      ENDIF
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

! Initialization of the particles (needs the solver, for nu and the
! I/O directories, and the box; see particles/inerpart_solver.f90 and
! icpart/icp_factory.f90)
      IF (dopart) THEN
         CALL particle%part_ctor('parameter.inp',fluid,workspace,part,part_nxt)
         icplist = init_icp_from_file('parameter.inp')
      ENDIF

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
      pind = 0                          ! index for the particles
      timet = tstep
      timep = pstep
      timec = cstep
      timef = 0
      ELSE
! If stat.ne.0 a previous run is continued
      ini = int((stat-1)*tstep) + 1
      tind = int(stat)
      pind = int((stat-1)*lgmult+1)
      timet = 0
      timep = 0
      timec = int(modulo(float(ini-1),float(cstep)))
      timef = int(modulo(float(ini-1),float(fstep)))
      ENDIF IC

! Initial conditions for the fields (iclist in &initial_conditions);
! field_nxt is the state that is advanced and written
      CALL init_allstates(iclist,fluid,field)
      CALL GState_update_to(field)        ! Device copies (no-op in host builds).
      CALL GState_copy(field_nxt,field)

      CALL init_forcing(forcemethod,fluid,force)
      CALL GState_update_to(force)
      CALL fluid%sync_device()            ! Auxiliary arrays of the solver (if any).
      fstatic = forcing_is_static(forcemethod)

! Initial conditions for the particles (icplist in
! &part_initial_conditions); part_nxt is the state that is advanced
      IF (dopart) THEN
         CALL init_allpstates(icplist,fluid,field,particle,part)
         IF (size(part(1)%rcomp) .ne. size(part_nxt(1)%rcomp)) THEN
            CALL GPState_resize(part_nxt,particle%partbuff_) ! We resize part_nxt
         ENDIF
         CALL GPState_update_to(part)     ! Device copies of the particle
         CALL particle%sync_device()      ! states and of the class arrays.
         CALL GPState_copy(part_nxt,part) ! We also update part_nxt.
      ENDIF

! Sets up the time stepper (&stepper, see steppers/stepper_factory.f90)
      IF (dopart) THEN
         CALL build_stepper_from_file('parameter.inp',stepper,workspace,fluid,particle)
      ELSE
         CALL build_stepper_from_file('parameter.inp',stepper,workspace,fluid)
      ENDIF

!
! Time integration scheme starts here. In offload builds the fields are
! worked on the device inside the time step (gdev_active set) and on
! their host copies elsewhere (for I/O and diagnostics): the host copies
! of the fields are refreshed before any output, and those of the
! auxiliary arrays of the solver before writing the fields; the forcing,
! which is computed on the host, is
! copied to the device in the steps in which it changes. If we are doing
! a benchmark, we measure cputime before starting.

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
! Refreshes the host copies of the fields if any output is due
         IF (((timet.eq.tstep).or.(timec.eq.cstep) &
              .or.(dopart.and.(timep.eq.pstep))).and.(bench.eq.0)) THEN
            CALL GState_update_from(field_nxt)
         ENDIF

! Every 'tstep' steps, stores the fields
! in binary files (with the auxiliary arrays of the solver,
! which are brought to the host only here)
 BIN:    IF ((timet.eq.tstep).and.(bench.eq.0)) THEN
            timet = 0
            tind = tind+1
            CALL fluid%sync_host()
            CALL fluid%write_states(field_nxt,planio)
         ENDIF BIN

! Every 'pstep' steps, stores the particle states
         IF (dopart) THEN
            IF ((timep.eq.pstep).and.(bench.eq.0)) THEN
               timep = 0
               pind = pind+1
               CALL GPState_update_from(part_nxt) ! host copies for the I/O
               CALL particle%write_pstate((t-1)*dt,fluid,field_nxt,part_nxt)
            ENDIF
         ENDIF

! Every 'cstep' steps, generates external files
! with global quantities.

         IF ((timec.eq.cstep).and.(bench.eq.0)) THEN
            timec = 0
            IF (.not.fstatic) CALL GState_update_from(force)
            CALL fluid%global(field_nxt,force,t)
         ENDIF

! Forcing update (nothing to do for a constant forcing)
         CALL update_forcing(forcemethod,fluid,force)
         IF (forcing_changed(forcemethod)) CALL GState_update_to(force)

! Time step (GHOST Canuto stepper): field keeps the state at the
! beginning of the step and field_nxt is advanced (and likewise
! part/part_nxt when running with particles)
         gdev_active = .TRUE.
         IF (dopart) THEN
            CALL GState_copy(field,field_nxt)
            CALL GPState_copy(part,part_nxt)
            CALL stepper%gstep((t-1)*dt,field,part,force,dt,field_nxt,part_nxt)
         ELSE
            CALL GState_copy(field,field_nxt)
            CALL stepper%gstep((t-1)*dt,field,force,dt,field_nxt)
         ENDIF
         gdev_active = .FALSE.

! Per-step particle diagnostics (particle_diag.txt: global count,
! z-range, deposition and reflection counters)
         IF (dopart.and.(bench.eq.0)) THEN
            CALL particle%write_diag(part_nxt,t,t*dt)
         ENDIF

         timet = timet+1
         timep = timep+1
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
! Particle benchmark (host-side timers; the totals are reliable,
! the per-phase GPU timers are not: read the kernel stats instead)
         IF (dopart) THEN
            rbal = GetLoadBal(particle)
            IF (myrank.eq.0) THEN
               INQUIRE( file='gpbenchmark.txt', exist=bbenchexist )
               OPEN(2,file='gpbenchmark.txt',position='append')
               IF ( .NOT. bbenchexist ) THEN
                  WRITE(2,*) '# maxparts nsteps nprocs nth rbal ',  &
                     'TSTEP TCOMM TSPLINE TTRANSP TDATAEX TINTERP TPUPDATE'
               ENDIF
               WRITE(2,*) maxparts,(step-ini+1),nprocs,nth,rbal,    &
                  GetTime(particle,GPTIME_STEP)   /(step-ini+1),    &
                  GetTime(particle,GPTIME_COMM)   /(step-ini+1),    &
                  GetTime(particle,GPTIME_SPLINE) /(step-ini+1),    &
                  GetTime(particle,GPTIME_TRANSP) /(step-ini+1),    &
                  GetTime(particle,GPTIME_DATAEX) /(step-ini+1),    &
                  GetTime(particle,GPTIME_INTERP) /(step-ini+1),    &
                  GetTime(particle,GPTIME_PUPDATE)/(step-ini+1)
               CLOSE(2)
            ENDIF
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

      IF (ALLOCATED(planfc%z%dir)) CALL gdev_free(C_LOC(planfc%z%dir))
      CALL fcgram_destroy_plan(planfc)

      CALL GState_dealloc(field)
      CALL GState_dealloc(field_nxt)
      CALL GState_dealloc(force)
      IF (dopart) THEN
         CALL GPState_dealloc(part)
         CALL GPState_dealloc(part_nxt)
      ENDIF

      CALL gfree( x )
      CALL gfree( y )
      CALL gfree( z )
      CALL gfree( kx )
      CALL gfree( ky )
      CALL gfree( kz )
      CALL gfree( kk2 )
      CALL gdev_free(C_LOC(khom))
      DEALLOCATE( khom )

      END PROGRAM SPECTER
