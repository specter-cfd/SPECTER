! =====================================================================
! NAME       : inerpart_solver.f90
! DESCRIPTION: Forms solver class for inertial particles in SPECTER,
!              computing:
!
!              dx/dt   = v_p
!              dv_p/dt = (cdrag/tau) * (u(x,t) - v_p) - grav * z_hat
!                        [ + Saffman lift + Brownian force ]
!
!              with cdrag = 1 for linear Stokes drag or the nonlinear
!              drag of Clift & Gauvin (1971) (a superset of the
!              Schiller-Naumann correction used in the old SPECTER
!              particles; the extra term is negligible for Re_p < 100):
!                 cdrag = 1 + 0.15 Re_p^0.687
!                       + (Re_p/57.14)/(1 + 4.25e4 Re_p^(-1.16)),
!                 Re_p^2 = (18 tau gamma/nu) |u - v_p|^2.
!              The Cunningham slip correction divides the drag by Cc
!              (the Stokes time tau itself is kept, so the particle
!              diameter derived from it stays the physical one).
!
!              Ported to the new paradigm from SPECTER_v1 (which hacked
!              the old GHOST lagpart classes), using GHOST's
!              inerpart_solver.f90 as the template so that the GPU
!              offloading is preserved. SPECTER specifics:
!
!              - z is bounded: the physical domain spans the grid
!                coords [0, nz-Cz-1] (walls at z=0 and z=Lz, grid
!                spacing dz = Lz/(nz-Cz-1), matching the solver grid
!                z(k) = (k-1) Lz/(nz-Cz-1)). The FC-Gram continuation
!                planes are non-physical: particles never enter them.
!              - The velocity (and vorticity) fields interpolated at
!                the particles are smoothed with fc_filter in Fourier
!                space before the transform to real space (the SPECTER
!                counterpart of GHOST's spectral truncation).
!              - Wall physics, applied in end_stage after each RK
!                stage (replaces GHOST's periodic wrap in z):
!                  iwall = 0: inelastic deposition (remove particles
!                             crossing a wall; NN exchange only, VDB
!                             falls back to elastic with a warning),
!                             with optional re-injection (doreinj) at
!                             a uniformly random position in the bulk
!                             so the global count is conserved and the
!                             deposition flux reaches a steady state;
!                  iwall = 1: elastic reflection (the normal velocity
!                             reverses; at intermediate RK stages the
!                             start-of-step state is mirrored about
!                             the wall, image method);
!                  iwall = 2: stick/trap (pin at the wall with zero
!                             velocity, keep counted).
!                A particle interacts with a wall when its CENTER is
!                within rpwall = d_p/2 (grid units) of it; gamma <= 0
!                gives the point-particle rule rpwall = 0.
!              - Saffman shear lift (dolift, needs gamma > 0):
!                a_L = clift * (u - v_p) x omega / sqrt(|omega|), with
!                omega = curl(u) interpolated at the particles.
!              - Brownian (thermal) force (dobrown): white-noise
!                velocity kick applied once per full time step (final
!                RK stage), dv_i = sqrt(2 brownkT dt/tau) xi_i with
!                xi_i ~ N(0,1) (fluctuation-dissipation balance of the
!                Langevin drag).
!              - Near-wall interpolation scheme (iwallinterp, see
!                gpsplineint_mod.F90): 0 = baseline cyclic spline over
!                the FC-periodic extended grid (default, validated),
!                1 = one-sided cubic, 2 = linear in the wall cells,
!                3 = natural-BC z-spline.
!              - Per-step diagnostics in particle_diag.txt
!                (write_diag): global count, z-range, cumulative
!                deposition and reflection counters per wall. The
!                first collogmax elastic reflections per wall are
!                logged in detail to collisions_rNNNN.txt.
!
!              State ordering is:
!                x1 (x2, x3), v1 (v2, v3)
!
!              State sector ids are:
!                POSITION (POSITION+1, POSITION+2)
!                VELOCITY (VELOCITY+1, VELOCITY+2)
!
! INPUT FILE : For psolver='inerpart', looks for a "&inerpart" namelist
!              with:
!              pidir   : changes class binary input  dir (default: status idir)
!              podir   : changes class binary output dir (default: status odir)
!              tau     : Stokes time (particle response time)
!              dograv  : .false.=no gravity [default], .true.=gravity
!              grav    : gravity acceleration (z-direction, positive
!                        downward, default=0)
!              donldrag: nonlinear drag (.false.=linear [default]);
!                        when true, nu is obtained from the PDE solver
!              gamma   : mass ratio rho_f/rho_p (default=0; must be >0
!                        for the particle diameter to be defined, which
!                        the wall contact gap, the lift and the
!                        nonlinear drag need)
!              cc      : Cunningham slip correction Cc >= 1 (drag is
!                        divided by Cc; default=1, off)
!              dolift  : Saffman shear lift (.false. [default])
!              dobrown : Brownian thermal force (.false. [default])
!              brownkT : Brownian intensity k_B*T/m_p (code units of
!                        velocity^2, default=0)
!              iwall   : wall mode: 0=inelastic deposition [default],
!                        1=elastic reflection, 2=stick/trap
!              doreinj : re-inject deposited particles in the bulk
!                        (.false. [default]; only with iwall=0 and NN)
!              iwallinterp: near-wall interpolation: 0=off [default],
!                        1=one-sided cubic, 2=linear, 3=natural-BC
!              collogmax: elastic walls: log the first collogmax
!                        reflections per wall in detail (default=0)
!              partlod : particle output level of detail (default=1):
!                         1: position (xlg), part veloc (vip), fluid veloc (vlg)
!                         2: Lagrangian vorticity (wlg),
!                            strain-rate tensor (s11,s12,s13,s22,s23)
!
! DATE       : ported to SPECTER 10/27/26 (from GHOST 04/02/26, BLE)
! =====================================================================

module inerpart_mod
  use particlebase_mod
  use gpstate_mod
  use gpselect
  use gmem
  use gdevice, only: gdev_active

  implicit none

  ! ================= Solver traits ===================================
  type, public  :: InerTraits
    integer       :: partlod  = 1       ! particle output level of detail
    integer       :: iwall    = 0       ! wall mode (0=deposit,1=elastic,2=stick)
    integer       :: iwallinterp = 0    ! near-wall interpolation scheme
    integer       :: collogmax = 0      ! no. logged reflections per wall
    real(kind=GP) :: tau      = 1.0_GP  ! Stokes time
    real(kind=GP) :: grav     = 0.0_GP  ! gravity (z, positive downward)
    real(kind=GP) :: gamma    = 0.0_GP  ! mass ratio rho_f/rho_p
    real(kind=GP) :: nu       = 0.0_GP  ! fluid kinematic viscosity
    real(kind=GP) :: cc       = 1.0_GP  ! Cunningham slip correction
    real(kind=GP) :: brownkT  = 0.0_GP  ! Brownian intensity k_B*T/m_p
    real(kind=GP) :: invtau   = 1.0_GP  ! precomputed 1/(Cc*tau)
    real(kind=GP) :: nld_rep  = 0.0_GP  ! precomputed NLD const (=18 tau gamma/nu)
    real(kind=GP) :: clift    = 0.0_GP  ! Saffman lift prefactor
    real(kind=GP) :: rpwall   = 0.0_GP  ! wall contact half-gap d_p/2 (grid units)
    logical       :: dograv   = .false. ! .false.=no gravity, .true.=grav
    logical       :: donldrag = .false. ! .false.=linear drag, .true.=NLD
    logical       :: dolift   = .false. ! Saffman shear lift
    logical       :: dobrown  = .false. ! Brownian thermal force
    logical       :: doreinj  = .false. ! re-inject deposited particles
  end type

  ! ================= Solver ==========================================
  type, extends(VelocParticleBase) :: Ipart
    ! Member data:
    type (InerTraits) :: traits_
    ! Interpolated fluid vorticity at the particles (Saffman lift);
    ! CONTIGUOUS pointers into the particle pool, like lvx_ etc.
    real(kind=GP), pointer, CONTIGUOUS, dimension(:) :: lwx_ => null()
    real(kind=GP), pointer, CONTIGUOUS, dimension(:) :: lwy_ => null()
    real(kind=GP), pointer, CONTIGUOUS, dimension(:) :: lwz_ => null()
    ! Deposition bookkeeping: idrem_ (host only) accumulates the global
    ! ids removed during the RK stages of the current step; they are
    ! re-injected (if doreinj) and flushed on the final stage.
    ! ndepbot_/ndeptop_ are cumulative per-wall deposition counters and
    ! ncolbot_/ncoltop_ cumulative elastic reflection counters, split
    ! by RK stage class (1 = intermediate stages, 2 = final stage); all
    ! local to this rank (write_diag reduces them).
    integer, allocatable, dimension(:) :: idrem_
    integer           :: nrem_      = 0
    integer           :: ndepbot_   = 0, ndeptop_ = 0
    integer           :: ncolbot_(2) = 0, ncoltop_(2) = 0
    ! Collision log (elastic walls): the first collogmax reflections
    ! per wall go to collisions_rNNNN.txt; particles logged at an
    ! intermediate stage are queued in idpend_ so their end-of-step
    ! state is logged too (follow-up line, stage 9).
    integer           :: nlogbot_   = 0, nlogtop_ = 0
    integer           :: logunit_   = 0, logstep_ = 0
    integer           :: npend_     = 0
    integer           :: idpend_(64) = 0
    integer           :: diagunit_  = 0
    logical           :: vdbwarn_   = .false.
  CONTAINS
    procedure, public :: init          =>          init_impl
    procedure, public :: dpdt          =>          dpdt_impl
    procedure, public :: feedback      =>      null_feedback
    procedure, public :: write_pstate  =>  write_pstate_impl
    procedure, public :: state_size    =>    state_size_impl
    procedure, public :: part_ctor     =>         Ipart_ctor
    procedure, public :: end_stage     =>    Ipart_end_stage
    procedure, public :: write_diag    =>   Ipart_write_diag
    procedure, public :: ResizeArrays  => Ipart_ResizeArrays
    final             :: Ipart_dtor
  end type Ipart

CONTAINS

  ! ===================================================================
  ! Solver initialization, this is where parameter files are read
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to initialize the solver.
  !! Reads the &inerpart namelist and sets sector indices.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_impl(this)
    use commtypes
    class      (Ipart), intent (inout) :: this
    real     (kind=GP)                 :: tau, grav, gamma, cc, brownkT
    integer                            :: ierr, partlod
    integer                            :: iwall, iwallinterp, collogmax
    logical                            :: dograv, donldrag, dolift, dobrown
    logical                            :: doreinj
    character(len=128)                 :: pidir, podir
    namelist/ inerpart / pidir,podir,partlod,tau,grav,gamma,dograv,donldrag
    namelist/ inerpart / cc,dolift,dobrown,brownkT,iwall,doreinj
    namelist/ inerpart / iwallinterp,collogmax

    this%POSITION = 1
    this%VELOCITY = this%POSITION + this%nc_
    ! Defaults
    pidir    = this%idir_ ! Set to the pde class idir at ctor
    podir    = this%odir_ ! Set to the pde class odir at ctor
    tau      = 1.0_GP
    grav     = 0.0_GP
    gamma    = 0.0_GP
    cc       = 1.0_GP
    brownkT  = 0.0_GP
    dograv   = .false.
    donldrag = .false.
    dolift   = .false.
    dobrown  = .false.
    doreinj  = .false.
    iwall    = 0
    iwallinterp = 0
    collogmax   = 0
    partlod  = 1
    if ( this%myrank_ .eq. 0 ) then
      open(1,file=this%infile_,status='unknown',form="formatted")
      read(1,NML=inerpart)
      close(1)
    endif
    call MPI_BCAST(pidir      ,128,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(podir      ,128,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(tau        ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(grav       ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(gamma      ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(cc         ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(brownkT    ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(dograv     ,1  ,MPI_LOGICAL  ,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(donldrag   ,1  ,MPI_LOGICAL  ,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(dolift     ,1  ,MPI_LOGICAL  ,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(dobrown    ,1  ,MPI_LOGICAL  ,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(doreinj    ,1  ,MPI_LOGICAL  ,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(iwall      ,1  ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(iwallinterp,1  ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(collogmax  ,1  ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(partlod    ,1  ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)

    this%idir_ = pidir ! If present in &inerpart, replaces the class default idir
    this%odir_ = podir ! If present in &inerpart, replaces the class default odir
    this%sstate_pos_ = 'xlg' ! state name of positions
    this%sstate_lag_ = 'vlg' ! state name of Lagrangian velocities
    this%sstate_vel_ = 'vip' ! state name of particles velocities
    this%traits_%        tau = tau
    this%traits_%       grav = grav
    this%traits_%      gamma = gamma
    this%traits_%         cc = cc
    this%traits_%    brownkT = brownkT
    this%traits_%     dograv = dograv
    this%traits_%   donldrag = donldrag
    this%traits_%     dolift = dolift
    this%traits_%    dobrown = dobrown
    this%traits_%    doreinj = doreinj
    this%traits_%      iwall = iwall
    this%traits_%iwallinterp = iwallinterp
    this%traits_%  collogmax = max(collogmax,0)
    this%traits_%     invtau = 1.0_GP / (cc*tau)
    this%traits_%    partlod = partlod
    if ( (cc.ne.1.0_GP) .and. (this%myrank_.eq.0) ) then
      write(*,'(A,F8.5,A,ES12.5)') ' Inerpart: Cunningham slip Cc = ', cc, &
            ', drag relaxation time Cc*tau = ', cc*tau
    endif
  end subroutine init_impl

  ! ===================================================================
  ! Computation of RHS, the solver equations are defined here
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute the rhs of the equations of motion
  !! of inertial particles:
  !!   dx/dt   = v_p
  !!   dv_p/dt = drag [ + lift ] - gravity
  !! The fields are smoothed with fc_filter in Fourier space
  !! before being transformed and interpolated.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE dpdt_impl(this, time, pde, fluidstate, pstate, dt, dpdtout)
    use equationbase_mod
    use fft
    IMPLICIT NONE
    class       (Ipart),             intent(inout) :: this
    class(EquationBase),             intent   (in) :: pde
    real      (kind=GP),             intent   (in) :: time, dt
    type   (GStateComp), target ,    intent   (in) :: fluidstate(:)
    type  (GPStateComp), target ,    intent   (in) :: pstate(:)
    type  (GPStateComp), target ,    intent(inout) :: dpdtout(:)
    complex   (KIND=GP), pointer, dimension(:,:,:) :: velc,vcx,vcy,vcz
    real      (KIND=GP), pointer, dimension(:,:,:) :: velr,tmp1,tmp2
    real      (kind=GP)                            :: rmp, invtau, grav
    real      (kind=GP)                            :: rep2_coef, clift
    real      (kind=GP), pointer, dimension(:)     :: dpx,dpy,dpz,dvx,dvy,dvz,pvx,pvy,pvz
    logical                                        :: bret

    CALL GTStart(this%htimers_(GPTIME_STEP))
    call this%workspace_%get_complex_tmp(velc,bret)
    call this%workspace_%get_real_tmp   (velr,bret)
    call this%workspace_%get_real_tmp   (tmp1,bret)
    call this%workspace_%get_real_tmp   (tmp2,bret)
    call this%AssignLagPos(pstate) ! We assign px_,py_,pz_ to the pstate

    invtau = this%traits_%invtau
    grav   = this%traits_%grav

    select type (pde)
    class is (VelocityBase)
      if (this%nc_ .ne. pde%nc_) then
        stop "Inerpart: # of components of the particles and pdes must be equal"
      endif

      ! IMPORTANT: pstate and dpdtout may alias the same array (some steppers
      ! can pass upout for both). We must be careful about ordering:
      !   1. Interpolate fluid velocity (and vorticity) with positions intact
      !   2. Write position RHS (overwrites pstate(POSITION), but we're done)
      !   3. Compute velocity RHS (reads particle velocity before overwrite)

      rmp = 1.0_GP/(real(this%nd_(1),kind=GP)*real(this%nd_(2),kind=GP)* &
                    real(this%nd_(3),kind=GP))
      vcx => fluidstate(pde%VELOCITY  )%ccomp
      vcy => fluidstate(pde%VELOCITY+1)%ccomp
      vcz => fluidstate(pde%VELOCITY+2)%ccomp

      ! Step 1: Interpolate fluid velocity to lvx_, lvy_, lvz_ (smoothed
      ! with fc_filter before the transform, as in the old SPECTER code)
      CALL ipart_copy_scal(vcx,velc,rmp)
      CALL fc_filter(velc)
      call fftp3d_complex_to_real(planfc,velc,velr,MPI_COMM_WORLD)
      call this%EulerToLag(this%lvx_,this%nparts_,velr,.true. ,tmp1,tmp2)
      CALL ipart_copy_scal(vcy,velc,rmp)
      CALL fc_filter(velc)
      call fftp3d_complex_to_real(planfc,velc,velr,MPI_COMM_WORLD)
      call this%EulerToLag(this%lvy_,this%nparts_,velr,.false.,tmp1,tmp2)
      CALL ipart_copy_scal(vcz,velc,rmp)
      CALL fc_filter(velc)
      call fftp3d_complex_to_real(planfc,velc,velr,MPI_COMM_WORLD)
      call this%EulerToLag(this%lvz_,this%nparts_,velr,.false.,tmp1,tmp2)

      ! Saffman lift: interpolate omega = curl(u) at the particles
      ! (same stencils as the velocity: doupdate=.false.). curlk with
      ! SPECTER's derivk handles the non-periodic z via FC-Gram, so the
      ! near-wall shear that drives the lift is captured correctly.
      if ( this%traits_%dolift .and. this%traits_%clift.gt.0.0_GP ) then
        CALL curlk(vcy,vcz,velc,1)               ! omega_x = dy(vz)-dz(vy)
        CALL ipart_scal(velc,rmp)
        CALL fc_filter(velc)
        call fftp3d_complex_to_real(planfc,velc,velr,MPI_COMM_WORLD)
        call this%EulerToLag(this%lwx_,this%nparts_,velr,.false.,tmp1,tmp2)
        CALL curlk(vcx,vcz,velc,2)               ! omega_y = dz(vx)-dx(vz)
        CALL ipart_scal(velc,rmp)
        CALL fc_filter(velc)
        call fftp3d_complex_to_real(planfc,velc,velr,MPI_COMM_WORLD)
        call this%EulerToLag(this%lwy_,this%nparts_,velr,.false.,tmp1,tmp2)
        CALL curlk(vcx,vcy,velc,3)               ! omega_z = dx(vy)-dy(vx)
        CALL ipart_scal(velc,rmp)
        CALL fc_filter(velc)
        call fftp3d_complex_to_real(planfc,velc,velr,MPI_COMM_WORLD)
        call this%EulerToLag(this%lwz_,this%nparts_,velr,.false.,tmp1,tmp2)
      endif

      rep2_coef = this%traits_%nld_rep
      clift     = this%traits_%clift
      ! Steps 2 and 3 in one kernel over the particles: position RHS
      ! dx/dt = v_p and velocity RHS dv_p/dt = cdrag/tau (u - v_p)
      ! [+ lift] - g. pstate and dpdtout may alias (some steppers pass
      ! upout for both): each particle reads its velocity before
      ! writing its RHS. Pointers and a kernel (ipart_rhs) are used to
      ! help offloading.
      dpx => dpdtout(this%POSITION  )%rcomp
      dpy => dpdtout(this%POSITION+1)%rcomp
      dpz => dpdtout(this%POSITION+2)%rcomp
      dvx => dpdtout(this%VELOCITY  )%rcomp
      dvy => dpdtout(this%VELOCITY+1)%rcomp
      dvz => dpdtout(this%VELOCITY+2)%rcomp
      pvx => pstate (this%VELOCITY  )%rcomp
      pvy => pstate (this%VELOCITY+1)%rcomp
      pvz => pstate (this%VELOCITY+2)%rcomp
      call ipart_rhs(this%nparts_,this%lvx_,this%lvy_,this%lvz_,          &
                     this%lwx_,this%lwy_,this%lwz_,pvx,pvy,pvz,           &
                     dpx,dpy,dpz,dvx,dvy,dvz,this%invdel_,invtau,grav,    &
                     rep2_coef,clift,this%traits_%donldrag,               &
                     this%traits_%dograv,                                 &
                     this%traits_%dolift .and. clift.gt.0.0_GP)
    class default
      stop "Inerpart: This solver does not support pdes without a velocity field"
    end select

    call this%workspace_%free_complex_tmp(velc)
    call this%workspace_%free_real_tmp   (velr)
    call this%workspace_%free_real_tmp   (tmp1)
    call this%workspace_%free_real_tmp   (tmp2)
    CALL GTAcc(this%htimers_(GPTIME_STEP))
  END SUBROUTINE dpdt_impl


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Functions to compute fluid and particle couplings
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine null_feedback(this, pstate, feedback)
    use gstate_mod
    class     (Ipart), intent   (in) :: this
    type(GPStateComp), intent   (in) :: pstate(:)
    type (GStateComp), intent(inout) :: feedback(:)
    return
  end subroutine null_feedback


  ! ===================================================================
  ! End of stage: Brownian kick, wall physics, re-injection, and the
  ! default synchronization/exchange of the base class
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : Ipart_end_stage
  !!  DESCRIPTION: SPECTER wall physics (replaces GHOST's periodic
  !!               wrap in z), applied after every RK stage before the
  !!               particle exchange:
  !!               - Brownian velocity kick (final stage only),
  !!               - wall interaction (stick / elastic / deposition),
  !!               - re-injection of the deposited particles (final
  !!                 stage only),
  !!               then the base end_stage (x-y periodicity, NN or VDB
  !!               exchange, buffer resizes). The physical z-range in
  !!               grid units is [wbot, wtop] with wbot = rpwall and
  !!               wtop = (nz-Cz-1) - rpwall (finite-size contact).
  !!  ARGUMENTS  :
  !!    upin    : state at the beginning of the step
  !!    upout   : state after the stage
  !!    stage   : RK stage (norder, ..., 1; 1 is the final stage)
  !!    norder  : number of RK stages
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Ipart_end_stage(this, upin, upout, stage, norder)
    use random
    use status, only: dt
    IMPLICIT NONE
    CLASS       (Ipart), INTENT(INOUT)         :: this
    TYPE  (GPStateComp), INTENT(INOUT), TARGET :: upin (:)
    TYPE  (GPStateComp), INTENT(INOUT), TARGET :: upout(:)
    INTEGER            , INTENT(IN)            :: stage, norder
    REAL(KIND=GP), pointer, dimension(:) :: pxo,pyo,pzo,pvxo,pvyo,pvzo
    REAL(KIND=GP), pointer, dimension(:) :: pxi,pyi,pzi,pvxi,pvyi,pvzi
    REAL(KIND=GP)   :: wbot,wtop,bamp,rpi,u1,u2
    INTEGER         :: j,istg,interm,ncb,nct,nb,nt,nrem,nkeep
    LOGICAL         :: dolog

    pxo  => upout(this%POSITION  )%rcomp
    pyo  => upout(this%POSITION+1)%rcomp
    pzo  => upout(this%POSITION+2)%rcomp
    pvxo => upout(this%VELOCITY  )%rcomp
    pvyo => upout(this%VELOCITY+1)%rcomp
    pvzo => upout(this%VELOCITY+2)%rcomp
    pxi  => upin (this%POSITION  )%rcomp
    pyi  => upin (this%POSITION+1)%rcomp
    pzi  => upin (this%POSITION+2)%rcomp
    pvxi => upin (this%VELOCITY  )%rcomp
    pvyi => upin (this%VELOCITY+1)%rcomp
    pvzi => upin (this%VELOCITY+2)%rcomp

    wbot = this%traits_%rpwall
    wtop = this%gext_(3) - this%traits_%rpwall

    ! Step counter for the collision log (self-contained: counts the
    ! steps taken in this run, starting at 1)
    IF ( stage .EQ. norder ) this%logstep_ = this%logstep_ + 1

    ! --- Brownian (thermal) force: white-noise velocity kick, operator-
    ! split from the deterministic RK and applied ONCE per full time
    ! step (the final RK stage). The kick std follows the fluctuation-
    ! dissipation balance of the Langevin drag, <dv_i^2> = 2(kT/m)dt/tau.
    ! The normal draws are generated on the host (per-rank stream, in
    ! particle order, deterministic) into the host copy of gptmp0_,
    ! uploaded, and added by a device kernel.
    IF ( (stage.EQ.1) .AND. this%traits_%dobrown .AND. &
         (this%traits_%brownkT.GT.0.0_GP) ) THEN
      rpi  = 4.0_GP*ATAN(1.0_GP)
      ! The kick uses the full time step dt (status module): it is
      ! applied once per step, at the final stage
      bamp = SQRT(2.0_GP*this%traits_%brownkT*dt*this%traits_%invtau)
      DO j = 1, this%nparts_
        CALL prandom_number(u1)
        CALL prandom_number(u2)
        u1 = MAX(u1,1.0e-30_GP)
        this%gptmp0_(1,j) = SQRT(-2.0_GP*LOG(u1))*COS(2.0_GP*rpi*u2)
        CALL prandom_number(u1)
        CALL prandom_number(u2)
        u1 = MAX(u1,1.0e-30_GP)
        this%gptmp0_(2,j) = SQRT(-2.0_GP*LOG(u1))*COS(2.0_GP*rpi*u2)
        CALL prandom_number(u1)
        CALL prandom_number(u2)
        u1 = MAX(u1,1.0e-30_GP)
        this%gptmp0_(3,j) = SQRT(-2.0_GP*LOG(u1))*COS(2.0_GP*rpi*u2)
      ENDDO
      CALL gupdate_to(this%gptmp0_)
      CALL ipart_brown(this%nparts_,bamp,this%gptmp0_,pvxo,pvyo,pvzo, &
                       this%traits_%iwall,wbot,wtop,pzo)
    ENDIF

    ! --- Wall physics --------------------------------------------------
    IF ( this%traits_%iwall .EQ. 2 ) THEN
      ! Stick / trap: pin to the wall with zero velocity, in upout AND
      ! in the start-of-step state upin, so the remaining stages (which
      ! rebuild from upin) keep the particle pinned. The particle stays
      ! in the active list (count conserved).
      CALL ipart_stick(this%nparts_,wbot,wtop,pxo,pyo,pzo,pvxo,pvyo,pvzo, &
                       pxi,pyi,pzi,pvxi,pvyi,pvzi)

    ELSE IF ( this%traits_%iwall .EQ. 1 ) THEN
      ! Elastic reflection: reverse the normal velocity. A reflection at
      ! an intermediate stage (stage > 1) also mirrors the start-of-step
      ! state upin about the wall (image method): the final stage
      ! rebuilds z and vz from upin, and without the mirror the
      ! end-of-step vz would still point into the wall. Reflections are
      ! counted per wall and stage class; tmpint_ flags the reflected
      ! particles (0 none, 1 bottom, 2 top) for the optional log.
      interm = 0
      istg   = 2
      IF ( stage.GT.1 ) THEN
        interm = 1
        istg   = 1
      ENDIF
      CALL ipart_reflect(this%nparts_,wbot,wtop,interm,pzo,pvzo,pzi,pvzi, &
                         this%tmpint_,ncb,nct)
      this%ncolbot_(istg) = this%ncolbot_(istg) + ncb
      this%ncoltop_(istg) = this%ncoltop_(istg) + nct
      dolog = (this%traits_%collogmax.GT.0) .AND. ((ncb.GT.0).OR.(nct.GT.0)) &
              .AND. ((this%nlogbot_.LT.this%traits_%collogmax) .OR.          &
                     (this%nlogtop_.LT.this%traits_%collogmax))
      IF ( dolog .OR. ((stage.EQ.1).AND.(this%npend_.GT.0)) ) THEN
        CALL Ipart_LogCol(this,upout,istg,wbot,wtop,dolog)
      ENDIF

    ELSE
      ! Inelastic deposition (default): remove the particles crossing a
      ! wall. Requires the NN exchange (in-place compaction); VDB cannot
      ! cleanly remove by global id, so it falls back to elastic
      ! reflection with a warning.
      IF ( this%iexchtype_ .EQ. GPEXCHTYPE_NN ) THEN
        CALL ipart_depflag(this%nparts_,wbot,wtop,pzo,this%tmpint_,nb,nt)
        IF ( nb+nt .GT. 0 ) THEN
          this%ndepbot_ = this%ndepbot_ + nb
          this%ndeptop_ = this%ndeptop_ + nt
          ! Record the removed global ids (used by the re-injection at
          ! the end of the step): compact the inverted flags, gather the
          ! ids at those indices and copy them to the host list.
          IF ( this%traits_%doreinj ) THEN
            CALL ipart_invflag(this%nparts_,this%tmpint_)
            CALL gpsel_compact(1,this%nparts_,this%tmpint_,this%iwrk_,nrem,0)
            CALL ipart_gather_i(nrem,this%iwrk_,this%id_,this%tmpint_)
            CALL gupdate_from(this%tmpint_)
            DO j = 1, nrem
              IF ( this%nrem_ .LT. SIZE(this%idrem_) ) THEN
                this%nrem_ = this%nrem_ + 1
                this%idrem_(this%nrem_) = this%tmpint_(j)
              ENDIF
            ENDDO
            ! Recompute the keep flags (tmpint_ was overwritten)
            CALL ipart_depflag(this%nparts_,wbot,wtop,pzo,this%tmpint_,nb,nt)
          ENDIF
          ! Deterministic compaction of all the per-particle arrays
          ! through the ascending keep list (lvx_ and tmpint_ are used
          ! as scratch; the Lagrangian velocity/vorticity temporaries
          ! are recomputed from scratch at the next dpdt call)
          CALL gpsel_compact(1,this%nparts_,this%tmpint_,this%iwrk_,nkeep,0)
          CALL ipart_compact_r(nkeep,this%iwrk_,pxo ,this%lvx_)
          CALL ipart_compact_r(nkeep,this%iwrk_,pyo ,this%lvx_)
          CALL ipart_compact_r(nkeep,this%iwrk_,pzo ,this%lvx_)
          CALL ipart_compact_r(nkeep,this%iwrk_,pvxo,this%lvx_)
          CALL ipart_compact_r(nkeep,this%iwrk_,pvyo,this%lvx_)
          CALL ipart_compact_r(nkeep,this%iwrk_,pvzo,this%lvx_)
          CALL ipart_compact_r(nkeep,this%iwrk_,pxi ,this%lvx_)
          CALL ipart_compact_r(nkeep,this%iwrk_,pyi ,this%lvx_)
          CALL ipart_compact_r(nkeep,this%iwrk_,pzi ,this%lvx_)
          CALL ipart_compact_r(nkeep,this%iwrk_,pvxi,this%lvx_)
          CALL ipart_compact_r(nkeep,this%iwrk_,pvyi,this%lvx_)
          CALL ipart_compact_r(nkeep,this%iwrk_,pvzi,this%lvx_)
          CALL ipart_compact_i(nkeep,this%iwrk_,this%id_,this%tmpint_)
          this%nparts_ = nkeep
        ENDIF
        ! Flush the deposition record once per full time step (final
        ! stage): re-inject the removed particles in the bulk if
        ! requested, then clear.
        IF ( stage .EQ. 1 ) THEN
          ! Ipart_ReInject is collective: every rank calls it, even the
          ! ones that deposited nothing this step
          IF ( this%traits_%doreinj ) THEN
            CALL Ipart_ReInject(this,upin,upout,wbot,wtop)
          ENDIF
          this%nrem_ = 0
        ENDIF
      ELSE
        IF ( (.NOT.this%vdbwarn_) .AND. (this%myrank_.EQ.0) ) THEN
          WRITE(*,'(A)') ' SPECTER WARNING (Inerpart end_stage): VDB exchange '//&
            'mode does not support inelastic deposition; applying elastic '  //&
            'reflection. Set iwall=1 or switch to NN exchange.'
          this%vdbwarn_ = .TRUE.
        ENDIF
        CALL ipart_reflect(this%nparts_,wbot,wtop,0,pzo,pvzo,pzi,pvzi, &
                           this%tmpint_,ncb,nct)
      ENDIF
    ENDIF

    ! --- Default synchronization (x-y periodicity, exchange, resizes) ---
    CALL end_stage(this, upin, upout, stage, norder)
  END SUBROUTINE Ipart_end_stage


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! The full time step of the current step, needed by the
  !! Brownian kick (the stage dt is not passed to end_stage).
  !! The Canuto stepper applies eff_dt = dt at the final stage,
  !! and the kick is applied only there, so the step dt is the
  !! status module dt.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function dtfull_(this) result(dtf)
    use status, only: dt
    class(Ipart), intent(in) :: this
    real(kind=GP)            :: dtf
    dtf = dt
  end function dtfull_


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : Ipart_ReInject
  !!  DESCRIPTION: Re-seeds the particles deposited during the current
  !!               time step at uniformly random positions in the bulk,
  !!               so the GLOBAL particle count is exactly conserved and
  !!               the deposition flux reaches a statistical steady
  !!               state. Each particle keeps its global id. Rank 0
  !!               draws every replacement position and broadcasts it;
  !!               each rank appends only the replacements that fall
  !!               inside its own slab (the NN exchange only moves
  !!               particles between neighboring slabs, so the particle
  !!               must be CREATED on its owner rank). Re-injected
  !!               particles start at rest. Runs on the host: the host
  !!               copies of the states are refreshed, appended to, and
  !!               uploaded (only on the ranks that append).
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Ipart_ReInject(this,upin,upout,wbot,wtop)
    use random
    use commtypes
    IMPLICIT NONE
    CLASS       (Ipart), INTENT(INOUT) :: this
    TYPE  (GPStateComp), INTENT(INOUT) :: upin (:)
    TYPE  (GPStateComp), INTENT(INOUT) :: upout(:)
    REAL(KIND=GP)      , INTENT(IN)    :: wbot,wtop
    INTEGER                                  :: i,j,ntot,nnew
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ncnt,ndsp,idg,itmp
    REAL(KIND=GP), ALLOCATABLE, DIMENSION(:) :: pg
    REAL(KIND=GP)                            :: r,zlo,zhi

    ALLOCATE(ncnt(this%nprocs_),ndsp(this%nprocs_))
    ! How many particles did each rank deposit this time step?
    CALL MPI_ALLGATHER(this%nrem_,1,MPI_INTEGER,ncnt,1,MPI_INTEGER, &
                       this%comm_,this%ierr_)
    ntot = 0
    DO i = 1, this%nprocs_
      ndsp(i) = ntot
      ntot    = ntot + ncnt(i)
    ENDDO
    IF ( ntot .EQ. 0 ) THEN
      DEALLOCATE(ncnt,ndsp)
      RETURN
    ENDIF

    ALLOCATE(idg(ntot),pg(3*ntot))
    ! Collect the global ids of every particle deposited this step:
    CALL MPI_ALLGATHERV(this%idrem_,this%nrem_,MPI_INTEGER,idg,ncnt,ndsp, &
                        MPI_INTEGER,this%comm_,this%ierr_)

    ! Rank 0 draws all replacement positions, then broadcasts them, so
    ! every rank agrees on where each id is re-injected.
    IF ( this%myrank_ .EQ. 0 ) THEN
      DO i = 1, ntot
        CALL prandom_number(r)
        pg(3*i-2) = r*this%gext_(1)             ! x: uniform, periodic
        CALL prandom_number(r)
        pg(3*i-1) = r*this%gext_(2)             ! y: uniform, periodic
        CALL prandom_number(r)
        pg(3*i  ) = wbot + r*(wtop-wbot)        ! z: uniform in the bulk
      ENDDO
    ENDIF
    CALL MPI_BCAST(pg,3*ntot,GC_REAL,0,this%comm_,this%ierr_)

    ! Each rank appends the replacements that land inside its own slab.
    ! Ownership matches the NN exchange test: lxbnds_(3,1) <= z < lxbnds_(3,2).
    zlo = this%lxbnds_(3,1)
    zhi = this%lxbnds_(3,2)
    nnew = 0
    DO i = 1, ntot
      IF ( pg(3*i).GE.zlo .AND. pg(3*i).LT.zhi ) nnew = nnew + 1
    ENDDO

    IF ( nnew .GT. 0 ) THEN
      ! Grow the buffers if the appended particles would overflow them
      IF ( this%nparts_+nnew .GT. this%partbuff_ ) THEN
        this%partbuff_ = this%partbuff_ + (1 + (this%nparts_+nnew - &
                         this%partbuff_)/this%partchunksize_)*this%partchunksize_
        CALL this%ResizeArrays(this%partbuff_,.true.)
        CALL GPState_resize(upin ,this%partbuff_)
        CALL GPState_resize(upout,this%partbuff_)
        IF ( SIZE(this%idrem_) .LT. this%partbuff_ ) THEN
          ! The removal record keeps its contents (it holds the ids
          ! deposited during this very step)
          ALLOCATE(itmp(this%partbuff_))
          itmp = 0
          itmp(1:SIZE(this%idrem_)) = this%idrem_
          CALL MOVE_ALLOC(itmp,this%idrem_)
        ENDIF
      ENDIF
      ! Append on the host: refresh the host copies, append, upload
      CALL GPState_update_from(upout)
      CALL GPState_update_from(upin)
      CALL gupdate_from(this%id_)
      DO i = 1, ntot
        IF ( pg(3*i).GE.zlo .AND. pg(3*i).LT.zhi ) THEN
          j = this%nparts_ + 1
          upout(this%POSITION  )%rcomp(j) = pg(3*i-2)
          upout(this%POSITION+1)%rcomp(j) = pg(3*i-1)
          upout(this%POSITION+2)%rcomp(j) = pg(3*i  )
          upout(this%VELOCITY  )%rcomp(j) = 0.0_GP
          upout(this%VELOCITY+1)%rcomp(j) = 0.0_GP
          upout(this%VELOCITY+2)%rcomp(j) = 0.0_GP
          upin (this%POSITION  )%rcomp(j) = pg(3*i-2)
          upin (this%POSITION+1)%rcomp(j) = pg(3*i-1)
          upin (this%POSITION+2)%rcomp(j) = pg(3*i  )
          upin (this%VELOCITY  )%rcomp(j) = 0.0_GP
          upin (this%VELOCITY+1)%rcomp(j) = 0.0_GP
          upin (this%VELOCITY+2)%rcomp(j) = 0.0_GP
          this%id_(j)  = idg(i)
          this%nparts_ = j
        ENDIF
      ENDDO
      CALL GPState_update_to(upout)
      CALL GPState_update_to(upin)
      CALL gupdate_to(this%id_)
    ENDIF
    DEALLOCATE(ncnt,ndsp,idg,pg)
  END SUBROUTINE Ipart_ReInject


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : Ipart_LogCol
  !!  DESCRIPTION: Detailed log of the elastic wall reflections of
  !!               this stage: appends one line per reflected particle
  !!               to this rank's collisions_rNNNN.txt (opened on
  !!               first use), physical units:
  !!                 step stage wall id x y z_before z_after
  !!                 vx vy vz_before vz_after
  !!               stage 1/2 = reflection applied in an intermediate/
  !!               final RK stage; stage 9 = end-of-step state of a
  !!               particle logged at stage 1 (follow-up, wall = -1).
  !!               The pre-reflection z and vz are reconstructed from
  !!               the reflected values (z_old = 2 w - z, vz_old =
  !!               -vz). Runs on the host: the host copies of the
  !!               flags, ids, positions and velocities are refreshed
  !!               (debugging aid; only when collogmax > 0 and a
  !!               reflection happened).
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Ipart_LogCol(this,upout,istg,wbot,wtop,dolog)
    IMPLICIT NONE
    CLASS       (Ipart), INTENT(INOUT) :: this
    TYPE  (GPStateComp), INTENT(INOUT) :: upout(:)
    INTEGER            , INTENT(IN)    :: istg
    REAL(KIND=GP)      , INTENT(IN)    :: wbot,wtop
    LOGICAL            , INTENT(IN)    :: dolog
    REAL(KIND=GP)   :: zold,vold,w
    INTEGER         :: j,k,iwl
    CHARACTER(len=32) :: fname

    ! Refresh the host copies this log reads
    CALL gupdate_from(this%tmpint_)
    CALL gupdate_from(this%id_)
    CALL gupdate_from(upout(this%POSITION  )%rcomp)
    CALL gupdate_from(upout(this%POSITION+1)%rcomp)
    CALL gupdate_from(upout(this%POSITION+2)%rcomp)
    CALL gupdate_from(upout(this%VELOCITY  )%rcomp)
    CALL gupdate_from(upout(this%VELOCITY+1)%rcomp)
    CALL gupdate_from(upout(this%VELOCITY+2)%rcomp)

    IF ( this%logunit_ .EQ. 0 ) THEN
      WRITE(fname,'(a,i4.4,a)') 'collisions_r',this%myrank_,'.txt'
      OPEN(NEWUNIT=this%logunit_,FILE=TRIM(fname),STATUS='replace',ACTION='write')
      WRITE(this%logunit_,'(a)') '# step stage wall id x y z_before z_after vx vy vz_before vz_after' // &
        '   (stage 1/2: reflection in intermediate/final RK stage; 9: end-of-step state of a stage-1 particle)'
    ENDIF

    IF ( dolog ) THEN
      DO j = 1, this%nparts_
        IF ( this%tmpint_(j) .EQ. 0 ) CYCLE
        iwl = this%tmpint_(j) - 1       ! 0 = bottom, 1 = top
        IF ( iwl.EQ.0 ) THEN
          IF ( this%nlogbot_ .GE. this%traits_%collogmax ) CYCLE
          this%nlogbot_ = this%nlogbot_ + 1
          w = wbot
        ELSE
          IF ( this%nlogtop_ .GE. this%traits_%collogmax ) CYCLE
          this%nlogtop_ = this%nlogtop_ + 1
          w = wtop
        ENDIF
        zold = 2.0_GP*w - upout(this%POSITION+2)%rcomp(j)
        vold = -upout(this%VELOCITY+2)%rcomp(j)
        WRITE(this%logunit_,'(i9,i3,i3,i10,8es16.8)') this%logstep_,istg,iwl, &
          this%id_(j),                                                        &
          upout(this%POSITION  )%rcomp(j)*this%delta_(1),                     &
          upout(this%POSITION+1)%rcomp(j)*this%delta_(2),                     &
          zold*this%delta_(3),                                                &
          upout(this%POSITION+2)%rcomp(j)*this%delta_(3),                     &
          upout(this%VELOCITY  )%rcomp(j),upout(this%VELOCITY+1)%rcomp(j),    &
          vold,upout(this%VELOCITY+2)%rcomp(j)
        FLUSH(this%logunit_)
        IF ( istg.EQ.1 .AND. this%npend_.LT.SIZE(this%idpend_) ) THEN
          this%npend_ = this%npend_ + 1
          this%idpend_(this%npend_) = this%id_(j)
        ENDIF
      ENDDO
    ENDIF

    ! End-of-step state of the particles logged at an intermediate
    ! stage (follow-up line; particles that changed MPI slab during the
    ! step are not found)
    IF ( istg.EQ.2 .AND. this%npend_.GT.0 ) THEN
      DO j = 1, this%nparts_
        DO k = 1, this%npend_
          IF ( this%id_(j) .EQ. this%idpend_(k) ) THEN
            WRITE(this%logunit_,'(i9,i3,i3,i10,8es16.8)') this%logstep_,9,-1, &
              this%id_(j),                                                    &
              upout(this%POSITION  )%rcomp(j)*this%delta_(1),                 &
              upout(this%POSITION+1)%rcomp(j)*this%delta_(2),                 &
              upout(this%POSITION+2)%rcomp(j)*this%delta_(3),                 &
              upout(this%POSITION+2)%rcomp(j)*this%delta_(3),                 &
              upout(this%VELOCITY  )%rcomp(j),upout(this%VELOCITY+1)%rcomp(j),&
              upout(this%VELOCITY+2)%rcomp(j),upout(this%VELOCITY+2)%rcomp(j)
            FLUSH(this%logunit_)
          ENDIF
        ENDDO
      ENDDO
      this%npend_ = 0
    ENDIF
  END SUBROUTINE Ipart_LogCol


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : Ipart_write_diag
  !!  DESCRIPTION: Per-step particle diagnostics to
  !!               particle_diag.txt: global active count, global
  !!               particle z-range (physical units), cumulative
  !!               deposition counts per wall, and cumulative elastic
  !!               reflection counts per wall and RK stage class.
  !!               Called by the main program at the end of every
  !!               step; only the z-positions are brought to the host.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Ipart_write_diag(this, pstate, t, time)
    use commtypes
    IMPLICIT NONE
    CLASS       (Ipart), INTENT(INOUT)         :: this
    TYPE  (GPStateComp), INTENT(INOUT), TARGET :: pstate(:)
    INTEGER            , INTENT(IN)            :: t
    REAL(KIND=GP)      , INTENT(IN)            :: time
    REAL(KIND=GP)   :: zl(2),zg(2),zmin,zmax
    INTEGER         :: j,np,ndb,ndt,ncb(2),nct(2)

    CALL gupdate_from(pstate(this%POSITION+2)%rcomp)
    zl(1) =  HUGE(1.0_GP)
    zl(2) = -HUGE(1.0_GP)
    DO j = 1, this%nparts_
      zl(1) = MIN(zl(1),pstate(this%POSITION+2)%rcomp(j))
      zl(2) = MAX(zl(2),pstate(this%POSITION+2)%rcomp(j))
    ENDDO
    zl(2) = -zl(2)
    CALL MPI_ALLREDUCE(zl,zg,2,GC_REAL,MPI_MIN,this%comm_,this%ierr_)
    zmin =  zg(1)*this%delta_(3)
    zmax = -zg(2)*this%delta_(3)
    CALL MPI_ALLREDUCE(this%nparts_ ,np ,1,MPI_INTEGER,MPI_SUM,this%comm_,this%ierr_)
    CALL MPI_ALLREDUCE(this%ndepbot_,ndb,1,MPI_INTEGER,MPI_SUM,this%comm_,this%ierr_)
    CALL MPI_ALLREDUCE(this%ndeptop_,ndt,1,MPI_INTEGER,MPI_SUM,this%comm_,this%ierr_)
    CALL MPI_ALLREDUCE(this%ncolbot_,ncb,2,MPI_INTEGER,MPI_SUM,this%comm_,this%ierr_)
    CALL MPI_ALLREDUCE(this%ncoltop_,nct,2,MPI_INTEGER,MPI_SUM,this%comm_,this%ierr_)
    IF ( this%myrank_ .EQ. 0 ) THEN
      IF ( this%diagunit_ .EQ. 0 ) THEN
        OPEN(NEWUNIT=this%diagunit_,FILE='particle_diag.txt',POSITION='append',&
             ACTION='write')
        WRITE(this%diagunit_,'(a)') '# step  t  nparts  zmin  zmax  '       // &
          'ndep_bot ndep_top  (cumulative deposition)  '                    // &
          'ncol_bot ncol_top  (cumulative elastic, final RK stage)  '       // &
          'ncol_bot_s1 ncol_top_s1  (intermediate stage)'
      ENDIF
      WRITE(this%diagunit_,'(i9,es16.8,i12,2es18.10,6i12)') t,time,np,       &
            zmin,zmax,ndb,ndt,ncb(2),nct(2),ncb(1),nct(1)
      FLUSH(this%diagunit_)
    ENDIF
  END SUBROUTINE Ipart_write_diag


  ! ===================================================================
  ! Internal routines: RHS and wall-physics kernels (module procedures
  ! with explicit-shape arrays; device kernels while gdev_active)
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! b = a*rmp (complex field copy with FFT normalization)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ipart_copy_scal(a,b,rmp)
    use grid
    use mpivars
    implicit none
    complex(kind=GP), intent (in), dimension(nz,ny,ista:iend) :: a
    complex(kind=GP), intent(out), dimension(nz,ny,ista:iend) :: b
    real   (kind=GP), intent (in)                             :: rmp
    integer :: i,j,k
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (k)
#endif
    do i = ista,iend
      do j = 1,ny
        do k = 1,nz
          b(k,j,i) = a(k,j,i)*rmp
        end do
      end do
    end do
  end subroutine ipart_copy_scal


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! a = a*rmp (in-place complex field scaling)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ipart_scal(a,rmp)
    use grid
    use mpivars
    implicit none
    complex(kind=GP), intent(inout), dimension(nz,ny,ista:iend) :: a
    real   (kind=GP), intent   (in)                             :: rmp
    integer :: i,j,k
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (k)
#endif
    do i = ista,iend
      do j = 1,ny
        do k = 1,nz
          a(k,j,i) = a(k,j,i)*rmp
        end do
      end do
    end do
  end subroutine ipart_scal


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! a = s*(a + b) (combination used for the strain components)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ipart_combine(a,b,s)
    use grid
    use mpivars
    implicit none
    complex(kind=GP), intent(inout), dimension(nz,ny,ista:iend) :: a
    complex(kind=GP), intent   (in), dimension(nz,ny,ista:iend) :: b
    real   (kind=GP), intent   (in)                             :: s
    integer :: i,j,k
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (k)
#endif
    do i = ista,iend
      do j = 1,ny
        do k = 1,nz
          a(k,j,i) = s*(a(k,j,i) + b(k,j,i))
        end do
      end do
    end do
  end subroutine ipart_combine


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Internal kernel to compute right-hand side of n particles
  !!   dx/dt = v_p/delta         (positions in grid units)
  !!   dv/dt = cdrag/tau (u - v_p) [+ lift] - grav z_hat
  !! with cdrag = 1 for linear Stokes drag, or the nonlinear
  !! drag of Clift & Gauvin (1971):
  !!   cdrag = 1 + 0.15 Re_p^0.687 + Re_p^2.16/(57.14 (Re_p^1.16+4.25e4))
  !!   Re_p^2 = rep2_coef |u - v_p|^2, rep2_coef = 18 tau gamma / nu
  !! and, with dolift, the Saffman shear lift
  !!   a_L = clift (u - v_p) x omega / sqrt(|omega|),
  !! symmetric between the two walls with no hand-coded sign: the
  !! mean spanwise vorticity reverses across the channel, so a
  !! wall-lagging particle is pushed away from whichever wall it
  !! sits near.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ipart_rhs(n,lvx,lvy,lvz,lwx,lwy,lwz,pvx,pvy,pvz,            &
                       dpx,dpy,dpz,dvx,dvy,dvz,invdel,invtau,grav,       &
                       rep2_coef,clift,donldrag,dograv,dolift)
    implicit none
    integer      , intent(in)    :: n
    real(kind=GP), intent(in)    :: lvx(n),lvy(n),lvz(n)
    real(kind=GP), intent(in)    :: lwx(n),lwy(n),lwz(n)
    real(kind=GP), intent(in)    :: pvx(n),pvy(n),pvz(n)
    real(kind=GP), intent(inout) :: dpx(n),dpy(n),dpz(n),dvx(n),dvy(n),dvz(n)
    real(kind=GP), intent(in)    :: invdel(3),invtau,grav,rep2_coef,clift
    logical      , intent(in)    :: donldrag,dograv,dolift
    real(kind=GP)                :: vx,vy,vz,dx,dy,dz,rep2,cdrag
    real(kind=GP)                :: tx,ty,tz,wmag,cf
    integer                      :: j
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active) &
!$omp   private(vx,vy,vz,dx,dy,dz,rep2,cdrag,tx,ty,tz,wmag,cf)
#else
!$omp parallel do private(vx,vy,vz,dx,dy,dz,rep2,cdrag,tx,ty,tz,wmag,cf)
#endif
    do j = 1,n
      vx = pvx(j); vy = pvy(j); vz = pvz(j)
      dpx(j) = vx*invdel(1)
      dpy(j) = vy*invdel(2)
      dpz(j) = vz*invdel(3)
      dx = lvx(j) - vx
      dy = lvy(j) - vy
      dz = lvz(j) - vz
      if ( donldrag ) then
        rep2 = rep2_coef * (dx**2 + dy**2 + dz**2)
        cdrag = 1.0_GP + 0.15_GP * rep2 ** 0.3435_GP                  &
              + rep2 ** 1.08_GP / (57.14_GP * (rep2 ** 0.58_GP + 4.25e4_GP))
        tx = dx * cdrag * invtau
        ty = dy * cdrag * invtau
        tz = dz * cdrag * invtau
      else
        tx = dx * invtau
        ty = dy * invtau
        tz = dz * invtau
      endif
      if ( dolift ) then
        wmag = SQRT(lwx(j)**2 + lwy(j)**2 + lwz(j)**2)
        if ( wmag .gt. 1.0e-12_GP ) then
          cf = clift/SQRT(wmag)
          tx = tx + cf*(dy*lwz(j) - dz*lwy(j))
          ty = ty + cf*(dz*lwx(j) - dx*lwz(j))
          tz = tz + cf*(dx*lwy(j) - dy*lwx(j))
        endif
      endif
      if ( dograv ) tz = tz - grav
      dvx(j) = tx
      dvy(j) = ty
      dvz(j) = tz
    end do
  end subroutine ipart_rhs


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Brownian kick: pv += bamp*g, skipping wall-trapped particles
  !! in stick mode
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ipart_brown(n,bamp,g,pvx,pvy,pvz,iwall,wbot,wtop,pz)
    implicit none
    integer      , intent(in)    :: n,iwall
    real(kind=GP), intent(in)    :: bamp,wbot,wtop
    real(kind=GP), intent(in)    :: g(3,n),pz(n)
    real(kind=GP), intent(inout) :: pvx(n),pvy(n),pvz(n)
    integer                      :: j
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active)
#else
!$omp parallel do
#endif
    do j = 1,n
      if ( (iwall.eq.2) .and. ((pz(j).le.wbot).or.(pz(j).ge.wtop)) ) cycle
      pvx(j) = pvx(j) + bamp*g(1,j)
      pvy(j) = pvy(j) + bamp*g(2,j)
      pvz(j) = pvz(j) + bamp*g(3,j)
    end do
  end subroutine ipart_brown


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Stick / trap wall mode: pin to the wall with zero velocity,
  !! in the current state and in the start-of-step state (so the
  !! remaining stages keep the particle frozen)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ipart_stick(n,wbot,wtop,pxo,pyo,pzo,pvxo,pvyo,pvzo, &
                         pxi,pyi,pzi,pvxi,pvyi,pvzi)
    implicit none
    integer      , intent(in)    :: n
    real(kind=GP), intent(in)    :: wbot,wtop
    real(kind=GP), intent(inout) :: pxo(n),pyo(n),pzo(n),pvxo(n),pvyo(n),pvzo(n)
    real(kind=GP), intent(inout) :: pxi(n),pyi(n),pzi(n),pvxi(n),pvyi(n),pvzi(n)
    integer                      :: j
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active)
#else
!$omp parallel do
#endif
    do j = 1,n
      if ( pzo(j) .le. wbot ) then
        pzo(j) = wbot
      else if ( pzo(j) .ge. wtop ) then
        pzo(j) = wtop
      else
        cycle
      endif
      pvxo(j) = 0.0_GP; pvyo(j) = 0.0_GP; pvzo(j) = 0.0_GP
      pxi (j) = pxo(j); pyi (j) = pyo(j); pzi (j) = pzo(j)
      pvxi(j) = 0.0_GP; pvyi(j) = 0.0_GP; pvzi(j) = 0.0_GP
    end do
  end subroutine ipart_stick


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Elastic reflection at the walls. At an intermediate stage
  !! (interm=1) the normal velocity reverses and the start-of-step
  !! state (pzi,pvzi) is mirrored about the wall (image method); at
  !! the final stage the velocity is pushed away from the wall
  !! (pvz = +|pvz| bottom, -|pvz| top). flag records the reflected
  !! particles (0 none, 1 bottom, 2 top), ncb/nct count them.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ipart_reflect(n,wbot,wtop,interm,pzo,pvzo,pzi,pvzi,flag,ncb,nct)
    implicit none
    integer      , intent(in)    :: n,interm
    real(kind=GP), intent(in)    :: wbot,wtop
    real(kind=GP), intent(inout) :: pzo(n),pvzo(n),pzi(n),pvzi(n)
    integer      , intent(inout) :: flag(n)
    integer      , intent(out)   :: ncb,nct
    integer                      :: j
    ncb = 0
    nct = 0
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active) reduction(+:ncb,nct)
#else
!$omp parallel do reduction(+:ncb,nct)
#endif
    do j = 1,n
      flag(j) = 0
      if ( pzo(j) .lt. wbot ) then
        pzo(j) = 2.0_GP*wbot - pzo(j)
        if ( interm .eq. 1 ) then
          pvzo(j) = -pvzo(j)
          pzi (j) = 2.0_GP*wbot - pzi(j)
          pvzi(j) = -pvzi(j)
        else
          pvzo(j) = ABS(pvzo(j))
        endif
        ncb = ncb + 1
        flag(j) = 1
      else if ( pzo(j) .ge. wtop ) then
        pzo(j) = 2.0_GP*wtop - pzo(j)
        if ( interm .eq. 1 ) then
          pvzo(j) = -pvzo(j)
          pzi (j) = 2.0_GP*wtop - pzi(j)
          pvzi(j) = -pvzi(j)
        else
          pvzo(j) = -ABS(pvzo(j))
        endif
        nct = nct + 1
        flag(j) = 2
      endif
    end do
  end subroutine ipart_reflect


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Deposition flags: flag=1 for the particles that stay (strictly
  !! inside the walls), 0 for the deposited ones; nb/nt count the
  !! removals per wall
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ipart_depflag(n,wbot,wtop,pz,flag,nb,nt)
    implicit none
    integer      , intent(in)    :: n
    real(kind=GP), intent(in)    :: wbot,wtop
    real(kind=GP), intent(in)    :: pz(n)
    integer      , intent(inout) :: flag(n)
    integer      , intent(out)   :: nb,nt
    integer                      :: j
    nb = 0
    nt = 0
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active) reduction(+:nb,nt)
#else
!$omp parallel do reduction(+:nb,nt)
#endif
    do j = 1,n
      if ( pz(j) .le. wbot ) then
        flag(j) = 0
        nb = nb + 1
      else if ( pz(j) .ge. wtop ) then
        flag(j) = 0
        nt = nt + 1
      else
        flag(j) = 1
      endif
    end do
  end subroutine ipart_depflag


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! flag = 1 - flag
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ipart_invflag(n,flag)
    implicit none
    integer, intent(in)    :: n
    integer, intent(inout) :: flag(n)
    integer                :: j
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active)
#else
!$omp parallel do
#endif
    do j = 1,n
      flag(j) = 1 - flag(j)
    end do
  end subroutine ipart_invflag


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! dst(k) = src(idx(k)) for k = 1..m (integer gather)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ipart_gather_i(m,idx,src,dst)
    implicit none
    integer, intent(in)    :: m
    integer, intent(in)    :: idx(*),src(*)
    integer, intent(inout) :: dst(*)
    integer                :: k
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active)
#else
!$omp parallel do
#endif
    do k = 1,m
      dst(k) = src(idx(k))
    end do
  end subroutine ipart_gather_i


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! In-place compaction of a real array through the ascending index
  !! list idx(1:m): a(k) = a(idx(k)), via the scratch array scr (the
  !! gather and the copy-back are separate kernels, so there are no
  !! read/write races between threads)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ipart_compact_r(m,idx,a,scr)
    implicit none
    integer      , intent(in)    :: m
    integer      , intent(in)    :: idx(*)
    real(kind=GP), intent(inout) :: a(*),scr(*)
    integer                      :: k
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active)
#else
!$omp parallel do
#endif
    do k = 1,m
      scr(k) = a(idx(k))
    end do
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active)
#else
!$omp parallel do
#endif
    do k = 1,m
      a(k) = scr(k)
    end do
  end subroutine ipart_compact_r


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! In-place compaction of an integer array through idx(1:m)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ipart_compact_i(m,idx,a,scr)
    implicit none
    integer, intent(in)    :: m
    integer, intent(in)    :: idx(*)
    integer, intent(inout) :: a(*),scr(*)
    integer                :: k
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active)
#else
!$omp parallel do
#endif
    do k = 1,m
      scr(k) = a(idx(k))
    end do
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active)
#else
!$omp parallel do
#endif
    do k = 1,m
      a(k) = scr(k)
    end do
  end subroutine ipart_compact_i


  ! ===================================================================
  ! Output methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute and write particle states.
  !! Writes: positions (xlg), fluid velocity (vlg),
  !!         inertial particle velocity (vip)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_pstate_impl(this, time, pde, fluidstate, pstate)
    use equationbase_mod
    use particlebase_mod
    use status
    use pstatus
    use fft
    class       (IPart),             intent(inout) :: this
    class(EquationBase),             intent   (in) :: pde
    type   (GStateComp), target ,    intent   (in) :: fluidstate(:)
    type  (GPStateComp), target ,    intent   (in) :: pstate(:)
    real      (kind=GP),             intent   (in) :: time
    complex   (kind=GP), pointer, dimension(:,:,:) :: velc, velc2
    complex   (kind=GP), pointer, dimension(:,:,:) :: vcx,vcy,vcz
    real      (kind=GP), pointer, dimension(:,:,:) :: velr,tmp1,tmp2
    real      (kind=GP)                            :: rmp
    logical                                        :: bret,wasdev

    call this%workspace_%get_complex_tmp(velc,bret)
    call this%workspace_%get_real_tmp   (velr,bret)
    call this%workspace_%get_real_tmp   (tmp1,bret)
    call this%workspace_%get_real_tmp   (tmp2,bret)
    call this%AssignLagPos(pstate)
    ! The interpolations run on the device copies (the fluid state is
    ! current there); the I/O routines switch to the host copies
    wasdev = gdev_active
    gdev_active = .TRUE.

    select type (pde)
    class is (VelocityBase)
      rmp = 1.0_GP/(real(this%nd_(1),kind=GP)*real(this%nd_(2),kind=GP)* &
                    real(this%nd_(3),kind=GP))
      vcx => fluidstate(pde%VELOCITY  )%ccomp
      vcy => fluidstate(pde%VELOCITY+1)%ccomp
      vcz => fluidstate(pde%VELOCITY+2)%ccomp
      ! Interpolate fluid velocity to particle positions (smoothed with
      ! fc_filter as in dpdt)
      CALL ipart_copy_scal(vcx,velc,rmp)
      CALL fc_filter(velc)
      CALL fftp3d_complex_to_real(planfc,velc,velr,MPI_COMM_WORLD)
      CALL this%EulerToLag(this%lvx_,this%nparts_,velr,.true. ,tmp1,tmp2)
      CALL ipart_copy_scal(vcy,velc,rmp)
      CALL fc_filter(velc)
      CALL fftp3d_complex_to_real(planfc,velc,velr,MPI_COMM_WORLD)
      CALL this%EulerToLag(this%lvy_,this%nparts_,velr,.false.,tmp1,tmp2)
      CALL ipart_copy_scal(vcz,velc,rmp)
      CALL fc_filter(velc)
      CALL fftp3d_complex_to_real(planfc,velc,velr,MPI_COMM_WORLD)
      CALL this%EulerToLag(this%lvz_,this%nparts_,velr,.false.,tmp1,tmp2)
      ! Write positions and Lagrangian fluid velocity
      WRITE(lgext,lgfmtext) pind
      CALL this%io_write_pdb(1,this%odir_,trim(this%sstate_pos_),lgext,time)
      CALL this%io_write_vec(1,this%odir_,trim(this%sstate_lag_),lgext,time)
      ! Write inertial particle velocity
      call gcopy(this%lvx_, pstate(this%VELOCITY  )%rcomp)
      call gcopy(this%lvy_, pstate(this%VELOCITY+1)%rcomp)
      call gcopy(this%lvz_, pstate(this%VELOCITY+2)%rcomp)
      CALL this%io_write_vec(1,this%odir_,trim(this%sstate_vel_),lgext,time)
! partlod >= 2: write Lagrangian vorticity and strain-rate tensor
      if ( this%traits_%partlod .ge. 2 ) then
! Write Lagrangian vorticity components
        CALL curlk(vcy,vcz,velc,1)
        CALL ipart_scal(velc,rmp)
        CALL fc_filter(velc)
        CALL fftp3d_complex_to_real(planfc,velc,velr,MPI_COMM_WORLD)
        CALL this%EulerToLag(this%lvx_,this%nparts_,velr,.false.,tmp1,tmp2)
        CALL curlk(vcx,vcz,velc,2)
        CALL ipart_scal(velc,rmp)
        CALL fc_filter(velc)
        CALL fftp3d_complex_to_real(planfc,velc,velr,MPI_COMM_WORLD)
        CALL this%EulerToLag(this%lvy_,this%nparts_,velr,.false.,tmp1,tmp2)
        CALL curlk(vcx,vcy,velc,3)
        CALL ipart_scal(velc,rmp)
        CALL fc_filter(velc)
        CALL fftp3d_complex_to_real(planfc,velc,velr,MPI_COMM_WORLD)
        CALL this%EulerToLag(this%lvz_,this%nparts_,velr,.false.,tmp1,tmp2)
        CALL this%io_write_vec(1,this%odir_,'wlg',lgext,time)
! Write strain-rate tensor components
        call this%workspace_%get_complex_tmp(velc2,bret)
        ! S11 = dv_x/dx
        CALL derivk(vcx, velc, 1)
        CALL ipart_scal(velc,rmp)
        CALL fc_filter(velc)
        CALL fftp3d_complex_to_real(planfc,velc,velr,MPI_COMM_WORLD)
        CALL this%io_write_euler(1,this%odir_,'s11',lgext,time,velr,.false.,tmp1,tmp2)
        ! S12 = 0.5*(dv_x/dy + dv_y/dx)
        CALL derivk(vcx, velc,  2)
        CALL derivk(vcy, velc2, 1)
        CALL ipart_combine(velc,velc2,0.5_GP*rmp)
        CALL fc_filter(velc)
        CALL fftp3d_complex_to_real(planfc,velc,velr,MPI_COMM_WORLD)
        CALL this%io_write_euler(1,this%odir_,'s12',lgext,time,velr,.false.,tmp1,tmp2)
        ! S13 = 0.5*(dv_x/dz + dv_z/dx)
        CALL derivk(vcx, velc,  3)
        CALL derivk(vcz, velc2, 1)
        CALL ipart_combine(velc,velc2,0.5_GP*rmp)
        CALL fc_filter(velc)
        CALL fftp3d_complex_to_real(planfc,velc,velr,MPI_COMM_WORLD)
        CALL this%io_write_euler(1,this%odir_,'s13',lgext,time,velr,.false.,tmp1,tmp2)
        ! S22 = dv_y/dy
        CALL derivk(vcy, velc, 2)
        CALL ipart_scal(velc,rmp)
        CALL fc_filter(velc)
        CALL fftp3d_complex_to_real(planfc,velc,velr,MPI_COMM_WORLD)
        CALL this%io_write_euler(1,this%odir_,'s22',lgext,time,velr,.false.,tmp1,tmp2)
        ! S23 = 0.5*(dv_y/dz + dv_z/dy)
        CALL derivk(vcy, velc,  3)
        CALL derivk(vcz, velc2, 2)
        CALL ipart_combine(velc,velc2,0.5_GP*rmp)
        CALL fc_filter(velc)
        CALL fftp3d_complex_to_real(planfc,velc,velr,MPI_COMM_WORLD)
        CALL this%io_write_euler(1,this%odir_,'s23',lgext,time,velr,.false.,tmp1,tmp2)
        call this%workspace_%free_complex_tmp(velc2)
      endif
    class default
      stop "inerpart: This solver does not support pdes without a velocity field"
    end select

    call this%workspace_%free_complex_tmp(velc)
    call this%workspace_%free_real_tmp   (velr)
    call this%workspace_%free_real_tmp   (tmp1)
    call this%workspace_%free_real_tmp   (tmp2)
    gdev_active = wasdev
  end subroutine write_pstate_impl


  ! ===================================================================
  ! Solver specific methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Constructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Ipart_ctor(this,infile, pde, workspace, pstate, pstate_cpy)
    USE equationbase_mod
    USE hd_mod,    ONLY: HDSolver
    USE var
    USE grid
    USE boxsize
    USE mpivars
    USE pstatus
    USE commtypes
    USE status
    USE random
    IMPLICIT NONE
    CLASS       (IPart), intent(inout)              :: this
    class(EquationBase), intent   (in)              :: pde
    type   (GWorkspace), intent(inout), target      :: workspace
    type  (GPStateComp), intent(inout), allocatable :: pstate(:), pstate_cpy(:)
    character   (len=*), intent   (in)              :: infile
    integer                                         :: tsta,tend,num_components
    integer                                         :: j,nzp,szreal
    real      (kind=GP)                             :: dpart
    logical                                         :: bret

    this%infile_      =  infile
    this%sclass_      =  'inerpart'
    this%idir_        =  pde%idir_ ! input  directory, same as in the pde class
    this%odir_        =  pde%odir_ ! output directory, same as in the pde class
    this%workspace_   => workspace
    this%hasfeedback_ = .false.
    this%nc_          = 3          ! fixed for now
    this%nparts_      = 0
    this%npartsm_     = 0
    this%nvdb_        = 0
    this%comm_        = MPI_COMM_WORLD
    this%maxparts_    = maxparts
    ! SPECTER geometry: x and y are periodic (spacing 2.pi.L/N); z is
    ! bounded, with Cz FC-Gram continuation points outside the physical
    ! domain [0,Lz]. The solver grid places nzp = nz-Cz physical points
    ! at z(k) = (k-1) Lz/(nzp-1), so the particle grid spacing in z is
    ! dz = Lz/(nzp-1), with the top wall at physical node nzp (grid
    ! coord nzp-1). The spline module must see the FULL extended grid
    ! (nd_(3) = nz): FC-Gram makes it periodic by construction, which
    ! the cyclic tridiagonal z-solve requires.
    IF ( (Cx.NE.0) .OR. (Cy.NE.0) ) THEN
      STOP 'Inerpart_ctor: only non-periodic z is supported (Cx = Cy = 0)'
    ENDIF
    IF ( Cz.EQ.0 ) THEN
      STOP 'Inerpart_ctor: SPECTER inertial particles need a non-periodic z (Cz > 0)'
    ENDIF
    nzp               = nz - Cz
    this%nd_(1)       = nx
    this%nd_(2)       = ny
    this%nd_(3)       = nz
    this%delta_(1)    = 2*pi*Lx/real(nx,kind=GP)
    this%delta_(2)    = 2*pi*Ly/real(ny,kind=GP)
    this%delta_(3)    = Lz/real(nzp-1,kind=GP)
    this%invdel_(1)   = real(nx,kind=GP)/(2*pi*Lx)
    this%invdel_(2)   = real(ny,kind=GP)/(2*pi*Ly)
    this%invdel_(3)   = real(nzp-1,kind=GP)/Lz
    this%seedfile_    = lgseedfile
    this%iinterp_     = 3          ! fixed for now
    this%itorder_     = 2
    this%intorder_    = max(intorder,1)
    this%iseed_       = 1000
    this%istep_       = 0
    this%iexchtype_   = ilgexchtype
    this%iouttype_    = ilgouttype
    this%bcollective_ = ilgcoll
    this%itimetype_   = GT_WTIME
    this%wrtunit_     = ilgwrtunit
    CALL this%SetRandSeed(seed)
    CALL MPI_COMM_SIZE(this%comm_,this%nprocs_,this%ierr_)
    CALL MPI_COMM_RANK(this%comm_,this%myrank_,this%ierr_)

    IF (this%iexchtype_.EQ.GPEXCHTYPE_VDB) THEN
      this%partbuff_ = maxparts
    ELSE IF (this%iexchtype_.EQ.GPEXCHTYPE_NN) THEN
      this%partbuff_      = 1 + (maxparts - 1)/this%nprocs_
      this%partchunksize_ = (this%partbuff_ + 9)/10
      this%partbuff_      =  this%partbuff_ + this%partchunksize_
      IF ((this%bcollective_.EQ.0).AND.(this%myrank_.EQ.0)) THEN
        this%partbuff_   = maxparts
      END IF
      this%stepcounter_ = 0
    END IF

    ! Initialize timers (get handles):
    DO j = 1, GPMAXTIMERS
      CALL GTInitHandle(this%htimers_(j),this%itimetype_)
      IF ( this%htimers_(j).EQ.GTNULLHANDLE ) THEN
        WRITE(*,*) 'Inerpart_ctor: Not enough timers available'
        STOP
      ENDIF
    ENDDO

    ! Initialize communicators
    CALL this%gpcomm_%GPartComm_ctor(GPCOMM_INTRFC_SF,this%partbuff_, &
         this%nd_,GPSI_NZGHOST,this%comm_,this%htimers_(GPTIME_COMM))
    CALL this%gpcomm_%Init()

    this%libnds_(1,1) = 1
    this%libnds_(1,2) = nx
    this%lxbnds_(1,1) = 0.0_GP
    this%lxbnds_(1,2) = real(nx,kind=GP)
    this%libnds_(2,1) = 1
    this%libnds_(2,2) = ny
    this%lxbnds_(2,1) = 0.0_GP
    this%lxbnds_(2,2) = real(ny,kind=GP)
    ! z local bounds: the full MPI slab (including the continuation
    ! zone if this rank holds part of it). The spline is valid over the
    ! whole extended domain; the particles are kept out of the
    ! continuation zone by the wall physics, not by the slab bounds.
    this%libnds_(3,1) = ksta
    this%libnds_(3,2) = kend
    this%lxbnds_(3,1) = real(ksta-1,kind=GP)
    this%lxbnds_(3,2) = real(kend-1,kind=GP) + 1.0_GP
    CALL range(1,nx,nprocs,myrank,tsta,tend)
    this%tibnds_(1,1) = 1
    this%tibnds_(1,2) = nz
    this%tibnds_(2,1) = 1
    this%tibnds_(2,2) = ny
    this%tibnds_(3,1) = tsta
    this%tibnds_(3,2) = tend
    ! Global extent in grid units: in z this is the TOP-WALL grid
    ! COORDINATE nzp-1 (node nzp sits at z=Lz), consistent with
    ! delta_(3) = Lz/(nzp-1); using nzp here stretched the particle
    ! z-coordinates by nzp/(nzp-1) and pushed the near-top-wall
    ! stencils into the continuation zone (see SPECTER_v1's
    ! KNOWN_LIMITATION_topwall_interp.md)
    this%gext_(1) = real(nx,kind=GP)
    this%gext_(2) = real(ny,kind=GP)
    this%gext_(3) = real(nzp-1,kind=GP)

    ! Call init (reads &inerpart namelist, sets POSITION/VELOCITY indices)
    call this%init()

    ! Coefficients that need the viscosity (from the PDE solver): the
    ! nonlinear drag Re_p constant, and the particle diameter
    ! d = sqrt(18 tau nu gamma) (from tau = rho_p d^2/(18 rho_f nu)),
    ! which sets the Saffman lift prefactor and the wall contact
    ! half-gap d/2 in grid units. These require a physical mass ratio
    ! (gamma > 0); otherwise they are disabled (point particles,
    ! linear Stokes drag).
    if ( this%traits_%donldrag .or. this%traits_%dolift .or. &
         (this%traits_%gamma .gt. 0.0_GP) ) then
      select type (pde)
        type is (HDSolver)
          this%traits_%nu = pde%traits_%nu
        class default
          stop "Inerpart_ctor: cannot get nu from this PDE solver"
      end select
    endif
    if ( this%traits_%gamma .gt. 0.0_GP ) then
      if ( this%traits_%nu .le. 0.0_GP ) then
        stop "Inerpart_ctor: gamma > 0 requires nu > 0 in the PDE solver"
      endif
      dpart = SQRT(18.0_GP*this%traits_%tau*this%traits_%nu*this%traits_%gamma)
      this%traits_%clift   = 1.615_GP*(6.0_GP/pi)*this%traits_%gamma* &
                             SQRT(this%traits_%nu)/dpart
      this%traits_%rpwall  = 0.5_GP*dpart*real(nzp-1,kind=GP)/Lz
      this%traits_%nld_rep = 18.0_GP*this%traits_%tau*this%traits_%gamma/ &
                             this%traits_%nu
    else
      this%traits_%clift   = 0.0_GP
      this%traits_%rpwall  = 0.0_GP   ! point-particle (center-on-plane)
      this%traits_%nld_rep = 0.0_GP
      if ( (this%traits_%dolift .or. this%traits_%donldrag) .and. &
           (this%myrank_.eq.0) ) then
        write(*,*) 'Inerpart_ctor: WARNING gamma<=0, Saffman lift and ', &
                   'nonlinear drag disabled (particle diameter undefined).'
      endif
    endif
    if ( this%traits_%donldrag .and. (this%traits_%nld_rep.le.0.0_GP) ) then
      this%traits_%donldrag = .false.
    endif
    if ( this%traits_%doreinj .and. (this%traits_%iwall.ne.0) .and. &
         (this%myrank_.eq.0) ) then
      write(*,*) 'Inerpart_ctor: WARNING doreinj has no effect unless ', &
                 'iwall=0 (particles are only removed by deposition).'
    endif
    if ( this%traits_%doreinj .and. (this%iexchtype_.ne.GPEXCHTYPE_NN) .and. &
         (this%myrank_.eq.0) ) then
      write(*,*) 'Inerpart_ctor: WARNING doreinj requires NN exchange; ', &
                 'no deposition occurs in VDB mode so nothing is re-injected.'
    endif

    ! Instantiate interp operation and select the near-wall
    ! interpolation scheme
    CALL this%intop_%GPSplineInt_ctor(3,this%nd_,this%libnds_,this%lxbnds_, &
         this%tibnds_,this%partbuff_,this%gpcomm_,&
         this%htimers_(GPTIME_DATAEX),this%htimers_(GPTIME_TRANSP))
    CALL this%intop_%SetWallInterp(this%traits_%iwallinterp,nzp)

    ! Allocate particle arrays
    CALL MPI_TYPE_SIZE(GC_REAL,szreal,this%ierr_)
    CALL galloc(this%id_    ,this%partbuff_)
    CALL galloc(this%tmpint_,this%partbuff_)
    CALL galloc(this%iwrk_  ,this%partbuff_)
    CALL galloc(this%ptmp0_ ,3,this%partbuff_)
    CALL galloc(this%gptmp0_,3,this%partbuff_)
    ALLOCATE(this%idrem_(this%partbuff_))
    this%idrem_ = 0
    num_components = this%state_size()  ! Returns 6
    CALL GPState_alloc(pstate    , num_components, this%partbuff_)
    CALL GPState_alloc(pstate_cpy, num_components, this%partbuff_)
    call this%workspace_%set_nparts(this%partbuff_)
    call this%workspace_%init_pcomp_arrays(this%partbuff_)
    call this%workspace_%get_pcomp_tmp(this%lvx_,bret)
    call this%workspace_%get_pcomp_tmp(this%lvy_,bret)
    call this%workspace_%get_pcomp_tmp(this%lvz_,bret)
    call this%workspace_%get_pcomp_tmp(this%lwx_,bret)
    call this%workspace_%get_pcomp_tmp(this%lwy_,bret)
    call this%workspace_%get_pcomp_tmp(this%lwz_,bret)
    IF ( this%iexchtype_.EQ.GPEXCHTYPE_VDB ) THEN
      CALL galloc(this%vdb_,3,this%partbuff_)
    ENDIF
  END SUBROUTINE Ipart_ctor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Destructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Ipart_dtor(this)
    IMPLICIT NONE
    TYPE(Ipart),INTENT(INOUT) :: this
    integer                   :: j
    CALL gfree(this%id_)
    CALL gfree(this%tmpint_)
    CALL gfree(this%iwrk_)
    CALL gfree(this%vdb_)
    CALL gfree(this%ptmp0_)
    CALL gfree(this%gptmp0_)
    IF ( ALLOCATED (this%idrem_) ) DEALLOCATE(this%idrem_)
    IF ( ASSOCIATED   (this%px_) ) NULLIFY       (this%px_)
    IF ( ASSOCIATED   (this%py_) ) NULLIFY       (this%py_)
    IF ( ASSOCIATED   (this%pz_) ) NULLIFY       (this%pz_)
    IF ( ASSOCIATED(this%lvx_) ) CALL this%workspace_%free_pcomp_tmp(this%lvx_)
    IF ( ASSOCIATED(this%lvy_) ) CALL this%workspace_%free_pcomp_tmp(this%lvy_)
    IF ( ASSOCIATED(this%lvz_) ) CALL this%workspace_%free_pcomp_tmp(this%lvz_)
    IF ( ASSOCIATED(this%lwx_) ) CALL this%workspace_%free_pcomp_tmp(this%lwx_)
    IF ( ASSOCIATED(this%lwy_) ) CALL this%workspace_%free_pcomp_tmp(this%lwy_)
    IF ( ASSOCIATED(this%lwz_) ) CALL this%workspace_%free_pcomp_tmp(this%lwz_)
    IF ( this%logunit_ .NE.0 ) CLOSE(this%logunit_)
    IF ( this%diagunit_.NE.0 ) CLOSE(this%diagunit_)
    DO j = 1, GPMAXTIMERS
      CALL GTFree(this%htimers_(j))
    ENDDO
  END SUBROUTINE Ipart_dtor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : Ipart_ResizeArrays
  !!  DESCRIPTION: As the base class, but also releases and
  !!               re-acquires the Lagrangian vorticity pointers
  !!               lwx_, lwy_, lwz_: the base ResizeArrays
  !!               reallocates every array of the particle pool
  !!               (and re-acquires its own lvx_, lvy_, lvz_), so
  !!               any other held pointer into the pool dangles
  !!               after the call.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Ipart_ResizeArrays(this,new_size,onlyinc,exc)
    IMPLICIT NONE
    CLASS(Ipart) ,INTENT(INOUT)          :: this
    INTEGER      ,INTENT(IN)             :: new_size
    LOGICAL      ,INTENT(IN)             :: onlyinc
    LOGICAL      ,INTENT(IN)   ,OPTIONAL :: exc
    LOGICAL                              :: bret,assocw(3)

    assocw = .false.
    IF ( ASSOCIATED(this%lwx_) ) THEN
      assocw(1) = .true.
      call this%workspace_%free_pcomp_tmp(this%lwx_)
    END IF
    IF ( ASSOCIATED(this%lwy_) ) THEN
      assocw(2) = .true.
      call this%workspace_%free_pcomp_tmp(this%lwy_)
    END IF
    IF ( ASSOCIATED(this%lwz_) ) THEN
      assocw(3) = .true.
      call this%workspace_%free_pcomp_tmp(this%lwz_)
    END IF
    CALL ResizeArrays(this,new_size,onlyinc,exc)
    IF (assocw(1)) call this%workspace_%get_pcomp_tmp(this%lwx_,bret)
    IF (assocw(2)) call this%workspace_%get_pcomp_tmp(this%lwy_,bret)
    IF (assocw(3)) call this%workspace_%get_pcomp_tmp(this%lwz_,bret)
  END SUBROUTINE Ipart_ResizeArrays


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute number of state members (equations)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE function state_size_impl(this) result(num)
    class(Ipart), intent(in) :: this
    integer                  :: num
    num = 2 * this%nc_          ! 3 positions + 3 velocities
  end function state_size_impl

end module inerpart_mod
