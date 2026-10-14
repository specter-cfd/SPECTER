! =====================================================================
! NAME       : hd_solver.f90
! DESCRIPTION: Forms class for incompressible HD solver, computing:
!
!              dv/dt + (v.Grad) v = - Grad p + nu Del^2 v + f
!
!              State ordering is:
!                vx, vy, vz
!
!              State sector ids are:
!                VELOCITY (VELOCITY+1, VELOCITY+2)
!
! Versión de SPECTER del hd_solver.f90 de GHOST (rama new_paradigm,
! afaf9ea), que se arma de a un método por vez. init lee &HD y
! &velbound, arma el plan de bordes, y reserva e inicializa la presión
! (init_pressure: cero en una corrida nueva, leída del archivo pr al
! continuar una); dudt calcula el lado derecho sin la presión y
! end_stage impone las paredes y proyecta después de cada etapa del RK
! (las llama el integrador, src/steppers/canuto_stepper.f90); global y
! write_states escriben las salidas. De &HD faltan la rotación, los
! escalares pasivos y spectlod de GHOST, que SPECTER no tiene.
!
! La presión pr_ (p' = p.dt, en el dominio (z,ky,kx)) es un arreglo
! auxiliar del solver, como los de los solvers cuánticos de GHOST: no
! es parte del estado que avanza el integrador, sobrevive entre etapas
! y pasos de tiempo, y la usan la proyección y la escritura.
!
! INPUT FILE : looks for a "&HD" namelist with:
!              fidir   : changes class binary input  dir (default: idir)
!              fodir   : changes class binary output dir (default: odir)
!              todir   : changes the class TXT output dir (default: '.')
!              nu      : fluid kinematic viscosity
!
!              and a "&velbound" namelist with:
!              vbcista : boundary condition in the i direction at the
!                        start of the domain ('periodic' or 'noslip')
!              vbciend : boundary condition in the i direction at the
!                        end of the domain
!              vizsta  : boundary velocity at z=0 in the i direction
!              vizend  : boundary velocity at z=Lz in the i direction
!
! DATE       : 11/30/25 (DLR)
! =====================================================================

module hd_mod
  use equationbase_mod
  use gstate_mod
  use boundary

  implicit none

  ! ================= Solver traits ===================================
  type, public  :: NHTraits
    real(kind=GP) :: nu = 0.0_GP  ! kinematic viscosity
  end type NHTraits

  ! ================= Solver ==========================================
  ! Define class:
  type, extends(VelocityBase) :: HDSolver
    ! Member data:
    logical           :: binit_  = .false. ! is initialized?
    type (NHTraits)   :: traits_
    type  (BCPLAN)    :: vplanbc_          ! plan de bordes de la velocidad
    real (kind=GP)    :: vzsta_(2)         ! (vx,vy) de la pared en z=0
    real (kind=GP)    :: vzend_(2)         ! (vx,vy) de la pared en z=Lz
    complex(kind=GP), allocatable, dimension(:,:,:) :: pr_ ! presión p'
  CONTAINS
    procedure, public :: init          =>          init_impl ! init method
    procedure, public :: dudt          =>          dudt_impl ! RHS method
    procedure, public :: end_stage     =>     end_stage_impl ! RK stage correction
    procedure, public :: global        =>        global_impl ! Writes global qtys
    procedure, public :: write_states  =>  write_states_impl ! Binary output
    procedure, public :: state_size    =>    state_size_impl ! state size
    procedure, public :: get_sstate    =>    get_sstate_impl ! get state name list
    procedure, public :: Solver_ctor   =>      HDSolver_ctor ! constructor
    final             :: HDSolver_dtor
  end type HDSolver

CONTAINS

  ! ===================================================================
  ! Solver initialization, this is where parameter files are read
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to initialize the solver
  !! Reads the &HD and &velbound namelists, sets up the
  !! boundary conditions, allocates and initializes the
  !! pressure and sets sector indices. Uses the FC-Gram plan,
  !! so it must be called after fcgram_create_plan.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_impl(this)
    use commtypes
    use status
    use fft
    class  (HDSolver), intent (inout) :: this

    ! Temporary data to read from namelists:
    character(len=15)  :: vbcxsta,vbcxend,vbcysta,vbcyend,vbczsta,vbczend
    real  (kind=GP)    :: vxzsta,vxzend,vyzsta,vyzend
    real  (kind=GP)    :: nu
    character(len=128) :: fidir,fodir,todir
    integer            :: ierr

    ! Required namelists:
    namelist / HD       / fidir,fodir,todir,nu
    namelist / velbound / vbcxsta,vbcxend,vbcysta,vbcyend,vbczsta,vbczend
    namelist / velbound / vxzsta,vxzend,vyzsta,vyzend

    call MPI_COMM_SIZE(MPI_COMM_WORLD,this%nprocs_,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,this%myrank_,ierr)

    ! Get I/O and trait variables from input file:
    fidir = idir ! Set the default to status idir
    fodir = odir ! Set the default to status odir
    todir = '.'  ! Set the default to the current dir
    nu    = 0.0_GP
    if ( this%myrank_ .eq. 0 ) then
      open(1,file=this%infile_,status='unknown',form="formatted")
      read(1,NML=HD)
      close(1)
    endif
    call MPI_BCAST(fidir,128,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(fodir,128,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(todir,128,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(nu   ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    this%traits_%nu = nu

    ! Boundary conditions for the velocity field
    vbcxsta = ""; vbcxend = ""
    vbcysta = ""; vbcyend = ""
    vbczsta = ""; vbczend = ""
    vxzsta = .0_GP; vxzend = .0_GP
    vyzsta = .0_GP; vyzend = .0_GP
    if ( this%myrank_ .eq. 0 ) then
      open(1,file=this%infile_,status='unknown',form="formatted")
      read(1,NML=velbound)
      close(1)
    endif
    call MPI_BCAST(vbcxsta,15,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(vbcxend,15,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(vbcysta,15,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(vbcyend,15,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(vbczsta,15,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(vbczend,15,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(vxzsta,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(vxzend,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(vyzsta,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(vyzend,1,GC_REAL,0,MPI_COMM_WORLD,ierr)

    ! Periodic directions get periodic boundary conditions
    if ( Cx .eq. 0 ) then
      vbcxsta = "periodic"
      vbcxend = "periodic"
    endif
    if ( Cy .eq. 0 ) then
      vbcysta = "periodic"
      vbcyend = "periodic"
    endif
    if ( Cz .eq. 0 ) then
      vbczsta = "periodic"
      vbczend = "periodic"
    endif
    call setup_bc(this%vplanbc_,planfc, &
         (/ vbcxsta,vbcxend,vbcysta,vbcyend,vbczsta,vbczend /),'v')
    this%vzsta_ = (/ vxzsta, vyzsta /)
    this%vzend_ = (/ vxzend, vyzend /)

    ! Pressure (auxiliary array of the solver)
    call galloc(this%pr_, nz, ny, ista, iend)

    ! Set I/O from inputfile data (if present in &HD, replaces the
    ! class defaults)
    this%idir_  = fidir
    this%odir_  = fodir
    this%todir_ = todir

    ! Initial pressure (uses idir_)
    call init_pressure(this)

    this%nd_      = 3                        ! 3d
    this%nc_      = this%nd_                 ! # field components
    this%VELOCITY = 1                        ! start of vel sector

    allocate(this%sstate_(this%state_size()))
    call this%get_sstate(this%sstate_)
    this%binit_ = .true.
  end subroutine init_impl


  ! ===================================================================
  ! Computation of RHS, the solver equations are defined here
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Computes the RHS of the momentum equation without the
  !! pressure, which is taken into account by end_stage:
  !!   dudt = nu Del^2 v - (v.Grad) v + f
  !! The nonlinear term is filtered (FC-Gram dealiasing). uin
  !! and dudt can be the same array (the stepper does that).
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dudt_impl(this, time, uin, uf, dt, dudt)
    use grid
    use mpivars
    implicit none

    class (HDSolver), intent(inout)             :: this
    real   (kind=GP), intent   (in)             :: time, dt
    type(GStateComp), intent(inout), target     :: uin(:),uf(:)
    type(GStateComp), intent(inout), target     :: dudt(:)
    complex(kind=GP), pointer, dimension(:,:,:) :: fx,fy,fz,vx,vy,vz
    complex(kind=GP), pointer, dimension(:,:,:) :: dx,dy,dz
    complex(kind=GP), pointer, dimension(:,:,:) :: C1,C2,C3,C4,C5,C6
    real   (kind=GP)                            :: nu
    integer                                     :: i,j,k
    logical                                     :: bret

    if ( .not. this%binit_ ) then
      stop 'HDSolver::dudt: Solver not initialized'
    endif

    nu = this%traits_%nu

    call this%workspace_%get_complex_tmp(C1,bret)
    call this%workspace_%get_complex_tmp(C2,bret)
    call this%workspace_%get_complex_tmp(C3,bret)
    call this%workspace_%get_complex_tmp(C4,bret)
    call this%workspace_%get_complex_tmp(C5,bret)
    call this%workspace_%get_complex_tmp(C6,bret)

    vx => uin (this%VELOCITY  )%ccomp
    vy => uin (this%VELOCITY+1)%ccomp
    vz => uin (this%VELOCITY+2)%ccomp
    fx => uf  (this%VELOCITY  )%ccomp
    fy => uf  (this%VELOCITY+1)%ccomp
    fz => uf  (this%VELOCITY+2)%ccomp
    dx => dudt(this%VELOCITY  )%ccomp
    dy => dudt(this%VELOCITY+1)%ccomp
    dz => dudt(this%VELOCITY+2)%ccomp

    call gradre(vx,vy,vz,C4,C5,C6)   ! (v.Grad) v
    call fc_filter(C4)               ! dealiasing of the
    call fc_filter(C5)               ! nonlinear term
    call fc_filter(C6)
    call laplak(vx,C1)               ! Del^2 v
    call laplak(vy,C2)
    call laplak(vz,C3)

!$omp parallel do collapse(2) private (k)
    do i = ista,iend
       do j = 1,ny
          do k = 1,nz
             dx(k,j,i) = nu*C1(k,j,i)-C4(k,j,i)+fx(k,j,i)
             dy(k,j,i) = nu*C2(k,j,i)-C5(k,j,i)+fy(k,j,i)
             dz(k,j,i) = nu*C3(k,j,i)-C6(k,j,i)+fz(k,j,i)
          end do
       end do
    end do

    call this%workspace_%free_complex_tmp(C1)
    call this%workspace_%free_complex_tmp(C2)
    call this%workspace_%free_complex_tmp(C3)
    call this%workspace_%free_complex_tmp(C4)
    call this%workspace_%free_complex_tmp(C5)
    call this%workspace_%free_complex_tmp(C6)
  end subroutine dudt_impl


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Correction after each RK stage o (norder, ..., 1): imposes
  !! the velocity of the walls and projects v onto solenoidal
  !! fields. The pressure p' of the previous stage is used for
  !! the walls, and the new one is left in pr_.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine end_stage_impl(this, uout, o, norder)
    use fft
    implicit none

    class (HDSolver), intent(inout)             :: this
    type(GStateComp), intent(inout), target     :: uout(:)
    integer         , intent   (in)             :: o, norder

    call v_imposebc_and_project(this%vplanbc_,planfc,               &
             uout(this%VELOCITY)%ccomp,uout(this%VELOCITY+1)%ccomp, &
             uout(this%VELOCITY+2)%ccomp,this%pr_,rki=o,rko=norder, &
             v_zsta=this%vzsta_,v_zend=this%vzend_)
  end subroutine end_stage_impl


  ! ===================================================================
  ! Computation of global quantities
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute and write global quantities
  !!   balance.txt          : energy, enstrophy and injection rate
  !!   helicity.txt         : kinetic helicity
  !!   noslip_diagnostic.txt: divergence and boundary conditions
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine global_impl(this, uin, uf, t)
    use status
    use fft
    class (HDSolver), intent(in)                :: this
    type(GStateComp), intent(in), target        :: uin(:), uf(:)
    integer         , intent(in)                :: t
    complex(kind=GP), pointer, dimension(:,:,:) :: fx,fy,fz,vx,vy,vz

    vx => uin(this%VELOCITY  )%ccomp
    vy => uin(this%VELOCITY+1)%ccomp
    vz => uin(this%VELOCITY+2)%ccomp
    fx => uf (this%VELOCITY  )%ccomp
    fy => uf (this%VELOCITY+1)%ccomp
    fz => uf (this%VELOCITY+2)%ccomp
    call hdcheck(vx,vy,vz,fx,fy,fz,t,dt,1)
    call vdiagnostic(this%vplanbc_,planfc,vx,vy,vz,t,dt)
  end subroutine global_impl


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Writes the states and the vorticity (write_states of
  !! EquationBase) and the pressure, transformed from p' to p
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_states_impl(this, uin, planio)
    use grid
    use iovar
    use status
    use filefmt
    use fft
    use commtypes
    implicit none

    class   (HDSolver), intent   (in)             :: this
    type  (GStateComp), intent(inout)             :: uin(:)
    type      (ioplan), intent   (in)             :: planio
    complex  (kind=GP), pointer, dimension(:,:,:) :: C1
    real     (kind=GP), pointer, dimension(:,:,:) :: R1
    real     (kind=GP)                            :: rmp
    integer                                       :: i,j,k
    logical                                       :: bret

    call write_states(this, uin, planio)
    call this%workspace_%get_complex_tmp(C1,bret)
    call this%workspace_%get_real_tmp   (R1,bret)
    rmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*dt)
!$omp parallel do collapse(2) private (k)
    do i = ista,iend
       do j = 1,ny
          do concurrent (k=1:nz)
            C1(k,j,i) = this%pr_(k,j,i)*rmp
          end do
       end do
    end do
    call fftp2d_complex_to_real_xy(planfc,C1,R1,MPI_COMM_WORLD)
    call io_write(1,this%odir_,'pr',ext,planio,R1)
    call this%workspace_%free_complex_tmp(C1)
    call this%workspace_%free_real_tmp   (R1)
  end subroutine write_states_impl


  ! ===================================================================
  ! Solver specific methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Presión inicial p' = p.dt: cero en una corrida nueva
  !! (stat=0); al continuar una corrida se lee del archivo pr
  !! numerado por stat, como la velocidad (ic/ic_velocity.f90)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_pressure(this)
    use grid
    use iovar
    use status
    use filefmt
    use fft
    use commtypes
    implicit none

    class (HDSolver), intent(inout)             :: this
    real   (kind=GP), pointer, dimension(:,:,:) :: R1
    integer                                     :: i,j,k
    logical                                     :: bret

    this%pr_ = 0.0_GP
    if ( stat .eq. 0 ) return
    call this%workspace_%get_real_tmp(R1,bret)
    WRITE(ext, fmtext) stat
    call io_read(1,this%idir_,'pr',ext,this%planio_,R1)
    call fftp2d_real_to_complex_xy(planfc,R1,this%pr_,MPI_COMM_WORLD)
!$omp parallel do collapse(2) private (k)
    do i = ista,iend
       do j = 1,ny
          do concurrent (k=1:nz-Cz)
             this%pr_(k,j,i) = this%pr_(k,j,i)*dt   ! p -> p'
          end do
       end do
    end do
    call this%workspace_%free_real_tmp(R1)
  end subroutine init_pressure


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Constructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine HDSolver_ctor(this, infile, workspace, plan)
    use iovar
    class  (HDSolver), intent(inout)         :: this
    type(GWorkspace) , intent(inout), target :: workspace
    type(ioplan)     , intent(inout), target :: plan
    character(len=*) , intent   (in)         :: infile
    this%infile_    =  infile    ! input file
    this%workspace_ => workspace
    this%planio_    => plan
    call this%init()
  end subroutine HDSolver_ctor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Destructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine HDSolver_dtor(this)
    use gmem
    type  (HDSolver), intent(inout) :: this
    if (associated(this%workspace_)) nullify(this%workspace_)
    if (associated(this%planio_))    nullify(this%planio_)
    if (allocated(this%sstate_))     deallocate(this%sstate_)
    if (allocated(this%pr_))         call gfree(this%pr_)
  end subroutine HDSolver_dtor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Get state variable names
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_sstate_impl(this, sstate)
    class  (HDSolver), intent   (in) :: this
    character (len=8), intent(inout) :: sstate(:)
    character(len=1)                 :: comp(3)
    integer                          :: j
    comp = ['x', 'y', 'z']
    do j = this%VELOCITY,this%VELOCITY+this%nc_-1
       sstate(j) = 'v' // comp(j-this%VELOCITY+1)
    enddo
  end subroutine get_sstate_impl


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute number of state members (equations)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE function state_size_impl(this) result(num)
    class(HDSolver), intent(in) :: this
    integer                     :: num
    num = this%nc_               ! # vel. components
  end function state_size_impl

end module hd_mod
