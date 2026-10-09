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
! afaf9ea), que se arma de a un método por vez. init lee nu de
! &velocity y &velbound, arma el plan de bordes y reserva la presión;
! dudt calcula el lado derecho sin la presión y end_stage impone las
! paredes y proyecta después de cada etapa del RK (las llama el
! integrador, src/steppers/canuto_stepper.f90); global y write_states
! escriben las salidas.
!
! La presión pr_ (p' = p.dt, en el dominio (z,ky,kx)) es un arreglo
! auxiliar del solver, como los de los solvers cuánticos de GHOST: no
! es parte del estado que avanza el integrador, sobrevive entre etapas
! y pasos de tiempo, y la usan la proyección y la escritura.
!
! INPUT FILE : looks for a "&velbound" namelist with:
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
    procedure, public :: spectra       =>       spectra_impl ! Writes spectra
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
  !! Reads the &velbound namelist, sets up the boundary
  !! conditions, allocates the pressure and sets sector
  !! indices. Uses the FC-Gram plan, so it must be called
  !! after fcgram_create_plan.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_impl(this)
    use commtypes
    use status
    use grid
    use fft
    use gmem
    use mpivars, only: ista, iend
    class  (HDSolver), intent (inout) :: this

    ! Temporary data to read from namelists:
    character(len=15) :: vbcxsta,vbcxend,vbcysta,vbcyend,vbczsta,vbczend
    real  (kind=GP)   :: vxzsta,vxzend,vyzsta,vyzend
    real  (kind=GP)   :: f0,u0,kdn,kup,nu
    real  (kind=GP)   :: fparam0,fparam1,fparam2,fparam3,fparam4
    real  (kind=GP)   :: fparam5,fparam6,fparam7,fparam8,fparam9
    real  (kind=GP)   :: vparam0,vparam1,vparam2,vparam3,vparam4
    real  (kind=GP)   :: vparam5,vparam6,vparam7,vparam8,vparam9
    integer           :: ierr

    ! Required namelists:
    namelist / velbound / vbcxsta,vbcxend,vbcysta,vbcyend,vbczsta,vbczend
    namelist / velbound / vxzsta,vxzend,vyzsta,vyzend
    namelist / velocity / f0,u0,kdn,kup,nu
    namelist / velocity / fparam0,fparam1,fparam2,fparam3,fparam4
    namelist / velocity / fparam5,fparam6,fparam7,fparam8,fparam9
    namelist / velocity / vparam0,vparam1,vparam2,vparam3,vparam4
    namelist / velocity / vparam5,vparam6,vparam7,vparam8,vparam9

    call MPI_COMM_SIZE(MPI_COMM_WORLD,this%nprocs_,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,this%myrank_,ierr)

    ! Viscosity, from the &velocity namelist (the rest of the list is
    ! for the initial condition and the forcing, still in the main)
    nu = 0.0_GP
    if ( this%myrank_ .eq. 0 ) then
      open(1,file=this%infile_,status='unknown',form="formatted")
      read(1,NML=velocity)
      close(1)
    endif
    call MPI_BCAST(nu,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
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

    ! I/O directories, as in GHOST without the &HD namelist
    this%idir_  = idir
    this%odir_  = odir
    this%todir_ = '.'

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

    if ( .not. this%binit_ ) then
      stop 'HDSolver::dudt: Solver not initialized'
    endif
    nu = this%traits_%nu

    call this%workspace_%get_complex_tmp(C1)
    call this%workspace_%get_complex_tmp(C2)
    call this%workspace_%get_complex_tmp(C3)
    call this%workspace_%get_complex_tmp(C4)
    call this%workspace_%get_complex_tmp(C5)
    call this%workspace_%get_complex_tmp(C6)

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
             uout(this%VELOCITY+2)%ccomp,this%pr_,rki=o,            &
             v_zsta=this%vzsta_,v_zend=this%vzend_)
  end subroutine end_stage_impl


  ! ===================================================================
  ! Computation of global quantities and spectra
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
  !! SPECTER no escribe espectros
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine spectra_impl(this, uin)
    use status
    class (HDSolver), intent(in)                :: this
    type(GStateComp), intent(in), target        :: uin(:)
    return
  end subroutine spectra_impl


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Writes the states and the vorticity (write_states_base)
  !! and the pressure, transformed from p' to p
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_states_impl(this, uin, planio)
    use grid
    use mpivars
    use commtypes
    use iovar
    use status
    use filefmt
    use fft
    implicit none

    class   (HDSolver), intent   (in)             :: this
    type  (GStateComp), intent(inout)             :: uin(:)
    type      (ioplan), intent   (in)             :: planio
    complex  (kind=GP), pointer, dimension(:,:,:) :: C1
    real     (kind=GP), pointer, dimension(:,:,:) :: R1
    real     (kind=GP)                            :: rmp
    integer                                       :: i,j,k

    call write_states_base(this, uin, planio)
    call this%workspace_%get_complex_tmp(C1)
    call this%workspace_%get_real_tmp   (R1)
    rmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*dt)
    do i = ista,iend
       do j = 1,ny
          do k = 1,nz
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
