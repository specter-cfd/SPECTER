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
! afaf9ea), que se arma de a un método por vez. Por ahora init lee
! &velbound y arma el plan de bordes de la velocidad; nu, el paso de
! RK, los globales y la escritura de binarios siguen en el main.
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

  ! ================= Solver ==========================================
  ! Define class:
  type, extends(VelocityBase) :: HDSolver
    ! Member data:
    logical           :: binit_  = .false. ! is initialized?
    type  (BCPLAN)    :: vplanbc_          ! plan de bordes de la velocidad
    real (kind=GP)    :: vzsta_(2)         ! (vx,vy) de la pared en z=0
    real (kind=GP)    :: vzend_(2)         ! (vx,vy) de la pared en z=Lz
  CONTAINS
    procedure, public :: init          =>          init_impl ! init method
    procedure, public :: dudt          =>          dudt_impl ! RHS method
    procedure, public :: global        =>        global_impl ! Writes global qtys
    procedure, public :: spectra       =>       spectra_impl ! Writes spectra
    procedure, public :: state_size    =>    state_size_impl ! state size
    procedure, public :: get_sstate    =>    get_sstate_impl ! get state name list
    procedure, public :: Solver_ctor   =>      HDSolver_ctor ! constructor
  end type HDSolver

CONTAINS

  ! ===================================================================
  ! Solver initialization, this is where parameter files are read
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to initialize the solver
  !! Reads the &velbound namelist, sets up the boundary
  !! conditions and sets sector indices. Uses the FC-Gram
  !! plan, so it must be called after fcgram_create_plan.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_impl(this)
    use commtypes
    use status
    use grid
    use fft
    class  (HDSolver), intent (inout) :: this

    ! Temporary data to read from namelists:
    character(len=15) :: vbcxsta,vbcxend,vbcysta,vbcyend,vbczsta,vbczend
    real  (kind=GP)   :: vxzsta,vxzend,vyzsta,vyzend
    integer           :: ierr

    ! Required namelists:
    namelist / velbound / vbcxsta,vbcxend,vbcysta,vbcyend,vbczsta,vbczend
    namelist / velbound / vxzsta,vxzend,vyzsta,vyzend

    call MPI_COMM_SIZE(MPI_COMM_WORLD,this%nprocs_,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,this%myrank_,ierr)

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
  !! Todavía no: el paso de RK sigue en el main
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dudt_impl(this, time, uin, uf, dt, dudt)
    class (HDSolver), intent(inout)             :: this
    real   (kind=GP), intent   (in)             :: time, dt
    type(GStateComp), intent(inout), target     :: uin(:),uf(:)
    type(GStateComp), intent(inout), target     :: dudt(:)
    stop 'HDSolver::dudt: todavía no implementado'
  end subroutine dudt_impl


  ! ===================================================================
  ! Computation of global quantities and spectra
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Todavía no: los globales siguen en el main
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine global_impl(this, uin, uf, t)
    use status
    class (HDSolver), intent(in)                :: this
    type(GStateComp), intent(in), target        :: uin(:), uf(:)
    integer         , intent(in)                :: t
    return
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
