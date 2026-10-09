! ===================================================================
! NAME       : equation_base.f90
! DESCRIPTION: Forms base class for all PDEs
!
! Versión de SPECTER del equation_base.f90 de GHOST (rama new_paradigm,
! afaf9ea). Quedan solo EquationBase y VelocityBase (sin escalares
! pasivos), que es lo que usa HD. Se agrega end_stage, que el integrador
! llama después de actualizar cada etapa; por defecto no hace nada, y
! los solvers con paredes lo usan para imponer los bordes y proyectar.
! write_states todavía no está: la escritura de binarios sigue en el
! main.
!
! DATE       : 11/29/25 (DLR)
! ===================================================================

module equationbase_mod
  use class_GWorkspace3D
  use gstate_mod
  use iovar
  implicit none

  private :: Solver_ctor_interface, init_interface
  private :: dudt_interface       , global_interface
  private :: spectra_interface    , state_size_interface

  ! ================= Base class for all PDEs =======================
  ! Define an abstract base class
  type, abstract :: EquationBase
      type(GWorkspace), pointer     :: workspace_ => null()
      type    (ioplan), pointer     :: planio_    => null()
      integer                       :: myrank_     ! MPI rank
      integer                       :: nprocs_     ! MPI rank
      character(len=8), allocatable :: sstate_(:)  ! state member nanes
      character(len=128)            :: infile_     ! config file name
      character(len=128)            :: odir_,idir_ ! internal class I/O directories
      character(len=128)            :: todir_      ! internal class TXT output dir
    contains
      procedure(Solver_ctor_interface), deferred :: Solver_ctor ! Constructor
      procedure(init_interface),        deferred :: init        ! init method
      procedure(dudt_interface),        deferred :: dudt        ! RHS method
      procedure(global_interface),      deferred :: global      ! Global qtys
      procedure(spectra_interface),     deferred :: spectra     ! Spectra
      procedure(state_size_interface),  deferred :: state_size  ! Number of states
      procedure, public                          :: end_stage   ! After each RK stage
      procedure, public                          :: sync_device ! Aux. arrays to the device
  end type EquationBase

  type, abstract, extends(EquationBase) :: VelocityBase
      integer :: VELOCITY    ! start of velocity sector
      integer :: nd_         ! problem dimension
      integer :: nc_         ! # vector field components
  end type VelocityBase

  abstract interface
     subroutine Solver_ctor_interface(this, infile, workspace, plan)
       use class_GWorkspace3D
       use iovar
       import :: EquationBase
       class(EquationBase), intent(inout)         :: this
       type(GWorkspace)   , intent(inout), target :: workspace
       type(ioplan)       , intent(inout), target :: plan
       character(len=*)   , intent   (in)         :: infile
     end subroutine Solver_ctor_interface

     subroutine init_interface(this)
       import :: EquationBase
       class (EquationBase), intent (inout) :: this
     end subroutine init_interface

     subroutine dudt_interface(this, time, uin, uf, dt, dudt)
       use gstate_mod
       import :: EquationBase
       class(EquationBase), intent(inout)         :: this
       real      (kind=GP), intent   (in)         :: time, dt
       type   (GStateComp), intent(inout), target :: uin(:),uf(:)
       type   (GStateComp), intent(inout), target :: dudt(:)
     end subroutine dudt_interface

     subroutine global_interface(this, uin, uf, t)
       use gstate_mod
       use fprecision
       use status
       import :: EquationBase
       class(EquationBase), intent(in)            :: this
       type   (GStateComp), intent(in), target    :: uin(:),uf(:)
       integer            , intent(in)            :: t
     end subroutine global_interface

     subroutine spectra_interface(this, uin)
       use gstate_mod
       use fprecision
       use filefmt
       use status
       import :: EquationBase
       class(EquationBase), intent(in)            :: this
       type   (GStateComp), intent(in), target    :: uin(:)
     end subroutine spectra_interface

     function state_size_interface(this) result(num)
       import :: EquationBase
       class(EquationBase), intent   (in)         :: this
       integer                                    :: num
     end function state_size_interface
  end interface

CONTAINS

  ! ===================================================================
  ! Concrete methods inherited by all solvers
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Lo llama el integrador después de actualizar cada etapa
  !! o (norder, ..., 1) del Runge-Kutta, con el estado ya
  !! actualizado en uout. Por defecto no hace nada.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine end_stage(this, uout, o, norder)
    class (EquationBase), intent(inout)         :: this
    type    (GStateComp), intent(inout), target :: uout(:)
    integer             , intent   (in)         :: o, norder
    return
  end subroutine end_stage


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Copies the auxiliary arrays of the solver (if any) to the
  !! device, after they are set by the initial conditions and the
  !! forcing. Default: the solver has no auxiliary arrays.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sync_device(this)
    class (EquationBase), intent(inout) :: this
    return
  end subroutine sync_device

end module equationbase_mod
