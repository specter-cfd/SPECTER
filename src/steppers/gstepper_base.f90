! ===================================================================
! NAME       : gstepper_base.f90
! DESCRIPTION: Forms base class for all time steppers
!
! Versión de SPECTER del gstepper_base.f90 de GHOST (rama new_paradigm,
! cb6943d), sin partículas: solo el paso de los campos (step), y el
! genérico gstep agrupa solo step.
!
! DATE       : 2/8/26 (DLR)
! ===================================================================

module gstepperbase_mod
  use equationbase_mod, only: EquationBase, GWorkspace
  use gstate_mod
  implicit none

  ! ================= Stepper traits ===================================
  type, public  :: GStepperTraits
    integer            :: itype       ! =1 Butcher type 
    integer            :: norder      ! Order
    integer            :: nstage      ! Number of stages
    integer            :: nstate  = 0 ! Number of field state components
    character(len=128) :: sname       ! Stepper name
  end type GStepperTraits

  ! ================= Base class for all steppers =======================
  ! Define an abstract base class
  type, abstract :: GStepperBase
      type   (GWorkspace), pointer   :: workspace_ => null()
      class(EquationBase), pointer   :: solver_    => null()
      character (len=128)            :: infile_    ! config file name
    contains
      procedure(GStepper_ctor_interface), deferred :: GStepper_ctor ! Constructor
      procedure         (init_interface), deferred :: init          ! init method
      procedure         (step_interface), deferred :: step          ! fluid step method
      generic              , public :: gstep => step                ! general method
  end type GStepperBase 
  
  abstract interface
    subroutine GStepper_ctor_interface(this, traits, workspace, solver)
      import :: GStepperBase, GStepperTraits, GWorkspace, EquationBase
      class (GStepperBase), intent(inout)         :: this
      type(GStepperTraits), intent(inout)         :: traits
      type    (GWorkspace), intent(inout), target :: workspace
      class (EquationBase), intent   (in), target :: solver
    end subroutine GStepper_ctor_interface

    subroutine init_interface(this, traits) 
      import :: GStepperBase, GStepperTraits
      class (GStepperBase), intent(inout) :: this
      type(GStepperTraits), intent   (in) :: traits
    end subroutine init_interface

    subroutine step_interface(this, time, uin, uf, dt, uout)
      import :: GStepperBase, GStateComp, GP
      class(GStepperBase), intent(inout) :: this
      real      (kind=GP), intent   (in) :: time, dt
      type   (GStateComp), intent(inout), target :: uin(:),uf(:)
      type   (GStateComp), intent(inout), target :: uout(:) 
    end subroutine step_interface
  end interface

end module gstepperbase_mod
