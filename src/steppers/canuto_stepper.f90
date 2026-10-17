! =====================================================================
! NAME       : canuto_stepper.f90
!
!              Performs time stepping using the modified RK method 
!              of Canuto et al., Spectral Methods in Fluid
!              Dynamics. This is a low storage time stepper
!              for a specifiable order, though it is strictly
!              speaking of full truncation order only for
!              norder = 2, or for linear ODEs/PDEs.
!
! Versión de SPECTER del canuto_stepper.f90 de GHOST (rama
! new_paradigm, cb6943d): solo el paso de los campos, y después de
! actualizar cada etapa llama a solver%end_stage, donde cada solver
! aplica su corrección (en HD con paredes, imponer los bordes y
! proyectar).
!
! DATE       : 2/8/26 (DLR)
! =====================================================================

module canuto_stepper_mod
  use gstepperbase_mod
  implicit none

  ! ================= Stepper ==========================================
  ! Define class:
  type, extends(GStepperBase) :: CanutoStepper
    ! Member data (extends GStepperBase member data)
    type(GStepperTraits)      :: traits_                     ! GStepper traits
  contains
    procedure, public :: init          => init_impl          ! initialize
    procedure, public :: step          => step_impl          ! take one field timestep
    procedure, public :: GStepper_ctor => CanutoStepper_ctor
    final             :: CanutoStepper_dtor
  end type CanutoStepper

contains

  ! ===================================================================
  ! Stepper-specific methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Constructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine CanutoStepper_ctor(this, traits, workspace, solver)
    class(CanutoStepper), intent(inout)         :: this
    type(GStepperTraits), intent(inout)         :: traits
    type    (GWorkspace), intent(inout), target :: workspace
    class (EquationBase), intent   (in), target :: solver
    this%workspace_ => workspace
    this%solver_    => solver
    if (.not. associated(this%workspace_)) then
      stop 'CanutoStepper::CanutoStepper_ctor: Worskpace not associated'
    endif
    call this%init(traits)
  end subroutine CanutoStepper_ctor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Destructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine CanutoStepper_dtor(this)
    type(CanutoStepper), intent(inout) :: this
    if (associated(this%workspace_)) nullify(this%workspace_)
    if (associated(this%solver_))    nullify(this%solver_)
  end subroutine CanutoStepper_dtor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to initialize the stepper
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_impl(this, traits)
    class(CanutoStepper), intent(inout) :: this
    type(GStepperTraits), intent   (in) :: traits
    this%traits_ = traits;
  end subroutine init_impl


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Implementation function to take one step of fields
  !! uin : state at the beginning of the step (not modified)
  !! uout: at input, a copy of uin; at output, the new state
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine step_impl(this, time, uin, uf, dt, uout)
    use gdevice, only: gdev_active
    implicit none

    class(CanutoStepper), intent(inout) :: this
    type    (GStateComp), intent(inout), target :: uin(:), uf(:), uout(:)
    real       (kind=GP), intent   (in) :: time, dt
    real       (kind=GP)                :: eff_dt
    integer                             :: i,j,k,o,ic
    complex    (kind=GP), pointer       :: po(:,:,:), pi(:,:,:)
       
    if ( size (uin) .ne. this%traits_%nstate &
     .or.size(uout) .ne. this%traits_%nstate  ) then
      stop 'CanutoStepper::step: Inconsistent input state'
    endif

    do o = this%traits_%norder,1,-1
      eff_dt = dt/real(o,kind=GP)
      ! We assume that at input, uout has a copy of uin
      call this%solver_%dudt(time, uout, uf, eff_dt, uout)
      ! The components are addressed through pointers: indexing
      ! uout(ic)%ccomp inside a target region does not work with flang
      do ic = 1,this%traits_%nstate
        po => uout(ic)%ccomp
        pi => uin (ic)%ccomp
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (k)
#endif
        do i = ista,iend
          do j = 1,ny
            do k = 1,nz
              po(k,j,i) = pi(k,j,i) + eff_dt*po(k,j,i)
            end do
          end do
        end do
      end do
      ! Correction of the stage by the solver (in HD with walls:
      ! boundary conditions and projection)
      call this%solver_%end_stage(uout, o, this%traits_%norder)
    end do
  end subroutine step_impl

end module canuto_stepper_mod
