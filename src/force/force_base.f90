! ===================================================================
! NAME       : force_base.f90
! DESCRIPTION: Forms base class for all forcing methods
!
! Versión de SPECTER del force_base.f90 de GHOST (rama new_paradigm,
! cb6943d), sin blend_forcing: usa saxpby_c de pseudospec_fluid, que
! SPECTER no tiene, y solo lo necesitan los métodos shuffle.
! DATE       : 01/16/26 (PDM)
! ===================================================================

module forcebase_mod
  USE fprecision
  
  IMPLICIT NONE

  ! ================= Base class for all forcings  ==================
  ! Define an abstract base class
  type, abstract :: forceBase
    contains
      procedure(init_GForce_interface), deferred   :: init_GForce
  end type forceBase

  type, abstract :: forceUpdt
      ! The old and new forcing states are CONTIGUOUS so that they are
      ! passed as is to the CONTIGUOUS dummies of gupdate_to: otherwise
      ! compilers may pass a copy (nvfortran 26.5 does), whose address
      ! is not in the registry of device copies of gdevice_mod.
      complex(kind=GP), pointer, contiguous, dimension(:,:,:) :: fxold_ => null()
      complex(kind=GP), pointer, contiguous, dimension(:,:,:) :: fyold_ => null()
      complex(kind=GP), pointer, contiguous, dimension(:,:,:) :: fzold_ => null()
      complex(kind=GP), pointer, contiguous, dimension(:,:,:) :: fxnew_ => null()
      complex(kind=GP), pointer, contiguous, dimension(:,:,:) :: fynew_ => null()
      complex(kind=GP), pointer, contiguous, dimension(:,:,:) :: fznew_ => null()
      logical                                      :: binit_
      ! Set by the update method when it modifies the host copy of the
      ! forcing state in the current step (the main program then copies
      ! the state to the device). Reset by update_forcing before each
      ! call. Methods that change the forcing every step (shuffle) do
      ! the intermediate updates on the device copies and leave the host
      ! copies untouched on those steps.
      logical                                      :: changed_ = .false.
    contains
      procedure(update_GForce_interface), deferred :: update_GForce
      ! Second pass of the update (after all the methods of the chain
      ! ran update_GForce): blend of the old and new states of the
      ! slowly evolving (shuffle) methods. No-op by default.
      procedure                                    :: blend_GForce => blend_GForce_noop
  end type forceUpdt

  ! Concrete data type to do chain operation of forcing methods
  type           :: forceChain
      class(forceBase), allocatable        :: force
      class(forceUpdt), allocatable        :: update
  end type forceChain
 
  abstract interface
    subroutine init_GForce_interface(this, solver, state)
      USE equationbase_mod
      USE gstate_mod
      import :: forceBase
      class   (forceBase),   intent   (in) :: this
      class(EquationBase),   intent(inout) :: solver
      type   (GStateComp),   intent(inout) :: state(:)
    end subroutine
  end interface

  abstract interface
    subroutine update_GForce_interface(this, force, solver, state)
      USE equationbase_mod
      USE gstate_mod
      import :: forceBase
      import :: forceUpdt
      class   (forceUpdt),   intent(inout) :: this
      class   (forceBase),   intent   (in) :: force
      class(EquationBase),   intent(inout) :: solver
      type   (GStateComp),   intent(inout) :: state(:)
    end subroutine
  end interface

CONTAINS
  
  ! ================= Chain operators for ICs =======================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Initializes all forcing states from a list
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_forcing(chain, solver, state)
    USE equationbase_mod
    USE gstate_mod
    implicit none
    type   (forceChain), intent   (in) :: chain(:)
    class(EquationBase), intent(inout) :: solver
    type   (GStateComp), intent(inout) :: state(:)
    integer                            :: i
    do i = 1,size(chain)
      call chain(i)%force%init_GForce(solver,state)
    end do
  end subroutine init_forcing

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Updates all forcing states from a list
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine update_forcing(chain, solver, state)
    USE equationbase_mod
    USE gstate_mod
    USE status
    implicit none
    type   (forceChain), intent(inout) :: chain(:)
    class(EquationBase), intent(inout) :: solver
    type   (GStateComp), intent(inout) :: state(:)
    integer                            :: i
    ! First pass: new random states, phase shifts, etc. (host code, at
    ! the correlation time). A method generating a state correlated with
    ! another forcing (random_fb with corr > 0) sees the newly generated
    ! state of the latter, as in the original code
    do i = 1,size(chain)
      if ( allocated(chain(i)%update) ) then
        chain(i)%update%changed_ = .false.
        call chain(i)%update%update_GForce(chain(i)%force,solver,state)
      endif
    end do
    ! Second pass: blends of the slowly evolving methods (every step)
    do i = 1,size(chain)
      if ( allocated(chain(i)%update) ) then
        call chain(i)%update%blend_GForce(solver,state)
      endif
    end do
    if (timef.eq.fstep) timef = 0 ! Updates state counters
    timef = timef + 1
  end subroutine update_forcing

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Default second pass of the update methods: nothing to do
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine blend_GForce_noop(this, solver, state)
    USE equationbase_mod
    USE gstate_mod
    implicit none
    class   (forceUpdt),   intent(inout) :: this
    class(EquationBase),   intent(inout) :: solver
    type   (GStateComp),   intent(inout) :: state(:)
  end subroutine blend_GForce_noop

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Logical function, true if any update method modified the
  !! host copy of the forcing state in the last call to
  !! update_forcing: the state must then be copied to the
  !! device (GState_update_to) before the time step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical function forcing_changed(chain) result(bchanged)
    implicit none
    type   (forceChain), intent(in) :: chain(:)
    integer                         :: i
    bchanged = .false.
    do i = 1,size(chain)
      if ( allocated(chain(i)%update) ) then
        if ( chain(i)%update%changed_ ) bchanged = .true.
      endif
    end do
  end function forcing_changed

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Logical function to keep track of static forcings:
  !! True if no forcing method has an update scheme, i.e. the
  !! forcing is constant in time (the arrays never change after
  !! init_forcing: the host copies are always current and the
  !! state never needs to be copied to or from a device again)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical function forcing_is_static(chain) result(bstatic)
    implicit none
    type   (forceChain), intent(in) :: chain(:)
    integer                         :: i
    bstatic = .true.
    do i = 1,size(chain)
      if ( allocated(chain(i)%update) ) bstatic = .false.
    end do
  end function forcing_is_static
  
end module forcebase_mod
