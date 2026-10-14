! ===================================================================
! NAME       : equation_base.f90
! DESCRIPTION: Forms base class for all PDEs
!
! Versión de SPECTER del equation_base.f90 de GHOST (rama new_paradigm,
! afaf9ea). Quedan solo EquationBase y VelocityBase (sin escalares
! pasivos), que es lo que usa HD. Se agrega end_stage, que el integrador
! llama después de actualizar cada etapa; por defecto no hace nada, y
! los solvers con paredes lo usan para imponer los bordes y proyectar.
! También se agrega sync_host, la contraparte de sync_device: copia al
! host, antes de las salidas, los arreglos auxiliares que el solver
! cambia en el dispositivo durante el paso (en HD, la presión).
! Los solvers con salidas propias (como la presión en HD) redefinen
! write_states y llaman antes al procedimiento write_states de este
! módulo, que escribe los estados y la vorticidad.
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
  private :: state_size_interface

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
      procedure(state_size_interface),  deferred :: state_size  ! Number of states
      procedure, public                          :: end_stage   ! After each RK stage
      procedure, public                          :: write_states
      procedure, public                          :: sync_device ! Aux. arrays to the device
      procedure, public                          :: sync_host   ! Aux. arrays to the host
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


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Copia al host los arreglos auxiliares del solver que
  !! cambian en el dispositivo durante el paso, antes de las
  !! salidas (propio de SPECTER). Por defecto el solver no
  !! tiene arreglos auxiliares.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sync_host(this)
    class (EquationBase), intent(inout) :: this
    return
  end subroutine sync_host


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Concrete method to write field states: every state
  !! component with its name and, for solvers with a velocity
  !! field and outs >= 1, the vorticity. The fields are
  !! written in the physical domain, through the FC-Gram plan.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_states(this, uin, planio)
    use grid
    use iovar
    use status
    use filefmt
    use fft
    use commtypes
    implicit none

    class (EquationBase), intent   (in)             :: this
    type    (GStateComp), intent(inout)             :: uin(:)
    type        (ioplan), intent   (in)             :: planio
    complex    (kind=GP), pointer, dimension(:,:,:) :: C1,C2,C3
    complex    (kind=GP), pointer, dimension(:,:,:) :: C4,C5,C6
    real       (kind=GP), pointer, dimension(:,:,:) :: R1,R2,R3
    real       (kind=GP)                            :: rmp
    integer                          :: i,j,k,state_size,nc
    logical                          :: bret

    WRITE(ext, fmtext) tind
    call this%workspace_%get_complex_tmp(C1,bret)
    call this%workspace_%get_real_tmp   (R1,bret)
    state_size = this%state_size()
    do nc = 1,state_size
      rmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
!$omp parallel do collapse(2) private (k)
      do i = ista,iend
         do j = 1,ny
            do concurrent (k=1:nz)
              C1(k,j,i) = uin(nc)%ccomp(k,j,i)*rmp
            end do
         end do
      end do
      call fftp3d_complex_to_real(planfc,C1,R1,MPI_COMM_WORLD)
      call io_write(1,this%odir_,trim(this%sstate_(nc)),ext,planio,R1)
    end do
    if ( outs .ge. 1) then
      call this%workspace_%get_complex_tmp(C2,bret)
      call this%workspace_%get_complex_tmp(C3,bret)
      call this%workspace_%get_complex_tmp(C4,bret)
      call this%workspace_%get_complex_tmp(C5,bret)
      call this%workspace_%get_complex_tmp(C6,bret)
      call this%workspace_%get_real_tmp   (R2,bret)
      call this%workspace_%get_real_tmp   (R3,bret)
      select type (this)
      class is (VelocityBase)
!$omp parallel do collapse(2) private (k)
        do i = ista,iend
           do j = 1,ny
              do concurrent (k=1:nz)
                C1(k,j,i) = uin(this%VELOCITY  )%ccomp(k,j,i)*rmp
                C2(k,j,i) = uin(this%VELOCITY+1)%ccomp(k,j,i)*rmp
                C3(k,j,i) = uin(this%VELOCITY+2)%ccomp(k,j,i)*rmp
              end do
           end do
        end do
        call curlk(C2,C3,C4,1)
        call curlk(C1,C3,C5,2)
        call curlk(C1,C2,C6,3)
        call fftp3d_complex_to_real(planfc,C4,R1,MPI_COMM_WORLD)
        call fftp3d_complex_to_real(planfc,C5,R2,MPI_COMM_WORLD)
        call fftp3d_complex_to_real(planfc,C6,R3,MPI_COMM_WORLD)
        call io_write(1,this%odir_,'wx',ext,planio,R1)
        call io_write(1,this%odir_,'wy',ext,planio,R2)
        call io_write(1,this%odir_,'wz',ext,planio,R3)
      end select
      call this%workspace_%free_complex_tmp(C2)
      call this%workspace_%free_complex_tmp(C3)
      call this%workspace_%free_complex_tmp(C4)
      call this%workspace_%free_complex_tmp(C5)
      call this%workspace_%free_complex_tmp(C6)
      call this%workspace_%free_real_tmp   (R2)
      call this%workspace_%free_real_tmp   (R3)
    endif
    call this%workspace_%free_complex_tmp(C1)
    call this%workspace_%free_real_tmp   (R1)
  end subroutine write_states

end module equationbase_mod
