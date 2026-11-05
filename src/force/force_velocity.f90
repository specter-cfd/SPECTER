! =====================================================================
! NAME       : force_velocity.f90
! DESCRIPTION: Forcing methods for all VelocityBase solver classes.
!              Each force needs an initial set up, and an update.
!              These forcing methods can be used to initialize and
!              update forcings on other solver classes, using
!              composition of operators, as long as no correlation
!              between the forces is needed.
!
! Forces avaliable:
!   null_fv    : Null forcing of the velocity field
!   uniform_fv : Uniform forcing along x (mean pressure gradient)
!   cbnoise_fv : Incompressible noise that satisfies no-slip
!
! Update methods available:
!   constant_fv: Constant forcing (no method)
!   shuffle_fv : Slowly evolving random phases
!
! SPECTER version of GHOST's force_velocity.f90 (new_paradigm branch).
! uniform_fv is the forcing of initialfv.f90 of SPECTER (M. Fontana),
! and cbnoise_fv with shuffle_fv the time dependent forcing of SPECTER
! by G. Marchetti.
!
! DATE       : 01/18/26 (PDM)
! =====================================================================

module force_velocity
  USE forcebase_mod

  IMPLICIT NONE

  ! ================= Forcing functions supported =====================
  type, extends(forceBase) :: forceNull_fv
    contains
      procedure :: init_GForce => init_nullfv
  end type forceNull_fv
  type, extends(forceBase) :: forceUniform_fv
    contains
      procedure :: init_GForce => init_uniformfv
  end type forceUniform_fv
  type, extends(forceBase) :: forceCbnoise_fv
    contains
      procedure :: init_GForce => init_cbnoisefv
  end type forceCbnoise_fv
! type, extends(forceBase) :: forceUserdef_fv
!   contains
!     procedure :: init_GForce => init_userdefv
! end type forceUserdef_fv

  ! ================= Update methods supported =======================
  type, extends(forceUpdt) :: shuffleupdt_fv
    contains
      procedure :: update_GForce => update_shufflefv
      procedure :: blend_GForce  => blend_shufflefv
  end type shuffleupdt_fv
! type, extends(forceUpdt) :: userupdt_fv
!   contains
!     procedure :: update_GForce => update_userupdtfv
! end type userupdt_fv

CONTAINS

  ! ===================================================================
  ! Forcing functions
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Null forcing
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_nullfv(this,solver,state)
    use gstate_mod
    use hd_mod
    use grid
    use mpivars
    implicit none

    class(forceNull_fv), intent   (in) :: this
    class(EquationBase), intent(inout) :: solver
    type   (GStateComp), intent(inout) :: state(:)
    integer                            :: i,j,k

    select type (solver)
    class is (VelocityBase)
!$omp parallel do collapse(2) private (k)
      DO i = ista,iend
         DO j = 1,ny
            DO CONCURRENT (k=1:nz)
              state(solver%VELOCITY  )%ccomp(k,j,i) = 0.0_GP
              state(solver%VELOCITY+1)%ccomp(k,j,i) = 0.0_GP
              state(solver%VELOCITY+2)%ccomp(k,j,i) = 0.0_GP
            END DO
         END DO
      END DO
    class default
      error stop "This solver does not support velocity forcing"
    end select
  end subroutine init_nullfv


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Uniform forcing along x (mean pressure gradient)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   f0: amplitude of the forcing
  subroutine init_uniformfv(this,solver,state)
    use gstate_mod
    use hd_mod
    use grid
    use mpivars
    use commtypes
    implicit none

    class(forceUniform_fv), intent   (in) :: this
    class(EquationBase)   , intent(inout) :: solver
    type   (GStateComp)   , intent(inout) :: state(:)
    real   (kind=GP)                      :: f0
    integer                               :: i,j,k

    namelist/ uniform_fv / f0
    select type (solver)
    class is (VelocityBase)
    ! Read parameters from a namelist in the input file
    if ( myrank .eq. 0 ) then
      open(1,file=solver%infile_,status='unknown',form="formatted")
      read(1,NML=uniform_fv)
      close(1)
    endif
    call mpi_bcast(f0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
!$omp parallel do collapse(2) private (k)
      DO i = ista,iend
         DO j = 1,ny
            DO CONCURRENT (k=1:nz)
              state(solver%VELOCITY  )%ccomp(k,j,i) = 0.0_GP
              state(solver%VELOCITY+1)%ccomp(k,j,i) = 0.0_GP
              state(solver%VELOCITY+2)%ccomp(k,j,i) = 0.0_GP
            END DO
         END DO
      END DO
    ! The (0,0,0) mode is in the task with ista = 1
    if ( ista .eq. 1 ) then
      state(solver%VELOCITY)%ccomp(1,1,1) = &
           f0*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
    endif
    class default
      error stop "This solver does not support velocity forcing"
    end select
  end subroutine init_uniformfv


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! 3D BC-satisfying, incompressible noise
  !! From Clever & Busse "Tertiary and quaternary solutions for
  !! plane Couette flow". C1 = dPhi/dz, C2 = Phi, C3 = Psi
  !! For amplitude decay, remember that is Phi~k^-3 => v_phi~k^-1
  !! and Psi~k^-2 => v_psi~k^-1.
  !! Derivatives in z direction are calculated analytically, to
  !! have a more accurate initial condition (specially the no-slip)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   f0    : energy of the forcing, <|f|^2> (see normvec)
  !   kdn   : minimum perpendicular wave number
  !   kup   : maximum perpendicular wave number
  !   nmin  : minimum Chandrasekhar-Reid harmonic
  !   nmax  : maximum Chandrasekhar-Reid harmonic
  !   rphase: random phase of Psi (.true.) or tied to the phase
  !           of Phi (.false., default)
  subroutine init_cbnoisefv(this,solver,state)
    use gstate_mod
    use hd_mod
    use grid
    use boxsize
    use kes
    use var
    use status
    use random
    use fft
    use mpivars
    use commtypes
    implicit none

    class(forceCbnoise_fv), intent   (in)       :: this
    class(EquationBase), intent(inout)          :: solver
    type   (GStateComp), intent(inout)          :: state(:)
    complex(kind=GP), pointer, dimension(:,:,:) :: C1,C2,C3,C4,C5,C6
    real   (kind=GP)                            :: f0,kdn,kup
    real   (kind=GP)                            :: rm1,rmp,rmq
    real   (kind=GP)                            :: rcr(8)
    integer                                     :: nmin,nmax
    logical                                     :: rphase
    integer                                     :: i,j,k,kk
    logical                                     :: bret

    namelist/ cbnoise_fv / f0,kdn,kup,nmin,nmax,rphase
    call solver%workspace_%get_complex_tmp(C1,bret)
    call solver%workspace_%get_complex_tmp(C2,bret)
    call solver%workspace_%get_complex_tmp(C3,bret)
    call solver%workspace_%get_complex_tmp(C4,bret)
    call solver%workspace_%get_complex_tmp(C5,bret)
    call solver%workspace_%get_complex_tmp(C6,bret)
    select type (solver)
    class is (VelocityBase)
    ! Read parameters from a namelist in the input file
    rphase = .false.
    if ( myrank .eq. 0 ) then
      open(1,file=solver%infile_,status='unknown',form="formatted")
      read(1,NML=cbnoise_fv)
      close(1)
    endif
    call mpi_bcast(f0    ,1,GC_REAL    ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kdn   ,1,GC_REAL    ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kup   ,1,GC_REAL    ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nmin  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nmax  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(rphase,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

    ! CR "wavenumbers"
    rcr = (/  4.73004074,  7.85320462, 10.99560784, 14.13716549, &
             17.27875966, 20.42035225, 23.56194490, 26.70353756 /)

    C1 = 0; C2 = 0; C3 = 0
    DO kk = nmin,nmax
    IF ( kk .le. 8 ) THEN
       rm1 = rcr(kk)
    ELSE IF ( MODULO(kk, 2) .eq. 1 ) THEN
       rm1 = (2*kk - 0.5)*pi
    ELSE
       rm1 = (2*kk + 0.5)*pi
    ENDIF

    IF (ista.eq.1) THEN
       DO j = 2,ny/2+1
          IF ((kk2(1,j,1).le.kup**2) .AND. (kk2(1,j,1).ge.kdn**2)) THEN
             rmp = 2*pi*randu(seed)                       ! Phase
             rmq = randu(seed)/sqrt(kk2(1,j,1)+rm1**2)**3 ! Amplitude

             ! Phi is in C2 and C1 is dPhi/dz
             DO k = 1, nz-Cz
                IF ( MODULO(kk, 2) .eq. 1 ) THEN
                C1(k,j,1) = C1(k,j,1) + rmq*rm1/Lz*(COS(rmp)+im*SIN(rmp))*(&
                      SINH(rm1*(z(k)/Lz-0.5))/COSH(0.5*rm1) +&
                      SIN(rm1*(z(k)/Lz-0.5))/COS(0.5*rm1))

                C2(k,j,1) = C2(k,j,1) + rmq*(COS(rmp)+im*SIN(rmp))*(&
                      COSH(rm1*(z(k)/Lz-0.5))/COSH(0.5*rm1) -&
                      COS(rm1*(z(k)/Lz-0.5))/COS(0.5*rm1))
                ELSE
                C1(k,j,1) = C1(k,j,1) + rmq*rm1/Lz*(COS(rmp)+im*SIN(rmp))*(&
                      COSH(rm1*(z(k)/Lz-0.5))/SINH(0.5*rm1) -&
                      COS(rm1*(z(k)/Lz-0.5))/SIN(0.5*rm1))

                C2(k,j,1) = C2(k,j,1) + rmq*(COS(rmp)+im*SIN(rmp))*(&
                      SINH(rm1*(z(k)/Lz-0.5))/SINH(0.5*rm1) -&
                      SIN(rm1*(z(k)/Lz-0.5))/SIN(0.5*rm1))
                ENDIF

                C1(k,ny-j+2,1) = conjg(C1(k,j,1))
                C2(k,ny-j+2,1) = conjg(C2(k,j,1))
             END DO

             ! Psi
             IF ( rphase ) THEN
                rmp = 2*pi*randu(seed)                     ! Phase
             ELSE IF ( MODULO(kk, 2) .eq. 0 ) THEN
                rmp = pi/2 + rmp
             ENDIF
             rmq = randu(seed)/sqrt(kk2(1,j,1)+kk**2*pi**2/Lz**2)**2 ! Amp
             DO k = 1, nz-Cz
                C3(k,j,1) = C3(k,j,1) + rmq*(COS(rmp)+im*SIN(rmp))*SIN(kk*pi*z(k)/Lz)
                C3(k,ny-j+2,1) = conjg(C3(k,j,1))
             END DO
          ENDIF
       ENDDO

       DO i = 2,iend
          DO j = 1,ny
             IF ((kk2(1,j,i).le.kup**2) .AND. (kk2(1,j,i).ge.kdn**2)) THEN
             rmp = 2*pi*randu(seed)                       ! Phase
             rmq = randu(seed)/sqrt(kk2(1,j,i)+rm1**2)**3 ! Amplitude

             ! Phi is in C2 and C1 is dPhi/dz
             DO k = 1, nz-Cz
                IF ( MODULO(kk, 2) .eq. 1 ) THEN
                C1(k,j,i) = C1(k,j,i) + rmq*rm1/Lz*(COS(rmp)+im*SIN(rmp))*(&
                      SINH(rm1*(z(k)/Lz-0.5))/COSH(0.5*rm1) +&
                      SIN(rm1*(z(k)/Lz-0.5))/COS(0.5*rm1))

                C2(k,j,i) = C2(k,j,i) + rmq*(COS(rmp)+im*SIN(rmp))*(&
                      COSH(rm1*(z(k)/Lz-0.5))/COSH(0.5*rm1) -&
                      COS(rm1*(z(k)/Lz-0.5))/COS(0.5*rm1))
                ELSE
                C1(k,j,i) = C1(k,j,i) + rmq*rm1/Lz*(COS(rmp)+im*SIN(rmp))*(&
                      COSH(rm1*(z(k)/Lz-0.5))/SINH(0.5*rm1) -&
                      COS(rm1*(z(k)/Lz-0.5))/SIN(0.5*rm1))

                C2(k,j,i) = C2(k,j,i) + rmq*(COS(rmp)+im*SIN(rmp))*(&
                      SINH(rm1*(z(k)/Lz-0.5))/SINH(0.5*rm1) -&
                      SIN(rm1*(z(k)/Lz-0.5))/SIN(0.5*rm1))
                ENDIF
             END DO

             ! Psi
             IF ( rphase ) THEN
                rmp = 2*pi*randu(seed)                     ! Phase
             ELSE IF ( MODULO(kk, 2) .eq. 0 ) THEN
                rmp = pi/2 + rmp
             ENDIF
             rmq = randu(seed)/sqrt(kk2(1,j,i)+kk**2*pi**2/Lz**2)**2 ! Amp
             DO k = 1, nz-Cz
                C3(k,j,i) = C3(k,j,i) + rmq*(COS(rmp)+im*SIN(rmp))*SIN(kk*pi*z(k)/Lz)
             END DO
             ENDIF
          END DO
       END DO
    ELSE
       DO i = ista,iend
          DO j = 1,ny
             IF ((kk2(1,j,i).le.kup**2) .AND. (kk2(1,j,i).ge.kdn**2)) THEN
             rmp = 2*pi*randu(seed)                        ! Phase
             rmq = randu(seed)/sqrt(kk2(1,j,i)+rm1**2)**3  ! Amplitude

             ! Phi is in C2 and C1 is dPhi/dz
             DO k = 1, nz-Cz
                IF ( MODULO(kk, 2) .eq. 1 ) THEN
                C1(k,j,i) = C1(k,j,i) + rmq*rm1/Lz*(COS(rmp)+im*SIN(rmp))*(&
                      SINH(rm1*(z(k)/Lz-0.5))/COSH(0.5*rm1) +&
                      SIN(rm1*(z(k)/Lz-0.5))/COS(0.5*rm1))

                C2(k,j,i) = C2(k,j,i) + rmq*(COS(rmp)+im*SIN(rmp))*(&
                      COSH(rm1*(z(k)/Lz-0.5))/COSH(0.5*rm1) -&
                      COS(rm1*(z(k)/Lz-0.5))/COS(0.5*rm1))
                ELSE
                C1(k,j,i) = C1(k,j,i) + rmq*rm1/Lz*(COS(rmp)+im*SIN(rmp))*(&
                      COSH(rm1*(z(k)/Lz-0.5))/SINH(0.5*rm1) -&
                      COS(rm1*(z(k)/Lz-0.5))/SIN(0.5*rm1))

                C2(k,j,i) = C2(k,j,i) + rmq*(COS(rmp)+im*SIN(rmp))*(&
                      SINH(rm1*(z(k)/Lz-0.5))/SINH(0.5*rm1) -&
                      SIN(rm1*(z(k)/Lz-0.5))/SIN(0.5*rm1))
                ENDIF
             END DO

             ! Psi
             IF ( rphase ) THEN
                rmp = 2*pi*randu(seed)                     ! Phase
             ELSE IF ( MODULO(kk, 2) .eq. 0 ) THEN
                rmp = pi/2 + rmp
             ENDIF
             rmq = randu(seed)/sqrt(kk2(1,j,i)+kk**2*pi**2/Lz**2)**2 ! Amp
             DO k = 1, nz-Cz
                C3(k,j,i) = C3(k,j,i) + rmq*(COS(rmp)+im*SIN(rmp))*SIN(kk*pi*z(k)/Lz)
             END DO
             ENDIF
          END DO
       END DO
    ENDIF
    END DO

    ! v normal
    CALL derivk(C2,C4,1)
    CALL derivk(C4,C5,1)
    CALL derivk(C2,C4,2)
    CALL derivk(C4,C6,2)
    state(solver%VELOCITY+2)%ccomp = -(C5+C6)                      !vz

    ! v tangential
    ! Poloidal
    CALL derivk(C1,C4,1)
    CALL derivk(C1,C5,2)
    state(solver%VELOCITY  )%ccomp = C4                            !vx
    state(solver%VELOCITY+1)%ccomp = C5                            !vy

    ! Toroidal
    CALL derivk(C3,C4,2)
    CALL derivk(C3,C5,1)
    state(solver%VELOCITY  )%ccomp = state(solver%VELOCITY  )%ccomp+C4
    state(solver%VELOCITY+1)%ccomp = state(solver%VELOCITY+1)%ccomp-C5

    ! 3D Fourier
    CALL fftp1d_real_to_complex_z(planfc, & !vx
                   state(solver%VELOCITY  )%ccomp,MPI_COMM_WORLD)
    CALL fftp1d_real_to_complex_z(planfc, & !vy
                   state(solver%VELOCITY+1)%ccomp,MPI_COMM_WORLD)
    CALL fftp1d_real_to_complex_z(planfc, & !vz
                   state(solver%VELOCITY+2)%ccomp,MPI_COMM_WORLD)

    ! Normalize
    CALL normvec(state(solver%VELOCITY  )%ccomp, &
                 state(solver%VELOCITY+1)%ccomp, &
                 state(solver%VELOCITY+2)%ccomp, &
                 f0,1)
    class default
      error stop "This solver does not support velocity forcing"
    end select
    call solver%workspace_%free_complex_tmp(C1)
    call solver%workspace_%free_complex_tmp(C2)
    call solver%workspace_%free_complex_tmp(C3)
    call solver%workspace_%free_complex_tmp(C4)
    call solver%workspace_%free_complex_tmp(C5)
    call solver%workspace_%free_complex_tmp(C6)
  end subroutine init_cbnoisefv


  ! ===================================================================
  ! Update methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Slowly evolving random phases. This update method
  !! requires randomnes in the generation of the forcing,
  !! i.e., it doesn't work with deterministic forcing as ABC.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine update_shufflefv(this, force, solver, state)
    use equationbase_mod
    use gmem
    use status
    use filefmt
    use fft
    use commtypes
    implicit none
    
    class(shuffleupdt_fv), intent(inout)             :: this
    class     (forceBase), intent   (in)             :: force
    class  (EquationBase), intent(inout)             :: solver
    type     (GStateComp), intent(inout)             :: state(:)
    complex     (kind=GP), pointer, dimension(:,:,:) :: C1,C2,C3
    real        (kind=GP), pointer, dimension(:,:,:) :: R1,R2,R3
    complex(kind=GP)                                 :: cdump
    real(kind=GP)                                    :: phase,rmp
    integer                                          :: i,j,k
    logical                                          :: bret

    ! If not initialized, initialize the method
    if ( .not. this%binit_ ) then
      call solver%workspace_%get_complex_tmp(this%fxold_,bret)
      call solver%workspace_%get_complex_tmp(this%fyold_,bret)
      call solver%workspace_%get_complex_tmp(this%fzold_,bret)
      call solver%workspace_%get_complex_tmp(this%fxnew_,bret)
      call solver%workspace_%get_complex_tmp(this%fynew_,bret)
      call solver%workspace_%get_complex_tmp(this%fznew_,bret)
      if ( stat .eq. 0 ) then ! Generate new forcing
         timef = fstep
      else                    ! Read previous forcing state
        write(ext, fmtext) tind
        call solver%workspace_%get_real_tmp(R1,bret)
        call solver%workspace_%get_real_tmp(R2,bret)
        call solver%workspace_%get_real_tmp(R3,bret)
        call io_read(1,solver%idir_,'fvxold',ext,solver%planio_,R1)
        call io_read(1,solver%idir_,'fvyold',ext,solver%planio_,R2)
        call io_read(1,solver%idir_,'fvzold',ext,solver%planio_,R3)
        call fftp3d_real_to_complex(planfc,R1,this%fxold_,MPI_COMM_WORLD)
        call fftp3d_real_to_complex(planfc,R2,this%fyold_,MPI_COMM_WORLD)
        call fftp3d_real_to_complex(planfc,R3,this%fzold_,MPI_COMM_WORLD)
        call io_read(1,solver%idir_,'fvxnew',ext,solver%planio_,R1)
        call io_read(1,solver%idir_,'fvynew',ext,solver%planio_,R2)
        call io_read(1,solver%idir_,'fvznew',ext,solver%planio_,R3)
        call fftp3d_real_to_complex(planfc,R1,this%fxnew_,MPI_COMM_WORLD)
        call fftp3d_real_to_complex(planfc,R2,this%fynew_,MPI_COMM_WORLD)
        call fftp3d_real_to_complex(planfc,R3,this%fznew_,MPI_COMM_WORLD)
        call solver%workspace_%free_real_tmp(R1)
        call solver%workspace_%free_real_tmp(R2)
        call solver%workspace_%free_real_tmp(R3)
        call gupdate_to(this%fxold_) ! Device copies (no-ops in host builds)
        call gupdate_to(this%fyold_)
        call gupdate_to(this%fzold_)
        call gupdate_to(this%fxnew_)
        call gupdate_to(this%fynew_)
        call gupdate_to(this%fznew_)
      endif    
      this%binit_ = .TRUE.
    endif
    select type (solver)
    class is (VelocityBase)
      ! Generate new forcing states when the correlation time is reached
      if (timef.eq.fstep) then
        ! Keeps a copy of the last forcing state (the host copies are
        ! refreshed first: the blends of the previous steps ran on the
        ! device copies in offload builds)
        call gupdate_from(state(solver%VELOCITY  )%ccomp)
        call gupdate_from(state(solver%VELOCITY+1)%ccomp)
        call gupdate_from(state(solver%VELOCITY+2)%ccomp)
        do i = ista,iend
          do j = 1,ny
            do k = 1,nz
              this%fxold_(k,j,i) = state(solver%VELOCITY  )%ccomp(k,j,i)
              this%fyold_(k,j,i) = state(solver%VELOCITY+1)%ccomp(k,j,i)
              this%fzold_(k,j,i) = state(solver%VELOCITY+2)%ccomp(k,j,i)
            end do
          end do
        end do
        call force%init_GForce(solver,state)     ! Creates a new forcing
        do i = ista,iend                ! Copies the new forcing to fnew
          do j = 1,ny
            do k = 1,nz
              this%fxnew_(k,j,i) = state(solver%VELOCITY  )%ccomp(k,j,i)
              this%fynew_(k,j,i) = state(solver%VELOCITY+1)%ccomp(k,j,i)
              this%fznew_(k,j,i) = state(solver%VELOCITY+2)%ccomp(k,j,i)
            end do
          end do
        end do
        ! Device copies of the old and new states for the blends of
        ! the next fstep steps (no-ops in host builds)
        call gupdate_to(this%fxold_)
        call gupdate_to(this%fyold_)
        call gupdate_to(this%fzold_)
        call gupdate_to(this%fxnew_)
        call gupdate_to(this%fynew_)
        call gupdate_to(this%fznew_)
        this%changed_ = .true. ! Host copies of the state modified
      endif
    class default
      error stop "This solver does not support velocity forcing"
    end select
    ! If we just wrote binary files, we save the force state
    if ((timet.eq.0).and.(bench.eq.0)) then
      write(ext, fmtext) tind
      rmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
      call solver%workspace_%get_complex_tmp(C1,bret)
      call solver%workspace_%get_complex_tmp(C2,bret)
      call solver%workspace_%get_complex_tmp(C3,bret)
      call solver%workspace_%get_real_tmp(R1,bret)
      call solver%workspace_%get_real_tmp(R2,bret)
      call solver%workspace_%get_real_tmp(R3,bret)
!$omp parallel do collapse(2) private (k)
      DO i = ista,iend
         DO j = 1,ny
            DO CONCURRENT (k=1:nz)
              C1(k,j,i) = this%fxold_(k,j,i)*rmp
              C2(k,j,i) = this%fyold_(k,j,i)*rmp
              C3(k,j,i) = this%fzold_(k,j,i)*rmp
            END DO
         END DO
      END DO
      call fftp3d_complex_to_real(planfc,C1,R1,MPI_COMM_WORLD)
      call fftp3d_complex_to_real(planfc,C2,R2,MPI_COMM_WORLD)
      call fftp3d_complex_to_real(planfc,C3,R3,MPI_COMM_WORLD)
      call io_write(1,solver%odir_,'fvxold',ext,solver%planio_,R1)
      call io_write(1,solver%odir_,'fvyold',ext,solver%planio_,R2)
      call io_write(1,solver%odir_,'fvzold',ext,solver%planio_,R3)
!$omp parallel do collapse(2) private (k)
      DO i = ista,iend
         DO j = 1,ny
            DO CONCURRENT (k=1:nz)
              C1(k,j,i) = this%fxnew_(k,j,i)*rmp
              C2(k,j,i) = this%fynew_(k,j,i)*rmp
              C3(k,j,i) = this%fznew_(k,j,i)*rmp
            END DO
         END DO
      END DO
      call fftp3d_complex_to_real(planfc,C1,R1,MPI_COMM_WORLD)
      call fftp3d_complex_to_real(planfc,C2,R2,MPI_COMM_WORLD)
      call fftp3d_complex_to_real(planfc,C3,R3,MPI_COMM_WORLD)
      call io_write(1,solver%odir_,'fvxnew',ext,solver%planio_,R1)
      call io_write(1,solver%odir_,'fvynew',ext,solver%planio_,R2)
      call io_write(1,solver%odir_,'fvznew',ext,solver%planio_,R3)
      call solver%workspace_%free_complex_tmp(C1)
      call solver%workspace_%free_complex_tmp(C2)
      call solver%workspace_%free_complex_tmp(C3)
      call solver%workspace_%free_real_tmp(R1)
      call solver%workspace_%free_real_tmp(R2)
      call solver%workspace_%free_real_tmp(R3)
    endif
  end subroutine update_shufflefv


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Second pass of the shuffle update: slowly updates the
  !! forcing at every time step, f = (1-rmp).fold + rmp.fnew.
  !! The weight is 1/fstep right after a new state is generated
  !! (timef = fstep) and reaches 1 after fstep steps. In the
  !! step of the generation the blend is done on the host
  !! copies (uploaded by the main program with the rest of the
  !! state), in the other steps on the device copies (offload
  !! builds; the host copies are left stale).
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine blend_shufflefv(this, solver, state)
    use equationbase_mod
    use status
    implicit none

    class(shuffleupdt_fv), intent(inout)             :: this
    class  (EquationBase), intent(inout)             :: solver
    type     (GStateComp), intent(inout)             :: state(:)
    real(kind=GP)                                    :: rmp

    select type (solver)
    class is (VelocityBase)
      rmp = float(mod(timef,fstep)+1)/float(fstep)
      call blend_forcing(state(solver%VELOCITY  )%ccomp,this%fxold_, &
                         this%fxnew_,rmp,timef.ne.fstep)
      call blend_forcing(state(solver%VELOCITY+1)%ccomp,this%fyold_, &
                         this%fynew_,rmp,timef.ne.fstep)
      call blend_forcing(state(solver%VELOCITY+2)%ccomp,this%fzold_, &
                         this%fznew_,rmp,timef.ne.fstep)
    class default
      error stop "This solver does not support velocity forcing"
    end select
  end subroutine blend_shufflefv

end module force_velocity
