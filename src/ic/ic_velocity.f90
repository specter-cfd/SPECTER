! =====================================================================
! NAME       : ic_velocity.f90
! DESCRIPTION: Initial conditions for all VelocityBase solver classes.
!              These initial conditions can be used to initialize
!              velocities in other solver classes, using composition of
!              initial condition operators, as long as no correlation
!              between the fields is needed.
!
! Initial conditions avaliable:
!   read_v       : Reads velocity from an input file numbered by stat
!   null_v       : Null velocity field
!   cbnoise_v    : Incompressible noise that satisfies no-slip at the
!                  walls (Clever & Busse)
!   couette_v    : Adds the laminar Couette profile to vx
!   poiseuille_v : Adds the laminar Poiseuille profile to vx
!
! Versión de SPECTER del ic_velocity.f90 de GHOST (rama new_paradigm,
! afaf9ea). read_v hace lo mismo que hacía el main de SPECTER al
! continuar una corrida. Respecto del de GHOST, transforma con planfc
! (el plan del FC-Gram) en vez de planrc, y con stat=1 escribe además
! la velocidad leída como la salida 0000. null_v es el de GHOST.
! cbnoise_v es el ruido de examples/initialv.f90 de SPECTER y, con
! rphase, el de examples/initialv.f90_random. couette_v y poiseuille_v
! suman su perfil a lo que ya hay, así que van después de otra
! condición inicial en iclist: "cbnoise_v;couette_v" es el
! initialv.f90_randcou del SPECTER adaptado, "cbnoise_v;poiseuille_v"
! el _randpois y "cbnoise_v;poiseuille_v;couette_v" el _randcpf.
! tg_v, abc_v y random_v de GHOST no están: suponen un dominio
! periódico en las tres direcciones.
!
! DATE       : 01/16/26 (PDM)
! =====================================================================

module ic_velocity
  USE icbase_mod

  IMPLICIT NONE

  ! ================= Initial conditions supported ====================
  type, extends(icBase) :: icRead_v
    contains
      procedure :: init_GState => init_readv
  end type icRead_v
  type, extends(icBase) :: icNull_v
    contains
      procedure :: init_GState => init_nullv
  end type icNull_v
  type, extends(icBase) :: icCbnoise_v
    contains
      procedure :: init_GState => init_cbnoisev
  end type icCbnoise_v
  type, extends(icBase) :: icCouette_v
    contains
      procedure :: init_GState => init_couettev
  end type icCouette_v
  type, extends(icBase) :: icPoiseuille_v
    contains
      procedure :: init_GState => init_poiseuillev
  end type icPoiseuille_v

CONTAINS

  ! ===================================================================
  ! Initial conditions
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Read the velocity field from restart files
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_readv(this,solver,state)
    use gstate_mod
    use hd_mod
    use iovar
    use status
    use filefmt
    use fft
    use commtypes
    implicit none

    class    (icRead_v), intent   (in)          :: this
    class(EquationBase), intent(inout)          :: solver
    type   (GStateComp), intent(inout)          :: state(:)
    complex(kind=GP), pointer, dimension(:,:,:) :: C1
    real   (kind=GP), pointer, dimension(:,:,:) :: R1
    real   (kind=GP)                            :: rmp
    integer                                     :: i,j,k,nc
    logical                                     :: bret

    if ((stat .eq. 0).and.(solver%myrank_ .eq. 0)) then
      stop 'Cannot read files if starting a new run with stat=0'
    endif
    call solver%workspace_%get_complex_tmp(C1,bret)
    call solver%workspace_%get_real_tmp(R1,bret)
    select type (solver)
    class is (VelocityBase)
      tind = int(stat)
      WRITE(ext, fmtext) tind
      do nc = solver%VELOCITY,solver%VELOCITY+solver%nc_-1
        call io_read(1,solver%idir_,trim(solver%sstate_(nc)),ext,solver%planio_,R1)
        call fftp3d_real_to_complex(planfc,R1,state(nc)%ccomp,MPI_COMM_WORLD)
      end do
      ! Con stat=1, la velocidad leída se escribe también como la salida 0000
      if ( stat .eq. 1 ) then
        rmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
        do nc = solver%VELOCITY,solver%VELOCITY+solver%nc_-1
!$omp parallel do collapse(2) private (k)
          do i = ista,iend
             do j = 1,ny
                do concurrent (k=1:nz)
                  C1(k,j,i) = state(nc)%ccomp(k,j,i)*rmp
                end do
             end do
          end do
          call fftp3d_complex_to_real(planfc,C1,R1,MPI_COMM_WORLD)
          call io_write(1,solver%odir_,trim(solver%sstate_(nc)),'0000',solver%planio_,R1)
        end do
      endif
    class default
      error stop "This solver does not support velocity field ICs"
    end select
    call solver%workspace_%free_complex_tmp(C1)
    call solver%workspace_%free_real_tmp(R1)
  end subroutine init_readv


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Null velocity
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_nullv(this,solver,state)
    use gstate_mod
    use hd_mod
    use grid
    use mpivars
    implicit none

    class    (icNull_v), intent   (in) :: this
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
      error stop "This solver does not support velocity field ICs"
    end select
  end subroutine init_nullv


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Ruido incompresible que cumple no deslizamiento en las
  !! paredes, de Clever & Busse ("Tertiary and quaternary
  !! solutions for plane Couette flow"). C1 = dPhi/dz, C2 = Phi
  !! (con las funciones de Chandrasekhar-Reid) y C3 = Psi (con
  !! senos); las derivadas en z se calculan analíticamente.
  !! Como Phi ~ k^-3 y Psi ~ k^-2, v_phi ~ v_psi ~ k^-1.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   u0    : amplitud de la velocidad
  !   kdn   : número de onda perpendicular mínimo
  !   kup   : número de onda perpendicular máximo
  !   nmin  : armónico de Chandrasekhar-Reid mínimo
  !   nmax  : armónico de Chandrasekhar-Reid máximo
  !   rphase: fase de Psi al azar (.true.) o atada a la de Phi
  !           (.false., por defecto)
  subroutine init_cbnoisev(this,solver,state)
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

    class (icCbnoise_v), intent   (in)          :: this
    class(EquationBase), intent(inout)          :: solver
    type   (GStateComp), intent(inout)          :: state(:)
    complex(kind=GP), pointer, dimension(:,:,:) :: C1,C2,C3,C4,C5,C6
    real   (kind=GP)                            :: u0,kdn,kup
    real   (kind=GP)                            :: rm1,rmp,rmq
    real   (kind=GP)                            :: rcr(8)
    integer                                     :: nmin,nmax
    logical                                     :: rphase
    integer                                     :: i,j,k,kk
    logical                                     :: bret

    namelist/ cbnoise_v / u0,kdn,kup,nmin,nmax,rphase
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
      read(1,NML=cbnoise_v)
      close(1)
    endif
    call mpi_bcast(u0    ,1,GC_REAL    ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kdn   ,1,GC_REAL    ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kup   ,1,GC_REAL    ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nmin  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nmax  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(rphase,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

    ! Primeros ocho "números de onda" de Chandrasekhar-Reid
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
                 u0,1)
    class default
      error stop "This solver does not support velocity field ICs"
    end select
    call solver%workspace_%free_complex_tmp(C1)
    call solver%workspace_%free_complex_tmp(C2)
    call solver%workspace_%free_complex_tmp(C3)
    call solver%workspace_%free_complex_tmp(C4)
    call solver%workspace_%free_complex_tmp(C5)
    call solver%workspace_%free_complex_tmp(C6)
  end subroutine init_cbnoisev


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Suma a vx el perfil laminar de Couette, lineal entre las
  !! velocidades de las paredes. Va después de otra condición
  !! inicial en iclist: "cbnoise_v;couette_v" es ruido más
  !! Couette y "null_v;couette_v", el laminar solo.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   vxzsta: velocidad de la pared en z=0
  !   vxzend: velocidad de la pared en z=Lz
  subroutine init_couettev(this,solver,state)
    use gstate_mod
    use hd_mod
    use grid
    use boxsize
    use fft
    use mpivars
    use commtypes
    implicit none

    class (icCouette_v), intent   (in)          :: this
    class(EquationBase), intent(inout)          :: solver
    type   (GStateComp), intent(inout)          :: state(:)
    real   (kind=GP), pointer, dimension(:,:,:) :: R1
    real   (kind=GP)                            :: vxzsta,vxzend
    integer                                     :: i,j,k
    logical                                     :: bret

    namelist/ couette_v / vxzsta,vxzend
    call solver%workspace_%get_real_tmp(R1,bret)
    select type (solver)
    class is (VelocityBase)
    ! Read parameters from a namelist in the input file
    if ( myrank .eq. 0 ) then
      open(1,file=solver%infile_,status='unknown',form="formatted")
      read(1,NML=couette_v)
      close(1)
    endif
    call mpi_bcast(vxzsta,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(vxzend,1,GC_REAL,0,MPI_COMM_WORLD,ierr)

    CALL fftp3d_complex_to_real(planfc,state(solver%VELOCITY)%ccomp, &
                   R1,MPI_COMM_WORLD)
    R1 = R1/nx/ny/nz
!$omp parallel do collapse(2) private (i)
    DO k = ksta,kend
       DO j = 1,ny
          DO CONCURRENT (i=1:nx)
             R1(i,j,k) = R1(i,j,k)+(vxzend-vxzsta)*z(k)/Lz+vxzsta
          END DO
       END DO
    END DO
    CALL fftp3d_real_to_complex(planfc,R1, & !vx
                   state(solver%VELOCITY)%ccomp,MPI_COMM_WORLD)
    class default
      error stop "This solver does not support velocity field ICs"
    end select
    call solver%workspace_%free_real_tmp(R1)
  end subroutine init_couettev


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Suma a vx el perfil laminar de Poiseuille del forzado
  !! uniforme f0 con paredes quietas, f0.z(Lz-z)/(2 nu), de
  !! máximo f0.Lz^2/(8 nu). Va después de otra condición
  !! inicial en iclist: "cbnoise_v;poiseuille_v".
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   f0: amplitud del forzado (la misma que en &uniform_fv)
  !   nu: viscosidad cinemática (la misma que en &HD)
  subroutine init_poiseuillev(this,solver,state)
    use gstate_mod
    use hd_mod
    use grid
    use boxsize
    use fft
    use mpivars
    use commtypes
    implicit none

    class(icPoiseuille_v), intent   (in)        :: this
    class(EquationBase), intent(inout)          :: solver
    type   (GStateComp), intent(inout)          :: state(:)
    real   (kind=GP), pointer, dimension(:,:,:) :: R1
    real   (kind=GP)                            :: f0,nu
    integer                                     :: i,j,k
    logical                                     :: bret

    namelist/ poiseuille_v / f0,nu
    call solver%workspace_%get_real_tmp(R1,bret)
    select type (solver)
    class is (VelocityBase)
    ! Read parameters from a namelist in the input file
    if ( myrank .eq. 0 ) then
      open(1,file=solver%infile_,status='unknown',form="formatted")
      read(1,NML=poiseuille_v)
      close(1)
    endif
    call mpi_bcast(f0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nu,1,GC_REAL,0,MPI_COMM_WORLD,ierr)

    CALL fftp3d_complex_to_real(planfc,state(solver%VELOCITY)%ccomp, &
                   R1,MPI_COMM_WORLD)
    R1 = R1/nx/ny/nz
!$omp parallel do collapse(2) private (i)
    DO k = ksta,kend
       DO j = 1,ny
          DO CONCURRENT (i=1:nx)
             R1(i,j,k) = R1(i,j,k)-f0*z(k)*(z(k)-Lz)/(2*nu)
          END DO
       END DO
    END DO
    CALL fftp3d_real_to_complex(planfc,R1, & !vx
                   state(solver%VELOCITY)%ccomp,MPI_COMM_WORLD)
    class default
      error stop "This solver does not support velocity field ICs"
    end select
    call solver%workspace_%free_real_tmp(R1)
  end subroutine init_poiseuillev

end module ic_velocity
