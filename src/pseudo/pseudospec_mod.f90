!=================================================================
! MODULES for 3D codes
!
! 2003 Pablo D. Mininni.
!      Department of Physics,
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar
!=================================================================

!=================================================================

  MODULE grid
      USE fprecision
!
! ni: number of points in the spatial grid in each direction,
! including the continuation points
      INTEGER :: nx
      INTEGER :: ny
      INTEGER :: nz
!
! Ci: number of continuation points, and oi: order at the
! boundaries, in each direction (Ci=0 in periodic directions)
      INTEGER :: Cx
      INTEGER :: Cy
      INTEGER :: Cz
      INTEGER :: ox
      INTEGER :: oy
      INTEGER :: oz
!
! Unidimensional grids
      REAL(KIND=GP), ALLOCATABLE, DIMENSION(:) :: x,y,z
      SAVE
    CONTAINS
!-----------------------------------------------------------------
      SUBROUTINE grid_init(infile_)
!
! Initializes the grid spatial resolution and the FC-Gram
! continuation. Ci and oi are zero (periodic direction) if
! they are not in the namelist.
      USE commtypes
      USE mpivars
      IMPLICIT NONE

      CHARACTER(len=*), INTENT(IN) :: infile_
      NAMELIST/ grid / nx,ny,nz,Cx,Cy,Cz,ox,oy,oz
      Cx = 0; Cy = 0; Cz = 0
      ox = 0; oy = 0; oz = 0
      IF ( myrank .EQ. 0 ) THEN
         OPEN(1,file=infile_,status='unknown',form="formatted")
         READ(1,NML=grid)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(nx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ny,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(Cx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(Cy,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(Cz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ox,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(oy,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(oz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      END SUBROUTINE grid_init

  END MODULE grid

!=================================================================

  MODULE filefmt
!
! Change the length of the string 'ext' to change the number
! of characters used to number the binary files and files with
! the spectra. The format fmtext should be consistent with the
! length of the string, e.g. for a 5 character extent
! len=7 and fmtext = '(i5.5)'.
      CHARACTER(len=4)      :: ext
      CHARACTER(len=6),SAVE :: fmtext = '(i4.4)'

  END MODULE filefmt
!=================================================================

  MODULE fft
      USE fcgram
      TYPE(FCPLAN), SAVE, TARGET :: planfc

  END MODULE fft
!=================================================================

  MODULE var
      USE fprecision
      REAL(KIND=GP), PARAMETER    :: pi = 3.14159265358979323846_GP
      COMPLEX(KIND=GP), PARAMETER :: im = (0.0_GP,1.0_GP)
      SAVE

  END MODULE var
!=================================================================

  MODULE kes
      USE fprecision
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:)     :: kx,ky,kz
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: kk2
      REAL(KIND=GP), ALLOCATABLE, TARGET, DIMENSION (:,:) :: khom
      SAVE
  END MODULE kes
!=================================================================

  MODULE boxsize
!
! Lx, Ly, Lz are the lengths of the sides of the physical (i.e.
! non-continued) box, either in units of 2.pi (if periodic
! direction) or 1 (non-periodic case). Dkx, Dky, Dkz are the
! spacings between wavenumbers in each direction.
! Como boxsize de GHOST, sin Dkk ni anis (SPECTER no escribe
! espectros); las grillas x, y, z siguen en grid y los números
! de onda en kes.
      USE fprecision
      REAL(KIND=GP)    :: Lx,Ly,Lz
      REAL(KIND=GP)    :: Dkx,Dky,Dkz
    CONTAINS
!-----------------------------------------------------------------
      SUBROUTINE box_init(infile_)
!
! Reads the box size from the namelist 'boxparams' and builds the
! grids in physical and Fourier space.
! For periodic diractions the length is in units of 2.pi (=1 gives a
! side of length 2.pi). For non-periodic directions it is measured
! in units of 1 (=1 gives a side length of 1).
!     Lx  : Length in x
!     Ly  : Length in y
!     Lz  : Length in z
      USE commtypes
      USE mpivars
      USE grid
      USE kes
      USE var
      USE gmem
!$    USE threads
      IMPLICIT NONE

      CHARACTER(len=*), INTENT(IN) :: infile_
      REAL(KIND=GP)                :: rmp
      INTEGER                      :: i,j,k
      NAMELIST / boxparams / Lx,Ly,Lz

      IF (myrank.eq.0) THEN
         OPEN(1,file=infile_,status='unknown',form="formatted")
         READ(1,NML=boxparams)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(Lx,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(Ly,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(Lz,1,GC_REAL,0,MPI_COMM_WORLD,ierr)

      ! Grillas y números de onda, con su copia en el dispositivo en
      ! las compilaciones para GPU (los kernels los leen de ahí), como
      ! en GHOST. galloc no tiene la forma (ny,ista:iend) de khom, así
      ! que su copia se registra directamente
      CALL galloc( x, nx )
      CALL galloc( y, ny )
      CALL galloc( z, nz )
      CALL galloc( kx, nx )
      CALL galloc( ky, ny )
      CALL galloc( kz, nz )
      CALL galloc( kk2, nz, ny, ista, iend )
      ALLOCATE( khom(ny,ista:iend) )
      CALL gdev_alloc(C_LOC(khom),INT(SIZE(khom),C_SIZE_T)*INT(STORAGE_SIZE(khom)/8,C_SIZE_T))

! Construct grids
      IF ( Cx .eq. 0 ) THEN
         ! Periodic in x
         rmp = Lx*2.0_GP*pi/nx
         Dkx = 1.0_GP/Lx
      ELSE
         ! Non-periodic in x
         rmp = Lx/(nx-Cx-1)
         Dkx = 2.0_GP*pi/(rmp*nx)
      ENDIF
      DO i=1,nx
        x(i) = rmp*(i-1)
      ENDDO

      IF ( Cy .eq. 0 ) THEN
         ! Periodic in y
         rmp = Ly*2.0_GP*pi/ny
         Dky = 1.0_GP/Ly
      ELSE
         ! Non-periodic in y
         rmp = Ly/(ny-Cy-1)
         Dky = 2.0_GP*pi/(rmp*ny)
      ENDIF
      DO j=1,ny
        y(j) = rmp*(j-1)
      ENDDO

      IF ( Cz .eq. 0 ) THEN
         ! Periodic in z
         rmp = Lz*2.0_GP*pi/nz
         Dkz = 1.0_GP/Lz
      ELSE
        ! Non-periodic in z
        rmp = Lz/(nz-Cz-1)
        Dkz = 2.0_GP*pi/(rmp*nz)
      ENDIF
      DO k=1,nz
        z(k) = rmp*(k-1)
      ENDDO
      pkend = min(nz-Cz,kend)  ! kend of physical domain

      ! Save grid
      IF ( myrank .eq. 0 ) THEN
         OPEN(1,file='x.txt',action='write')
         OPEN(2,file='y.txt',action='write')
         OPEN(3,file='z.txt',action='write')
         DO i = 1,nx-Cx
            WRITE(1,FMT='(1P E23.15)') x(i)
         ENDDO
         DO i = 1,ny-Cy
            WRITE(2,FMT='(1P E23.15)') y(i)
         ENDDO
         DO i = 1,nz-Cz
            WRITE(3,FMT='(1P E23.15)') z(i)
         ENDDO
         CLOSE(1)
         CLOSE(2)
         CLOSE(3)
      ENDIF

! Construct grids for wavenumber domain
      DO i = 1,nx/2
         kx(i) = real(i-1,kind=GP)
         kx(i+nx/2) = real(i-nx/2-1,kind=GP)
      END DO
      IF (ny.eq.1) THEN
        ky(1) = 0.0_GP
      ENDIF
      DO j = 1,ny/2
         ky(j) = real(j-1,kind=GP)
         ky(j+ny/2) = real(j-ny/2-1,kind=GP)
      END DO
      DO k = 1,nz/2
         kz(k) = real(k-1,kind=GP)
         kz(k+nz/2) = real(k-nz/2-1,kind=GP)
      END DO
      kx = kx*Dkx
      ky = ky*Dky
      kz = kz*Dkz

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               kk2(k,j,i) = kx(i)**2+ky(j)**2+kz(k)**2
            END DO
            khom(j,i) = SQRT(kk2(1,j,i))
         END DO
      END DO

      ! Copias en el dispositivo (no hacen nada en las compilaciones
      ! para CPU)
      CALL gupdate_to( x )
      CALL gupdate_to( y )
      CALL gupdate_to( z )
      CALL gupdate_to( kx )
      CALL gupdate_to( ky )
      CALL gupdate_to( kz )
      CALL gupdate_to( kk2 )
      CALL gdev_update_to(C_LOC(khom),INT(SIZE(khom),C_SIZE_T)*INT(STORAGE_SIZE(khom)/8,C_SIZE_T))
      END SUBROUTINE box_init

  END MODULE boxsize
!=================================================================

  MODULE random
      USE var
      USE fprecision
      CONTAINS
!-----------------------------------------------------------------
!-----------------------------------------------------------------

       REAL(KIND=GP) FUNCTION randu(idum)
!
! Uniform distributed random numbers between -1 and
! 1. The seed idum must be between 0 and the value
! of mask

       INTEGER, PARAMETER :: iq=127773,ir=2836,mask=123459876
       INTEGER, PARAMETER :: ia=16807,im=2147483647
       INTEGER            :: k,idum
       REAL(KIND=GP), PARAMETER :: am=1./im

       idum = ieor(idum,mask)
       k = idum/iq
       idum = ia*(idum-k*iq)-ir*k
       IF (idum.lt.0) idum = idum+im
       randu = am*idum
       randu = (randu-.5)*2
       idum = ieor(idum,mask)
       END FUNCTION randu
!-----------------------------------------------------------------
!-----------------------------------------------------------------

       REAL(KIND=GP) FUNCTION randn(idum)
!
! Normally distributed random numbers with zero mean
! and unit variance. The seed idum must be between 0
! and the value of mask in randu.

       REAL(KIND=GP)      :: v1,v2,ran1
       REAL(KIND=GP)      :: fac,rsq
       REAL(KIND=GP),SAVE :: gset
       INTEGER       :: idum
       INTEGER, SAVE :: iset

       IF ((iset.ne.0).or.(iset.ne.1)) iset=0
       IF (idum.lt.0) iset=0
       IF (iset.eq.0) THEN
          rsq = 2.
          DO WHILE ((rsq.ge.1.).or.(rsq.eq.0.))
             v1 = randu(idum)
             v2 = randu(idum)
             rsq = v1**2+v2**2
          END DO
          fac = sqrt(-2.*log(rsq)/rsq)
          gset = v1*fac
          randn = v2*fac
          iset = 1
       ELSE
          randn = gset
          iset = 0
       ENDIF
       END FUNCTION randn
!-----------------------------------------------------------------
!-----------------------------------------------------------------

       SUBROUTINE randn_cmplx(x,y,idum)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Normally distributed complex random numbers with
! unit variance. The seed idum must be between 0
! and the value of mask in randu.
!-----------------------------------------------------------------

       REAL(KIND=GP), INTENT (OUT) :: x,y
       REAL(KIND=GP)       :: fac,phi
       INTEGER, INTENT(IN) :: idum

       x = randu(idum)
       y = randu(idum)
       do while (y.eq.0)
          y = randu(idum)
       end do
       phi = 2.*pi*abs(x)
       fac = sqrt(-2.*log(abs(y)))
       x = fac*cos(phi)
       y = fac*sin(phi)

       END SUBROUTINE randn_cmplx
!-----------------------------------------------------------------
!-----------------------------------------------------------------

       SUBROUTINE prandom_seed(iseed)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes a unique seed for each MPI task. Note: there is a limit
! of 2^(size(integer)) for the number of tasks that can have a
! unique seed. iseed must be an integer seed, and iseed+myrank
! must not be the value of mask in randu (randu stays at that
! value, and all the words of the seed would be equal).
! (Fixed as in GHOST, rama new_paradigm: the old array constructor
! had size(seed)+1 elements, an out-of-bounds write, and seeds with
! many small or correlated words give biased first numbers.)
!
! Parameters
!     iseed : integer seed
!-----------------------------------------------------------------
       USE fprecision
       USE commtypes
       USE mpivars

       INTEGER, INTENT(INOUT)   :: iseed
       INTEGER, ALLOCATABLE     :: iseed1(:)
       INTEGER                  :: j,k,idum
       REAL(KIND=GP)            :: r

       iseed     = mod(iseed+myrank,abs(huge(0)-iseed)-1)
       CALL random_seed(size=k)
       ALLOCATE (iseed1(k))
! Fills the seed with successive states of randu, since some
! compilers (nvfortran) do not scramble the seed, and seeds with
! many small or correlated words give biased first numbers. The
! first states are discarded because the seeds of adjacent ranks
! start close.
       idum = iseed
       DO j = 1,3
          r = randu(idum)
       END DO
       DO j = 1,k
          r = randu(idum)
          iseed1(j) = idum
       END DO
       CALL random_seed(put=iseed1)
       DEALLOCATE(iseed1)

       END SUBROUTINE prandom_seed
!-----------------------------------------------------------------
!-----------------------------------------------------------------

       SUBROUTINE prandom_number(r)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Generates a random number for this MPI task on the interval [0,1].
! Seed for this must be generated by using prandom_seed call.
!
! Parameters
!     iseed : integer seed
!-----------------------------------------------------------------
       USE fprecision
       USE commtypes
       USE mpivars

       REAL(KIND=GP), INTENT(OUT)   :: r

       CALL random_number(r)

       END SUBROUTINE prandom_number
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  END MODULE random
!=================================================================
