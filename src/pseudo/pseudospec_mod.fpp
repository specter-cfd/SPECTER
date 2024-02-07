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
! Number of points in the spatial grid in each direction
      INTEGER, PARAMETER :: nx = NX_
      INTEGER, PARAMETER :: ny = NY_
      INTEGER, PARAMETER :: nz = NZ_

!
! Lx, Ly, Lz, are the lengths of the sides of the physical (i.e.
! non-continued) box, either in units of 2.pi (if periodic
! direction) or 1 (non-periodic case)
      REAL(KIND=GP) :: Lx,Ly,Lz

!
! Global number of continuation points and derivative order.
      INTEGER, PARAMETER :: Cx = CX_
      INTEGER, PARAMETER :: Cy = CY_
      INTEGER, PARAMETER :: Cz = CZ_
      INTEGER, PARAMETER :: ox = OX_
      INTEGER, PARAMETER :: oy = OY_
      INTEGER, PARAMETER :: oz = OZ_


! Unidimensional grids
      REAL(KIND=GP), ALLOCATABLE, DIMENSION(:) :: x,y,z
      SAVE

  END MODULE grid

!=================================================================

  MODULE order
!
! ord: number of iterations in the Runge-Kutta method
      INTEGER :: ord = ORD_
      SAVE

  END MODULE order
!=================================================================

  MODULE filefmt
!
! Change the length of the string 'ext' to change the number
! of characters used to number the binary files and files with
! the spectra. The format fmtext should be consistent with the
! length of the string, e.g. for a 5 character extent
! len=7 and fmtext = '(i5.5)'.
      CHARACTER(len=5)      :: ext
      CHARACTER(len=7),SAVE :: fmtext = '(i5.5)'

  END MODULE filefmt
!=================================================================

  MODULE fft
      USE fcgram
      TYPE(FCPLAN), SAVE :: planfc

  END MODULE fft
!=================================================================

!  MODULE pseudo_boundaries
!      USE boundary

!  END MODULE pseudo_boundaries
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
      REAL(KIND=GP)    :: Dkx,Dky,Dkz
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:)     :: kx,ky,kz
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: kk2
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:)   :: khom
      SAVE
  END MODULE kes
!=================================================================

  MODULE newtmod
      USE fprecision
      INTEGER :: n_dim_1d 
      INTEGER :: n, n_max 
      INTEGER :: step_evol
      INTEGER :: t_rk
      REAL(KIND=GP) :: epsilon = 0.0001_GP, d_sx = 0.0001_GP, d_sy = 0.0000_GP, dT = 0.0001_GP !TODO: Wont assigning values here overwrite the updated values in iteration on every call of newtmod?
      REAL(KIND=GP) :: proj_f, proj_x, proj_y
      REAL(KIND=GP) :: tol
      SAVE

  END MODULE newtmod

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
! Computes unique seed for each MPI task. Note: there is a limit
! of 2^(size(integer)) for number of tasks that can have a unique
! seed
!
! Parameters
!     iseed : integer seed
!-----------------------------------------------------------------
       USE fprecision
       USE commtypes
       USE mpivars

       INTEGER, INTENT(INOUT)   :: iseed
       INTEGER, ALLOCATABLE     :: iseed1(:)
       INTEGER                  :: j,k

       iseed     = mod(iseed+myrank,abs(huge(0)-iseed)-1)
       CALL random_seed(size=k)
       ALLOCATE (iseed1(k), source = iseed*[(j, j=0, k)])
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
