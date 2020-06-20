!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Subroutines to compute spatial derivatives and nonlinear
! terms in Navier-Stokes equations in 3D using a 
! pseudospectral-FCGram method. You should use the FFTPLANS
! and MPIVARS modules (see the file 'fftp_mod.f90') in each
! program that calls any of the subroutines in this file.
!
! NOTATION: index 'i' is 'x'
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2019 Mauro Fontana & Pablo D. Mininni.
!      Department of Physics,
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mfontana@df.uba.ar
!=================================================================

!*****************************************************************
      SUBROUTINE derivk(a,b,dir)
!-----------------------------------------------------------------
!
! Three-dimensional derivative of the matrix 'a'
!
! Parameters
!     a  : input matrix
!     b  : at the output contains the derivative da/dk_dir
!     dir: =1 derivative in the x-direction
!          =2 derivative in the y-direction
!          =3 derivative in the z-direction
!
      USE kes
      USE var
      USE grid
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: b
      INTEGER, INTENT(IN) :: dir
      INTEGER             :: i,j,k

!
! Derivative in the x-direction
!
      IF (dir.eq.1) THEN
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  b(k,j,i) = im*kx(i)*a(k,j,i)
               END DO
            END DO
         END DO
!
! Derivative in the y-direction
!
      ELSE IF (dir.eq.2) THEN
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  b(k,j,i) = im*ky(j)*a(k,j,i)
               END DO
            END DO
         END DO
!
! Derivative in the z-direction
!
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  b(k,j,i) = im*kz(k)*a(k,j,i)
               END DO
            END DO
         END DO
      ENDIF

      RETURN
      END SUBROUTINE derivk

!*****************************************************************
      SUBROUTINE laplak(a,b)
!-----------------------------------------------------------------
!
! Three-dimensional Laplacian of the matrix 'a'
!
! Parameters
!     a: input matrix
!     b: at the output contains the Laplacian
!
      USE kes
      USE grid
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: b
      INTEGER :: i,j,k

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               b(k,j,i) = -kk2(k,j,i)*a(k,j,i)
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE laplak

!*****************************************************************
      SUBROUTINE curlk(a,b,c,dir)
!-----------------------------------------------------------------
!
! Computes the curl of the vector field A in Fourier
! space. The needed components of the field A are
! given by the matrixes a and b, and the order must
! follow the right hand convention [(x,y), (x,z) or (y,z)].
!
! Parameters
!     a  : input matrix
!     b  : input matrix
!     c  : at the output contains curl(A)_dir
!     dir: =1 computes the x-component
!          =2 computes the y-component
!          =3 computes the z-component
!
      USE fprecision
      USE grid
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a,b
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: c
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: c1,c2
      INTEGER, INTENT(IN) :: dir
      INTEGER             :: i,j,k

!
! Computes the x-component
!
      IF (dir.eq.1) THEN
         CALL derivk(a,c1,3)
         CALL derivk(b,c2,2)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  c(k,j,i) = c2(k,j,i)-c1(k,j,i)
               END DO
            END DO
         END DO
!
! Computes the y-component
!
      ELSE IF (dir.eq.2) THEN
         CALL derivk(a,c1,3)
         CALL derivk(b,c2,1)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  c(k,j,i) = c1(k,j,i)-c2(k,j,i)
               END DO
            END DO
         END DO
!
! Computes the z-component
!
      ELSE
         CALL derivk(a,c1,2)
         CALL derivk(b,c2,1)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  c(k,j,i) = c2(k,j,i)-c1(k,j,i)
               END DO
            END DO
         END DO
      ENDIF

      RETURN
      END SUBROUTINE curlk

!*****************************************************************
      SUBROUTINE gradre(a,b,c,d,e,f)
!-----------------------------------------------------------------
!
! Three-dimensional inner product A.grad(A) in
! real space. The components of the field A are
! given by the matrixes a, b and c.
!
! Parameters
!     a: input matrix in the x-direction
!     b: input matrix in the y-direction
!     c: input matrix in the z-direction
!     d: product (A.grad)A_x in Fourier space [output]
!     e: product (A.grad)A_y in Fourier space [output]
!     f: product (A.grad)A_z in Fourier space [output]
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: d,e,f
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: c1,c2
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: c3,c4
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: r1,r2
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: r3,r4
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: rx,ry,rz
      REAL(KIND=GP)    :: tmp
      INTEGER :: i,j,k

!
! Computes (A_x.dx)A_dir
!
      c1 = a
      CALL derivk(a,c2,1)
      CALL derivk(b,c3,1)
      CALL derivk(c,c4,1)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c3,r3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c4,r4,MPI_COMM_WORLD)

!$omp parallel do if (pkend-ksta.ge.nth) private (j,i)
      DO k = ksta,pkend
!$omp parallel do if (pkend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               rx(i,j,k) = r1(i,j,k)*r2(i,j,k)
               ry(i,j,k) = r1(i,j,k)*r3(i,j,k)
               rz(i,j,k) = r1(i,j,k)*r4(i,j,k)
            END DO
         END DO
      END DO
!
! Computes (A_y.dy)A_dir
!
      c1 = b
      CALL derivk(a,c2,2)
      CALL derivk(b,c3,2)
      CALL derivk(c,c4,2)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c3,r3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c4,r4,MPI_COMM_WORLD)

!$omp parallel do if (pkend-ksta.ge.nth) private (j,i)
      DO k = ksta,pkend
!$omp parallel do if (pkend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               rx(i,j,k) = rx(i,j,k)+r1(i,j,k)*r2(i,j,k)
               ry(i,j,k) = ry(i,j,k)+r1(i,j,k)*r3(i,j,k)
               rz(i,j,k) = rz(i,j,k)+r1(i,j,k)*r4(i,j,k)
            END DO
         END DO
      END DO
!
! Computes (A_z.dz)A_dir
!
      c1 = c
      CALL derivk(a,c2,3)
      CALL derivk(b,c3,3)
      CALL derivk(c,c4,3)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c3,r3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c4,r4,MPI_COMM_WORLD)

      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
!$omp parallel do if (pkend-ksta.ge.nth) private (j,i)
      DO k = ksta,pkend
!$omp parallel do if (pkend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               rx(i,j,k) = (rx(i,j,k)+r1(i,j,k)*r2(i,j,k))*tmp
               ry(i,j,k) = (ry(i,j,k)+r1(i,j,k)*r3(i,j,k))*tmp
               rz(i,j,k) = (rz(i,j,k)+r1(i,j,k)*r4(i,j,k))*tmp
            END DO
         END DO
      END DO

      CALL fftp3d_real_to_complex(planrc,rx,d,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,ry,e,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,rz,f,MPI_COMM_WORLD)

      RETURN
      END SUBROUTINE gradre

!*****************************************************************
      SUBROUTINE prodre(a,b,c,d,e,f)
!-----------------------------------------------------------------
!
! Three-dimensional cross product curl(A)xA in
! real space. The components of the field A are
! given by the matrixes a, b and c
!
! Parameters
!     a  : input matrix in the x-direction
!     b  : input matrix in the y-direction
!     c  : input matrix in the z-direction
!     d  : product [curl(A)xA]_x in Fourier space [output]
!     e  : product [curl(A)xA]_y in Fourier space [output]
!     f  : product [curl(A)xA]_z in Fourier space [output]
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: d,e,f
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: r1,r2
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: r3,r4
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: r5,r6
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: r7
      REAL(KIND=GP)    :: tmp
      INTEGER :: i,j,k

!
! Computes curl(A)
!
      CALL curlk(b,c,d,1)
      CALL curlk(a,c,e,2)
      CALL curlk(a,b,f,3)
      CALL fftp3d_complex_to_real(plancr,d,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,e,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,f,r3,MPI_COMM_WORLD)
!
! Computes A
!
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               d(k,j,i) = a(k,j,i)
               e(k,j,i) = b(k,j,i)
               f(k,j,i) = c(k,j,i)
            END DO
         END DO
      END DO
      CALL fftp3d_complex_to_real(plancr,d,r4,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,e,r5,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,f,r6,MPI_COMM_WORLD)
!
! Computes curl(A)xA
!
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
!$omp parallel do if (pkend-ksta.ge.nth) private (j,i)
      DO k = ksta,pkend
!$omp parallel do if (pkend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               r7(i,j,k) = (r2(i,j,k)*r6(i,j,k)-r5(i,j,k)*r3(i,j,k))*tmp
               r3(i,j,k) = (r3(i,j,k)*r4(i,j,k)-r6(i,j,k)*r1(i,j,k))*tmp
               r1(i,j,k) = (r1(i,j,k)*r5(i,j,k)-r4(i,j,k)*r2(i,j,k))*tmp
            END DO
         END DO
      END DO

      CALL fftp3d_real_to_complex(planrc,r7,d,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,r3,e,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,r1,f,MPI_COMM_WORLD)

      RETURN
      END SUBROUTINE prodre

!*****************************************************************
      SUBROUTINE pressure(a,b,c,d)
!-----------------------------------------------------------------
!     Solves the Poisson equation laplacian(p') = div(v) with 
!     boundary  conditions such that v_z - deriv(p,z) = 0 on the
!     surface. Then applies v = v - grad(p'). The pressure
!     in the mixed (z,ky,kx) domain is returned.
!     p' = dt*p for convenience.
!
! Parameters
!     a, b, c : x,y,z components of v, respectively (kz,ky,kx) [INOUT]
!     d       : solution to the Poisson equation p' (z,ky,kx)  [OUT]
!
!TODO benchmark applying p_i and p_h separatedly in x and y components
! vs applying them together at the end. In z it is necessary to apply
! them in order to generate the boundary condition.
      USE grid
      USE var
      USE kes
      USE fcgram
      USE commtypes
      USE mpivars
      USE fft
!$    USE threads

      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend)   :: d
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)           :: C1,C2,C3

      REAL(KIND=GP)    :: tmp
      INTEGER          :: i,j,k

      ! Solve Poisson
      CALL poisson(a,b,c,d)

      ! Apply inhomogeneous solution
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i=ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j=1,ny
            DO k=1,nz
               a(k,j,i) = a(k,j,i) - im*kx(i)*d(k,j,i)
               b(k,j,i) = b(k,j,i) - im*ky(j)*d(k,j,i)
               c(k,j,i) = c(k,j,i) - im*kz(k)*d(k,j,i)
            ENDDO
         ENDDO
      ENDDO

      ! Boundary condition
      tmp = 1.0_GP/nz
      C1 = c
      CALL fftp1d_complex_to_real_z(plancr,C1,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i=ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j=1,ny
            DO k=1,nz
               C1(k,j,i) = C1(k,j,i)*tmp
            ENDDO
         ENDDO
      ENDDO

      ! Solve Laplace equation and construct p_h(z,ky,kx)
      CALL laplace(C1(1,:,:),C1(nz-Cz,:,:),C2,C3)

      ! Total pressure in (z,ky,kx) domain
      CALL fftp1d_complex_to_real_z(plancr,d,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i=ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j=1,ny
            DO k=1,nz
               d(k,j,i) = d(k,j,i)*tmp + C2(k,j,i)
            ENDDO
         ENDDO
      ENDDO

      ! Hom. pressure and normal derivative in (kz,ky,kx) domain
      CALL fftp1d_real_to_complex_z(planrc,C2,MPI_COMM_WORLD)
      CALL fftp1d_real_to_complex_z(planrc,C3,MPI_COMM_WORLD)

      ! Apply homogeneous solution
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i=ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j=1,ny
            DO k=1,nz
               a(k,j,i) = a(k,j,i) - im*kx(i)*C2(k,j,i)
               b(k,j,i) = b(k,j,i) - im*ky(j)*C2(k,j,i)
               c(k,j,i) = c(k,j,i) - C3(k,j,i)
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE pressure

!*****************************************************************
      SUBROUTINE poisson(a,b,c,d)
!-----------------------------------------------------------------
!
! Solves the inhomogeneous part of the Poisson equation
! laplacian(p') = div(v). It takes the components of v as input
! matrixes and p'_inhom in Fourier space. Note that the boundary
! conditions must be enforced via the homogeneous solution a
! a posteriori.
!
! Parameters
!     a  : input matrix in the x-direction
!     b  : input matrix in the y-direction
!     c  : input matrix in the z-direction
!     d  : at the output contains the result in Fourier space
!
      USE fprecision
      USE var
      USE kes
      USE grid
      USE mpivars
      USE commtypes
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: d
      INTEGER             :: i,j,k

!NOTE: This generates a NaN, but shouldn't propagate to relevant
!variables. Written this way is easier for the compiler to vectorize
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               d(k,j,i) = -im*(kx(i)*a(k,j,i) + ky(j)*b(k,j,i) + &
                               kz(k)*c(k,j,i))/kk2(k,j,i)
            END DO
          END DO
      END DO
      IF (ista.eq.1) d(1,1,1) = 0.0_GP 


      RETURN
      END SUBROUTINE poisson

!*****************************************************************
      SUBROUTINE laplace(bb,tb,a,b)
!-----------------------------------------------------------------
! Returns the gradient of phi in (z,ky,kx) coordinates, where  phi
! is a Laplace equation solution satisfing dphi/dz|z=0 = bb and
! dphi/dz|z=L = tb.

! Parameters
!     bb : input matrix for the bottom boundary condition
!     tb : input matrix for the top boundary condition
!     a  : at the output contains the result in the mixed
!          (z,ky,kx) space

      USE kes
      USE var
      USE grid
      USE fcgram
      USE mpivars
      USE fft
!$    USE threads


      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(ny,ista:iend)     :: bb,tb
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: a,b

      COMPLEX(KIND=GP)             :: coef1, coef2
      REAL(KIND=GP)                :: tmp

      INTEGER                      :: i,j,k

      IF (ista.eq.1) THEN
!$omp parallel do
         DO k=1,nz-Cz
            a(k,1,1) = real(bb(1,1), kind=GP)*z(k)
            b(k,1,1) = real(bb(1,1), kind=GP)
         ENDDO
!$omp parallel do private (k,tmp,coef1,coef2)
         DO j = 2,ny
            tmp = 1.0_GP/(khom(j,1)*(1-exp(-2*khom(j,1)*Lz)))
            coef1 = (tb(j,1)-bb(j,1)*exp(-khom(j,1)*Lz))*tmp
            coef2 = - (bb(j,1)-tb(j,1)*exp(-khom(j,1)*Lz))*tmp
            DO k=1,nz-Cz
               a(k,j,1) = coef1*exp(khom(j,1)*(z(k)-Lz)) + &
                          coef2*exp(-khom(j,1)*z(k))
               b(k,j,1) = khom(j,1)*(coef1*exp(khom(j,1)*(z(k)-Lz)) - &
                          coef2*exp(-khom(j,1)*z(k)))
            END DO
         END DO
!$omp parallel do if (iend-ista.ge.nth) private (j,k,tmp,coef1,coef2)
         DO i = 2,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,tmp,coef1,coef2)
            DO j = 1,ny
               tmp = 1.0_GP/(khom(j,i)*(1-exp(-2*khom(j,i)*Lz)))
               coef1 = (tb(j,i)-bb(j,i)*exp(-khom(j,i)*Lz))*tmp
               coef2 = - (bb(j,i)-tb(j,i)*exp(-khom(j,i)*Lz))*tmp
               DO k=1,nz-Cz
                  a(k,j,i) = coef1*exp(khom(j,i)*(z(k)-Lz)) + &
                             coef2*exp(-khom(j,i)*z(k))
                  b(k,j,i) = khom(j,i)*(coef1*exp(khom(j,i)*(z(k)-Lz)) - &
                             coef2*exp(-khom(j,i)*z(k)))
               END DO
            END DO
         END DO
      ELSE 
!$omp parallel do if (iend-ista.ge.nth) private (j,k,tmp,coef1,coef2)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,tmp,coef1,coef2)
            DO j = 1,ny
               tmp = 1.0_GP/(khom(j,i)*(1-exp(-2*khom(j,i)*Lz)))
               coef1 = (tb(j,i)-bb(j,i)*exp(-khom(j,i)*Lz))*tmp
               coef2 = - (bb(j,i)-tb(j,i)*exp(-khom(j,i)*Lz))*tmp
               DO k=1,nz-Cz
                  a(k,j,i) = coef1*exp(khom(j,i)*(z(k)-Lz)) + &
                             coef2*exp(-khom(j,i)*z(k))
                  b(k,j,i) = khom(j,i)*(coef1*exp(khom(j,i)*(z(k)-Lz)) - &
                             coef2*exp(-khom(j,i)*z(k)))
               END DO
            END DO
         END DO
      ENDIF

      RETURN
      END SUBROUTINE laplace

!*****************************************************************
      SUBROUTINE energy(a,b,c,d,kin)
!-----------------------------------------------------------------
!
! Computes the mean energy of a vector field utilizing Parseval's
! theorem to compute sum(u^2) in (kx,ky,z) domain.
! The output is only valid in the first node.
!
! Parameters
!     a  : input matrix in the x-direction
!     b  : input matrix in the y-direction
!     c  : input matrix in the z-direction
!     d  : at the output contains the energy
!     kin: =0 computes the magnetic energy
!          =1 computes the kinetic energy
!          =2 computes the magnetic enstrophy
!
      USE fprecision
      USE commtypes
      USE grid
      USE fcgram
      USE mpivars
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      INTEGER, INTENT(IN)           :: kin
      DOUBLE PRECISION, INTENT(OUT) :: d

      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: C1
      REAL(KIND=GP), DIMENSION(nz,ny,ista:iend)    :: R1

      DOUBLE PRECISION  :: dloc
      REAL(KIND=GP)     :: tmp
      INTEGER           :: i,j,k


      dloc = 0.0D0
      ! One nx*ny from parsevals, an nx*ny*(nz-Cz) from averaging over
      ! physical domain and a nz**2 due to fft normalization.
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2/ &
            real((nz-Cz),kind=GP)

! Computes the energy of field (a,b,c)
      IF (kin .eq. 1) THEN
         C1 = a
      ELSE IF (kin.eq.0) THEN
         CALL curlk(b,c,c1,1)
      ENDIF
      CALL fftp1d_complex_to_real_z(plancr,c1,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
      DO j = 1,ny
      DO k = 1,nz
         R1(k,j,i) = real(C1(k,j,i), kind=GP)**2+aimag(C1(k,j,i))**2
      ENDDO
      ENDDO
      ENDDO

      IF (kin .eq. 1) THEN
         C1 = b
      ELSE IF (kin.eq.0) THEN
         CALL curlk(a,c,c1,2)
      ENDIF
      c1 = b
      CALL fftp1d_complex_to_real_z(plancr,c1,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
      DO j = 1,ny
      DO k = 1,nz
         R1(k,j,i) = R1(k,j,i) + real(C1(k,j,i), kind=GP)**2 +&
                        aimag(C1(k,j,i))**2
      ENDDO
      ENDDO
      ENDDO

      IF (kin .eq. 1) THEN
         C1 = c
      ELSE IF (kin.eq.0) THEN
         CALL curlk(a,b,c1,3)
      ENDIF
      CALL fftp1d_complex_to_real_z(plancr,c1,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
      DO j = 1,ny
      DO k = 1,nz
         R1(k,j,i) = R1(k,j,i) + real(C1(k,j,i), kind=GP)**2 +&
                        aimag(C1(k,j,i))**2
      ENDDO
      ENDDO
      ENDDO

      IF (ista.eq.1) THEN
!$omp parallel do private (k) reduction(+:dloc)
         DO j = 1,ny
            DO k = 1,nz-Cz
               dloc = dloc + R1(k,j,1)*tmp
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k) reduction(+:dloc)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k) reduction(+:dloc)
            DO j = 1,ny
               DO k = 1,nz-Cz
                  dloc = dloc + 2*R1(k,j,i)*tmp
               END DO
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k) reduction(+:dloc)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k) reduction(+:dloc)
           DO j = 1,ny
              DO k = 1,nz-Cz
                 dloc = dloc + 2*R1(k,j,i)*tmp
              END DO
           END DO
        END DO
      ENDIF

! Computes the reduction between nodes
      CALL MPI_REDUCE(dloc,d,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE energy

!*****************************************************************
      SUBROUTINE helicity(a,b,c,d)
!-----------------------------------------------------------------
!
! Computes the mean helicity of a vector field.
! The output is only valid in the first node.
!
! Parameters
!     a: input matrix in the x-direction
!     b: input matrix in the y-direction
!     c: input matrix in the z-direction
!     d: at the output contains the helicity
!
      USE fprecision
      USE commtypes
      USE fft
      USE grid
      USE fcgram
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      DOUBLE PRECISION, INTENT(OUT) :: d

      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)             :: c1,c2

      DOUBLE PRECISION :: dloc
      REAL(KIND=GP)    :: tmp
      INTEGER          :: i,j,k

      dloc = 0.0D0
      ! One nx*ny from parsevals, an nx*ny*(nz-Cz) from averaging over
      ! physical domain and a nz**2 due to fft normalization.
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2 / &
            real(nz-Cz,kind=GP)

      c1 = a
      CALL curlk(b,c,c2,1)
      CALL fftp1d_complex_to_real_z(plancr,c1,MPI_COMM_WORLD)
      CALL fftp1d_complex_to_real_z(plancr,c2,MPI_COMM_WORLD)
      IF (ista.eq.1) THEN
!$omp parallel do private (k) reduction(+:dloc)
         DO j = 1,ny
            DO k = 1,nz-Cz
                 dloc = dloc+real(c1(k,j,1)*conjg(c2(k,j,1)))*tmp
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k) reduction(+:dloc)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k) reduction(+:dloc)
            DO j = 1,ny
               DO k = 1,nz-Cz
                  dloc = dloc+2*real(c1(k,j,i)*conjg(c2(k,j,i)))*tmp
               END DO
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k) reduction(+:dloc)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k) reduction(+:dloc)
            DO j = 1,ny
               DO k = 1,nz-Cz
                  dloc = dloc+2*real(c1(k,j,i)*conjg(c2(k,j,i)))*tmp
               END DO
            END DO
         END DO
      ENDIF

      C1 = b
      CALL curlk(a,c,c2,2)
      CALL fftp1d_complex_to_real_z(plancr,c1,MPI_COMM_WORLD)
      CALL fftp1d_complex_to_real_z(plancr,c2,MPI_COMM_WORLD)
      IF (ista.eq.1) THEN
!$omp parallel do private (k) reduction(+:dloc)
         DO j = 1,ny
            DO k = 1,nz-Cz
               dloc = dloc+real(c1(k,j,1)*conjg(c2(k,j,1)))*tmp
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k) reduction(+:dloc)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k) reduction(+:dloc)
            DO j = 1,ny
               DO k = 1,nz-Cz
                  dloc = dloc+2*real(c1(k,j,i)*conjg(c2(k,j,i)))*tmp
               END DO
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k) reduction(+:dloc)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k) reduction(+:dloc)
            DO j = 1,ny
               DO k = 1,nz-Cz
                  dloc = dloc+2*real(c1(k,j,i)*conjg(c2(k,j,i)))*tmp
               END DO
            END DO
         END DO
      ENDIF

      C1 = c
      CALL curlk(a,b,c2,3)
      CALL fftp1d_complex_to_real_z(plancr,c1,MPI_COMM_WORLD)
      CALL fftp1d_complex_to_real_z(plancr,c2,MPI_COMM_WORLD)
      IF (ista.eq.1) THEN
!$omp parallel do private (k) reduction(+:dloc)
         DO j = 1,ny
            DO k = 1,nz-Cz
               dloc = dloc+real(c1(k,j,1)*conjg(c2(k,j,1)))*tmp
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k) reduction(+:dloc)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k) reduction(+:dloc)
            DO j = 1,ny
               DO k = 1,nz-Cz
                  dloc = dloc+2*real(c1(k,j,i)*conjg(c2(k,j,i)))*tmp
               END DO
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k) reduction(+:dloc)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k) reduction(+:dloc)
            DO j = 1,ny
               DO k = 1,nz-Cz
                  dloc = dloc+2*real(c1(k,j,i)*conjg(c2(k,j,i)))*tmp
               END DO
            END DO
         END DO
      ENDIF

      CALL MPI_REDUCE(dloc,d,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE helicity

!*****************************************************************
      SUBROUTINE cross(a,b,c,d,e,f,g,kin)
!-----------------------------------------------------------------
!
! Computes the cross helicity or averaged inner
! product of two vector fields. The output is
! only valid in the first node.
!
! Parameters
!     a  : first field x-component
!     b  : first field y-component
!     c  : first field z-component
!     d  : second field x-component
!     e  : second field y-component
!     f  : second field z-component
!     g  : at the output contains the inner product
!     kin: =0 computes the inner product of the curls
!          =1 computes the inner product of the fields
!
      USE fprecision
      USE commtypes
      USE kes
      USE fft
      USE grid
      USE fcgram
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: d,e,f
      DOUBLE PRECISION, INTENT(OUT) :: g

      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)             :: c1,c2
      DOUBLE PRECISION, DIMENSION(nz,ny,ista:iend)             :: r1

      DOUBLE PRECISION              :: gloc
      REAL(KIND=GP)                 :: tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j,k

      gloc = 0.0D0
      ! One nx*ny from parsevals, an nx*ny*(nz-Cz) from averaging over
      ! physical domain and a nz**2 due to fft normalization.
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2/ &
            real(nz-Cz, kind=GP)

! Computes the averaged inner product between the fields
      IF (kin.eq.1) THEN
         c1 = a
         c2 = d
         CALL fftp1d_complex_to_real_z(plancr,c1,MPI_COMM_WORLD)
         CALL fftp1d_complex_to_real_z(plancr,c2,MPI_COMM_WORLD)
         r1 = real(c1*conjg(c2))

         c1 = b
         c2 = e
         CALL fftp1d_complex_to_real_z(plancr,c1,MPI_COMM_WORLD)
         CALL fftp1d_complex_to_real_z(plancr,c2,MPI_COMM_WORLD)
         r1 = r1 + real(c1*conjg(c2))

         c1 = c
         c2 = f
         CALL fftp1d_complex_to_real_z(plancr,c1,MPI_COMM_WORLD)
         CALL fftp1d_complex_to_real_z(plancr,c2,MPI_COMM_WORLD)
         r1 = r1 + real(c1*conjg(c2))

         IF (ista.eq.1) THEN
!$omp parallel do private (k) reduction(+:gloc)
            DO j = 1,ny
               DO k = 1,nz-Cz
                  gloc = gloc+r1(k,j,1)*tmp
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k) reduction(+:gloc)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k) reduction(+:gloc)
               DO j = 1,ny
                  DO k = 1,nz-Cz
                     gloc = gloc+2*r1(k,j,i)*tmp
                  END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k) reduction(+:gloc)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k) reduction(+:gloc)
               DO j = 1,ny
                  DO k = 1,nz-Cz
                     gloc = gloc+2*r1(k,j,i)*tmp
                  END DO
               END DO
            END DO
         ENDIF
      ENDIF

! Computes the reduction between nodes
      CALL MPI_REDUCE(gloc,g,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                         MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE cross

!*****************************************************************
      SUBROUTINE hdcheck(a,b,c,d,e,f,t,dt,hel)
!-----------------------------------------------------------------
!
!  Computes the kinetic energy and helicity of the fluid.
!
! Output files contain:
! 'balance.txt':  time, <v^2>, <omega^2>, mechanic injection rate
! 'helicity.txt': time, kinetic helicity
!
! Parameters
!     a  : velocity field in the x-direction
!     b  : velocity field in the y-direction
!     c  : velocity field in the z-direction
!     d  : force in the x-direction
!     e  : force in the y-direction
!     f  : force in the z-direction
!     t  : number of time steps made
!     dt : time step
!     hel: =0 skips kinetic helicity computation
!          =1 computes the kinetic helicity
!
      USE fprecision
      USE commtypes
      USE grid
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: d,e,f
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)          :: c1,c2,c3
      DOUBLE PRECISION    :: eng,ens,pot,khe
      REAL(KIND=GP)       :: dt
      INTEGER, INTENT(IN) :: hel
      INTEGER, INTENT(IN) :: t
      INTEGER             :: i,j,k

!
! Computes the mean energy, enstrophy, and kinetic helicity
!
      CALL energy(a,b,c,eng,1)
      CALL energy(a,b,c,ens,0)
      IF (hel.eq.1) THEN
         CALL helicity(a,b,c,khe)
      ENDIF
!
! Computes the energy injection rate
!
      CALL cross(a,b,c,d,e,f,pot,1)
!
! Creates external files to store the results
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='balance.txt',position='append')
         WRITE(1,10) (t-1)*dt,eng,ens,pot
   10    FORMAT( 1P E13.6,2E23.16,E24.16 )
         CLOSE(1)
         IF (hel.eq.1) THEN
            OPEN(1,file='helicity.txt',position='append')
            WRITE(1,FMT='(1P E13.6,E24.16)') (t-1)*dt,khe
            CLOSE(1)
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE hdcheck

!*****************************************************************
      SUBROUTINE maxabs(a,b,c,d,kin)
!-----------------------------------------------------------------
!
! Computes the maximum absolute value of the field 
! vorticity, current density, or of the original field. 
! The output is only valid in the first node.
!
! Parameters
!     a  : field x-component
!     b  : field y-component
!     c  : field z-component
!     d  : at the output contains the maximum value
!     kin: =0 computes the maximum of vorticity
!          =1 computes the maximum of current density
!          =2 computes the maximum of the field
!
      USE fprecision
      USE commtypes
      USE fft
      USE grid
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)          :: c1,c2,c3
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)             :: r1,r2,r3
      REAL(KIND=GP), INTENT(OUT)   :: d
      REAL(KIND=GP)                :: dloc
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j,k

      IF (kin.eq.0) THEN
         CALL curlk(b,c,c1,1)
         CALL curlk(a,c,c2,2)
         CALL curlk(a,b,c3,3)
      ELSE IF (kin.eq.1) THEN
         CALL laplak(a,c1)
         CALL laplak(b,c2)
         CALL laplak(c,c3)
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  c1(k,j,i) = a(k,j,i)
                  c2(k,j,i) = b(k,j,i)
                  c3(k,j,i) = c(k,j,i)
               END DO
            END DO
         END DO
      ENDIF
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c3,r3,MPI_COMM_WORLD)
      dloc = 0.0_GP
!$omp parallel do if (pkend-ksta.ge.nth) private (j,i) reduction(max:dloc)
      DO k = ksta,pkend
!$omp parallel do if (pkend-ksta.lt.nth) private (i) reduction(max:dloc)
         DO j = 1,ny
            DO i = 1,nx
               dloc = max(dloc, &
                          sqrt(r1(i,j,k)**2+r2(i,j,k)**2+r3(i,j,k)**2))
            END DO
         END DO
      END DO
      dloc = dloc/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
      CALL MPI_REDUCE(dloc,d,1,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE maxabs


!*****************************************************************
      SUBROUTINE fc_filter(a)
!-----------------------------------------------------------------
!
! Filters a field in Fourier space using a cube with exponential 
! tails using the formula a*exp(-alpha(2 k_i/ N_i)**2p).
!
! Parameters
!     a    : input/output array

      USE fprecision
      USE kes
      USE mpivars
      USE grid
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: a
      DOUBLE PRECISION, PARAMETER  :: alpha=16*log(10d0), p=50d0
      INTEGER :: i,j,k

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
      DO j = 1,ny
      DO k = 1,nz
         a(k,j,i) = a(k,j,i)*exp(-alpha*(2*kx(i)/nx/Dkx)**(2*p))* &
                    exp(-alpha*(2*ky(j)/ny/Dky)**(2*p))* &
                    exp(-alpha*(2*kz(k)/nz/Dkz)**(2*p))
      END DO
      END DO
      END DO

      RETURN
      END SUBROUTINE fc_filter

!*****************************************************************
      SUBROUTINE divergence(a,b,c,d)
!-----------------------------------------------------------------
!
! Computes the mean square divergence of a vector field.
! The output is only valid in the first node.
!
! Parameters
!     a: input matrix in the x-direction
!     b: input matrix in the y-direction
!     c: input matrix in the z-direction
!     d: at the output contains the divergence
!
      USE var
      USE grid
      USE fcgram
      USE mpivars
      USE commtypes
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      DOUBLE PRECISION, INTENT(OUT) :: d

      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: C1,C2
      REAL(KIND=GP), DIMENSION(nz,ny,ista:iend)    :: R1

      DOUBLE PRECISION  :: dloc
      REAL(KIND=GP)     :: tmp

      INTEGER :: i,j,k

      ! One nx*ny from parsevals, an nx*ny*(nz-Cz) from averaging over
      ! physical domain and a nz**2 due to fft normalization.
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2/ &
            real(nz-Cz,kind=GP)
      d = 0d0
      dloc = 0d0

      CALL derivk(a,C1,1)
      C2 = C1

      CALL derivk(b,C1,2)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
      DO j = 1,ny
      DO k = 1,nz
         C2(k,j,i) = C2(k,j,i) + C1(k,j,i)
      ENDDO
      ENDDO
      ENDDO

      CALL derivk(c,C1,3)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
      DO j = 1,ny
      DO k = 1,nz
         C2(k,j,i) = C2(k,j,i) + C1(k,j,i)
      ENDDO
      ENDDO
      ENDDO

      CALL fftp1d_complex_to_real_z(plancr,C2,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
      DO j = 1,ny
      DO k = 1,nz
         R1(k,j,i) = real(c2(k,j,i), kind=GP)**2 + aimag(c2(k,j,i))**2
      ENDDO
      ENDDO
      ENDDO

      IF (ista.eq.1) THEN
!$omp parallel do private (k) reduction(+:dloc)
         DO j = 1,ny
            DO k = 1,nz-Cz
               dloc = dloc + R1(k,j,1)*tmp
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k) reduction(+:dloc)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k) reduction(+:dloc)
            DO j = 1,ny
               DO k = 1,nz-Cz
                  dloc = dloc + 2*R1(k,j,i)*tmp
               END DO
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k) reduction(+:dloc)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k) reduction(+:dloc)
           DO j = 1,ny
              DO k = 1,nz-Cz
                 dloc = dloc + 2*R1(k,j,i)*tmp
              END DO
           END DO
        END DO
      ENDIF

      CALL MPI_REDUCE(dloc,d,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
          MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE divergence

!*****************************************************************
      SUBROUTINE bouncheck(a,b,c,d,e,f,g)
!-----------------------------------------------------------------
!
! Computes the mean squared values of the tangential and normal
! components of a vector field at the boundaries.
! The output is only valid in the first node.
!
! Parameters
!     a: input matrix in the x-direction
!     b: input matrix in the y-direction
!     c: input matrix in the z-direction
!     d: at the output contains the mean squared value of the
!          tangential component at z=0.
!     d: at the output contains the mean squared value of the
!          tangential component at z=Lz.
!     d: at the output contains the mean squared value of the
!          normal component at z=0.
!     d: at the output contains the mean squared value of the
!          normal component at z=Lz.
!
      USE var
      USE grid
      USE fcgram
      USE mpivars
      USE commtypes
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      DOUBLE PRECISION, INTENT(OUT) :: d,e,f,g

      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)  :: C1
      REAL(KIND=GP), DIMENSION(ny,ista:iend)        :: R1,R2,R3,R4

      DOUBLE PRECISION   :: dloc,eloc,floc,gloc
      REAL(KIND=GP)      :: tmp

      INTEGER  :: i,j,k

      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2

      d = 0d00; e = 0d0; f = 0d0; g = 0d0
      dloc = 0d0; eloc = 0d0; floc = 0d0; gloc = 0d0

      ! Slip
      C1 = a
      CALL fftp1d_complex_to_real_z(plancr,C1,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth)
      DO j = 1,ny
         R1(j,i) = real(C1(1,j,i), kind=GP)**2 + aimag(c1(1,j,i))**2
         R2(j,i) = real(C1(nz-Cz,j,i), kind=GP)**2 +&
                     aimag(c1(nz-Cz,j,i))**2
      ENDDO
      ENDDO

      C1 = b
      CALL fftp1d_complex_to_real_z(plancr,C1,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth)
      DO j = 1,ny
         R1(j,i) = R1(j,i) + real(C1(1,j,i), kind=GP)**2 +&
                      aimag(c1(1,j,i))**2
         R2(j,i) = R2(j,i) + real(C1(nz-Cz,j,i), kind=GP)**2 +&
                     aimag(c1(nz-Cz,j,i))**2
      ENDDO
      ENDDO

      ! Normal
      C1 = c
      CALL fftp1d_complex_to_real_z(plancr,C1,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth)
      DO j = 1,ny
         R3(j,i) = real(C1(1,j,i), kind=GP)**2 + aimag(c1(1,j,i))**2
         R4(j,i) = real(C1(nz-Cz,j,i), kind=GP)**2 +&
                     aimag(c1(nz-Cz,j,i))**2
      ENDDO
      ENDDO

      IF (ista.eq.1) THEN
!$omp parallel do reduction(+:dloc,eloc,floc,gloc)
         DO j = 1,ny
               dloc = dloc + R1(j,1)*tmp
               eloc = eloc + R2(j,1)*tmp
               floc = floc + R3(j,1)*tmp
               gloc = gloc + R4(j,1)*tmp
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j) &
!$omp& reduction(+:dloc,eloc,floc,gloc)
         DO i = 2,iend
            DO j = 1,ny
                  dloc = dloc + 2*R1(j,i)*tmp
                  eloc = eloc + 2*R2(j,i)*tmp
                  floc = floc + 2*R3(j,i)*tmp
                  gloc = gloc + 2*R4(j,i)*tmp
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j) &
!$omp& reduction(+:dloc,eloc,floc,gloc)
         DO i = ista,iend
           DO j = 1,ny
               dloc = dloc + 2*R1(j,i)*tmp
               eloc = eloc + 2*R2(j,i)*tmp
               floc = floc + 2*R3(j,i)*tmp
               gloc = gloc + 2*R4(j,i)*tmp
           END DO
        END DO
      ENDIF

      CALL MPI_REDUCE(dloc,d,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
          MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(eloc,e,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
          MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(floc,f,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
          MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(gloc,g,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
          MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE bouncheck

!*****************************************************************
      SUBROUTINE diagnostic(a,b,c,t,dt)
!-----------------------------------------------------------------
!
!  Computes the mean squared divergence in the bulk of the fluid
!  as well as the mean squared slip velocity v_s and normal
! velocity v_n  at the boundaries.
!
! Output files contain:
! 'diagnostic.txt': time, <|div(v)|^2>, <|v_s|^2>|z=0,
!                   <|v_s|^2>|z=Lz, <|v_n|^2>|z=0, <|v_n|^2>|z=Lz
!
! Parameters
!     a  : velocity field in the x-direction
!     b  : velocity field in the y-direction
!     c  : velocity field in the z-direction
!     t  : number of time steps made
!     dt : time step
!
      USE fprecision
      USE grid
      USE mpivars

      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      REAL(KIND=GP), INTENT(IN)        :: dt
      INTEGER, INTENT(IN)              :: t
      DOUBLE PRECISION                 :: tmp,tmq,tmr,tms,tm1

      CALL divergence(a,b,c,tm1)
      CALL bouncheck(a,b,c,tmp,tmq,tmr,tms)

      IF ( myrank .eq. 0) THEN
          OPEN(1,file='diagnostic.txt',position='append')
          WRITE(1,FMT='(1P 2E13.6, 2E25.17, 2E13.6)') &
              (t-1)*dt,tm1,tmp,tmq,tmr,tms
          CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE diagnostic

!*****************************************************************
      SUBROUTINE normvec(a,b,c,d,kin)
!-----------------------------------------------------------------
!
! Normalizes a vector field.
!
! Parameters
!     a  : input matrix in the x-direction
!     b  : input matrix in the y-direction
!     c  : input matrix in the z-direction
!     d  : input specifying the desired value
!     kin: =0 normalize to specified magnetic energy
!          =1 normalize to specified kinetic energy
!          =2 normalize to specified magnetic enstrophy
!

      USE fprecision
      USE grid
      USE commtypes
      USE mpivars
!$    USE threads

      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: a,b,c
      REAL(KIND=GP), INTENT(IN) :: d
      INTEGER, INTENT(IN)       :: kin

      DOUBLE PRECISION  :: tmp
      REAL(KIND=GP)     :: rmp
      INTEGER :: i,j,k

      CALL energy(a,b,c,tmp,kin)
      CALL MPI_BCAST(tmp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      rmp = sqrt(d/tmp)

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
      DO j = 1,ny
      DO k = 1,nz
         a(k,j,i) = a(k,j,i)*rmp
         b(k,j,i) = b(k,j,i)*rmp
         c(k,j,i) = c(k,j,i)*rmp
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE normvec
