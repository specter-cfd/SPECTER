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

!******************************************************************
!  TODO: benchmark populating with zeros after nz-Cz and summing up
! to nz to see if it helps vectorization. (Energy and related
! subroutines)
!******************************************************************

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
      CALL fftp3d_complex_to_real(planfc,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(planfc,c2,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(planfc,c3,r3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(planfc,c4,r4,MPI_COMM_WORLD)

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
      CALL fftp3d_complex_to_real(planfc,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(planfc,c2,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(planfc,c3,r3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(planfc,c4,r4,MPI_COMM_WORLD)

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
      CALL fftp3d_complex_to_real(planfc,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(planfc,c2,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(planfc,c3,r3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(planfc,c4,r4,MPI_COMM_WORLD)

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

      CALL fftp3d_real_to_complex(planfc,rx,d,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planfc,ry,e,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planfc,rz,f,MPI_COMM_WORLD)

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
      CALL fftp3d_complex_to_real(planfc,d,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(planfc,e,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(planfc,f,r3,MPI_COMM_WORLD)
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
      CALL fftp3d_complex_to_real(planfc,d,r4,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(planfc,e,r5,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(planfc,f,r6,MPI_COMM_WORLD)
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

      CALL fftp3d_real_to_complex(planfc,r7,d,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planfc,r3,e,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planfc,r1,f,MPI_COMM_WORLD)

      RETURN
      END SUBROUTINE prodre

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
!     kin: =0 computes the energy of the curl of the field
!          =1 computes the energy of the field
!          =2 computes the energy of the curl^2 of the field
!
      USE fprecision
      USE commtypes
      USE grid
      USE fft
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      INTEGER, INTENT(IN)           :: kin
      DOUBLE PRECISION, INTENT(OUT) :: d

      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: C1,C2,C3,C4
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
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  C1(k,j,i) = a(k,j,i)
               END DO
            END DO
         END DO
         CALL fftp1d_complex_to_real_z(planfc,C1,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  R1(k,j,i) = real(C1(k,j,i), kind=GP)**2+aimag(C1(k,j,i))**2
               ENDDO
            ENDDO
         ENDDO

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  C1(k,j,i) = b(k,j,i)
               END DO
            END DO
         END DO
         CALL fftp1d_complex_to_real_z(planfc,C1,MPI_COMM_WORLD)
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

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  C1(k,j,i) = c(k,j,i)
               END DO
            END DO
         END DO
         CALL fftp1d_complex_to_real_z(planfc,C1,MPI_COMM_WORLD)
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


      ! Compute the energy of the curl of the field (a,b,c)
      ELSE IF (kin .eq. 0) THEN
         CALL curlk(b,c,C1,1)
         CALL fftp1d_complex_to_real_z(planfc,C1,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  R1(k,j,i) = real(C1(k,j,i), kind=GP)**2+aimag(C1(k,j,i))**2
               ENDDO
            ENDDO
         ENDDO

         CALL curlk(a,c,C1,2)
         CALL fftp1d_complex_to_real_z(planfc,C1,MPI_COMM_WORLD)
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

         CALL curlk(a,c,C1,2)
         CALL fftp1d_complex_to_real_z(planfc,C1,MPI_COMM_WORLD)
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


      ! Compute the energy of the curl^2 of the field (a,b,c)
      ELSE IF (kin .eq. 2) THEN
         CALL curlk(b,c,C1,1)
         CALL curlk(a,c,C2,2)
         CALL curlk(a,b,C3,3)
         
         CALL curlk(C2,C3,C4,1)
         CALL curlk(C1,C3,C3,2)
         CALL curlk(C1,C2,C1,3)

         CALL fftp1d_complex_to_real_z(planfc,C4,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  R1(k,j,i) = real(C4(k,j,i), kind=GP)**2+aimag(C4(k,j,i))**2
               ENDDO
            ENDDO
         ENDDO

         CALL fftp1d_complex_to_real_z(planfc,C3,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  R1(k,j,i) = R1(k,j,i) + real(C3(k,j,i), kind=GP)**2 +&
                                 aimag(C3(k,j,i))**2
               ENDDO
            ENDDO
         ENDDO

         CALL fftp1d_complex_to_real_z(planfc,C1,MPI_COMM_WORLD)
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
      ENDIF


      ! Compute mean value
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
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      DOUBLE PRECISION, INTENT(OUT) :: d

      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)             :: c1,c2
      REAL(KIND=GP), DIMENSION(nz,ny,ista:iend)                :: r1

      DOUBLE PRECISION :: dloc
      REAL(KIND=GP)    :: tmp
      INTEGER          :: i,j,k

      dloc = 0.0D0
      ! One nx*ny from parsevals, an nx*ny*(nz-Cz) from averaging over
      ! physical domain and a nz**2 due to fft normalization.
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2 / &
            real(nz-Cz,kind=GP)

      ! Compute the pointwise helicity
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               c1(k,j,i) = a(k,j,i)
            END DO
         END DO
      END DO
      CALL curlk(b,c,c2,1)
      CALL fftp1d_complex_to_real_z(planfc,c1,MPI_COMM_WORLD)
      CALL fftp1d_complex_to_real_z(planfc,c2,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               r1(k,j,i) = real(c1(k,j,i)*conjg(c2(k,j,i)))
            END DO
         END DO
      END DO

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               c1(k,j,i) = b(k,j,i)
            END DO
         END DO
      END DO
      CALL curlk(a,c,c2,2)
      CALL fftp1d_complex_to_real_z(planfc,c1,MPI_COMM_WORLD)
      CALL fftp1d_complex_to_real_z(planfc,c2,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               r1(k,j,i) = r1(k,j,i) + real(c1(k,j,i)*conjg(c2(k,j,i)))
            END DO
         END DO
      END DO

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               c1(k,j,i) = c(k,j,i)
            END DO
         END DO
      END DO
      CALL curlk(a,b,c2,3)
      CALL fftp1d_complex_to_real_z(planfc,c1,MPI_COMM_WORLD)
      CALL fftp1d_complex_to_real_z(planfc,c2,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               r1(k,j,i) = r1(k,j,i) + real(c1(k,j,i)*conjg(c2(k,j,i)))
            END DO
         END DO
      END DO

      ! Compute average values
      IF (ista.eq.1) THEN
!$omp parallel do private (k) reduction(+:dloc)
         DO j = 1,ny
            DO k = 1,nz-Cz
                 dloc = dloc+r1(k,j,1)*tmp
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k) reduction(+:dloc)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k) reduction(+:dloc)
            DO j = 1,ny
               DO k = 1,nz-Cz
                  dloc = dloc+2*r1(k,j,i)*tmp
               END DO
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k) reduction(+:dloc)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k) reduction(+:dloc)
            DO j = 1,ny
               DO k = 1,nz-Cz
                  dloc = dloc+2*r1(k,j,i)*tmp
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

      ! Computes the pointwise inner product between the fields
      IF (kin.eq.1) THEN
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  c1(k,j,i) = a(k,j,i)
                  c2(k,j,i) = d(k,j,i)
               END DO
            END DO
         END DO
      ELSE
         CALL curlk(b,c,c1,1)
         CALL curlk(e,f,c2,1)
      ENDIF 
      CALL fftp1d_complex_to_real_z(planfc,c1,MPI_COMM_WORLD)
      CALL fftp1d_complex_to_real_z(planfc,c2,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               r1(k,j,i) = real(c1(k,j,i)*conjg(c2(k,j,i)))
            END DO
         END DO
      END DO

      IF (kin.eq.1) THEN
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  c1(k,j,i) = b(k,j,i)
                  c2(k,j,i) = e(k,j,i)
               END DO
            END DO
         END DO
      ELSE
         CALL curlk(a,c,c1,2)
         CALL curlk(d,f,c2,2)
      ENDIF
      CALL fftp1d_complex_to_real_z(planfc,c1,MPI_COMM_WORLD)
      CALL fftp1d_complex_to_real_z(planfc,c2,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               r1(k,j,i) = r1(k,j,i) + real(c1(k,j,i)*conjg(c2(k,j,i)))
            END DO
         END DO
      END DO

      IF (kin.eq.1) THEN
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  c1(k,j,i) = c(k,j,i)
                  c2(k,j,i) = f(k,j,i)
               END DO
            END DO
         END DO
      ELSE
         CALL curlk(a,b,c1,3)
         CALL curlk(d,e,c2,3)
      ENDIF
      CALL fftp1d_complex_to_real_z(planfc,c1,MPI_COMM_WORLD)
      CALL fftp1d_complex_to_real_z(planfc,c2,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               r1(k,j,i) = r1(k,j,i) + real(c1(k,j,i)*conjg(c2(k,j,i)))
            END DO
         END DO
      END DO

      ! Compute average value
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

! Computes the mean energy, enstrophy, and kinetic helicity
      CALL energy(a,b,c,eng,1)
      CALL energy(a,b,c,ens,0)
      IF (hel.eq.1) THEN
         CALL helicity(a,b,c,khe)
      ENDIF

! Computes the energy injection rate
      CALL cross(a,b,c,d,e,f,pot,1)

! Creates external files to store the results
      IF (myrank.eq.0) THEN
         OPEN(1,file='balance.txt',position='append')
   10       FORMAT( 1P E13.6,2E23.16,E24.16 )
            WRITE(1,10) (t-1)*dt,eng,ens,pot
         CLOSE(1)
         IF (hel.eq.1) THEN
            OPEN(1,file='helicity.txt',position='append')
   11          FORMAT( 1P E13.6,E24.16 )
               WRITE(1,11) (t-1)*dt,khe
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
      CALL fftp3d_complex_to_real(planfc,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(planfc,c2,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(planfc,c3,r3,MPI_COMM_WORLD)
      dloc = 0.0_GP
!$omp parallel do if (pkend-ksta.ge.nth) private (j,i) reduction(max:dloc)
      DO k = ksta,pkend
!$omp parallel do if (pkend-ksta.lt.nth) private (i) reduction(max:dloc)
         DO j = 1,ny
            DO i = 1,nx
               dloc = max(dloc, r1(i,j,k)**2+r2(i,j,k)**2+r3(i,j,k)**2)
            END DO
         END DO
      END DO
      dloc = sqrt(dloc)/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
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
      USE fft
      USE mpivars
      USE commtypes
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

      ! Get pointwise divergence
      CALL derivk(a,C1,1)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               C2(k,j,i) = C1(k,j,i)
            ENDDO
         ENDDO
      ENDDO

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

      CALL fftp1d_complex_to_real_z(planfc,C2,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               R1(k,j,i) = real(c2(k,j,i),kind=GP)**2 + aimag(c2(k,j,i))**2
            ENDDO
         ENDDO
      ENDDO

      ! Compute mean value
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
