!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Extra subroutines to compute the passive/active scalar 
! spectrum, transfer function, and associated global quantities 
! in the HD, MHD, Hall-MHD, and Boussinesq equations when a 
! passive or active scalar is present. You should use the 
! FCPLAN and MPIVARS modules (see the files 'fftp_mod.f90') in 
! each program that calls any of the subroutines in this file. 
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2009 Pablo D. Mininni.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar 
!=================================================================

!*****************************************************************
      SUBROUTINE advect(a,b,c,d,e)
!-----------------------------------------------------------------
!
! Three-dimensional inner product A.grad(B) in 
! real space. The components of the field A are 
! given by the arrays a, b and c, B is a scalar 
! quantity given by d.
!
! Parameters
!     a: input matrix in the x-direction
!     b: input matrix in the y-direction
!     c: input matrix in the z-direction
!     d: input matrix with the scalar
!     e: product (A.grad)B in Fourier space [output]
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a,b
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: c,d
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: e
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: c1,c2
      REAL(KIND=GP),    DIMENSION(nx,ny,ksta:kend) :: r1,r2
      REAL(KIND=GP),    DIMENSION(nx,ny,ksta:kend) :: r3
      REAL(KIND=GP)    :: tmp
      INTEGER :: i,j,k

!
! Computes (A_x.dx)B
!
      c1 = a
      CALL derivk(d,c2,1)
      CALL fftp3d_complex_to_real(planfc,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(planfc,c2,r2,MPI_COMM_WORLD)

!$omp parallel do if (pkend-ksta.ge.nth) private (j,i)
      DO k = ksta,pkend
!$omp parallel do if (pkend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               r3(i,j,k) = r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO
!
! Computes (A_y.dy)B
!
      c1 = b
      CALL derivk(d,c2,2)
      CALL fftp3d_complex_to_real(planfc,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(planfc,c2,r2,MPI_COMM_WORLD)

!$omp parallel do if (pkend-ksta.ge.nth) private (j,i)
      DO k = ksta,pkend
!$omp parallel do if (pkend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               r3(i,j,k) = r3(i,j,k)+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO
!
! Computes (A_z.dz)B
!
      c1 = c
      CALL derivk(d,c2,3)
      CALL fftp3d_complex_to_real(planfc,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(planfc,c2,r2,MPI_COMM_WORLD)

! Add and normalize at the same time
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
!$omp parallel do if (pkend-ksta.ge.nth) private (j,i)
      DO k = ksta,pkend
!$omp parallel do if (pkend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               r3(i,j,k) = (r3(i,j,k)+r1(i,j,k)*r2(i,j,k))*tmp
            END DO
         END DO
      END DO

      CALL fftp3d_real_to_complex(planfc,r3,e,MPI_COMM_WORLD)

      RETURN
      END SUBROUTINE advect

!*****************************************************************
      SUBROUTINE variance(a,b,kin)
!-----------------------------------------------------------------
!
! Computes the mean variance of the scalar.
! The output is only valid in the first node.
!
! Parameters
!     a  : input matrix with the scalar
!     b  : at the output contains the variance
!     kin: =0 computes the variance of k^2 times the scalar
!          =1 computes the variance of the scalar
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a
      DOUBLE PRECISION, INTENT(OUT) :: b

      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: at
      DOUBLE PRECISION              :: bloc
      REAL(KIND=GP)                 :: tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j,k

      bloc = 0.0D0
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2/ &
            real(nz-Cz,kind=GP)

      IF (kin .eq. 1) THEN
         at = a
      ELSE
         at = kk2*a
      ENDIF
      
      CALL fftp1d_complex_to_real_z(planfc,at,MPI_COMM_WORLD)

!
! Computes the variance
!
      IF (ista.eq.1) THEN
!$omp parallel do private (k) reduction(+:bloc)
         DO j = 1,ny
            DO k = 1,nz-Cz
               bloc = bloc+tmp*abs(at(k,j,1))**2
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k) reduction(+:bloc)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k) reduction(+:bloc)
            DO j = 1,ny
               DO k = 1,nz-Cz
                  bloc = bloc+2*tmp*abs(at(k,j,i))**2
               END DO
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k) reduction(+:bloc)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k) reduction(+:bloc)
            DO j = 1,ny
               DO k = 1,nz-Cz
                  bloc = bloc+2*tmp*abs(at(k,j,i))**2
               END DO
            END DO
         END DO
      ENDIF
!
! Computes the reduction between nodes
!
      CALL MPI_REDUCE(bloc,b,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE variance

!*****************************************************************
      SUBROUTINE product(a,b,c)
!-----------------------------------------------------------------
!
! Computes the integral of the product of two scalars. 
! The output is only valid in the first node.
!
! Parameters
!     a  : first scalar
!     b  : second scalar
!     c  : at the output contains the product
!
      USE fprecision
      USE commtypes
      USE grid
      USE mpivars
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      DOUBLE PRECISION, INTENT(OUT) :: c

      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: at,bt
      DOUBLE PRECISION              :: cloc
      REAL(KIND=GP)                 :: tmp
      INTEGER             :: i,j,k

      cloc = 0.0D0
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2/ &
            real((nz-Cz),kind=GP)

      at = a; bt = b
      CALL fftp1d_complex_to_real_z(planfc,at,MPI_COMM_WORLD)
      CALL fftp1d_complex_to_real_z(planfc,bt,MPI_COMM_WORLD)

!
! Computes the averaged inner product between the fields
!
      IF (ista.eq.1) THEN
!$omp parallel do private (k) reduction(+:cloc)
         DO j = 1,ny
            DO k = 1,nz-Cz
               cloc = cloc+tmp*real(at(k,j,1)*conjg(bt(k,j,1)))
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k) reduction(+:cloc)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k) reduction(+:cloc)
            DO j = 1,ny
               DO k = 1,nz-Cz
                  cloc = cloc+2*tmp*real(at(k,j,i)*conjg(bt(k,j,i)))
               END DO
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k) reduction(+:cloc)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k) reduction(+:cloc)
            DO j = 1,ny
               DO k = 1,nz-Cz
                  cloc = cloc+2*tmp*real(at(k,j,i)*conjg(bt(k,j,i)))
               END DO
            END DO
         END DO
      ENDIF
!
! Computes the reduction between nodes
!
      CALL MPI_REDUCE(cloc,c,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                         MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE product

!*****************************************************************
      SUBROUTINE pscheck(a,b,t,dt)
!-----------------------------------------------------------------
!
! Consistency check for the conservation of energy, 
! helicity, and null divergency of the velocity field
!
! Output file contains:
! 'scalar.txt':  time, <theta^2>, <|grad(theta)|^2>, injection rate
!
! Parameters
!     a : scalar concentration
!     b : source of the scalar
!     t : number of time steps made
!     dt: time step
!
      USE fprecision
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      DOUBLE PRECISION    :: eng,ens,pot
      REAL(KIND=GP), INTENT(IN)    :: dt
      INTEGER, INTENT(IN) :: t
      INTEGER             :: i,j,k

!
! Computes the variance and the variance of k^2 times the scalar
!
      CALL variance(a,eng,1)
      CALL variance(a,ens,0)
!
! Computes the scalar injection rate
!
      CALL product(a,b,pot)
!
! Creates external files to store the results
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='scalar.txt',position='append')
         WRITE(1,10) (t-1)*dt,eng,ens,pot
   10    FORMAT( 1P E13.6,2E22.14,E23.14 )
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE pscheck

!*****************************************************************
      SUBROUTINE normsca(a,b,kin)
!-----------------------------------------------------------------
!
! Normalizes a scalar field.
!
! Parameters
!     a  : input matrix in the x-direction
!     b  : input specifying the desired value
!     kin: =0 normalizes to specified variance of k^2 times the
!             scalar
!          =1 normalizes to specified variance of the scalar
!

      USE fprecision
      USE grid
      USE commtypes
      USE mpivars
!$    USE threads

      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: a
      REAL(KIND=GP), INTENT(IN) :: b
      INTEGER, INTENT(IN)       :: kin

      DOUBLE PRECISION  :: tmp
      REAL(KIND=GP)     :: rmp
      INTEGER :: i,j,k

      CALL variance(a,tmp,kin)
      CALL MPI_BCAST(tmp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      rmp = sqrt(b/tmp)

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
      DO j = 1,ny
      DO k = 1,nz
         a(k,j,i) = a(k,j,i)*rmp
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE normsca
