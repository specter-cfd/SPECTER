!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Extra subroutines to compute spatial derivatives and 
! nonlinear terms in MHD and equations in 3D using a 
! pseudo-spectral method. You should use the FFTPLANS and 
! MPIVARS modules (see the file 'fftp_mod.f90') in each 
! program that calls any of the subroutines in this file. 
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2003 Pablo D. Mininni.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar
!=================================================================

!*****************************************************************
      SUBROUTINE vector(a,b,c,d,e,f,c1,c2,c3)
!-----------------------------------------------------------------
!
! Computes the product AxB in real space. The components 
! of the vector fields A and B are given by the matrixes 
! a,b,c,d,e and f, following the right hand convention.
!
! Parameters
!     a: input matrix with A_x
!     b: input matrix with A_y
!     c: input matrix with A_z
!     d: input matrix with B_x
!     e: input matrix with B_y
!     f: input matrix with B_z
!    c1: at the output contains (AxB)_x in Fourier space
!    c2: at the output contains (AxB)_y in Fourier space
!    c3: at the output contains (AxB)_z in Fourier space
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: d,e,f
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: c1,c2,c3
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend) :: r1,r2
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend) :: r3,r4
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend) :: r5,r6
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend) :: r7
      REAL(KIND=GP)    :: tmp
      INTEGER :: i,j,k

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
      CALL fftp3d_complex_to_real(planfc,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(planfc,c2,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(planfc,c3,r3,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               c1(k,j,i) = d(k,j,i)
               c2(k,j,i) = e(k,j,i)
               c3(k,j,i) = f(k,j,i)
            END DO
         END DO
      END DO
      CALL fftp3d_complex_to_real(planfc,c1,r4,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(planfc,c2,r5,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(planfc,c3,r6,MPI_COMM_WORLD)

      tmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
!$omp parallel do if (iend-ista.ge.nth) private (j,i)
      DO k = ksta,pkend
!$omp parallel do if (iend-ista.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               r7(i,j,k) = (r2(i,j,k)*r6(i,j,k)-r5(i,j,k)*r3(i,j,k))*tmp
               r3(i,j,k) = (r3(i,j,k)*r4(i,j,k)-r6(i,j,k)*r1(i,j,k))*tmp
               r1(i,j,k) = (r1(i,j,k)*r5(i,j,k)-r4(i,j,k)*r2(i,j,k))*tmp
            END DO
         END DO
      END DO

      CALL fftp3d_real_to_complex(planfc,r7,c1,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planfc,r3,c2,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planfc,r1,c3,MPI_COMM_WORLD)

      RETURN
      END SUBROUTINE vector


!*****************************************************************
      SUBROUTINE mhdcheck(a,b,c,ma,mb,mc,t,dt,hel,crs)
!-----------------------------------------------------------------
!
! Consistency check for the conservation of the total
! energy, helicity, and null divergency of the fields
!
! Output files contain:
! 'balance.txt':  time, <v^2>+<b^2>, <omega^2>, <j^2>
! 'energy.txt':   time, <v^2>, <b^2>
! 'helicity.txt': time, kinetic helicity, magnetic helicity
! 'cross.txt':    time, <v.b>, <a^2>, generalized helicity [OPT]
!
! Parameters
!     a  : velocity field in the x-direction
!     b  : velocity field in the y-direction
!     c  : velocity field in the z-direction
!     ma : vector potential in the x-direction
!     mb : vector potential in the y-direction
!     mc : vector potential in the z-direction
!     t  : number of time steps made
!     dt : time step
!     hel: =0 skips helicity computation
!          =1 computes the helicity
!     crs: =0 skips cross helicity computation
!          =1 computes cross helicity
!
      USE fprecision
      USE commtypes
      USE grid
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: ma,mb,mc
      REAL(KIND=GP), INTENT(IN)    :: dt
      INTEGER, INTENT(IN) :: t
      INTEGER, INTENT(IN) :: hel,crs

      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: c1,c2,c3
      DOUBLE PRECISION    :: engk,engm,eng,ens
      DOUBLE PRECISION    :: divk,divm,asq,crh
      DOUBLE PRECISION    :: helk,helm,cur,tmp
      DOUBLE PRECISION    :: helg, g1sq, g2sq
      REAL(KIND=GP)                :: tmq
      INTEGER             :: i,j,k

      tmp = 0.0D0
      tmq = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2

!
! Computes the mean energy, enstrophy 
! and square current, and the kinetic 
! and magnetic helicity
!
      CALL energy(a,b,c,engk,1)
      CALL energy(a,b,c,ens,0)
      CALL energy(ma,mb,mc,engm,0)
      CALL energy(ma,mb,mc,cur,2)
      eng = engk+engm
      IF (hel.eq.1) THEN
         CALL helicity(a,b,c,helk)
         CALL helicity(ma,mb,mc,helm)
      ENDIF
!
! Computes the square vector potential, the 
! cross helicity
!
      IF (crs .eq. 1) THEN
         CALL energy(ma,mb,mc,asq,1)

         CALL derivk(ma,c1,1)
         CALL derivk(mb,c2,2)
         CALL derivk(mc,c3,3)

         CALL cross(a,b,c,c1,c2,c3,crh,1)
      ENDIF
!
! Creates external files to store the results
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='balance.txt',position='append')
   10       FORMAT( 1P E13.6,3E23.16 )
            WRITE(1,10) (t-1)*dt,eng,ens,cur
         CLOSE(1)
         OPEN(1,file='energy.txt',position='append')
   11       FORMAT( 1P E13.6,2E23.16 )
            WRITE(1,11) (t-1)*dt,engk,engm
         CLOSE(1)
         IF (hel.eq.1) THEN
            OPEN(1,file='helicity.txt',position='append')
   12          FORMAT( 1P E13.6,2E24.16 )
               WRITE(1,12) (t-1)*dt,helk,helm
            CLOSE(1)
         ENDIF
         IF (crs.eq.1) THEN
            OPEN(1,file='cross.txt',position='append')
   13          FORMAT( 1P E13.6,E23.16,E24.16 )
               WRITE(1,13) (t-1)*dt,crh,asq
            CLOSE(1)
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE mhdcheck
