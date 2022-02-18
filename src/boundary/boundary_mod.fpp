!======================================================================
! SPECTER boundary module
!
!
! 2020 Mauro Fontana. DF-UBA.
!======================================================================
! Definitions for conditional compilation
#include "../specter.h"

MODULE boundary
      USE fprecision

      IMPLICIT NONE

      TYPE BCPLAN
         INTEGER  :: bcxsta,bcxend,bcysta,bcyend,bczsta,bczend
      ENDTYPE BCPLAN

      CONTAINS

      INCLUDE 'vboundary.f90'
#ifdef SCALAR_
      INCLUDE 'sboundary.f90'
#endif
#ifdef MAGFIELD_
      INCLUDE 'bboundary.f90'
#endif

!***************************************************************************
      SUBROUTINE setup_bc(this,planfc,bckind,field,vplan)
!---------------------------------------------------------------------------
!  Add field and boundary conditions to the plan
!  ARGUMENTS:
!    this   : 'this' instance of a BCPLAN. [INOUT]
!    fcplan : FCPLAN that will be used for the boundary conditions. [INOUT]
!    bckind : Array of strings with the kind of boundary condition
!             for each boundary. [INOUT]
!    field  : Field whose boundary conditions must be set up. Options
!             include 'v', 'b', 'th'.
!---------------------------------------------------------------------------
      USE fcgram

      IMPLICIT NONE

      TYPE(BCPLAN), INTENT(INOUT)          :: this
      TYPE(FCPLAN), INTENT(INOUT)          :: planfc
      CHARACTER(len=15), INTENT(IN)        :: bckind(6)
      CHARACTER(len=*), INTENT(IN)         :: field

      TYPE(BCPLAN), INTENT(IN), OPTIONAL   :: vplan

      IF ( field .eq. 'v' ) THEN
         CALL v_setup(this,planfc,bckind)
#ifdef SCALAR_
      ELSEIF ( field .eq. 's' ) THEN
         CALL s_setup(this,planfc,bckind)
#endif
#ifdef MAGFIELD_
      ELSEIF ( field .eq. 'b' ) THEN
         IF (PRESENT(vplan)) THEN
            CALL b_setup(this,planfc,bckind,vplan=vplan)
         ELSE
            CALL b_setup(this,planfc,bckind)
         ENDIF
#endif
      ENDIF

      RETURN
      ENDSUBROUTINE setup_bc

!**********************************************************************
      SUBROUTINE goto_domain_w_boundaries(planbc,planfc,a,b,c)
!----------------------------------------------------------------------
!  Transform given arrays from (kz,ky,kx) domain to a representation
! that encompasses all non-periodic boundaries. The interior points of
! the array normalized by the appropriate FFT factor.
!  ARGUMENTS:
!    planbc  : Plan for the BC of the target field. [IN]
!    planfc  : FCPLAN to use. [IN]
!         a  : target array. [IN] 
!         b  : target array. (OPTIONAL) [IN] 
!         c  : target array. (OPTIONAL) [IN] 
!----------------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE commtypes
      USE fcgram
!$    USE threads
      IMPLICIT NONE

      TYPE(BCPLAN), INTENT(IN)  :: planbc
      TYPE(FCPLAN), INTENT(IN)  :: planfc
      
      COMPLEX(kind=GP),INTENT(INOUT),DIMENSION(planfc%z%n,planfc%y%n,ista:iend)          :: a
      COMPLEX(kind=GP),OPTIONAL,INTENT(INOUT),DIMENSION(planfc%z%n,planfc%y%n,ista:iend) :: b,c

      REAL(KIND=GP)  :: tmp
      INTEGER        :: i,j,k

      ! Note: borders don't require normalization as they will be overwritten
      ! however normalizing them anyway helps vectorization.

      IF ( planbc%bcysta .gt. -1 .AND. planbc%bcyend .gt. -1 ) THEN
         !TODO transform to (z,y,kx)
         STOP
      ELSEIF ( planbc%bczsta .gt. -1 .AND. planbc%bczend .gt. -1 ) THEN
         tmp = 1.0_GP/real(planfc%z%n, kind=GP)
         CALL fftp1d_complex_to_real_z(planfc,a,MPI_COMM_WORLD)
         ! Normalize
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,planfc%y%n
               DO k = 1,planfc%z%n-planfc%z%C
                  a(k,j,i) = a(k,j,i)*tmp
               ENDDO
            ENDDO
         ENDDO

         IF ( PRESENT(b) ) THEN
            CALL fftp1d_complex_to_real_z(planfc,b,MPI_COMM_WORLD)
            ! Normalize
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
               DO j = 1,planfc%y%n
                  DO k = 1,planfc%z%n-planfc%z%C
                     b(k,j,i) = b(k,j,i)*tmp
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

         IF ( PRESENT(c) ) THEN
            CALL fftp1d_complex_to_real_z(planfc,c,MPI_COMM_WORLD)
            ! Normalize
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
               DO j = 1,planfc%y%n
                  DO k = 1,planfc%z%n-planfc%z%C
                     c(k,j,i) = c(k,j,i)*tmp
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE goto_domain_w_boundaries 

!**********************************************************************
      SUBROUTINE goto_3d_fourier(planbc,planfc,a,b,c)
!----------------------------------------------------------------------
!  Transform given arrays from the minimum representation that
! encompasses all non-periodic boundaries to the 3D wavenumber domain.
!  ARGUMENTS:
!    planbc  : Plan for the BC of the target field. [IN]
!    planfc  : FCPLAN to use. [IN]
!         a  : target array. [IN] 
!         b  : target array. (OPTIONAL) [IN] 
!         c  : target array. (OPTIONAL) [IN] 
!----------------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE commtypes
      USE fcgram
!$    USE threads
      IMPLICIT NONE

      TYPE(BCPLAN), INTENT(IN)  :: planbc
      TYPE(FCPLAN), INTENT(IN)  :: planfc
      
      COMPLEX(kind=GP),INTENT(INOUT),DIMENSION(planfc%z%n,planfc%y%n,ista:iend)          :: a
      COMPLEX(kind=GP),OPTIONAL,INTENT(INOUT),DIMENSION(planfc%z%n,planfc%y%n,ista:iend) :: b,c


      IF ( planbc%bcysta .gt. -1 .AND. planbc%bcyend .gt. -1 ) THEN
         !TODO transform (z,y,kx) to (kz,ky,kx)
         STOP   
      ELSEIF ( planbc%bczsta .gt. -1 .AND. planbc%bczend .gt. -1 ) THEN
         CALL fftp1d_real_to_complex_z(planfc,a,MPI_COMM_WORLD)

         IF ( PRESENT(b) ) THEN
            CALL fftp1d_real_to_complex_z(planfc,b,MPI_COMM_WORLD)
         ENDIF

         IF ( PRESENT(c) ) THEN
            CALL fftp1d_real_to_complex_z(planfc,c,MPI_COMM_WORLD)
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE goto_3d_fourier 

!***********************************************************************
      SUBROUTINE sol_project(a,b,c,d,bctarget,bczsta,bczend)
!-----------------------------------------------------------------------
!    Projects a given vector field onto a solenoidal space. For this it 
!   solves the Poisson equation laplacian(d) = div(v) with a,b,c being
!   the x,y,z components of the vector field v. It allows for different
!   combinations of homogeneous boundary conditions for the potential 
!   and also Dirichlet boundary conditions for the normal component
!   of the projected field.
!   After the computation the potential solution is applied to project
!   v' = v - grad(d). The potential d is returned in the mixed (z,ky,kx)
!   domain.
!   NOTE: Robin boundary conditions assume the proportionality constant
!   to be sqrt(kx**2+ky**2).
!   !TODO update documentation about p' = dt*p for convenience.
!
! Parameters
!     a, b, c  : x,y,z components of v, respectively (kz,ky,kx) [INOUT]
!     d        : solution to the Poisson equation p  (z,ky,kx)    [OUT]
!     bctarget : 0 if BC should be applied to d                    [IN]
!                1 if BC should be satisfied by c (after projection).
!                Only Dirichlet supported for c
!     bczsta   : kind of boundary condition at z=0                 [IN]
!              * =0 for homogeneous Dirichlet BC
!              * =2 for homogeneous Robin BC
!     bczend   : same as bczsta but for z=Lz                       [IN]
!              * =0 for homogeneous Dirichlet BC
!              * =2 for homogeneous Robin BC
!
      USE grid
      USE var
      USE kes
      USE fft
      USE commtypes
      USE mpivars
!$    USE threads

      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend)   :: d
      INTEGER, INTENT(IN) :: bctarget,bczsta,bczend

      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)                :: C1,C2,C3
      COMPLEX(KIND=GP), DIMENSION(2,ny,ista:iend)                 :: bc

      REAL(KIND=GP)    :: tmp
      INTEGER          :: i,j,k

      ! Get inhomogeneous solution and apply it
      ! ---------------------------------------
      CALL poisson_inhomogeneous(a,b,c,d)

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

      ! Boundary condition for Laplace
      ! ------------------------------
      tmp = 1.0_GP/nz
      IF (bctarget .eq. 0) THEN
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i=ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j=1,ny
               DO k=1,nz
                  C1(k,j,i) = d(k,j,i)*tmp
               ENDDO
            ENDDO
         ENDDO
      ELSEIF (bctarget .eq. 1) THEN
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i=ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j=1,ny
               DO k=1,nz
                  C1(k,j,i) = c(k,j,i)*tmp
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      IF ( bczsta .eq. 2 .OR. bczend .eq. 2 ) THEN  !For Robin
          CALL derivk(C1,C2,3)
          CALL fftp1d_complex_to_real_z(planfc,C2,MPI_COMM_WORLD)
      ENDIF
      CALL fftp1d_complex_to_real_z(planfc,C1,MPI_COMM_WORLD)

      IF (bczsta .eq. 0) THEN
         IF (bctarget .eq. 0) THEN
!$omp parallel do private (j)
            DO i=ista,iend
               DO j=1,ny
                  bc(1,j,i) = -C1(1,j,i)
               ENDDO
            ENDDO
         ELSEIF (bctarget .eq. 1) THEN
!$omp parallel do private (j)
            DO i=ista,iend
               DO j=1,ny
                  bc(1,j,i) = C1(1,j,i)
               ENDDO
            ENDDO
         ENDIF
      ELSEIF (bczsta .eq. 2) THEN
!$omp parallel do private (j)
         DO i=ista,iend
            DO j=1,ny
               ! C2 positive because external solution for z<0 goes like
               ! exp(khom*z) instead of exp(-khom*z)
               bc(1,j,i) = C2(1,j,i) - khom(j,i)*C1(1,j,i)
            ENDDO
         ENDDO
      ELSE
          IF ( myrank .eq. 0) PRINT*, "[ERROR] Unsupported BC kind in call "&
              "to sol_project. Aborting..."
         CALL MPI_FINALIZE(MPI_COMM_WORLD, ierr)
         STOP
      ENDIF

      IF (bczend .eq. 0) THEN
         IF (bctarget .eq. 0) THEN
!$omp parallel do private (j)
            DO i=ista,iend
               DO j=1,ny
                  bc(2,j,i) = -C1(nz-Cz,j,i)
               ENDDO
            ENDDO
         ELSEIF (bctarget .eq. 1) THEN
!$omp parallel do private (j)
            DO i=ista,iend
               DO j=1,ny
                  bc(2,j,i) = C1(nz-Cz,j,i)
               ENDDO
            ENDDO
         ENDIF
      ELSEIF (bczend .eq. 2) THEN
!$omp parallel do private (j)
         DO i=ista,iend
            DO j=1,ny
               bc(2,j,i) = -(C2(nz-Cz,j,i) + khom(j,i)*C1(nz-Cz,j,i))
            ENDDO
         ENDDO
      ELSE
          IF ( myrank .eq. 0) PRINT*, "[ERROR] Unsupported BC kind in call "&
              "to sol_project. Aborting..."
         CALL MPI_FINALIZE(MPI_COMM_WORLD, ierr)
         STOP
      ENDIF

      ! Solve Laplace equation and construct d_hom(z,ky,kx)
      ! --------------------------------------------------
      ! Dirichlet in field means Neumann condition on potential
      CALL laplace_z(bc,C2,C3,bczsta+bctarget,bczend+bctarget)

      ! Get d = d_inh + d_hom in (z,ky,kx) domain
      IF (bctarget .eq. 0) THEN
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i=ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j=1,ny
               DO k=1,nz
                  d(k,j,i) = C1(k,j,i) + C2(k,j,i)
               ENDDO
            ENDDO
         ENDDO
      ELSEIF (bctarget .eq. 1) THEN
         CALL fftp1d_complex_to_real_z(planfc,d,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i=ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j=1,ny
               DO k=1,nz
                  d(k,j,i) = d(k,j,i)*tmp + C2(k,j,i)
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      ! Remove harmonic contribution from field
      ! ---------------------------------------
      CALL fftp1d_real_to_complex_z(planfc,C2,MPI_COMM_WORLD)
      CALL fftp1d_real_to_complex_z(planfc,C3,MPI_COMM_WORLD)

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
      END SUBROUTINE sol_project

!*****************************************************************
      SUBROUTINE poisson_inhomogeneous(a,b,c,d)
!-----------------------------------------------------------------
!   Solves the inhomogeneous part of the Poisson equation
!   laplacian(d') = div(v), with a,b,c the x,y,z components of v.
!   It takes the components of v as input matrices and returns 
!   d'_inhom in Fourier space. Note that the boundary
!   conditions must be enforced via the homogeneous solution a
!   a posteriori.
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
      END SUBROUTINE poisson_inhomogeneous

!******************************************************************
      SUBROUTINE laplace_z(bc,a,b,bczsta,bczend)
!------------------------------------------------------------------
!  Returns a in (z,ky,kx) coordinates, where a is a solution to the
! Laplace equation. It assumes periodic boundary conditions in x,y
! and supports Dirichlet, Neumann and Robin boundary conditions.
! For the Robin case the proportionality constant is assumed to be
! sqrt(kx**2+ky**2).
!
! Parameters
!     bc : (2,ky,kx)input matrix for the boundary conditions for
!           a, with          shape (2,ky,kx).                  [IN]
!     a  : at the output contains the result in the mixed
!          (z,ky,kx) space.                                   [OUT]
!     b  : at the output contains the normal derivative
!          of the result in the mixed (z,ky,kx) space.        [OUT]
!     bczsta  : kind of boundary condition at z=0              [IN]
!             * =0 for homogeneous Dirichlet BC
!             * =1 for homogeneous Neumann BC
!             * =2 for homogeneous Robin BC
!     bczend  : kind of boundary condition at z=Lz             [IN]
!             * =0 for homogeneous Dirichlet BC
!             * =1 for homogeneous Neumann BC
!             * =2 for homogeneous Robin BC
!

      USE kes
      USE var
      USE grid
      USE fft
      USE mpivars
      USE commtypes
!$    USE threads


      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(2,ny,ista:iend)    :: bc
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend)  :: a,b
      INTEGER, INTENT(IN) :: bczsta,bczend

      COMPLEX(KIND=GP), DIMENSION(2,ny,ista:iend)    :: coef

      REAL(KIND=GP)                :: tmp

      INTEGER                      :: i,j,k


      ! Get coefficients
      ! ----------------
      ! Pure Dirichlet
      IF ( bczsta .eq. bczend .AND. bczsta .eq. 0 ) THEN
      IF (ista.eq.1) THEN
            coef(1,1,1) = (bc(2,1,1)-bc(1,1,1))/Lz
            coef(2,1,1) = bc(1,1,1)
!$omp parallel do private (k,tmp)
         DO j = 2,ny
            tmp   = 1.0_GP/(1 - exp(-2*khom(j,1)*Lz))
            coef(1,j,1) = (bc(2,j,1) - bc(1,j,1)*exp(-khom(j,1)*Lz))*tmp
            coef(2,j,1) = (bc(1,j,1) - bc(2,j,1)*exp(-khom(j,1)*Lz))*tmp
         END DO
!$omp parallel do if (iend-ista.ge.nth) private (j,k,tmp)
         DO i = 2,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,tmp)
            DO j = 1,ny
               tmp   = 1.0_GP/(1 - exp(-2*khom(j,i)*Lz))
               coef(1,j,i) = (bc(2,j,i) - bc(1,j,i)*exp(-khom(j,i)*Lz))*tmp
               coef(2,j,i) = (bc(1,j,i) - bc(2,j,i)*exp(-khom(j,i)*Lz))*tmp
            END DO
         END DO
      ELSE 
!$omp parallel do if (iend-ista.ge.nth) private (j,k,tmp)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,tmp)
            DO j = 1,ny
                 tmp   = 1.0_GP/(1 - exp(-2*khom(j,i)*Lz))
                 coef(1,j,i) = (bc(2,j,i) - bc(1,j,i)*exp(-khom(j,i)*Lz))*tmp
                 coef(2,j,i) = (bc(1,j,i) - bc(2,j,i)*exp(-khom(j,i)*Lz))*tmp
            END DO
         END DO
      ENDIF  !End Dirichlet

      ! Pure Neumann
      ELSEIF ( bczsta .eq. bczend .AND. bczsta .eq. 1 ) THEN
      IF (ista.eq.1) THEN
         coef(1,1,1) = bc(1,1,1)
         coef(2,1,1) = 0.0_GP
!$omp parallel do private (k,tmp)
         DO j = 2,ny
            tmp   = 1.0_GP/(khom(j,1)*(1-exp(-2*khom(j,1)*Lz)))
            coef(1,j,1) = ( bc(2,j,1) - bc(1,j,1)*exp(-khom(j,1)*Lz))*tmp
            coef(2,j,1) = (-bc(1,j,1) + bc(2,j,1)*exp(-khom(j,1)*Lz))*tmp
         ENDDO
!$omp parallel do if (iend-ista.ge.nth) private (j,k,tmp)
         DO i = 2,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,tmp)
            DO j = 1,ny
               tmp   = 1.0_GP/(khom(j,i)*(1-exp(-2*khom(j,i)*Lz)))
               coef(1,j,i) = ( bc(2,j,i) - bc(1,j,i)*exp(-khom(j,i)*Lz) )*tmp
               coef(2,j,i) = (-bc(1,j,i) + bc(2,j,i)*exp(-khom(j,i)*Lz) )*tmp
            END DO
         END DO
      ELSE 
!$omp parallel do if (iend-ista.ge.nth) private (j,k,tmp)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,tmp)
            DO j = 1,ny
               tmp   = 1.0_GP/(khom(j,i)*(1-exp(-2*khom(j,i)*Lz)))
               coef(1,j,i) = ( bc(2,j,i) - bc(1,j,i)*exp(-khom(j,i)*Lz) )*tmp
               coef(2,j,i) = (-bc(1,j,i) + bc(2,j,i)*exp(-khom(j,i)*Lz) )*tmp
            END DO
         END DO
      ENDIF  ! End Neumann

      ! Pure Robin
      ELSEIF ( bczsta .eq. bczend .AND. bczsta .eq. 2 ) THEN
      IF (ista.eq.1) THEN
         coef(1,1,1) = bc(1,1,1)
         coef(2,1,1) = 0.0_GP
!$omp parallel do private (k,tmp)
         DO j = 2,ny
            tmp   = 1.0_GP/(2*khom(j,1))
            coef(1,j,1) = bc(2,j,1)*tmp
            coef(2,j,1) = bc(1,j,1)*tmp
         END DO
!$omp parallel do if (iend-ista.ge.nth) private (j,k,tmp)
         DO i = 2,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,tmp)
            DO j = 1,ny
               tmp   = 1.0_GP/(2*khom(j,i))
               coef(1,j,i) = bc(2,j,i)*tmp
               coef(2,j,i) = bc(1,j,i)*tmp
            END DO
         END DO
      ELSE 
!$omp parallel do if (iend-ista.ge.nth) private (j,k,tmp)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,tmp)
            DO j = 1,ny
               tmp   = 1.0_GP/(2*khom(j,i))
               coef(1,j,i) = bc(2,j,i)*tmp
               coef(2,j,i) = bc(1,j,i)*tmp
           END DO
         END DO
      ENDIF  !End Robin

      ! Dirichlet bottom - Robin top
      ELSEIF ( bczsta .eq. 0 .AND. bczend .eq. 2 ) THEN
      IF (ista.eq.1) THEN
         coef(1,1,1) = bc(2,1,1)
         coef(2,1,1) = bc(1,1,1)
!$omp parallel do private (k,tmp)
         DO j = 2,ny
            tmp   = 1.0_GP/(2*khom(j,1))
            coef(1,j,1) = bc(2,j,1)*tmp
            coef(2,j,1) = (bc(1,j,1)*2*khom(j,1) - bc(2,j,1)*exp(-khom(j,1)*Lz))*tmp
         END DO
!$omp parallel do if (iend-ista.ge.nth) private (j,k,tmp)
         DO i = 2,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,tmp)
            DO j = 1,ny
               tmp   = 1.0_GP/(2*khom(j,i))
               coef(1,j,i) = bc(2,j,i)*tmp
               coef(2,j,i) = (bc(1,j,i)*2*khom(j,i) - bc(2,j,i)*exp(-khom(j,i)*Lz))*tmp
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,tmp)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,tmp)
            DO j = 1,ny
               tmp   = 1.0_GP/(2*khom(j,i))
               coef(1,j,i) = bc(2,j,i)*tmp
               coef(2,j,i) = (bc(1,j,i)*2*khom(j,i) - bc(2,j,i)*exp(-khom(j,i)*Lz))*tmp
           END DO
         END DO
      ENDIF  !End Dirichlet bottom - Robin top

      ELSE  ! Unsupported combination of BC
         IF ( myrank .eq. 0) PRINT*, "[ERROR] Unsupported BC combination in call "&
              "to laplace_z. Aborting..."
         CALL MPI_FINALIZE(MPI_COMM_WORLD, ierr)
         STOP
      ENDIF  !End coefficients

      ! Construct solution and normal derivative
      ! ----------------------------------------
      IF (ista.eq.1) THEN
!$omp parallel do
         DO k=1,nz
            a(k,1,1) = real(coef(1,1,1), kind=GP)*z(k) + real(coef(2,1,1), kind=GP)
            b(k,1,1) = real(coef(1,1,1), kind=GP)
         ENDDO
!$omp parallel do private (k)
         DO j = 2,ny
            DO k=1,nz
               a(k,j,1) = coef(1,j,1)*exp(khom(j,1)*(z(k)-Lz)) + &
                          coef(2,j,1)*exp(-khom(j,1)*z(k))
               b(k,j,1) = khom(j,1)*(coef(1,j,1)*exp(khom(j,1)*(z(k)-Lz)) - &
                          coef(2,j,1)*exp(-khom(j,1)*z(k)))
            END DO
         END DO
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = 2,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k=1,nz
                  a(k,j,i) = coef(1,j,i)*exp(khom(j,i)*(z(k)-Lz)) + &
                             coef(2,j,i)*exp(-khom(j,i)*z(k))
                  b(k,j,i) = khom(j,i)*(coef(1,j,i)*exp(khom(j,i)*(z(k)-Lz)) - &
                             coef(2,j,i)*exp(-khom(j,i)*z(k)))
               END DO
            END DO
         END DO
      ELSE 
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k=1,nz
                  a(k,j,i) = coef(1,j,i)*exp(khom(j,i)*(z(k)-Lz)) + &
                             coef(2,j,i)*exp(-khom(j,i)*z(k))
                  b(k,j,i) = khom(j,i)*(coef(1,j,i)*exp(khom(j,i)*(z(k)-Lz)) - &
                             coef(2,j,i)*exp(-khom(j,i)*z(k)))
               END DO
            END DO
         END DO
      ENDIF

      RETURN
      END SUBROUTINE laplace_z

!**********************************************************************
      SUBROUTINE bouncheck_z(planfc,bot,top,a,b)
!----------------------------------------------------------------------
!
!  Computes the mean values of the a**2 (or a**2+b*2 if the 
! latter is present) at z=0 and z=Lz.
!  The output is only valid in the first node.
!
! Parameters
!     bot: at the output contains the mean squared value of the
!            at z=0. [OUT]
!     top: at the output contains the mean squared value of the
!            at z=Lz. [OUT]
!     a:   the array whose mean squared value is to be computed. [IN]
!     b:   if present, the mean value of a**2+b**2 is 
!            computed. [IN] (OPTIONAL)
!
      USE var
      USE grid
      USE fcgram
      USE mpivars
      USE commtypes
!$    USE threads
      IMPLICIT NONE

      TYPE(FCPLAN), INTENT(IN)  :: planfc

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend)           :: a
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend), OPTIONAL :: b

      DOUBLE PRECISION, INTENT(OUT) :: bot,top

      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)  :: C1
      REAL(KIND=GP), DIMENSION(ny,ista:iend)        :: R1,R2

      DOUBLE PRECISION   :: botloc,toploc
      REAL(KIND=GP)      :: tmp

      INTEGER  :: i,j,k

      tmp = (1.0_GP/(planfc%x%n*planfc%y%n*planfc%z%n))**2

      bot    = 0d0; top    = 0d0
      botloc = 0d0; toploc = 0d0

      ! Get quantity at boundaries
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               C1(k,j,i) = a(k,j,i)
            ENDDO
         ENDDO
      ENDDO
      CALL fftp1d_complex_to_real_z(planfc,C1,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth)
         DO j = 1,ny
            R1(j,i) = real(C1(1,j,i), kind=GP)**2 + aimag(c1(1,j,i))**2
            R2(j,i) = real(C1(nz-Cz,j,i), kind=GP)**2 +&
                        aimag(c1(nz-Cz,j,i))**2
         ENDDO
      ENDDO

      IF ( PRESENT(b) ) THEN
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  C1(k,j,i) = b(k,j,i)
               ENDDO
            ENDDO
         ENDDO
         CALL fftp1d_complex_to_real_z(planfc,C1,MPI_COMM_WORLD)
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
      ENDIF

      ! Get mean value
      IF (ista.eq.1) THEN
!$omp parallel do reduction(+:dloc,eloc,floc,gloc)
         DO j = 1,ny
               botloc = botloc + R1(j,1)*tmp
               toploc = toploc + R2(j,1)*tmp
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j) &
!$omp& reduction(+:dloc,eloc,floc,gloc)
         DO i = 2,iend
            DO j = 1,ny
                  botloc = botloc + 2*R1(j,i)*tmp
                  toploc = toploc + 2*R2(j,i)*tmp
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j) &
!$omp& reduction(+:dloc,eloc,floc,gloc)
         DO i = ista,iend
           DO j = 1,ny
               botloc = botloc + 2*R1(j,i)*tmp
               toploc = toploc + 2*R2(j,i)*tmp
           END DO
        END DO
      ENDIF

      CALL MPI_REDUCE(botloc,bot,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
          MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(toploc,top,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
          MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE bouncheck_z

!***************************************************************************
      FUNCTION preprocess(string) RESULT(output)
!---------------------------------------------------------------------------
!  Returns a string preprocessed and ready for comparison with another
! alphabetic lower case string.
!---------------------------------------------------------------------------
      IMPLICIT NONE
         CHARACTER(len=*), INTENT(in) :: string
         CHARACTER(len=len(string))   :: output

         output = lowcase(trim(adjustl(string)))

      END FUNCTION preprocess

!***************************************************************************
      FUNCTION lowcase(string) RESULT(lower)
!---------------------------------------------------------------------------
!  Returns a lower case version of a string. Assumes all characters are 
!  ASCII encoded. Based on: https://www.star.le.ac.uk/%7ecgp/fortran.html
!---------------------------------------------------------------------------
         IMPLICIT NONE
         CHARACTER(len=*), INTENT(in) :: string
         CHARACTER(len=len(string))   :: lower
         INTEGER :: j

         DO j = 1,len(string)
         IF(string(j:j) >= "A" .AND. string(j:j) <= "Z") then
            lower(j:j) = achar(iachar(string(j:j)) + 32)
         ELSE
            lower(j:j) = string(j:j)
         END IF
         END DO
      END FUNCTION lowcase

ENDMODULE boundary
