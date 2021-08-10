!======================================================================
! SPECTER boundary conditions for the scalar field.
!
!  -1 : Periodic boundary conditions.
!   0 : Constant concentration (or temperature) boundary conditions.
!
!  References:
!  - Constant: https://doi.org/10.1016/j.cpc.2020.107482 
!
! 2020 Mauro Fontana. DF-UBA.
!======================================================================

!***********************************************************************
      FUNCTION s_parsebc(str) RESULT(bc)
!-----------------------------------------------------------------------
!  Function to map a given string value denoting a BC type into an 
!  integer for internal module usage.
!-----------------------------------------------------------------------
      USE mpivars
      USE commtypes
      IMPLICIT NONE

      CHARACTER(len=*), INTENT(IN) :: str
      INTEGER                      :: bc

      IF ( str .eq. 'constant' ) THEN
         bc = 0
      ELSEIF ( str .eq. 'periodic' ) THEN
         bc = -1
      ELSE
         IF ( myrank .eq. 0 ) PRINT*, "[ERROR] Unknown boundary ",& 
             "condition type", str, " Aborting..."
         CALL MPI_FINALIZE(MPI_COMM_WORLD, ierr) 
         STOP
      ENDIF
      RETURN
      ENDFUNCTION s_parsebc

!***************************************************************************
      SUBROUTINE s_setup(this,planfc,bckind)
!---------------------------------------------------------------------------
!  Add velocity field and boundary conditions to the boundary plan.
!  ARGUMENTS:
!    this   : 'this' instance of a BCPLAN. [INOUT]
!    fcplan : FCPLAN that will be used for the boundary conditions. [INOUT]
!    bckind : Array of strings with the kind of boundary condition
!             for each boundary. [INOUT]
!---------------------------------------------------------------------------
      USE fcgram
 
      IMPLICIT NONE
 
      TYPE(BCPLAN), INTENT(INOUT)   :: this
      TYPE(FCPLAN), INTENT(INOUT)   :: planfc
      CHARACTER(len=15), INTENT(IN) :: bckind(6)

      this%bcxsta = s_parsebc(preprocess(bckind(1)))
      this%bcxend = s_parsebc(preprocess(bckind(2)))
      this%bcysta = s_parsebc(preprocess(bckind(3)))
      this%bcyend = s_parsebc(preprocess(bckind(4)))
      this%bczsta = s_parsebc(preprocess(bckind(5)))
      this%bczend = s_parsebc(preprocess(bckind(6)))

      END SUBROUTINE s_setup

!***********************************************************************
      SUBROUTINE s_imposebc(planbc,planfc,th)
!-----------------------------------------------------------------------
!   Subroutine to apply boundary conditions to the velocity field and
!  project the field into a solenoidal space, for the latter a Poisson
!  equation with appropriate boundary conditions is solved.
!  ARGUMENTS:
!    planbc : BCPLAN to use during the imposition. [IN]
!    planfc : FCPLAN to use during the imposition. [IN]
!    th     : Array containing the scala field field in the
!              (kz,ky,kx) domain. [INOUT]
!-----------------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE commtypes
      USE fcgram

      IMPLICIT NONE

      TYPE(BCPLAN), INTENT(IN)  :: planbc
      TYPE(FCPLAN), INTENT(IN)  :: planfc
      COMPLEX(KIND=GP),INTENT(INOUT),DIMENSION(planfc%z%n,planfc%y%n,ista:iend) :: th

!TODO Create interface for the case of non-periodic BC in X
!(in that case pr should be real)

      ! Check BC in Y direction
      IF ( planbc%bcysta .ne. -1 .OR. planbc%bcyend .ne. -1) THEN
         IF ( myrank .eq. 0 ) PRINT*, "[ERROR] Non-periodic boundary ", &
             "conditions in Y direction not supported yet. Aborting..."
         CALL MPI_FINALIZE(MPI_COMM_WORLD, ierr) 
         STOP
      ENDIF
      ! Check BC in Z direction
      IF ( planbc%bczsta .ne. 0 .OR. planbc%bczend .ne. 0 ) THEN
         IF ( myrank .eq. 0 ) PRINT*, "[ERROR] Unsupported boundary ", &
             "conditions in Z direction. Aborting..."
         CALL MPI_FINALIZE(MPI_COMM_WORLD, ierr) 
         STOP
      ENDIF

      CALL goto_domain_w_boundaries(planbc,planfc,th)

      ! Z=0
      IF (planbc%bczsta .eq. 0) CALL s_constant_z(planfc,th,0)

      ! Z=Lz
      IF (planbc%bczend .eq. 0) CALL s_constant_z(planfc,th,1)

      ! Back to 3D Fourier
      CALL goto_3d_fourier(planbc,planfc,th)

      RETURN
      END SUBROUTINE s_imposebc
      
!***********************************************************************
      SUBROUTINE s_constant_z(planfc,th,pos)
!-----------------------------------------------------------------------
!  Subroutine to apply non-slip boundary conditions in top and bottom
!  boundaries.
!  ARGUMENTS:
!    fcplan : FCPLAN to use during the imposition. [IN]
!    th     : Array containing the scalar field in the
!              (kz,ky,kx) domain. [INOUT]
!    pos    : 0 for bottom boundary, 1 for top boundary.  [IN]
!-----------------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE commtypes
      USE fcgram
      USE kes
      USE var
      USE order
!$    USE threads

      IMPLICIT NONE

      TYPE(FCPLAN), INTENT(IN)  :: planfc
      COMPLEX(KIND=GP),INTENT(INOUT),DIMENSION(planfc%z%n,planfc%y%n,ista:iend) :: th
      INTEGER, INTENT(IN)  :: pos

      REAL(KIND=GP) :: tmp
      INTEGER       :: ind
      INTEGER       :: i,j

      IF (pos .eq. 0) THEN
         ind = 1
      ELSEIF ( pos .eq. 1 ) THEN
         ind = planfc%z%n - planfc%z%c
      ENDIF

 !$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
 !$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,planfc%y%n
            th(ind,j,i) = real(0.0, kind=GP)  !TODO get value at call time
         ENDDO
      ENDDO

      END SUBROUTINE s_constant_z

!*****************************************************************
      SUBROUTINE sdiagnostic(planbc,planfc,a,t,dt)
!-----------------------------------------------------------------
!  Computes the accuracy of the boundary conditions for the scalar.
!
! Output files contain:
! 'isothermal_diagnostic.txt': time, <|th|^2>|z=0, <|th|^2>|z=Lz
!
! Parameters
!     a  : scalar field 
!     t  : number of time steps made
!     dt : time step
!
      USE fprecision
      USE fcgram
      USE grid
      USE mpivars

      IMPLICIT NONE

      TYPE(BCPLAN), INTENT(IN)  :: planbc
      TYPE(FCPLAN), INTENT(IN)  :: planfc

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a
      REAL(KIND=GP), INTENT(IN)        :: dt
      INTEGER, INTENT(IN)              :: t

      DOUBLE PRECISION                 :: tmp,tmq

      IF ( (planbc%bcxsta .eq. 0) .OR. (planbc%bcxend .eq. 0) .OR. &
           (planbc%bcysta .eq. 0) .OR. (planbc%bcyend .eq. 0) .OR. &
           (planbc%bczsta .eq. 0) .OR. (planbc%bczend .eq. 0) ) THEN

         CALL bouncheck_z(planfc,tmp,tmq,a)
         IF ( myrank .eq. 0) THEN
            OPEN(1,file='scalar_constant_diagnostic.txt',position='append')
10             FORMAT( 1P 3E13.6 )
               WRITE(1,10) (t-1)*dt,tmp,tmq
            CLOSE(1)
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE sdiagnostic
