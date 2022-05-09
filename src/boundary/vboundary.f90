!======================================================================
! SPECTER boundary conditions for the velocity field.
!
!  -1 : Periodic boundary conditions.
!   0 : Non-slip boundary conditions.
!
!  References:
!  - No-Slip: https://doi.org/10.1016/j.cpc.2020.107482 
!
! 2020 Mauro Fontana. DF-UBA.
!======================================================================

!***********************************************************************
      FUNCTION v_parsebc(str) RESULT(bc)
!-----------------------------------------------------------------------
!  Function to map a given string value denoting a BC type into an 
!  integer for internal module usage.
!-----------------------------------------------------------------------
      USE mpivars
      USE commtypes
      IMPLICIT NONE

      CHARACTER(len=*), INTENT(IN) :: str
      INTEGER                      :: bc

      IF ( str .eq. 'noslip' ) THEN
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
      ENDFUNCTION v_parsebc

!***************************************************************************
      SUBROUTINE v_setup(this,planfc,bckind)
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

      this%bcxsta = v_parsebc(preprocess(bckind(1)))
      this%bcxend = v_parsebc(preprocess(bckind(2)))
      this%bcysta = v_parsebc(preprocess(bckind(3)))
      this%bcyend = v_parsebc(preprocess(bckind(4)))
      this%bczsta = v_parsebc(preprocess(bckind(5)))
      this%bczend = v_parsebc(preprocess(bckind(6)))

      END SUBROUTINE v_setup

!***********************************************************************
      SUBROUTINE v_imposebc_and_project(planbc,planfc,vx,vy,vz,pr,rki, &
                                        v_zsta,v_zend)
!-----------------------------------------------------------------------
!   Subroutine to apply boundary conditions to the velocity field and
!  project the field into a solenoidal space, for the latter a Poisson
!  equation with appropriate boundary conditions is solved.
!  ARGUMENTS:
!    planbc : BCPLAN to use during the imposition. [IN]
!    planfc : FCPLAN to use during the imposition. [IN]
!    vi     : Array containing the ith component of the velocity field
!              (kz,ky,kx) domain. [INOUT]
!    pr     : Array containing the pressure in the (z,ky,kx) domain
!             or (z,y,kx) domain if y is non-periodic. [IN]
!    rki    : Current iteration of Runge-Kutta scheme. (OPTIONAL) [IN]
!    v_zsta : x and y components of v at z=0.  (OPTIONAL) [IN]
!    v_zend : x and y components of v at z=Lz. (OPTIONAL) [IN]
!-----------------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE commtypes
      USE fcgram

      IMPLICIT NONE

      TYPE(BCPLAN), INTENT(IN)  :: planbc
      TYPE(FCPLAN), INTENT(IN)  :: planfc
      COMPLEX(KIND=GP),INTENT(INOUT),DIMENSION(planfc%z%n,planfc%y%n,ista:iend) :: vx,vy,vz,pr

      REAL(KIND=GP), INTENT(IN), DIMENSION(2), OPTIONAL :: v_zsta, v_zend
      INTEGER, INTENT(IN), OPTIONAL       :: rki

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

      CALL goto_domain_w_boundaries(planbc,planfc,vx,vy)

      IF ( (planbc%bczsta .eq. 0 .OR. planbc%bczend .eq. 0) .AND. &
            .NOT. PRESENT(rki) ) THEN
         IF ( myrank .eq. 0 ) PRINT*, "[ERROR] Non-slip BC require `rki` ", &
             "keyword argument in call to imposebc_and_project. Aborting..."
         CALL MPI_FINALIZE(MPI_COMM_WORLD, ierr) 
         STOP
      ENDIF

      ! Z=0
      IF (planbc%bczsta .eq. 0) THEN
          IF( .NOT. PRESENT(v_zsta) ) THEN
             CALL noslip_z(planfc,rki,vx,vy,vz,pr,(/.0_GP,.0_GP/),0)
          ELSE
             CALL noslip_z(planfc,rki,vx,vy,vz,pr,v_zsta,0)
          ENDIF
      ENDIF

      ! Z=Lz
      IF (planbc%bczend .eq. 0) THEN
          IF( .NOT. PRESENT(v_zend) ) THEN
             CALL noslip_z(planfc,rki,vx,vy,vz,pr,(/.0_GP,.0_GP/),1)
          ELSE
             CALL noslip_z(planfc,rki,vx,vy,vz,pr,v_zend,1)
          ENDIF
      ENDIF

      ! Back to 3D Fourier
      CALL goto_3d_fourier(planbc,planfc,vx,vy)

      ! Now project
      CALL sol_project(vx,vy,vz,pr,1,0,0)

      RETURN
      END SUBROUTINE v_imposebc_and_project
      
!***********************************************************************
      SUBROUTINE noslip_z(planfc,o,vx,vy,vz,pr,vbound,pos)
!-----------------------------------------------------------------------
!  Subroutine to apply non-slip boundary conditions in top and bottom
!  boundaries.
!  ARGUMENTS:
!    fcplan : FCPLAN to use during the imposition. [IN]
!    o      : Current iteration of the Runge-Kutta scheme. [IN]
!    vi     : Array containing the ith component of the velocity field
!              (kz,ky,kx) domain. [INOUT]
!    pr     : Array containing the pressure in the (z,ky,kx) domain
!             or (z,y,kx) domain if y is non-periodic. [IN]
!    vbound : Tangential velocity of the boundary. [IN]
!    pos    : 0 for bottom boundary, 1 for top boundary. [IN]
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
      COMPLEX(KIND=GP),INTENT(INOUT),DIMENSION(planfc%z%n,planfc%y%n,ista:iend) :: vx,vy,vz
      COMPLEX(KIND=GP),INTENT(IN),DIMENSION(planfc%z%n,planfc%y%n,ista:iend)    :: pr
      REAL(KIND=GP), INTENT(IN), DIMENSION(2) :: vbound
      INTEGER, INTENT(IN)  :: o,pos

      REAL(KIND=GP) :: tmp
      INTEGER       :: ind
      INTEGER       :: i,j

      IF (pos .eq. 0) THEN
         ind = 1
      ELSEIF ( pos .eq. 1 ) THEN
         ind = planfc%z%n - planfc%z%c
      ENDIF

      tmp = 1/real(o, kind=GP)
      IF (o .ne. ord) tmp = real(o+1, kind=GP)*tmp

 !$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
 !$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,planfc%y%n
            vx(ind,j,i) = im*kx(i)*pr(ind,j,i)*tmp
            vy(ind,j,i) = im*ky(j)*pr(ind,j,i)*tmp
         ENDDO
      ENDDO

      ! Add velocity of the boundary
      IF ( ista .eq. 1 ) vx(ind,1,1) = planfc%x%n * planfc%y%n * vbound(1)
      IF ( ista .eq. 1 ) vy(ind,1,1) = planfc%x%n * planfc%y%n * vbound(2)

      END SUBROUTINE noslip_z

!*****************************************************************
      SUBROUTINE vdiagnostic(planbc,planfc,a,b,c,t,dt)
!-----------------------------------------------------------------
!  Computes the mean squared divergence of the velocity field
! inside the domain.
!  Also printed are the accuracy of the boundary conditions. Note
! that if mixed boundary conditions are employed, only certain 
! columns of each file will contain the relevant boundary information.
!
! Output files contain:
! 'noslip_diagnostic.txt': time, <|div(v)|^2>, <|v_t|^2>|z=0,
!                     <|v_t|^2>|z=Lz, <|v_n|^2>|z=0, <|v_n|^2>|z=Lz
!
!    where in both cases n and t denote the normal and tangential
!    part of the vector field, respectively
!
! Parameters
!     a  : velocity field in the x-direction
!     b  : velocity field in the y-direction
!     c  : velocity field in the z-direction
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

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      REAL(KIND=GP), INTENT(IN)        :: dt
      INTEGER, INTENT(IN)              :: t

      DOUBLE PRECISION                 :: tmp,tmq,tmr,tms,tm1

      CALL divergence(a,b,c,tm1)

      IF ( (planbc%bcxsta .eq. 0) .OR. (planbc%bcxend .eq. 0) .OR. &
           (planbc%bcysta .eq. 0) .OR. (planbc%bcyend .eq. 0) .OR. &
           (planbc%bczsta .eq. 0) .OR. (planbc%bczend .eq. 0) ) THEN

         CALL bouncheck_z(planfc,tmp,tmq,a,b)
         CALL bouncheck_z(planfc,tmr,tms,c)
         IF ( myrank .eq. 0) THEN
            OPEN(1,file='noslip_diagnostic.txt',position='append')
10             FORMAT( 1P 7E13.6 )
               WRITE(1,10) (t-1)*dt,tm1,tmp,tmq,tmr,tms
            CLOSE(1)
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE vdiagnostic
