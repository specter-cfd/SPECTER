!======================================================================
! SPECTER boundary conditions for the magnetic vector potential.
!
!  -1 : Periodic boundary conditions.
!   0 : Conducting boundary conditions.
!   1 : Insulating boundary conditions.
!
!  References:
!  - TBD
!
! 2020 Mauro Fontana. DF-UBA.
!======================================================================

!***********************************************************************
      FUNCTION b_parsebc(str) RESULT(bc)
!-----------------------------------------------------------------------
!  Function to map a given string value denoting a BC type into an 
!  integer for internal module usage.
!-----------------------------------------------------------------------
      USE mpivars
      USE commtypes
      IMPLICIT NONE

      CHARACTER(len=*), INTENT(IN) :: str
      INTEGER                      :: bc

      IF ( str .eq. 'conducting' ) THEN
         bc = 0
      ELSEIF ( str .eq. 'vacuum' ) THEN
         bc = 1
      ELSEIF ( str .eq. 'periodic' ) THEN
         bc = -1
      ELSE
         IF ( myrank .eq. 0 ) PRINT*, "[ERROR] Unknown boundary ",& 
             "condition type", str, " Aborting..."
         CALL MPI_FINALIZE(MPI_COMM_WORLD, ierr) 
         STOP
      ENDIF
      RETURN
      ENDFUNCTION b_parsebc

!***************************************************************************
      SUBROUTINE b_setup(this,planfc,bckind,vplan)
!---------------------------------------------------------------------------
!  Add vector potential and boundary conditions to the boundary plan.
!  ARGUMENTS:
!    this   : 'this' instance of a BCPLAN. [INOUT]
!    fcplan : FCPLAN that will be used for the boundary conditions. [INOUT]
!    bckind : Array of strings with the kind of boundary condition
!             for each boundary. [INOUT]
!---------------------------------------------------------------------------
      USE fcgram
      USE mpivars
      USE commtypes
 
      IMPLICIT NONE
 
      TYPE(BCPLAN), INTENT(INOUT)   :: this
      TYPE(FCPLAN), INTENT(INOUT)   :: planfc
      CHARACTER(len=15), INTENT(IN) :: bckind(6)

      TYPE(BCPLAN), INTENT(IN), OPTIONAL   :: vplan


      this%bcxsta = b_parsebc(preprocess(bckind(1)))
      this%bcxend = b_parsebc(preprocess(bckind(2)))
      this%bcysta = b_parsebc(preprocess(bckind(3)))
      this%bcyend = b_parsebc(preprocess(bckind(4)))
      this%bczsta = b_parsebc(preprocess(bckind(5)))
      this%bczend = b_parsebc(preprocess(bckind(6)))

      IF (PRESENT(vplan)) THEN
         !TODO Add possibility of non-slip with moving boundary
         IF (vplan%bczsta .ne. 0 .AND. vplan%bczend .ne. 0) THEN  
            IF ( myrank .eq. 0 ) PRINT*, "[ERROR] Conducting and insulating ", &
                "boundary conditions with non-steady walls not supported yet. ", &
                "Aborting..."
            CALL MPI_FINALIZE(MPI_COMM_WORLD, ierr) 
            STOP
         ENDIF
      ENDIF

      IF ( this%bczsta.eq.0 .OR. this%bczend.eq.0) THEN
         IF (.NOT. is_neumann_loaded(planfc%z,1) ) THEN
            CALL load_neumann_tables(planfc%z,planfc%tdir,1)
         ENDIF
         IF (.NOT. is_neumann_loaded(planfc%z,2) ) THEN
            CALL load_neumann_tables(planfc%z,planfc%tdir,2)
         ENDIF
      ELSEIF ( this%bczsta.eq.1 .OR. this%bczend.eq.1) THEN
         IF (.NOT. is_neumann_loaded(planfc%z,1) ) THEN
            CALL load_neumann_tables(planfc%z,planfc%tdir,1)
         ENDIF
      ENDIF

      END SUBROUTINE b_setup


!***********************************************************************
      SUBROUTINE a_imposebc_and_project(planbc,planfc,ax,ay,az,ph)
!-----------------------------------------------------------------------
!  Subroutine to apply boundary conditions to the magnetic vector
! potential and project to a solenoidal space, for the latter a Poisson
! equation with appropriate boundary conditions is solved. For the
! conducting case, appropriate intermediate boundary conditions are
! also imposed.
!
!  ARGUMENTS:
!    planbc : BCPLAN to use during the imposition. [IN]
!    planfc : FCPLAN to use during the imposition. [IN]
!    ai     : Array containing the ith component of the magnetic vector
!             potential in the (kz,ky,kx) domain. [INOUT]
!    ph     : Array containing the electrostatic potential in the
!             (z,ky,kx) domain or (z,y,kx) domain if y is 
!             non-periodic. [IN]
!-----------------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE commtypes
      USE fcgram

      IMPLICIT NONE

      TYPE(BCPLAN), INTENT(IN)  :: planbc
      TYPE(FCPLAN), INTENT(IN)  :: planfc
      COMPLEX(KIND=GP),INTENT(INOUT),DIMENSION(planfc%z%n,planfc%y%n,ista:iend) :: ax,ay,az,ph

!TODO Create interface for the case of non-periodic BC in X
!(in that case ph should be real)

      ! Check BC in Y direction
      IF ( planbc%bcysta .ne. -1 .OR. planbc%bcyend .ne. -1) THEN
         IF ( myrank .eq. 0 ) PRINT*, "[ERROR] Non-periodic boundary ", &
             "conditions in Y direction not supported yet. Aborting..."
         CALL MPI_FINALIZE(MPI_COMM_WORLD, ierr) 
         STOP
      ENDIF
      ! Check BC in Z direction
      IF ( (planbc%bczsta .lt. 0 .OR. planbc%bczsta .gt. 1) .OR. &
           (planbc%bczend .lt. 0 .OR. planbc%bczend .gt. 1) ) THEN
         IF ( myrank .eq. 0 ) PRINT*, "[ERROR] Unsupported boundary ", &
             "conditions in Z direction. Aborting..."
         CALL MPI_FINALIZE(MPI_COMM_WORLD, ierr) 
         STOP
      ENDIF

      ! Set zero mode for Az
      IF (myrank .eq. 0) THEN
         az(1,1,1) = 0.0_GP
      ENDIF

      ! Intermediate boundary conditions
      IF (planbc%bczsta .eq. 0 .OR. planbc%bczend .eq. 0) THEN
         ! Go to nearest domain with boudaries only for conducting walls
         CALL goto_domain_w_boundaries(planbc,planfc,ax,ay)

         IF (planbc%bczsta .eq. 0) CALL int_conducting_z(planfc,ax,ay,0)
         IF (planbc%bczend .eq. 0) CALL int_conducting_z(planfc,ax,ay,1)

         ! Back to 3D Fourier
         CALL goto_3d_fourier(planbc,planfc,ax,ay)
      ENDIF

      ! Project using appropriate BC
      CALL sol_project(ax,ay,az,ph,0,2*planbc%bczsta,2*planbc%bczend)


      ! Go to nearest domain with boudaries
      CALL goto_domain_w_boundaries(planbc,planfc,ax,ay,az)

      ! Z=0
      IF (planbc%bczsta .eq. 0) THEN
         CALL conducting_z(planfc,ax,ay,az,0)
      ELSEIF (planbc%bczsta .eq. 1) THEN
         CALL insulating_z(planfc,ax,ay,az,0)
      ENDIF

      ! Z=Lz
      IF (planbc%bczend .eq. 0) THEN
         CALL conducting_z(planfc,ax,ay,az,1)
      ELSEIF (planbc%bczend .eq. 1) THEN
         CALL insulating_z(planfc,ax,ay,az,1)
      ENDIF

      ! Back to 3D Fourier
      CALL goto_3d_fourier(planbc,planfc,ax,ay,az)

      RETURN
      END SUBROUTINE a_imposebc_and_project

!***********************************************************************
      SUBROUTINE int_conducting_z(planfc,ax,ay,pos)
!-----------------------------------------------------------------------
!   Subroutine to apply intermediate conducting boundary conditions in
!  top and bottom boundaries.
!  ARGUMENTS:
!    fcplan : FCPLAN to use during the imposition. [IN]
!    ai     : Array containing the ith component of the vector potential
!              (kz,ky,kx) domain. [INOUT]
!    pos    : 0 for bottom boundary, 1 for top boundary.
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
      COMPLEX(KIND=GP),INTENT(INOUT),DIMENSION(planfc%z%n,planfc%y%n,ista:iend) :: ax,ay

      INTEGER, INTENT(IN)  :: pos

      INTEGER  :: ind
      INTEGER  :: i,j

      IF (pos .eq. 0) THEN
         ind = 1
      ELSEIF ( pos .eq. 1 ) THEN
         ind = planfc%z%n - planfc%z%c
      ENDIF

 !$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
 !$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,planfc%y%n
            ax(ind,j,i) = 0.0_GP 
            ay(ind,j,i) = 0.0_GP
          ENDDO
      ENDDO
 
      END SUBROUTINE int_conducting_z

!***********************************************************************
      SUBROUTINE conducting_z(planfc,ax,ay,az,pos)
!-----------------------------------------------------------------------
!  Subroutine to apply conducting boundary conditions in top and bottom
!  boundaries.
!  ARGUMENTS:
!    fcplan : FCPLAN to use during the imposition. [IN]
!    ai     : Array containing the ith component of the vector potential
!              (kz,ky,kx) domain. [INOUT]
!    pos    : 0 for bottom boundary, 1 for top boundary.
!    inter  : 0 for final BC, 1 for intermediate BC
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
      COMPLEX(KIND=GP),INTENT(INOUT),DIMENSION(planfc%z%n,planfc%y%n,ista:iend) :: ax,ay,az

      INTEGER, INTENT(IN)  :: pos

      INTEGER  :: ind
      INTEGER  :: i,j

      IF (pos .eq. 0) THEN
         ind = 1
      ELSEIF ( pos .eq. 1 ) THEN
         ind = planfc%z%n - planfc%z%c
      ENDIF

 !$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
 !$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,planfc%y%n
            ax(ind,j,i) = 0.0_GP 
            ay(ind,j,i) = 0.0_GP
            az(ind,j,i) = 0.0_GP
         ENDDO
      ENDDO

      ! Second order homogeneous Neumann for a_para
      ! and first order for a_perp
      CALL neumann_reconstruct(planfc,ax,5+pos,2)
      CALL neumann_reconstruct(planfc,ay,5+pos,2)
      CALL neumann_reconstruct(planfc,az,5+pos,1)

      END SUBROUTINE conducting_z


!***********************************************************************
      SUBROUTINE insulating_z(planfc,ax,ay,az,pos)
!-----------------------------------------------------------------------
!  Subroutine to apply insulating boundary conditions in top and bottom
!  boundaries.
!  ARGUMENTS:
!    fcplan : FCPLAN to use during the imposition. [IN]
!    ai     : Array containing the ith component of the vector potential
!              (kz,ky,kx) domain. [INOUT]
!    pos    : 0 for bottom boundary, 1 for top boundary.
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
      COMPLEX(KIND=GP),INTENT(INOUT),DIMENSION(planfc%z%n,planfc%y%n,ista:iend) :: ax,ay,az
      INTEGER, INTENT(IN)  :: pos

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
            ax(ind,j,i) = 0.0_GP 
            ay(ind,j,i) = 0.0_GP
            az(ind,j,i) = 0.0_GP
         ENDDO
      ENDDO

      ! Second order homogeneous Neumann for a_para
      ! and first order for a_perp
      CALL robin_reconstruct(planfc,ax,5+pos,khom)
      CALL robin_reconstruct(planfc,ay,5+pos,khom)
      CALL robin_reconstruct(planfc,az,5+pos,khom)

      END SUBROUTINE insulating_z


!*****************************************************************
      SUBROUTINE bdiagnostic(planbc,planfc,a,b,c,t,dt)
!-----------------------------------------------------------------
!  Computes the mean squared divergence of both the magnetic field
! and the vector potential inside the domain.
!  Also printed are the accuracy of the boundary conditions. Note
! that if mixed boundary conditions are employed, only certain 
! columns of each file will contain the relevant boundary information.
!
! Output files contain:
! 'conducting_diagnostic.txt': time, <|div(a)|^2>, <|div(b)|^2>, <|j_t|^2>|z=0,
!                   <|j_t|^2>|z=Lz, <|b_n|^2>|z=0, <|b_n|^2>|z=Lz
!
! 'vacuum_diagnostic.txt': time, <|div(a)|^2>, <|div(b)|^2>, <|bc_t|^2>|z=0,
!                   <|bc_t|^2>|z=Lz, <|bc_n|^2>|z=0, <|bc_n|^2>|z=Lz
!
!    where in both cases n and t denote the normal and tangential
!    part of the vector field, respectively. In the vacuum case 
!    bc = da/dn + sqrt(kx**2+ky**2)*a. 
!
! Parameters
!     a   : magnetic vector potential in the x-direction
!     b   : magnetic vector potential in the y-direction
!     c   : magnetic vector potential in the z-direction
!     t   : number of time steps made
!     dt  : time step
!
      USE fprecision
      USE grid
      USE fcgram
      USE mpivars

      IMPLICIT NONE

      TYPE(BCPLAN), INTENT(IN)  :: planbc
      TYPE(FCPLAN), INTENT(IN)  :: planfc

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      REAL(KIND=GP), INTENT(IN)        :: dt
      INTEGER, INTENT(IN)              :: t

      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: c1,c2,c3,c4

      DOUBLE PRECISION                 :: tmp,tmq,tmr,tms,tm1,tm2

      ! Get magnetic field
      CALL curlk(b,c,c1,1)
      CALL curlk(a,c,c2,2)
      CALL curlk(a,b,c3,3)

      CALL divergence(a,b,c,tm1)
      CALL divergence(c1,c2,c3,tm2)


10    FORMAT( 1P 7E13.6 )
      IF ( (planbc%bcxsta .eq. 0) .OR. (planbc%bcxend .eq. 0) .OR. &
           (planbc%bcysta .eq. 0) .OR. (planbc%bcyend .eq. 0) .OR. &
           (planbc%bczsta .eq. 0) .OR. (planbc%bczend .eq. 0) ) THEN
         CALL bouncheck_z(planfc,tmr,tms,c3)

         ! Current
         CALL curlk(c2,c3,c4,1)
         CALL curlk(c1,c3,c2,2)
         CALL bouncheck_z(planfc,tmp,tmq,c4,c2)

         IF ( myrank .eq. 0) THEN
            OPEN(1,file='conducting_diagnostic.txt',position='append')
               WRITE(1,10) (t-1)*dt,tm1,tm2,tmp,tmq,tmr,tms
            CLOSE(1)
         ENDIF

      ELSEIF ( (planbc%bcxsta .eq. 1) .OR. (planbc%bcxend .eq. 1) .OR. &
               (planbc%bcysta .eq. 1) .OR. (planbc%bcyend .eq. 1) .OR. &
               (planbc%bczsta .eq. 1) .OR. (planbc%bczend .eq. 1) ) THEN
         CALL robcheck(planfc,a,b,c,tmp,tmq,tmr,tms)
         IF ( myrank .eq. 0) THEN
            OPEN(1,file='vacuum_diagnostic.txt',position='append')
               WRITE(1,10) (t-1)*dt,tm1,tm2,tmp,tmq,tmr,tms
            CLOSE(1)
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE bdiagnostic


!*****************************************************************
      SUBROUTINE robcheck(planfc,a,b,c,d,e,f,g)
!-----------------------------------------------------------------
!
!  Computes the mean squared values of error in the vacuum Robin
! boudary condition, both for the tangential and the normal
! component of the field.
! The output is only valid in the first node.
!
! Parameters:
!     a: x compoenent of the vector potential 
!     b: y compoenent of the vector potential 
!     c: z compoenent of the vector potential 
!     d: at the output contains the mean squared value of the
!          error in the tangential component at z=0.
!     e: at the output contains the mean squared value of the
!          error in the tangential component at z=Lz.
!     f: at the output contains the mean squared value of the
!          error in the normal component at z=0.
!     g: at the output contains the mean squared value of the
!          error in the normal component at z=Lz.
!
      USE var
      USE grid
      USE fcgram
      USE mpivars
      USE commtypes
      USE kes
!$    USE threads
      IMPLICIT NONE

      TYPE(FCPLAN), INTENT(IN)  :: planfc

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      DOUBLE PRECISION, INTENT(OUT)    :: d,e,f,g

      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)  :: C1,C2
      REAL(KIND=GP), DIMENSION(ny,ista:iend)        :: R1,R2,R3,R4

      DOUBLE PRECISION   :: dloc,eloc,floc,gloc
      REAL(KIND=GP)      :: tmp

      INTEGER  :: i,j,k

      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2

      d = 0d00; e = 0d0; f = 0d0; g = 0d0
      dloc = 0d0; eloc = 0d0; floc = 0d0; gloc = 0d0

      ! Ax
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               C1(k,j,i) = a(k,j,i)
            ENDDO
         ENDDO
      ENDDO
      CALL derivk(C1,C2,3)
      CALL fftp1d_complex_to_real_z(planfc,C1,MPI_COMM_WORLD)
      CALL fftp1d_complex_to_real_z(planfc,C2,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth)
         DO j = 1,ny
            R1(j,i) = real(-C2(1,j,i)+khom(j,i)*C1(1,j,i), kind=GP)**2 +&
                      aimag(-C2(1,j,i)+khom(j,i)*C1(1,j,i))**2

            R2(j,i) = real(C2(nz-Cz,j,i)+khom(j,i)*C1(nz-Cz,j,i),kind=GP)**2 +&
                      aimag(C2(nz-Cz,j,i)+khom(j,i)*C1(nz-Cz,j,i))**2
         ENDDO
      ENDDO

      ! Ay
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               C1(k,j,i) = b(k,j,i)
            ENDDO
         ENDDO
      ENDDO
      CALL derivk(C1,C2,3)
      CALL fftp1d_complex_to_real_z(planfc,C1,MPI_COMM_WORLD)
      CALL fftp1d_complex_to_real_z(planfc,C2,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth)
         DO j = 1,ny
            R1(j,i) = R1(j,i) + real(-C2(1,j,i)+khom(j,i)*C1(1,j,i), kind=GP)**2 +&
                      aimag(-C2(1,j,i)+khom(j,i)*C1(1,j,i))**2

            R2(j,i) = R2(j,i) + real(C2(nz-Cz,j,i)+khom(j,i)*C1(nz-Cz,j,i),kind=GP)**2 +&
                      aimag(C2(nz-Cz,j,i)+khom(j,i)*C1(nz-Cz,j,i))**2
         ENDDO
      ENDDO

      ! Az
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               C1(k,j,i) = c(k,j,i)
            ENDDO
         ENDDO
      ENDDO
      CALL derivk(C1,C2,3)
      CALL fftp1d_complex_to_real_z(planfc,C1,MPI_COMM_WORLD)
      CALL fftp1d_complex_to_real_z(planfc,C2,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth)
         DO j = 1,ny
            R3(j,i) = real(-C2(1,j,i)+khom(j,i)*C1(1,j,i), kind=GP)**2 +&
                      aimag(-C2(1,j,i)+khom(j,i)*C1(1,j,i))**2

            R4(j,i) = real(C2(nz-Cz,j,i)+khom(j,i)*C1(nz-Cz,j,i),kind=GP)**2 +&
                      aimag(C2(nz-Cz,j,i)+khom(j,i)*C1(nz-Cz,j,i))**2
         ENDDO
      ENDDO

      !
      ! Compute average
      !
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
      END SUBROUTINE robcheck
