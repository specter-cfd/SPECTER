!======================================================================
! SPECTER fcgram module
!
!   Implements helper datatypes and subroutines used to enforce
! different kinds of boundary conditions using an FC-Gram method.
!   Module-wide notation: boundaries 1, 2, 3, 4, 5, 6 represent
! boundaries located at x=0, x=Lx, y=0, y=Ly, z=0, z=Lz, respectively.
!
! TODO Add support for reconstructions in x direction. This would
!      require interfaces to subroutines that operate on real arrays.
!      (Continuation in x should result in real array)
!
! 2020 Mauro Fontana. DF-UBA.
!======================================================================

MODULE fcgram
      USE fprecision
      USE fftplans

      IMPLICIT NONE

      TYPE FCCOORD
         ! Tables used for Dirichlet periodic extension and for
         ! reconstructions using Neumann or second normal derivative
         ! at the boundaries
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: dir
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)    :: neu,neu2
         ! Grid spacing
         REAL(KIND=GP)  :: dx
         ! Total num of grid, continuation and matching points
         INTEGER        :: N,C,d 
      END TYPE FCCOORD

      TYPE FCPLAN
         TYPE(FCCOORD)      :: x,y,z
         TYPE(FFTPLAN)      :: planrc, plancr
         CHARACTER(len=100) :: tdir
      END TYPE FCPLAN

      ! Timers for reconstruction subroutines. The periodic
      ! extension is benchmarked in as part of the FFT.
      INTEGER, SAVE             :: hneu,hrob
      DOUBLE PRECISION, SAVE    :: neutime = 0.0
      DOUBLE PRECISION, SAVE    :: robtime = 0.0


      CONTAINS
!**********************************************************************
      SUBROUTINE fcgram_create_plan(this,N,C,o,dx,tdir,fftflags)
!----------------------------------------------------------------------
!  Main explicit constructor
!  ARGUMENTS:
!    this     : 'this' instance constructor. [OUT]
!    N        : A 3-tuple containing the total number of points (i.e. 
!             the array size) in each direction. [IN]
!    C        : A 3-tuple containing the number of continuation points.
!             (Cx,Cy,Cz). [IN]
!    o        : A 3-tuple containing the accuracy order desired in each
!             direction (ox,oy,oz). [IN]
!    dx       : A 3-tuple containing the spacings (dx,dy,dz). [IN]
!    tdir     : Directory containing the FC-Gram tables. [IN]
!    fftflags : Flags for FFT plans creation (see fftp/fftp.fpp)
!----------------------------------------------------------------------
      USE commtypes
      USE fprecision
      USE mpivars
      USE gtimer

      IMPLICIT NONE

      TYPE(FCPLAN), INTENT(OUT)    :: this
      INTEGER, INTENT(IN)          :: N(3), C(3), o(3), fftflags
      REAL(KIND=GP), INTENT(IN)    :: dx(3)
      CHARACTER(len=*), INTENT(IN) :: tdir

      this%x%dx = dx(1)
      this%y%dx = dx(2)
      this%z%dx = dx(3)
  
      this%x%N  = N(1)
      this%y%N  = N(2)
      this%z%N  = N(3)

      this%x%C  = C(1)
      this%y%C  = C(2)
      this%z%C  = C(3)

      this%x%d  = o(1)
      this%y%d  = o(2)
      this%z%d  = o(3)

      IF (this%y%C .le. 0) THEN
         CALL fftp3d_create_plan(this%planrc,N,-1,fftflags,0)
         CALL fftp3d_create_plan(this%plancr,N,1,fftflags,0)
      ELSE
         CALL fftp3d_create_plan(this%planrc,N,-1,fftflags,1)
         CALL fftp3d_create_plan(this%plancr,N,1,fftflags,1)
      ENDIF


      this%tdir = tdir

      CALL GTStart(hneu,GT_WTIME)
      CALL GTStart(hrob,GT_WTIME)

      ! Load Dirichlet tables for each direction
      ! X
      IF (this%x%C .eq. 0 .AND. this%x%d .eq. 0) THEN
         IF (myrank .eq. 0) PRINT*, "[INFO] Tables for x dimension ", &
            "not loaded as it is a periodic direction."
      ELSE IF (this%x%C .gt. 0 .AND. this%x%d .gt. 0) THEN
         CALL load_dirichlet_tables(this%x, tdir)
      ELSE 
         IF (myrank .eq. 0) PRINT*, "[ERROR] Mismatch in continuation",&
            " or matching points in x direction. Aborting..."
         CALL MPI_FINALIZE(MPI_COMM_WORLD, ierr)
         STOP
      ENDIF
      ! Y
      IF (this%y%C .eq. 0 .AND. this%y%d .eq. 0) THEN
         IF (myrank .eq. 0) PRINT*, "[INFO] Tables for y dimension ", &
            "not loaded as it is a periodic direction."
      ELSE IF (this%y%C .gt. 0 .AND. this%y%d .gt. 0) THEN
         CALL load_dirichlet_tables(this%y, tdir)
      ELSE 
         IF (myrank .eq. 0) PRINT*, "[ERROR] Mismatch in continuation",&
            " or matching points in y direction. Aborting..."
         CALL MPI_FINALIZE(MPI_COMM_WORLD, ierr)
         STOP
      ENDIF
      ! Z
      IF (this%z%C .eq. 0 .AND. this%z%d .eq. 0) THEN
         IF (myrank .eq. 0) PRINT*, "[INFO] Tables for z dimension ", &
            "not loaded as it is a periodic direction."
      ELSE IF (this%z%C .gt. 0 .AND. this%z%d .gt. 0) THEN
         CALL load_dirichlet_tables(this%z, tdir)
      ELSE 
         IF (myrank .eq. 0) PRINT*, "[ERROR] Mismatch in continuation",&
            " or matching points in z direction. Aborting..."
         CALL MPI_FINALIZE(MPI_COMM_WORLD, ierr)
         STOP
      ENDIF

      END SUBROUTINE fcgram_create_plan

!*****************************************************************
      SUBROUTINE fcgram_destroy_plan(this)
!-----------------------------------------------------------------
!  Main explicit destructor
!  ARGUMENTS:
!    this    : 'this' fc_plan instance. [INOUT]
!----------------------------------------------------------------------
      USE gtimer

      IMPLICIT NONE
      TYPE(FCPLAN), INTENT(INOUT) :: this

      IF ( ALLOCATED(this%x%dir) ) DEALLOCATE(this%x%dir)
      IF ( ALLOCATED(this%y%dir) ) DEALLOCATE(this%y%dir)
      IF ( ALLOCATED(this%z%dir) ) DEALLOCATE(this%z%dir)

      IF ( ALLOCATED(this%x%neu) ) DEALLOCATE(this%x%neu)
      IF ( ALLOCATED(this%y%neu) ) DEALLOCATE(this%y%neu)
      IF ( ALLOCATED(this%z%neu) ) DEALLOCATE(this%z%neu)

      IF ( ALLOCATED(this%x%neu2) ) DEALLOCATE(this%x%neu2)
      IF ( ALLOCATED(this%y%neu2) ) DEALLOCATE(this%y%neu2)
      IF ( ALLOCATED(this%z%neu2) ) DEALLOCATE(this%z%neu2)

      CALL fftp3d_destroy_plan(this%plancr)
      CALL fftp3d_destroy_plan(this%planrc)

      CALL GTFree(hneu)
      CALL GTFree(hrob)

      END SUBROUTINE fcgram_destroy_plan


!*****************************************************************
      SUBROUTINE load_dirichlet_tables(coord,tdir)
!-----------------------------------------------------------------
! Helper function to populate Dirichlet tables in a coord
! member of a given fcplan datatype.
! ARGUMENTS:
!     coord : A coord type variable. At the output contains
!             the tables required to compute a periodic extension
!             with Dirichlet boundary conditions. [INOUT]
!     tdir  : directory containing the FC-Gram tables [IN]
!-----------------------------------------------------------------
      USE commtypes
      USE mpivars

      IMPLICIT NONE

      TYPE(FCCOORD), INTENT(INOUT) :: coord
      CHARACTER(len=*), INTENT(IN) :: tdir

      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: A, Q

      INTEGER          :: i,j
      LOGICAL          :: fexists

      CHARACTER(len=256) :: astr,qstr
      CHARACTER(len=5)   :: cstr,dstr

      ALLOCATE( A(coord%C, coord%d), Q(coord%d, coord%d) )

      ! Check for multiple allocations and allocate
      IF ( ALLOCATED(coord%dir) ) THEN
         IF (myrank .eq. 0) PRINT*, "[WARNING] Multiple allocations ", &
             "attempted at load_dir_tables. Skipping..."
      ELSE
          ALLOCATE( coord%dir(coord%C,coord%d) )
      ENDIF


      ! Placeholder strings for numeric and compose quantities
      WRITE(cstr,'(I5)') coord%C
      WRITE(dstr,'(I5)') coord%d

      astr = trim(tdir) // '/A' // trim(adjustl(cstr)) // '-' // &
              trim(adjustl(dstr)) //  '.dat'
      qstr = trim(tdir) // '/Q' // trim(adjustl(dstr)) // '.dat'

      ! Try to load A and Q matrices
      INQUIRE( FILE=trim(astr), EXIST=fexists )
      IF ( fexists ) THEN
         OPEN(10, FILE=trim(astr), FORM='unformatted', ACCESS='stream', &
              ACTION='read')
            READ(10) A
         CLOSE(10)
      ELSE
         IF (myrank .eq. 0) PRINT*, "Could not find table ", trim(astr)
         CALL MPI_FINALIZE(MPI_COMM_WORLD, ierr)
         STOP
      ENDIF

      INQUIRE( FILE=trim(qstr), EXIST=fexists )
      IF ( fexists ) THEN
         OPEN(10, FILE=trim(qstr), FORM='unformatted', ACCESS='stream', &
              ACTION='read')
            READ(10) Q
         CLOSE(10)
      ELSE
         IF (myrank .eq. 0) PRINT*, "Could not find table ", trim(qstr)
         CALL MPI_FINALIZE(MPI_COMM_WORLD, ierr)
         STOP
      ENDIF

      ! Construct continuation operators
      Q = TRANSPOSE(Q)
      coord%dir = MATMUL(A,Q)

      DEALLOCATE( A,Q )

      RETURN
      END SUBROUTINE load_dirichlet_tables


!*****************************************************************
      SUBROUTINE load_neumann_tables(coord,tdir,ord)
!-----------------------------------------------------------------
! Helper function to populate the Neumann reconstructors for a
! coord member of a given fcplan datatype.
! ARGUMENTS:
!     coord : A coord type variable. At the output contains
!             the tables required to compute a periodic extension
!             with normal derivative boundary conditions. [INOUT]
!     tdir  : directory containing the FC-Gram tables. [IN]
!     ord   : Order. 1 for Neumann reconstructor, 2 for second
!             normal derivative reconstructor. [IN]
!-----------------------------------------------------------------

      USE commtypes
      USE mpivars

      IMPLICIT NONE

      INTEGER, INTENT(IN)          :: ord
      CHARACTER(len=*), INTENT(IN) :: tdir
      TYPE(FCCOORD), INTENT(INOUT) :: coord

      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: Q, Qn

      DOUBLE PRECISION   :: dxp
      INTEGER            :: i,j
      LOGICAL            :: fexists
      CHARACTER(len=256) :: qstr,qnstr
      CHARACTER(len=5)   :: dstr

      ALLOCATE( Q(coord%d,coord%d), Qn(coord%d,coord%d) )

      ! Check for multiple allocations and allocate
      IF ( ord .eq. 1 ) THEN
         IF ( ALLOCATED(coord%neu) ) THEN
            IF (myrank .eq. 0) PRINT*, "[WARNING] Multiple allocations ", &
                "attempted at load_neu_tables. Skipping..."
            RETURN
         ELSE
             ALLOCATE( coord%neu(coord%d) )
         ENDIF
      ELSEIF ( ord .eq. 2 ) THEN
         IF ( ALLOCATED(coord%neu2) ) THEN
            IF (myrank .eq. 0) PRINT*, "[WARNING] Multiple allocations ", &
                "attempted at load_neu_tables. Skipping..."
            RETURN
         ELSE
             ALLOCATE( coord%neu2(coord%d) )
         ENDIF
      ELSE
         IF (myrank .eq. 0) PRINT*, "[WARNING] load_neu_tables called with ", &
             "undefined order. Skipping..."
         RETURN
      ENDIF



      WRITE(dstr,'(I5)') coord%d
      qstr = trim(tdir) // '/Q' // trim(adjustl(dstr)) // '.dat'
      IF (ord .eq. 1) THEN
         qnstr = trim(tdir) // '/Q1n' // trim(adjustl(dstr)) // '.dat'
      ELSEIF (ord .eq. 2) THEN
         qnstr = trim(tdir) // '/Q2n' // trim(adjustl(dstr)) // '.dat'
      ENDIF

      ! Try to load Q and Qn matrices
      INQUIRE( FILE=trim(qstr), EXIST=fexists )
      IF ( fexists ) THEN
         OPEN(10, FILE=trim(qstr), FORM='unformatted', ACCESS='stream', &
              ACTION='read')
            READ(10) Q
         CLOSE(10)
      ELSE
         IF (myrank .eq. 0) PRINT*, "Could not find table ", trim(qstr)
         CALL MPI_FINALIZE(MPI_COMM_WORLD, ierr)
         STOP
      ENDIF

      INQUIRE( FILE=trim(qnstr), EXIST=fexists )
      IF ( fexists ) THEN
         OPEN(10, FILE=trim(qnstr), FORM='unformatted', ACCESS='stream', &
              ACTION='read')
            READ(10) dxp
            READ(10) Qn
         CLOSE(10)
      ELSE
         IF (myrank .eq. 0) PRINT*, "Could not find table ", trim(qnstr)
         CALL MPI_FINALIZE(MPI_COMM_WORLD, ierr)
         STOP
      ENDIF

      ! To reconstruct the last value I only need the last row
      ! of the inverse. dx/dxp accounts for different grid spacings
      ! between current grid and the grid used to generate the tables
      IF ( ord .eq. 1 ) THEN
         coord%neu           = MATMUL(Q(coord%d,:), TRANSPOSE(Qn))
         coord%neu(coord%d)  = coord%neu(coord%d) * coord%dx/dxp
      ELSEIF ( ord .eq. 2 ) THEN
         coord%neu2          = MATMUL(Q(coord%d,:), TRANSPOSE(Qn))
         coord%neu2(coord%d) = coord%neu2(coord%d) * (coord%dx/dxp)**2
      ENDIF

      DEALLOCATE( Q, Qn )

      END SUBROUTINE load_neumann_tables

!***********************************************************************
      SUBROUTINE neumann_reconstruct(plan,f,boun,ord)
!-----------------------------------------------------------------------
! Reconstructs function values at the boundary from given its
! ordth-normal derivative at the endpoint d^(ord)_n f = g. The
! array f must contain the values of g at the desired endpoint.
! ARGUMENTS:
!     plan : fcplan to use for the reconstruction.
!     f    : Array whose values must be reconstructed. [INOUT]
!     boun : Boundary to which the reconstruction should 
!            be performed. 1 for x=0, 2 for x=Lx ... 6 for z=Lz [IN]
!     ord  : Order of the derivative. 1 for Neumann, 2 for second
!            normal derivative. [IN]
!-----------------------------------------------------------------------
      USE commtypes
      USE fprecision
      USE mpivars
      USE gtimer
!$    USE threads

      IMPLICIT NONE

      TYPE(FCPLAN), INTENT(IN)                                                :: plan
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(plan%z%N,plan%y%N,ista:iend) :: f
      INTEGER, INTENT(IN)                                                     :: boun, ord

      INTEGER        :: Cy,Cz,dy,dz,ny,nz
      INTEGER        :: i,j,k

      CALL GTStart(hneu)

      ny = plan%y%N
      nz = plan%z%N

      Cy = plan%y%C
      Cz = plan%z%C
      dy = plan%y%d
      dz = plan%z%d

      ! ...........
      ! Y direction
      ! ...........
      IF ( boun .eq. 3 .AND. ord .eq. 1 .AND. ALLOCATED(plan%y%neu) ) THEN
!$omp paralleldo if (iend-ista.ge.nth) private (j,k,tmp)
         DO i=ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,tmp)
            DO j=1,nz-Cz
               f(j,1,i) = plan%y%neu(dy)*f(j,1,i)
               DO k=1,dy-1
                  f(j,1,i) = f(j,1,i) + plan%y%neu(k)*f(j,dy+1-k,i)
               ENDDO
            ENDDO
         ENDDO
      ELSEIF ( boun .eq. 3 .AND. ord .eq. 2 .AND. ALLOCATED(plan%y%neu) ) THEN
!$omp paralleldo if (iend-ista.ge.nth) private (j,k,tmp)
         DO i=ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,tmp)
            DO j=1,nz-Cz
               f(j,1,i) = plan%y%neu2(dy)*f(j,1,i)
               DO k=1,dy-1
                  f(j,1,i) = f(j,1,i) + plan%y%neu2(k)*f(j,dy+1-k,i)
               ENDDO
            ENDDO
         ENDDO
      ELSEIF ( boun .eq. 4 .AND. ord .eq. 1 .AND. ALLOCATED(plan%y%neu) ) THEN
!$omp paralleldo if (iend-ista.ge.nth) private (j,k,tmp)
         DO i=ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,tmp)
            DO j=1,nz-Cz
               f(j,ny-Cy,i) = plan%y%neu(dy)*f(j,ny-Cy,i)
               DO k=1,dy-1
                  f(j,ny-Cy,i) = f(j,ny-Cy,i) + plan%y%neu(k)*f(j,ny-Cy-dy+k,i)
               ENDDO
            ENDDO
         ENDDO
      ELSEIF ( boun .eq. 4 .AND. ord .eq. 2 .AND. ALLOCATED(plan%y%neu) ) THEN
!$omp paralleldo if (iend-ista.ge.nth) private (j,k,tmp)
         DO i=ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,tmp)
            DO j=1,nz-Cz
               f(j,ny-Cy,i) = plan%y%neu2(dy)*f(j,ny-Cy,i)
               DO k=1,dy-1
                  f(j,ny-Cy,i) = f(j,ny-Cy,i) + plan%y%neu2(k)*f(j,ny-Cy-dy+k,i)
               ENDDO
            ENDDO
         ENDDO
      ! ...........
      ! Z direction
      ! ...........
      ELSEIF ( boun .eq. 5 .AND. ord .eq. 1 .AND. ALLOCATED(plan%z%neu) ) THEN
!$omp paralleldo if (iend-ista.ge.nth) private (j,k,tmp)
         DO i=ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,tmp)
            DO j=1,ny
               f(1,j,i) = plan%z%neu(dz)*f(1,j,i)
               DO k=1,dz-1
                  f(1,j,i) = f(1,j,i) + plan%z%neu(k)*f(dz+1-k,j,i)
               ENDDO
            ENDDO
         ENDDO
      ELSEIF ( boun .eq. 5 .AND. ord .eq. 2 .AND. ALLOCATED(plan%z%neu2) ) THEN
!$omp paralleldo if (iend-ista.ge.nth) private (j,k,tmp)
         DO i=ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,tmp)
            DO j=1,ny
               f(1,j,i) = plan%z%neu2(dz)*f(1,j,i)
               DO k=1,dz-1
                  f(1,j,i) = f(1,j,i) + plan%z%neu2(k)*f(dz+1-k,j,i)
               ENDDO
            ENDDO
         ENDDO
      ELSEIF ( boun .eq. 6 .AND. ord .eq. 1 .AND. ALLOCATED(plan%z%neu) ) THEN
!$omp paralleldo if (iend-ista.ge.nth) private (j,k,tmp)
         DO i=ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,tmp)
            DO j=1,ny
               f(nz-Cz,j,i) = plan%z%neu(dz)*f(nz-Cz,j,i)
               DO k=1,dz-1
                  f(nz-Cz,j,i) = f(nz-Cz,j,i) + plan%z%neu(k)*f(nz-Cz-dz+k,j,i)
               ENDDO
            ENDDO
         ENDDO
      ELSEIF ( boun .eq. 6 .AND. ord .eq. 2 .AND. ALLOCATED(plan%z%neu2) ) THEN
!$omp paralleldo if (iend-ista.ge.nth) private (j,k,tmp)
         DO i=ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,tmp)
            DO j=1,ny
               f(nz-Cz,j,i) = plan%z%neu2(dz)*f(nz-Cz,j,i)
               DO k=1,dz-1
                  f(nz-Cz,j,i) = f(nz-Cz,j,i) + plan%z%neu2(k)*f(nz-Cz-dz+k,j,i)
               ENDDO
            ENDDO
         ENDDO
      ELSE
         IF ( myrank .eq. 0) PRINT*, "[ERROR] Neumann reconstruction not ", &
             "performed. Unsupported derivative order, boundary surface or ", &
             "unallocated table. Aborting..."
         CALL MPI_FINALIZE(MPI_COMM_WORLD, ierr)
         STOP
      ENDIF

      CALL GTStop(hneu); neutime = neutime + GTGetTime(hneu); 


      END SUBROUTINE neumann_reconstruct

!**********************************************************************
      SUBROUTINE robin_reconstruct(plan,f,boun,a)
!----------------------------------------------------------------------
! Reconstructs function values at the boundary from the Robin
! boundary condition f' + a*f = g. The array f must contain the
! values of g at the desired endpoint. This version of the 
! subroutine supports the case in which a is a (real) function
! of the space/wavenumber domain.
! ARGUMENTS:
!     plan : fcplan to use for the reconstruction.
!     a    : Array whose values must be reconstructed. [INOUT]
!     boun : Boundary to which the reconstruction should 
!            be performed. 1 for x=0, 2 for x=Lx ... 6 for z=Lz [IN]
!     a    : Array with the coefficients for the linear combination
!            between the derivative and the function. [IN]
!
! NOTE: Deferring shape of a allows to consider both the y and z
! cases (where a has shapes nz,ista:iend and ny,ista:iend,
! respectively) in the same subroutine.
! TODO: Create a separate subroutine and an interface for the
! constant a case.
!-----------------------------------------------------------------------
      USE commtypes
      USE fprecision
      USE mpivars
      USE gtimer
!$    USE threads

      IMPLICIT NONE

      TYPE(FCPLAN), INTENT(IN)                                                :: plan
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(plan%z%N,plan%y%N,ista:iend) :: f
      REAL(KIND=GP), INTENT(INOUT), DIMENSION(:,:)                            :: a
      INTEGER, INTENT(IN)                                                     :: boun

      INTEGER          :: Cy,Cz,dy,dz,ny,nz
      INTEGER          :: i,j,k

      CALL GTStart(hrob)

      ny = plan%y%N
      nz = plan%z%N

      Cy = plan%y%C
      Cz = plan%z%C
      dy = plan%y%d
      dz = plan%z%d

      IF (boun .eq. 3 .OR. boun .eq. 4) THEN
         IF ( .NOT. ALLOCATED(plan%y%neu) ) THEN
            IF (myrank.eq.0) PRINT*, "[ERROR] Neumann table not loaded in", &
                "robin_reconstruct. Aborting..."
            CALL MPI_FINALIZE(MPI_COMM_WORLD, ierr)
            STOP
         ENDIF
      ELSE IF (boun .eq. 5 .OR. boun .eq. 6) THEN
         IF ( .NOT. ALLOCATED(plan%z%neu) ) THEN
            IF (myrank.eq.0) PRINT*, "[ERROR] Neumann table not loaded in", &
                "robin_reconstruct. Aborting..."
            CALL MPI_FINALIZE(MPI_COMM_WORLD, ierr)
            STOP
         ENDIF
      ENDIF

      ! ...........
      ! Y direction
      ! ...........
      IF ( boun .eq. 3 ) THEN
!$omp paralleldo if (iend-ista.ge.nth) private (j,k,tmp)
         DO i=ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,tmp)
            DO j=1,nz-Cz
               f(j,1,i) = plan%y%neu(dy)*f(j,1,i)
               DO k=1,dy-1
                  f(j,1,i) = f(j,1,i) + plan%y%neu(k)*f(j,dy+1-k,i)
               ENDDO
               !ista+1 required due to deferred shape
               f(j,1,i) = f(j,1,i)/(a(j,i-ista+1)*plan%y%neu(dy)+1)
            ENDDO
         ENDDO
      ELSEIF ( boun .eq. 4 ) THEN
!$omp paralleldo if (iend-ista.ge.nth) private (j,k,tmp)
         DO i=ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,tmp)
            DO j=1,nz-Cz
               f(j,ny-Cy,i) = plan%y%neu(dy)*f(j,ny-Cy,i)
               DO k=1,dy-1
                  f(j,ny-Cy,i) = f(j,ny-Cy,i) + plan%y%neu(k)*f(j,ny-Cy-dy+k,i)
               ENDDO
               !ista+1 required due to deferred shape
               f(j,ny-Cy,i) = f(j,ny-Cy,i)/(a(j,i-ista+1)*plan%y%neu(dy)+1)  
            ENDDO
         ENDDO
      ! ...........
      ! Z direction
      ! ...........
      ELSEIF ( boun .eq. 5 ) THEN
!$omp paralleldo if (iend-ista.ge.nth) private (j,k,tmp)
         DO i=ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,tmp)
            DO j=1,ny
               f(1,j,i) = plan%z%neu(dz)*f(1,j,i)
               DO k=1,dz-1
                  f(1,j,i) = f(1,j,i) + plan%z%neu(k)*f(dz+1-k,j,i)
               ENDDO
               !ista+1 required due to deferred shape
               f(1,j,i) = f(1,j,i)/(a(j,i-ista+1)*plan%z%neu(dz)+1)
            ENDDO
         ENDDO
      ELSEIF ( boun .eq. 6 ) THEN
!$omp paralleldo if (iend-ista.ge.nth) private (j,k,tmp)
         DO i=ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,tmp)
            DO j=1,ny
               f(nz-Cz,j,i) = plan%z%neu(dz)*f(nz-Cz,j,i)
               DO k=1,dz-1
                  f(nz-Cz,j,i) = f(nz-Cz,j,i) + plan%z%neu(k)*f(nz-Cz-dz+k,j,i)
               ENDDO
               !ista+1 required due to deferred shape
               f(nz-Cz,j,i) = f(nz-Cz,j,i)/(a(j,i-ista+1)*plan%z%neu(dz)+1)  
            ENDDO
         ENDDO
      ELSE
         IF ( myrank .eq. 0) PRINT*, "[ERROR] Robin reconstruction not ", &
             "performed. Wrong boundary specified. Aborting..."
         CALL MPI_FINALIZE(MPI_COMM_WORLD, ierr)
         STOP
      ENDIF

      CALL GTStop(hrob); robtime = robtime + GTGetTime(hrob); 

      END SUBROUTINE robin_reconstruct

!***************************************************************************
      FUNCTION is_neumann_loaded(coord,ord) RESULT(output)
!---------------------------------------------------------------------------
!  Checks if the corresponding Neumann table of order ord has been loaded at
! the specified border.
!---------------------------------------------------------------------------
      USE commtypes
      USE mpivars

      IMPLICIT NONE
      TYPE(FCCOORD), INTENT(IN) :: coord
      INTEGER, INTENT(IN)       :: ord
      LOGICAL                   :: output
 
      IF (ord .eq. 1) THEN
         output = ALLOCATED(coord%neu)
      ELSE IF (ord .eq. 2) THEN
         output = ALLOCATED(coord%neu2)
      ELSE
         IF ( myrank .eq. 0) PRINT*, "[ERROR] Failed to check status ", &
         "of Neumann tables, unrecognized order. Aborting..."
         CALL MPI_FINALIZE(MPI_COMM_WORLD, ierr)
         STOP
      ENDIF

      END FUNCTION is_neumann_loaded

ENDMODULE fcgram
