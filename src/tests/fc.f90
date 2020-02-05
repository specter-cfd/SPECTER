C1=0;C2=0;C3=0;C4=0;C5=0;C6=0
R1=0;R2=0;R3=0

! Target function
DO k=ksta,pkend
DO j=1,ny
DO i=1,nx
   R1(i,j,k) = sin(4*x(i))*cos(8*y(j))*sin(6*z(k))
ENDDO
ENDDO
ENDDO

! -----------------------------------------
! Test FFTP
! -----------------------------------------
CALL fftp3d_real_to_complex(planrc,R1,C1,MPI_COMM_WORLD)
CALL fftp2d_real_to_complex_xy(planrc,R1,C2,MPI_COMM_WORLD)
C3 = C2
CALL fftp1d_real_to_complex_z(planrc,C3,MPI_COMM_WORLD)

rmp = MAXVAL(ABS(C3-C1))
CALL MPI_REDUCE(rmp,rmq,1,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) PRINT*,  "Maximum difference between 3D transform",&
       " and 2D+1D transform:", char(10), rmq

C4 = C1
CALL fftp3d_complex_to_real(plancr,C4,R2,MPI_COMM_WORLD)
C5 = C3
CALL fftp1d_complex_to_real_z(plancr,C5,MPI_COMM_WORLD)
C6 = C5
CALL fftp2d_complex_to_real_xy(plancr,C6,R3,MPI_COMM_WORLD)

rmp = MAXVAL(ABS(R3-R2))
CALL MPI_REDUCE(rmp,rmq,1,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) PRINT*, "Maximum difference between 3D ",&
        "antitransform and 1D+2D antitransform:", char(10), rmq
rmp = MAXVAL(ABS(C2(1:nz-Cz,:,:)-C5(1:nz-Cz,:,:)/nz))
CALL MPI_REDUCE(rmp,rmq,1,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) PRINT*, "Maximum difference between 3D ",&
        "transform + 1D antitransform vs 2D transform:", char(10), rmq

R2 = R2/nx/ny/nz
rmp = MAXVAL(ABS(R1(:,:,ksta:pkend)-R2(:,:,ksta:pkend)))
CALL MPI_REDUCE(rmp,rmq,1,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) PRINT*, "Maximum error in original data vs 3D ", &
        "transform + 3D antitransform:", char(10), rmq


! -----------------------------------------
! Diagnostic derivatives
! -----------------------------------------
! Numerical derivatives
CALL derivk3(C1,C4,1); CALL derivk3(C1,C5,2); CALL derivk3(C1,C6,3)

! ...........
! X direction
! ...........
! Exact derivative
DO k=ksta,pkend
DO j=1,ny
DO i=1,nx
   R1(i,j,k) = 4*cos(4*x(i))*cos(8*y(j))*sin(6*z(k))
ENDDO
ENDDO
ENDDO
! Print diagnostic
IF (myrank .eq. 0) PRINT*, CHAR(10), "-----------------------", CHAR(10)
CALL fftp3d_complex_to_real(plancr,C4,R2,MPI_COMM_WORLD)
R2 = R2/nx/ny/nz
rmp = MAXVAL(ABS(R2(:,:,ksta:pkend)-R1(:,:,ksta:pkend)))
CALL MPI_REDUCE(rmp,rmq,1,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) PRINT*, "Maximum error in x derivative:", char(10), rmq
! Averages
rmp = SUM(ABS(R2(:,:,ksta:pkend)-R1(:,:,ksta:pkend)))
CALL MPI_REDUCE(rmp,rmq,1,GC_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) PRINT*, "Avg. error in x derivative:", char(10),&
    rmq/nx/ny/(nz-Cz)

! ...........
! Y direction
! ...........
! Exact derivative
DO k=ksta,pkend
DO j=1,ny
DO i=1,nx
   R1(i,j,k) = sin(4*x(i))*(-8*sin(8*y(j)))*sin(6*z(k))
ENDDO
ENDDO
ENDDO
! Print diagnostic
IF (myrank .eq. 0) PRINT*, CHAR(10), "-----------------------", CHAR(10)
CALL fftp3d_complex_to_real(plancr,C5,R2,MPI_COMM_WORLD)
R2 = R2/nx/ny/nz
rmp = MAXVAL(ABS(R2(:,:,ksta:pkend)-R1(:,:,ksta:pkend)))
CALL MPI_REDUCE(rmp,rmq,1,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) PRINT*, "Maximum error in y derivative:", char(10), rmq
! Averages
rmp = SUM(ABS(R2(:,:,ksta:pkend)-R1(:,:,ksta:pkend)))
CALL MPI_REDUCE(rmp,rmq,1,GC_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) PRINT*, "Avg. error in y derivative:", char(10),&
    rmq/nx/ny/(nz-Cz)

! ...........
! Z direction
! ...........
! Exact derivative
DO k=ksta,pkend
DO j=1,ny
DO i=1,nx
   R1(i,j,k) = sin(4*x(i))*cos(8*y(j))*6*cos(6*z(k))
ENDDO
ENDDO
ENDDO
! Print diagnostic
IF (myrank .eq. 0) PRINT*, CHAR(10), "-----------------------", CHAR(10)
CALL fftp3d_complex_to_real(plancr,C6,R2,MPI_COMM_WORLD)
R2 = R2/nx/ny/nz
rmp = MAXVAL(ABS(R2(:,:,ksta:pkend)-R1(:,:,ksta:pkend)))
CALL MPI_REDUCE(rmp,rmq,1,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) PRINT*, "Maximum error in z derivative:", char(10), rmq
! Averages
rmp = SUM(ABS(R2(:,:,ksta:pkend)-R1(:,:,ksta:pkend)))
CALL MPI_REDUCE(rmp,rmq,1,GC_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) PRINT*, "Avg. error in z derivative:", char(10),&
    rmq/nx/ny/(nz-Cz)
