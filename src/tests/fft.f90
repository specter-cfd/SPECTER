CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) THEN
PRINT*, "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"
PRINT*, "-*-*-*-*-*-*-*-*-*-*-*-*-*-  fft.f90 -* -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"
PRINT*, "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"
ENDIF

C1=0;C2=0;C3=0;C4=0;C5=0;C6=0
R1=0;R2=0;R3=0

! Target function
DO k=ksta,pkend
DO j=1,ny-Cy
DO i=1,nx
   R1(i,j,k) = sin(4*x(i))*cos(8*y(j))*sin(6*z(k))
ENDDO
ENDDO
ENDDO

CALL fftp2d_real_to_complex_xy(planfc,R1,C1,MPI_COMM_WORLD)
C3 = C1
CALL fftp1d_real_to_complex_z(planfc,C3,MPI_COMM_WORLD)

CALL fftp1d_real_to_complex_x(planfc,R1,C2,MPI_COMM_WORLD)
C4 = C2
CALL fftp2d_real_to_complex_yz(planfc,C4,MPI_COMM_WORLD)

rmp = MAXVAL(ABS(C3-C4))
CALL MPI_REDUCE(rmp,rmq,1,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) PRINT*, "Maximum difference between XY->Z transform ", &
        "and X-> YZ transform:", char(10), rmq

rmp = MAXVAL(ABS(C3))
CALL MPI_REDUCE(rmp,rmq,1,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
rmp = MAXVAL(ABS(C4))
CALL MPI_REDUCE(rmp,rms,1,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)

IF (Cx.eq.0 .AND. Cy.eq.0 .AND. Cz.eq.0) THEN
IF(iend.ge.5 .AND. ista.le.5) PRINT*, "Abs. value expected vs. found at: ", &
    "k=(4,8,6)", nx*ny*nz/8, "[XY->Z]",ABS(C3(7,9,5)), "[X->YZ]",ABS(C4(7,9,5))
IF (myrank .eq. 0) PRINT*, "Maximum value found: [XY->Z]", rmq, &
    "[X->YZ]", rms
ENDIF

CALL fftp1d_complex_to_real_z(planfc,C3,MPI_COMM_WORLD)
rmp = MAXVAL(ABS(C1(1:nz-Cz,:,:)-C3(1:nz-Cz,:,:)/nz))
CALL MPI_REDUCE(rmp,rmq,1,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) PRINT*, "Maximum difference between Z transform ", &
        "and Z anti-transform:", char(10), rmq

CALL fftp2d_complex_to_real_xy(planfc,C1,R2,MPI_COMM_WORLD)
rmp = MAXVAL(ABS(R1(:,:,ksta:pkend)-R2(:,:,ksta:pkend)/nx/ny))
CALL MPI_REDUCE(rmp,rmq,1,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) PRINT*, "Maximum difference between XY transform ", &
        "and XY anti-transform:", char(10), rmq

CALL fftp2d_complex_to_real_yz(planfc,C4,MPI_COMM_WORLD)
rmp = MAXVAL(ABS(C2(1:nz-Cz,1:ny-Cy,:)-C4(1:nz-Cz,1:ny-Cy,:)/nz/ny))
CALL MPI_REDUCE(rmp,rmq,1,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) PRINT*, "Maximum difference between YZ transform ", &
        "and YZ anti-transform:", char(10), rmq
!
CALL fftp1d_complex_to_real_x(planfc,C2,R2,MPI_COMM_WORLD)
rmp = MAXVAL(ABS(R1(:,:,ksta:pkend)-R2(:,:,ksta:pkend)/nx))
CALL MPI_REDUCE(rmp,rmq,1,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) PRINT*, "Maximum difference between X transform ", &
        "and X anti-transform:", char(10), rmq

IF (myrank .eq. 0) &
PRINT*, "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"
