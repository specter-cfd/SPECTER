CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) THEN
PRINT*, "-*-*-*-*-*-*-*-*-*-*-*-*  fc_dirichlet.f90  -*-*-*-*-*-*-*-*-*-*-*-*-*"
PRINT*, "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"
ENDIF

C1=0;C2=0;C3=0;C4=0;C5=0;C6=0
R1=0;R2=0;R3=0

! Target function
DO k=ksta,kend
DO j=1,ny
DO i=1,nx
   R1(i,j,k) = sin(4*x(i))*cos(8*y(j))*sin(6*z(k))
ENDDO
ENDDO
ENDDO

! Get continued field
CALL fftp3d_real_to_complex(planfc,R1,C1,MPI_COMM_WORLD)
! Numerical derivatives
CALL derivk(C1,C4,1); CALL derivk(C1,C5,2); CALL derivk(C1,C6,3)

! ...........
! X direction
! ...........
! Exact derivative
DO k=ksta,kend
DO j=1,ny
DO i=1,nx
   R2(i,j,k) = 4*cos(4*x(i))*cos(8*y(j))*sin(6*z(k))
ENDDO
ENDDO
ENDDO
! Print diagnostic
CALL fftp3d_complex_to_real(planfc,C4,R3,MPI_COMM_WORLD)
R3 = R3/nx/ny/nz
rmp = MAXVAL(ABS(R3(1:nx-Cx,1:ny-Cy,ksta:pkend)-R2(1:nx-Cx,1:ny-Cy,ksta:pkend)))
CALL MPI_REDUCE(rmp,rmq,1,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
rmp = SUM(ABS(R3(1:nx-Cx,1:ny-Cy,ksta:pkend)-R2(1:nx-Cx,1:ny-Cy,ksta:pkend)))
CALL MPI_REDUCE(rmp,rms,1,GC_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) PRINT*, "Error in x derivative:"
IF (myrank .eq. 0) PRINT*, "Maximum: ", rmq, &
    "    Average: ", rms/(nx-Cx)/(ny-Cy)/(nz-Cz)

! ...........
! Y direction
! ...........
! Exact derivative
DO k=ksta,pkend
DO j=1,ny
DO i=1,nx
   R2(i,j,k) = sin(4*x(i))*(-8*sin(8*y(j)))*sin(6*z(k))
ENDDO
ENDDO
ENDDO
! Print diagnostic
CALL fftp3d_complex_to_real(planfc,C5,R3,MPI_COMM_WORLD)
R3 = R3/nx/ny/nz
rmp = MAXVAL(ABS(R3(1:nx-Cx,1:ny-Cy,ksta:pkend)-R2(1:nx-Cx,1:ny-Cy,ksta:pkend)))
CALL MPI_REDUCE(rmp,rmq,1,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
rmp = SUM(ABS(R3(1:nx-Cx,1:ny-Cy,ksta:pkend)-R2(1:nx-Cx,1:ny-Cy,ksta:pkend)))
CALL MPI_REDUCE(rmp,rms,1,GC_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) PRINT*, "Error in y derivative:"
IF (myrank .eq. 0) PRINT*, "Maximum: ", rmq, &
    "    Average: ", rms/(nx-Cx)/(ny-Cy)/(nz-Cz)

! ...........
! Z direction
! ...........
! Exact derivative
DO k=ksta,pkend
DO j=1,ny
DO i=1,nx
   R2(i,j,k) = sin(4*x(i))*cos(8*y(j))*6*cos(6*z(k))
ENDDO
ENDDO
ENDDO

! Print diagnostic
CALL fftp3d_complex_to_real(planfc,C6,R3,MPI_COMM_WORLD)
R3 = R3/nx/ny/nz
rmp = MAXVAL(ABS(R3(1:nx-Cx,1:ny-Cy,ksta:pkend)-R2(1:nx-Cx,1:ny-Cy,ksta:pkend)))
CALL MPI_REDUCE(rmp,rmq,1,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
rmp = SUM(ABS(R3(1:nx-Cx,1:ny-Cy,ksta:pkend)-R2(1:nx-Cx,1:ny-Cy,ksta:pkend)))
CALL MPI_REDUCE(rmp,rms,1,GC_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) PRINT*, "Error in z derivative:"
IF (myrank .eq. 0) PRINT*, "Maximum: ", rmq, &
    "    Average: ", rms/(nx-Cx)/(ny-Cy)/(nz-Cz)

IF (myrank .eq. 0) &
PRINT*, "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"
