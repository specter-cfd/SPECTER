! TODO Add tests for x direction when properly
! implemented in the fcgram module.
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) THEN
PRINT*, "-*-*-*-*-*-*-*-*-*-*-*-*-*  fc_robin.f90  -*-*-*-*-*-*-*-*-*-*-*-*-*-*"
PRINT*, "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"
ENDIF
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

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

! Exact derivative in y direction
DO k=ksta,pkend
DO j=1,ny
DO i=1,nx
   R2(i,j,k) = -8*sin(4*x(i))*sin(8*y(j))*sin(6*z(k))
ENDDO
ENDDO
ENDDO

! Exact derivative in z direction
DO k=ksta,pkend
DO j=1,ny
DO i=1,nx
   R3(i,j,k) = 6*sin(4*x(i))*cos(8*y(j))*cos(6*z(k))
ENDDO
ENDDO
ENDDO

IF ( Cy .gt. 0 ) THEN
   ! Field and derivatives in (z,y,kx) domain
   CALL fftp1d_real_to_complex_x(planfc,R1,C1,MPI_COMM_WORLD)
   CALL fftp1d_real_to_complex_x(planfc,R2,C2,MPI_COMM_WORLD)
   CALL fftp1d_real_to_complex_x(planfc,R3,C3,MPI_COMM_WORLD)

   ! Generate sqrt(kx**2+kz**2)
   DO i=ista,iend
      DO k=ksta,kend
         R2(k,1,i) = sqrt(kx(i)**2 + kz(k)**2)
      ENDDO
   ENDDO

   ! Populate Robin BC at endpoints
   C1(:,1,:)     = -C2(:,1,:) + R2(:,1,:)*C1(:,1,:) ! d/dy = -d/dn
   C1(:,ny-Cy,:) = C2(:,ny-Cy,:) + R2(:,1,:)*C1(:,ny-Cy,:)

   ! Reconstruct Y boundaries
   IF (.NOT. ALLOCATED(planfc%y%neu)) THEN
       CALL load_neumann_tables(planfc%y,planfc%tdir,1)
   ENDIF
   CALL robin_reconstruct(planfc,C1,3,R2(:,1,:))
   CALL robin_reconstruct(planfc,C1,4,R2(:,1,:))
ELSE
   CALL fftp2d_real_to_complex_xy(planfc,R1,C1,MPI_COMM_WORLD)
   CALL fftp2d_real_to_complex_xy(planfc,R3,C3,MPI_COMM_WORLD)
ENDIF

! Populate Robin BC at endpoints
C1(1,:,:)     = -C3(1,:,:) + khom(:,:)*C1(1,:,:) ! d/dz = -d/dn
C1(nz-Cz,:,:) = C3(nz-Cz,:,:) + khom(:,:)*C1(nz-Cz,:,:)
! Load tables and reconstruct
IF (.NOT. ALLOCATED(planfc%z%neu)) THEN
    CALL load_neumann_tables(planfc%z,planfc%tdir,1)
ENDIF
CALL robin_reconstruct(planfc,C1,5,khom)
CALL robin_reconstruct(planfc,C1,6,khom)

! -----------------------------
! Print errors in real space
! -----------------------------
IF ( Cy .gt. 0 ) THEN
   C4 = C1
   CALL fftp1d_complex_to_real_x(planfc,C4,R3,MPI_COMM_WORLD)
   R3 = R3/nx
   rmp = MAXVAL(ABS(R3(1:nx-Cx,1,:)-R1(1:nx-Cx,1,:)))
   CALL MPI_REDUCE(rmp,rmq,1,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
   rmp = MAXVAL(ABS(R3(1:nx-Cx,ny-Cy,:)-R1(1:nx-Cx,ny-Cy,:)))
   CALL MPI_REDUCE(rmp,rms,1,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
   IF ( myrank .eq. 0 ) THEN
      PRINT*, "Maximum reconstruction error at y=0 : ", rmq
      PRINT*, "Maximum reconstruction error at y=Ly: ", rms
   ENDIF
ELSE
   C4 = C1
   CALL fftp2d_complex_to_real_xy(planfc,C4,R3,MPI_COMM_WORLD)
   R3 = R3/nx/ny
ENDIF
IF ( myrank .eq. 0 ) THEN
   rm1 = MAXVAL(ABS(R2(:,1:ny-Cy,1)-R1(:,1:ny-Cy,1)))
   PRINT*, "Maximum reconstruction error at z=0 : ", rm1
ENDIF
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
IF ( ksta .le. nz-Cz .AND. nz-Cz .le. kend ) THEN
   rm2 = MAXVAL(ABS(R3(:,1:ny-Cy,nz-Cz)-R1(:,1:ny-Cy,nz-Cz)))
   PRINT*, "Maximum reconstruction error at z=Lz: ", rm2
ENDIF
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

! -------------------------------------------
! Print derivative errors in real space
! -------------------------------------------

! Get exact d/dy again
! Exact derivative in y direction
DO k=ksta,pkend
DO j=1,ny
DO i=1,nx
   R2(i,j,k) = -8*sin(4*x(i))*sin(8*y(j))*sin(6*z(k))
ENDDO
ENDDO
ENDDO

IF ( Cy .gt. 0 ) THEN
   ! Get analytical z-derivative again 
   CALL fftp1d_complex_to_real_x(planfc,C3,R3,MPI_COMM_WORLD)
   R3 = R3/nx

   CALL fftp2d_real_to_complex_yz(planfc,C1,MPI_COMM_WORLD)
   CALL derivk(C1,C2,2)

   ! Print diagnostic
   CALL fftp3d_complex_to_real(planfc,C2,R1,MPI_COMM_WORLD)
   R1 = R1/nx/ny/nz
   rmp = MAXVAL(ABS(R1(1:nx-Cx,1:ny-Cy,ksta:pkend)-R2(1:nx-Cx,1:ny-Cy,ksta:pkend)))
   CALL MPI_REDUCE(rmp,rmq,1,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
   rmp = SUM(ABS(R1(1:nx-Cx,1:ny-Cy,ksta:pkend)-R2(1:nx-Cx,1:ny-Cy,ksta:pkend)))
   CALL MPI_REDUCE(rmp,rms,1,GC_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)

   CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
   IF (myrank .eq. 0) PRINT*, "Error in y derivative:"
   IF (myrank .eq. 0) PRINT*, "Maximum: ", rmq, &
       "    Average: ", rms/(nx-Cx)/(ny-Cy)/(nz-Cz)
   CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

   CALL derivk(C1,C3,3) ! 3D Fourier z-derivative
ELSE
   ! Get analytical z-derivative again 
   CALL fftp2d_complex_to_real_xy(planfc,C3,R3,MPI_COMM_WORLD)
   R3 = R3/ny/ny

   CALL fftp1d_real_to_complex_z(planfc,C1,MPI_COMM_WORLD)
   CALL derivk(C1,C3,3) ! 3D Fourier z-derivative
ENDIF

! Print diagnostic
CALL fftp3d_complex_to_real(planfc,C3,R1,MPI_COMM_WORLD)
R1 = R1/nx/ny/nz
rmp = MAXVAL(ABS(R1(1:nx-Cx,1:ny-Cy,ksta:pkend)-R3(1:nx-Cx,1:ny-Cy,ksta:pkend)))
CALL MPI_REDUCE(rmp,rmq,1,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
rmp = SUM(ABS(R1(1:nx-Cx,1:ny-Cy,ksta:pkend)-R3(1:nx-Cx,1:ny-Cy,ksta:pkend)))
CALL MPI_REDUCE(rmp,rms,1,GC_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) PRINT*, "Error in z derivative:"
IF (myrank .eq. 0) PRINT*, "Maximum: ", rmq, &
    "    Average: ", rms/(nx-Cx)/(ny-Cy)/(nz-Cz)
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

IF (myrank .eq. 0) &
PRINT*, "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"
