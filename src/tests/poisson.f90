!TODO Add proper tests when laplace_xyz and laplace_yz are implemented.

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) THEN
PRINT*, "-*-*-*-*-*-*-*-*-*-*-*-*-*- poisson.f90 *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-"
PRINT*, "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"
ENDIF
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

vx=0;vy=0;vz=0;pr=0
C1=0;C2=0;C3=0;C4=0;C5=0;C6=0
R1=0;R2=0;R3=0

! Target field
! Constant are added to test the zero modes of Poisson
DO k=ksta,pkend
DO j=1,ny
DO i=1,nx
   R1(i,j,k) = sin(4*x(i))*cos(8*y(j))*exp(.4*z(k)/Lz) !+ 1000
   R2(i,j,k) = sin(4*x(i))*cos(8*y(j))*exp(.4*z(k)/Lz) !+ 2000
   R3(i,j,k) = sin(4*x(i))*cos(8*y(j))*exp(.4*z(k)/Lz) !+ 3000
ENDDO
ENDDO
ENDDO

! ----------------------------------------------
! Numeric divergence
! -----------------------------------------------
CALL fftp3d_real_to_complex(planfc,R1,C1,MPI_COMM_WORLD)
CALL fftp3d_real_to_complex(planfc,R2,C2,MPI_COMM_WORLD)
CALL fftp3d_real_to_complex(planfc,R3,C3,MPI_COMM_WORLD)

! -----------------------------------------
! Exact divergence
! -----------------------------------------
DO k=ksta,pkend
DO j=1,ny
DO i=1,nx
   R1(i,j,k) = exp(.4*z(k)/Lz)*(4*cos(4*x(i))*cos(8*y(j)) - 8*sin(4*x(i))*sin(8*y(j)) &
           + .4*sin(4*x(i))*cos(8*y(j))/Lz)
ENDDO
ENDDO
ENDDO
CALL fftp3d_real_to_complex(planfc,R1,C4,MPI_COMM_WORLD)

! ----------------------------------------------
! Compare
! ----------------------------------------------
CALL energy(C4,C4,C4,tmp,1)
tmp = tmp/dble(3.0)

CALL divergence(C1,C2,C3,tmq)

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) PRINT*, "Mean squared divergence (analytical): ", &
        char(10), tmp

IF (myrank .eq. 0) PRINT*, "Relative difference between exact and ", & 
        "numerical mean divergence is:", char(10), ABS(tmp-tmq)/tmp
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
! -------------------------------------------------------------------
! Now, test Poisson for velocity field (homogeneous condition on vz)
! -------------------------------------------------------------------
IF (myrank .eq. 0) PRINT*, "--------------------------"
IF (myrank .eq. 0) PRINT*, "Checking pressure solver: "
IF (myrank .eq. 0) PRINT*, "--------------------------"
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

C4=C1; C5=C2; C6=C3   ! Save the field
CALL sol_project(C1,C2,C3,pr,1,0,0)

CALL divergence(C1,C2,C3,tmq)
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) PRINT*, "Mean squared divergence ", &
  "after Poisson:", char(10), tmq
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

CALL bouncheck_z(planfc,tmp,tmq,C1,C2)
CALL bouncheck_z(planfc,tmr,tms,C3)

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) PRINT*, "Mean squared normal velocity at borders:",&
        char(10), tmr, tms, "(after Laplace)"
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) PRINT*, "--------------------------"
IF (myrank .eq. 0) PRINT*, "Checking conductor solver:"
IF (myrank .eq. 0) PRINT*, "--------------------------"
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

C1=C4; C2=C5; C3=C6   ! Restore the field
CALL sol_project(C1,C2,C3,pr,0,0,0)

CALL divergence(C1,C2,C3,tmq)
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) PRINT*, "Mean squared divergence ", &
  "after Poisson:", char(10), tmq
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

CALL fftp1d_real_to_complex_z(planfc,pr,MPI_COMM_WORLD)

CALL bouncheck_z(planfc,tmr,tms,pr)
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) PRINT*, "Mean squared electrostatic potential at borders:",&
        char(10), tmr/3, tms/3, "(after Laplace)"
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) PRINT*, "--------------------------"
IF (myrank .eq. 0) PRINT*, "Checking insulator solver:"
IF (myrank .eq. 0) PRINT*, "--------------------------"
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

C1=C4; C2=C5; C3=C6   ! Restore the field
CALL sol_project(C1,C2,C3,pr,0,2,2)

CALL divergence(C1,C2,C3,tmq)
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) PRINT*, "Mean squared divergence ", &
  "after Poisson:", char(10), tmq
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
CALL fftp1d_real_to_complex_z(planfc,pr,MPI_COMM_WORLD)

CALL robcheck(planfc,pr,pr,pr,tmp,tmq,tmr,tms)
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) PRINT*, "Mean squared error in Robin boundary for the ", &
        "electrostatic potential:", char(10), tmr/3, tms/3, "(after Laplace)"
    PRINT*, "(This error is higher than the others cases because it involves ",&
        "estimating a derivative with FC..."
    PRINT*, "should be of the FC order for this resolution)"
IF (myrank .eq. 0) &
PRINT*, "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

