C1=0;C2=0;C3=0;C4=0;C5=0;C6=0
R1=0;R2=0;R3=0

! Target field
DO k=ksta,pkend
DO j=1,ny
DO i=1,nx
   R1(i,j,k) = sin(4*x(i))*cos(8*y(j))*exp(.4*z(k)/Lz)
   R2(i,j,k) = sin(4*x(i))*cos(8*y(j))*exp(.4*z(k)/Lz)
   R3(i,j,k) = sin(4*x(i))*cos(8*y(j))*exp(.4*z(k)/Lz)
ENDDO
ENDDO
ENDDO

! ----------------------------------------------
! Numeric divergence
! -----------------------------------------------
CALL fftp3d_real_to_complex(planrc,R1,C1,MPI_COMM_WORLD)
CALL fftp3d_real_to_complex(planrc,R2,C2,MPI_COMM_WORLD)
CALL fftp3d_real_to_complex(planrc,R3,C3,MPI_COMM_WORLD)

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
CALL fftp3d_real_to_complex(planrc,R1,C4,MPI_COMM_WORLD)

! ----------------------------------------------
! Compare
! ----------------------------------------------
CALL energy(C4,C4,C4,tmp,1)
tmp = tmp/dble(3.0)

CALL divergence(C1,C2,C3,tmq)

IF (myrank .eq. 0) PRINT*, "Mean squared divergence (analytical): ", &
        char(10), tmp

IF (myrank .eq. 0) PRINT*, "Relative difference between exact and ", & 
        "numerical mean divergence is:", char(10), ABS(tmp-tmq)/tmp

! ----------------------------------------------
! Now, test Poisson
! ----------------------------------------------

CALL pressure(C1,C2,C3,C4,C5,C6)

CALL energy(C1,C2,C3,tmp,1)
CALL divergence(C1,C2,C3,tmq)
IF (myrank .eq. 0) PRINT*, "Mean squared divergence ", &
  "after Poisson:", char(10), tmq!/tmp, "   (energy = ", tmp, ")"

! ----------------------------------------------
! Check flux condition
! ----------------------------------------------
CALL bouncheck(C1,C2,C3,tmp,tmq,tmr,tms)
!IF (myrank .eq. 0) PRINT*, "Mean squared slip velocity at borders:",&
!        char(10), rmp, rmq
IF (myrank .eq. 0) PRINT*, "Mean squared normal velocity at borders:",&
        char(10), tmr, tms, "(after Laplace)"
