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
R2 = R1
R3 = R1

tmp = SUM(DBLE(R1**2 + R2**2 + R3**2)/nx/ny/(nz-Cz))

CALL MPI_REDUCE(tmp,tmq,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

CALL fftp3d_real_to_complex(planrc,R1,C1,MPI_COMM_WORLD)
CALL fftp3d_real_to_complex(planrc,R2,C2,MPI_COMM_WORLD)
CALL fftp3d_real_to_complex(planrc,R3,C3,MPI_COMM_WORLD)

CALL energy(C1,C2,C3,tmp,1)

IF ( myrank .eq. 0 ) PRINT*, "Relative difference between mean energy ", &
        "in real and kx,ky,z domains is:", char(10), &
        ABS(tmq-tmp)/ABS(tmq)

IF ( myrank .eq. 0 ) PRINT*, tmq, tmp
