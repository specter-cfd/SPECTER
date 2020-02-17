! Step 2 of Runge-Kutta for the HD equations
! Computes the nonlinear terms and evolves the equations in dt/o
         rmp = 1/real(o, kind=GP)

         ! Non-linear term
         CALL gradre(vx,vy,vz,C4,C5,C6)

         ! Dealias non-linear 
         CALL fc_filter(C4)
         CALL fc_filter(C5)
         CALL fc_filter(C6)

         ! Laplacian
         CALL laplak(vx,vx)
         CALL laplak(vy,vy)
         CALL laplak(vz,vz)

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
         DO k = 1,nz
            ! Perform oth-RK step
            vx(k,j,i) = C1(k,j,i) + dt*(nu*vx(k,j,i)-C4(k,j,i)+&
                fx(k,j,i))*rmp
            vy(k,j,i) = C2(k,j,i) + dt*(nu*vy(k,j,i)-C5(k,j,i)+&
                fy(k,j,i))*rmp
            vz(k,j,i) = C3(k,j,i) + dt*(nu*vz(k,j,i)-C6(k,j,i)+&
                fz(k,j,i))*rmp
         ENDDO
         ENDDO
         ENDDO

         ! Apply boundary conditions for intermediate velocity field
         CALL fftp1d_complex_to_real_z(plancr,vx,MPI_COMM_WORLD)
         CALL fftp1d_complex_to_real_z(plancr,vy,MPI_COMM_WORLD)

         IF (o .ne. ord) rmp = real(o+1, kind=GP)*rmp
         rmq = 1/real(nz, kind=GP)

         ! Normalize IFFT and apply boundary conditions at the same time
         ! BC are only for tangential velocity, normal velocity
         ! is enforced via the pressure. Details in: !TODO paper
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               vx(1,j,i) = im*kx(i)*pr(1,j,i)*rmp
               vy(1,j,i) = im*ky(j)*pr(1,j,i)*rmp
            ENDDO
         ENDDO
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 2,nz-Cz-1
                  vx(k,j,i) = vx(k,j,i)*rmq
                  vy(k,j,i) = vy(k,j,i)*rmq
               ENDDO
            ENDDO
         ENDDO
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               vx(nz-Cz,j,i) = im*kx(i)*pr(nz-Cz,j,i)*rmp
               vy(nz-Cz,j,i) = im*ky(j)*pr(nz-Cz,j,i)*rmp
            ENDDO
         ENDDO

         ! Add slip of the boundaries
         IF ( ista .eq. 1) vx(1,1,1)    =nx*ny*vxbot
         IF ( ista .eq. 1) vx(nz-Cz,1,1)=nx*ny*vxtop
         IF ( ista .eq. 1) vy(1,1,1)    =nx*ny*vybot
         IF ( ista .eq. 1) vy(nz-Cz,1,1)=nx*ny*vytop


         ! Return to (kz,ky,kx) domain
         CALL fftp1d_real_to_complex_z(planrc,vx,MPI_COMM_WORLD)
         CALL fftp1d_real_to_complex_z(planrc,vy,MPI_COMM_WORLD)
         
         ! Project into solenoidal space with null flux
         CALL pressure(vx,vy,vz,pr)
