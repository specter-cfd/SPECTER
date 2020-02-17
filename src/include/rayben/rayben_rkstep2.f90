! Step 2 of Runge-Kutta for the HD equations
! Computes the nonlinear terms and evolves the equations in dt/o
         rmp = 1/real(o, kind=GP)

         ! Non-linear terms
         CALL gradre(vx,vy,vz,C4,C5,C6)
         CALL advect(vx,vy,vz,th,C8)

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  C6(k,j,i) = C6(k,j,i) - xmom*th(k,j,i)  ! Bouyancy
                  C8(k,j,i) = C8(k,j,i) - xtemp*vz(k,j,i) ! Heat current
               ENDDO
            ENDDO
         ENDDO

         ! Dealias non-linear 
         CALL fc_filter(C4)
         CALL fc_filter(C5)
         CALL fc_filter(C6)
         CALL fc_filter(C8)

         ! Laplacian
         CALL laplak(vx,vx)
         CALL laplak(vy,vy)
         CALL laplak(vz,vz)
         CALL laplak(th,th)

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
         DO k = 1,nz
            ! Perform half oth-RK step
            vx(k,j,i) = C1(k,j,i) + dt*(nu*vx(k,j,i)-C4(k,j,i)+&
                fx(k,j,i))*rmp
            vy(k,j,i) = C2(k,j,i) + dt*(nu*vy(k,j,i)-C5(k,j,i)+&
                fy(k,j,i))*rmp
            vz(k,j,i) = C3(k,j,i) + dt*(nu*vz(k,j,i)-C6(k,j,i)+&
                fz(k,j,i))*rmp
            th(k,j,i) = C7(k,j,i) + dt*(kappa*th(k,j,i)-C8(k,j,i)+&
                fs(k,j,i))*rmp
         ENDDO
         ENDDO
         ENDDO

         ! Apply boundary conditions for intermediate velocity field
         CALL fftp1d_complex_to_real_z(plancr,vx,MPI_COMM_WORLD)
         CALL fftp1d_complex_to_real_z(plancr,vy,MPI_COMM_WORLD)
         CALL fftp1d_complex_to_real_z(plancr,vz,MPI_COMM_WORLD)
         CALL fftp1d_complex_to_real_z(plancr,th,MPI_COMM_WORLD)

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
               vz(1,j,i) = vz(1,j,i)*rmq
               th(1,j,i) = 0.0_GP 
            ENDDO
         ENDDO
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 2,nz-Cz-1
                  vx(k,j,i) = vx(k,j,i)*rmq
                  vy(k,j,i) = vy(k,j,i)*rmq
                  vz(k,j,i) = vz(k,j,i)*rmq
                  th(k,j,i) = th(k,j,i)*rmq
               ENDDO
            ENDDO
         ENDDO
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               vx(nz-Cz,j,i) = im*kx(i)*pr(nz-Cz,j,i)*rmp
               vy(nz-Cz,j,i) = im*ky(j)*pr(nz-Cz,j,i)*rmp
               vz(nz-Cz,j,i) = vz(nz-Cz,j,i)*rmq
               th(nz-Cz,j,i) = 0.0_GP
            ENDDO
         ENDDO

         ! Add slip and mean temperature perturbation at the boundaries
         IF ( ista .eq. 1) vx(1,1,1)    =nx*ny*vxbot
         IF ( ista .eq. 1) vx(nz-Cz,1,1)=nx*ny*vxtop
         IF ( ista .eq. 1) vy(1,1,1)    =nx*ny*vybot
         IF ( ista .eq. 1) vy(nz-Cz,1,1)=nx*ny*vytop
         IF ( ista .eq. 1) th(1,1,1)    =nx*ny*thbot
         IF ( ista .eq. 1) th(nz-Cz,1,1)=nx*ny*thtop

         ! Return to (kz,ky,kx) domain
         CALL fftp1d_real_to_complex_z(planrc,vx,MPI_COMM_WORLD)
         CALL fftp1d_real_to_complex_z(planrc,vy,MPI_COMM_WORLD)
         CALL fftp1d_real_to_complex_z(planrc,vz,MPI_COMM_WORLD)
         CALL fftp1d_real_to_complex_z(planrc,th,MPI_COMM_WORLD)

         CALL fc_filter(vx); CALL fc_filter(vy); CALL fc_filter(vz)
         CALL fc_filter(th)

         ! Project into solenoidal space with null flux
         ! (second half of rkstep)
         CALL pressure(vx,vy,vz,pr)
