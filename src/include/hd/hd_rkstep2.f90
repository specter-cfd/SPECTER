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

         ! Apply boundary conditions and project onto solenoidal space
         CALL v_imposebc_and_project(vplanbc,planfc,vx,vy,vz,pr,rki=o,
                 v_zsta=(/ vxzsta, vysta/), vz_end=(/ vxzend, vyzend/))
