! Step 2 of Runge-Kutta for the HD equations
! Computes the nonlinear terms and evolves the equations in dt/o
         rmp = 1/real(o, kind=GP)

         ! Magnetic field
         CALL curlk(ay,az,C12,1)
         CALL curlk(ax,az,C13,2)
         CALL curlk(ax,ay,C14,3)

         ! Add background magnetic field
         IF (myrank .eq. 0) THEN
            C12(1,1,1) = bx0*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
            C13(1,1,1) = by0*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
            C14(1,1,1) = bz0*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
         ENDIF

         ! Current density
         CALL curlk(C13,C14,ax,1)
         CALL curlk(C12,C14,ay,2)
         CALL curlk(C12,C13,az,3)

         ! Todo explore reducing 3 global arrays by calling a subroutine that
         ! computes gradre - lorentz (and hence C15,C16,C17 end up being local
         ! to that subroutine). In this case another subroutine 'emf' could
         ! compute the emf replacing the arrays for b.

         ! Advective term
!         CALL gradre(vx,vy,vz,C4,C5,C6)
         CALL prodre(vx,vy,vz,C4,C5,C6)

         ! Lorentz force
         CALL vector(ax,ay,az,C12,C13,C14,C15,C16,C17)

         ! Non-linear term
         C4 = C4 - C15
         C5 = C5 - C16
         C6 = C6 - C17

         ! Dealias non-linear
         CALL fc_filter(C4)
         CALL fc_filter(C5)
         CALL fc_filter(C6)

         ! Electromotive force
         CALL vector(vx,vy,vz,C12,C13,C14,C15,C16,C17)
        
         ! Dealias EMF
         CALL fc_filter(C15)
         CALL fc_filter(C16)
         CALL fc_filter(C17)

         ! lap(v)
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
            ax(k,j,i) = C9 (k,j,i) + dt*(-mu*ax(k,j,i)+C15(k,j,i)+&
                mx(k,j,i))*rmp
            ay(k,j,i) = C10(k,j,i) + dt*(-mu*ay(k,j,i)+C16(k,j,i)+&
                my(k,j,i))*rmp
            az(k,j,i) = C11(k,j,i) + dt*(-mu*az(k,j,i)+C17(k,j,i)+&
                mz(k,j,i))*rmp
         ENDDO
         ENDDO
         ENDDO


         ! Apply BC and project into solenoidal space
         ! NOTE: Magnetic boundary conditions don't work with moving
         ! boundaries yet
         CALL v_imposebc_and_project(vplanbc,planfc,vx,vy,vz,pr,rki=o)
         CALL a_imposebc_and_project(bplanbc,planfc,ax,ay,az,ph)
