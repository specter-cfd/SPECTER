!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!Evolve system in T time!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

step_evol = ini + NINT(T_guess/dt)

DO t_rk = ini,step_evol
   ! Runge-Kutta step 1
   ! Copies the fields into auxiliary arrays

   !$omp parallel do if (iend-ista.ge.nth) private (j,k)
   DO i = ista,iend
   !$omp parallel do if (iend-ista.lt.nth) private (k)
      DO j = 1,ny
         DO k = 1,nz

            INCLUDE 'include/bouss/bouss_rkstep1.f90'

         END DO
      END DO
   END DO

   ! Runge-Kutta step 2
   ! Evolves the system in time

   DO o = ord,1,-1
      INCLUDE 'include/bouss/bouss_rkstep1.f90'

   END DO
END DO 
