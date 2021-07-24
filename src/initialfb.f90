! Initial condition for the velocity.
! This file contains the expression used for the initial
! velocity field. You can use temporary real arrays R1-R3
! of size (1:nx,1:ny,ksta:kend) and temporary complex arrays
! C1-C6 of size (1:nz,1:ny,ista:iend) to do intermediate
! computations. The variable u0 should control the global
! amplitude of the velocity, and variables vparam0-9 can be
! used to control the amplitudes of individual terms. At the
! end, the three components of the velocity in spectral
! space should be stored in the arrays vx, vy, and vz.

! 3D BC-satisfying, incompressible noise
! From Clever & Busse "Tertiary and quaternary solutions for plane
! Couette flow". C1 = dPhi/dz, C2 = Phi, C3 = Psi
! For amplitude decay, remember that is Phi~k^-3 => v_phi~k^-1
! and Psi~k^-2 => v_psi~k^-1.

! Derivatives in z direction are calculated analytically, to have a
! more accurate initial condition (specially the no-slip)

! aparam0 and aparam1 are the minimum and maximum harmonic in the
! z direction kdn and kup are the miniumum and maximum
! perpendicular wavenumber.

      DO kk = int(mparam2), int(mparam3)

      IF (ista.eq.1) THEN
!$omp parallel do private (k,rmp,rmq,rms)
         DO j = 2,ny/2+1
            IF ((kk2(1,j,1).le.mparam1**2).and.(kk2(1,j,1).ge.mparam0**2)) THEN
               !Psi
               rmp = 2*pi*randu(seed)  !phase
               DO k = 1, nz-Cz
                  rmq = SIN(kk*pi*z(k)/Lz)/sqrt(kk2(1,j,1)+kk**2*pi**2/Lz**2)**2
                  C3(k,j,1) = (COS(rmp)+im*SIN(rmp))*rmq
                  C3(k,ny-j+2,1) = conjg(C3(k,j,1))
               END DO
            ENDIF
         ENDDO
         
!$omp parallel do if (iend-ista-1.ge.nth) private (j,k,rmp,rmq,rms)
         DO i = 2,iend
!$omp parallel do if (iend-ista-1.lt.nth) private (k,rmp,rmq,rms)
            DO j = 1,ny
               IF ((kk2(1,j,i).le.mparam1**2).and.(kk2(1,j,i).ge.mparam0**2)) THEN
               ! Psi
               rmp = 2*pi*randu(seed)  !phase
               DO k = 1, nz-Cz
                  rmq = SIN(kk*pi*z(k)/Lz)/sqrt(kk2(1,j,i)+kk**2*pi**2/Lz**2)**2
                  C3(k,j,i) = (COS(rmp)+im*SIN(rmp))*rmq
               END DO
               ENDIF
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,rmp,rmq,rms)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,rmp,rmq,rms)
            DO j = 1,ny
               IF ((kk2(1,j,i).le.mparam1**2).and.(kk2(1,j,i).ge.mparam0**2)) THEN
               ! Psi
               rmp = 2*pi*randu(seed)  !phase
               DO k = 1, nz-Cz
                  rmq = SIN(kk*pi*z(k)/Lz)/sqrt(kk2(1,j,i)+kk**2*pi**2/Lz**2)**2
                  C3(k,j,i) = (COS(rmp)+im*SIN(rmp))*rmq
               END DO
               ENDIF
            END DO
        END DO
      ENDIF
      END DO



      ! Toroidal 
      CALL derivk(C3,mx,2)
      CALL derivk(-C3,my,1)
      mz = 0

      ! 3D Fourier
      CALL fftp1d_real_to_complex_z(planfc,mx,MPI_COMM_WORLD)
      CALL fftp1d_real_to_complex_z(planfc,my,MPI_COMM_WORLD)
      CALL fftp1d_real_to_complex_z(planfc,mz,MPI_COMM_WORLD)

      ! Normalize
      CALL normvec(mx,my,mz,m0,0)

      mx = 0
      my = 0
      mz = 0
