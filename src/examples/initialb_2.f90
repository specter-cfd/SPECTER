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

! vparam0 and vparam1 are the minimum and maximum Chandrasekhar-Reid
! harmonics to use. kdn and kup are the miniumum and maximum
! perpendicular wavenumber.

      DO kk = int(aparam0), int(aparam1)
      ! CR "wavenumbers"
      IF ( kk .eq. 1 ) THEN
         rm1 = 4.73004074
         rm2 = 7.85320462
      ELSE IF ( kk .eq. 2 ) THEN
         rm1 = 10.99560784 
         rm2 = 14.13716549
      ELSE IF ( kk .eq. 3 ) THEN
         rm1 = 17.27875966
         rm2 = 20.42035225
      ELSE IF ( kk .eq. 4 ) THEN
         rm1 = 23.56194490
         rm2 = 26.70353756
      ELSE
         rm1 = (2*kk - 0.5)*pi
         rm2 = (2*kk + 0.5)*pi
      ENDIF


      IF (ista.eq.1) THEN
!$omp parallel do private (k,rmp,rmq,rms)
         DO j = 2,ny/2+1
            IF ((kk2(1,j,1).le.kup**2).and.(kk2(1,j,1).ge.kdn**2)) THEN
               rmp = 2*pi*randu(seed)  !phase

               ! Phi
               DO k = 1, nz-Cz
                  ! Amplitude of even and odd contributions
                  rmq = COS(rm1)/sqrt(kk2(1,j,1)+rm1**2)**4
                  rms = SIN(rm2)/sqrt(kk2(1,j,1)+rm2**2)**4

                  ! Even contribution first, then odd contribution
                  C1(k,j,1) = rmq*rm1/Lz*(COS(rmp)+im*SIN(rmp))*(&
                        SINH(rm1*(z(k)/Lz-0.5))/COSH(0.5*rm1) +&
                        SIN(rm1*(z(k)/Lz-0.5))/COS(0.5*rm1))
                  C1(k,j,1) = C1(k,j,1) +  rms*rm2/Lz*(COS(rmp)+im*SIN(rmp))*(&
                        COSH(rm2*(z(k)/Lz-0.5))/SINH(0.5*rm2) -&
                        COS(rm2*(z(k)/Lz-0.5))/SIN(0.5*rm2))

                  C2(k,j,1) = rmq*(COS(rmp)+im*SIN(rmp))*(&
                        COSH(rm1*(z(k)/Lz-0.5))/COSH(0.5*rm1) -&
                        COS(rm1*(z(k)/Lz-0.5))/COS(0.5*rm1))
                  C2(k,j,1) = C2(k,j,1) + rms*(COS(rmp)+im*SIN(rmp))*(&
                        SINH(rm2*(z(k)/Lz-0.5))/SINH(0.5*rm2) -&
                        SIN(rm2*(z(k)/Lz-0.5))/SIN(0.5*rm2))
                 
                  C1(k,ny-j+2,1) = conjg(C1(k,j,1))
                  C2(k,ny-j+2,1) = conjg(C2(k,j,1))
               END DO

               !Psi
               rmp = 2*pi*randu(seed)
               DO k = 1, nz-Cz
                  rmq = SIN(kk*pi*z(k)/Lz)/sqrt(kk2(1,j,1)+kk**2*pi**2/Lz**2)**3
                  C3(k,j,1) = (COS(rmp)+im*SIN(rmp))*rmq
                  C3(k,ny-j+2,1) = conjg(C3(k,j,1))
               END DO
            ENDIF
         ENDDO
         
!$omp parallel do if (iend-ista-1.ge.nth) private (j,k,rmp,rmq,rms)
         DO i = 2,iend
!$omp parallel do if (iend-ista-1.lt.nth) private (k,rmp,rmq,rms)
            DO j = 1,ny
               IF ((kk2(1,j,i).le.kup**2).and.(kk2(1,j,i).ge.kdn**2)) THEN
               rmp = 2*pi*randu(seed)  !phase

               ! Phi
               DO k = 1, nz-Cz
                  ! Amplitude of even and odd contributions
                  rmq = COS(rm1)/sqrt(kk2(1,j,i)+rm1**2)**4
                  rms = SIN(rm2)/sqrt(kk2(1,j,i)+rm2**2)**4

                  ! Even contribution first, then odd contribution
                  C1(k,j,i) = rmq*rm1/Lz*(COS(rmp)+im*SIN(rmp))*(&
                        SINH(rm1*(z(k)/Lz-0.5))/COSH(0.5*rm1) +&
                        SIN(rm1*(z(k)/Lz-0.5))/COS(0.5*rm1))
                  C1(k,j,i) = C1(k,j,i) +  rms*rm2/Lz*(COS(rmp)+im*SIN(rmp))*(&
                        COSH(rm2*(z(k)/Lz-0.5))/SINH(0.5*rm2) -&
                        COS(rm2*(z(k)/Lz-0.5))/SIN(0.5*rm2))

                  C2(k,j,i) = rmq*(COS(rmp)+im*SIN(rmp))*(&
                        COSH(rm1*(z(k)/Lz-0.5))/COSH(0.5*rm1) -&
                        COS(rm1*(z(k)/Lz-0.5))/COS(0.5*rm1))
                  C2(k,j,i) = C2(k,j,i) + rms*(COS(rmp)+im*SIN(rmp))*(&
                        SINH(rm2*(z(k)/Lz-0.5))/SINH(0.5*rm2) -&
                        SIN(rm2*(z(k)/Lz-0.5))/SIN(0.5*rm2))
               END DO

               ! Psi
               rmp = 2*pi*randu(seed)
               DO k = 1, nz-Cz
                  rmq = SIN(kk*pi*z(k)/Lz)/sqrt(kk2(1,j,i)+kk**2*pi**2/Lz**2)**3
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
               IF ((kk2(1,j,i).le.kup**2).and.(kk2(1,j,i).ge.kdn**2)) THEN
               rmp = 2*pi*randu(seed)  !phase

               ! Phi
               DO k = 1, nz-Cz
                  ! Amplitude of even and odd contributions
                  rmq = COS(rm1)/sqrt(kk2(1,j,i)+rm1**2)**4
                  rms = SIN(rm2)/sqrt(kk2(1,j,i)+rm2**2)**4

                  ! Even contribution first, then odd contribution
                  C1(k,j,i) = rmq*rm1/Lz*(COS(rmp)+im*SIN(rmp))*(&
                        SINH(rm1*(z(k)/Lz-0.5))/COSH(0.5*rm1) +&
                        SIN(rm1*(z(k)/Lz-0.5))/COS(0.5*rm1))
                  C1(k,j,i) = C1(k,j,i) +  rms*rm2/Lz*(COS(rmp)+im*SIN(rmp))*(&
                        COSH(rm2*(z(k)/Lz-0.5))/SINH(0.5*rm2) -&
                        COS(rm2*(z(k)/Lz-0.5))/SIN(0.5*rm2))

                  C2(k,j,i) = rmq*(COS(rmp)+im*SIN(rmp))*(&
                        COSH(rm1*(z(k)/Lz-0.5))/COSH(0.5*rm1) -&
                        COS(rm1*(z(k)/Lz-0.5))/COS(0.5*rm1))
                  C2(k,j,i) = C2(k,j,i) + rms*(COS(rmp)+im*SIN(rmp))*(&
                        SINH(rm2*(z(k)/Lz-0.5))/SINH(0.5*rm2) -&
                        SIN(rm2*(z(k)/Lz-0.5))/SIN(0.5*rm2))
               END DO

               ! Psi
               rmp = 2*pi*randu(seed)
               DO k = 1, nz-Cz
                  rmq = SIN(kk*pi*z(k)/Lz)/sqrt(kk2(1,j,i)+kk**2*pi**2/Lz**2)**3
                  C3(k,j,i) = (COS(rmp)+im*SIN(rmp))*rmq
               END DO
               ENDIF
            END DO
        END DO
      ENDIF
      END DO


!      ! v normal
!      CALL derivk(C2,C4,1)
!      CALL derivk(C4,C5,1)
!      CALL derivk(C2,C4,2)
!      CALL derivk(C4,C6,2)
!      az = -(C5+C6)
!
!      ! v tangential
!      ! Poloidal
!      CALL derivk(C1,C4,1)
!      CALL derivk(C1,C5,2)
!      ax = C4
!      ay = C5

!      ! Toroidal 
!      CALL derivk(C3,C4,2)
!      CALL derivk(C3,C5,1)
!      ax = ax+C4
!      ay = ay-C5

      CALL derivk(

      ! 3D Fourier
      CALL fftp1d_real_to_complex_z(planfc,ax,MPI_COMM_WORLD)
      CALL fftp1d_real_to_complex_z(planfc,ay,MPI_COMM_WORLD)
      CALL fftp1d_real_to_complex_z(planfc,az,MPI_COMM_WORLD)

      ! Normalize
      CALL normvec(ax,ay,az,a0,0)
