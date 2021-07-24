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

      C1 = 0; C2 = 0; C3 = 0
      R1(1,1,ksta) =  4.73004074
      R1(2,1,ksta) =  7.85320462
      R1(3,1,ksta) = 10.99560784 
      R1(4,1,ksta) = 14.13716549 
      R1(5,1,ksta) = 17.27875966 
      R1(6,1,ksta) = 20.42035225 
      R1(7,1,ksta) = 23.56194490 
      R1(8,1,ksta) = 26.70353756 

      DO kk = int(vparam0), int(vparam1)
      ! CR "wavenumbers"
      IF ( kk .le. 8 ) THEN
         rm1 = R1(kk,1,ksta)
      ELSE IF ( MODULO(kk, 2) .eq. 1 ) THEN
         rm1 = (2*kk - 0.5)*pi
      ELSE
         rm1 = (2*kk + 0.5)*pi
      ENDIF

      IF (ista.eq.1) THEN
!$omp parallel do private (k,rmp,rmq)
          DO j = 1,ny/2+1
            IF ((kk2(1,j,1).le.kup**2) .AND. (kk2(1,j,1).ge.kdn**2)) THEN
               rmp = 2*pi*randu(seed)                       ! Phase
               rmq = randu(seed)/sqrt(kk2(1,j,1)+rm1**2)**3 ! Amplitude

               ! Phi is in C2 and C1 is dPhi/dz
               DO k = 1, nz-Cz
                  IF ( MODULO(kk, 2) .eq. 1 ) THEN
                  C1(k,j,1) = C1(k,j,1) + rmq*rm1/Lz*(COS(rmp)+im*SIN(rmp))*(&
                        SINH(rm1*(z(k)/Lz-0.5))/COSH(0.5*rm1) +&
                        SIN(rm1*(z(k)/Lz-0.5))/COS(0.5*rm1))

                  C2(k,j,1) = C2(k,j,1) + rmq*(COS(rmp)+im*SIN(rmp))*(&
                        COSH(rm1*(z(k)/Lz-0.5))/COSH(0.5*rm1) -&
                        COS(rm1*(z(k)/Lz-0.5))/COS(0.5*rm1))
                  ELSE
                  C1(k,j,1) = C1(k,j,1) + rmq*rm1/Lz*(COS(rmp)+im*SIN(rmp))*(&
                        COSH(rm1*(z(k)/Lz-0.5))/SINH(0.5*rm1) -&
                        COS(rm1*(z(k)/Lz-0.5))/SIN(0.5*rm1))

                  C2(k,j,1) = C2(k,j,1) + rmq*(COS(rmp)+im*SIN(rmp))*(&
                        SINH(rm1*(z(k)/Lz-0.5))/SINH(0.5*rm1) -&
                        SIN(rm1*(z(k)/Lz-0.5))/SIN(0.5*rm1))
                  ENDIF
                 
                  C1(k,ny-j+2,1) = conjg(C1(k,j,1))
                  C2(k,ny-j+2,1) = conjg(C2(k,j,1))
               END DO
               
               !Psi
               IF ( MODULO(kk, 2) .eq. 1 ) THEN
                   rmp = rmp 
               ELSE
                   rmp = pi/2 + rmp
               ENDIF
               rmq = randu(seed)/sqrt(kk2(1,j,1)+kk**2*pi**2/Lz**2)**2 ! Amp
               DO k = 1, nz-Cz
                  C3(k,j,1) = C3(k,j,1) + rmq*(COS(rmp)+im*SIN(rmp))*SIN(kk*pi*z(k)/Lz)
                  C3(k,ny-j+2,1) = conjg(C3(k,j,1))
               END DO
            ENDIF
         ENDDO
         
!$omp parallel do if (iend-ista-1.ge.nth) private (j,k,rmp,rmq)
         DO i = 2,iend
!$omp parallel do if (iend-ista-1.lt.nth) private (k,rmp,rmq)
            DO j = 1,ny
               IF ((kk2(1,j,i).le.kup**2) .AND. (kk2(1,j,i).ge.kdn**2)) THEN
               rmp = 2*pi*randu(seed)                       ! Phase
               rmq = randu(seed)/sqrt(kk2(1,j,i)+rm1**2)**3 ! Amplitude

               ! Phi is in C2 and C1 is dPhi/dz
               DO k = 1, nz-Cz
                  IF ( MODULO(kk, 2) .eq. 1 ) THEN
                  C1(k,j,i) = C1(k,j,i) + rmq*rm1/Lz*(COS(rmp)+im*SIN(rmp))*(&
                        SINH(rm1*(z(k)/Lz-0.5))/COSH(0.5*rm1) +&
                        SIN(rm1*(z(k)/Lz-0.5))/COS(0.5*rm1))

                  C2(k,j,i) = C2(k,j,i) + rmq*(COS(rmp)+im*SIN(rmp))*(&
                        COSH(rm1*(z(k)/Lz-0.5))/COSH(0.5*rm1) -&
                        COS(rm1*(z(k)/Lz-0.5))/COS(0.5*rm1))
                  ELSE
                  C1(k,j,i) = C1(k,j,i) + rmq*rm1/Lz*(COS(rmp)+im*SIN(rmp))*(&
                        COSH(rm1*(z(k)/Lz-0.5))/SINH(0.5*rm1) -&
                        COS(rm1*(z(k)/Lz-0.5))/SIN(0.5*rm1))

                  C2(k,j,i) = C2(k,j,i) + rmq*(COS(rmp)+im*SIN(rmp))*(&
                        SINH(rm1*(z(k)/Lz-0.5))/SINH(0.5*rm1) -&
                        SIN(rm1*(z(k)/Lz-0.5))/SIN(0.5*rm1))
                  ENDIF
               END DO

               IF ( MODULO(kk, 2) .eq. 1 ) THEN
                   rmp = rmp 
               ELSE
                   rmp = pi/2 + rmp
               ENDIF
               rmq = randu(seed)/sqrt(kk2(1,j,i)+kk**2*pi**2/Lz**2)**2 ! Amp
               DO k = 1, nz-Cz
                  C3(k,j,i) = C3(k,j,i) + rmq*(COS(rmp)+im*SIN(rmp))*SIN(kk*pi*z(k)/Lz)
               END DO
               ENDIF
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,rmp,rmq)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,rmp,rmq)
            DO j = 1,ny
               IF ((kk2(1,j,i).le.kup**2) .AND. (kk2(1,j,i).ge.kdn**2)) THEN
               rmp = 2*pi*randu(seed)                        ! Phase
               rmq = randu(seed)/sqrt(kk2(1,j,i)+rm1**2)**3  ! Amplitude

               ! Phi is in C2 and C1 is dPhi/dz
               DO k = 1, nz-Cz
                  IF ( MODULO(kk, 2) .eq. 1 ) THEN
                  C1(k,j,i) = C1(k,j,i) + rmq*rm1/Lz*(COS(rmp)+im*SIN(rmp))*(&
                        SINH(rm1*(z(k)/Lz-0.5))/COSH(0.5*rm1) +&
                        SIN(rm1*(z(k)/Lz-0.5))/COS(0.5*rm1))

                  C2(k,j,i) = C2(k,j,i) + rmq*(COS(rmp)+im*SIN(rmp))*(&
                        COSH(rm1*(z(k)/Lz-0.5))/COSH(0.5*rm1) -&
                        COS(rm1*(z(k)/Lz-0.5))/COS(0.5*rm1))
                  ELSE
                  C1(k,j,i) = C1(k,j,i) + rmq*rm1/Lz*(COS(rmp)+im*SIN(rmp))*(&
                        COSH(rm1*(z(k)/Lz-0.5))/SINH(0.5*rm1) -&
                        COS(rm1*(z(k)/Lz-0.5))/SIN(0.5*rm1))

                  C2(k,j,i) = C2(k,j,i) + rmq*(COS(rmp)+im*SIN(rmp))*(&
                        SINH(rm1*(z(k)/Lz-0.5))/SINH(0.5*rm1) -&
                        SIN(rm1*(z(k)/Lz-0.5))/SIN(0.5*rm1))
                  ENDIF
               END DO

               ! Psi
               IF ( MODULO(kk, 2) .eq. 1 ) THEN
                   rmp = rmp 
               ELSE
                   rmp = pi/2 + rmp
               ENDIF
               rmq = randu(seed)/sqrt(kk2(1,j,i)+kk**2*pi**2/Lz**2)**2 ! Amp
               DO k = 1, nz-Cz
                  C3(k,j,i) = C3(k,j,i) + rmq*(COS(rmp)+im*SIN(rmp))*SIN(kk*pi*z(k)/Lz)
               END DO
               ENDIF
            END DO
        END DO
      ENDIF
      END DO

      ! v normal
      CALL derivk(C2,C4,1)
      CALL derivk(C4,C5,1)
      CALL derivk(C2,C4,2)
      CALL derivk(C4,C6,2)
      vz = -(C5+C6)

      ! v tangential
      ! Poloidal
      CALL derivk(C1,C4,1)
      CALL derivk(C1,C5,2)
      vx = C4
      vy = C5

      ! Toroidal 
      CALL derivk(C3,C4,2)
      CALL derivk(C3,C5,1)
      vx = vx+C4
      vy = vy-C5

      ! 3D Fourier
      CALL fftp1d_real_to_complex_z(planfc,vx,MPI_COMM_WORLD)
      CALL fftp1d_real_to_complex_z(planfc,vy,MPI_COMM_WORLD)
      CALL fftp1d_real_to_complex_z(planfc,vz,MPI_COMM_WORLD)

      ! Normalize
      CALL normvec(vx,vy,vz,u0,1)
