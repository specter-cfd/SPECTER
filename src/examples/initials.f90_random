! Initial condition for the scalar field.
! This file contains the expression used for the initial
! scalar field. You can use temporary real arrays R1-R3
! of size (1:nx,1:ny,ksta:kend) and temporary complex arrays
! C1-C8 of size (1:nz,1:ny,ista:iend) to do intermediate
! computations. The variable c0 should control the global
! amplitude of the scalar, and variables cparam0-9 can be
! used to control the amplitudes of individual terms. At the
! end, the scalar field in spectral space should be stored
! in the array th.


      ! Noise that vanishes at z=0,Lz using a sin expansion
      ! and decays as 1/k

      ! Min and max wavenumbers indices in z
      DO kk = int(cparam0), int(cparam1)

      IF (ista.eq.1) THEN
!$omp parallel do if (iend-ista.lt.nth) private (k,rmp,rmq)
         DO j = 2,ny/2+1
            IF ((kk2(1,j,1).le.skup**2).and.(kk2(1,j,1).ge.skdn**2)) THEN
               rmq = 2*pi*randu(seed)
               rmp = 1.0_GP/sqrt(kk2(kk,j,1))
               DO k = 1, nz-Cz
                  ! Even contribution first, then odd contribution
                  th(k,j,1) = rmp*(COS(rmq)+im*SIN(rmq))*&
                              SIN(2*kk*pi/Lz*z(k))
                  th(k,j,1) = th(k,j,1) + rmp*(COS(rmq)+im*SIN(rmq))*&
                              SIN((2*kk-1)*pi/Lz*z(k))
                  th(k,ny-j+2,1) = conjg(th(k,j,1))
               END DO
            ENDIF
         ENDDO
!$omp parallel do if (iend-ista.ge.nth) private (j,k,rmp,rmq)
         DO i = 2,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,rmp,rmq)
            DO j = 1,ny
               IF ((kk2(1,j,i).le.skup**2).and.(kk2(1,j,i).ge.skdn**2)) THEN
                  rmq = 2*pi*randu(seed)
                  rmp = 1.0_GP/sqrt(kk2(kk,j,i))
                  DO k=1,nz-Cz
                  ! Even contribution first, then odd contribution
                  th(k,j,i) = rmp*(COS(rmq)+im*SIN(rmq))*&
                              SIN(2*kk*pi/Lz*z(k))
                  th(k,j,i) = th(k,j,i) + rmp*(COS(rmq)+im*SIN(rmq))*&
                              SIN((2*kk-1)*pi/Lz*z(k))
                  END DO
               ENDIF
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,rmp,rmq)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,rmp,rmq)
            DO j = 1,ny
               IF ((kk2(1,j,i).le.skup**2).and.(kk2(1,j,i).ge.skdn**2)) THEN
                  rmq = 2*pi*randu(seed)
                  rmp = 1.0_GP/sqrt(kk2(kk,j,i))
                  DO k = 1, nz-Cz
                  ! Even contribution first, then odd contribution
                  th(k,j,i) = rmp*(COS(rmq)+im*SIN(rmq))*&
                              SIN(2*kk*pi/Lz*z(k))
                  th(k,j,i) = th(k,j,i) + rmp*(COS(rmq)+im*SIN(rmq))*&
                              SIN((2*kk-1)*pi/Lz*z(k))
                  END DO
               ENDIF
            END DO
        END DO
      ENDIF
      END DO

      ! 3D Fourier
      CALL fftp1d_real_to_complex_z(planfc,th,MPI_COMM_WORLD)

      ! Normalize
      CALL normsca(th,c0,1)
