! External forcing for the scalar
! This file contains the expression used for the external
! scalar forcing. You can use temporary real arrays R1-R3
! of size (1:nx,1:ny,ksta:kend) and temporary complex arrays
! C1-C8 of size (1:nz,1:ny,ista:iend) to do intermediate
! computations. The variable s0 should control the global
! amplitude of the forcing, and variables sparam0-9 can be
! used to control the amplitudes of individual terms. At the
! end, the forcing in spectral in the arrays fs.

! Null mechanical forcing
      DO i = ista,iend
         DO j = 1,ny
            DO k = 1,nz
               fs(k,j,i) = 0.
            END DO
         END DO
      END DO
