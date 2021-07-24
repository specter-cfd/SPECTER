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

      fx = 0.0_GP
      fy = 0.0_GP
      fz = 0.0_GP

      IF (myrank .eq. 0) THEN
         fx(1,1,1) = f0*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
      ENDIF
