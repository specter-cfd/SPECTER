!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Subroutines to do Newton and stabilized biconjugate gradient
! methods to prepare data for the GPE solver. Can be easily
! modified to prepare initial conditions close to a fixed
! point for any solver. You should use the FFTPLANS and
! MPIVARS modules (see the file 'fftp_mod.f90') in each 
! program that calls any of the subroutines in this file. 
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2015 P.D. Mininni and M.E. Brachet
!
! 17 May 2018: Support for elongated box (N.Muller & P.D.Mininni) 
!=================================================================



!*****************************************************************
      SUBROUTINE ThreeTo1D(a,b,c,d,xvec1d)
!-----------------------------------------------------------------
!
! Copies the data into a 1D complex array to do the Newton method.
!
! Parameters
!     a     : x component
!     b     : y component
!     c     : z component
!     d     : thermal component
!     xvec1d  : 1D vector
!
      USE fprecision
      USE mpivars
      USE newtmod
      USE grid
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c,d
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(n_dim_1d)          :: xvec1D
      INTEGER             :: i,j,k
      INTEGER             :: offset1,offset2

!$omp parallel do if ((iend-ista).ge.nth) private(j,k,offset1,offset2)
      DO i = ista,iend
         offset1 = 4*(i-ista)*ny*nz
!$omp parallel do if ((iend-ista).lt.nth) private(k,offset2)
         DO j = 1,ny
            offset2 = offset1 + 4*(j-1)*nz
            DO k = 1,nz
               xvec1D(1+4*(k-1)+offset2) = a(k,j,i)
               xvec1D(2+4*(k-1)+offset2) = b(k,j,i)
               xvec1D(3+4*(k-1)+offset2) = c(k,j,i)
               xvec1D(4+4*(k-1)+offset2) = d(k,j,i)
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE ThreeTo1D



!*****************************************************************
      SUBROUTINE OneTo3D(xvec1d,a,b,c,d)
!-----------------------------------------------------------------
!
! Copies the data back into 3D complex arrays, after doing the 
! Newton method.
!
! Parameters
!     xvec1d  : 1D vector
!     a     : x component
!     b     : y component
!     c     : z component
!     d     : thermal component
!
      USE fprecision
      USE mpivars
      USE newtmod
      USE grid
      USE var
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: a,b,c,d
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n_dim_1d)            :: xvec1D
      INTEGER             :: i,j,k
      INTEGER             :: offset1,offset2

!$omp parallel do if ((iend-ista).ge.nth) private(j,k,offset1,offset2)
      DO i = ista,iend
         offset1 = 4*(i-ista)*ny*nz
!$omp parallel do if ((iend-ista).lt.nth) private(k,offset2)
         DO j = 1,ny
            offset2 = offset1+4*(j-1)*nz
            DO k = 1,nz
               a(k,j,i) = xvec1D(1+4*(k-1)+offset2)
               b(k,j,i) = xvec1D(2+4*(k-1)+offset2)
               c(k,j,i) = xvec1D(3+4*(k-1)+offset2)
               d(k,j,i) = xvec1D(4+4*(k-1)+offset2)
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE OneTo3D

!*****************************************************************
      SUBROUTINE Traslation(xvec1D, xvec1D_out, direc, d)
!-----------------------------------------------------------------
!
! Computes the matrix element -T_direc Y0.
! Where the operator T_direc is an infinitesimal  generator of the group
! of translations along direc. And direc = 1 is in x direction and
! direc = 2 is in y direction.
! Therefore, if direc=1: -T1 Y0(i,j,k) = Y0(i,j,k) * im * i / Lambdax 
! where Lambdax = x_dim / 2pi = 1
! 
! Parameters
!     Y  : 1D vector
!     direc  : direction in which to translate
!     Y_out  : 1D vector

      USE fprecision
      USE commtypes
      USE newtmod
      USE mpivars
      USE grid
      USE kes
      USE var
   !$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n_dim_1d) :: xvec1D 
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(n_dim_1d) :: xvec1D_out 
      INTEGER, INTENT(IN) :: direc !if direc = 1 traslates in x, if =2 traslates in y
      INTEGER :: i, j, k
      REAL(KIND=GP), INTENT(IN) :: d
      INTEGER :: offset1,offset2


      IF (direc.eq.1) THEN
            !$omp parallel do if ((iend-ista).ge.nth) private(j,k,offset1,offset2)
            DO i = ista,iend
            offset1 = 4*(i-ista)*ny*nz
                  !$omp parallel do if ((iend-ista).lt.nth) private(k,offset2)
                  DO j = 1,ny
                        offset2 = offset1 + 4*(j-1)*nz
                        DO k = 1,nz
                        xvec1D_out(1+4*(k-1)+offset2) = xvec1D(1+4*(k-1)+offset2) * EXP( im * kx(i)*d) 
                        xvec1D_out(2+4*(k-1)+offset2) = xvec1D(2+4*(k-1)+offset2) * EXP( im * kx(i)*d)
                        xvec1D_out(3+4*(k-1)+offset2) = xvec1D(3+4*(k-1)+offset2) * EXP( im * kx(i)*d)
                        xvec1D_out(4+4*(k-1)+offset2) = xvec1D(4+4*(k-1)+offset2) * EXP( im * kx(i)*d)
                        END DO
                  END DO
            END DO
      ELSE 
            !$omp parallel do if ((iend-ista).ge.nth) private(j,k,offset1,offset2)
            DO i = ista,iend
                  offset1 = 4*(i-ista)*ny*nz
      !$omp parallel do if ((iend-ista).lt.nth) private(k,offset2)
                  DO j = 1,ny
                  offset2 = offset1 + 4*(j-1)*nz
                  DO k = 1,nz
                        xvec1D_out(1+4*(k-1)+offset2) = xvec1D(1+4*(k-1)+offset2) * EXP(im * ky(j) * d)
                        xvec1D_out(2+4*(k-1)+offset2) = xvec1D(2+4*(k-1)+offset2) * EXP(im * ky(j) * d)
                        xvec1D_out(3+4*(k-1)+offset2) = xvec1D(3+4*(k-1)+offset2) * EXP(im * ky(j) * d)
                        xvec1D_out(4+4*(k-1)+offset2) = xvec1D(4+4*(k-1)+offset2) * EXP(im * ky(j) * d)
                  END DO
                  END DO
            END DO
      ENDIF
      RETURN 
      END SUBROUTINE Traslation


!*****************************************************************
     SUBROUTINE f_Y_RB(Y0, f_Y0)
!-----------------------------------------------------------------
!
! Computes the matrix element f(Y0) in Rayleigh Benard flow.
! Where f(Y0) = Y0_dot is the time derivative of the state Y0
! given by the flow governing equations.
! Parameters
!     Y0  : 1D complex vector
!     f_Y0  : 1D complex vector

      USE fprecision
      USE commtypes
      USE newtmod
      USE mpivars
      USE grid
   !$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n_dim_1d) :: Y0 
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(n_dim_1d) :: f_Y0
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: vx, vy, vz, th
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: fx, fy, fz, fs
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: C4, C5, C6, C8
      REAL(KIND=GP)    :: xmom,xtemp,nu,kappa
      INTEGER :: i, j, k
      INTEGER :: offset1,offset2

      CALL OneTo3D(Y0, vx, vy, vz, th)

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

      !$omp parallel do if ((iend-ista).ge.nth) private(j,k,offset1,offset2)
      DO i = ista,iend
      offset1 = 4*(i-ista)*ny*nz
            !$omp parallel do if ((iend-ista).lt.nth) private(k,offset2)
            DO j = 1,ny
                  offset2 = offset1 + 4*(j-1)*nz
                  DO k = 1,nz
                  f_Y0(1+4*(k-1)+offset2) = nu*vx(k,j,i)-C4(k,j,i)+fx(k,j,i)
                  f_Y0(2+4*(k-1)+offset2) = nu*vy(k,j,i)-C5(k,j,i)+fy(k,j,i)
                  f_Y0(3+4*(k-1)+offset2) = nu*vz(k,j,i)-C6(k,j,i)+fz(k,j,i)
                  f_Y0(4+4*(k-1)+offset2) = kappa*th(k,j,i)-C8(k,j,i)+fs(k,j,i)
                  END DO
            END DO
      END DO
      RETURN 
      END SUBROUTINE f_Y_RB

!*****************************************************************
     SUBROUTINE Perturb(X0, X_pert, c)
!-----------------------------------------------------------------
!
      USE fprecision
      USE commtypes
      USE newtmod
      USE mpivars
      USE grid
   !$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n_dim_1d) :: X0
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(n_dim_1d) :: X_pert
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n_dim_1d) :: c
      INTEGER :: i

!$omp parallel do
      DO i = 1,n_dim_1d
         X_pert(i) = X0(i) + epsilon * c(i)
      ENDDO

      ! !$omp parallel do if ((iend-ista).ge.nth) private(j,k,offset1,offset2)
      ! DO i = ista,iend
      ! offset1 = 4*(i-ista)*ny*nz
      !       !$omp parallel do if ((iend-ista).lt.nth) private(k,offset2)
      !       DO j = 1,ny
      !             offset2 = offset1 + 4*(j-1)*nz
      !             DO k = 1,nz
      !             X_pert(1+4*(k-1)+offset2) = X0(1+4*(k-1)+offset2) + epsilon * c(1+4*(k-1)+offset2)
      !             X_pert(2+4*(k-1)+offset2) = X0(2+4*(k-1)+offset2) + epsilon * c(2+4*(k-1)+offset2)
      !             X_pert(3+4*(k-1)+offset2) = X0(3+4*(k-1)+offset2) + epsilon * c(3+4*(k-1)+offset2)
      !             X_pert(4+4*(k-1)+offset2) = X0(4+4*(k-1)+offset2) + epsilon * c(4+4*(k-1)+offset2)
      !             END DO
      !       END DO
      ! END DO

      RETURN 
      END SUBROUTINE Perturb


!*****************************************************************
     SUBROUTINE X_fin_diff(X_partial_dif, X_evol, Y_1d, sx, sy)
!-----------------------------------------------------------------

      USE fprecision
      USE commtypes
      USE newtmod
      USE mpivars
      USE grid
      USE kes
      USE var
   !$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(n_dim_1d) :: X_partial_dif
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n_dim_1d) :: X_evol
      COMPLEX(KIND=GP), DIMENSION(n_dim_1d) :: X_evol_shift
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n_dim_1d) :: Y_1d
      REAL(KIND=GP), INTENT(IN)    :: sx, sy
      INTEGER :: i, j, k
      INTEGER :: offset1,offset2

      CALL Traslation(X_evol, X_evol_shift, 1, sx) !traslado en sx, sy dados por el guess inicial
      CALL Traslation(X_evol_shift, X_evol_shift, 2, sy) !COMENTARIO: faltaría definir sx, sy en algún lado


      !$omp parallel do if ((iend-ista).ge.nth) private(j,k,offset1,offset2)
      DO i = ista,iend
      offset1 = 4*(i-ista)*ny*nz
            !$omp parallel do if ((iend-ista).lt.nth) private(k,offset2)
            DO j = 1,ny
                  offset2 = offset1 + 4*(j-1)*nz
                  DO k = 1,nz
                  X_partial_dif(1+4*(k-1)+offset2) = (X_evol_shift(1+4*(k-1)+offset2) - Y_1d(1+4*(k-1)+offset2))/epsilon
                  X_partial_dif(2+4*(k-1)+offset2) = (X_evol_shift(2+4*(k-1)+offset2) - Y_1d(2+4*(k-1)+offset2))/epsilon
                  X_partial_dif(3+4*(k-1)+offset2) = (X_evol_shift(3+4*(k-1)+offset2) - Y_1d(3+4*(k-1)+offset2))/epsilon
                  X_partial_dif(4+4*(k-1)+offset2) = (X_evol_shift(4+4*(k-1)+offset2) - Y_1d(4+4*(k-1)+offset2))/epsilon
                  END DO
            END DO
      END DO

      RETURN 
      END SUBROUTINE X_fin_diff




!*****************************************************************
      SUBROUTINE scal(u1,u2,s)
!-----------------------------------------------------------------
!
! Routine to compute the reduced scalar product of two 1D
! vectors in double precision (even if GP=SINGLE).
!
! Parameters
!     u1      : First 1D vector
!     u2      : Second 1D vector
!     s       : at the output contais the reduced scalar product
!     n_dim_1d: size of the 1D vectors
!
      USE fprecision
      USE commtypes
      USE newtmod
      USE mpivars
      USE grid
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(OUT) :: s
      DOUBLE PRECISION              :: stemp,tmp
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n_dim_1d) :: u1,u2
      INTEGER :: i

      stemp = 0.0D0
      tmp = 1.0_GP/  &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
!$omp parallel do reduction(+:stemp)
      DO i = 1,n_dim_1d
         stemp = stemp+CONJG(u1(i))*u2(i)*tmp !no me acuerdo cual va conjugado
      ENDDO
      CALL MPI_ALLREDUCE(stemp,s,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                         MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE scal

!*****************************************************************
     SUBROUTINE CalculateProjection(c, X0, proj_f, proj_x, proj_y)
!-----------------------------------------------------------------
!
! Removes the projection of c = dX0 in the directions with traslational
! symmetries such as x and y. It also removes the components of dX0 
! along the direction of f(X0) such that the correction doesnt follow
! the original orbit.

      USE fprecision
      USE commtypes
      USE newtmod
      USE mpivars
      USE grid
      USE kes
      USE var
   !$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n_dim_1d) :: c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n_dim_1d) :: X0
      COMPLEX(KIND=GP), DIMENSION(n_dim_1d) :: Xtras_x
      COMPLEX(KIND=GP), DIMENSION(n_dim_1d) :: Xtras_y
      COMPLEX(KIND=GP), DIMENSION(n_dim_1d) :: f_X0
      DOUBLE PRECISION, INTENT(OUT) :: proj_f, proj_x, proj_y
      INTEGER :: i, j, k
      INTEGER :: offset1,offset2

      CALL f_Y_RB(X0, f_X0)
      CALL scal(c, f_X0, proj_f)
            
      !$omp parallel do if ((iend-ista).ge.nth) private(j,k,offset1,offset2)
      DO i = ista,iend
      offset1 = 4*(i-ista)*ny*nz
            !$omp parallel do if ((iend-ista).lt.nth) private(k,offset2)
            DO j = 1,ny
                  offset2 = offset1 + 4*(j-1)*nz
                  DO k = 1,nz
                  Xtras_x(1+4*(k-1)+offset2) = X0(1+4*(k-1)+offset2) *  im * kx(i) / Lx 
                  Xtras_x(2+4*(k-1)+offset2) = X0(2+4*(k-1)+offset2) *  im * kx(i) / Lx
                  Xtras_x(3+4*(k-1)+offset2) = X0(3+4*(k-1)+offset2) *  im * kx(i) / Lx
                  Xtras_x(4+4*(k-1)+offset2) = X0(4+4*(k-1)+offset2) *  im * kx(i) / Lx
                  END DO
            END DO
      END DO
      
      !$omp parallel do if ((iend-ista).ge.nth) private(j,k,offset1,offset2)
      DO i = ista,iend
      offset1 = 4*(i-ista)*ny*nz
            !$omp parallel do if ((iend-ista).lt.nth) private(k,offset2)
            DO j = 1,ny
                  offset2 = offset1 + 4*(j-1)*nz
                  DO k = 1,nz
                  Xtras_y(1+4*(k-1)+offset2) = X0(1+4*(k-1)+offset2) *  im * ky(j) / Ly 
                  Xtras_y(2+4*(k-1)+offset2) = X0(2+4*(k-1)+offset2) *  im * ky(j) / Ly
                  Xtras_y(3+4*(k-1)+offset2) = X0(3+4*(k-1)+offset2) *  im * ky(j) / Ly
                  Xtras_y(4+4*(k-1)+offset2) = X0(4+4*(k-1)+offset2) *  im * ky(j) / Ly
                  END DO
            END DO
      END DO

      CALL scal(c, Xtras_x, proj_x)
      CALL scal(c, Xtras_y, proj_y)


      !COMENTARIO: lo siguiente incluye sacar la proyección, pero puede no ser necesario
      ! CALL f_Y_RB(X0, f_X0)
      ! CALL scal(f_X0, f_X0, norm_f_X0)!COMENTARIO: chequear si scal funciona como está o hay que modificarla
      ! CALL scal(c, f_X0, s)

      ! !$omp parallel do if ((iend-ista).ge.nth) private(j,k,offset1,offset2)
      ! DO i = ista,iend
      ! offset1 = 4*(i-ista)*ny*nz
      !       !$omp parallel do if ((iend-ista).lt.nth) private(k,offset2)
      !       DO j = 1,ny
      !             offset2 = offset1 + 4*(j-1)*nz
      !             DO k = 1,nz
      !             c(1+4*(k-1)+offset2) = c(1+4*(k-1)+offset2) - s * f_X0(1+4*(k-1)+offset2)/norm_fX0
      !             c(2+4*(k-1)+offset2) = c(2+4*(k-1)+offset2) - s * f_X0(2+4*(k-1)+offset2)/norm_fX0
      !             c(3+4*(k-1)+offset2) = c(3+4*(k-1)+offset2) - s * f_X0(3+4*(k-1)+offset2)/norm_fX0
      !             c(4+4*(k-1)+offset2) = c(4+4*(k-1)+offset2) - s * f_X0(4+4*(k-1)+offset2)/norm_fX0
      !             END DO
      !       END DO
      ! END DO

      ! !Me faltaría sacar los comp. en x e y.

      RETURN 
      END SUBROUTINE CalculateProjection



! !*****************************************************************
!      SUBROUTINE Evol_T(X0, X_evol, T_guess)
! !-----------------------------------------------------------------
! !
! ! Evolves the 1d vector in a time given by T_guess. The 1d vector
! ! contains the vx, vy, vz, th fields.
! !
! ! Parameters
! !     X0  : 1D complex vector
! !     X_evol  : 1D complex vector
! !     T_guess  : real


!       USE fprecision
!       USE commtypes
!       USE newtmod
!       USE mpivars
!       USE grid
!       USE kes
!       USE var
!       USE order
!       USE pseudo_boundaries !COMENTARIO: no deja usar este módulo, seguro porque no está en pseudomod
!       USE fft
!    !$    USE threads
!       IMPLICIT NONE

!       COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n_dim_1d) :: X0
      ! COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(n_dim_1d) :: X_evol
      ! COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: vx, vy, vz, th, pr
      ! COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: fx, fy, fz, fs
      ! COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: C1, C2, C3, C4, C5, C6, C7, C8 
      ! REAL(KIND=GP),  DIMENSION (nz,ny,ista:iend)    :: R1
      ! REAL(KIND=GP), INTENT(IN)    :: T_guess !COMENTARIO: Tiene que tener GP o al darselo yo en el parameter.inp debería ser real nomás?
      ! REAL(KIND=GP)   :: dt,rmp
      ! REAL(KIND=GP)   :: xmom,xtemp,nu,kappa
      ! INTEGER   :: t, ini, step_rk
      ! INTEGER :: i, j, k
      ! INTEGER :: offset1,offset2
      ! INTEGER :: rki, o 

!       CALL OneTo3D(X0, vx, vy, vz, th)

!       step_rk = ini + NINT(T_guess/dt)

!       RK : DO t = ini,step_rk
!       ! Runge-Kutta step 1
! ! Copies the fields into auxiliary arrays

! !$omp parallel do if (iend-ista.ge.nth) private (j,k)
!       DO i = ista,iend
! !$omp parallel do if (iend-ista.lt.nth) private (k)
!       DO j = 1,ny
!       DO k = 1,nz

!       INCLUDE RKSTEP1_

!       END DO
!       END DO
!       END DO

! ! Runge-Kutta step 2
! ! Evolves the system in time

!       DO o = ord,1,-1
!       INCLUDE RKSTEP2_

!       END DO
!       END DO RK

      ! CALL ThreeTo1D(vx, vy, vz, th, X_evol)

      ! RETURN 
      ! END SUBROUTINE Evol_T




! !*****************************************************************
!       SUBROUTINE newton(xvec1D, sx, sy)
! !-----------------------------------------------------------------
! !no desarrollada todavía

!       USE fprecision
!       USE commtypes
!       USE newtmod
!       USE mpivars
!       USE grid
!       USE kes
!       USE var
!    !$    USE threads
!       IMPLICIT NONE

!       COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(n_dim_1d) :: xvec1D 
!       REAL, INTENT(IN) :: sx, sy !if direc = 1 traslates in x, if =2 traslates in y
!       INTEGER :: i, j, k
      
!       !Primero transformo a Y0.
!       CALL Traslation(xvec1d, 1, sx) !traslado en sx, sy dados por el guess inicial
!       CALL Traslation(xvec1d, 2, sy) !COMENTARIO: faltaría definir sx, sy en algún lado





!       RETURN 
!       END SUBROUTINE newton

! !*****************************************************************
!       SUBROUTINE Comp_T1_X_conj(X0, X0_out)
! !-----------------------------------------------------------------
! !
! ! Computes the matrix element (T1 X0)*, where * indicates
! ! the complex conjugate.
! ! 
! ! Parameters
! !     X  : 1D complex vector
! !     X_out  : 1D complex vector
!       USE fprecision
!       USE commtypes
!       USE newtmod
!       USE mpivars
!       USE grid
!    !$    USE threads
!       IMPLICIT NONE

!       COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n_dim_1d) :: X
!       COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n_dim_1d) :: X_out
!       INTEGER :: i

!       CALL(Comp_T1_Y(X,-X_out))
!       ! Es válido hacer eso? 

!       ! Ahora me queda conjugar el vector (trasponerlo se hará cuando se haga el producto)
! !$omp parallel do reduction(+:stemp)
!       DO i = 1,n_dim_1d
!          X_out(i) = conjg(X_out(i))
!       ENDDO

!       RETURN
!       END SUBROUTINE Comp_T1_X_conj

!*****************************************************************
!      SUBROUTINE Comp_direc_deriv(X, Y, Y_out)
!-----------------------------------------------------------------
!
! Computes the matrix element e^(-s_x T1) d/dX0 (X(T,X0))
! To do so, one approximates using finite differences like so:
! (exp(-s_x T1) X(T, X0+e c) - Y0)/e
! where e is an epsilon, and c is the direction in which to take the derivative,

! Therefore, one must evolve the system starting from X0+ e c in time T



! !*****************************************************************
!       SUBROUTINE gtrhs(xguess,rhs,cflow,dt,vx,vy,vz,vsq,R1,C3,C4,C5,C6)
! !-----------------------------------------------------------------
! !
! ! Computes -rhs of the ARGL equation.
! !
! ! Parameters
! !     xguess      : 1D vector state
! !     rhs         : rhs of the ARGL equation
! !     cflow       : =1 if mean constant velocity is present
! !     dt          : time step
! !     vx,vy,vz,vsq: arrays with the velocity field
! !     R1,C1-C6    : auxiliary arrays
! !
!       USE fprecision
!       USE mpivars
!       USE newtmod
!       USE grid
!       USE hbar
!       USE kes
!       USE ali
! !$    USE threads
!       IMPLICIT NONE
      
!       COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend)    :: vx,vy,vz
!       COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: C3,C4
!       COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: C5,C6
!       COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)              :: zre,zim
!       REAL(KIND=GP), INTENT(IN), DIMENSION(nx,ny,ksta:kend)     :: vsq
!       REAL(KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend)  :: R1
!       REAL(KIND=GP), INTENT(IN), DIMENSION(n_dim_1d)  :: xguess
!       REAL(KIND=GP), INTENT(OUT), DIMENSION(n_dim_1d) :: rhs
!       REAL(KIND=GP), INTENT(IN) :: dt
!       REAL(KIND=GP)       :: rmp,rmq
!       INTEGER, INTENT(IN) :: cflow
!       INTEGER             :: i,j,k

!       CALL OneTo3D(xguess,zre,zim)
      
! ! Here we basically have the contents of "include/argl_rkstep2.f90".
! ! We use Newton when we don't have counterflow, so we have no forcing
! ! function. However we still allow for counterflow.neq.0 for the case 
! ! in which we have a mean constant velocity (e.g., to advect a ring).
! ! The equations are:
! !    rhs.re = dt.(omegag.zre - beta.|z|^2 zre - |v^2|.zre/(4.alpha) +
! !             + v.grad(zim) - alpha.k^2.zre)/(1+alpha.k^2.dt)
! !    rhs.im = dt.(omegag.zim - beta.|z|^2 zim - |v^2|.zim/(4.alpha) -
! !             - v.grad(zre) - alpha.k^2.zim)/(1+alpha.k^2.dt)

!       CALL squareabs(zre,zim,R1,1)
!       rmq = real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)*omegag/beta
!       IF (cflow.eq.0) THEN ! If not doing counterflow we have |v^2|
! !$omp parallel do if (kend-ksta.ge.nth) private (j,i)
!          DO k = ksta,kend
! !$omp parallel do if (kend-ksta.lt.nth) private (i)
!             DO j = 1,ny
!                DO i = 1,nx
!                   R1(i,j,k) = rmq-R1(i,j,k)-vsq(i,j,k)
!                END DO
!             END DO
!          END DO
!       ELSE
! !$omp parallel do if (kend-ksta.ge.nth) private (j,i)
!          DO k = ksta,kend
! !$omp parallel do if (kend-ksta.lt.nth) private (i)
!             DO j = 1,ny
!                DO i = 1,nx
!                   R1(i,j,k) = rmq-R1(i,j,k)
!                END DO
!             END DO
!          END DO
!       ENDIF
!       CALL nonlgpe(R1,zre,C3)
!       CALL nonlgpe(R1,zim,C4)
!       CALL advect3(vx,vy,vz,zre,C5) ! -v.grad(zre)
!       CALL advect3(vx,vy,vz,zim,C6) ! -v.grad(zim)
! !$omp parallel do if (iend-ista.ge.nth) private (j,k)
!       DO i = ista,iend
! !$omp parallel do if (iend-ista.lt.nth) private (k)
!          DO j = 1,ny
!             DO k = 1,nz
!                IF (kn2(k,j,i).le.kmax) THEN
!                   rmp = dt/(1.0_GP+alpha*kk2(k,j,i)*dt)
!                   zre(k,j,i) = -(beta*C3(k,j,i)-C6(k,j,i)-    &
!                                alpha*kk2(k,j,i)*zre(k,j,i))*rmp
!                   zim(k,j,i) = -(beta*C4(k,j,i)+C5(k,j,i)-    &
!                                alpha*kk2(k,j,i)*zim(k,j,i))*rmp
!                ELSE
!                   zre(k,j,i) = 0.0_GP
!                   zim(k,j,i) = 0.0_GP
!                ENDIF
!             END DO
!          END DO
!       END DO
      
!       CALL ThreeTo1D(zre,zim,rhs)
!       RETURN
!       END SUBROUTINE gtrhs

! !*****************************************************************
!       SUBROUTINE lin(xguess,from,to,cflow,dt,vx,vy,vz,vsq,R1,C3,C4,C5,C6)
! !-----------------------------------------------------------------
! !
! ! Computes the linearized rhs of the ARGL equation:
! ! lin = dt.(omega.dz -|v|^2.dz/(4.alpha) - i v.grad(dz) - beta. 
! !       (2 |z|^2.dz + z^2.dz*) - alpha k^2 dz)/(1 + alpha k^2 dt)
! !
! ! Parameters
! !     xguess      : 1D vector state
! !     from        : perturbation to the 1D vector state
! !     to          : linearized rhs
! !     cflow       : =1 if mean constant velocity is present
! !     dt          : time step
! !     vx,vy,vz,vsq: arrays with the velocity field
! !     R1,C1-C6    : auxiliary arrays
! !
!       USE fprecision
!       USE commtypes
!       USE mpivars
!       USE newtmod
!       USE grid
!       USE hbar
!       USE kes
!       USE ali
!       USE fft
! !$    USE threads
!       IMPLICIT NONE
      
!       COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend)    :: vx,vy,vz
!       COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: C3,C4
!       COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: C5,C6
!       COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: C1,C2
!       COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)             :: zre,zim
!       COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)             :: zreg,zimg
!       REAL(KIND=GP), INTENT(IN), DIMENSION(nx,ny,ksta:kend)    :: vsq
!       REAL(KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend) :: R1
!       REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)       :: R2
!       REAL(KIND=GP), INTENT(IN), DIMENSION(n_dim_1d)  :: xguess,from
!       REAL(KIND=GP), INTENT(OUT), DIMENSION(n_dim_1d) :: to
!       REAL(KIND=GP), INTENT(IN) :: dt
!       REAL(KIND=GP)       :: rmp,rmq
!       INTEGER, INTENT(IN) :: cflow
!       INTEGER             :: i,j,k

!       CALL OneTo3D(xguess,zre,zim)
!       CALL OneTo3D(from,zreg,zimg)

!       CALL squareabs(zre,zim,R1,1)
!       rmq = real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)*omegag/beta
!       IF (cflow.eq.0) THEN ! If not doing counterflow we have |v^2|
! !$omp parallel do if (kend-ksta.ge.nth) private (j,i)
!          DO k = ksta,kend
! !$omp parallel do if (kend-ksta.lt.nth) private (i)
!             DO j = 1,ny
!                DO i = 1,nx
!                   R1(i,j,k) = rmq-2*R1(i,j,k)-vsq(i,j,k)
!                END DO
!             END DO
!          END DO
!       ELSE
! !$omp parallel do if (kend-ksta.ge.nth) private (j,i)
!          DO k = ksta,kend
! !$omp parallel do if (kend-ksta.lt.nth) private (i)
!             DO j = 1,ny
!                DO i = 1,nx
!                   R1(i,j,k) = rmq-2*R1(i,j,k)
!                END DO
!             END DO
!          END DO
!       ENDIF
!       CALL nonlgpe(R1,zreg,C3)
!       CALL nonlgpe(R1,zimg,C4)
!       CALL advect3(vx,vy,vz,zreg,C5)        ! -v.grad(zreg)
!       CALL advect3(vx,vy,vz,zimg,C6)        ! -v.grad(zimg)
!       CALL zsquare(zre,zim,R1,R2,1)
!       CALL pertgpe(R1,R2,zreg,zimg,zre,zim) ! (z^2).zreg*
! !$omp parallel do if (iend-ista.ge.nth) private (j,k)
!       DO i = ista,iend
! !$omp parallel do if (iend-ista.lt.nth) private (k)
!          DO j = 1,ny
!             DO k = 1,nz
!                IF (kn2(k,j,i).le.kmax) THEN
!                   rmp = dt/(1.0_GP+alpha*kk2(k,j,i)*dt)
!                   zre(k,j,i) = (beta*(C3(k,j,i)-zre(k,j,i))-C6(k,j,i)- &
!                                alpha*kk2(k,j,i)*zreg(k,j,i))*rmp
!                   zim(k,j,i) = (beta*(C4(k,j,i)-zim(k,j,i))+C5(k,j,i)- &
!                                alpha*kk2(k,j,i)*zimg(k,j,i))*rmp
!                ELSE
!                   zre(k,j,i) = 0.0_GP
!                   zim(k,j,i) = 0.0_GP
!                ENDIF
!             END DO
!          END DO
!       END DO

!       CALL ThreeTo1D(zre,zim,to)

!       RETURN
!       END SUBROUTINE lin

! !*****************************************************************
!       SUBROUTINE newton(zre,zim,cflow,dt,vx,vy,vz,vsq,R1,C3,C4,C5,C6)
! !-----------------------------------------------------------------
! !
! ! Newton's method to find a steady state by solving:
! !   eulstep(x)-x = 0.
! ! Starting from xguess, one iterates:
! !   x_{n+1} = x_{n} + deltax
! ! with
! !   lin deltax = rhs
! ! where
! !   rhs = -(eulstep(x_{n})-x_{n})
! ! and
! !   lin = d/dx (eulstep(x)-x) evaluated at x_{n}.
! !
! ! Parameters
! !     zre         : Real part of the wavefunction
! !     zim         : Imaginary part of the wavefunction
! !     cflow       : =1 if mean constant velocity is present
! !     dt          : time step
! !     vx,vy,vz,vsq: arrays with the velocity field
! !     R1,R2,C3-C6 : auxiliary arrays
! !     outs        : controls the amount of output
! !
!       USE fprecision
!       USE mpivars
!       USE newtmod
!       USE grid
!       USE var
! !$    USE threads
!       IMPLICIT NONE

!       COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend)    :: vx,vy,vz
!       COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: zre,zim
!       COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: C3,C4
!       COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: C5,C6
!       DOUBLE PRECISION    :: errnewt,ss
!       REAL(KIND=GP), INTENT(IN), DIMENSION(nx,ny,ksta:kend)    :: vsq
!       REAL(KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend) :: R1
!       REAL(KIND=GP), ALLOCATABLE, DIMENSION(:) :: xguess,rhs
!       REAL(KIND=GP), ALLOCATABLE, DIMENSION(:) :: deltax
!       REAL(KIND=GP), INTENT(IN) :: dt
!       INTEGER, INTENT(IN) :: cflow
!       INTEGER             :: i,j,k,inewt
      
!       ALLOCATE( xguess(1:n_dim_1d) )
!       ALLOCATE( rhs(1:n_dim_1d) )
!       ALLOCATE( deltax(1:n_dim_1d) )
!       CALL ThreeTo1D(zre,zim,xguess)

! ! The Newton loop
!       DO inewt=1,iter_max_newt

! ! Compute actual error
!          CALL gtrhs(xguess,rhs,cflow,dt,vx,vy,vz,vsq,R1,C3,C4,C5,C6)
!          CALL scal(rhs,rhs,ss)
!          errnewt = sqrt(ss)
!          IF (myrank.eq.0) &
!             PRINT *,' NEWTON LOOP, iter= ',inewt,' err= ',errnewt
!          IF (errnewt.lt.tol_newt) then
!             IF (myrank.eq.0) print *, 'Newton converged!'
!             CALL OneTo3D(xguess,zre,zim)
!             DEALLOCATE ( rhs,deltax )
!             RETURN
!          ELSE
!             CALL bicgstab(xguess,rhs,deltax,cflow,dt,errnewt,vx,vy,vz, &
!                           vsq,R1,C3,C4,C5,C6)
!             DO i = 1,n_dim_1d
!                xguess(i)=xguess(i)+deltax(i)
!             ENDDO
!          ENDIF
         
!       ENDDO
      
!       IF (myrank.eq.0) &
!            PRINT *,'Newton failed to converge in ',iter_max_newt,' steps!'
!       CALL OneTo3D(xguess,zre,zim)
!       DEALLOCATE ( rhs,deltax )

!       RETURN
!     END SUBROUTINE newton

! !*****************************************************************
!       SUBROUTINE scal(u1,u2,s)
! !-----------------------------------------------------------------
! !
! ! Routine to compute the reduced scalar product of two 1D
! ! vectors in double precision (even if GP=SINGLE).
! !
! ! Parameters
! !     u1      : First 1D vector
! !     u2      : Second 1D vector
! !     s       : at the output contais the reduced scalar product
! !     n_dim_1d: size of the 1D vectors
! !
!       USE fprecision
!       USE commtypes
!       USE newtmod
!       USE mpivars
!       USE grid
! !$    USE threads
!       IMPLICIT NONE

!       DOUBLE PRECISION, INTENT(OUT) :: s
!       DOUBLE PRECISION              :: stemp,tmp
!       REAL(KIND=GP), INTENT(IN), DIMENSION(n_dim_1d) :: u1,u2
!       INTEGER :: i

!       stemp = 0.0D0
!       tmp = 1.0_GP/  &
!             (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
! !$omp parallel do reduction(+:stemp)
!       DO i = 1,n_dim_1d
!          stemp = stemp+u1(i)*u2(i)*tmp
!       ENDDO
!       CALL MPI_ALLREDUCE(stemp,s,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
!                          MPI_COMM_WORLD,ierr)

!       RETURN
!       END SUBROUTINE scal

! !*****************************************************************
!       SUBROUTINE zsquare(a,b,ra,rb,dealias)
! !-----------------------------------------------------------------
! !
! ! Square of the complex wavefunction Z. Note the output is not
! ! normalized (i.e., not divided by N^3).
! !
! ! Parameters
! !     a : real part of the wavefunction in Fourier space
! !     b : imaginary part of the wavefunction in Fourier space
! !     ra: real part of Z^2 (= a^2-b^2) [output]
! !     rb: imag part of Z^2 (= 2 a.b  ) [output]
! !     dealias: =0 does not dealias the result
! !              =1 dealiases the result

!       USE fprecision
!       USE commtypes
!       USE mpivars
!       USE newtmod
!       USE grid
!       USE kes
!       USE ali
!       USE fft
! !$    USE threads
!       IMPLICIT NONE

!       COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
!       COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)             :: c1,c2
!       REAL(KIND=GP), INTENT(OUT), DIMENSION(nx,ny,ksta:kend)   :: ra,rb
!       REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: r1
!       REAL(KIND=GP)       :: rmp
!       INTEGER, INTENT(IN) :: dealias
!       INTEGER :: i,j,k

! !
! ! Computes real and imaginary parts
! !
! !$omp parallel do if (iend-ista.ge.nth) private (j,k)
!       DO i = ista,iend
! !$omp parallel do if (iend-ista.lt.nth) private (k)
!          DO j = 1,ny
!             DO k = 1,nz
!                c1(k,j,i) = a(k,j,i)
!                c2(k,j,i) = b(k,j,i)
!             END DO
!          END DO
!       END DO
!       CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
!       CALL fftp3d_complex_to_real(plancr,c2,rb,MPI_COMM_WORLD)
!       rmp = 1.0_GP/  &
!             (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
! !$omp parallel do if (kend-ksta.ge.nth) private (j,i)
!       DO k = ksta,kend
! !$omp parallel do if (kend-ksta.lt.nth) private (i)
!          DO j = 1,ny
!             DO i = 1,nx
!                ra(i,j,k) = (r1(i,j,k)**2-rb(i,j,k)**2)*rmp
!                rb(i,j,k) = 2*r1(i,j,k)*rb(i,j,k)*rmp
!             END DO
!          END DO
!       END DO
! !
! ! Dealiases the result and returns to real space
! !
!       IF (dealias.eq.1) THEN
!          CALL fftp3d_real_to_complex(planrc,ra,c1,MPI_COMM_WORLD)
!          CALL fftp3d_real_to_complex(planrc,rb,c2,MPI_COMM_WORLD)
! !$omp parallel do if (iend-ista.ge.nth) private (j,k)
!          DO i = ista,iend
! !$omp parallel do if (iend-ista.lt.nth) private (k)
!             DO j = 1,ny
!                DO k = 1,nz
!                   IF (kn2(k,j,i).gt.kmax) THEN
!                      c1(k,j,i) = 0.0_GP
!                      c2(k,j,i) = 0.0_GP
!                   ENDIF
!                END DO
!             END DO
!          END DO
!          CALL fftp3d_complex_to_real(plancr,c1,ra,MPI_COMM_WORLD)
!          CALL fftp3d_complex_to_real(plancr,c2,rb,MPI_COMM_WORLD)
!       ENDIF

!       RETURN
!       END SUBROUTINE zsquare

! !*****************************************************************
!       SUBROUTINE pertgpe(ra,rb,ca,cb,a,b)
! !-----------------------------------------------------------------
! !
! ! Computes one of the linearized terms of Z.|Z|^2 in real space:
! ! pertgpe = Z^2.dZ* = (Zre^2 - Zim^2 + 2i Zre.Zim)*(dZre - i dZim)
! !
! ! Parameters
! !     ra : input matrix with real(Z^2) in real space (not normalized)
! !     rb : input matrix with imag(Z^2) in real space (not normalized)
! !     ca : input matrix with real(dZ) in Fourier space
! !     cb : input matrix with imag(dZ) in Fourier space
! !     a  : real part of pertgpe in Fourier space [output]
! !     b  : imag part of pertgpe in Fourier space [output]
! !
!       USE fprecision
!       USE commtypes
!       USE mpivars
!       USE newtmod
!       USE grid
!       USE fft
! !$    USE threads
!       IMPLICIT NONE

!       COMPLEX(KIND=GP), INTENT(IN),  DIMENSION(nz,ny,ista:iend) :: ca,cb
!       COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: a,b
!       COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: c1,c2
!       REAL(KIND=GP), INTENT(IN), DIMENSION(nx,ny,ksta:kend)     :: ra,rb
!       REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: r1,r2,r3
!       REAL(KIND=GP)    :: rmp
!       INTEGER :: i,j,k

! !
! ! Computes real and imaginary parts
! !
! !$omp parallel do if (iend-ista.ge.nth) private (j,k)
!       DO i = ista,iend
! !$omp parallel do if (iend-ista.lt.nth) private (k)
!          DO j = 1,ny
!             DO k = 1,nz
!                c1(k,j,i) = ca(k,j,i)
!                c2(k,j,i) = cb(k,j,i)
!             END DO
!          END DO
!       END DO
!       CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
!       CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
!       rmp = 1.0_GP/ &
!             (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
! !$omp parallel do if (kend-ksta.ge.nth) private (j,i)
!       DO k = ksta,kend
! !$omp parallel do if (kend-ksta.lt.nth) private (i)
!          DO j = 1,ny
!             DO i = 1,nx
!                r3(i,j,k) = (ra(i,j,k)*r1(i,j,k)+rb(i,j,k)*r2(i,j,k))*rmp
!                r1(i,j,k) = (rb(i,j,k)*r1(i,j,k)-ra(i,j,k)*r2(i,j,k))*rmp
!             END DO
!          END DO
!       END DO
!       CALL fftp3d_real_to_complex(planrc,r3,a,MPI_COMM_WORLD)
!       CALL fftp3d_real_to_complex(planrc,r1,b,MPI_COMM_WORLD)

!       RETURN
!       END SUBROUTINE pertgpe
      
! !*****************************************************************
!       SUBROUTINE bicgstab(xguess,rhs,xnew,cflow,dt,errnewt,vx,vy,vz, &
!                           vsq,R1,C3,C4,C5,C6)
! !-----------------------------------------------------------------
! !
! ! Stabilized biconjugate gradient method.
! ! Solves lin xnew = rhs
! ! xnew and rhs are vectors of dimension n
! ! the action of linear operator lin on foo
! ! goo = lin(foo) is given by the call
! ! call lin(xguess,foo,goo,n)
! ! and the scalar product ss= foo1 . foo2
! ! is given by the call
! ! call scal(foo1,foo2,ss,n).
! ! The program exits when the following condition is met:
! ! sqrt((lin(xnew)-rhs).(lin(xnew)-rhs)) < tol
! !
! ! Parameters
! !     xguess: 1D vector state
! !     rhs   : rhs of the ARGL equation
! !     xnew  : correction to the 1D vector state
! !     vx,vy,vz,vsq : arrays with the velocity field
! !     R1,C3-C6     : auxiliary arrays
! !
!       USE fprecision
!       USE commtypes
!       USE mpivars
!       USE newtmod
!       USE grid
! !$    USE threads
!       IMPLICIT NONE

!       COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend)    :: vx,vy,vz
!       COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: C3,C4
!       COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: C5,C6
!       DOUBLE PRECISION, INTENT(IN) :: errnewt
!       DOUBLE PRECISION    :: ss,tt,ts,rhonew
!       DOUBLE PRECISION    :: omeganew,omegaold,rhoold
!       DOUBLE PRECISION    :: alpha,beta,err
      
!       REAL(KIND=GP), INTENT(IN), DIMENSION(nx,ny,ksta:kend)    :: vsq
!       REAL(KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend) :: R1
!       REAL(KIND=GP), INTENT(IN), DIMENSION(n_dim_1d)    :: xguess
!       REAL(KIND=GP), INTENT(INOUT), DIMENSION(n_dim_1d) :: rhs
!       REAL(KIND=GP), INTENT(OUT), DIMENSION(n_dim_1d)   :: xnew
!       REAL(KIND=GP), INTENT(IN) :: dt
!       REAL(KIND=GP), ALLOCATABLE, DIMENSION(:) :: pold,pnew,xold,rold,rnew
!       REAL(KIND=GP), ALLOCATABLE, DIMENSION(:) :: rhat0,vold,vnew,s,t
!       REAL(KIND=GP)       :: tol
!       INTEGER, INTENT(IN) :: cflow
!       INTEGER             :: i,ilp
     
!       ALLOCATE ( pold(1:n_dim_1d) )
!       ALLOCATE ( pnew(1:n_dim_1d) )
!       ALLOCATE ( xold(1:n_dim_1d) )
!       ALLOCATE ( rold(1:n_dim_1d) )
!       ALLOCATE ( rnew(1:n_dim_1d) )
!       ALLOCATE ( rhat0(1:n_dim_1d) )
!       ALLOCATE ( vold(1:n_dim_1d) )
!       ALLOCATE ( vnew(1:n_dim_1d) )
!       ALLOCATE ( s(1:n_dim_1d) )
!       ALLOCATE ( t(1:n_dim_1d) )
   
! ! Set absolute error goal
!       tol = tolbicg_rel*errnewt

! ! Set xold to zero
! !$omp parallel do 
!       DO i = 1,n_dim_1d
!          xold(i) = 0.0_GP
!       ENDDO
      
!       CALL lin(xguess,xold,rold,cflow,dt,vx,vy,vz,vsq,R1,C3,C4,C5,C6)
!       CALL gtrhs(xguess,rhs,cflow,dt,vx,vy,vz,vsq,R1,C3,C4,C5,C6)
! !$omp parallel do
!       DO i = 1,n_dim_1d
!          rold(i) = rhs(i)-rold(i)
!       ENDDO
! !$omp parallel do 
!       DO i = 1,n_dim_1d
!          rhat0(i) = rold(i)
!       ENDDO
!       rhoold = 1.0D0
!       alpha = 1.0D0
!       omegaold = 1.0D0
! !$omp parallel do
!       DO i = 1,n_dim_1d
!          vold(i) = 0.0_GP
!          pold(i) =0.0_GP
!       ENDDO
      
! ! This is the main bicgstab loop
!       DO ilp = 1,iter_max_bicg
!          CALL scal(rhat0,rold,rhonew)
!          beta = (rhonew/rhoold)*(alpha/omegaold)
         
! !$omp parallel do 
!          DO i = 1,n_dim_1d
!             pnew(i) = rold(i) + beta*(pold(i) - omegaold*vold(i))
!          ENDDO

!          CALL lin(xguess,pnew,vnew,cflow,dt,vx,vy,vz,vsq,R1,C3,C4,C5,C6)
!          CALL scal(rhat0,vnew,ss)
!          alpha = rhonew/ss

! !$omp parallel do
!          DO i = 1,n_dim_1d
!             s(i) = rold(i) - alpha*vnew(i)
!          ENDDO

!          CALL lin(xguess,s,t,cflow,dt,vx,vy,vz,vsq,R1,C3,C4,C5,C6)
!          call scal(s,s,ss)
!          err = sqrt(ss)
!          IF (myrank.eq.0) PRINT *, ' BiCgStab LOOP, iter= ',ilp,' err= ',err

!          IF (err.lt.tol) THEN
! !$omp parallel do
!             DO i=1,n_dim_1d
!                xnew(i) = xold(i) + alpha*pnew(i)
!             ENDDO
!             DEALLOCATE ( pold,pnew,xold,rold,rnew,rhat0,vold,vnew,s,t )
!             RETURN
!          ELSE
!             CALL scal(t,s,ts)
!             CALL scal(t,t,tt)
!             omeganew = ts/tt
! !$omp parallel do
!             DO i=1,n_dim_1d
!                xnew(i) = xold(i) + alpha*pnew(i) + omeganew*s(i)
!                rnew(i) = s(i) - omeganew*t(i)
!                xold(i) = xnew(i)
!                rold(i) = rnew(i)
!                vold(i) = vnew(i)
!                pold(i) = pnew(i)
!             ENDDO
!             rhoold = rhonew
!             omegaold = omeganew
!          ENDIF

!       ENDDO

! ! We're out of the loop without converging
!       IF (myrank.eq.0) PRINT *, ' BiCgStab FAILED TO CONVERGE IN ',iter_max_bicg,' ITERATIONS'
!       DEALLOCATE ( pold,pnew,xold,rold,rnew,rhat0,vold,vnew,s,t )

!       RETURN
!       END SUBROUTINE bicgstab
