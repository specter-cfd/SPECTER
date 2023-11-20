!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!Newton Solver!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Save initial field in 1d variable X0
CALL ThreeTo1D(vx, vy, vz, th, X0) 
! CALL Evol_T(X0, X_evol, T_guess)

!performs evolution of vx, vy, vz, th in time T
INCLUDE 'include/bouss/bouss_evol_T.f90' 

!Transforms to a 1d variable X_evol
CALL ThreeTo1D(vx, vy, vz, th, X_evol)
! IF (myrank.eq.0) THEN
!    WRITE(*,*) 'Pase X_evol a 3d'
! ENDIF

!Traslate in the guessed shifts:
CALL Traslation(X_evol, Y0, 1, sx) 
CALL Traslation(Y0, Y0, 2, sy) 

!Calculate initial direction for directional derivative
!$omp parallel do
   DO i = 1,n_dim_1d
      dX0(i) = X0(i) - Y0(i)
   ENDDO

!Perform small shifts:
CALL Traslation(Y0, Y_shift, 1, d_sx) 
CALL Traslation(Y_shift, Y_shift, 2, d_sy) 

!Calculate f(Y)
CALL f_Y_RB(Y0, f_Y)!COMENTARIO: Faltaría multiplicar todos los elementos por dT

!Calculate the directional derivative
CALL Perturb(X0, X_pert, dX0)
CALL OneTo3D(X_pert, vx, vy, vz, th)
 !CALL Evol_T(X_pert, X_evol, T_guess)

INCLUDE 'include/bouss/bouss_evol_T.f90' 

!Transforms to a 1d variable X_evol
CALL ThreeTo1D(vx, vy, vz, th, X_pert_evol) 

!Calculates the directional derivative term
CALL X_fin_diff(X_partial_dif, X_pert_evol, Y0, sx, sy)