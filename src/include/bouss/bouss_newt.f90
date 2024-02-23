!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!Newton Solver!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!! Algorithm for solving near recurrences of Rayleigh Benard flow.

!! A solution for a recurrence is such that F(X) = 0, where X = (Omega_0, s_x, T) being Omega_0: velocity and thermal  
!! field at time t_0, s_x the guessed shift for relative periodic orbit, and T the guessed period of the orbit.  

!! F(X) is the difference between the shifted velocity and thermal field at time t+T (Omega_T) and Omega_0 

!! Each Newton step involves solving for s_k the equation J(X_k)*s_k = -F(X_k) (1) and then computing X_k+1 = X_k+s_k.
!! J(X_k) is the Jacobian of F(X)

!!Explain extra 2 equations -> AX = b

!! We solve AX = b by means of GMRES:



GMRES: DO n = 1,n_max
IF (myrank.eq.0) THEN
print *, 'Iteration:',n
ENDIF

!!!First we obtain the RHS of (1)

! On the first iteration we save X0 from the guessed fields
   IF (n.eq.1) THEN    ! In later iterations dX is updated by the Arnoldi algorithm.
      !Save initial field in 1d variable X0
      CALL ThreeTo1D(X0, vx, vy, vz, th) 

      !performs evolution of vx, vy, vz, th in time T
      INCLUDE 'include/bouss/bouss_evol_T.f90' 

      !Transforms to a 1d variable X_evol
      CALL ThreeTo1D(X_evol, vx, vy, vz, th)
      CALL Norm(aux_norm,X_evol)
      ! IF (myrank.eq.0) THEN
      ! print *, 'norm_X_evol=', aux_norm
      ! ENDIF
      
      !Traslate in the guessed shifts:
      CALL Translation(Y0, X_evol, 1, sx) 
      CALL Translation(Y0, Y0, 2, sy) 

      CALL Norm(aux_norm,Y0)
      ! IF (myrank.eq.0) THEN
      ! print *, 'norm_Y0=', aux_norm
      ! ENDIF

      !Calculate initial direction for directional derivative
      !$omp parallel do
         DO i = 1,n_dim_1d
            dX(i) = X0(i) - Y0(i)
         ENDDO

      CALL Norm(b_norm,dX)
      CALL Norm(aux_norm,X0)
      ! IF (myrank.eq.0) THEN
      ! print *, 'norm_b=', b_norm,'norm_X0=', aux_norm
      ! ENDIF

   ENDIF
   
   !!! We then compute the LHS of (1) term by term:

   !Calculate the directional derivative
   CALL Perturb(X_pert,epsilon,X0,dX )
   CALL OneTo3D(vx, vy, vz, th, X_pert)

   INCLUDE 'include/bouss/bouss_evol_T.f90' 


   !Transforms to a 1d variable X_evol
   CALL ThreeTo1D(X_pert_evol, vx, vy, vz, th) 

   !Calculates the directional derivative term    
   CALL X_fin_diff(X_partial_dif, X_pert_evol, Y0, sx, sy, epsilon)

   !Computes terms of translation generator:
   CALL Shift_term(Y_shift_x, Y0, 1, d_sx) 
   CALL Shift_term(Y_shift_y, Y0, 2, d_sy) 

   !TODO: Check if extra computation needed for the traslation terms

   !Calculate f(Y)
   CALL f_Y_RB(f_Y, Y0, dT)

   !!! Now for the rest of the first column of A:

   !Calculates the projection along the shifted directions and along the direction of flow
   CALL CalculateProjection(proj_f, proj_x, proj_y, dX, X0)
   IF (myrank.eq.0) THEN
   print *, 'Calculated Projection'
   ENDIF

   ! Can now form r_n = b - A*X_n

   CALL Form_Res(Res, Res_aux, dX, X_partial_dif, f_Y, Y_shift_x, Y_shift_y, proj_f, proj_x, proj_y, n)

   IF (myrank.eq.0) THEN
   print *, 'Formed Res'
   ENDIF


   CALL Arnoldi_step(Res, Res_aux, Q, Q_aux, H, res_norm, n)
   IF (myrank.eq.0) THEN
   print *, 'Performed Arnoldi_step'
   ENDIF

   CALL Update_values(dX, d_sx, d_sy, dT_guess, Res, Res_aux,n)
   IF (myrank.eq.0) THEN
   print *, 'Updated values'
   ENDIF

   CALL Givens_rotation(H, cs, sn, n)
   IF (myrank.eq.0) THEN
   print *, 'Performed Givens rotation'
   ENDIF

   CALL Update_error(beta, cs, sn, e, b_norm, res_norm, n)
   IF (myrank.eq.0) THEN
   print *, 'Updated error. e(n)=',e(n)
   ENDIF

   IF (e(n)<tol) THEN
      n_max = n
      EXIT GMRES
   !Check if it will lead to syntax error
   ENDIF

END DO GMRES

IF (myrank.eq.0) THEN
print *, 'Exited GMRES'
ENDIF

! mu = 0.0_GP

! DO n_hook = 1, n_hook_max

!    CALL Hookstep_transform(H, mu, n)
!    !TODO: Check if n=n
!    CALL Backpropagation(y, H, beta, n_max)

!    CALL Norm(y, norm_y)

!    IF (norm_y.gt.Delta) THEN
!       mu = mu + 0.01_GP
!    ELSE
!       EXIT
!    END IF

! END DO

CALL Backpropagation(y_sol, H, beta, n_max)
IF (myrank.eq.0) THEN
print *, 'Did Backpropagation',e(n)
ENDIF

!TODO: calculate m and check again if n=n
CALL Update_X(vx, vy, vz, th, X0, Q, y_sol, n)
IF (myrank.eq.0) THEN
print *, 'Updated X'
ENDIF

T_guess = T_guess + dT_guess
sx = sx + d_sx
sy = sy + d_sy