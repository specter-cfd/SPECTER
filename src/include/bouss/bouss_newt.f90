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

!!!First we obtain the RHS of (1)

! On the first iteration we save X0 from the guessed fields
   IF (n.eq.1) THEN    ! In later iterations dX is updated by the Arnoldi algorithm.
      !Save initial field in 1d variable X0
      CALL ThreeTo1D(vx, vy, vz, th, X0) 


      !performs evolution of vx, vy, vz, th in time T
      INCLUDE 'include/bouss/bouss_evol_T.f90' 

      !Transforms to a 1d variable X_evol
      CALL ThreeTo1D(vx, vy, vz, th, X_evol)

      !Traslate in the guessed shifts:
      CALL Traslation(X_evol, Y0, 1, sx) 
      CALL Traslation(Y0, Y0, 2, sy) 

      !Calculate initial direction for directional derivative
      !$omp parallel do
         DO i = 1,n_dim_1d
            dX(i) = X0(i) - Y0(i)
         ENDDO
      
      CALL Norm(dX, b_norm)
   ENDIF
   
   !!! We then compute the LHS of (1) term by term:

   !Calculate the directional derivative
   CALL Perturb(X0, X_pert, dX, epsilon)
   CALL OneTo3D(X_pert, vx, vy, vz, th)
   !CALL Evol_T(X_pert, X_evol, T_guess)

   INCLUDE 'include/bouss/bouss_evol_T.f90' 

   !Transforms to a 1d variable X_evol
   CALL ThreeTo1D(vx, vy, vz, th, X_pert_evol) 

   !Calculates the directional derivative term    
   CALL X_fin_diff(X_partial_dif, X_pert_evol, Y0, sx, sy, epsilon)

   !Computes infinitesi:
   CALL Shift_term(Y0, Y_shift_x, 1, d_sx) 
   CALL Shift_term(Y0, Y_shift_y, 2, d_sy) 
   !TODO: Check if extra computation needed for the traslation terms

   !Calculate f(Y)
   CALL f_Y_RB(Y0, f_Y, dT)

   !!! Now for the rest of the first column of A:

   !Calculates the projection along the shifted directions and along the direction of flow
   CALL CalculateProjection(dX, X0, proj_f, proj_x, proj_y)

   ! Can now form r_n = b - A*X_n

   CALL Form_Res(Res, dX, X_partial_diff, proj_f, proj_x, proj_y, f_Y, Y_shift_x, Y_shift_y, n)

   CALL Arnoldi_step(Res, Q, H, n)

   CALL Update_values(Res, dX, d_sx, d_sy, dT)

   CALL Givens_rotation(H, cs, sn, n)

   CALL Update_error(beta, cs, sn, e, b_norm, n)

   IF e(n)<tol THEN
      n_max = n
      EXIT GMRES
   !Check if it will lead to syntax error

END DO GMRES

mu = 0.0_GP

DO n_hook = 1, n_hook_max

   CALL Hookstep_transform(H, mu, n)
   !TODO: Check if n=n
   CALL Backpropagation(H, beta, n_max, y)

   CALL Norm(y, norm_y)

   IF (norm_y.gt.Delta) THEN
      mu = mu + 0.01_GP
   ELSE
      EXIT
   END IF

END DO

!TODO: calculate m and check again if n=n
CALL Update_X(X0, Q, y, vx, vy, vz, th, n, m)


