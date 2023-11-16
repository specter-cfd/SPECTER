!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!Newton Solver!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Save initial field in 1d variable X0
CALL ThreeTo1D(vx, vy, vz, th, X0) 
! CALL Evol_T(X0, X_evol, T_guess)

!performs evolution of vx, vy, vz, th in time T
INCLUDE 'include/bouss/bouss_evol_T.f90' 

tind = tind + NINT(T_guess/(tstep*dt))

WRITE(ext, fmtext) tind


rmp = 1.0_GP/ &
	          (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
               DO j = 1,ny
                  DO k = 1,nz
                     C1(k,j,i) = vx(k,j,i)*rmp
                     C2(k,j,i) = vy(k,j,i)*rmp
                     C3(k,j,i) = vz(k,j,i)*rmp
                  END DO
               END DO
            END DO
            CALL fftp3d_complex_to_real(planfc,C1,R1,MPI_COMM_WORLD)
            CALL fftp3d_complex_to_real(planfc,C2,R2,MPI_COMM_WORLD)
            CALL fftp3d_complex_to_real(planfc,C3,R3,MPI_COMM_WORLD)
            IF (myrank.eq.0) THEN
               WRITE(*,*) 'Estoy a punto de escribir los campos'
            ENDIF
            CALL io_write(1,odir,'vx',ext,planio,R1)
            IF (myrank.eq.0) THEN
            WRITE(*,*) 'Ya escribí vx'
            ENDIF
            CALL io_write(1,odir,'vy',ext,planio,R2)
            CALL io_write(1,odir,'vz',ext,planio,R3)
            IF (myrank.eq.0) THEN
               WRITE(*,*) 'Ya escribí los 3 campos'
            ENDIF

! !Transforms to a 1d variable X_evol
CALL ThreeTo1D(vx, vy, vz, th, X_evol)

IF (myrank.eq.0) THEN
   WRITE(*,*) 'Ya pase evol a 3d'
ENDIF

! !Traslate in the guessed shifts:
CALL Traslation(X_evol, Y0, 1, sx) 
CALL Traslation(Y0, Y0, 2, sy) 
IF (myrank.eq.0) THEN
   WRITE(*,*) 'Ya trasladé en sx y sy'
ENDIF


!$omp parallel do
   DO i = 1,n_dim_1d
      dX0(i) = X0(i) - Y0(i)
   ENDDO
IF (myrank.eq.0) THEN
   WRITE(*,*) 'Ya calcule dX0'
ENDIF


!Perform small shifts:
CALL Traslation(Y0, Y_shift, 1, d_sx) 
CALL Traslation(Y_shift, Y_shift, 2, d_sy) 

IF (myrank.eq.0) THEN
   WRITE(*,*) 'Ya hice pequenos shifts'
ENDIF


!Calculate f(Y)
CALL f_Y_RB(Y0, f_Y)!COMENTARIO: Faltaría multiplicar todos los elementos por dT
IF (myrank.eq.0) THEN
   WRITE(*,*) 'Ya calcule f_Y'
ENDIF

!Calculate the directional derivative
CALL Perturb(X0, X_pert, dX0)
CALL OneTo3D(X_pert, vx, vy, vz, th)
!CALL Evol_T(X_pert, X_evol, T_guess)

IF (myrank.eq.0) THEN
   WRITE(*,*) 'Estoy por evolucionar sistema perturbado'
ENDIF

INCLUDE 'include/bouss/bouss_evol_T.f90' 
IF (myrank.eq.0) THEN
   WRITE(*,*) 'Termine de evolucionar sist perturbado'
ENDIF

! !Transforms to a 1d variable X_evol
CALL ThreeTo1D(vx, vy, vz, th, X_pert_evol) 


CALL X_fin_diff(X_partial_dif, X_pert_evol, Y0, sx, sy)
IF (myrank.eq.0) THEN
   WRITE(*,*) 'Calcule X_fin_diff'
ENDIF
