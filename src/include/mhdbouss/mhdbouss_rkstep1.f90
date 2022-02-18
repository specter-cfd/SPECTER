! Step 1 of Runge-Kutta for the Boussinesq MHD equations
! Copies v_x,y,z into the auxiliary arrays C1-C3
! Copies b_x,y,z into the auxiliary arrays C9-C11
! Copies th      into the auxiliary array  C7

         C1(k,j,i)  = vx(k,j,i)
         C2(k,j,i)  = vy(k,j,i)
         C3(k,j,i)  = vz(k,j,i)
         C7(k,j,i)  = th(k,j,i)
         C9(k,j,i)  = ax(k,j,i)
         C10(k,j,i) = ay(k,j,i)
         C11(k,j,i) = az(k,j,i)
