         ! Save energy, enstrophy and bulk injection rate
         CALL hdcheck(vx,vy,vz,fx,fy,fz,t,dt,1,0)

         ! Print diagnostic quantities
         ! Mean squared divergence, <vx^2+vy^2> at z=0,Lz
         ! and <vz^2> at z=0,Lz
         CALL diagnostic(vx,vy,vz,t,dt)
