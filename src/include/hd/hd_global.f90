         ! Save energy, enstrophy and bulk injection rate
         CALL hdcheck(vx,vy,vz,fx,fy,fz,t,dt,1,0)

         ! Print diagnostic quantities
         CALL vdiagnostic(vplanbc,planfc,vx,vy,vz,t,dt)
