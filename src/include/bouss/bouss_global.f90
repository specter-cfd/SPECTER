         ! Save velocity field energy, enstrophy and bulk injection rate
         CALL hdcheck(vx,vy,vz,fx,fy,fz,t,dt,1,0)

         ! Save scalar field energy, squared gradient and injection rate
         CALL pscheck(th,fs,t,dt)

         ! Print diagnostic quantities
         CALL vdiagnostic(vplanbc,planfc,vx,vy,vz,t,dt)
         CALL sdiagnostic(splanbc,planfc,th,t,dt)
