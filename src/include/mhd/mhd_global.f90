         ! Save total energy, enstrophy and bulk injection rate
         CALL mhdcheck(vx,vy,vz,ax,ay,az,t,dt,1,1)

         ! Print diagnostic quantities
         CALL vdiagnostic(vplanbc,planfc,vx,vy,vz,t,dt)
         CALL bdiagnostic(bplanbc,planfc,ax,ay,az,t,dt)

!         CALL maxabs(vx,vy,vz,rmp,2)
!         CALL maxabs(ax,ay,az,rmq,0)
!         CALL maxabs(ax,ay,az,rm1,2)
!         CALL maxabs(ax,ay,az,rm2,1)
!
!         IF ( myrank .eq. 0 ) THEN
!            OPEN(1,file='maximum.txt',position='append')
!            WRITE(1,FMT='(5E13.6)') (t-1)*dt,rmp,rmq,rm1,rm2
!            CLOSE(1)
!         ENDIF
