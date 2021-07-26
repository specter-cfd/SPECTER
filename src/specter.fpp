!=================================================================
      PROGRAM SPECTER
!=================================================================
! SPECTER code: Special PEriodic Continuation Turbulence solvER
!
! SPECTER performs numerical integration for a set of PDEs usually
! found in fluid dynamics. A pseudo-spectral method is used to
! compute spatial derivatives, coupled with an adjustable order
! Runge-Kutta time stepping. When non-periodic directions are
! specified, it utilizes FC-Gram continuations to compute Fourier
! transforms with minimal overhead and considerable accuracy.
! To compile, you need the FFTW and MPI libraries installed on
! your system. More info in the attached README file or in 
! https://github.com/mfontanaar/SPECTER.
!
! Notation: index 'i' moves in the 'x' direction (or kx)
!           index 'j' moves in the 'y' direction (or ky)
!           index 'k' moves in the 'z' direction (or kz)
!
! Conditional compilation options:
!           HD         Hydrodynamic solver
!           BOUSS      Boussinesq solver
!           ROTBOUSS   Boussinesq solver with rotation
!           MHD        Magnetohydrodinamic solver
!
! 2019 Mauro Fontana
!      Department of Physics,
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mfontana@df.uba.ar
!
! Feb 2020: Initial release.
! Jul 2021: Added MHD support.
!
! References:
! - Fontana M, Bruno OP, Mininni PD, Dmitruk P; Comp. Phys. Comm. 256,
!     107482 (2020) -HD,BOUSS-
! - Fontana M, Mininni PD, Bruno OP, Dmitruk P; arXiv: 2107.07077 (2021) -MHD-
! - Rosenberg DL, Mininni PD, Reddy R, Pouquet A.; Atmosphere 11, 178 (2020)
! - Mininni PD, Rosenberg DL, Reddy R, Pouquet A.; P. Comp. 37, 123 (2011)
! - Bruno OP, Lyon M; Jour. Comp. Phys. 229, 2009 (2010)
!
! Definitions for conditional compilation
#include "specter.h"

!
! Modules
      USE fprecision
      USE commtypes
      USE mpivars
      USE filefmt
      USE iovar
      USE grid
      USE fcgram
      USE boundary
      USE fft
      USE var
      USE kes
      USE order
      USE random
      USE threads
      USE gtimer

#if defined(DEF_GHOST_CUDA_)
      USE, INTRINSIC :: iso_c_binding
      USE cuda_bindings
      USE cutypes
#endif

      IMPLICIT NONE

      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: vx,vy,vz
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: fx,fy,fz
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: pr
#ifdef SCALAR_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: th,fs
#endif
#ifdef MAGFIELD_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: ax,ay,az
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: mx,my,mz
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: ph
#endif

!
! Temporal data storage arrays

      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C1,C2,C3
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C4,C5,C6
#ifdef SCALAR_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C7,C8
#endif
#ifdef MAGFIELD_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C9,C10,C11
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C12,C13,C14
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C15,C16,C17
#endif
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:)    :: R1,R2,R3

!
! Variables for temporal storage
      COMPLEX(KIND=GP) :: cdump,jdump
      COMPLEX(KIND=GP) :: cdumq,jdumq
      COMPLEX(KIND=GP) :: cdumr,jdumr
      REAL(KIND=GP)    :: rmp,rmq,rms
      REAL(KIND=GP)    :: rmt,rm1,rm2
      DOUBLE PRECISION :: tmp,tmq,tmr,tms
      DOUBLE PRECISION :: eps,epm
!$    DOUBLE PRECISION, EXTERNAL :: omp_get_wtime


!
! Simulation, forcing and initial conditions parameters
      REAL(KIND=GP)    :: dt,nu,mu
      REAL(KIND=GP)    :: kup,kdn
      REAL(KIND=GP)    :: dump
      REAL(KIND=GP)    :: stat
      REAL(KIND=GP)    :: f0,u0
      REAL(KIND=GP)    :: vxzsta,vxzend,vyzsta,vyzend
      REAL(KIND=GP)    :: fparam0,fparam1,fparam2,fparam3,fparam4
      REAL(KIND=GP)    :: fparam5,fparam6,fparam7,fparam8,fparam9
      REAL(KIND=GP)    :: vparam0,vparam1,vparam2,vparam3,vparam4
      REAL(KIND=GP)    :: vparam5,vparam6,vparam7,vparam8,vparam9
#ifdef BOUSSINESQ_
      REAL(KIND=GP)    :: gama,xmom,xtemp
#endif
#ifdef SCALAR_
      REAL(KIND=GP)    :: kappa
      REAL(KIND=GP)    :: skup,skdn
      REAL(KIND=GP)    :: c0,s0
      REAL(KIND=GP)    :: szsta,szend
      REAL(KIND=GP)    :: sparam0,sparam1,sparam2,sparam3,sparam4
      REAL(KIND=GP)    :: sparam5,sparam6,sparam7,sparam8,sparam9
      REAL(KIND=GP)    :: cparam0,cparam1,cparam2,cparam3,cparam4
      REAL(KIND=GP)    :: cparam5,cparam6,cparam7,cparam8,cparam9
#endif
#ifdef MAGFIELD_
      REAL(KIND=GP)    :: mkup,mkdn
      REAL(KIND=GP)    :: m0,a0
      REAL(KIND=GP)    :: bx0,by0,bz0
      REAL(KIND=GP)    :: mparam0,mparam1,mparam2,mparam3,mparam4
      REAL(KIND=GP)    :: mparam5,mparam6,mparam7,mparam8,mparam9
      REAL(KIND=GP)    :: aparam0,aparam1,aparam2,aparam3,aparam4
      REAL(KIND=GP)    :: aparam5,aparam6,aparam7,aparam8,aparam9
#endif

#ifdef ROTATION_
      REAL(KIND=GP)    :: omegax,omegay,omegaz
#endif

!
! Book keeping variables, options holders and clocks 
      INTEGER :: ini,step
      INTEGER :: tstep,cstep
      INTEGER :: bench,outs
      INTEGER :: seed
      INTEGER :: mult
      INTEGER :: t,o
      INTEGER :: i,j,k
      INTEGER :: ki,kj,kk
      INTEGER :: tind
      INTEGER :: timet,timec
      INTEGER :: ihcpu1,ihcpu2
      INTEGER :: ihomp1,ihomp2
      INTEGER :: ihwtm1,ihwtm2
!$    INTEGER, EXTERNAL     :: omp_get_max_threads

#ifdef MAGFIELD_
      INTEGER :: dyna
#endif

      TYPE(BCPLAN)          :: vplanbc
#ifdef SCALAR_
      TYPE(BCPLAN)          :: splanbc
#endif
#ifdef MAGFIELD_
      TYPE(BCPLAN)          :: bplanbc
#endif


#if defined(DEF_GHOST_CUDA_)
      INTEGER :: idevice, iret, ncuda, ngcuda, ppn
      TYPE(cudaDevicePropG) :: devprop
#endif
      TYPE(IOPLAN)          :: planio
      CHARACTER(len=100)    :: odir,idir,tdir

! Strings to hold the kind of boundary condition
      CHARACTER(len=15)     :: vbcxsta,vbcxend,vbcysta,vbcyend,vbczsta,vbczend
#ifdef SCALAR_
      CHARACTER(len=15)     :: sbcxsta,sbcxend,sbcysta,sbcyend,sbczsta,sbczend
#endif
#ifdef MAGFIELD_
      CHARACTER(len=15)     :: bbcxsta,bbcxend,bbcysta,bbcyend,bbczsta,bbczend
#endif

      LOGICAL               :: bbenchexist


!
! Namelists for the input files
      NAMELIST / status / idir,odir,tdir,stat,mult,bench,outs
      NAMELIST / parameter / dt,step,tstep,cstep,seed
      NAMELIST / boxparams / Lx,Ly,Lz

      NAMELIST / velocity / f0,u0,kdn,kup,nu
      NAMELIST / velocity / fparam0,fparam1,fparam2,fparam3,fparam5
      NAMELIST / velocity / fparam5,fparam6,fparam7,fparam8,fparam9
      NAMELIST / velocity / vparam0,vparam1,vparam2,vparam3,vparam4
      NAMELIST / velocity / vparam5,vparam6,vparam7,vparam8,vparam9

      NAMELIST / velbound / vbcxsta,vbcxend,vbcysta,vbcyend,vbczsta,vbczend
      NAMELIST / velbound / vxzsta,vxzend,vyzsta,vyzend
#ifdef BOUSSINESQ_
      NAMELIST / boussinesq / gama,xmom,xtemp
#endif
#ifdef SCALAR_
      NAMELIST / scalar   / c0,s0,skdn,skup,kappa
      NAMELIST / scalar   / sparam0,sparam1,sparam2,sparam3,sparam4
      NAMELIST / scalar   / sparam5,sparam6,sparam7,sparam8,sparam9
      NAMELIST / scalar   / cparam0,cparam1,cparam2,cparam3,cparam4
      NAMELIST / scalar   / cparam5,cparam6,cparam7,sparam8,sparam9

      NAMELIST / scabound / sbcxsta,sbcxend,sbcysta,sbcyend,sbczsta,sbczend
      NAMELIST / scabound / szsta,szend
#endif
#ifdef MAGFIELD_
      NAMELIST / magfield / m0,a0,mkdn,mkup,mu
      NAMELIST / magfield / mparam0,mparam1,mparam2,mparam3,mparam4
      NAMELIST / magfield / mparam5,mparam6,mparam7,aparam8,aparam9
      NAMELIST / magfield / aparam0,aparam1,aparam2,aparam3,aparam4
      NAMELIST / magfield / aparam5,aparam6,aparam7,aparam8,aparam9

      NAMELIST / magbound / bbcxsta,bbcxend,bbcysta,bbcyend,bbczsta,bbczend
      NAMELIST / dynamo   / dyna
      NAMELIST / uniformb / bx0,by0,bz0
#endif
#ifdef ROTATION_
      NAMELIST / rotation / omegax,omegay,omegaz
#endif

! Initializes the MPI and I/O libraries
      CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED,provided,ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

! NOTE: On systems with a single GPU per node (e.g., Titan)
!       we remove the following block. But on systems with
!       multiple devices per node, this will have to be
!       considered carefully, and possibly re-def'd:
#if defined(DEF_GHOST_CUDA_)
#if defined(CUDA_BIND_LINUX_)
! Initializes CUDA for Linux-based systems. This is a call to an
! NVIDIA-developed intermediary code that gets the GPU dev. no.
! by looking in cpu_info and finding the device that resides on
! its PCI bus:
      iret = cudaGetDeviceCount(ncuda)  ! diagnostic , for now
      CALL MPI_REDUCE(ncuda,ngcuda,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      ppn     = G_PPN_
      idevice = -1
      iret    = setaffinity_for_nvidia(myrank,ppn,idevice)
      iret    = cudaSetDevice(idevice);
#endif
#if defined(CUDA_BIND_DEVICE_)
! Initializes CUDA by selecting device. The list of devices can
! be changed by modifying the env. variable CUDA_VISIBLE_DEVICES:
      iret    = cudaGetDeviceCount(ncuda)
      idevice = mod(myrank,ncuda)
      iret    = cudaSetDevice(idevice);
      IF ( iret .EQ. cudaErrorInvalidDevice ) THEN
        WRITE(*,*)'MAIN: Invalid CUDA device selected: ', &
        idevice, '; myrank=',myrank, '; NUM_CUDA_DEV=',ncuda
        STOP
      ENDIF
      CALL cudaGetDeviceProperties(devprop,idevice)
      IF ( devprop%major .GT. 999 ) THEN
        WRITE(*,*)'MAIN: CUDA device emulation not allowed!'
        STOP
      ENDIF
      IF ( nstreams .GT. 1 .AND. devprop%deviceOverlap .EQ. 0 ) THEN
        WRITE(*,*)'MAIN: Async transfer and computation overlap not supported!'
        STOP
      ENDIF
      iret = cudaGetDevice(idevice)
#endif
#endif

! Create partitions across tasks
     CALL range(1,nx/2+1,nprocs,myrank,ista,iend)
     CALL range(1,nz,nprocs,myrank,ksta,kend)
 ! Initialize array IO
      CALL io_init(myrank,(/nx-Cx,ny-Cy,nz-Cz/),ksta,kend,planio)


!
! Allocates memory for distributed blocks
      ALLOCATE( vx(nz,ny,ista:iend), vy(nz,ny,ista:iend), vz(nz,ny,ista:iend) )
      ALLOCATE( pr(nz,ny,ista:iend) )
      ALLOCATE( fx(nz,ny,ista:iend), fy(nz,ny,ista:iend), fz(nz,ny,ista:iend) )

      ALLOCATE( C1(nz,ny,ista:iend), C2(nz,ny,ista:iend), C3(nz,ny,ista:iend) )
      ALLOCATE( C4(nz,ny,ista:iend), C5(nz,ny,ista:iend), C6(nz,ny,ista:iend) )

      ALLOCATE( R1(nx,ny,ksta:kend), R2(nx,ny,ksta:kend), R3(nx,ny,ksta:kend) )

      ALLOCATE( x(nx), y(ny), z(nz) )
      ALLOCATE( kx(nx), ky(ny), kz(nz) )
      ALLOCATE( kk2(nz,ny,ista:iend), khom(ny,ista:iend) )

#ifdef SCALAR_
      ALLOCATE( th(nz,ny,ista:iend), fs(nz,ny,ista:iend) )
      ALLOCATE( C7(nz,ny,ista:iend), C8(nz,ny,ista:iend) )
#endif

#ifdef MAGFIELD_
      ALLOCATE( ax (nz,ny,ista:iend), ay (nz,ny,ista:iend), az (nz,ny,ista:iend) )
      ALLOCATE( ph(nz,ny,ista:iend) )
      ALLOCATE( mx (nz,ny,ista:iend), my (nz,ny,ista:iend), mz (nz,ny,ista:iend) )

      ALLOCATE( C9 (nz,ny,ista:iend), C10(nz,ny,ista:iend), C11(nz,ny,ista:iend) )
      ALLOCATE( C12(nz,ny,ista:iend), C13(nz,ny,ista:iend), C14(nz,ny,ista:iend) )
      ALLOCATE( C15(nz,ny,ista:iend), C16(nz,ny,ista:iend), C17(nz,ny,ista:iend) )
#endif


!
! The following lines read the file 'parameter.inp'
! Reads general configuration flags from the namelist
! 'status' on the external file 'parameter.inp'
!     idir : directory for unformatted input
!     odir : directory for unformatted output
!     tdir : directory for FC-Gram tables
!     stat : = 0 starts a new run
!            OR  gives the number of the file used to continue a run
!     mult : time step multiplier
!     bench: = 0 production run
!            = 1 benchmark run (no I/O)
!            = 2 higher level benchmark run (+time to create plans)
!     outs : = 0 writes velocity [and vector potential (MAGFIELD_)]
!            = 1 writes vorticity [and magnetic field (MAGFIELD_)]
!            = 2 writes current density (MAGFIELD_)

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=status)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(idir,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(odir,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(tdir,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(stat,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mult,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(bench,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(outs,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

! Reads parameters that will be used to control the
! time integration from the namelist 'parameter' on
! the external file 'parameter.inp'
!     dt   : time step size
!     step : total number of time steps to compute
!     tstep: number of steps between binary output
!     cstep: number of steps between output of global quantities
!     seed : seed for the random number generator

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=parameter)
         CLOSE(1)
         dt = dt/real(mult,kind=GP)
         step = step*mult
         tstep = tstep*mult
         cstep = cstep*mult
      ENDIF
      CALL MPI_BCAST(dt,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(step,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(tstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(seed,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)


! Reads parameters to set the box size
! For periodic diractions the length is in units of 2.pi (=1 gives a
! side of length 2.pi). For non-periodic directions it is measured
! in units of 1 (=1 gives a side length of 1).
!     Lx  : Length in x
!     Ly  : Length in y
!     Lz  : Length in z

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=boxparams)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(Lx,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(Ly,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(Lz,1,GC_REAL,0,MPI_COMM_WORLD,ierr)


! Reads parameters for the velocity field from the 
! namelist 'velocity' on the external file 'parameter.inp' 
!     f0    : amplitude of the mechanical forcing
!     u0    : amplitude of the initial velocity field
!     kdn   : minimum wave number in v/mechanical forcing
!     kup   : maximum wave number in v/mechanical forcing
!     nu    : kinematic viscosity
!     fparam0-9 : ten real numbers to control properties of 
!            the mechanical forcing
!     vparam0-9 : ten real numbers to control properties of
!            the initial conditions for the velocity field

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=velocity)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(f0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(u0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(kdn,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(kup,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nu,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(fparam0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(fparam1,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(fparam2,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(fparam3,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(fparam4,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(fparam5,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(fparam6,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(fparam7,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(fparam8,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(fparam9,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(vparam0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(vparam1,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(vparam2,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(vparam3,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(vparam4,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(vparam5,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(vparam6,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(vparam7,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(vparam8,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(vparam9,1,GC_REAL,0,MPI_COMM_WORLD,ierr)

! Reads parameters for the velocity field boundary conditions
! namelist 'velbound' on the external file 'parameter.inp' 
!     vbcista : boundary condition in the i direction at the
!              start of the domain.
!     vbciend : boundary condition in the i direction at the
!              start of the domain.
!     vizsta  : boundary velocity at z=0 in the i direction
!     vizend  : boundary velocity at z=Lz in the i direction
!
      vbcxsta = ""; vbcxend = ""
      vbcysta = ""; vbcyend = ""
      vbczsta = ""; vbczend = ""

      vxzsta = .0_GP; vxzend = .0_GP
      vyzsta = .0_GP; vyzend = .0_GP

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=velbound)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(vbcxsta,15,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(vbcxend,15,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(vbcysta,15,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(vbcyend,15,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(vbczsta,15,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(vbczend,15,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(vxzsta,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(vxzend,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(vyzsta,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(vyzend,1,GC_REAL,0,MPI_COMM_WORLD,ierr)


#ifdef BOUSSINESQ_
! Reads parameters specifically for Boussinesq solver from the 
! namelist 'boussinesq' on the external file 'parameter.inp'
!     gama  : Bouyancy amplitude (positive definite)
!     xmom  : multiplies bouyancy term in momentum equation
!     xtemp : multiplies temperature-current term in 
!             temperature/density equation
      xmom  = 1.0
      xtemp = 1.0
      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=boussinesq)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(gama ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(xmom ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(xtemp,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      xmom  = xmom  * gama
      xtemp = xtemp * gama
#endif

#ifdef SCALAR_
! Reads parameters for the passive/active scalar from the 
! namelist 'scalar' on the external file 'parameter.inp'
!     s0   : amplitude of the passive scalar source
!     c0   : amplitude of the initial concentration
!     skdn : minimum wave number in concentration/source
!     skup : maximum wave number in concentration/source
!     kappa: diffusivity
!     szsta: value of the scalar at z=Lz
!     szend: value of the scalar at z=0 
!     sparam0-9 : ten real numbers to control properties of 
!            the source
!     cparam0-9 : ten real numbers to control properties of
!            the initial concentration
      szsta = 0.0_GP
      szend = 0.0_GP
      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=scalar)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(s0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(skdn,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(skup,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(kappa,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(szsta,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(szend,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sparam0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sparam1,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sparam2,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sparam3,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sparam4,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sparam5,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sparam6,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sparam7,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sparam8,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sparam9,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cparam0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cparam1,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cparam2,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cparam3,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cparam4,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cparam5,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cparam6,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cparam7,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cparam8,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cparam9,1,GC_REAL,0,MPI_COMM_WORLD,ierr)

! Reads parameters for the scalar field boundary conditions
! namelist 'scabound' on the external file 'parameter.inp' 
!     bbcista : boundary condition in the i direction at the
!              start of the domain.
!     bbciend : boundary condition in the i direction at the
!              start of the domain.
      sbcxsta = ""; sbcxend = ""
      sbcysta = ""; sbcyend = ""
      sbczsta = ""; sbczend = ""

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=scabound)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(sbcxsta,15,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sbcxend,15,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sbcysta,15,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sbcyend,15,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sbczsta,15,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sbczend,15,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
#endif

#ifdef MAGFIELD_
! Reads general configuration flags for runs with 
! magnetic fields from the namelist 'dynamo' on 
! the external file 'parameter.inp'
!     dyna : = 0 when stat=0 generates initial v and B (MAGFIELD_)
!            = 1 when stat.ne.0 imports v and generates B (MAGFIELD_) 

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=dynamo)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(dyna,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

! Reads parameters for the magnetic field from the 
! namelist 'magfield' on the external file 'parameter.inp' 
!     m0   : amplitude of the electromotive forcing
!     a0   : amplitude of the initial vector potential
!     mkdn : minimum wave number in B/electromotive forcing
!     mkup : maximum wave number in B/electromotive forcing
!     mu   : magnetic diffusivity
!     mparam0-9 : ten real numbers to control properties of 
!            the electromotive forcing
!     aparam0-9 : ten real numbers to control properties of
!            the initial conditions for the magnetic field

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=magfield)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(m0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(a0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mkdn,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mkup,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mu,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mparam0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mparam1,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mparam2,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mparam3,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mparam4,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mparam5,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mparam6,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mparam7,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mparam8,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mparam9,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(aparam0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(aparam1,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(aparam2,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(aparam3,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(aparam4,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(aparam5,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(aparam6,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(aparam7,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(aparam8,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(aparam9,1,GC_REAL,0,MPI_COMM_WORLD,ierr)

! Reads parameters for the magnetic field boundary conditions
! namelist 'magbound' on the external file 'parameter.inp' 
!     bbcista : boundary condition in the i direction at the
!              start of the domain.
!     bbciend : boundary condition in the i direction at the
!              start of the domain.
      bbcxsta = ""; bbcxend = ""
      bbcysta = ""; bbcyend = ""
      bbczsta = ""; bbczend = ""

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=magbound)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(bbcxsta,15,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(bbcxend,15,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(bbcysta,15,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(bbcyend,15,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(bbczsta,15,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(bbczend,15,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

! Reads parameters for runs with a uniform magnetic 
! field from the namelist 'uniformb' on the external 
! file 'parameter.inp' 
!     bx0: uniform magnetic field in x
!     by0: uniform magnetic field in y
!     bz0: uniform magnetic field in z
      bx0 = 0.0_GP; by0 = 0.0_GP; bz0 = 0.0_GP
      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=uniformb)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(bx0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(by0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(bz0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
#endif

#ifdef ROTATION_
! Reads parameters for runs with rotation from the 
! namelist 'rotation' on the external file 'parameter.inp'
!     omegax: amplitude of the uniform rotation along x
!     omegay: amplitude of the uniform rotation along y
!     omegaz: amplitude of the uniform rotation along z

      omegax = 0.0_GP; omegay = 0.0_GP; omegaz = 0.0_GP
      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=rotation)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(omegax,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(omegay,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(omegaz,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
#endif


!
! Construct grids
      IF ( Cx .eq. 0 ) THEN
         ! Periodic in x
         rmp = Lx*2.0_GP*pi/nx
         Dkx = 1.0_GP/Lx
         vbcxsta = "periodic"
         vbcxend = "periodic"
#ifdef MAGFIELD_
         bbcxsta = "periodic"
         bbcxend = "periodic"
#endif
      ELSE
         ! Non-periodic in x
         rmp = Lx/(nx-Cx-1)
         Dkx = 2.0_GP*pi/(rmp*nx)
      ENDIF
      DO i=1,nx
        x(i) = rmp*(i-1)
      ENDDO

      IF ( Cy .eq. 0 ) THEN
         ! Periodic in y
         rmp = Ly*2.0_GP*pi/ny
         Dky = 1.0_GP/Ly
         vbcysta = "periodic"
         vbcyend = "periodic"
#ifdef MAGFIELD_
         bbcysta = "periodic"
         bbcyend = "periodic"
#endif
      ELSE
         ! Non-periodic in y
         rmp = Ly/(ny-Cy-1)
         Dky = 2.0_GP*pi/(rmp*ny)
      ENDIF
      DO j=1,ny
        y(j) = rmp*(j-1)
      ENDDO
      IF ( Cz .eq. 0 ) THEN
         ! Periodic in z
         rmp = Lz*2.0_GP*pi/nz
         Dkz = 1.0_GP/Lz
         vbczsta = "periodic"
         vbczend = "periodic"
#ifdef MAGFIELD_
         bbczsta = "periodic"
         bbczend = "periocic"
#endif
      ELSE
        ! Non-periodic in z
        rmp = Lz/(nz-Cz-1)
        Dkz = 2.0_GP*pi/(rmp*nz)
      ENDIF
      DO k=1,nz
        z(k) = rmp*(k-1)
      ENDDO
      pkend = min(nz-Cz,kend)  ! kend of physical domain

      ! Save grid
      IF ( myrank .eq. 0 ) THEN
         OPEN(1,file='x.txt',action='write')
         OPEN(2,file='y.txt',action='write')
         OPEN(3,file='z.txt',action='write')
         DO i = 1,nx-Cx
            WRITE(1,FMT='(1P E23.15)') x(i)
         ENDDO
         DO i = 1,ny-Cy
            WRITE(2,FMT='(1P E23.15)') y(i)
         ENDDO
         DO i = 1,nz-Cz
            WRITE(3,FMT='(1P E23.15)') z(i)
         ENDDO
         CLOSE(1)
         CLOSE(2)
         CLOSE(3)
      ENDIF

! Construct grids for wavenumber domain
      DO i = 1,nx/2
         kx(i) = real(i-1,kind=GP)
         kx(i+nx/2) = real(i-nx/2-1,kind=GP)
      END DO
      DO j = 1,ny/2
         ky(j) = real(j-1,kind=GP)
         ky(j+ny/2) = real(j-ny/2-1,kind=GP)
      END DO
      DO k = 1,nz/2
         kz(k) = real(k-1,kind=GP)
         kz(k+nz/2) = real(k-nz/2-1,kind=GP)
      END DO
      kx = kx*Dkx
      ky = ky*Dky
      kz = kz*Dkz

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               kk2(k,j,i) = kx(i)**2+ky(j)**2+kz(k)**2
            END DO
            khom(j,i) = SQRT(kk2(1,j,i))
         END DO
      END DO

 ! Initializes the FFT library. This must be done at
 ! this stage as it requires the variable "bench" to
 ! be properly initialized.
 ! Use FFTW_ESTIMATE or FFTW_MEASURE in short runs
 ! Use FFTW_PATIENT or FFTW_EXHAUSTIVE in long runs
 ! FFTW 2.x only supports FFTW_ESTIMATE or FFTW_MEASURE

      nth = 1
!$    nth = omp_get_max_threads()
#if !defined(DEF_GHOST_CUDA_)
!$    CALL fftp3d_init_threads(ierr)
#endif
      IF (bench.eq.2) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL GTStart(ihcpu2,GT_CPUTIME)
         CALL GTStart(ihomp2,GT_OMPTIME)
         CALL GTStart(ihwtm2,GT_WTIME)
      ENDIF

      CALL fcgram_create_plan(planfc, (/nx,ny,nz/), (/Cx,Cy,Cz/), &
                              (/ox,oy,oz/), (/x(2),y(2),z(2)/), &
                              tdir, FFTW_ESTIMATE)
      IF (bench.eq.2) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL GTStop(ihcpu2)
         CALL GTStop(ihomp2)
         CALL GTStop(ihwtm2)
      ENDIF

#ifndef TESTS_  /* If building tests jump to after Runge-Kutta*/
! Set up boundary conditions
      CALL setup_bc(vplanbc,planfc, &
          (/ vbcxsta,vbcxend,vbcysta,vbcyend,vbczsta,vbczend /),'v')
#ifdef SCALAR_
      CALL setup_bc(splanbc,planfc, &
          (/ sbcxsta,sbcxend,sbcysta,sbcyend,sbczsta,sbczend /),'s')
#endif
#ifdef MAGFIELD_
      CALL setup_bc(bplanbc,planfc, &
          (/ bbcxsta,bbcxend,bbcysta,bbcyend,bbczsta,bbczend/),'b', &
          vplan=vplanbc)
#endif
!
! Sets the external forcing
      INCLUDE 'initialfv.f90'           ! mechanical forcing
#ifdef SCALAR_
      INCLUDE 'initialfs.f90'           ! scalar source
#endif
#ifdef MAGFIELD_
      INCLUDE 'initialfb.f90'           ! electromotive forcing
#endif

! If stat=0 we start a new run.
! Generates initial conditions for the fields and particles.
 IC : IF (stat.eq.0) THEN
      ini = 1
      tind = 0                          ! index for the binaries
      timet = tstep
      timec = cstep
      INCLUDE 'initialv.f90'            ! initial velocity
      ! Zero pressure gradient
      DO i=ista,iend
      DO j=1,ny
      DO k=1,nz-Cz
         pr(k,j,i)  = 0._GP
      ENDDO
      ENDDO
      ENDDO
#ifdef SCALAR_
      INCLUDE 'initials.f90'            ! initial scalar density
#endif
#ifdef MAGFIELD_
      INCLUDE 'initialb.f90'            ! initial vector potential
      ! Zero electric potential
      DO i=ista,iend
      DO j=1,ny
      DO k=1,nz-Cz
         ph(k,j,i)  = 0._GP
      ENDDO
      ENDDO
      ENDDO

#endif
      ELSE
! If stat.ne.0 a previous run is continued
      ini = int((stat-1)*tstep) + 1
      tind = int(stat)
      WRITE(ext, fmtext) tind
      timet = 0
      timec = int(modulo(float(ini-1),float(cstep)))

      CALL io_read(1,idir,'vx',ext,planio,R1)
      CALL io_read(1,idir,'vy',ext,planio,R2)
      CALL io_read(1,idir,'vz',ext,planio,R3)
      CALL fftp3d_real_to_complex(planfc,R1,vx,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planfc,R2,vy,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planfc,R3,vz,MPI_COMM_WORLD)

      CALL io_read(1,idir,'pr',ext,planio,R1)
      CALL fftp2d_real_to_complex_xy(planfc,R1,pr,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz-Cz
               pr(k,j,i) = pr(k,j,i)*dt   ! Convert to p'
             END DO
          END DO
      END DO

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

      IF (stat .eq. 1) THEN
         CALL fftp3d_complex_to_real(planfc,C1,R1,MPI_COMM_WORLD)
         CALL fftp3d_complex_to_real(planfc,C2,R2,MPI_COMM_WORLD)
         CALL fftp3d_complex_to_real(planfc,C3,R3,MPI_COMM_WORLD)
         CALL io_write(1,odir,'vx','0000',planio,R1)
         CALL io_write(1,odir,'vy','0000',planio,R2)
         CALL io_write(1,odir,'vz','0000',planio,R3)
      ENDIF


#ifdef SCALAR_
      CALL io_read(1,idir,'th',ext,planio,R1)
      CALL fftp3d_real_to_complex(planfc,R1,th,MPI_COMM_WORLD)
#endif
#ifdef MAGFIELD_
 DYN: IF (dyna.eq.0) THEN
         CALL io_read(1,idir,'ax',ext,planio,R1)
         CALL io_read(1,idir,'ay',ext,planio,R2)
         CALL io_read(1,idir,'az',ext,planio,R3)
         CALL fftp3d_real_to_complex(planfc,R1,ax,MPI_COMM_WORLD)
         CALL fftp3d_real_to_complex(planfc,R2,ay,MPI_COMM_WORLD)
         CALL fftp3d_real_to_complex(planfc,R3,az,MPI_COMM_WORLD)

         CALL io_read(1,idir,'ph',ext,planio,R1)
         CALL fftp2d_real_to_complex_xy(planfc,R1,ph,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz-Cz
                  ph(k,j,i) = ph(k,j,i)*dt   ! Convert to p'
                END DO
             END DO
         END DO
      ELSE
         INCLUDE 'initialb.f90'      ! initial vector potential
         ! Zero electric potential
         DO i=ista,iend
         DO j=1,ny
         DO k=1,nz-Cz
            ph(k,j,i)  = 0._GP
         ENDDO
         ENDDO
         ENDDO

         ini = 1                     ! resets all counters (the
         tind = 0                    ! dynamo run starts at t=0)
         timet = tstep
         timec = cstep
      ENDIF DYN
#endif
      ENDIF IC

!
! Time integration scheme starts here.
! Does ord iterations of Runge-Kutta. If 
! we are doing a benchmark, we measure 
! cputime before starting.

      IF (bench.eq.1) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL GTStart(ihcpu1,GT_CPUTIME)
         CALL GTStart(ihomp1,GT_OMPTIME)
         CALL GTStart(ihwtm1,GT_WTIME)
! Must reset FFT internal timers after initialization:
         comtime = 0.0D0
         ffttime = 0.0D0
         tratime = 0.0D0
#if defined(DEF_GHOST_CUDA_)
         memtime = 0.0D0
         tottime = 0.0D0
#endif

      ENDIF

 RK : DO t = ini,step
! Every 'tstep' steps, stores the fields 
! in binary files
 BIN:    IF ((timet.eq.tstep).and.(bench.eq.0)) THEN
            timet = 0
            tind = tind+1
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
            IF (outs.ge.1) THEN
               CALL curlk(C2,C3,C4,1)
               CALL curlk(C1,C3,C5,2)
               CALL curlk(C1,C2,C6,3)
               CALL fftp3d_complex_to_real(planfc,C4,R1,MPI_COMM_WORLD)
               CALL fftp3d_complex_to_real(planfc,C5,R2,MPI_COMM_WORLD)
               CALL fftp3d_complex_to_real(planfc,C6,R3,MPI_COMM_WORLD)
               CALL io_write(1,odir,'wx',ext,planio,R1)
               CALL io_write(1,odir,'wy',ext,planio,R2)
               CALL io_write(1,odir,'wz',ext,planio,R3)
            ENDIF
            CALL fftp3d_complex_to_real(planfc,C1,R1,MPI_COMM_WORLD)
            CALL fftp3d_complex_to_real(planfc,C2,R2,MPI_COMM_WORLD)
            CALL fftp3d_complex_to_real(planfc,C3,R3,MPI_COMM_WORLD)
            CALL io_write(1,odir,'vx',ext,planio,R1)
            CALL io_write(1,odir,'vy',ext,planio,R2)
            CALL io_write(1,odir,'vz',ext,planio,R3)

            ! Pressure. Transform from p' to p to save
            rmp = 1.0_GP/ &
	          (real(nx,kind=GP)*real(ny,kind=GP)*dt)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
               DO j = 1,ny
                  DO k = 1,nz
                     C1(k,j,i) = pr(k,j,i)*rmp
                  END DO
               END DO
            END DO
            CALL fftp2d_complex_to_real_xy(planfc,C1,R1,MPI_COMM_WORLD)
            CALL io_write(1,odir,'pr',ext,planio,R1)

#ifdef SCALAR_
            rmp = 1.0_GP/ &
	          (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
               DO j = 1,ny
                  DO k = 1,nz
                     C1(k,j,i) = th(k,j,i)*rmp
                  END DO
               END DO
            END DO
            CALL fftp3d_complex_to_real(planfc,C1,R1,MPI_COMM_WORLD)
            CALL io_write(1,odir,'th',ext,planio,R1)
#endif
#ifdef MAGFIELD_
            rmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
               DO j = 1,ny
                  DO k = 1,nz
                     C1(k,j,i) = ax(k,j,i)*rmp
                     C2(k,j,i) = ay(k,j,i)*rmp
                     C3(k,j,i) = az(k,j,i)*rmp
                  END DO
               END DO
            END DO
            IF (outs.ge.1) THEN
               CALL curlk(C2,C3,C4,1)
               CALL curlk(C1,C3,C5,2)
               CALL curlk(C1,C2,C6,3)
               CALL fftp3d_complex_to_real(planfc,C4,R1,MPI_COMM_WORLD)
               CALL fftp3d_complex_to_real(planfc,C5,R2,MPI_COMM_WORLD)
               CALL fftp3d_complex_to_real(planfc,C6,R3,MPI_COMM_WORLD)
               CALL io_write(1,odir,'bx',ext,planio,R1)
               CALL io_write(1,odir,'by',ext,planio,R2)
               CALL io_write(1,odir,'bz',ext,planio,R3)
            ENDIF
            IF (outs.eq.2) THEN
               CALL laplak(C1,C4)
               CALL laplak(C2,C5)
               CALL laplak(C3,C6)
               CALL fftp3d_complex_to_real(planfc,C4,R1,MPI_COMM_WORLD)
               CALL fftp3d_complex_to_real(planfc,C5,R2,MPI_COMM_WORLD)
               CALL fftp3d_complex_to_real(planfc,C6,R3,MPI_COMM_WORLD)
               CALL io_write(1,odir,'jx',ext,planio,R1)
               CALL io_write(1,odir,'jy',ext,planio,R2)
               CALL io_write(1,odir,'jz',ext,planio,R3)
            ENDIF
            CALL fftp3d_complex_to_real(planfc,C1,R1,MPI_COMM_WORLD)
            CALL fftp3d_complex_to_real(planfc,C2,R2,MPI_COMM_WORLD)
            CALL fftp3d_complex_to_real(planfc,C3,R3,MPI_COMM_WORLD)
            CALL io_write(1,odir,'ax',ext,planio,R1)
            CALL io_write(1,odir,'ay',ext,planio,R2)
            CALL io_write(1,odir,'az',ext,planio,R3)

            ! Electric potential. Transform from phi' to phi to save
            rmp = 1.0_GP/ &
	          (real(nx,kind=GP)*real(ny,kind=GP)*dt)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
               DO j = 1,ny
                  DO k = 1,nz
                     C1(k,j,i) = ph(k,j,i)*rmp
                  END DO
               END DO
            END DO
            CALL fftp2d_complex_to_real_xy(planfc,C1,R1,MPI_COMM_WORLD)
            CALL io_write(1,odir,'ph',ext,planio,R1)

#endif
         ENDIF BIN

! Every 'cstep' steps, generates external files 
! with global quantities.

         IF ((timec.eq.cstep).and.(bench.eq.0)) THEN
            timec = 0
            INCLUDE GLOBALOUTPUT_
         ENDIF

! Runge-Kutta step 1
! Copies the fields into auxiliary arrays

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
         DO k = 1,nz

            INCLUDE RKSTEP1_

         END DO
         END DO
         END DO

! Runge-Kutta step 2
! Evolves the system in time

         DO o = ord,1,-1

         INCLUDE RKSTEP2_
      
         END DO

         timet = timet+1
         timec = timec+1
      END DO RK
!
! End of Runge-Kutta
#else
!      INCLUDE 'tests/playground.f90'
!      INCLUDE 'tests/fft.f90'
      INCLUDE 'tests/fc_dirichlet.f90'
      INCLUDE 'tests/fc_neumann.f90'
      INCLUDE 'tests/fc_neumann2.f90'
#ifdef MAGFIELD_
      INCLUDE 'tests/fc_robin.f90'
#endif
      INCLUDE 'tests/energy.f90'
      INCLUDE 'tests/poisson.f90'
#endif /*End of if TESTS_*/

! Computes the benchmark
      IF (bench.gt.0) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL GTStop(ihcpu1)
         CALL GTStop(ihomp1)
         CALL GTStop(ihwtm1)
         IF (myrank.eq.0) THEN
            INQUIRE( file='benchmark.txt', exist=bbenchexist )
            OPEN(1,file='benchmark.txt',position='append')
#if defined(DEF_GHOST_CUDA_)
            IF ( .NOT. bbenchexist ) THEN
               WRITE(1,*) &
	       '# nx ny nz nsteps nprocs nth nstrm TCPU TOMP TWTIME TFFT TTRA',&
           'TCOM TMEM TASS TCONT TNEU TROB TTOT'
            ENDIF
            WRITE(1,*) nx,ny,nz,(step-ini+1),nprocs,nth, &
                       nstreams                        , &
                       GTGetTime(ihcpu1)/(step-ini+1)  , &
                       GTGetTime(ihomp1)/(step-ini+1)  , &
                       GTGetTime(ihwtm1)/(step-ini+1)  , &
                       ffttime/(step-ini+1), tratime/(step-ini+1), &
                       comtime/(step-ini+1), memtime/(step-ini+1), &
                       conttime/(step-ini+1), neutime/(step-ini+1),&
                       robtime/(step-ini+1), tottime/(step-ini+1)
#else
            IF ( .NOT. bbenchexist ) THEN
               WRITE(1,*) &
	       '# nx ny nz nsteps nprocs nth TCPU TOMP TWTIME TFFT TTRA TCOM',&
           'TCONT TNEU TROB TTOT'
            ENDIF
            WRITE(1,*) nx,ny,nz,(step-ini+1),nprocs,nth, &
                       GTGetTime(ihcpu1)/(step-ini+1),   &
                       GTGetTime(ihomp1)/(step-ini+1),   &
                       GTGetTime(ihwtm1)/(step-ini+1),   &
                       ffttime/(step-ini+1), tratime/(step-ini+1), &
                       comtime/(step-ini+1), conttime/(step-ini+1),&
                       neutime/(step-ini+1), robtime/(step-ini+1), &
                       tottime/(step-ini+1)
#endif
            IF (bench.eq.2) THEN
               WRITE(1,*) 'FFTW: Create_plan = ',      &
                       GTGetTime(ihcpu2)/(step-ini+1), &
                       GTGetTime(ihomp2)/(step-ini+1), &
                       GTGetTime(ihwtm2)/(step-ini+1)
            ENDIF
            CLOSE(1)
         ENDIF
      ENDIF

      CALL MPI_FINALIZE(ierr)

      CALL GTFree(ihcpu1)
      CALL GTFree(ihomp1)
      CALL GTFree(ihwtm1)
      CALL GTFree(ihcpu2)
      CALL GTFree(ihomp2)
      CALL GTFree(ihwtm2)

      CALL fcgram_destroy_plan(planfc)

      DEALLOCATE( vx,vy,vz )
      DEALLOCATE( pr )
      DEALLOCATE( fx,fy,fz )

      DEALLOCATE( C1,C2,C3 )
      DEALLOCATE( C4,C5,C6 )

      DEALLOCATE( R1,R2,R3 )

      DEALLOCATE( x,y,z )
      DEALLOCATE( kx,ky,kz )
      DEALLOCATE( kk2,khom )

#ifdef SCALAR_
      DEALLOCATE( th,fs )
      DEALLOCATE( C7,C8 )
#endif
#ifdef MAGFIELD_
      DEALLOCATE(  ax,  ay, az  )
      DEALLOCATE(  ph )
      DEALLOCATE(  mx,  my, mz  )
      DEALLOCATE(  C9, C10, C11 )
      DEALLOCATE( C12, C13, C14 )
      DEALLOCATE( C15, C16, C17 )
#endif

      END PROGRAM SPECTER
