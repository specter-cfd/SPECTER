!=================================================================
! GHOST code: Geophysical High Order Suite for Turbulence
!
! Header file with definitions for conditional compilation.
! Main code is in main3D.fpp. See that file for more details.
!
! 2015 Pablo D. Mininni.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar
!=================================================================

! The lines below define preprocessor variables for each solver.
! New solvers can be created by adding new definitions here, 
! adding the necessary solver files in 'include', and listing
! all objects in 'Makefile'.

! Fluid solvers

#ifdef CHANNEL_SOL
#define INCLUDEFNAME_ 'channel/channel_
#endif

#ifdef RAYBEN_SOL
#define BOUSSINESQ_
#define SCALAR_
#define INCLUDEFNAME_ 'rayben/rayben_
#endif

#ifdef ROTCON_SOL
#define BOUSSINESQ_
#define SCALAR_
#define ROTATION_
#define INCLUDEFNAME_ 'rotcon/rotcon_
#endif

! Do not edit below this line!
! Builds the names of all files to include for each solver
#define STRINGIFY(a) a
#define SOLVERCHECK_ STRINGIFY(INCLUDEFNAME_)validate.f90'
#define GLOBALOUTPUT_ STRINGIFY(INCLUDEFNAME_)global.f90'
#define RKSTEP1_ STRINGIFY(INCLUDEFNAME_)rkstep1.f90'
#define RKSTEP2_ STRINGIFY(INCLUDEFNAME_)rkstep2.f90'
