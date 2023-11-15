!=================================================================
! SPECTER code: Special PEriodic Continuation Turbulence solvER
!
! Header file with definitions for conditional compilation.
! Main code is in specter.fpp. See that file for more details.
!
! 2015 Pablo D. Mininni.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar
!
! 2020 Mauro Fontana. Refactored from GHOST to usage in SPECTER.
!=================================================================

! The lines below define preprocessor variables for each solver.
! New solvers can be created by adding new definitions here, 
! adding the necessary solver files in 'include', and listing
! all objects in 'Makefile'.

! Fluid solvers

#ifdef HD_SOL
#define INCLUDEFNAME_ 'hd/hd_
#endif

#ifdef BOUSS_SOL
#define BOUSSINESQ_
#define SCALAR_
#define INCLUDEFNAME_ 'bouss/bouss_
#endif

#ifdef ROTBOUSS_SOL
#define BOUSSINESQ_
#define SCALAR_
#define ROTATION_
#define INCLUDEFNAME_ 'rotbouss/rotbouss_
#endif

#ifdef MHD_SOL
#define MAGFIELD_
#define INCLUDEFNAME_ 'mhd/mhd_
#endif

#ifdef MHDBOUSS_SOL
#define BOUSSINESQ_
#define SCALAR_
#define MAGFIELD_
#define INCLUDEFNAME_ 'mhdbouss/mhdbouss_
#endif


! Do not edit below this line!
! Builds the names of all files to include for each solver
#define STRINGIFY(a) a
#define SOLVERCHECK_ STRINGIFY(INCLUDEFNAME_)validate.f90'
#define GLOBALOUTPUT_ STRINGIFY(INCLUDEFNAME_)global.f90'
#define RKSTEP1_ STRINGIFY(INCLUDEFNAME_)rkstep1.f90'
#define RKSTEP2_ STRINGIFY(INCLUDEFNAME_)rkstep2.f90'

#ifdef BOUSS_SOL
#define BOUSS_NEWT_ STRINGIFY(INCLUDEFNAME_)newt.f90'
!#define BOUSS_EVOL_T_ STRINGIFY(INCLUDEFNAME_)evol_T.f90'
#endif