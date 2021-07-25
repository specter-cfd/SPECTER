Solvers supported by the main code and tools for analysis
=========================================================
The main code is built using `make` or `make main`. The following 
solvers are currently supported (pass `SOLVER=` to the make command 
to build a solver other than the default):

   HD         Hydrodynamic solver
   BOUSS      Boussinesq solver
   ROTBOUSS   Boussinesq solver with rotation
   MHD        Magnetodydrodynamic solver

Boundary conditions supported
=============================
- Velocity field: 
    * `noslip`      No slip boundary conditions.

- Magnetic field:
    * `conducting`  Perfectly conducting boundary conditions (requires
      stationary walls and `noslip` boundary conditions).
    * `vacuum`      Vacuum boundary conditions.

- Scalar:
    * `isothermal`  Uniform temperature boundary conditions.

Analysis tools
==============
   BOOTS      Spectral interpolation tool to take binary outputs from
              resolution NXxNYxNZ to NX'xNY'xNZ'.
