Solvers supported by the main code and tools for analysis
=========================================================
The main code is built using CMake (see `build/README`). The solver
is chosen at run time with the `pdes` namelist in `parameter.inp`.
The following solvers are currently supported:

   HD         Hydrodynamic solver

Boundary conditions supported
=============================
- Velocity field: 
    * `noslip`      No slip boundary conditions.

Analysis tools
==============
   BOOTS      Spectral interpolation tool to take binary outputs from
              resolution NXxNYxNZ to NX'xNY'xNZ'.
