SPECTER code: Special PEriodic Continuation Turbulence solvER
==============================================================

The code and different solvers provided here numerically integrate
fluids equations in 3 dimensions with non-periodic boundary conditions.
The solver is chosen at run time; for now the only one is the
hydrodynamic solver `HD`. A 
pseudo-spectral method is used to compute spatial derivatives, 
while a Runge-Kutta method of adjustable order is used to evolve the 
system in the time domain. For non-periodic boundary conditions, an 
FC-Gram continuation is performed before transforming real fields. To 
compile, you need the FFTW library installed on your system and some 
flavor of MPI. Note that only serial transforms in the FFTW package 
are used. The parallel part of the FFT is handled by SPECTER explicitly.

1 - Installing FFTW
-------------------
The parallel FFT based on `FFTW` is in the directory `fftp`. Only `FFTW`
versions `>= 3.0` are supported, as the API of `FFTW 2.x` is incompatible
with that of `FFTW 3.x`.

To compile `FFTW 3.x` (freely available at http://www.fftw.org/), `cd` to
the directory containing the source and run

```
./configure
make
make install
make distclean
./configure --enable-float
make
make install
```

The first pass installs the double precision libraries, and the
second pass installs the single precision libraries. By default,
the libraries are installed in `/usr/local/lib`, `/usr/local/man`,
etc. You can specify an installation prefix other than `/usr/local`
(e.g. the user's home directory) by giving `configure` the option
`--prefix=PATH`.

2 - Compiling the code
----------------------
`SPECTER` is compiled with CMake, out of the source tree. The options
are set in the file `UserConfig.cmake`, or passed to `cmake` with `-D`
(the values passed with `-D` are kept in the cache of the build
directory). From the root of the repository:

```
mkdir build && cd build
cmake .. -DFFTW3_ROOT=/path/to/fftw -DPRECISION=SINGLE
make
```

The binary is placed in `bin/` and named like the solver (`HD`),
together with `BOOTS` (section 5).

### 2.1. Code options:

The variable `PRECISION` (`SINGLE` or `DOUBLE`) sets the precision
of the floating point numbers, and `FFTW3_ROOT` points to the
installation of FFTW.

The grid, the FC-Gram continuation, the solver, the initial
conditions, the forcing and the order of the Runge-Kutta method are
chosen at run time (see section 3.1).

The remaining options are optional and the defaults give good
performance in most computers (don't change the values if you are
unsure about the meaning of these variables). The variables

```
IKIND  = 8
CSIZE  = 32
NSTRIP = 1
```

set the variables `ikind`, `csize`, and `nstrip` used in the
parallel FFTs. `IKIND` is the size of C pointers. For optimal
results, it should be 4 in 32-bits machines, and 8 in 64-bits
machines. `CSIZE` is used by the parallel FFT to do cache-friendly
transpositions. The optimal value for a particular machine can be
adjusted doing benchmarks, but rule of thumb values are 8 if the L1
cache is smaller or equal to 64 kb, and 16 if the L1 cache is
larger than 64 kb (24 or 32 may be best in some processors with
large cache). The variable `NSTRIP` controls strip mining during
the transpositions. The default value is 1 and in that case no
strip mining is done.

Finally, if required, the users can change the number of digits
used to label the binary files. The variables in the module
`filefmt` (in the file `pseudo/pseudospec_mod.f90`) control the
length of the time label used to name these files.

### 2.2. Parallelization:

`SPECTER` uses 1D domain decomposition with pure MPI parallelization.
The options for hybrid OpenMP-MPI parallelization and for GPUs in
`UserConfig.cmake` come from `GHOST` and have no effect yet.

### 2.3. Structure of the code:

As in `GHOST`, the main program `specter.F90` builds the solver, the
initial conditions, the forcing and the time stepper with their
factories, from the lists in the input file:

- `pdes/`     : the solver (`hd_solver.f90`) and its factory
- `ic/`       : the initial conditions
- `force/`    : the forcing
- `steppers/` : the time stepper
- `boundary/` : the boundary conditions of the velocity

The pseudo-spectral routines are in `pseudo/`, and the parallel FFT
with the FC-Gram continuation in `fftp/`.


3 - Running the code
--------------------
The codes can be executed using `mpiexec` or the correspondent vendor
equivalent. The first node then reads the parameters for the run
from the input file, and passes the values to the rest of the
nodes.

### 3.1. Input file:

The input file for all the solvers is named `parameter.inp`. The
file is separated into several `lists`, which have the following
syntax:

```
! General comments
&listname1
variable1 = 1.00 ! Comment
variable2 = 2.00 ! Another comment
/

! More comments
&listname2
variable3 = 3.00 ! Comment
variable4 = 4.00 ! Another comment
/
...
```

The number of lists required depends on the solver. An example with
the list of all variables for all the solvers is given in the
directory `examples`. Note all variables are not needed for all the
solvers. The order of the lists in the file `parameter.inp`, and
the order of the variables inside each list, are not important. All
solvers need a list 'grid' with the grid and the FC-Gram
continuation, as in `GHOST`:

```
&grid
nx = 256             ! grid points in x
ny = 128             ! grid points in y
nz = 128             ! grid points in z
Cz = 25              ! continuation points in z
oz = 5               ! order at the boundaries in z
/
```

The variables in this list are
- `nx`, `ny`, `nz`: number of grid points in each direction,
    including the continuation points.
- `Cx`, `Cy`, `Cz`: number of continuation points in each
    direction. The physical domain has `nx-Cx` points in x, and
    `Cx = 0` (the default) makes x a periodic direction.
- `ox`, `oy`, `oz`: order of the approximation at the boundaries in
    each non-periodic direction (0, the default, in periodic
    directions).

The FC-Gram tables for the chosen `C` and `o` must exist in `tdir`
(see section 4). All solvers also need a list 'status' with variables

```
&status
idir = "/ptmp/run"   ! read files from idir
odir = "/ptmp/run"   ! write files to odir
tdir = "/ptmp/run"   ! load FC-Gram tables from tdir
stat = 0             ! status
mult = 1             ! timestep multiplier
bench = 0            ! benchmark
outs = 0             ! output
dt = 1.0e-3          ! time step
step = 10001         ! number of steps
tstep = 500          ! binary output
cstep = 10           ! global quantities
seed = 1000          ! random seed
/
```

The variables in this list are
- `idir`  : directory for binary input
- `odir`  : directory for binary output
- `tdir`  : directory for the FC-Gram tables
- `status`: the number of the last binary file available if a
        simulation is restarted. If zero, a new simulation is done.
- `mult`  : time step multiplier. Useful for restarting a run with
    a different time step but keeping the same cadence when saving
    to files. mult = 2 results in an actual timestep of half what 
    is defined in the dt parameter.
- `benchmark`: if 1, performs a benchmark run.
- `output`: controls the amount of output written in the binary files.
- `dt`  : sets the time step, which must satisfy the CFL condition
      `dt <= d/u`, where `u` is a characteristic speed, and `d` is
      the smallest of `dx,dy,dz` (the grid spacings).
- `step`: controls the total number of time steps the code will do.
- `tstep` and `cstep`: number of time steps done between exporting
      global quantities (`cstep`), and binary files (`tstep`).
- `seed`: seed for the random number generator.

# 
The solver is chosen with the list 'pdes', as in `GHOST`:

```
&pdes
solver = "HD"     ! the only solver available
/
```

# 
The HD solver reads its parameters from the list 'HD', as in
`GHOST`:

```
&HD
nu = 1.e-3        ! kinematic viscosity
/
```

where
- `nu`   : kinematic viscosity.
- `fidir`, `fodir` and `todir` (optional): directories for the binary
    input, the binary output and the text output of the solver. By
    default they are `idir`, `odir` and the current directory.

The initial conditions and the forcing have their own lists (see
below); the list 'velocity' of older versions is no longer read.

# 
The time stepper is set by the list 'stepper', as in `GHOST`:

```
&stepper
sname = "TRADITIONAL"   ! Runge-Kutta of Canuto et al.
itype = 1
norder = 2              ! order of the Runge-Kutta method
nstage = 2              ! number of stages (= norder)
/
```

The only stepper available is `TRADITIONAL`, a low storage
Runge-Kutta with `norder` stages.

# 
The initial velocity field is set by the list 'initial_conditions':

```
&initial_conditions
iclist = "cbnoise_v;couette_v"
/
```

`iclist` holds one or more initial conditions separated by `;`,
applied in order (see `ic/ic_velocity.f90`):
- `read_v`: reads vx, vy and vz from `idir`, numbered by `stat`
    (to continue a run).
- `null_v`: null velocity field.
- `cbnoise_v`: incompressible noise that satisfies no-slip at the
    walls (Clever & Busse), with parameters in the list `cbnoise_v`
    (`u0`, `kdn`, `kup`, `nmin`, `nmax` and `rphase`).
- `couette_v`: adds to vx the laminar Couette profile between the
    wall velocities `vxzsta` and `vxzend` (list `couette_v`).
- `poiseuille_v`: adds to vx the laminar Poiseuille profile of the
    forcing `f0` with viscosity `nu` (list `poiseuille_v`).

The profiles add to the previous field, so they go after another
initial condition: `"cbnoise_v;couette_v"` is noise plus Couette,
and `"null_v;couette_v"` the laminar flow alone.

# 
The mechanical forcing is set by the list 'forcing':

```
&forcing
forces = "uniform_fv,constant_fv"
cort = 1000.
/
```

`forces` holds one or more pairs `forcing,update` separated by `;`
(see `force/force_velocity.f90`). The forcings are
- `null_fv`: null forcing.
- `uniform_fv`: uniform forcing `f0` along x, the mean pressure
    gradient of a Poiseuille flow (list `uniform_fv`, with `f0`).

and the only update method is `constant_fv` (forcing constant in
time). `cort` is the correlation time of the update methods.

# 
And then the boundary conditions are specified with

```
&velbound
vbczsta = "noslip"
vbczend = "noslip"
vxzsta  = 0.0
vyzsta  = 0.0
vxzend  = 1.0
vyzend  = 1.0
/
```

where
- `vbczsta` and `vbczend` : kind of boundary conditions at `z=0` (`zsta`)
     and `z=Lz` (`zend`). No-slip in the example.
- `vizsta`and `vizend`    : average value of the `i`-th component of the
     velocity field at `z=0` (`zsta`) and `z=Lz` (`zend`).

The only boundary condition available in a non-periodic direction is
`noslip`; directions without continuation points are periodic. See
the file `examples/parameter.inp` for a complete input file.

### 3.2 Output files:

When the code is not running in benchmark mode, it writes several
files with global quantities and components of the fields.
Global quantities is written to text files in the directory where
the binary was executed, while binary files are written to the
directory defined in the `parameter.inp` file.

Binary files are named by the field name, component, and a number
that labels the time of the snapshot of the field. As an example,
the x-component of the velocity field in a hydrodynamic simulation
is stored as `vx.NNN.out`, where NNN is proportional to the time at
which the snapshot was taken. The actual time can be computed as
`t = dt*tstep*(NNN-1)`.

The code also writes text files with quantities integrated over all
volume (`helicity.txt`, `balance.txt`, etc.).  In a hydrodynamic run,
the file `balance.txt` has four columns with the time (`t`), two times
the energy (`2.E`), two times the enstrophy (`2.Omega`), and the energy
injection rate (`eps`).

The binary files of the HD solver are the velocity (`vx`, `vy`, `vz`),
the vorticity (`wx`, `wy`, `wz`, if `outs` is 1 or larger) and the
pressure (`pr`). When a run is continued with `stat=1`, the velocity
read is also written with the number `0000`. The file
`noslip_diagnostic.txt` checks the divergence and the boundary
conditions. More details are in the comments of `pdes/hd_solver.f90`.


4 - FC-Gram continuations
-------------------------
FC-Gram continuations employ a series of precomputed tables to compute
a periodic extension of any given data and remove Gibbs ringing in the
Fourier transform. These tables (and, hence, the continuation itself)
depend on two parameters, the length of the extension in number of
points (`C`) and the order of approximation at the boudary (`O`). The
values of `C` and `O` for each direction are read at run time from
the list 'grid' (section 3.1), and the appropriate tables must exist
in the tdir directory, named `AC-O.dat` and `QO.dat`, with `C` and
`O` having the meaning stated above. A set of frequently
used tables can be found in the tables directory, together with a
`README.info` file for more information. Note that for boundary
conditions which involve the i-th derivative the tables `QinO.dat` are
required.


5 - Regridding (BOOTS)
----------------------
Together with the solver, the build generates `BOOTS`, the 'bootstrap
regridding' procedure of `GHOST`. It reads fields in physical space
on an old grid and writes them on a larger grid, prolongating them
in Fourier space; the output files are the input files with the new
size appended (`_PNNNNN-NNNNN-NNNNN`, with the physical points). The
new grid is read from the list 'grid' of `boots.inp` (section 3.1),
and the old grid and the files to process from the list 'regrid'.
In the non-periodic direction the old data is extended with an
FC-Gram continuation of `Czt` points, chosen so that the extended
periods of both grids match, and of the order `ozt` given in the
list 'order' (if zero, the order of the new grid). The tables for
that continuation must exist in `tdir`; the set in `tables/boots`
covers some common cases. If the grid changes in x, the number of
MPI tasks must be 1 or of the form 2^i nx_new/nx_old with i > 0. An
example input file can be found in `examples/boots.inp`.

6 - References
--------------
1. The Navier-Stokes solver is explained in:  Fontana, M., Bruno, O. P.,
Mininni P. D. & Dmitruk P.; Fourier continuation method for
incompressible fluids with boundaries. Computer Physics Comm. 256,
107482 (2020).
https://doi.org/10.3390/atmos11020178.

2. The MHD solver is detailed in: Fontana M., Mininni P. D., Bruno O. P.
& Dmitruk P.; Vector potential-based MHD solver for non-periodic flows
using Fourier continuation expansions. (2021).
https://arxiv.org/abs/2107.07077

3. The FC-Gram technique is explained in: Bruno, O. P., & Lyon, M.
(2010). High-order unconditionally stable FC-AD solvers for general
smooth domains I. Basic elements. Journal of Computational Physics,
229(6), 2009–2033 (2009).
https://doi.org/10.1016/j.jcp.2009.11.020.

4. The hybrid parallelization scheme is described in: Mininni, P. D.,
Rosenberg, D., Reddy, R., & Pouquet, A. (2011). A hybrid MPI-OpenMP
scheme for scalable parallel pseudospectral computations for fluid
turbulence. Parallel Computing, 37(6–7), 316–326 (2011).
https://doi.org/10.1016/j.parco.2011.05.004.

5. CUDA paralellization is described in: Rosenberg, D., Mininni, P. D.,
Reddy, R., & Pouquet, A. (2018). GPU parallelization of a hybrid
pseudospectral fluid turbulence framework using CUDA (2020).
https://doi.org/10.3390/atmos11020178.

6 - Citing `SPECTER`
--------------------
We ask users to cite this paper when referencing `SPECTER`.
- Fontana, M., Bruno, O. P., Mininni P. D. & Dmitruk P.; Fourier
continuation method for incompressible fluids with boundaries.
https://doi.org/10.1016/j.cpc.2020.107482.
