SPECTER code: Special PEriodic Continuation Turbulence solvER
==============================================================

The code and different solvers provided here numerically integrate
fluids equations in 3 dimensions with non-periodic boundary conditions.
Different solvers are provided (check the Makefile and the comments in
`specter.fpp` for the current list of supported solvers). A 
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
### 2.1. Setting code options:

The linear resolution is set before compilation in the file
Makefile.in. The lines

```
NX       = 64
NY       = 64
NZ       = 64
```

set the number of grid points (including continuation regions) in
each direction. Then the lines

```
CX       = 0
CY       = 0
CZ       = 25
```

set the number of continuation points to use in each direction.
For periodic directions the desired variable should be set to zero
(i.e. CX=0 implies that X is periodic). Similarly, the section

```
OX       = 0
OY       = 0
OZ       = 5
```

defines the order of approximation in each non-periodic direction.
Should be set to zero for periodic directions.

Also, the variable

```
ORD      = 2
```

controls the number of iterations in the Runge-Kutta method.

The set of PDEs and boundary conditions are defined by the variable

```
SOLVER   = BOUSS
```

See the file `SOLVERS_AND_BOUNDARY_CONDITIONS.info` for a list of the
implemented solvers.

Finally the variable

```
PRECISION = SINGLE
```

adjust the usage of single or double precision for the representaton
of floating point numbers in the code.

The remaining options in this section are optional and the defaults
give good performance in most computers. But if needed, you can set
three more variables for optimization purposes in the file
Makefile.in (again, the default values work fine in most computers,
don't change the values if you are unsure about the meaning of
these variables). The lines

```
IKIND    = 8
CSIZE    = 32
NSTRIP   = 1
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
used to label spectrum, transfer and binary files. The variables
in the module `filefmt` (in the file `pseudo/pseudospec3D_mod.fpp`)
control the length of the time label used to name these files.

### 2.2 Compiling the code:

Initial conditions for the fields and the external forcings are
defined in the files `initial*.f90`. Examples of these files can
be found in the directory `examples`.

Note that the name of all files in the directory `examples` finish
with `_example-name`. To use any of these files, you must copy the
file to the `src` directory and remove the `_example-name` from
the name of the file. As an example, to use null initial
conditions for the velocity field, from the `src` directory you
must do

```
cp examples/initialv.f90_null initialv.f90
```

Check the Makefile to see all paths and variables are correct. This
includes the variables explained above as well as the relevant
compilers, MPI libraries and their respective flags.

The path to the FFTW libraries can be changed by editting the
variable `FFTWDIR` in `Makefile.in`. In addition, the Makefile
supports `make clean`, `make dist`, and extra options to make code
for data analysis after the runs are finished.

### 2.3 Parallelization models

By default, `SPECTER` uses 1D domain decomposition with pure MPI
parallelization. This behavior can be changed with two settings
in the Makefile.in file, or equivalently with two variables that
can be passed to make at compile time.

#### 2.3.1 Hybrid parallelization

The code provides support for OpenMP-MPI hybridization, for use in
supercomputers with multi-core processors and for simulations using
a large number of cores. When the variable

```
P_HYBRID=yes
```

is defined in Makefile.in, a hybrid parallelization model is used.
Then, the number of MPI jobs and the number of threads in each job
can be set independently at run time. To define the number
of threads (in this example, set to 4 threads) in a Bourne-like
shell set the environment variable

```
export OMP_NUM_THREADS=4
```

and in a C-like shell do

```
setenv OMP_NUM_THREADS 4
```

For the hybrid code to scale properly, processor affinity and/or
binding are crucial. How to set processor affinity and binding is
platform dependent; check the platform documentation where you are
running for more information.

#### 2.3.2 CUDA support

SPECTER has experimental support for GPU-based computation of FFTs
using CUDA. To use this option, the NVIDIA CUDA Toolkit (the CUDA
compilers plus the GPU-accelerated math libraries) must be
installed in the system. Paths to the CUDA compiler and libraries
must be declared in the file Makefile.in (check all the variables
defined in the section "CUDA compilers" of Makefile.in).

To enable this option, pass the following variable to make at
compile time:

```
P_CUDA=yes
```

Only one GPU can be used per MPI task. In systems with multiple
GPUs, the GPU binded to a given MPI task should be the one with
the fastest access to the CPU in which the MPI task is running.
How to set affinity between CPUs and GPUs is platform dependent;
check the platform documentation where you are running for more
information.

### 2.4 About the solvers:
The main program `specter.fpp` includes files for the different
solvers at compile time from the directory `include`. The files
in `include` are named `solver/solver_component.f90` where `solver`
is the name of the solver (in lowercase) and `component` indicates
what action is done by the file. Most common components are

`global`   : contains routines to compute global quantities
`rkstep1`  : first step of Runge-Kutta
`rkstep2`  : second step of Runge-Kutta (contains the equations)


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
solvers need at least a list 'status' with variables

```
&status
idir = "/ptmp/run"   ! read files from idir
odir = "/ptmp/run"   ! write files to odir
tdir = "/ptmp/run"   ! load FC-Gram tables from tdir
stat = 0             ! status
mult = 1             ! timestep multiplier
bench = 0            ! benchmark
outs = 0             ! output
/
```

The variables in this list are
- `idir`  : directory for binary input
- `odir`  : directory for binary output
- `tdir`  : directory for binary output
- `status`: the number of the last binary file available if a
        simulation is restarted. If zero, a new simulation is done.
- `mult`  : time step multiplier. Useful for restarting a run with
    a different time step but keeping the same cadence when saving
    to files. mult = 2 results in an actual timestep of half what 
    is defined in the dt parameter.
- `benchmark`: if 1, performs a benchmark run.
- `output`: controls the amount of output written in the binary files.

# 
All runs also read a list with parameters for the integration,
named `parameter`

```
&parameter
dt = 1.0e-3      ! time step
step = 10001     ! number of steps
tstep = 500      ! binary output
cstep = 10       ! global quantities
seed = 1000      ! random seed
/
```

The variables in this list are
- `dt`  : sets the time step, which must satisfy the CFL condition
      `dt <= d/u`, where `u` is a characteristic speed, and `d` is
      the smallest of `dx,dy,dz` (the grid spacings).
- `step`: controls the total number of time steps the code will do.
- `tstep` and `cstep`: number of time steps done between exporting
      global quantities (`cstep`), and binary files (`tstep`).
- `seed`: seed for the random number generator.

# 
If solvers with a velocity field are used, `SPECTER` reads initial
parameters for the velocity field from the list 'velocity':

```
&velocity
f0 = 0.37714265   ! amplitude of mechanical forcing
u0 = 1.00         ! amplitude of initial velocity field
kdn = 2.00        ! minimum wavenumber
kup = 4.00        ! maximum wavenumber
nu = 1.e-3        ! viscosity
vibot = 0.0       ! value of vi at z=0 (i can be x or y)
vitop = 0.0       ! value of vi at z=Lz (i can be x or y)
fparam0 = 0.90    ! free parameter for the forcing
vparam0 = 0.90    ! free parameter for the velocity
/
```

These variables are:
- `f0`  : amplitude of the external force.
- `u0`  : amplitude of the initial velocity field.
- `kdn` and `kup`: lower and upper bounds in Fourier space for the
      external force and/or the initial velocity field.
- `nu`  : kinematic viscosity
- `fparam0-9`: variables that can be used to adjust the amplitude of
    individual terms in the expressions of the external force. See the
    files `initialf*` in the directory `examples` for more details.
- `vparam0-9`: idem for the initial velocity field.

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

A list of the supported boundary conditions for a given field is given in
`SOLVERS_AND_BOUNDARY_CONDITIONS.info`

Similar lists are used for other solvers See the file
`examples/parameter.inp` and `specter.fpp` for more details.

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

A description of the contents of all output files for a specific
solver is automatically generated at compile time, named
`README_output.infor`, and placed in the same directory as the
executable file (by default, `../bin`). More details about the output
files can be found in the comments in the source file (the text in
`README_output.info` indicates the specific source file in which the
comments for a given output file can be found).

Examples of simple scripts to read and plot output from `SPECTER`
using Python can be found in the directory `../contrib`. These
examples also include some explanantions on the contents in the
output files.


4 - FC-Gram continuations
-------------------------
FC-Gram continuations employ a series of precomputed tables to compute
a periodic extension of any given data and remove Gibbs ringing in the
Fourier transform. These tables (and, hence, the continuation itself)
depend on two parameters, the length of the extension in number of
points (`C`) and the order of approximation at the boudary (`O`). The
values of `C` and `O` for each direction are defined at compilation
time as described above in section 2. At runtime, the appropriate
tables must exist in the tdir directory, named `AC-O.dat` and `QO.dat`,
with `C` and `O` having the meaning stated above. A set of frequently
used tables can be found in the tables directory, together with a
`README.info` file for more information. Note that for boundary
conditions which involve the i-th derivative the tables `QinO.dat` are
required.


5 - References
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
