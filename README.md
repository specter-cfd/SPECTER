![](./specter.svg)

# SPECTER

----
**SPECTER** (**S**pecial **Pe**riodic **C**ontinuation **T**urbulence Solv**er**) is a **Fortran** program (with some **C** bindings) that performs numerical integration in 3D space for a set of PDEs usually found in fluid dynamics. 

A pseudo-spectral method is used to compute spatial derivatives, coupled with an adjustable order Runge-Kutta time stepping. When non-periodic directions are specified, **SPECTER** utilizes FC-Gram continuations to compute Fourier transforms with minimal overhead and considerable accuracy. At its present state, **SPECTER** has support for a single non-periodic direction. Support for two or more non-periodic directions is planned.

**SPECTER** is programmed to scale. It can be run in your workstation as well as in top notch supercomputers. With this in mind, the code is parallelized employing **MPI**, **OpenMP** and **CUDA**.

Note: **CUDA** support is currently experimental and not provided publicly yet. If you are interested in CUDA support please contact the author.

----
## Getting SPECTER
If you have git installed in your system, you can run 
```
git clone https://github.com/mfontanaar/SPECTER
```
Alternatively, you can fetch it using `curl` or `wget` with
```
curl/wget https://github.com/mfontanaar/SPECTER
```
If you are using a web browser, you can also get it via [this link](https://github.com/mfontanaar/SPECTER/archive/master.zip).

## Prequisites
**SPECTER** needs **Fortran** and **C** compilers to be present in the system, and an installation of both **FFTW** 3.x and an **MPI** flavor (`mpich`, `openmpi`, etc.) is required. Additionally, for compiling with **CUDA** support the [nVidia Cuda Compiler](https://developer.nvidia.com/cuda-llvm-compiler) is needed.

## Compiling
Compiling **SPECTER** is as simple as setting appropriate paths to the compilers, and the `fftw` and `mpi` libraries in the file `src/Makefile.in`. Parameters like simulation resolution, solver to employ, time-stepping order and floating-point precision can also be specified in the same file. The code is then compiled by running 
```
make
```
while in the `src/` directory. After running `make`, the output binary should be placed in the `bin/` directory, and named like the chosen solver.

## Running SPECTER
After the compilation, a binary file with the name of the selected solver is placed in `bin/` directory. In that same directory, a file named `parameter.inp` can be found, where you can specify several important parameters like magnitude of the time step, the number of steps to perform, how often to save output information. Physical parameters like the forcing amplitude, the mmagnitude of the initial condition or the viscosity are also specified in `bin/parameter.inp/`. A full list of the parameters with their description can be found in the same file, as well as in `src/README.txt`. After setting up `parameters.inp` you can run the code with
```
./SOLVER
```
If you are running **SPECTER** in a cluster that uses a queuing system like `TORQUE` or `Slurm` you must follow the instructions provided by the cluster administrator to add `SOLVER` to the queue.

## Forcing and initial conditions
The expressions for the forcing fields and the initial conditions must be specified at compile time. They are defined in the files `initialv.f90` or `initialfv.f90` for the initial velocity field and the initial mechanical forcing, respectively. The same convention is used for the scalar or other fields. A set of commonly used forcings and initial conditions can be found in the `src/examples/` directory.

## Getting FC-Gram tables
For non-periodic directions, **SPECTER** uses tables to compute appropriate periodic extensions. A collection of commonly employed tables can be found in the `tables/` directory. However, a supplementary **Fortran** program to generate these tables is available at [https://github.com/mfontanaar/fctables](https://github.com/mfontanaar/fctables).

## More information
To get the best out of **SPECTER** we recommend reading the information provided in `src/README.txt` for a detailed explanation of each compilation and runtime parameter.

## References
- Fontana, M., Bruno, O. P., Mininni P. D. & Dmitruk P.; *Fourier continuation method for incompressible fluids with boundaries*. arXiv: [2002.01392](https://arxiv.org/abs/2002.01392).
- Bruno, O. P. & Lyon, M: *High-order unconditionally stable FC-AD solvers for general smooth domains I. Basic elements*. JCP.  229, (2010). DOI: [10.1016/j.jcp.2009.11.020](https://doi.org/10.1016/j.jcp.2009.11.020).
- Mininni P. D., Rosenberg D. L., Reddy R. & Pouquet A.; *A hybrid MPI–OpenMP scheme for scalable parallel pseudospectral computations for fluid turbulence*. P. Comp. 37, 123 (2011). DOI: [10.1016/j.parco.2011.05.004](https://doi.org/10.1016/j.parco.2011.05.004).
- Rosenberg D. L., Mininni P. D., Reddy R. & Pouquet A.; ArXiV (2018). *GPU parallelization of a hybrid pseudospectral fluid turbulence framework using CUDA*. arXiv: [arXiv:1808.01309](https://arxiv.org/abs/1808.01309).

## Citation
If you use **SPECTER** for a publication, we kindly ask you to cite the following article **Fontana, M., Bruno, O. P., Mininni P. D. & Dmitruk P.; *Fourier continuation method for incompressible fluids with boundaries*. arXiv: [arXiv:2002.01392](https://arxiv.org/abs/2002.01392)**.

## Authors
- Mauro Fontana - Physics Department at University of Buenos Aires.
- Pablo D. Mininni - Physics Department at University of Buenos Aires.
