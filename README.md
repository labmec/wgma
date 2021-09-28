# wgma

wgma is a library for modal analysis of waveguides using the 
[NeoPZ](https://github.com/labmec/neopz) C++ library.

Along with the library, a few examples available in this repository illustrate the
functionalities and capabilities of the library. All examples are based on 
[a FEM formulation](http://labmec.github.io/neopz/material/availablemats.html#modal-analysis-of-waveguides) using a HCurl-conforming (H1-conforming) approximation space for the
transverse (longitudinal) field components.

The wgma library is currently organised in the following namespaces
- `wgma` : main namespace. Contains the `WGAnalysis` class for managing the modal analysis of waveguides
- `gmeshtools`: auxiliary routines for creating structured curved meshes
- `cmeshtools`: auxiliary routines for dealing with computational meshes containing multiple approximation spaces
- `bctype`: `enum` class with commonly used electromagnetic boundary conditions
- `pmltype`: `enum` class for identifying direction of attenuation of a given PML region
- `slepc`: a few handlers for solving the eigensystem using the EPS module of SLEPc (usage is optional)

## requirements
- A C++ 17 compiler
- [CMake](https://cmake.org/download/) 3.14.0+
- A NeoPZ install configured with MKL support
### optional
- [SLEPc](https://slepc.upv.es) (tested with 3.15.2)

If the SLEPc solver is used, NeoPZ need not have been configured with MKL.
In order to use this package with SLEPc, one must configure both PETSc and SLEPc 
with complex scalar types. 

*Note*: not all SLEPc configurations are available. PRs are welcome.

The following is merely a suggestion on how to configure these libraries:


#### PETSc
```sh
./configure --with-scalar-type=complex --with-debugging=0 
  --with-blaslapack-dir=$MKL_ROOT_DIR 
  --with-mpi-dir=$MPI_ROOT_DIR 
  --download-scalapack=yes --download-mumps=yes
```

#### SLEPc
```sh
./configure  --with-arpack=1 --with-arpack-dir=$ARPACK_BUILD_DIR
  --with-arpack-lib="-L$MKL_LIB_DIR -larpack -lparpack 
  -lmkl_gf_lp64 -lmkl_gnu_thread" --with-feast=0
```

## available examples

### wr90

Analysis of the WR90 waveguide. 

Illustrates
- how to use the wgma library
- how to generate simple meshes in NeoPZ
- how to use the `TPZKrylovSolver`

### stepfiber

Analysis of a step-index optical fiber. 

Illustrates 
- usage of PML
- usage of the `wgma::slepc::EPSHandler` for using SLEPc solvers
- creation of more complex geometries in NeoPZ
