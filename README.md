# wgma

wgma is a library for analysis of optical waveguides using the 
[NeoPZ](https://github.com/labmec/neopz) C++ library.

**notice:** currently, wgma must be linked against the `develop` version of NeoPZ.

Along with the library, a few examples available in this repository illustrate the
functionalities and capabilities of the library. 

For the modal analysis of waveguides, all examples are based on 
[a FEM formulation](http://labmec.github.io/neopz/material/availablemats.html#modal-analysis-of-waveguides) 
using a HCurl-conforming (H1-conforming) approximation space for the
transverse (longitudinal) field components.

For the scattering analysis of planar waveguides, the implementation follows from
[Tsuji,Koshiba, 2002](http://opg.optica.org/jlt/abstract.cfm?URI=jlt-20-3-463).

The wgma library is currently organised in the following namespaces

- `wganalysis` : namespace for modal analysis of waveguides. 
  Contains the `Wgma2D` class for managing the modal analysis
of waveguides with 2D cross-section.
- `gmeshtools`: auxiliary routines for creating structured curved meshes, reading `.msh`
files and directional *h*-refinement
- `cmeshtools`: auxiliary routines for dealing with computational meshes and PML
- `bctype`: `enum` class with commonly used electromagnetic boundary conditions
- `pmltype`: `enum` class for identifying direction of attenuation of a given PML region
- `modetype`: `enum` class for distinguishing between TE/TM modes in the analysis of
planar waveguides
- `slepc`: a few handlers for solving the eigensystem using the EPS module of SLEPc (usage is optional)

## requirements
- A C++ 17 compiler
- [CMake](https://cmake.org/download/) 3.14.0+
- A NeoPZ install configured with MKL support (**note**: `develop` version)
### optional
- [SLEPc](https://slepc.upv.es) (tested with 3.15.2)
- [gmsh](https://gmsh.info) (formats msh3 and msh4 are supported. tested only with msh4)

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

### ribwg

Analysis of a rib waveguide.

Illustrates
- reading a `.msh` file from gmsh
- directional refinement (useful for singularities)


### ecf

Analysis of an Exposed-Core-Fiber

Illustrates
- performance in a real-world scenario
- how to generate dispersion curves and export to `.csv`

### planar_wg

Scattering analysis of a planar slab waveguide with a discontinuity

Illustrates
- how to analyse 2D planar waveguides for a given excitation source
- how to analyse 2D planar waveguide discontinuities using port truncation with PMLs
- how to prescribe custom sources

### pcwg

Modal analysis of a 2D photonic crystal followed by scattering analysis

Illustrates
- how to easily generate `.msh` meshes in python allowing for exact representation 
of curved geometries
- how PMLs can be set up for periodic domains
- how to perform more complex analysis and transfer the solution between different 
problems