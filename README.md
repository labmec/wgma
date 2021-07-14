# wgma

wgma is a library for modal analysis of waveguides using the 
[NeoPZ](https://github.com/labmec/neopz) C++ library.

Along with the library, a few examples available in this repository illustrate the
functionalities and capabilities of the library. All examples are based on 
[a FEM formulation](http://labmec.github.io/neopz/material/availablemats.html#modal-analysis-of-waveguides) using a HCurl-conforming (H1-conforming) approximation space for the
transverse (longitudinal) field components.

The wgma library currently consists of
- `gmeshtools`: auxiliary routines for creating structured curved meshes
- `cmeshtools`: auxiliary routines for dealing with computational meshes containing multiple approximation spaces
- `bctype`: `enum` class with commonly used electromagnetic boundary conditions

## requirements
- A C++ 17 compiler
- [CMake](https://cmake.org/download/) 3.13.0+
- A NeoPZ install configured with LAPACK support
