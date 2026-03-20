# Dense_and_sparse_matrix_operations

Overview
--------
This project contains a C++ implementation of dense and sparse linear algebra
operations developed for a finite-volume-method assignment. It includes custom
classes for vectors, dense matrices, sparse matrices, dense linear systems, and
sparse linear systems, together with diagnostics for measuring computational
performance.

The project benchmarks:
- vector addition, subtraction, scalar multiplication, and norms
- dense-matrix addition, subtraction, and scalar multiplication
- sparse-matrix addition, subtraction, and scalar multiplication
- dense and sparse matrix-vector multiplication
- dense and sparse residual evaluation
- nominal and optimised multiplication strategies

The code also supports:
- reading and writing Matrix Market (.mtx) files
- generation of CSV outputs for post-processing and visualisation
- automated diagnostics across multiple matrix sizes

How to run
----------
Run the project from the root directory using:

    ./run.sh

This is the intended way to build and execute the code.

What run.sh does
----------------
The script:
1. removes any existing build directory,
2. configures the project with CMake,
3. builds the executable,
4. enters the build directory,
5. runs the executable `dsm_app`.

Requirements
------------
The project requires:
- a C++17-compatible compiler
- CMake 3.16 or later
- BLAS / OpenBLAS
- a Unix-like shell environment

Project structure
-----------------
Main folders:
- LinearVector/
  Vector class and vector operations

- DenseOperations/
  Dense matrix and dense linear-system classes

- SparseOperations/
  Sparse address, sparse matrix, and sparse linear-system classes

- Mesh/
  Face-addressed mesh and mesh-reading utilities used to generate sparse
  connectivity patterns from finite-volume meshes

- Files/
  Input data used by the diagnostics
  - MTXFiles/ contains Matrix Market test matrices organised by size
  - MeshFiles/ contains mesh-based input data

Main files:
- main.cpp
  Runs the full set of diagnostics
- Diagnostics.hpp
  Contains the timing and benchmark routines
- PathConfig.hpp
  Handles project input/output paths
- run.sh
  Build-and-run script
- CMakeLists.txt
  Build configuration

Output
------
During execution, the program creates an Output/ directory in the build folder.
The diagnostics write CSV files containing timing data and sample outputs for
post-processing.

Notes
-----
- The recommended execution method is:

      ./run.sh

- The project links against BLAS for the optimised dense matrix-vector
  multiplication benchmark.
- The mesh-based sparse operations use owner-neighbour connectivity derived
  from face-addressed finite-volume meshes.
