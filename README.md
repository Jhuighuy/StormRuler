<!--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-->
# StormRuler🦜 — A very high order CFD solver
<!--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-->
**StormRuler** is a very high order multidimensional CFD solver, 
written in Fortran 2018 and C11.

<!----------------------------------------------------------------->
## 🏗Compiling
<!----------------------------------------------------------------->

Supported compilers:
* _GCC_ version _9.0_ and more recent 
  (tested on _10.3.0_).
* _Intel Classic compilers_ version _19.0_ and more recent
  (tested on _2021.3.0_).
<!--* _AMD AOCC_ version _3.1.0_ and more recent
  (tested on _3.1.0_).
* _PGI Compilers_ (from _NVIDIA HPC SDK_) version 21 and more recent 
  (tested on _21.07_).
* _NAG Fortran Compiler_ version 7.0 and more recent
  (tested on _7.0 build 7048_).-->

macOS with Intel compilers:
```zsh
export LIBRARY_PATH="$LIBRARY_PATH:/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib"
```

<!----------------------------------------------------------------->
## 🌀Equations solved
<!----------------------------------------------------------------->
**StormRuler** features support of the various set of the
partial differential equations, including:
* 💧 Cahn-Hilliard equation,
* 🌊 Incompressible Navier-Stokes equations,
* 🌪 _(planned)_ Сompressible Navier-Stokes/Euler equations,
* ...

<!----------------------------------------------------------------->
## 🌐Numerical methods
<!----------------------------------------------------------------->
The heart of the **StormRuler** is the ✨Finite Difference Method✨.

<!----------------------------------------------------------------->
## 🌈Algebraic solvers
<!----------------------------------------------------------------->
In order to implement the high performance implicit schemes,
several algebraic problems, like systems of linear and nonlinear
equations, have to be solved.

For the sake of convenience, all auxiliary solvers are implemented 
in the _matrix-free_ manner: no assembled matrix is required to find 
a solution of the algebraic problem, only the matrix-vector product 
function is used.

Although most of the problems can be solved in the matrix-free 
manner using the Krylov subspace iterative solver, in some 
cases an assembled matrix be required to construct a suitable 
preconditioner or utilize a direct solver.
**StormRuler** extracts a matrix using the matrix-vector 
product function automatically, using the _column coloring algorithm_ 
in order to minimize an amount of the extraction matrix-vector 
products.

**StormRuler** contains:
- 🛸 Matrix-free nonlinear solvers:
  * Newton-Raphson solver 
    (`Newton`, for the general nonlinear problems),
  * Jacobian-Free Newton-Krylov solver 
    (`JFNK`, for the general nonlinear problems);

- 🏎 Matrix-free linear iterative solvers:
  * Conjugate Gradients solver 
    (`CG`, for the _definite symmetric_ linear problems),
  * _(planned)_ Flexible Conjugate Gradients solver 
    (`FCG`, for the _definite symmetric_ linear problems with),
  * Biconjugate Gradients (stabilized) solver
    (`BiCGStab`, for the general _non-singular_ linear problems),
  * Minimal Residual solver
    (`MINRES`, for the indefinite _symmetric_ linear problems),
  * Generalized Minimal Residual method solver
    (`GMRES`, for the general linear problems),
  * _(planned)_ Flexible Generalized Minimal Residual method solver
    (`FGMRES`, for the general linear problems),
  * _(planned)_ Transpose-free Quasi-Minimal Residual solver
    (`TFQMR`, for the general linear problems);

- 🚜 Matrix-free linear iterative least squares solvers:
  * Least squares-QR solver
    (`LSQR`, for the general linear least squares problems),
  * Least squares-MINRES solver
    (`LSMR`, for the general linear least squares problems);

- 🚂 Linear direct solvers (embedded into the matrix-free environment):
  * _(planned)_ MKL Direct Sparse Solver 
    (`DSS_MKL`, from [MKL DSS](https://intel.ly/37N95pe)).
  * _(planned)_ PARDISO direct solver
    (`PARDISO`/`PARDISO_MKL`),
  * _(planned)_ SuperLU direct solver
    (`SuperLU`);

- 🚨 Matrix-based preconditioners (embedded into the matrix-free environment):
  * _(planned)_ Jacobi preconditioner
    (`Jacobi`),
  * LU Symmetric Gauss-Seidel (LU-SGS) preconditioner
    (`LU-SGS` for any linear problems),
  * _(planned)_ Incomplete Cholesky preconditioner,
    (`ICHOL`)
  * _(planned)_ Incomplete LU preconditioner,
    (`ILU0`, `ILUT`),
  * _(planned)_ Approximate Inverse preconditioner
    (`AINV`),
  * _(planned)_ SPAI preconditioner,
    (`SPAI`);

- Matrix-free preconditioners:
  * _(planned)_ Polynomial preconditioner
    (`Poly`),
  * _(planned)_ Krylov preconditioner
    (`Krylov`),

<!----------------------------------------------------------------->
## 🛤Road map
<!----------------------------------------------------------------->

Legend:
- 🧸 — _easy problem_,
- 🪓 — _medium complexity problem_,
- 🚬 — _hardcore feature/problem_,
- 🦜 — _unidentified complexity_,
- 💄 — _requires creativity_,
- 🧻 — _refactoring required_,
- 🐏 — _system programming skills required_,
- 🐞 — _our bug, fix required_,
- 🪦 — _compiler bug, workaround required_.

* C/C++ API:
  - [x] 🐏 Pure C API,
  - [ ] 🐏 Higher-level C++ API,
  - [ ] 🚬🐏 Python API.

* General architecture:
  - [ ] 🦜🐞 Segfaults,
  - [x] 🐞 with `type(c_ptr), value :: env` fixed globally for ifort,
  - [ ] 🚬 GPU support,
  - [ ] 🚬🚬 MPI support.

* Options system:
  - [ ] 💄 ???

* Mesh:
  - [ ] 🧸 Move kernel runners away from mesh,
  - [ ] 🧸 BC kernels,
  - [x] 🪓 Reimplement mesh generation with support for the varous DnQm models.
  - [ ] 🪓 Generate nodes,
  - [ ] 🪓 Generate faces,
  - [ ] 🧸 Redesigned VTK output,
  - [ ] 🪓 Move VTK output away from mesh.

* Mesh ordering:
  - [ ] 🦜🧻 Some C/C++ API for mesh loading,
  - [x] 🪓 Cache-friendly cell sorting: Hilbert Sort,
  - [x] 🪓 Cache-friendly cell sorting: METIS,
  - [ ] 🐞 Something looks broken..
  - [ ] 🚬 Better cell ordering quality functional, 
  - [ ] 🧸 Functional-based unified API for cell ordering,
  - [ ] 🪓 BC cells sorting and better BCs parallelization.

* AMR/cut cell:
  - [ ] 🪓 Block mesh (pre MPI),
  - [ ] 🚬🚬 Non-conforming multilevel mesh,
  - [ ] 🚬🚬🚬 AMR...
  - [ ] 🚬🚬🚬 Cut cell methods...

* GMG:
  - [ ] 🚬 Mesh coarsening and refinement (pre GMG),
  - [ ] 🚬🚬 V-cycle GMG,
  - [ ] 🚬🚬 F-cycle GMG,
  - [ ] 🚬🚬 W-cycle GMG.

* LBM:
  - [ ] 🪓 Correct streaming operator.
  - [ ] 🧸 SRT collision operator,
  - [ ] 🪓 MRT collision operator,
  - [ ] 🪓 Free boundary LBM boundary conditions,
  - [ ] 🪓 Bounce-back LBM boundary conditions,
  - [ ] 🦜 Color gradients, ...

* New differential operators and boundary conditions:
  - [ ] 🧸 Variable weight Laplace operator with 4+ order approx.,
  - [ ] 🦜 Tensor weight Laplace operator with 4+ order approx.,
  - [ ] 🦜 High order convection approx.,
  - [x] 🪓 Cylindrical symmetry for 2D domains,
  - [ ] 💄 More special boundary conditions,
  - [ ] 🪓 Godunov/WENO linear convection operator,
  - [ ] 🪓 Godunov/WENO nonlinear convection operator,
  - [ ] 🪓 Riemann solvers, Euler equations...

* Linear iterative solvers:
  - [ ] 🧻 Clean-up unified solver to use conjugate MatVec,
  - [ ] 🧻 Convergence parameters in C/C++ API,
  - [ ] 🧻 Non-uniform solver on higher-level,
  - [ ] 💄 Some better residual monitor,
  - [x] 🚬 `GMRES` solver implementation,
  - [ ] 🪓 Preconditioned `GMRES` implementation (right preconditioned?),
  - [ ] 🪓 `TFQMR` solver implementation.
  - [ ] 🪓 `FCG` solver implementation,
  - [ ] 🪓 `FGMRES` solver implementation.

* Matrix extraction:
  - [x] 🧸 CSR matrix class, CSR matvec,
  - [ ] 🧸 CSC matrix class, fast CSR-CSC tranpositions,
  - [x] 🧸 CSR Extraction with prescribed coloring,
  - [ ] 🧸 CSC Extraction with prescribed coloring,
  - [ ] 🧸 Fill matrix diagonal function.
  - [ ] 🧸 Extract matrix diagonal function.
  - [ ] 🧸 Extract matrix row as a sparse vector function.
  - [ ] 🧸 Sparse-sparse approximate AXPY. 
  - [ ] 🧸 Matrix symmetrization.
  - [ ] 🧸 Partial matrix-vector products in DL, DU modes.
  - [ ] 🦜 Block extraction with prescribed coloring,
  - [x] 🪓 Bandwidth-based column coloring problem,
  - [x] 🪓 Portrait-based column coloring problem,
  - [ ] 🚬 Some more optimal column coloring algorthms..

* Preconditioning:
  - [ ] 🦜🧻 Refactor unified and C/C++ API for preconditioning,
  - [ ] 🦜🧻 Add user-defined preconditioner in C/C++ API,
  - [x] 🧻 Refactor precondtioner from function pointer to class,
  - [ ] 🦜 `Jacobi` preconditioner,
  - [x] 🪓 `LU_SGS` preconditioner,
  - [x] 🪓 MKL-based `ILU0_MKL` preconditioners,
  - [ ] 🪓 MKL-based `ILUT_MKL` preconditioners,
  - [ ] 🚬 `ILU`/`ICHOL` preconditioners,
  - [ ] 🚬 Static `SPAI` preconditioner,
  - [ ] 🚬🚬 Dynamic `SPAI` preconditioner,
  - [ ] 🚬 'Left' `SPAI` preconditioner,
  - [ ] 🦜 `AINV` preconditioner,
  - [ ] 🦜 polynomial preconditioner.

- Direct solvers:
  - [x] 🚬 Optimized partial matrix operations with MKL-comparable performance,
  - [ ] 🪓 Dense direct solver,
  - [ ] 🪓 Sparse-approximate direct solver,
  - [x] 🧸 Sequential triangular solvers,
  - [x] 🚬 Parallel DAG-based triangular solvers,
  - [ ] 🚬🚬 Parallel block diagonal extraction-based triangular solvers,
  - [ ] 🦜 Built-in direct solver,
  - [ ] 🦜 Direct solvers (`MKL_DSS`, `PARDISO`, `SuperLU`).

* Nonlinear solvers:
  - [x] 🧸 Newton-Raphson solver,
  - [ ] 💄 Better API for the exact Newton-Raphson solver, 
  - [ ] 🦜 Relaxed Newton solver,
  - [x] 🧸 Jacobian-Free Newton-Krylov solver,
  - [x] 🧻 Optimized first order JFNK,
  - [x] 🧸 Select an epsilon in the first order JFNK,
  - [ ] 🦜 Nonlinear preconditioning..
