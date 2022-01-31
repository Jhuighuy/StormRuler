<!--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-->
# StormRuler🦜 — A very high order CFD solver
<!--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-->
**StormRuler** is a very high order multidimensional CFD solver, 
written in C++20 and Fortran 2018.

<!----------------------------------------------------------------->
## 🏗 Compiling
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
## 🌀 Equations solved
<!----------------------------------------------------------------->
**StormRuler** features support of the various set of the
partial differential equations, including:
* 💧 Cahn-Hilliard equation,
* 🌊 Incompressible Navier-Stokes equations,
* 🌪 _(planned)_ Сompressible Navier-Stokes/Euler equations,
* ...

<!----------------------------------------------------------------->
## 🌐 Numerical methods
<!----------------------------------------------------------------->
The heart of the **StormRuler** is the ✨Finite Difference Method✨.

<!----------------------------------------------------------------->
## 🌈 Algebraic solvers
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
- 🛸 Nonlinear solvers:
  * Newton solver 
    (`Newton`, for the general nonlinear problems),
  * Jacobian-Free Newton-Krylov solver 
    (`JFNK`, for the general nonlinear problems);

- 🏎 Linear iterative solvers:
  * _(planned)_ Richardson iteration solver
    (`Richardson`, for the _general_ linear problems);
  * Conjugate Gradients solver 
    (`CG`, for the _definite symmetric_ linear problems),
  * _(planned)_ Flexible Conjugate Gradients solver 
    (`FCG`, for the _definite symmetric_ linear problems
     with _flexible preconditioning_),
  * Minimal Residual solver
    (`MINRES`, for the _indefinite symmetric_ linear problems),
  * Conjugate Gradients Squared solver 
    (`CGS`, for the _general_ linear problems),
  * Biconjugate Gradients (stabilized) solver
    (`BiCGStab`/`BiCGStab(l)`, for the _general_ linear problems),
  * Transpose-Free Quasi-Minimal Residual solver
    (`TFQMR`/`TFQMR(1)`, for the _general_ linear problems),
  * Induced Dimension Reduction solver
    (`IDR(s)`, for the _general_ linear problems);
  * Generalized Minimal Residual solver
    (`GMRES`, for the _general_ linear problems),
  * Flexible Generalized Minimal Residual solver
    (`FGMRES`, for the _general_ linear problems
     with _flexible preconditioning_),
  * _(planned)_ Loose (Flexible) Generalized Minimal Residual solver,
    (`LGMRES`/`LFGMRES`, for the _general_ linear problems),

- 🚜 Linear iterative least squares solvers:
  * Least squares-QR solver
    (`LSQR`, for the linear _least squares_ problems),
  * Least squares-MINRES solver
    (`LSMR`, for the linear _least squares_ problems);

- 🚂 Linear direct solvers:
  * _(planned)_ PARDISO direct solver
    (`PARDISO`/`PARDISO_MKL`),
  * _(planned)_ SuperLU direct solver
    (`SuperLU`);

- 🚨 Preconditioners:
  * _(planned)_ Diagonal preconditioner
    (`Jacobi`, for any problems),
  * Symmetric Gauss-Seidel preconditioner
    (`LU_SGS`, for any linear problems),
  * _(planned)_ Incomplete Cholesky preconditioner
    (`IC0`, `ICT`, for _definite symmetric_ problems),
  * Incomplete LU preconditioner
    (`ILU0`, `ILUT`, for _unsymmetric_ problems, _currently requires MKL_),
  * _(planned)_ Approximate Inverse preconditioner
    (`AINV0`, `AINV`, for _symmetric_ problems),
  * _(planned)_ SPAI preconditioner
    (`SPAI0`, `SPAI`, for _unsymmetric_ problems);
  * Chebyshev Polynomial preconditioner
    (`Cheby`, for _definite symmetric_ problems),
  * _(planned)_ Polynomial preconditioner
    (`Poly`, ...),
  * _(planned)_ Krylov preconditioner
    (`Krylov`, ...),

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
  - [ ] 🚬🐏 Python API.

* Image:
  - [ ] 🦜💄 Image API.

* Mesh ordering:
  - [ ] 🦜🧻 Some C/C++ API for mesh loading,
  - [x] 🪓 Cache-friendly cell sorting: Hilbert Sort,
  - [x] 🪓 Cache-friendly cell sorting: METIS,
  - [ ] 🐞 Something looks broken..
  - [ ] 🚬 Better cell ordering quality functional, 
  - [ ] 🧸 Functional-based unified API for cell ordering,
  - [ ] 🪓 BC cells sorting and better BCs parallelization.

* IO:
  - [ ] 🧻 Refactor IO lists into some more intersting API,
  - [ ] 🪓 Move Neato output away from mesh,
  - [x] 🪓 Move VTK output away from mesh,
  - [x] 🚬 Redesigned VTK output (as `.vti`),
  - [x] 🦜 ZLib compression,
  - [ ] 🧻 Refactor the compression headers,
  - [ ] 🪓 Parallel VTK output (as `.pvti`),
  - [ ] 🚬 Multiblock VTK output (as `.vtm`),
  - [ ] 🦜💄 Better image IO.

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
  - [ ] 💄 Some better residual monitor,
  - [x] 🧸 Report true residual in `CG`,
  - [ ] 🪓 Eigenvalue estimates in `CG`,
  - [ ] 🪓 `FCG` solver implementation,
  - [ ] 🧻 Fix unpreconditioned `MINRES` solver,
  - [ ] 🧻 `MINRES` uses too many vectors (10), can be implemented with 7.
  - [ ] 🦜 `CGS` solver penaltization,
  - [ ] 🦜 `CGS` breakdown detection,
  - [ ] 🦜 `BiCGStab` breakdown detection,
  - [x] 🚬 `BiCGStabL` solver implementation,
  - [x] 🪓 `TFQMR` solver implementation,
  - [x] 🪓 Preconditioned `TFQMR` solver,
  - [ ] 🪓 More accurate `TFQMR` residual estimator,
  - [x] 🪓 `TFQMR1` solver implementation,
  - [ ] 🦜 `TFQMR`/`TFQMR1` breakdown detection,
  - [x] 🚬 `GMRES` solver implementation,
  - [x] 🪓 Preconditioned `GMRES` implementation,
  - [x] 🪓 Right preconditioned `FGMRES` implementation,
  - [ ] 🦜 `GMRES`/`FGMRES` breakdown detection,
  - [ ] 🪓 Eigenvalue estimates in `GMRES`,
  - [ ] 🧻 Refactor `LSQR` and `LSMR` to support the non-square operators,

* Matrix operations and extraction:
  - [x] 🧸 CSR matrix class, CSR matvec,
  - [x] 🧸 CSR Extraction with prescribed coloring,
  - [ ] 🦜 Block extraction with prescribed coloring,
  - [x] 🪓 Bandwidth-based column coloring problem,
  - [x] 🪓 Portrait-based column coloring problem,
  - [ ] 🦜 Diagonal part extraction,
  - [ ] 🦜 Triangular part extraction,
  - [ ] 🚬 Some more optimal column coloring algorthms..

* Preconditioning:
  - [x] 🪓 `LU_SGS` preconditioner,
  - [x] 🪓 MKL-based `ILU0` preconditioners,
  - [x] 🪓 MKL-based `ILUT` preconditioners,
  - [ ] 🚬 Custom `ILU0`/`IC0` preconditioners,
  - [ ] 🚬 Custom `ILUT`/`ICT` preconditioners,
  - [x] 🧸 Chebyshev polynomial preconditioner,
  - [ ] 🦜 Some other polynomial preconditioners,
<!--
  - [ ] 🦜 `Jacobi` preconditioner,
  - [ ] 🚬 `SPAI0` preconditioner,
  - [ ] 🚬🚬 `SPAI` preconditioner,
  - [ ] 🚬 'Left' `SPAI` preconditioner,
  - [ ] 🦜 `AINV0` preconditioner,
  - [ ] 🦜 `AINV` preconditioner,
  - [ ] 🦜 Krylov preconditioner.
-->

- Direct solvers:
  - [x] 🚬 Optimized partial matrix operations with MKL-comparable performance,
  - [x] 🧸 Sequential triangular solvers,
  - [x] 🚬 Parallel DAG-based triangular solvers,
  - [ ] 🚬🚬 Parallel block diagonal extraction-based triangular solvers,
<!--
  - [ ] 🪓 Dense direct solver,
  - [ ] 🪓 Sparse-approximate direct solver,
  - [ ] 🦜 Built-in direct solver,
  - [ ] 🦜 Direct solvers (`PARDISO`, `SuperLU`).
-->

* Nonlinear solvers:
  - [x] 🧸 Newton-Raphson solver,
  - [ ] 💄 Better API for the exact Newton-Raphson solver, 
  - [ ] 🦜 Relaxed Newton solver,
  - [x] 🧸 Jacobian-Free Newton-Krylov solver,
  - [x] 🧻 Optimized first order JFNK,
  - [x] 🧸 Select an epsilon in the first order JFNK,
  - [ ] 🦜 Nonlinear preconditioning..
