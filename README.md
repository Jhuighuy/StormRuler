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

<!--For the sake of convenience, all auxiliary solvers are implemented 
in the matrix-free manner: no assembled matrix is required to find 
a solution of the algebraic problem, only the matrix-vector product 
function is used.

Although most of the problems can be solved in the matrix-free 
manner using the Krylov subspace iterative solver, in some 
pathological cases an assembled matrix be required to 
construct a suitable preconditioner or utilize a direct solver.
**StormRuler** reconstructs a matrix using the matrix-vector 
product function automatically, using the 
_graph coloring based-algorithm_ in order to minimize an 
amount of the matrix-vector products required to construct it.-->

Several preconditioners are available in pure matrix-free context,
although some of them are not parallelized.

**StormRuler** contains:
- 🛸 Matrix-free nonlinear solvers:
  * Newton-Raphson solver 
    (`Newton`, for the general nonlinear problems),
  * Jacobian-Free Newton-Krylov solver 
    (`JFNK`, for the general nonlinear problems);

- 🏎 Matrix-free linear iterative solvers:
  * Conjugate Gradients solver 
    (`CG`, for the _definite symmetric_ linear problems),
  * Biconjugate Gradients (stabilized) solver
    (`BiCGSTAB`, for the general _non-singular_ linear problems),
  * Minimal Residual solver
    (`MINRES`, for the indefinite _symmetric_ linear problems),
  * Generalized Minimal Residual method solver
    (`GMRES`, for the general linear problems),
  * _(planned)_ Quasi-Minimal Residual solver
    (`QMR`, for the general linear problems);
  * _(planned)_ Transpose-free Quasi-Minimal Residual solver
    (`TFQMR`, for the general linear problems);

<!--
- 🚂 Linear direct solvers (embedded into the matrix-free environment):
  * MKL Direct Sparse Solver 
    (`DSS_MKL`, from [MKL DSS](https://intel.ly/37N95pe)).-->

- 🚜 Matrix-free linear iterative least squares solvers:
  * Least squares-QR solver
    (`LSQR`, for the general linear least squares problems),
  * Least squares-MINRES solver
    (`LSMR`, for the general linear least squares problems);

- 🚨 Matrix-free preconditioners:
  * Block Jacobi preconditioner
    (`Jacobi`),
  * Block LU-SGS preconditioner
    (`LU_SGS`),
  * _(planned)_ Block SPAI preconditioner.

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

* Parameters system:
  - [ ] 💄 ???

* Mesh:
  - [ ] 🧸 Move kernel runners away from mesh,
  - [ ] 🧸 BC kernels,
  - [ ] 🪓 Refactor mesh generation in rectangle/cube,
  - [ ] 🪓 Refactor mesh generation with image,
  - [ ] 🪓 Generate nodes,
  - [ ] 🪓 Generate faces,
  - [ ] 🧸 Redesigned VTK output,
  - [ ] 🚬 Better cell ordering quality functional, 
  - [x] 🪓 Cache-friendly cell sorting: Hilbert Sort,
  - [x] 🪓 Cache-friendly cell sorting: METIS,
  - [ ] 🧸 Unified API for cell sorting,
  - [ ] 🪓 BC cells sorting and better BCs parallelization,
  - [ ] 🚬 Mesh coarsening and refinement (pre GMG),
  - [ ] 🪓 Block mesh (pre MPI),
  - [ ] 🚬🚬 Non-conforming multilevel mesh,
  - [ ] 🚬🚬🚬 AMR,
  - [ ] 🚬🚬🚬 Cut cell methods.

* New differential operators and boundary conditions:
  - [ ] 🧸 Variable weight Laplace operator with 4+ order approx.,
  - [ ] 🦜 Tensor weight Laplace operator with 4+ order approx.,
  - [ ] 🦜 High order convection approx.,
  - [x] 🪓 Cylindrical symmetry for 2D domains,
  - [ ] 💄 More special boundary conditions,
  - [ ] 🪓 Godunov/WENO linear convection operator,
  - [ ] 🪓 Godunov/WENO nonlinear convection operator,
  - [ ] 🪓 Riemann solvers, Euler equations...

* Linear solvers:
  - [ ] 🧻 Clean-up unified solver to use conjugate MatVec,
  - [ ] 🧻 Convergence parameters in C/C++ API,
  - [ ] 💄 Some better residual monitor,
  - [x] 🚬 `GMRES` solver implementation,
  - [ ] 🪓 Preconditioned `GMRES` implementation (right preconditioned?),
  - [ ] 🪓 `QMR` solver implementation,
  - [ ] 🪓 `TFQMR` solver implementation.

* Matrix reconstruction:
  - [ ] 🪓 Matrix portrait construction,
  - [ ] 🚬🚬🚬 Graph-coloring problem,
  - [ ] 🪓 `ILU`/`ICHOL` preconditioners,
  - [ ] 🪓 `SPAI` preconditioner,
  - [ ] 🦜 Direct solver (`MKL_DSS`, `PARDISO`, `SuperLU`).

* Preconditioning:
  - [x] 🧻 Refactor precondtioner from function pointer to class,
  - [ ] 🧸 Add user-defined preconditioner in C/C++ API,
  - [ ] 🦜 Fix block case diagonal extraction,
  - [ ] 🦜 Block tridiagonal preconditioner,
  - [ ] 🚬🚬🚬 Matrix-free SPAI preconditioner,
  - [ ] 🚬 V-cycle GMG,
  - [ ] 🚬🚬 F-cycle GMG,
  - [ ] 🚬🚬🚬 W-cycle GMG.

* Nonlinear solvers:
  - [x] 🧸 Newton-Raphson solver,
  - [ ] 💄 Better API for the exact Newton-Raphson solver, 
  - [ ] 🦜 Relaxed Newton solver,
  - [x] 🧸 Jacobian-Free Newton-Raphson solver,
  - [ ] 🧻 Optimized first order JFNK,
  - [ ] 🧸 Select an epsilon in the first order JFNK,
  - [ ] 🧸 Second order JFNK solver,
  - [ ] 🦜 Nonlinear preconditioning..
