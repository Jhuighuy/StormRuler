<!--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-->
# StormRuler🦜 — A very high order CFD solver
<!--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-->
**StormRuler** is a very high order multidimensional CFD solver, 
written in Fortran 2018 and C11.

<!----------------------------------------------------------------->
## Compiling
<!----------------------------------------------------------------->

Supported compilers:
* _GCC_ version _9.0_ and more recent 
  (tested on _10.3.0_).
* _Intel Classic compilers_ version _19.0_ and more recent
  (tested on _2021.3.0_).
* _AMD AOCC_ version _3.1.0_ and more recent
  (tested on _3.1.0_).
* _PGI Compilers_ (from _NVIDIA HPC SDK_) version 21 and more recent 
  (tested on _21.07_).
* _NAG Fortran Compiler_ version 7.0 and more recent
  (tested on _7.0 build 7048_).

macOS with Intel compilers:
```zsh
export LIBRARY_PATH="$LIBRARY_PATH:/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib"
```

<!----------------------------------------------------------------->
## Equations solved
<!----------------------------------------------------------------->
**StormRuler** features support of the various set of the
partial differential equations, including:
* Incompressible Navier-Stokes equations,
* Cahn-Hilliard equation,
* _(planned)_ Streamfunction-Vorticity formulation for the
  Incompressible Navier-Stokes equations in 2D,
* Incompressible Navier-Stokes-Cahn-Hilliard equations,
* _(planned)_ Streamfunction-Vorticity formulation for the
  Incompressible Navier-Stokes-Cahn-Hilliard equations in 2D.

<!----------------------------------------------------------------->
## Numerical methods
<!----------------------------------------------------------------->
The heart of the **StormRuler** is the Finite Difference Method.

<!----------------------------------------------------------------->
## Auxiliary solvers
<!----------------------------------------------------------------->
In order to implement the high performance implicit schemes,
several auxiliary problems, like systems of linear and nonlinear
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
- Matrix-free nonlinear solvers:
  * Newton-Raphson solver 
    (`Newton`, for the general nonlinear problems),
  * Jacobian-Free Newton-Krylov solver 
    (`JFNK`, for the general nonlinear problems);

- Matrix-free linear iterative solvers:
  * Chebyshev Semi-Iterative solver
    (`Chebyshev`, for the _definite symmetric_ linear problems),
  * Conjugate Gradients solver 
    (`CG`, for the _definite symmetric_ linear problems),
  * _(planned)_ Hybrid Chebyshev-CG solver,
    (`ConjugateChebyshev`, for the _definite symmetric_ linear problems,
     _theoretically **the fastest** existing SPD-solver 
      for the repeated-computation_),
  * Biconjugate Gradients (stabilized) solver
    (`BiCGSTAB`, for the general _non-singular_ linear problems),
  * Minimal Residual solver
    (`MINRES`, for the indefinite _symmetric_ linear problems),
  * _(planned)_ Generalized Minimal Residual method solver
    (`GMRES`, for the general linear problems);

- Matrix-free linear iterative least squares solvers:
  * Least squares-QR solver:
    (`LSQR`, for the general linear least squares problems),
  * Least squares-MINRES solver:
    (`LSMR`, for the general linear least squares problems);

- Matrix-free eigensolvers:
  * _(planned)_ Power iteration eigensolver
    (`PowerIteration`, for the general problems),
  * _(planned)_ Lanczos eigensolver
    (`Lanczos`, for the _symmetric_ problems,
     coupled with optional Lanczos-CG linear solver),
  * _(planned)_ Arnoldi eigensolver
    (`Arnoldi`, for the problems);

- Matrix-free preconditioners:
  * Block Jacobi preconditioner
    (`Jacobi`),
  * Block LU-SGS preconditioner
    (`LU_SGS`),
  * _(planned)_ Block SPAI preconditioner.

- Some very specific solvers:
  * _(planned)_ Symmetric tridiagonal eigensolvers
    (`bisection`, `divide-and-conquer`, `QR`,
     for the symmetric tridiagonal matrices, ),
  * LAPACK Symmetric tridiagonal eigensolvers
    (for the symmetric tridiagonal matrices, from [LAPACK](https://bit.ly/3yWg8qM)),

<!--
- Linear direct solvers (embedded into the matrix-free environment):
  * MKL Direct Sparse Solver 
    (`DSS_MKL`, from [MKL DSS](https://intel.ly/37N95pe)).-->

<!----------------------------------------------------------------->
## Road map
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
  - [ ] 💄 Untyped C API,
  - [ ] 💄 Redesign the C API,
  - [x] 💄 Unified Fortran solver API,
  - [x] 💄 Unified Fortran solver RCI API,
  - [ ] 🐏 PThread and semaphores for Windows,
  - [x] 🚬🪦 Unified RCI C solver API,
  - [ ] 🚬 MATLAB API,
  - [ ] 🚬 Lua API.

* Mesh:
  - [ ] 🪓 Nodes,
  - [ ] 🪓 Faces,
  - [ ] 🧸 Redesigned VTK output,
  - [ ] 🪓 Cache-friendly cell sorting,
  - [ ] 🚬 GMG,
  - [ ] 🚬 Non-conforming multilevel mesh,
  - [ ] 🚬🚬 Cut cell methods.

* General architecture:
  - [ ] 🦜🐞 Segfaults,
  - [ ] 🚬🪦 Custom threading (without OpenMP on Fortran level),
  - [ ] 🚬 GPU support,
  - [ ] 🚬🚬 MPI support.

* New differential operators and boundary conditions:
  - [ ] 🧸 Variable weight Laplace operator with 4+ order approx.,
  - [ ] 🦜 Tensor weight Laplace operator with 4+ order approx.,
  - [ ] 🦜 High order convection approx.,
  - [x] 🪓 Cylindrical symmetry for 2D domains,
  - [ ] 💄 More special boundary conditions,
  - [ ] 🪓 Godunov/WENO linear convection operator,
  - [ ] 🪓 Godunov/WENO nonlinear convection operator,
  - [ ] 🪓 Riemann solvers, Euler equations...

* Matrix reconstruction:
  - [ ] 🚬 Matrix portrait construction,
  - [ ] 🚬🚬🚬 Graph-coloring problem,
  - [ ] 🚬 Some matrix API...

* New linear solvers/eigensolvers:
  - [ ] 🚬 `GMRES` solver implementation.
  - [ ] 🧸 Recover existing eigensolvers.
  - [ ] 🧸 `PowerIterations` eigensolver,
  - [ ] 🪓 `Lanczos` hybrid eigensolver/solver,
  - [ ] 🪓 `Arnoldi` hybrid eigensolver/solver,
  - [ ] 🪓 `ConjugateChebyshev` solver,
  - [ ] 🦜 Custom eigensolvers for the tridiagonal Hermitian matrices,
  - [ ] 🧸 LAPACK eigensolvers for the Hessenberg matrices,
  - [ ] 🦜 Custom eigensolvers for the Hessenberg matrices,
  - [ ] 🦜 Direct solver (`MKL_DSS`, `PARDISO`, `SuperLU`)...

* Nonlinear solvers:
  - [x] 🧸 Newton-Raphson solver,
  - [ ] 🦜 Relaxed Newton solver,
  - [x] 🧸 Jacobian-Free Newton-Raphson solver,
  - [ ] 🧸 Select an epsilon in JFNK,
  - [ ] 🧻 Refactor JFNK without using the general Newton solver,
  - [ ] 🧸 Second order JFNK solver.

* Preconditioning:
  - [ ] 🦜 Fix block case diagonal extraction,
  - [ ] 🦜 Block tridiagonal preconditioner,
  - [ ] 🚬 Matrix-free SPAI preconditioner,
  - [ ] 🚬 Matrix-based preconditioning,
  - [ ] 🦜 Nnlinear preconditioner.

* Add support for complex numbers:
  - [x] 🪓 Support for complex numbers on BLAS level in Fortran,
  - [ ] 🪓 Support for complex linear solvers,
  - [ ] 🪓 Support for complex differential operators and boundary condition. 

* Symbolic arithmetics:
  - [ ] 🧸 Recover existing symbolic interface.
  - [ ] 🪓 Implement symbolic drivers for the linear case.
  - [ ] 🪓 Implement symbolic drivers for the nonlinear case.
  - [ ] 🦜 Direct/reverse auto-differentiation..
