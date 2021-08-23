<!--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-->
# ♂StormRuler♂ — A very high order CFD solver
<!--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-->
**StormRuler** is a very high order multidimensional CFD solver, 
written in Fortran 2008 and C++17.

Supported compilers:
* _GCC_ version _9.0_ and more recent 
  (tested on _10.3.0_).
* _Intel Classic compilers_ version _19.0_ and more recent
  (tested on _2021.3.0_).
* _PGI Compilers_ (from _NVIDIA HPC SDK_) version 21 and more recent 
  (tested on _21.07_).
* _NAG Fortran Compiler_ version 7.0 and more recent
  (tested on _7.0 build 7048_).

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

<!--
For the sake of convenience, all auxiliary solvers are implemented 
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

**StormRuler** contains:
- Matrix-free Linear iterative solvers: 
  * Conjugate Gradients solver 
    (`CG`, for the definite symmetric linear problems),
  * MKL Conjugate Gradients solver 
    (`CG_MKL`, from [MKL RCI ISS](https://intel.ly/3s4XF9F)),
  * _(planned)_ MINRES solver
    (for the indefinite symmetric linear problems),
  * Biconjugate Gradients (stabilized) solver
    (`BiCGSTAB`, for the general linear problems),
  * _(planned)_ Generalized minimal residual method solver
    (`GMRES`, for the general SLAE),
  * MKL Flexible generalized minimal residual method solver
    (`FGMRES_MKL`, from [MKL RCI ISS](https://intel.ly/3s4XF9F)).

- Matrix-free preconditioners:
  * Block Jacobi preconditioner,
  * Block LU-SGS preconditioner.

<!--
- Linear direct solvers (embedded into the matrix-free environment):
  * MKL Direct Sparse Solver 
    (`DSS_MKL`, from [MKL DSS](https://intel.ly/37N95pe)).-->

- Nonlinear solvers:
  * _(planned)_ Newton-Raphson solver 
    (for the general nonlinear problems),
  * _(planned)_ Jacobian-Free Newton-Raphson solver 
    (for the general nonlinear problems),
