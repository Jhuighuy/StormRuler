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
**StormRuler** contains:
* Conjugate Gradients solver 
  (_CG_, for the definite symmetric linear problems),
* Conjugate Gradients solver 
  (_CG_MKL_, from [MKL RCI ISS](https://intel.ly/3s4XF9F)),
* _(planned)_ MINRES solver
  (for the indefinite symmetric linear problems),
* Biconjugate Gradients (stabilized) solver
  (_BiCGSTAB_, for the general linear problems),
* _(planned)_ Generalized minimal residual method solver
  (_GMRES_, for the general SLAE),
* Flexible generalized minimal residual method solver
  (_FGMRES_MKL_, from [MKL RCI ISS](https://intel.ly/3s4XF9F)),
* _(planned)_ Newton-Raphson solver 
  (for the general nonlinear problems),
* _(planned)_ Jacobian-Free Newton-Raphson solver 
  (for the general nonlinear problems),
