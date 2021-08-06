# StormRuler
**StormRuler** is a very high order multidimensional CFD solver, 
written in modern Fortran.

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
