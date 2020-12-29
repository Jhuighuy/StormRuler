# StormRuler
**StormRuler** is a multidimensional very high order CFD solver, 
written in modern Fortran.

<!----------------------------------------------------------------->
## Equations solved
**StormRuler** features support of the various set of the
partial differential equations, including:
* Incompressible Navier-Stokes equations
  $$
    \frac{\partial\vec{v}}{\partial t} +
    (\vec{v}\cdot\nabla)\vec{v} + 
    \frac{1}{\rho}\nabla p = 
    \nu\Delta\vec{v},\quad
    \nabla\cdot\vec{v} = 0,
  $$
* Cahn-Hilliard equation
  $$
    \frac{\partial c}{\partial t} =
    \nabla\cdot M(c)\nabla\left(
      f(c) - \varepsilon^{2}\Delta c
    \right),
  $$
* _(planned)_ Streamfunction-Vorticity formulation for the
  Incompressible Navier-Stokes equations in 2D,
* Incompressible Navier-Stokes-Cahn-Hilliard equations
  $$
    \frac{\partial\vec{v}}{\partial t} +
    (\vec{v}\cdot\nabla)\vec{v} + 
    \frac{1}{\rho}\nabla p + \frac{1}{\rho}c\nabla w = 
    \nu\Delta\vec{v},\quad
    \nabla\cdot\vec{v} = 0, \\
    \frac{\partial c}{\partial t} +
    (\vec{v}\cdot\nabla)c =
    \nabla\cdot M(c)\nabla w,\quad
    w = f(c) - \varepsilon^{2}\Delta c,
  $$
* _(planned)_ Streamfunction-Vorticity formulation for the
  Incompressible Navier-Stokes-Cahn-Hilliard equations in 2D.
<!----------------------------------------------------------------->

<!----------------------------------------------------------------->
## Numerical methods
The heart of the **StormRuler** is the Finite Difference Method.
<!----------------------------------------------------------------->

<!----------------------------------------------------------------->
## Auxiliary solvers
In order to implement the high performance implicit schemes,
several auxiliary problems, like systems of linear and nonlinear
equations, have to be solved.
**StormRuler** contains:
* Conjugate Gradients solver 
  (for the definite symmetric linear problems),
* Biconjugate Gradients (stabilized) solver
  (for the general linear problems),
* _(planned)_ MINRES solver
  (for the indefinite symmetric linear problems),
* _(planned)_ GMRES solver
  (for the general SLAE),
* _(planned)_ Newton-Raphson solver 
  (for the general nonlinear problems),
* _(planned)_ Jacobian-Free Newton-Raphson solver 
  (for the general nonlinear problems),
<!----------------------------------------------------------------->
