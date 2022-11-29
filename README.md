<!--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-->
# StormRuler🦜 — A very high order FVM framework
<!--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-->

[![CodeFactor](https://www.codefactor.io/repository/github/jhuighuy/stormruler/badge)](https://www.codefactor.io/repository/github/jhuighuy/stormruler)

**StormRuler** is a FVM-based multidimensional partial 
differential equations solving framework, written in C++23.

<!----------------------------------------------------------------------------->
## 🌀 Equations solved
<!----------------------------------------------------------------------------->

**StormRuler** can be used to solve various partial differential equations, 
including:
* 🌊 Incompressible Navier-Stokes equations,
* 🌪 _(planned)_ Сompressible Navier-Stokes/Euler equations,
* 💧 _(planned)_ Cahn-Hilliard equation,
* ...

<!----------------------------------------------------------------------------->
## 🌐 Numerical methods
<!----------------------------------------------------------------------------->
The heart of the **StormRuler** is the _✨Finite Volume Method✨_.

<!----------------------------------------------------------------------------->
## 🌈 Algebraic solvers
<!----------------------------------------------------------------------------->
In order to implement the high performance implicit schemes, several algebraic 
problems, like systems of linear and nonlinear equations, have to be solved.

For the sake of convenience, all auxiliary solvers are implemented  in the 
_matrix-free_ manner: no assembled matrix is required to find a solution of the 
algebraic problem, only the matrix-vector product function is used.

Although most of the problems can be solved in the matrix-free manner using the 
Krylov subspace iterative solver, in some cases an assembled matrix be required 
to construct a suitable preconditioner or utilize a direct solver.
**StormRuler** extracts a matrix from the matrix-vector product function 
automatically using the sophisticated template metaprogramming techniques.

### Iterative solvers:
| Name                   | Problem type                | Flexible | Status   |
|------------------------|-----------------------------|----------|----------|
| **Richardson**         | General Square Nonsingular  | Yes      | Complete |
| **Broyden**            | General Square Nonsingular  | No       | Planned  |
| **Newton**             | General Square Nonsingular  | Yes      | Complete |
| **JFNK**               | General Square Nonsingular  | No       | Complete |
| **CG**                 | Linear Definite Symmetric   | No       | Complete |
| **FCG**                | Linear Definite Symmetric   | Yes      | Planned  |
| **MINRES**             | Linear Indefinite Symmetric | No       | Planned  |
| **CGS**                | Linear Square Nonsingular   | No       | Complete |
| **BiCGStab**           | Linear Square Nonsingular   | No       | Complete |
| **BiCGStab(l)**        | Linear Square Nonsingular   | No       | Complete |
| **TFQMR**              | Linear Square Nonsingular   | No       | Complete |
| **TFQMR(1)**           | Linear Square Nonsingular   | No       | Complete |
| **IDR(s)**             | Linear Square Nonsingular   | No       | Complete |
| **GMRES**              | Linear Square               | No       | Complete |
| **FGMRES**             | Linear Square               | Yes      | Complete |
| **LGMRES**             | Linear Square               | No       | Planned  |
| **LFGMRES**            | Linear Square               | Yes      | Planned  |
| **LSQR**               | Linear Rectangular          | No       | Planned  |
| **LSMR**               | Linear Rectangular          | No       | Planned  |

### Preconditioners:
| Name                   | Problem type                | Flexible | Status   |
|------------------------|-----------------------------|----------|----------|
| Diagonal               |                             | No       | Planned  |
| Symmetric Gauss-Seidel |                             | No       | Planned  |
| Incomplete Cholesky    |                             | No       | Planned  |
| Incomplete LU          |                             | No       | Planned  |
| AINV                   |                             | No       | Planned  |
| SPAI                   |                             | No       | Planned  |
| Krylov                 |                             | Yes      | Planned  |

<!----------------------------------------------------------------------------->
## 🏗 Compiling
<!----------------------------------------------------------------------------->

_To be written..._
