<!--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-->
# StormRuler🦜 — A very high order FVM framework
<!--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-->

[![SonarCloud](https://sonarcloud.io/images/project_badges/sonarcloud-white.svg)](https://sonarcloud.io/summary/new_code?id=Jhuighuy_StormRuler)

[![SonarCloud](https://github.com/Jhuighuy/StormRuler/actions/workflows/analysis-sonar.yml/badge.svg)](https://github.com/Jhuighuy/StormRuler/actions/workflows/analysis-sonar.yml)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=Jhuighuy_StormRuler&metric=alert_status)](https://sonarcloud.io/summary/new_code?id=Jhuighuy_StormRuler)
[![Coverage](https://sonarcloud.io/api/project_badges/measure?project=Jhuighuy_StormRuler&metric=coverage)](https://sonarcloud.io/summary/new_code?id=Jhuighuy_StormRuler)
[![Lines of Code](https://sonarcloud.io/api/project_badges/measure?project=Jhuighuy_StormRuler&metric=ncloc)](https://sonarcloud.io/summary/new_code?id=Jhuighuy_StormRuler)
[![codecov](https://codecov.io/github/Jhuighuy/StormRuler/branch/main/graph/badge.svg?token=GUSIUDW3G0)](https://codecov.io/github/Jhuighuy/StormRuler)
[![CodeQL](https://github.com/Jhuighuy/StormRuler/actions/workflows/analysis-codeql.yml/badge.svg)](https://github.com/Jhuighuy/StormRuler/actions/workflows/analysis-codeql.yml)
[![Coverity Scan Build Status](https://scan.coverity.com/projects/27159/badge.svg)](https://scan.coverity.com/projects/jhuighuy-stormruler)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/e7a26478673f403aa32f41a7c2a86d8d)](https://www.codacy.com/gh/Jhuighuy/StormRuler/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=Jhuighuy/StormRuler&amp;utm_campaign=Badge_Grade)
[![CodeFactor](https://www.codefactor.io/repository/github/jhuighuy/stormruler/badge)](https://www.codefactor.io/repository/github/jhuighuy/stormruler)
[![OpenSSF Best Practices](https://bestpractices.coreinfrastructure.org/projects/6812/badge)](https://bestpractices.coreinfrastructure.org/projects/6812)

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
_To be written..._

<!----------------------------------------------------------------------------->
## 🌐 Numerical methods
<!----------------------------------------------------------------------------->

The heart of the **StormRuler** is the _✨Finite Volume Method✨_.
_To be written..._

<!----------------------------------------------------------------------------->
## 🌈 Algebra
<!----------------------------------------------------------------------------->

_To be written..._

### Iterative solvers:
| Name                    | Problem type                 | Flexible | Status   |
|:------------------------|:-----------------------------|:--------:|:--------:|
| **Richardson**          | General Square Nonsingular   | Yes      | ✅       |
| **Broyden**             | General Square Nonsingular   | No       | Planned  |
| **Newton**              | General Square Nonsingular   | Yes      | ✅       |
| **JFNK**                | General Square Nonsingular   | No       | ✅       |
| **CG**                  | Linear Definite Symmetric    | No       | ✅       |
| **FCG**                 | Linear Definite Symmetric    | Yes      | Planned  |
| **MINRES**              | Linear Indefinite Symmetric  | No       | Planned  |
| **CGS**                 | Linear Square Nonsingular    | No       | ✅       |
| **BiCGStab**            | Linear Square Nonsingular    | No       | ✅       |
| **BiCGStab(l)**         | Linear Square Nonsingular    | No       | ✅       |
| **TFQMR**               | Linear Square Nonsingular    | No       | ✅       |
| **TFQMR(1)**            | Linear Square Nonsingular    | No       | ✅       |
| **IDR(s)**              | Linear Square Nonsingular    | No       | ✅       |
| **GMRES**               | Linear Square                | No       | ✅       |
| **FGMRES**              | Linear Square                | Yes      | ✅       |
| **LGMRES**              | Linear Square                | No       | Planned  |
| **LFGMRES**             | Linear Square                | Yes      | Planned  |
| **LSQR**                | Linear Rectangular           | No       | Planned  |
| **LSMR**                | Linear Rectangular           | No       | Planned  |

### Preconditioners:
| Name                        | Problem type             | Flexible | Status   |
|:----------------------------|:-------------------------|:--------:|:--------:|
| **Block Diagonal**          | Square Nonsingular       | No       | Planned  |
| **Symmetric Gauss-Seidel**  | Square Nonsingular       | No       | Planned  |
| **Incomplete Cholesky**     | Definite Symmetric       | No       | Planned  |
| **Incomplete LU**           | Square Nonsingular       | No       | Planned  |
| **Incomplete QR**           | Rectangular              | No       | Planned  |
| **AINV**                    | Definite Symmetric       | No       | Planned  |
| **SPAI**                    | Square Nonsingular       | No       | Planned  |
| **AMG**                     | Square Nonsingular       | No       | Planned  |
| **Krylov**                  | Square Nonsingular       | Yes      | Planned  |

<!----------------------------------------------------------------------------->
## 🏗 Compiling
<!----------------------------------------------------------------------------->

_To be written..._

| Compiler               | Linux           | macOS           | Windows         |
|:----------------------:|:---------------:|:---------------:|:---------------:|
| **GCC** 12.1+          | ✅              | ✅              | ✅              |
| **Clang** 16.0+        | Partial         | Planned         | Planned         |
| **Intel LLVM**         | Planned         |                 | Planned         |
| **MSVC** 19.34+        |                 |                 | ✅              |
