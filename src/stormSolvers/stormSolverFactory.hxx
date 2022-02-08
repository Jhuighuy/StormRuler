/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// Copyright (C) 2022 Oleg Butakov
///
/// Permission is hereby granted, free of charge, to any person
/// obtaining a copy of this software and associated documentation
/// files (the "Software"), to deal in the Software without
/// restriction, including without limitation the rights  to use,
/// copy, modify, merge, publish, distribute, sublicense, and/or
/// sell copies of the Software, and to permit persons to whom the
/// Software is furnished to do so, subject to the following
/// conditions:
///
/// The above copyright notice and this permission notice shall be
/// included in all copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
/// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
/// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
/// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
/// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
/// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
/// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
/// OTHER DEALINGS IN THE SOFTWARE.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
#ifndef _STORM_SOLVER_FACTORY_HXX_
#define _STORM_SOLVER_FACTORY_HXX_

#include <string_view>
#include <stdexcept>

#include <stormBase.hxx>

#include <stormSolvers/stormSolver.hxx>
#include <stormSolvers/stormSolverCg.hxx>
#include <stormSolvers/stormSolverMinres.hxx>
#include <stormSolvers/stormSolverCgs.hxx>
#include <stormSolvers/stormSolverBiCgStab.hxx>
#include <stormSolvers/stormSolverTfqmr.hxx>
#include <stormSolvers/stormSolverIdrs.hxx>
#include <stormSolvers/stormSolverGmres.hxx>
//#include <stormSolvers/stormSolverLsqr.hxx>

#include <stormSolvers/stormSolverRichardson.hxx>
#include <stormSolvers/stormSolverNewton.hxx>

_STORM_NAMESPACE_BEGIN_

/// ----------------------------------------------------------------- ///
/// @brief Solver types.
/// ----------------------------------------------------------------- ///
namespace SolverType {
  
  /// @brief @c CG iterative solver.
  static std::string_view const Cg = "CG";
  /// @brief @c FCG iterative solver.
  static std::string_view const Fcg = "FCG";
  /// @brief @c MINRES iterative solver.
  static std::string_view const Minres = "MINRES";
  /// @brief @c CGS iterative solver.
  static std::string_view const Cgs = "CGS";
  /// @brief @c BiCGStab iterative solver.
  static std::string_view const BiCgStab = "BiCgStab";
  /// @brief @c BiCGStab(l) iterative solver.
  static std::string_view const BiCgStabL = "BiCgStab(l)";
  /// @brief @c TFQMR iterative solver.
  static std::string_view const Tfqmr = "TFQMR";
  /// @brief @c TFQMR(1) iterative solver.
  static std::string_view const Tfqmr1 = "TFQMR(1)";
  /// @brief @c IDR(s) iterative solver.
  static std::string_view const Idrs = "IDR(s)";
  /// @brief @c GMRES iterative solver.
  static std::string_view const Gmres = "GMRES";
  /// @brief @c FGMRES iterative solver.
  static std::string_view const Fgmres = "FGMRES";
  /// @brief @c LGMRES iterative solver.
  static std::string_view const Lgmres = "LGMRES";
  /// @brief @c LFGMRES iterative solver.
  static std::string_view const Lfgmres = "LFGMRES";
  
  /// @brief @c LSQR iterative solver.
  static std::string_view const Lsqr = "LSQR";
  /// @brief @c LSMR iterative solver.
  static std::string_view const Lsmr = "LSMR";

  /// @brief @c Richardson iterative solver.
  static std::string_view const Richarson = "Richardson";
  /// @brief @c Broyden iterative solver.
  static std::string_view const Broyden = "Broyden";
  /// @brief @c Newton iterative solver.
  static std::string_view const Newton = "Newton";
  /// @brief @c JFNK iterative solver.
  static std::string_view const Jfnk = "JFNK";

} // namespace SolverType

/// ----------------------------------------------------------------- ///
/// @brief Make iterative solver of the specified type.
/// ----------------------------------------------------------------- ///
template<class InVector, class OutVector = InVector>
std::unique_ptr<IterativeSolver<InVector, OutVector>>
    MakeIterativeSolver(std::string_view const& solverType) {

  // ----------------------
  // Try the square Krylov subspace solver first:
  // ----------------------
  if constexpr (std::is_same_v<InVector, OutVector>) {
    if (solverType == SolverType::Cg) {
      return std::make_unique<CgSolver<InVector>>();
    }
    if (solverType == SolverType::Fcg) {
      //return std::make_unique<FcgSolver<InVector>>();
    }
    if (solverType == SolverType::Minres) {
      return std::make_unique<MinresSolver<InVector>>();
    }
    if (solverType == SolverType::Cgs) {
      return std::make_unique<CgsSolver<InVector>>();
    }
    if (solverType == SolverType::BiCgStab) {
      return std::make_unique<BiCgStabSolver<InVector>>();
    }
    if (solverType == SolverType::BiCgStabL) {
      return std::make_unique<BiCGStabLSolver<InVector>>();
    }
    if (solverType == SolverType::Tfqmr) {
      return std::make_unique<TfqmrSolver<InVector>>();
    }
    if (solverType == SolverType::Tfqmr1) {
      return std::make_unique<Tfqmr1Solver<InVector>>();
    }
    if (solverType == SolverType::Idrs) {
      return std::make_unique<IdrsSolver<InVector>>();
    }
    if (solverType == SolverType::Gmres) {
      return std::make_unique<GmresSolver<InVector>>();
    }
    if (solverType == SolverType::Fgmres) {
      return std::make_unique<FgmresSolver<InVector>>();
    }
    if (solverType == SolverType::Lgmres) {
      //return std::make_unique<LgmresSolver<InVector>>();
    }
    if (solverType == SolverType::Lfgmres) {
      //return std::make_unique<LfgmresSolver<InVector>>();
    }
  }

  // ----------------------
  // Try the Krylov subspace least squares solvers next:
  // ----------------------
  if (solverType == SolverType::Lsqr) {
    //return std::make_unique<LsqrSolver<InVector, OutVector>>();
  }
  if (solverType == SolverType::Lsmr) {
    //return std::make_unique<LsmrSolver<InVector, OutVector>>();
  }

  // ----------------------
  // Finally, try the other general solvers:
  // ----------------------
  if constexpr (std::is_same_v<InVector, OutVector>) {
    if (solverType == SolverType::Richarson) {
      return std::make_unique<RichardsonSolver<InVector>>();
    }
    if (solverType == SolverType::Broyden) {
      //return std::make_unique<BroydenSolver<InVector>>();
    }
    if (solverType == SolverType::Newton) {
      //return std::make_unique<NewtonSolver<InVector>>();
    }
    if (solverType == SolverType::Jfnk) {
      return std::make_unique<JfnkSolver<InVector>>();
    }
  }

  throw std::invalid_argument("Invalid iterative solver type specified.");

} // MakeSolver<...>

_STORM_NAMESPACE_END_

#endif // ifndef _STORM_SOLVER_FACTORY_HXX_
