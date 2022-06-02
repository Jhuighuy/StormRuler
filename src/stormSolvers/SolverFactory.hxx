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

#pragma once

#include <stdexcept>

#include <stormBase.hxx>
#include <stormUtils/Enum.hxx>

#include <stormSolvers/Solver.hxx>
#include <stormSolvers/SolverCg.hxx>
//#include <stormSolvers/SolverMinres.hxx>
#include <stormSolvers/SolverBiCgStab.hxx>
#include <stormSolvers/SolverCgs.hxx>
#include <stormSolvers/SolverGmres.hxx>
#include <stormSolvers/SolverIdrs.hxx>
#include <stormSolvers/SolverTfqmr.hxx>
//#include <stormSolvers/SolverLsqr.hxx>

#include <stormSolvers/SolverNewton.hxx>
#include <stormSolvers/SolverRichardson.hxx>

namespace Storm {

/// ----------------------------------------------------------------- ///
/// @brief Solver types.
/// ----------------------------------------------------------------- ///
class SolverType final : public Enum<SolverType> {
  // clang-format off
  StormEnum_(SolverType)

  /// @brief Default solver.
  StormEnumValue_(Default)

  /// @brief @c CG iterative solver.
  StormEnumValue_(Cg, "CG")

  /// @brief @c FCG iterative solver.
  StormEnumValue_(Fcg, "FCG")

  /// @brief @c MINRES iterative solver.
  StormEnumValue_(Minres, "MINRES")

  /// @brief @c CGS iterative solver.
  StormEnumValue_(Cgs, "CGS")

  /// @brief @c BiCGStab iterative solver.
  StormEnumValue_(BiCgStab, "BiCgStab")

  /// @brief @c BiCGStab(l) iterative solver.
  StormEnumValue_(BiCgStabL, "BiCgStab(l)")

  /// @brief @c TFQMR iterative solver.
  StormEnumValue_(Tfqmr, "TFQMR")

  /// @brief @c TFQMR(1) iterative solver.
  StormEnumValue_(Tfqmr1, "TFQMR(1)")

  /// @brief @c IDR(s) iterative solver.
  StormEnumValue_(Idrs, "IDR(s)")

  /// @brief @c GMRES iterative solver.
  StormEnumValue_(Gmres, "GMRES")

  /// @brief @c FGMRES iterative solver.
  StormEnumValue_(Fgmres, "FGMRES")

  /// @brief @c LGMRES iterative solver.
  StormEnumValue_(Lgmres, "LGMRES")

  /// @brief @c LFGMRES iterative solver.
  StormEnumValue_(Lfgmres, "LFGMRES")

  /// @brief @c LSQR iterative solver.
  StormEnumValue_(Lsqr, "LSQR")

  /// @brief @c LSMR iterative solver.
  StormEnumValue_(Lsmr, "LSMR")

  /// @brief @c Richardson iterative solver.
  StormEnumValue_(Richarson)

  /// @brief @c Broyden iterative solver.
  StormEnumValue_(Broyden)

  /// @brief @c Newton iterative solver.
  StormEnumValue_(Newton)

  /// @brief @c JFNK iterative solver.
  StormEnumValue_(Jfnk, "JFNK")

  // clang-format on

}; // class SolverType

/// ----------------------------------------------------------------- ///
/// @brief Make iterative solver of the specified type.
/// ----------------------------------------------------------------- ///
template<vector_like InVector, vector_like OutVector = InVector>
auto MakeIterativeSolver(SolverType solverType = SolverType::Default)
    -> std::unique_ptr<IterativeSolver<InVector, OutVector>> {
  // Try the square Krylov subspace solver first:
  if constexpr (std::is_same_v<InVector, OutVector>) {
    if (solverType == SolverType::Cg) {
      return std::make_unique<CgSolver<InVector>>();
    }
    if (solverType == SolverType::Fcg) {
      // return std::make_unique<FcgSolver<InVector>>();
    }
    if (solverType == SolverType::Minres) {
      // return std::make_unique<MinresSolver<InVector>>();
    }
    if (solverType == SolverType::Cgs) {
      return std::make_unique<CgsSolver<InVector>>();
    }
    if (solverType == SolverType::BiCgStab) {
      return std::make_unique<BiCgStabSolver<InVector>>();
    }
    if (solverType == SolverType::BiCgStabL) {
      return std::make_unique<BiCgStabLSolver<InVector>>();
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
    if (solverType == SolverType::Default || solverType == SolverType::Gmres) {
      // Note: GMRES is the default square solver.
      return std::make_unique<GmresSolver<InVector>>();
    }
    if (solverType == SolverType::Fgmres) {
      return std::make_unique<FgmresSolver<InVector>>();
    }
    if (solverType == SolverType::Lgmres) {
      // return std::make_unique<LgmresSolver<InVector>>();
    }
    if (solverType == SolverType::Lfgmres) {
      // return std::make_unique<LfgmresSolver<InVector>>();
    }
  }

  // Next, try the other general solvers:
  if constexpr (std::is_same_v<InVector, OutVector>) {
    if (solverType == SolverType::Richarson) {
      return std::make_unique<RichardsonSolver<InVector>>();
    }
    if (solverType == SolverType::Broyden) {
      // return std::make_unique<BroydenSolver<InVector>>();
    }
    if (solverType == SolverType::Newton) {
      return std::make_unique<NewtonSolver<InVector>>();
    }
    if (solverType == SolverType::Jfnk) {
      return std::make_unique<JfnkSolver<InVector>>();
    }
  }

  // Finally, try the Krylov subspace least squares solvers:
  if (solverType == SolverType::Lsqr) {
    // return std::make_unique<LsqrSolver<InVector, OutVector>>();
  }
  if (solverType == SolverType::Default || solverType == SolverType::Lsmr) {
    // Note: LSMR is the default rectangular solver.
    // return std::make_unique<LsmrSolver<InVector, OutVector>>();
  }

  throw std::invalid_argument("Invalid iterative solver type specified.");

} // MakeIterativeSolver

} // namespace Storm
