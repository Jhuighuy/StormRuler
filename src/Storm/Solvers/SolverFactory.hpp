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

#include <concepts>
#include <memory>
#include <stdexcept>

#include <Storm/Base.hpp>

#include <Storm/Utils/Enum.hpp>

#include <Storm/Solvers/Solver.hpp>
// clang-format off
#include <Storm/Solvers/SolverCg.hpp>
#include <Storm/Solvers/SolverMinres.hpp>
#include <Storm/Solvers/SolverCgs.hpp>
#include <Storm/Solvers/SolverBiCgStab.hpp>
#include <Storm/Solvers/SolverTfqmr.hpp>
#include <Storm/Solvers/SolverIdrs.hpp>
#include <Storm/Solvers/SolverGmres.hpp>
#include <Storm/Solvers/SolverLsqr.hpp>
#include <Storm/Solvers/SolverRichardson.hpp>
#include <Storm/Solvers/SolverNewton.hpp>
// clang-format on

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

/// @brief Make iterative solver of the specified @p solver_type.
template<VectorLike InVector, VectorLike OutVector = InVector>
auto MakeIterativeSolver(SolverType solver_type = SolverType::Default)
    -> std::unique_ptr<IterativeSolver<InVector, OutVector>> {
  // Try the Krylov subspace square solver first:
  if constexpr (std::same_as<InVector, OutVector>) {
    if (solver_type == SolverType::Cg) {
      return std::make_unique<CgSolver<InVector>>();
    }
    if (solver_type == SolverType::Fcg) {
      // return std::make_unique<FcgSolver<InVector>>();
    }
    if (solver_type == SolverType::Minres) {
      // return std::make_unique<MinresSolver<InVector>>();
    }
    if (solver_type == SolverType::Cgs) {
      return std::make_unique<CgsSolver<InVector>>();
    }
    if (solver_type == SolverType::BiCgStab) {
      return std::make_unique<BiCgStabSolver<InVector>>();
    }
    if (solver_type == SolverType::BiCgStabL) {
      return std::make_unique<BiCgStabLSolver<InVector>>();
    }
    if (solver_type == SolverType::Tfqmr) {
      return std::make_unique<TfqmrSolver<InVector>>();
    }
    if (solver_type == SolverType::Tfqmr1) {
      return std::make_unique<Tfqmr1Solver<InVector>>();
    }
    if (solver_type == SolverType::Idrs) {
      return std::make_unique<IdrsSolver<InVector>>();
    }
    if (solver_type == SolverType::Default ||
        solver_type == SolverType::Gmres) {
      // Note: GMRES is the default square solver.
      return std::make_unique<GmresSolver<InVector>>();
    }
    if (solver_type == SolverType::Fgmres) {
      return std::make_unique<FgmresSolver<InVector>>();
    }
    if (solver_type == SolverType::Lgmres) {
      // return std::make_unique<LgmresSolver<InVector>>();
    }
    if (solver_type == SolverType::Lfgmres) {
      // return std::make_unique<LfgmresSolver<InVector>>();
    }
  }

  // Next, try the other general square solvers:
  if constexpr (std::same_as<InVector, OutVector>) {
    if (solver_type == SolverType::Richarson) {
      return std::make_unique<RichardsonSolver<InVector>>();
    }
    if (solver_type == SolverType::Broyden) {
      // return std::make_unique<BroydenSolver<InVector>>();
    }
    if (solver_type == SolverType::Newton) {
      return std::make_unique<NewtonSolver<InVector>>();
    }
    if (solver_type == SolverType::Jfnk) {
      return std::make_unique<JfnkSolver<InVector>>();
    }
  }

  // Finally, try the Krylov subspace least squares solvers:
  if (solver_type == SolverType::Lsqr) {
    // return std::make_unique<LsqrSolver<InVector, OutVector>>();
  }
  if (solver_type == SolverType::Default || solver_type == SolverType::Lsmr) {
    // Note: LSMR is the default rectangular solver.
    // return std::make_unique<LsmrSolver<InVector, OutVector>>();
  }

  throw std::invalid_argument("Invalid iterative solver type specified.");

} // MakeIterativeSolver

} // namespace Storm
