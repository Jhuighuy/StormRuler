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

#include <stormSolvers/stormSolver.hxx>
#include <stormSolvers/stormSolverRichardson.hxx>
#include <stormSolvers/stormSolverCg.hxx>
#include <stormSolvers/stormSolverMinres.hxx>
#include <stormSolvers/stormSolverCgs.hxx>
#include <stormSolvers/stormSolverBiCgStab.hxx>
#include <stormSolvers/stormSolverTfqmr.hxx>
#include <stormSolvers/stormSolverIdrs.hxx>
#include <stormSolvers/stormSolverGmres.hxx>
//#include <stormSolvers/stormSolverLsqr.hxx>

#include <stormSolvers/stormSolverNewton.hxx>

/// Krylov-subspace solver types.
/// @{
#define STORM_KSP_CG        "CG"
#define STORM_KSP_FCG       "FCG"
#define STORM_KSP_MINRES    "MINRES"
#define STORM_KSP_CGS       "CGS"
#define STORM_KSP_BiCGStab  "BiCGStab"
#define STORM_KSP_BiCGStabL "BiCGStab(l)"
#define STORM_KSP_TFQMR     "TFQMR"
#define STORM_KSP_TFQMR1    "TFQMR(1)"
#define STORM_KSP_IDRs      "IDR(s)"
#define STORM_KSP_GMRES     "GMRES"
#define STORM_KSP_FGMRES    "FGMRES"
#define STORM_KSP_LGMRES    "LGMRES"
#define STORM_KSP_LFGMRES   "LFGMRES"
#define STORM_KSP_LSQR      "LSQR"
#define STORM_KSP_LSMR      "LSMR"
/// @}

/// @{
#define STORM_Richardson "Richardson"
#define STORM_Newton     "Newton"
#define STORM_JFNK       "JFNK"
/// @}

/// ----------------------------------------------------------------- ///
/// @brief Make iterative solver of the specified type.
/// ----------------------------------------------------------------- ///
/** @{ */
template<class Vector>
std::unique_ptr<stormIterativeSolver<Vector>>
    stormMakeIterativeSolver(std::string_view const& solverType) {

  // ----------------------
  // Try the Krylov subspace solver first:
  // ----------------------
  if (solverType == STORM_KSP_CG) {

    return std::make_unique<stormCgSolver<Vector>>();

  } else if (solverType == STORM_KSP_FCG) {

    /// @todo Implement me!
    //return std::make_unique<stormFcgSolver<Vector>>();

  } else if (solverType == STORM_KSP_MINRES) {

    return std::make_unique<stormMinresSolver<Vector>>();

  } else if (solverType == STORM_KSP_CGS) {

    return std::make_unique<stormCgsSolver<Vector>>();

  } else if (solverType == STORM_KSP_BiCGStab) {

    return std::make_unique<stormBiCgStabSolver<Vector>>();

  } else if (solverType == STORM_KSP_BiCGStabL) {

    return std::make_unique<stormBiCGStabLSolver<Vector>>();

  } else if (solverType == STORM_KSP_TFQMR) {

    return std::make_unique<stormTfqmrSolver<Vector>>();

  } else if (solverType == STORM_KSP_TFQMR1) {

    return std::make_unique<stormTfqmr1Solver<Vector>>();

  } else if (solverType == STORM_KSP_IDRs) {

    return std::make_unique<stormIdrsSolver<Vector>>();

  } else if (solverType == STORM_KSP_GMRES) {

    return std::make_unique<stormGmresSolver<Vector>>();

  } else if (solverType == STORM_KSP_FGMRES) {

    return std::make_unique<stormFgmresSolver<Vector>>();

  } else if (solverType == STORM_KSP_LGMRES) {

    /// @todo Implement me!
    //return std::make_unique<stormLgmresSolver<Vector>>();

  } else if (solverType == STORM_KSP_LFGMRES) {

    /// @todo Implement me!
    //return std::make_unique<stormLfgmresSolver<Vector>>();

  }

  // ----------------------
  // Try the Krylov subspace least squares solvers next:
  // ----------------------
  if (solverType == STORM_KSP_LSQR) {

    //return std::make_unique<stormLsqrSolver<Vector>>();

  } else if (solverType == STORM_KSP_LSMR) {

    //return std::make_unique<stormLsmrSolver<Vector>>();

  }

  // ----------------------
  // Finally, try the other general solvers:
  // ----------------------
  if (solverType == STORM_Richardson) {

    return std::make_unique<stormRichardsonSolver<Vector>>();

  } else if (solverType == STORM_Newton) {

    //return std::make_unique<stormNewtonSolver<Vector>>();

  } else if (solverType == STORM_JFNK) {

    return std::make_unique<stormJfnkSolver<Vector>>();

  }

  throw std::invalid_argument("Invalid iterative solver type specified.");

} // stormMakeSolver<...>
template<class InVector, class OutVector>
std::unique_ptr<stormIterativeSolver<InVector, OutVector>>
    stormMakeIterativeSolver(std::string_view const& solverType) {

  // ----------------------
  // Try the square solvers first:
  // ----------------------
  if constexpr (std::is_same_v<InVector, OutVector>) {

    return stormMakeIterativeSolver<InVector>(solverType);

  }

  // ----------------------
  // Finally, try the Krylov subspace least squares solvers:
  // ----------------------
  if (solverType == STORM_KSP_LSQR) {

    //return std::make_unique<stormLsqrSolver<InVector, OutVector>>();

  } else if (solverType == STORM_KSP_LSMR) {

    //return std::make_unique<stormLsmrSolver<InVector, OutVector>>();

  }

  throw std::invalid_argument("Invalid iterative rectangular solver type specified.");

} // stormMakeSolver<...>
/** @} */

#endif // ifndef _STORM_SOLVER_FACTORY_HXX_
