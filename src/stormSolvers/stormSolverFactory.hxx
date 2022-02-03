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

#include <cstring>
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
#define STORM_KSP_Richardson "Richardson"   /// @todo Implement me!
#define STORM_KSP_CG         "CG"
#define STORM_KSP_FCG        "FCG"          /// @todo Implement me!
#define STORM_KSP_MINRES     "MINRES"
#define STORM_KSP_CGS        "CGS"
#define STORM_KSP_BiCGStab   "BiCGStab"
#define STORM_KSP_BiCGStabL  "BiCGStab(l)"
#define STORM_KSP_TFQMR      "TFQMR"
#define STORM_KSP_TFQMR1     "TFQMR(1)"
#define STORM_KSP_IDRs       "IDR(s)"
#define STORM_KSP_GMRES      "GMRES"
#define STORM_KSP_FGMRES     "FGMRES"
#define STORM_KSP_LGMRES     "LGMRES"       /// @todo Implement me!
#define STORM_KSP_LFGMRES    "LFGMRES"      /// @todo Implement me!
#define STORM_KSP_LSQR       "LSQR"
#define STORM_KSP_LSMR       "LSMR"
/// @}

/// @{
#define STORM_NEWTON "Newton"
#define STORM_JFNK   "JFNK"
/// @}

#define _STORM_MAKE_KSP_SOLVER_(solverType) \
  if (std::strcmp(solverType, STORM_KSP_Richardson) == 0) { \
    \
    return std::make_unique<stormRichardsonSolver<Vector>>(); \
    \
  } else if (std::strcmp(solverType, STORM_KSP_CG) == 0) { \
    \
    return std::make_unique<stormCgSolver<Vector>>(); \
    \
  } else if (std::strcmp(solverType, STORM_KSP_MINRES) == 0) { \
    \
    return std::make_unique<stormMinresSolver<Vector>>(); \
    \
  } else if (std::strcmp(solverType, STORM_KSP_CGS) == 0) { \
    \
    return std::make_unique<stormCgsSolver<Vector>>(); \
    \
  } else if (std::strcmp(solverType, STORM_KSP_BiCGStab) == 0) { \
    \
    return std::make_unique<stormBiCgStabSolver<Vector>>(); \
    \
  } else if (std::strcmp(solverType, STORM_KSP_BiCGStabL) == 0) { \
    \
    return std::make_unique<stormBiCGStabLSolver<Vector>>(); \
    \
  } else if (std::strcmp(solverType, STORM_KSP_TFQMR) == 0) { \
    \
    return std::make_unique<stormTfqmrSolver<Vector>>(); \
    \
  } else if (std::strcmp(solverType, STORM_KSP_TFQMR1) == 0) { \
    \
    return std::make_unique<stormTfqmr1Solver<Vector>>(); \
    \
  } else if (std::strcmp(solverType, STORM_KSP_IDRs) == 0) { \
    \
    return std::make_unique<stormIdrsSolver<Vector>>(); \
    \
  } else if (std::strcmp(solverType, STORM_KSP_GMRES) == 0) { \
    \
    return std::make_unique<stormGmresSolver<Vector>>(); \
    \
  } else if (std::strcmp(solverType, STORM_KSP_FGMRES) == 0) { \
    \
    return std::make_unique<stormFgmresSolver<Vector>>(); \
    \
  } /*else if (std::strcmp(solverType, STORM_KSP_LSQR) == 0) { \
    \
    return std::make_unique<stormLsqrSolver<Vector>>(); \
    \
  } else if (std::strcmp(solverType, STORM_KSP_LSMR) == 0) { \
    \
    return std::make_unique<stormLsmrSolver<Vector>>(); \
    \
  } */

#define _STORM_MAKE_RECT_KSP_SOLVER_(solverType) \
  /*if (std::strcmp(solverType, STORM_KSP_LSQR) == 0) { \
    \
    return std::make_unique<stormLsqrSolver<InVector, OutVector>>(); \
    \
  } else if (std::strcmp(solverType, STORM_KSP_LSMR) == 0) { \
    \
    return std::make_unique<stormLsmrSolver<InVector, OutVector>>(); \
    \
  }*/

/// ----------------------------------------------------------------- ///
/// @brief Make iterative solver of the specified type.
/// ----------------------------------------------------------------- ///
/** @{ */
template<class Vector>
std::unique_ptr<stormIterativeSolver<Vector>>
    stormMakeIterativeSolver(stormString_t solverType) {

  _STORM_MAKE_KSP_SOLVER_(solverType);

  throw std::invalid_argument("Invalid iterative solver type specified.");

} // stormMakeSolver<...>
template<class InVector, class OutVector>
std::unique_ptr<stormIterativeSolver<InVector, OutVector>>
    stormMakeIterativeSolver(stormString_t solverType) {

  if constexpr (std::is_same_v<InVector, OutVector>) {
    return stormMakeIterativeSolver<InVector>(solverType);
  }

  _STORM_MAKE_RECT_KSP_SOLVER_(solverType);

  throw std::invalid_argument("Invalid iterative rectangular solver type specified.");

} // stormMakeSolver<...>
/** @} */

#endif // ifndef _STORM_SOLVER_FACTORY_HXX_
