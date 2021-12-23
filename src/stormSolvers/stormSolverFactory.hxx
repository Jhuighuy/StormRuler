/// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ///
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
/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ///
#ifndef _STORM_SOLVER_FACTORY_HXX_
#define _STORM_SOLVER_FACTORY_HXX_

#include <stormSolvers/stormSolver.hxx>
#include <stormSolvers/stormSolverCg.hxx>
#include <stormSolvers/stormSolverMinres.hxx>
#include <stormSolvers/stormSolverLsqr.hxx>

#include <cstring>
#include <stdexcept>

/// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ///
/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ///

/// Krylov-subspace solver types.
/// @{
#define STORM_CG       "CG"
#define STORM_BiCGStab "BiCGStab"
#define STORM_MINRES   "MINRES"
#define STORM_GMRES    "GMRES"
#define STORM_TFQMR    "TFQMR"
#define STORM_LSQR     "LSQR"
#define STORM_LSMR     "LSMR"
/// @}

/// ----------------------------------------------------------------- ///
/// @brief Make solver of the specified type.
/// ----------------------------------------------------------------- ///
/** @{ */
template<class tArray, class... tArgs>
std::unique_ptr<stormSolver<tArray>> 
    stormMakeSolver(stormString_t solverType, tArgs&&... args) {

  if (std::strcmp(solverType, STORM_CG) == 0) {

    return std::make_unique<stormCgSolver<tArray>>(std::forward<tArgs>(args)...);

  } else if (std::strcmp(solverType, STORM_BiCGStab) == 0) {

    return std::make_unique<stormBiCgStabSolver<tArray>>(std::forward<tArgs>(args)...);

  } else if (std::strcmp(solverType, STORM_MINRES) == 0) {

    return std::make_unique<stormMinresSolver<tArray>>(std::forward<tArgs>(args)...);

  } else if (std::strcmp(solverType, STORM_GMRES) == 0) {

    return std::make_unique<stormGmresSolver<tArray>>(std::forward<tArgs>(args)...);

  } else if (std::strcmp(solverType, STORM_TFQMR) == 0) {

    return std::make_unique<stormTfqmrSolver<tArray>>(std::forward<tArgs>(args)...);

  } else if (std::strcmp(solverType, STORM_LSQR) == 0) {

    return std::make_unique<stormLsqrSolver<tArray>>(std::forward<tArgs>(args)...);

  } else if (std::strcmp(solverType, STORM_LSMR) == 0) {

    return std::make_unique<stormLsmrSolver<tArray>>(std::forward<tArgs>(args)...);

  } 

  throw std::invalid_argument("Invalid solver type specified.");

} // stormMakeSolver<...>
template<class tInArray, class tOutArray, class... tArgs>
std::unique_ptr<stormSolver<tInArray, tOutArray>> 
    stormMakeSolver(stormString_t solverType, tArgs&&... args) {

  if constexpr (std::is_same_v<tInArray, tOutArray>) {

    return stormMakeSolver<tInArray>(solverType, std::forward<tArgs>(args)...);

  } else if (std::strcmp(solverType, STORM_LSQR) == 0) {

    return std::make_unique<stormLsqrSolver<tInArray, tOutArray>>(std::forward<tArgs>(args)...);

  } else if (std::strcmp(solverType, STORM_LSMR) == 0) {

    return std::make_unique<stormLsmrSolver<tInArray, tOutArray>>(std::forward<tArgs>(args)...);

  } 

  throw std::invalid_argument("Invalid rectangular solver type specified.");

} // stormMakeSolver<...>
/** @} */

/// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ///
/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ///

#endif // ifndef _STORM_SOLVER_FACTORY_HXX_
