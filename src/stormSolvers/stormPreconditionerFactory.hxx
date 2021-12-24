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
#ifndef _STORM_PRECONDITIONER_FACTORY_
#define _STORM_PRECONDITIONER_FACTORY_

#include <cstring>
#include <stdexcept>

#include <stormSolvers/stormPreconditioner.hxx>
#include <stormSolvers/stormPreconditionerPolynomial.hxx>

/// Precondtioner types.
/// @{
#define STORM_PRE_NONE   ""
#define STORM_PRE_ID     "ID"
#define STORM_PRE_JACOBI "JACOBI"
#define STORM_PRE_LU_SGS "LU_SGS"
#define STORM_PRE_IC0    "IC0"
#define STORM_PRE_ICT    "ICT"
#define STORM_PRE_ILU0   "ILU0"
#define STORM_PRE_ILUT   "ILUT"
#define STORM_PRE_AINV0  "AINV0"
#define STORM_PRE_AINV   "AINV"
#define STORM_PRE_SPAI0  "SPAI0"
#define STORM_PRE_SPAI   "SPAI"
#define STORM_PRE_CHEBY  "Cheby"
#define STORM_PRE_KRYLOV "Krylov"
/// @}

/// ----------------------------------------------------------------- ///
/// @brief Make preconditioner of the specified type.
/// ----------------------------------------------------------------- ///
template<class tArray>
std::unique_ptr<stormPreconditioner<tArray>>
    stormMakePreconditioner(stormString_t preType) {

  if (std::strcmp(preType, STORM_PRE_NONE) == 0) {

    return nullptr;

  } else if (std::strcmp(preType, STORM_PRE_ID) == 0) {

    return std::make_unique<stormIdentityPreconditioner<tArray>>();

  } else if (std::strcmp(preType, STORM_PRE_JACOBI) == 0) {

    _STORM_NOT_IMPLEMENTED_();

  } else if (std::strcmp(preType, STORM_PRE_LU_SGS) == 0) {

    _STORM_NOT_IMPLEMENTED_();

  } else if (std::strcmp(preType, STORM_PRE_IC0) == 0) {

    _STORM_NOT_IMPLEMENTED_();

  } else if (std::strcmp(preType, STORM_PRE_ICT) == 0) {

    _STORM_NOT_IMPLEMENTED_();

  } else if (std::strcmp(preType, STORM_PRE_ILU0) == 0) {

    _STORM_NOT_IMPLEMENTED_();

  } else if (std::strcmp(preType, STORM_PRE_ILUT) == 0) {

    _STORM_NOT_IMPLEMENTED_();

  } else if (std::strcmp(preType, STORM_PRE_AINV0) == 0) {

    _STORM_NOT_IMPLEMENTED_();

  } else if (std::strcmp(preType, STORM_PRE_AINV) == 0) {

    _STORM_NOT_IMPLEMENTED_();

  } else if (std::strcmp(preType, STORM_PRE_SPAI0) == 0) {

    _STORM_NOT_IMPLEMENTED_();

  } else if (std::strcmp(preType, STORM_PRE_SPAI) == 0) {

    _STORM_NOT_IMPLEMENTED_();

  } else if (std::strcmp(preType, STORM_PRE_CHEBY) == 0) {

    return std::make_unique<stormChebyshevPreconditioner<tArray>>();

  } else if (std::strcmp(preType, STORM_PRE_KRYLOV) == 0) {

    _STORM_NOT_IMPLEMENTED_();

  }

  throw std::invalid_argument("Invalid preconditioner type specified.");

} // stormMakePreconditioner<...>

#endif // _STORM_PRECONDITIONER_FACTORY_
