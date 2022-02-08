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

#include <string_view>
#include <stdexcept>

#include <stormBase.hxx>

#include <stormSolvers/stormPreconditioner.hxx>
#include <stormSolvers/stormPreconditionerPolynomial.hxx>

_STORM_NAMESPACE_BEGIN_

/// ----------------------------------------------------------------- ///
/// @brief Precondtioner types.
/// ----------------------------------------------------------------- ///
namespace PreconditionerType {

  /// @brief No preconditioning.
  static std::string_view const None = "";

  /// @brief Identity preconditioner.
  static std::string_view const Identity = "Identity";

  /// @brief @c Jacobi preconditioner.
  static std::string_view const Jacobi = "Jacobi";

  /// @brief @c SGS preconditioner.
  static std::string_view const Sgs = "Sgs";

  /// @brief @c IC(0) preconditioner.
  static std::string_view const Ic0 = "Ic0";

  /// @brief @c IC(t) preconditioner.
  static std::string_view const Ict = "Ict";

  /// @brief @c ILU(0) preconditioner.
  static std::string_view const Ilu0 = "Ilu0";

  /// @brief @c ILU(t) preconditioner.
  static std::string_view const Ilut = "Ilut";

  /// @brief @c ILQ(0) preconditioner.
  static std::string_view const Ilq0 = "Ilq0";

  /// @brief @c ILQ(t) preconditioner.
  static std::string_view const Ilqt = "Ilqt";

  /// @brief @c AINV(0) preconditioner.
  static std::string_view const Ainv0 = "Ainv0";

  /// @brief @c AINV preconditioner.
  static std::string_view const Ainv = "Ainv";

  /// @brief @c SPAI(0) preconditioner.
  static std::string_view const Spai0 = "Spai0";

  /// @brief @c SPAI preconditioner.
  static std::string_view const Spai = "Spai";

  /// @brief @c Chebyshev polynomial preconditioner.
  static std::string_view const Chebyshev = "Chebyshev";

  /// @brief @c Krylov preconditioner.
  static std::string_view const Krylov = "Krylov";

} // namespace PreconditionerType

/// ----------------------------------------------------------------- ///
/// @brief Make preconditioner of the specified type.
/// ----------------------------------------------------------------- ///
template<class Vector>
std::unique_ptr<Preconditioner<Vector>>
    MakePreconditioner(std::string_view const& preType = {}) {

  if (preType.empty() || preType == PreconditionerType::None) {
    return nullptr;
  }
  if (preType == PreconditionerType::Identity) {
    return std::make_unique<IdentityPreconditioner<Vector>>();
  }
  if (preType == PreconditionerType::Jacobi) {
    //return std::make_unique<JacobiPreconditioner<Vector>>();
  }
  if (preType == PreconditionerType::Sgs) {
    //return std::make_unique<SgsPreconditioner<Vector>>();
  }
  if (preType == PreconditionerType::Ic0) {
    //return std::make_unique<Ic0Preconditioner<Vector>>();
  }
  if (preType == PreconditionerType::Ict) {
    //return std::make_unique<IctPreconditioner<Vector>>();
  }
  if (preType == PreconditionerType::Ilu0) {
    //return std::make_unique<Ilu0Preconditioner<Vector>>();
  }
  if (preType == PreconditionerType::Ilut) {
    //return std::make_unique<IlutPreconditioner<Vector>>();
  }
  if (preType == PreconditionerType::Ilq0) {
    //return std::make_unique<Ilq0Preconditioner<Vector>>();
  }
  if (preType == PreconditionerType::Ilqt) {
    //return std::make_unique<IlqtPreconditioner<Vector>>();
  }
  if (preType == PreconditionerType::Ainv0) {
    //return std::make_unique<Ainv0Preconditioner<Vector>>();
  }
  if (preType == PreconditionerType::Ainv) {
    //return std::make_unique<AinvPreconditioner<Vector>>();
  }
  if (preType == PreconditionerType::Spai0) {
    //return std::make_unique<Spai0Preconditioner<Vector>>();
  }
  if (preType == PreconditionerType::Spai) {
    //return std::make_unique<SpaiPreconditioner<Vector>>();
  }
  if (preType == PreconditionerType::Chebyshev) {
    return std::make_unique<ChebyshevPreconditioner<Vector>>();
  }
  if (preType == PreconditionerType::Krylov) {
    //return std::make_unique<KrylovPreconditioner<Vector>>();
  }

  throw std::invalid_argument("Invalid preconditioner type specified.");

} // MakePreconditioner<...>

_STORM_NAMESPACE_END_
