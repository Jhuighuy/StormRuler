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

#include <stormBase.hxx>

namespace Storm {

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Set of the operations for given vector type.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Vector>
class VectorOperations {
public:

  /// @brief Compute a dot product of @p xVec and @p yVec.
  static auto Dot(Vector const& xVec, Vector const& yVec) = delete;

  /// @brief Compute a norm of @p xVec.
  static auto Norm2(Vector const& xVec) = delete;

  /// @brief Swap @p xVec and @p yVec.
  static void Swap(Vector& xVec, Vector& yVec) = delete;

  /// @brief Compute @p xVec = @p yVec.
  static void Set(Vector& xVec, Vector const& yVec) = delete;

  /// @brief Randomly fill the @p xVec with value @p a.
  static void Fill(Vector& xVec, auto a) = delete;

  /// @brief Randomly fill the @p xVec.
  static void RandFill(Vector& xVec) = delete;

  /// @brief Compute @p xVec *= @p a.
  static void ScaleAssign(Vector& xVec, auto a);

}; // class VectorOperations

/// @brief A type of a dot product of the two vectors.
template<class Vector>
using DotType = decltype(VectorOperations<Vector>::Dot(std::declval<Vector>(),
                                                       std::declval<Vector>()));

// clang-format off

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Vector concept that supports the level 1 BLAS operations.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Vector>
concept VectorLike = 
  /// @pre Require the vector type to be assignable.
  requires(Vector& targetVec, Vector const& sourceVector, 
           bool copyContents) {
    { targetVec.Assign(sourceVector) };
    { targetVec.Assign(sourceVector, copyContents) };
  } &&
  /// @pre Require the dot product operation. 
  requires {
    typename DotType<Vector>;
  } &&
  requires(Vector const& xVec) {
    { VectorOperations<Vector>::Norm2(xVec) } -> std::convertible_to<real_t>;
  } &&
  /// @pre Require the fill and random fill operations. 
  requires(Vector& xVec, DotType<Vector> a) {
    VectorOperations<Vector>::Fill(xVec, a);
    VectorOperations<Vector>::RandFill(xVec);
  } &&
  /// @pre Require the addition operation. 
  requires(Vector& xVec, Vector const& yVec, DotType<Vector> a) {
    VectorOperations<Vector>::ScaleAssign(xVec, a);
  };

// clang-format on

namespace Blas {

/// @brief Compute a dot product of @p xVec and @p yVec.
template<VectorLike Vector>
auto Dot(Vector const& xVec, Vector const& yVec) {
  return VectorOperations<Vector>::Dot(xVec, yVec);
}

/// @brief Compute a norm of @p xVec.
template<VectorLike Vector>
real_t Norm2(Vector const& xVec) {
  return VectorOperations<Vector>::Norm2(xVec);
}

/// @brief Randomly fill the @p xVec with value @p a.
template<VectorLike Vector>
void Fill(Vector& xVec, auto a) {
  VectorOperations<Vector>::Fill(xVec, a);
}

/// @brief Randomly fill the @p xVec.
template<VectorLike Vector>
void RandFill(Vector& xVec) {
  VectorOperations<Vector>::RandFill(xVec);
}

/// @brief Compute @p xVec *= @p a.
template<VectorLike Vector>
void ScaleAssign(Vector& xVec, auto a) {
  VectorOperations<Vector>::ScaleAssign(xVec, a);
}

/// @brief Compute @p xVec = @p a * @p yVec.
template<VectorLike Vector>
void Scale(Vector& xVec, Vector const& yVec, auto a) {
  VectorOperations<Vector>::Scale(xVec, yVec, a);
}

} // namespace Blas

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Operator-like concept.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Operator, class InVector, class OutVector = InVector>
concept operator_like = requires(Operator& anyOp, OutVector& yVec,
                                 InVector const& xVec) {
  anyOp(yVec, xVec);
};

} // namespace Storm
