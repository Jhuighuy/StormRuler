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

  /// @brief Compute a dot product of @p x_vec and @p y_vec.
  static auto Dot(Vector const& x_vec, Vector const& y_vec) = delete;

  /// @brief Compute a norm of @p x_vec.
  static auto Norm2(Vector const& x_vec) = delete;

  /// @brief Swap @p x_vec and @p y_vec.
  static void Swap(Vector& x_vec, Vector& y_vec) = delete;

  /// @brief Compute @p x_vec = @p y_vec.
  static void Set(Vector& x_vec, Vector const& y_vec) = delete;

  /// @brief Randomly fill the @p x_vec with value @p a.
  static void Fill(Vector& x_vec, auto a) = delete;

  /// @brief Randomly fill the @p x_vec.
  static void RandFill(Vector& x_vec) = delete;

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
    { targetVec.assign(sourceVector) };
    { targetVec.assign(sourceVector, copyContents) };
  } &&
  /// @pre Require the dot product operation. 
  requires {
    typename DotType<Vector>;
  } &&
  requires(Vector const& x_vec) {
    { VectorOperations<Vector>::Norm2(x_vec) } -> std::convertible_to<real_t>;
  } &&
  /// @pre Require the fill and random fill operations. 
  requires(Vector& x_vec, DotType<Vector> a) {
    VectorOperations<Vector>::Fill(x_vec, a);
    VectorOperations<Vector>::RandFill(x_vec);
  };

// clang-format on

namespace Blas {

/// @brief Compute a dot product of @p x_vec and @p y_vec.
template<VectorLike Vector>
auto Dot(Vector const& x_vec, Vector const& y_vec) {
  return VectorOperations<Vector>::Dot(x_vec, y_vec);
}

/// @brief Compute a norm of @p x_vec.
template<VectorLike Vector>
real_t Norm2(Vector const& x_vec) {
  return VectorOperations<Vector>::Norm2(x_vec);
}

/// @brief Randomly fill the @p x_vec with value @p a.
template<VectorLike Vector>
void Fill(Vector& x_vec, auto a) {
  VectorOperations<Vector>::Fill(x_vec, a);
}

/// @brief Randomly fill the @p x_vec.
template<VectorLike Vector>
void RandFill(Vector& x_vec) {
  VectorOperations<Vector>::RandFill(x_vec);
}

} // namespace Blas

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Operator-like concept.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Operator, class InVector, class OutVector = InVector>
concept operator_like = requires(Operator& any_op, OutVector& y_vec,
                                 InVector const& x_vec) {
  any_op(y_vec, x_vec);
};

} // namespace Storm
