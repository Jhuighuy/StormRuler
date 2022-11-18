// Copyright (C) 2020 - 2023 Oleg Butakov
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR Allocator PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT
// SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
// OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

// Automatically included with `Matrix.hpp`.
#pragma once

#include <Storm/Base.hpp>

#include <Storm/Utils/Math.hpp>

#include <Storm/Bittern/Matrix.hpp>

#include <concepts>
#include <functional>
#include <tuple>
#include <type_traits>
#include <utility>

namespace Storm {

/// @brief Element-wise product of function to matrices view.
template<std::copy_constructible Func, matrix_view... Matrices>
  requires std::is_object_v<Func> && (sizeof...(Matrices) >= 1) &&
           std::regular_invocable<Func, matrix_element_decltype_t<Matrices>...>
class MapMatrixView final :
    public MatrixViewInterface<MapMatrixView<Func, Matrices...>> {
private:

  STORM_NO_UNIQUE_ADDRESS_ Func func_;
  STORM_NO_UNIQUE_ADDRESS_ std::tuple<Matrices...> mats_;

public:

  /// @brief Construct a map view.
  constexpr MapMatrixView(Func func, Matrices... mats)
      : func_{std::move(func)}, mats_{std::move(mats)...} {
#if !STORM_COMPILER_MSVC_
    STORM_ASSERT_( //
        std::apply(
            [](const auto& first_mat, const auto&... rest_mats) {
              return ((first_mat.shape() == rest_mats.shape()) && ...);
            },
            mats_),
        "Shapes of the matrix arguments are mismatched.");
#endif
  }

  /// @copydoc MatrixViewInterface::shape
  [[nodiscard]] constexpr auto shape() const noexcept {
    return std::get<0>(mats_).shape();
  }

  /// @copydoc MatrixViewInterface::operator()
  [[nodiscard]] constexpr decltype(auto)
  operator()(size_t row_index, size_t col_index) const noexcept {
    STORM_ASSERT_(shape().in_range(row_index, col_index),
                  "Indices are out of range!");
    return std::apply(
        [&](const Matrices&... mats) {
          return func_(mats(row_index, col_index)...);
        },
        mats_);
  }

}; // class MapMatrixView

template<class Func, class... Matrices>
MapMatrixView(Func, Matrices&&...)
    -> MapMatrixView<Func, forward_as_matrix_view_t<Matrices>...>;

/// @brief Make a element-wise product of function @p func
///   to matrices @p mats.
template<class Func, viewable_matrix... Matrices>
  requires std::regular_invocable<Func, matrix_element_decltype_t<Matrices>...>
[[nodiscard]] constexpr auto map(Func&& func, Matrices&&... mats) {
  return MapMatrixView(std::forward<Func>(func),
                       std::forward<Matrices>(mats)...);
}

// -----------------------------------------------------------------------------

/// @brief Logically negate the boolean matrix @p mat.
template<viewable_matrix Matrix>
  requires bool_matrix<Matrix>
[[nodiscard]] constexpr auto operator!(Matrix&& mat) {
  return MapMatrixView(std::logical_not{}, std::forward<Matrix>(mat));
}

/// @brief Element-wise logically "and"
/// the boolean matrices @p mat1 and @p mat2.
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
  requires bool_matrix<Matrix1> && bool_matrix<Matrix2>
[[nodiscard]] constexpr auto operator&&(Matrix1&& mat1, Matrix2&& mat2) {
  return MapMatrixView(std::logical_and{}, //
                       std::forward<Matrix1>(mat1),
                       std::forward<Matrix2>(mat2));
}

/// @brief Element-wise logically "or"
/// the boolean matrices @p mat1 and @p mat2.
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
  requires bool_matrix<Matrix1> && bool_matrix<Matrix2>
[[nodiscard]] constexpr auto operator||(Matrix1&& mat1, Matrix2&& mat2) {
  return MapMatrixView(std::logical_or{}, //
                       std::forward<Matrix1>(mat1),
                       std::forward<Matrix2>(mat2));
}

// -----------------------------------------------------------------------------

/// @brief Element-wise compare the matrices @p mat1 and @p mat2.
/// @{
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
[[nodiscard]] constexpr auto operator==(Matrix1&& mat1, Matrix2&& mat2) {
  return MapMatrixView(std::equal_to{}, //
                       std::forward<Matrix1>(mat1),
                       std::forward<Matrix2>(mat2));
}
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
[[nodiscard]] constexpr auto operator!=(Matrix1&& mat1, Matrix2&& mat2) {
  return MapMatrixView(std::not_equal_to{}, //
                       std::forward<Matrix1>(mat1),
                       std::forward<Matrix2>(mat2));
}
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
[[nodiscard]] constexpr auto operator<(Matrix1&& mat1, Matrix2&& mat2) {
  return MapMatrixView(std::less{}, //
                       std::forward<Matrix1>(mat1),
                       std::forward<Matrix2>(mat2));
}
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
[[nodiscard]] constexpr auto operator<=(Matrix1&& mat1, Matrix2&& mat2) {
  return MapMatrixView(std::less_equal{}, //
                       std::forward<Matrix1>(mat1),
                       std::forward<Matrix2>(mat2));
}
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
[[nodiscard]] constexpr auto operator>(Matrix1&& mat1, Matrix2&& mat2) {
  return MapMatrixView(std::greater{}, //
                       std::forward<Matrix1>(mat1),
                       std::forward<Matrix2>(mat2));
}
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
[[nodiscard]] constexpr auto operator>=(Matrix1&& mat1, Matrix2&& mat2) {
  return MapMatrixView(std::greater_equal{}, //
                       std::forward<Matrix1>(mat1),
                       std::forward<Matrix2>(mat2));
}
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
[[nodiscard]] constexpr auto operator<=>(Matrix1&& mat1, Matrix2&& mat2) {
  return MapMatrixView(std::compare_three_way{}, //
                       std::forward<Matrix1>(mat1),
                       std::forward<Matrix2>(mat2));
}
/// @}

/// @brief Element-wise compare the matrices @p mat1 and @p mat2
/// to be approximately equal.
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
  requires numeric_matrix<Matrix1> && numeric_matrix<Matrix2>
[[nodiscard]] constexpr auto approx_eq(Matrix1&& mat1, //
                                       Matrix2&& mat2, real_t tolerance) {
  STORM_ASSERT_(tolerance > 0.0, "Negative comparison tolerance!");
  return MapMatrixView(
      [=]<class Elem1, class Elem2>(Elem1&& elem1, Elem2&& elem2) {
        return abs(std::forward<Elem1>(elem1) - //
                   std::forward<Elem2>(elem2)) <= tolerance;
      },
      std::forward<Matrix1>(mat1), std::forward<Matrix2>(mat2));
}

// -----------------------------------------------------------------------------

/// @brief "+" the matrix @p mat.
template<viewable_matrix Matrix>
  requires numeric_matrix<Matrix>
[[nodiscard]] constexpr auto operator+(Matrix&& mat) {
  // Since operator "+" does nothing on artimetic types,
  // simply forward matrices with artimetic elements as matrix views.
  return forward_as_matrix_view(std::forward<Matrix>(mat));
}

/// @brief Negate the matrix @p mat.
template<viewable_matrix Matrix>
  requires numeric_matrix<Matrix>
[[nodiscard]] constexpr auto operator-(Matrix&& mat) {
  return MapMatrixView(std::negate{}, std::forward<Matrix>(mat));
}

/// @brief Multiply the matrix @p mat by a scalar @p scal.
/// @{
template<viewable_matrix Matrix, std::copyable Scalar>
  requires (numeric_matrix<Matrix> && !matrix<Scalar>)
[[nodiscard]] constexpr auto operator*(Matrix&& mat, Scalar scal) {
  return MapMatrixView(
      [scal = std::move(scal)]<class Elem>(Elem&& elem) noexcept {
        return std::forward<Elem>(elem) * scal;
      },
      std::forward<Matrix>(mat));
}
template<std::copyable Scalar, viewable_matrix Matrix>
  requires (!matrix<Scalar> && numeric_matrix<Matrix>)
[[nodiscard]] constexpr auto operator*(Scalar scal, Matrix&& mat) {
  return MapMatrixView(
      [scal = std::move(scal)]<class Elem>(Elem&& elem) noexcept {
        return scal * std::forward<Elem>(elem);
      },
      std::forward<Matrix>(mat));
}
/// @}

/// @brief Divide the matrix @p mat by a scalar @p scal.
template<viewable_matrix Matrix, std::copyable Scalar>
  requires (numeric_matrix<Matrix> && !matrix<Scalar>)
[[nodiscard]] constexpr auto operator/(Matrix&& mat, Scalar scal) {
  return MapMatrixView(
      [scal = std::move(scal)]<class Elem>(Elem&& elem) noexcept {
        return std::forward<Elem>(elem) / scal;
      },
      std::forward<Matrix>(mat));
}

/// @brief Divide the scalar @p scal by a matrix @p mat.
template<std::copyable Scalar, viewable_matrix Matrix>
  requires (!matrix<Scalar> && numeric_matrix<Matrix>)
[[nodiscard]] constexpr auto operator/(Scalar scal, Matrix&& mat) {
  return MapMatrixView(
      [scal = std::move(scal)]<class Elem>(Elem&& elem) noexcept {
        return scal / std::forward<Elem>(elem);
      },
      std::forward<Matrix>(mat));
}

/// @brief Add the matrices @p mat1 and @p mat2.
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
  requires numeric_matrix<Matrix1> && numeric_matrix<Matrix2>
[[nodiscard]] constexpr auto operator+(Matrix1&& mat1, Matrix2&& mat2) {
  return MapMatrixView(std::plus{}, //
                       std::forward<Matrix1>(mat1),
                       std::forward<Matrix2>(mat2));
}

/// @brief Subtract the matrices @p mat1 and @p mat2.
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
  requires numeric_matrix<Matrix1> && numeric_matrix<Matrix2>
[[nodiscard]] constexpr auto operator-(Matrix1&& mat1, Matrix2&& mat2) {
  return MapMatrixView(std::minus{}, //
                       std::forward<Matrix1>(mat1),
                       std::forward<Matrix2>(mat2));
}

/// @brief Element-wise multiply the matrices @p mat1 and @p mat2.
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
  requires numeric_matrix<Matrix1> && numeric_matrix<Matrix2>
[[nodiscard]] constexpr auto operator*(Matrix1&& mat1, Matrix2&& mat2) {
  return MapMatrixView(std::multiplies{}, //
                       std::forward<Matrix1>(mat1),
                       std::forward<Matrix2>(mat2));
}

/// @brief Element-wise divide the matrices @p mat1 and @p mat2.
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
  requires numeric_matrix<Matrix1> && numeric_matrix<Matrix2>
[[nodiscard]] constexpr auto operator/(Matrix1&& mat1, Matrix2&& mat2) {
  return MapMatrixView(std::divides{}, //
                       std::forward<Matrix1>(mat1),
                       std::forward<Matrix2>(mat2));
}

/// @brief Cast the @p mat elements to another type.
/// @tparam To Type, the elements would be cased to.
template<class To, viewable_matrix Matrix>
  requires std::convertible_to<matrix_element_decltype_t<Matrix>, To>
[[nodiscard]] constexpr auto matrix_cast(Matrix&& mat) {
  return MapMatrixView(
      []<class Elem>(Elem&& elem) noexcept {
        return static_cast<To>(std::forward<Elem>(elem));
      },
      std::forward<Matrix>(mat));
}

// -----------------------------------------------------------------------------

#define MAKE_UNARY_MATRIX_FUNC_(func)               \
  template<viewable_matrix Matrix>                  \
    requires numeric_matrix<Matrix>                 \
  [[nodiscard]] constexpr auto func(Matrix&& mat) { \
    return MapMatrixView(                           \
        []<class Elem>(Elem&& elem) noexcept {      \
          return func(std::forward<Elem>(elem));    \
        },                                          \
        std::forward<Matrix>(mat));                 \
  }

#define MAKE_BINARY_MATRIX_FUNC_(func)                                         \
  template<viewable_matrix Matrix, std::copyable Scalar>                       \
    requires (numeric_matrix<Matrix> && !matrix<Scalar>)                       \
  [[nodiscard]] constexpr auto func(Matrix&& mat, Scalar scal) {               \
    return MapMatrixView(                                                      \
        [scal = std::move(scal)]<class Elem>(Elem&& elem) noexcept {           \
          return func(std::forward<Elem>(elem), scal);                         \
        },                                                                     \
        std::forward<Matrix>(mat));                                            \
  }                                                                            \
  template<std::copyable Scalar, viewable_matrix Matrix>                       \
    requires (!matrix<Scalar> && numeric_matrix<Matrix>)                       \
  [[nodiscard]] constexpr auto func(Scalar scal, Matrix&& mat) {               \
    return MapMatrixView(                                                      \
        [scal = std::move(scal)]<class Elem>(Elem&& elem) noexcept {           \
          return func(scal, std::forward<Elem>(elem));                         \
        },                                                                     \
        std::forward<Matrix>(mat));                                            \
  }                                                                            \
  template<viewable_matrix Matrix1, viewable_matrix Matrix2>                   \
    requires numeric_matrix<Matrix1> && numeric_matrix<Matrix2>                \
  [[nodiscard]] constexpr auto func(Matrix1&& mat1, Matrix2&& mat2) {          \
    return MapMatrixView(                                                      \
        []<class Elem1, class Elem2>(Elem1&& elem1, Elem2&& elem2) noexcept {  \
          return func(std::forward<Elem1>(elem1), std::forward<Elem2>(elem2)); \
        },                                                                     \
        std::forward<Matrix1>(mat1), std::forward<Matrix2>(mat2));             \
  }

/// @brief Element-wise absolute value of the matrix.
MAKE_UNARY_MATRIX_FUNC_(abs)
/// @brief Element-wise sign of the matrix.
MAKE_UNARY_MATRIX_FUNC_(sign)
/// @brief Element-wise minimum of the matrices.
MAKE_BINARY_MATRIX_FUNC_(min)
/// @brief Element-wise maximum of the matrices.
MAKE_BINARY_MATRIX_FUNC_(max)

MAKE_BINARY_MATRIX_FUNC_(pow)
/// @brief Element-wise square root.
MAKE_UNARY_MATRIX_FUNC_(sqrt)
/// @brief Element-wise cube root.
MAKE_UNARY_MATRIX_FUNC_(cbrt)

/// @brief Element-wise exponent.
MAKE_UNARY_MATRIX_FUNC_(exp)
/// @brief Element-wise exponent base 2.
MAKE_UNARY_MATRIX_FUNC_(exp2)
/// @brief Element-wise logarithm.
MAKE_UNARY_MATRIX_FUNC_(log)
/// @brief Element-wise logarithm base 2.
MAKE_UNARY_MATRIX_FUNC_(log2)
/// @brief Element-wise logarithm base 10.
MAKE_UNARY_MATRIX_FUNC_(log10)

/// @brief Element-wise sine.
MAKE_UNARY_MATRIX_FUNC_(sin)
/// @brief Element-wise cosine.
MAKE_UNARY_MATRIX_FUNC_(cos)
/// @brief Element-wise tangent.
MAKE_UNARY_MATRIX_FUNC_(tan)
/// @brief Element-wise inverse sine.
MAKE_UNARY_MATRIX_FUNC_(asin)
/// @brief Element-wise inverse cosine.
MAKE_UNARY_MATRIX_FUNC_(acos)
/// @brief Element-wise inverse tangent.
MAKE_UNARY_MATRIX_FUNC_(atan)
/// @brief Element-wise inverse tangent.
MAKE_BINARY_MATRIX_FUNC_(atan2)

/// @brief Element-wise hyperbolic sine.
MAKE_UNARY_MATRIX_FUNC_(sinh)
/// @brief Element-wise hyperbolic cosine.
MAKE_UNARY_MATRIX_FUNC_(cosh)
/// @brief Element-wise hyperbolic tangent.
MAKE_UNARY_MATRIX_FUNC_(tanh)
/// @brief Element-wise hyperbolic inverse sine.
MAKE_UNARY_MATRIX_FUNC_(asinh)
/// @brief Element-wise hyperbolic inverse cosine.
MAKE_UNARY_MATRIX_FUNC_(acosh)
/// @brief Element-wise hyperbolic inverse tangent.
MAKE_UNARY_MATRIX_FUNC_(atanh)

#undef MAKE_UNARY_MATRIX_FUNC_
#undef MAKE_BINARY_MATRIX_FUNC_

} // namespace Storm
