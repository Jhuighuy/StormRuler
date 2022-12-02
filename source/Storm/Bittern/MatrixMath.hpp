// Copyright (C) 2020-2023 Oleg Butakov
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

#pragma once

#include <Storm/Base.hpp>

#include <Storm/Bittern/Functors.hpp>
#include <Storm/Bittern/Math.hpp>
#include <Storm/Bittern/Matrix.hpp>
#include <Storm/Bittern/MatrixAlgorithms.hpp>
#include <Storm/Bittern/MatrixView.hpp>

#include <concepts>
#include <tuple>
#include <type_traits>
#include <utility>

namespace Storm
{

// -----------------------------------------------------------------------------

/// @brief Element-wise apply function to the matrices.
template<std::copy_constructible Func, matrix_view... Matrices>
  requires std::is_object_v<Func> && (sizeof...(Matrices) >= 1) &&
           std::regular_invocable<Func, matrix_element_t<Matrices>...>
class MapMatrixView final :
    public MatrixViewInterface<MapMatrixView<Func, Matrices...>>
{
private:

  STORM_NO_UNIQUE_ADDRESS_ Func func_;
  STORM_NO_UNIQUE_ADDRESS_ std::tuple<Matrices...> mats_;

public:

  /// @brief Construct a map view.
  constexpr MapMatrixView(Func func, Matrices... mats)
      : func_{std::move(func)}, mats_{std::move(mats)...}
  {
    std::apply(
        [](const auto& first_mat, const auto&... rest_mats) {
          STORM_ASSERT_((first_mat.shape() == rest_mats.shape()) && ...,
                        "Shapes of the matrix arguments are mismatched.");
        },
        mats_);
  }

  /// @copydoc MatrixViewInterface::shape
  constexpr auto shape() const noexcept
  {
    return std::get<0>(mats_).shape();
  }

  /// @copydoc MatrixViewInterface::operator()
  template<class... Indices>
    requires compatible_matrix_indices_v<MapMatrixView, Indices...>
  constexpr auto operator()(Indices... indices) const noexcept
  {
    STORM_ASSERT_(in_range(shape(), indices...), "Indices are out of range!");
    auto compute_element = [this, &indices...](const Matrices&... mats) {
      return func_(mats(indices...)...);
    };
    return std::apply(compute_element, mats_);
  }

}; // class MapMatrixView

template<class Func, class... Matrices>
MapMatrixView(Func, Matrices&&...)
    -> MapMatrixView<Func, forward_as_matrix_view_t<Matrices>...>;

/// @brief Element-wise apply function @p func to the matrices @p mats.
template<class Func, viewable_matrix... Matrices>
  requires compatible_matrices_v<Matrices...> &&
           std::regular_invocable<Func, matrix_element_t<Matrices>...>
constexpr auto map(Func&& func, Matrices&&... mats)
{
  return MapMatrixView{std::forward<Func>(func),
                       std::forward<Matrices>(mats)...};
}

// -----------------------------------------------------------------------------

/// @brief Cast the @p mat elements to another type.
/// @tparam To Type, the elements would be cased to.
template<class To, viewable_matrix Matrix>
  requires std::convertible_to<matrix_element_t<Matrix>, To>
constexpr auto matrix_cast(Matrix&& mat)
{
  return map(Cast<To>{}, std::forward<Matrix>(mat));
}

// -----------------------------------------------------------------------------

/// @brief Logically negate the boolean matrix @p mat.
template<viewable_matrix Matrix>
  requires bool_matrix<Matrix>
constexpr auto operator!(Matrix&& mat)
{
  return map(Not{}, std::forward<Matrix>(mat));
}

/// @brief Element-wise logically "and"
/// the boolean matrices @p mat1 and @p mat2.
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
  requires compatible_matrices_v<Matrix1, Matrix2> && //
           bool_matrix<Matrix1> && bool_matrix<Matrix2>
constexpr auto matrix_and(Matrix1&& mat1, Matrix2&& mat2)
{
  return map(And{}, //
             std::forward<Matrix1>(mat1), std::forward<Matrix2>(mat2));
}

/// @brief Element-wise logically "or"
/// the boolean matrices @p mat1 and @p mat2.
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
  requires compatible_matrices_v<Matrix1, Matrix2> && //
           bool_matrix<Matrix1> && bool_matrix<Matrix2>
constexpr auto matrix_or(Matrix1&& mat1, Matrix2&& mat2)
{
  return map(Or{}, //
             std::forward<Matrix1>(mat1), std::forward<Matrix2>(mat2));
}

/// @brief Element-wise merge the matrices @p then_mat and @p else_mat
/// based on the condition matrix @p cond_mat.
template<viewable_matrix CondMatrix, //
         viewable_matrix ThenMatrix, viewable_matrix ElseMatrix>
  requires compatible_matrices_v<CondMatrix, ThenMatrix, ElseMatrix> &&
           bool_matrix<CondMatrix> &&
           std::common_with<matrix_element_t<ThenMatrix>,
                            matrix_element_t<ElseMatrix>>
constexpr auto merge(CondMatrix&& cond_mat, //
                     ThenMatrix&& then_mat, ElseMatrix&& else_mat)
{
  return map(Merge{}, //
             std::forward<CondMatrix>(cond_mat),
             std::forward<ThenMatrix>(then_mat),
             std::forward<ElseMatrix>(else_mat));
}

// -----------------------------------------------------------------------------

/// @brief Element-wise compare the matrices @p mat1 and @p mat2.
/// @{
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
  requires compatible_matrices_v<Matrix1, Matrix2>
constexpr auto operator==(Matrix1&& mat1, Matrix2&& mat2)
{
  return map(Equal{}, //
             std::forward<Matrix1>(mat1), std::forward<Matrix2>(mat2));
}
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
  requires compatible_matrices_v<Matrix1, Matrix2>
constexpr auto operator!=(Matrix1&& mat1, Matrix2&& mat2)
{
  return map(NotEqual{}, //
             std::forward<Matrix1>(mat1), std::forward<Matrix2>(mat2));
}
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
  requires compatible_matrices_v<Matrix1, Matrix2>
constexpr auto operator<(Matrix1&& mat1, Matrix2&& mat2)
{
  return map(Less{}, //
             std::forward<Matrix1>(mat1), std::forward<Matrix2>(mat2));
}
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
  requires compatible_matrices_v<Matrix1, Matrix2>
constexpr auto operator<=(Matrix1&& mat1, Matrix2&& mat2)
{
  return map(LessEqual{}, //
             std::forward<Matrix1>(mat1), std::forward<Matrix2>(mat2));
}
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
  requires compatible_matrices_v<Matrix1, Matrix2>
constexpr auto operator>(Matrix1&& mat1, Matrix2&& mat2)
{
  return map(Greater{}, //
             std::forward<Matrix1>(mat1), std::forward<Matrix2>(mat2));
}
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
  requires compatible_matrices_v<Matrix1, Matrix2>
constexpr auto operator>=(Matrix1&& mat1, Matrix2&& mat2)
{
  return map(GreaterEqual{}, //
             std::forward<Matrix1>(mat1), std::forward<Matrix2>(mat2));
}
/// @}

/// @brief Element-wise compare the matrices @p mat1 and @p mat2
/// to be approximately equal (within tolerance @p tolerance ).
/// @{
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
  requires compatible_matrices_v<Matrix1, Matrix2> && //
           numeric_matrix<Matrix1> && numeric_matrix<Matrix2>
constexpr auto approx_equal(Matrix1&& mat1, Matrix2&& mat2)
{
  return map(ApproxEqual{}, //
             std::forward<Matrix1>(mat1), std::forward<Matrix2>(mat2));
}
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
  requires compatible_matrices_v<Matrix1, Matrix2> && //
           numeric_matrix<Matrix1> && numeric_matrix<Matrix2>
constexpr auto approx_equal(Matrix1&& mat1, Matrix2&& mat2,
                            long double tolerance)
{
  STORM_ASSERT_(tolerance > 0.0l, "Negative tolerance!");
  return map(ApproxEqual{tolerance}, //
             std::forward<Matrix1>(mat1), std::forward<Matrix2>(mat2));
}
/// @}

/// @brief Element-wise minimum of the matrices.
template<viewable_matrix... Matrices>
  requires compatible_matrices_v<Matrices...>
constexpr auto min(Matrices&&... matrices)
{
  return map(Min{}, std::forward<Matrices>(matrices)...);
}

/// @brief Element-wise maximum of the matrices.
template<viewable_matrix... Matrices>
  requires compatible_matrices_v<Matrices...>
constexpr auto max(Matrices&&... matrices)
{
  return map(Max{}, std::forward<Matrices>(matrices)...);
}

// -----------------------------------------------------------------------------

/// @brief "+" the matrix @p mat.
template<viewable_matrix Matrix>
  requires numeric_matrix<Matrix>
constexpr auto operator+(Matrix&& mat)
{
  // Since operator "+" does nothing on artimetic types,
  // simply forward matrices with artimetic elements as matrix views.
  return forward_as_matrix_view(std::forward<Matrix>(mat));
}

/// @brief Negate the matrix @p mat.
template<viewable_matrix Matrix>
  requires numeric_matrix<Matrix>
constexpr auto operator-(Matrix&& mat)
{
  return map(Negate{}, std::forward<Matrix>(mat));
}

/// @brief Multiply the matrix @p mat by a scalar @p scal.
/// @{
template<std::copyable Scalar, viewable_matrix Matrix>
  requires numeric_type<Scalar> && numeric_matrix<Matrix>
constexpr auto operator*(Scalar scal, Matrix&& mat)
{
  return map(BindFirst{Multiply{}, std::move(scal)}, //
             std::forward<Matrix>(mat));
}
template<viewable_matrix Matrix, std::copyable Scalar>
  requires numeric_matrix<Matrix> && numeric_type<Scalar>
constexpr auto operator*(Matrix&& mat, Scalar scal)
{
  return map(BindLast{Multiply{}, std::move(scal)}, //
             std::forward<Matrix>(mat));
}
/// @}
/// @brief Multiply-assign the matrix @p mat by a scalar @p scal.
template<output_matrix Matrix, std::copyable Scalar>
  requires numeric_matrix<Matrix> && numeric_type<Scalar>
constexpr decltype(auto) operator*=(Matrix&& mat, Scalar scal)
{
  return assign(BindLast{MultiplyAssign{}, std::move(scal)},
                std::forward<Matrix>(mat));
}

/// @brief Divide the scalar @p scal by a matrix @p mat.
template<std::copyable Scalar, viewable_matrix Matrix>
  requires numeric_type<Scalar> && numeric_matrix<Matrix>
constexpr auto operator/(Scalar scal, Matrix&& mat)
{
  return map(BindFirst{Divide{}, std::move(scal)}, //
             std::forward<Matrix>(mat));
}
/// @brief Divide the matrix @p mat by a scalar @p scal.
template<viewable_matrix Matrix, std::copyable Scalar>
  requires numeric_matrix<Matrix> && numeric_type<Scalar>
constexpr auto operator/(Matrix&& mat, Scalar scal)
{
  return map(BindLast{Divide{}, std::move(scal)}, //
             std::forward<Matrix>(mat));
}
/// @brief Divide-assign the matrix @p mat by a scalar @p scal.
template<output_matrix Matrix, std::copyable Scalar>
  requires numeric_matrix<Matrix> && numeric_type<Scalar>
constexpr decltype(auto) operator/=(Matrix&& mat, Scalar scal)
{
  return assign(BindLast{DivideAssign{}, std::move(scal)}, //
                std::forward<Matrix>(mat));
}

/// @brief Add the matrices @p mat1 and @p mat2.
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
  requires compatible_matrices_v<Matrix1, Matrix2> && //
           numeric_matrix<Matrix1> && numeric_matrix<Matrix2>
constexpr auto operator+(Matrix1&& mat1, Matrix2&& mat2)
{
  return map(Add{}, //
             std::forward<Matrix1>(mat1), std::forward<Matrix2>(mat2));
}
/// @brief Add-assign the matrices @p mat1 and @p mat2.
template<output_matrix Matrix1, matrix Matrix2>
  requires compatible_matrices_v<Matrix1, Matrix2> && //
           numeric_matrix<Matrix1> && numeric_matrix<Matrix2>
constexpr decltype(auto) operator+=(Matrix1&& mat1, Matrix2&& mat2)
{
  return assign(AddAssign{}, //
                std::forward<Matrix1>(mat1), std::forward<Matrix2>(mat2));
}

/// @brief Subtract the matrices @p mat1 and @p mat2.
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
  requires compatible_matrices_v<Matrix1, Matrix2> && //
           numeric_matrix<Matrix1> && numeric_matrix<Matrix2>
constexpr auto operator-(Matrix1&& mat1, Matrix2&& mat2)
{
  return map(Subtract{}, //
             std::forward<Matrix1>(mat1), std::forward<Matrix2>(mat2));
}
/// @brief Subtract-assign the matrices @p mat1 and @p mat2.
template<output_matrix Matrix1, matrix Matrix2>
  requires compatible_matrices_v<Matrix1, Matrix2> && //
           numeric_matrix<Matrix1> && numeric_matrix<Matrix2>
constexpr decltype(auto) operator-=(Matrix1&& mat1, Matrix2&& mat2)
{
  return assign(SubtractAssign{}, //
                std::forward<Matrix1>(mat1), std::forward<Matrix2>(mat2));
}

/// @brief Element-wise multiply the matrices @p mat1 and @p mat2.
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
  requires compatible_matrices_v<Matrix1, Matrix2> && //
           numeric_matrix<Matrix1> && numeric_matrix<Matrix2>
constexpr auto operator*(Matrix1&& mat1, Matrix2&& mat2)
{
  return map(Multiply{}, //
             std::forward<Matrix1>(mat1), std::forward<Matrix2>(mat2));
}
/// @brief Element-wise multiply-assign the matrices @p mat1 and @p mat2.
template<output_matrix Matrix1, matrix Matrix2>
  requires compatible_matrices_v<Matrix1, Matrix2> && //
           numeric_matrix<Matrix1> && numeric_matrix<Matrix2>
constexpr decltype(auto) operator*=(Matrix1&& mat1, Matrix2&& mat2)
{
  return assign(MultiplyAssign{}, //
                std::forward<Matrix1>(mat1), std::forward<Matrix2>(mat2));
}

/// @brief Element-wise divide the matrices @p mat1 and @p mat2.
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
  requires compatible_matrices_v<Matrix1, Matrix2> && //
           numeric_matrix<Matrix1> && numeric_matrix<Matrix2>
constexpr auto operator/(Matrix1&& mat1, Matrix2&& mat2)
{
  return map(Divide{}, //
             std::forward<Matrix1>(mat1), std::forward<Matrix2>(mat2));
}
/// @brief Element-wise divide-assign the matrices @p mat1 and @p mat2.
template<output_matrix Matrix1, matrix Matrix2>
  requires compatible_matrices_v<Matrix1, Matrix2> && //
           numeric_matrix<Matrix1> && numeric_matrix<Matrix2>
constexpr decltype(auto) operator/=(Matrix1&& mat1, Matrix2&& mat2)
{
  return assign(DivideAssign{}, //
                std::forward<Matrix1>(mat1), std::forward<Matrix2>(mat2));
}

/// @brief Normalize the matrix @p mat (divide by it's norm).
template<viewable_matrix Matrix>
  requires numeric_matrix<Matrix>
constexpr auto normalize(Matrix&& mat)
{
  auto mat_view = forward_as_matrix_view(std::forward<Matrix>(mat));
  const auto mat_norm = norm_2(mat_view);
  return safe_divide(1.0, mat_norm) * std::move(mat_view);
}

// -----------------------------------------------------------------------------

#define MAKE_UNARY_MATRIX_FUNC_(func)                   \
  template<viewable_matrix Matrix>                      \
    requires numeric_matrix<Matrix>                     \
  constexpr auto func(Matrix&& mat)                     \
  {                                                     \
    auto func_object = []<class Arg>(Arg&& arg) {       \
      return func(std::forward<Arg>(arg));              \
    };                                                  \
    return map(func_object, std::forward<Matrix>(mat)); \
  }

#define MAKE_BINARY_MATRIX_FUNC_(func)                                        \
  template<std::copyable Scalar, viewable_matrix Matrix>                      \
    requires numeric_type<Scalar> && numeric_matrix<Matrix>                   \
  constexpr auto func(Scalar scal, Matrix&& mat)                              \
  {                                                                           \
    auto func_object = []<class Arg1, class Arg2>(Arg1&& arg1, Arg2&& arg2) { \
      return func(std::forward<Arg1>(arg1), std::forward<Arg2>(arg2));        \
    };                                                                        \
    return map(BindFirst{func_object, scal}, std::forward<Matrix>(mat));      \
  }                                                                           \
  template<viewable_matrix Matrix, std::copyable Scalar>                      \
    requires numeric_matrix<Matrix> && numeric_type<Scalar>                   \
  constexpr auto func(Matrix&& mat, Scalar scal)                              \
  {                                                                           \
    auto func_object = []<class Arg1, class Arg2>(Arg1&& arg1, Arg2&& arg2) { \
      return func(std::forward<Arg1>(arg1), std::forward<Arg2>(arg2));        \
    };                                                                        \
    return map(BindLast{func_object, scal}, std::forward<Matrix>(mat));       \
  }                                                                           \
  template<viewable_matrix Matrix1, viewable_matrix Matrix2>                  \
    requires numeric_matrix<Matrix1> && numeric_matrix<Matrix2>               \
  constexpr auto func(Matrix1&& mat1, Matrix2&& mat2)                         \
  {                                                                           \
    auto func_object = []<class Arg1, class Arg2>(Arg1&& arg1, Arg2&& arg2) { \
      return func(std::forward<Arg1>(arg1), std::forward<Arg2>(arg2));        \
    };                                                                        \
    return map(func_object, /**/                                              \
               std::forward<Matrix1>(mat1), std::forward<Matrix2>(mat2));     \
  }

/// @brief Element-wise absolute value of the matrix.
MAKE_UNARY_MATRIX_FUNC_(abs)
/// @brief Element-wise sign of the matrix.
MAKE_UNARY_MATRIX_FUNC_(sign)

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

// -----------------------------------------------------------------------------

} // namespace Storm
