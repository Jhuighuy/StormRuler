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

#include <Storm/Bittern/Functions.hpp>
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

  static_assert(std::copyable<Func>, "Boxing is not implemented yet!");

  STORM_NO_UNIQUE_ADDRESS_ Func func_;
  STORM_NO_UNIQUE_ADDRESS_ std::tuple<Matrices...> mats_;

public:

  /// @brief Construct a map view.
  constexpr explicit MapMatrixView(Func func, Matrices... mats)
      : func_{std::move(func)}, mats_{std::move(mats)...}
  {
    std::apply(
        [](const auto& first_mat, const auto&... rest_mats) {
          STORM_ASSERT_((first_mat.shape() == rest_mats.shape()) && ...,
                        "Shapes of the matrix arguments are mismatched.");
        },
        mats_);
  }

  /// @brief Get the matrix shape.
  constexpr auto shape() const noexcept
  {
    return std::get<0>(mats_).shape();
  }

  /// @brief Get the matrix element at @p indices.
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
    -> MapMatrixView<Func, matrix_view_t<Matrices>...>;

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
constexpr auto cast_matrix(Matrix&& mat)
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

/// @brief Element-wise logically AND the boolean matrices @p mats.
template<viewable_matrix... Matrices>
  requires compatible_matrices_v<Matrices...> && (... && bool_matrix<Matrices>)
constexpr auto matrix_and(Matrices&&... mats)
{
  return map(And{}, std::forward<Matrices>(mats)...);
}

/// @brief Element-wise logically OR the boolean matrices @p mats.
template<viewable_matrix... Matrices>
  requires compatible_matrices_v<Matrices...> && (... && bool_matrix<Matrices>)
constexpr auto matrix_or(Matrices&&... mats)
{
  return map(Or{}, std::forward<Matrices>(mats)...);
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
  return map(Identity{}, std::forward<Matrix>(mat));
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

/// @brief Add the matrices @p mat1 and @p mat2.
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
  requires compatible_matrices_v<Matrix1, Matrix2> && //
           numeric_matrix<Matrix1> && numeric_matrix<Matrix2>
constexpr auto operator+(Matrix1&& mat1, Matrix2&& mat2)
{
  return map(Add{}, std::forward<Matrix1>(mat1), std::forward<Matrix2>(mat2));
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

/// @brief Element-wise multiply the matrices @p mat1 and @p mat2.
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
  requires compatible_matrices_v<Matrix1, Matrix2> && //
           numeric_matrix<Matrix1> && numeric_matrix<Matrix2>
constexpr auto operator*(Matrix1&& mat1, Matrix2&& mat2)
{
  return map(Multiply{}, //
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

/// @brief Normalize the matrix @p mat (divide by it's norm).
template<viewable_matrix Matrix>
  requires numeric_matrix<Matrix>
constexpr auto normalize(Matrix&& mat)
{
  auto mat_view = make_matrix_view(std::forward<Matrix>(mat));
  const auto mat_norm = norm_2(mat_view);
  return safe_divide(1.0_dp, mat_norm) * std::move(mat_view);
}

// -----------------------------------------------------------------------------

/// @brief Get the element-wise real component of the matrix @p mat.
template<viewable_matrix Matrix>
  requires numeric_matrix<Matrix>
constexpr auto real(Matrix&& mat)
{
  return map(Real{}, std::forward<Matrix>(mat));
}

/// @brief Get the element-wise imaginary component of the matrix @p mat.
template<viewable_matrix Matrix>
  requires numeric_matrix<Matrix>
constexpr auto imag(Matrix&& mat)
{
  return map(Imag{}, std::forward<Matrix>(mat));
}

/// @brief Element-wise conjugate the matrix @p mat.
template<viewable_matrix Matrix>
  requires numeric_matrix<Matrix>
constexpr auto conj(Matrix&& mat)
{
  return map(Conj{}, std::forward<Matrix>(mat));
}

// -----------------------------------------------------------------------------

/// @brief Compute the element-wise absolute value of the matrix @p mat.
template<viewable_matrix Matrix>
  requires numeric_matrix<Matrix>
constexpr auto abs(Matrix&& mat)
{
  return map(Abs{}, std::forward<Matrix>(mat));
}

/// @brief Compute the element-wise sign of the matrix @p mat.
template<viewable_matrix Matrix>
  requires real_matrix<Matrix> || integer_matrix<Matrix>
constexpr auto sign(Matrix&& mat)
{
  return map(Sign{}, std::forward<Matrix>(mat));
}

// -----------------------------------------------------------------------------

/// @brief Compute the element-wise square root of the matrix @p mat.
template<viewable_matrix Matrix>
  requires numeric_matrix<Matrix>
constexpr auto sqrt(Matrix&& mat)
{
  return map(Sqrt{}, std::forward<Matrix>(mat));
}

/// @brief Compute the element-wise cube root of the matrix @p mat.
template<viewable_matrix Matrix>
  requires numeric_matrix<Matrix>
constexpr auto cbrt(Matrix&& mat)
{
  return map(Cbrt{}, std::forward<Matrix>(mat));
}

template<std::copyable Scalar, viewable_matrix Matrix>
  requires numeric_type<Scalar> && numeric_matrix<Matrix>
constexpr auto pow(Scalar scal, Matrix&& mat)
{
  return map(BindFirst{Pow{}, scal}, std::forward<Matrix>(mat));
}
/// @brief Element-wise raise matrix @p to the power @p power.
template<viewable_matrix Matrix, std::copyable Scalar>
  requires numeric_matrix<Matrix> && numeric_type<Scalar>
constexpr auto pow(Matrix&& mat, Scalar power)
{
  return map(BindLast{Pow{}, power}, std::forward<Matrix>(mat));
}
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
  requires numeric_matrix<Matrix1> && numeric_matrix<Matrix2>
constexpr auto pow(Matrix1&& mat1, Matrix2&& mat2)
{
  return map(Pow{}, std::forward<Matrix1>(mat1), std::forward<Matrix2>(mat2));
}

// -----------------------------------------------------------------------------

/// @brief Compute the element-wise exponent of the matrix @p mat.
template<viewable_matrix Matrix>
  requires numeric_matrix<Matrix>
constexpr auto exp(Matrix&& mat)
{
  return map(Exp{}, std::forward<Matrix>(mat));
}

/// @brief Compute the element-wise exponent (base 2) of the matrix @p mat.
template<viewable_matrix Matrix>
  requires numeric_matrix<Matrix>
constexpr auto exp2(Matrix&& mat)
{
  return map(Exp2{}, std::forward<Matrix>(mat));
}

/// @brief Compute the element-wise logarithm of the matrix @p mat.
template<viewable_matrix Matrix>
  requires numeric_matrix<Matrix>
constexpr auto log(Matrix&& mat)
{
  return map(Log{}, std::forward<Matrix>(mat));
}

/// @brief Compute the element-wise logarithm (base 2) of the matrix @p mat.
template<viewable_matrix Matrix>
  requires numeric_matrix<Matrix>
constexpr auto log2(Matrix&& mat)
{
  return map(Log2{}, std::forward<Matrix>(mat));
}

/// @brief Compute the element-wise logarithm (base 10) of the matrix @p mat.
template<viewable_matrix Matrix>
  requires numeric_matrix<Matrix>
constexpr auto log10(Matrix&& mat)
{
  return map(Log10{}, std::forward<Matrix>(mat));
}

// -----------------------------------------------------------------------------

/// @brief Compute the element-wise sine of the matrix @p mat.
template<viewable_matrix Matrix>
  requires numeric_matrix<Matrix>
constexpr auto sin(Matrix&& mat)
{
  return map(Sin{}, std::forward<Matrix>(mat));
}

/// @brief Compute the element-wise cosine of the matrix @p mat.
template<viewable_matrix Matrix>
  requires numeric_matrix<Matrix>
constexpr auto cos(Matrix&& mat)
{
  return map(Cos{}, std::forward<Matrix>(mat));
}

/// @brief Compute the element-wise tangent of the matrix @p mat.
template<viewable_matrix Matrix>
  requires numeric_matrix<Matrix>
constexpr auto tan(Matrix&& mat)
{
  return map(Tan{}, std::forward<Matrix>(mat));
}

/// @brief Compute the element-wise inverse sine of the matrix @p mat.
template<viewable_matrix Matrix>
  requires numeric_matrix<Matrix>
constexpr auto asin(Matrix&& mat)
{
  return map(Asin{}, std::forward<Matrix>(mat));
}

/// @brief Compute the element-wise inverse cosine of the matrix @p mat.
template<viewable_matrix Matrix>
  requires numeric_matrix<Matrix>
constexpr auto acos(Matrix&& mat)
{
  return map(Acos{}, std::forward<Matrix>(mat));
}

/// @brief Compute the element-wise inverse tangent of the matrix @p mat.
template<viewable_matrix Matrix>
  requires numeric_matrix<Matrix>
constexpr auto atan(Matrix&& mat)
{
  return map(Atan{}, std::forward<Matrix>(mat));
}

// -----------------------------------------------------------------------------

/// @brief Compute the element-wise hyperbolic sine of the matrix @p mat.
template<viewable_matrix Matrix>
  requires numeric_matrix<Matrix>
constexpr auto sinh(Matrix&& mat)
{
  return map(Sinh{}, std::forward<Matrix>(mat));
}

/// @brief Compute the element-wise hyperbolic cosine of the matrix @p mat.
template<viewable_matrix Matrix>
  requires numeric_matrix<Matrix>
constexpr auto cosh(Matrix&& mat)
{
  return map(Cosh{}, std::forward<Matrix>(mat));
}

/// @brief Compute the element-wise hyperbolic tangent of the matrix @p mat.
template<viewable_matrix Matrix>
  requires numeric_matrix<Matrix>
constexpr auto tanh(Matrix&& mat)
{
  return map(Tanh{}, std::forward<Matrix>(mat));
}

/// @brief Compute the element-wise inverse hyperbolic sine of the
/// matrix @p mat.
template<viewable_matrix Matrix>
  requires numeric_matrix<Matrix>
constexpr auto asinh(Matrix&& mat)
{
  return map(Asinh{}, std::forward<Matrix>(mat));
}

/// @brief Compute the element-wise inverse hyperbolic cosine of the
/// matrix @p mat.
template<viewable_matrix Matrix>
  requires numeric_matrix<Matrix>
constexpr auto acosh(Matrix&& mat)
{
  return map(Acosh{}, std::forward<Matrix>(mat));
}

/// @brief Compute the element-wise inverse hyperbolic tangent of the
/// matrix @p mat.
template<viewable_matrix Matrix>
  requires numeric_matrix<Matrix>
constexpr auto atanh(Matrix&& mat)
{
  return map(Atanh{}, std::forward<Matrix>(mat));
}

// -----------------------------------------------------------------------------

} // namespace Storm
