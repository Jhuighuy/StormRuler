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
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <Storm/Base.hpp>

#include <Storm/Crow/ConceptUtils.hpp>
#include <Storm/Crow/FunctionalUtils.hpp>
#include <Storm/Crow/MathUtils.hpp>

#include <Storm/Bittern/Matrix.hpp>
#include <Storm/Bittern/Shape.hpp>

#include <concepts>
#include <functional>
#include <limits>
#include <ostream>
#include <random>
#include <utility>

namespace Storm {

// -----------------------------------------------------------------------------

/// @brief Print a @p mat.
/// @todo This is too trivial implementation, we need something fancier :)
std::ostream& operator<<(std::ostream& out, const matrix auto& mat) {
  for (size_t row_index = 0; row_index < num_rows(mat); ++row_index) {
    out << "( ";
    for (size_t col_index = 0; col_index < num_cols(mat); ++col_index) {
      out << mat(row_index, col_index) << " ";
    }
    out << ")" << std::endl;
  }
  return out;
}

// -----------------------------------------------------------------------------

template<class Func, matrix Matrix, matrix... RestMatrices>
  requires compatible_matrices_v<Matrix, RestMatrices...> &&
           std::invocable<Func, //
                          matrix_element_ref_t<Matrix>,
                          matrix_element_ref_t<RestMatrices>...>
constexpr void matrix_for_each(Func func, //
                               Matrix&& mat,
                               RestMatrices&&... rest_mats) noexcept {
  STORM_ASSERT((mat.shape() == rest_mats.shape()) && ...,
               "Matrix shapes do not match!");
  if constexpr (matrix_r<Matrix, 1>) {
    for (size_t row_index = 0; row_index < num_rows(mat); ++row_index) {
      func(mat(row_index), rest_mats(row_index)...);
    }
  } else {
    /// @todo Below:
    static_assert(matrix_r<Matrix, 2>, "General rank is not implemented yet!");
    for (size_t row_index = 0; row_index < num_rows(mat); ++row_index) {
      for (size_t col_index = 0; col_index < num_cols(mat); ++col_index) {
        func(mat(row_index, col_index), rest_mats(row_index, col_index)...);
      }
    }
  }
}

// -----------------------------------------------------------------------------

/// @brief Move the matrix @p mat elements into the matrix @p out_mat.
template<output_matrix OutMatrix, matrix Matrix>
  requires assignable_matrix<OutMatrix, Matrix>
constexpr OutMatrix&& move_elements(OutMatrix&& out_mat,
                                    Matrix&& mat) noexcept {
  matrix_for_each(MoveAssign{}, out_mat, std::forward<Matrix>(mat));
  return std::forward<OutMatrix>(out_mat);
}

/// @brief Copy the matrix @p mat elements to the matrix @p out_mat.
template<output_matrix OutMatrix, matrix Matrix>
  requires assignable_matrix<OutMatrix, Matrix>
constexpr OutMatrix&& copy_elements(OutMatrix&& out_mat,
                                    Matrix&& mat) noexcept {
  matrix_for_each(CopyAssign{}, out_mat, std::forward<Matrix>(mat));
  return std::forward<OutMatrix>(out_mat);
}

/// @brief Evaluate the matrix @p out_mat elements with
/// action function @p action_func and arguments @p mats.
template<class ActionFunc, output_matrix OutMatrix, matrix... Matrices>
  requires compatible_matrices_v<OutMatrix, Matrices...> &&
           std::regular_invocable<ActionFunc, //
                                  matrix_element_ref_t<OutMatrix>,
                                  matrix_element_ref_t<Matrices>...>
constexpr OutMatrix&& eval_elements(ActionFunc action_func, //
                                    OutMatrix&& out_mat,
                                    Matrices&&... mats) noexcept {
  matrix_for_each(std::move(action_func), //
                  out_mat, std::forward<Matrices>(mats)...);
  return std::forward<OutMatrix>(out_mat);
}

/// @brief Assign the matrices.
/// @todo To be removed!
template<output_matrix OutMatrix, matrix Matrix>
constexpr OutMatrix& operator<<=(OutMatrix&& out_mat, Matrix&& mat) noexcept {
  return copy_elements(std::forward<OutMatrix>(out_mat),
                       std::forward<Matrix>(mat));
}

// -----------------------------------------------------------------------------

/// @brief Fill the matrix @p out_mat with a scalar @p scal.
template<matrix OutMatrix, std::copyable Scalar>
  requires std::assignable_from<matrix_element_ref_t<OutMatrix>, Scalar>
constexpr OutMatrix& fill(OutMatrix& out_mat, Scalar scal) {
  return eval_elements(
      [scal = std::move(scal)](auto& out_elem) noexcept { out_elem = scal; },
      out_mat);
}

/// @brief Fill the matrix @p out_mat elements with
/// the uniformly-distributed random numbers.
/// @note This is a sequential operation!
template<output_matrix Matrix>
STORM_CPP23_CONSTEXPR Matrix&
fill_randomly(Matrix&& out_mat, //
              matrix_element_t<Matrix> min = 0,
              matrix_element_t<Matrix> max = 1) noexcept {
  static std::mt19937_64 random_engine{}; // NOSONAR
  std::uniform_real_distribution distribution{min, max};
  for (size_t row_index = 0; row_index < num_rows(out_mat); ++row_index) {
    for (size_t col_index = 0; col_index < num_cols(out_mat); ++col_index) {
      out_mat(row_index, col_index) = distribution(random_engine);
    }
  }
  return out_mat;
}

// -----------------------------------------------------------------------------

/// @brief Reduce the matrix @p mat coefficients to a single value.
/// @param init Initial reduction value.
/// @param reduce_func Reduction function.
/// @todo Restrictions!
/// @{
template<class Value, class ReduceFunc, matrix_r<1> Matrix>
constexpr auto reduce(Value init, ReduceFunc reduce_func, Matrix&& mat) {
  for (size_t row_index = 0; row_index < num_rows(mat); ++row_index) {
    init = reduce_func(init, mat(row_index));
  }
  return init;
}
template<class Value, class ReduceFunc, matrix_r<2> Matrix>
constexpr auto reduce(Value init, ReduceFunc reduce_func, Matrix&& mat) {
  for (size_t row_index = 0; row_index < num_rows(mat); ++row_index) {
    for (size_t col_index = 0; col_index < num_cols(mat); ++col_index) {
      init = reduce_func(init, mat(row_index, col_index));
    }
  }
  return init;
}

template<class Value, class ReduceFunc, class Func, //
         matrix_r<1> Matrix, matrix... RestMatrices>
  requires compatible_matrices_v<Matrix, RestMatrices...>
constexpr auto reduce(Value init, ReduceFunc reduce_func, Func func,
                      Matrix&& mat, RestMatrices&&... mats) noexcept {
  STORM_ASSERT((mat.shape() == mats.shape()) && ...,
               "Matrix shapes doesn't match!");
  for (size_t row_index = 0; row_index < num_rows(mat); ++row_index) {
    init = reduce_func(init, func(mat(row_index), mats(row_index)...));
  }
  return init;
}
template<class Value, class ReduceFunc, class Func, //
         matrix_r<2> Matrix, matrix... RestMatrices>
  requires compatible_matrices_v<Matrix, RestMatrices...>
constexpr auto reduce(Value init, ReduceFunc reduce_func, Func func,
                      Matrix&& mat, RestMatrices&&... mats) noexcept {
  STORM_ASSERT((mat.shape() == mats.shape()) && ...,
               "Matrix shapes doesn't match!");
  for (size_t row_index = 0; row_index < num_rows(mat); ++row_index) {
    for (size_t col_index = 0; col_index < num_cols(mat); ++col_index) {
      init = reduce_func(
          init, func(mat(row_index, col_index), mats(row_index, col_index)...));
    }
  }
  return init;
}
/// @}

// -----------------------------------------------------------------------------

/// @brief Sum the matrix @p mat elements.
template<matrix Matrix>
constexpr auto sum(Matrix&& mat) {
  return reduce(matrix_element_t<Matrix>{0.0}, Add{},
                std::forward<Matrix>(mat));
}

// -----------------------------------------------------------------------------

/// @brief Check if all the boolean matrix @p mat elements are true.
template<matrix Matrix>
constexpr auto all(Matrix&& mat) {
  return reduce(true, LogicalAnd{}, std::forward<Matrix>(mat));
}

/// @brief Check if any of the boolean matrix @p mat elements is true.
template<matrix Matrix>
constexpr auto any(Matrix&& mat) {
  return reduce(false, LogicalOr{}, std::forward<Matrix>(mat));
}

// -----------------------------------------------------------------------------

/// @brief Minimum matrix @p mat element.
template<matrix Matrix>
  requires std::totally_ordered<matrix_element_t<Matrix>>
constexpr auto min_element(Matrix&& mat) {
  using Elem = matrix_element_t<Matrix>;
  return reduce(std::numeric_limits<Elem>::max(), Min{},
                std::forward<Matrix>(mat));
}

/// @brief Maximum matrix @p mat element.
template<matrix Matrix>
  requires std::totally_ordered<matrix_element_t<Matrix>>
constexpr auto max_element(Matrix&& mat) {
  using Elem = matrix_element_t<Matrix>;
  return reduce(std::numeric_limits<Elem>::lowest(), Max{},
                std::forward<Matrix>(mat));
}

// -----------------------------------------------------------------------------

/// @brief Element-wise matrix @p mat @f$ L_{1} @f$-norm.
template<matrix Matrix>
constexpr auto norm_1(Matrix&& mat) {
  using Result = std::invoke_result_t<Abs, matrix_element_t<Matrix>>;
  return reduce(Result{}, Add{}, Abs{}, std::forward<Matrix>(mat));
}

/// @brief Element-wise matrix @p mat @f$ L_{2} @f$-norm.
/// @{
template<matrix Matrix>
constexpr auto norm_2_2(Matrix&& mat) {
  using Result = std::invoke_result_t<Abs, matrix_element_t<Matrix>>;
  return reduce(Result{}, Add{}, AbsSquared{}, std::forward<Matrix>(mat));
}
template<matrix Matrix>
constexpr auto norm_2(Matrix&& mat) {
  return sqrt(norm_2_2(std::forward<Matrix>(mat)));
}
/// @}

/// @brief Element-wise matrix @p mat @f$ L_{p} @f$-norm.
/// @{
template<matrix Matrix, scalar Number>
constexpr auto norm_p_p(Matrix&& mat, Number p) {
  STORM_ASSERT(p > 0, "Invalid p-norm parameter!");
  using Result = std::invoke_result_t<Abs, matrix_element_t<Matrix>>;
  return reduce(Result{}, Add{}, Compose{Abs{}, BindLast{Pow{}, std::move(p)}},
                std::forward<Matrix>(mat));
}
template<matrix Matrix, scalar Number>
constexpr auto norm_p(Matrix&& mat, Number p) {
  /// @todo We should have something better than `1.0 / p`.
  return pow(norm_p_p(std::forward<Matrix>(mat), p), 1.0 / p);
}
/// @}

/// @brief Element-wise matrix @p mat @f$ L_{\infty} @f$-norm.
template<matrix Matrix>
constexpr auto norm_inf(Matrix&& mat) {
  using Result = std::invoke_result_t<Abs, matrix_element_t<Matrix>>;
  return reduce(Result{}, Max{}, Abs{}, std::forward<Matrix>(mat));
}

/// @copydoc norm_2_2
template<matrix Matrix>
constexpr auto length_2(Matrix&& mat) {
  return norm_2_2(std::forward<Matrix>(mat));
}
/// @copydoc norm_2
template<matrix Matrix>
constexpr auto length(Matrix&& mat) {
  return norm_2(std::forward<Matrix>(mat));
}

// -----------------------------------------------------------------------------

/// @brief Element-wise dot product of the matrices @p mat1 and @p mat2.
template<matrix Matrix1, matrix Matrix2>
  requires compatible_matrices_v<Matrix1, Matrix2>
constexpr auto dot_product(Matrix1&& mat1, Matrix2&& mat2) noexcept {
  using Result = std::invoke_result_t< //
      DotProduct, matrix_element_t<Matrix1>, matrix_element_t<Matrix2>>;
  return reduce(Result{}, Add{}, DotProduct{}, //
                std::forward<Matrix1>(mat1), std::forward<Matrix2>(mat2));
}

// -----------------------------------------------------------------------------

} // namespace Storm
