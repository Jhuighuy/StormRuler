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

#include <Storm/Bittern/Math.hpp>

#include <algorithm>
#include <concepts>
#include <limits>
#include <type_traits>
#include <utility>

namespace Storm {

/// @brief Matrix shape: a pair of the row and column extents.
class MatrixShape final {
private:

  size_t num_rows_, num_cols_;

public:

  /// @brief Construct the matrix shape.
  constexpr MatrixShape(size_t num_rows, size_t num_cols) noexcept
      : num_rows_{num_rows}, num_cols_{num_cols} {
    STORM_ASSERT_(num_rows_ != 0 && num_cols_ != 0,
                  "Matrix shape cannot be zero!");
  }

  /// @brief Number of the matrix rows.
  constexpr size_t num_rows() const noexcept {
    return num_rows_;
  }

  /// @brief Number of the matrix columns.
  constexpr size_t num_cols() const noexcept {
    return num_cols_;
  }

  /// @brief Comparison operator.
  constexpr auto operator<=>(const MatrixShape&) const = default;

  /// @brief Check if the @p row_index and @p col_index are in range.
  constexpr bool in_range(size_t row_index, size_t col_index) const noexcept {
    return row_index < num_rows_ && col_index < num_cols_;
  }

}; // class MatrixShape

/// @brief Matrix: has matrix shape and two subscripts.
template<class Matrix>
concept matrix = requires(Matrix& mat) {
                   { mat.shape() } -> std::same_as<MatrixShape>;
                   { mat(size_t{}, size_t{}) };
                 };

/// @brief Number of the matrix rows.
template<matrix Matrix>
constexpr size_t num_rows(Matrix&& mat) noexcept {
  return mat.shape().num_rows();
}

/// @brief Number of the matrix columns.
template<matrix Matrix>
constexpr size_t num_cols(Matrix&& mat) noexcept {
  return mat.shape().num_cols();
}

// -----------------------------------------------------------------------------

/// @brief Matrix element type, as is.
template<matrix Matrix>
using matrix_element_decltype_t =
    decltype(std::declval<Matrix>()(size_t{}, size_t{}));

/// @brief Matrix element type.
template<matrix Matrix>
using matrix_element_t = std::remove_cvref_t<matrix_element_decltype_t<Matrix>>;

/// @brief Matrix element reference type.
template<matrix Matrix>
  requires std::is_lvalue_reference_v<matrix_element_decltype_t<Matrix>>
using matrix_element_ref_t = matrix_element_decltype_t<Matrix>;

/// @brief Matrix with assignable elements.
template<class Matrix>
concept output_matrix = matrix<Matrix>; /// @todo Implement me!

/// @brief Matrix with matrix elements (block matrix).
template<class Matrix>
concept block_matrix = matrix<Matrix> && matrix<matrix_element_t<Matrix>>;

template<matrix Matrix>
struct matrix_inner_element {
  using type = matrix_element_t<Matrix>;
};
template<block_matrix Matrix>
struct matrix_inner_element<Matrix> {
  using type = typename matrix_inner_element<matrix_element_t<Matrix>>::type;
};

/// @brief Matrix inner element type.
/// It is the regular element type for the regular matrices, and
/// the element type of a block if the matrix is a block matrix.
template<matrix Matrix>
using matrix_inner_element_t = typename matrix_inner_element<Matrix>::type;

/// @brief Matrix with boolean elements.
template<class Matrix>
concept bool_matrix =
    matrix<Matrix> && bool_type<matrix_inner_element_t<Matrix>>;

/// @brief Matrix with integral elements.
/// We do not have any special integer matrix operations,
/// integer matrices are defined for completeness.
template<class Matrix>
concept integer_matrix =
    matrix<Matrix> && integer_type<matrix_inner_element_t<Matrix>>;

/// @brief Matrix with real (floating-point) elements.
template<class Matrix>
concept real_matrix =
    matrix<Matrix> && real_type<matrix_inner_element_t<Matrix>>;

/// @brief Matrix with complex (floating-point) elements.
template<class Matrix>
concept complex_matrix =
    matrix<Matrix> && complex_type<matrix_inner_element_t<Matrix>>;

/// @brief Matrix with real or complex (floating-point) elements.
template<class Matrix>
concept real_or_complex_matrix = real_matrix<Matrix> || complex_matrix<Matrix>;

/// @brief Matrix with numerical elements.
template<class Matrix>
concept numeric_matrix =
    integer_matrix<Matrix> || real_or_complex_matrix<Matrix>;

// -----------------------------------------------------------------------------

/// @brief Treat matrices of the specified type as the sparse matrices.
template<class Real>
inline constexpr bool enable_sparse_matrix_v = false;

/// @brief Sparse matrix.
template<class Matrix>
concept sparse_matrix = matrix<Matrix> && enable_sparse_matrix_v<Matrix>;

// -----------------------------------------------------------------------------

/// @brief Reduce the matrix @p mat coefficients to a single value.
/// @param init Initial reduction value.
/// @param reduce_func Reduction function.
/// @todo Restrictions!
/// @{
template<class Value, class ReduceFunc, matrix Matrix>
constexpr auto reduce(Value init, ReduceFunc reduce_func, Matrix&& mat) {
  for (size_t row_index = 0; row_index < num_rows(mat); ++row_index) {
    for (size_t col_index = 0; col_index < num_cols(mat); ++col_index) {
      if constexpr (!block_matrix<Matrix>) {
        init = reduce_func(init, mat(row_index, col_index));
      } else {
        init = reduce(init, reduce_func, mat(row_index, col_index));
      }
    }
  }
  return init;
}
template<class Value, class ReduceFunc, class Func, //
         matrix Matrix, matrix... RestMatrices>
constexpr auto reduce(Value init, ReduceFunc reduce_func, Func func,
                      Matrix&& mat, RestMatrices&&... mats) noexcept {
  STORM_ASSERT_((mat.shape() == mats.shape()) && ...,
                "Matrix shapes doesn't match!");
  for (size_t row_index = 0; row_index < num_rows(mat); ++row_index) {
    for (size_t col_index = 0; col_index < num_cols(mat); ++col_index) {
      if constexpr (!block_matrix<Matrix>) {
        init = reduce_func(init, func(mat(row_index, col_index),
                                      mats(row_index, col_index)...));
      } else {
        init = reduce(init, reduce_func, func, mat(row_index, col_index),
                      mats(row_index, col_index)...);
      }
    }
  }
  return init;
}
/// @}

/// @brief Sum the matrix @p mat elements.
template<real_or_complex_matrix Matrix>
constexpr auto sum(Matrix&& mat) {
  return reduce(matrix_inner_element_t<Matrix>{0.0}, std::plus{},
                std::forward<Matrix>(mat));
}

/// @brief Check if all the boolean matrix @p mat elements are true.
template<bool_matrix Matrix>
  requires std::same_as<matrix_element_t<Matrix>, bool>
constexpr auto all(Matrix&& mat) {
  return reduce(true, std::logical_and{}, std::forward<Matrix>(mat));
}
/// @brief Check if any of the boolean matrix @p mat elements is true.
template<bool_matrix Matrix>
  requires std::same_as<matrix_element_t<Matrix>, bool>
constexpr auto any(Matrix&& mat) {
  return reduce(false, std::logical_or{}, std::forward<Matrix>(mat));
}

/// @brief Minimum matrix @p mat element.
template<matrix Matrix>
  requires std::totally_ordered<matrix_element_t<Matrix>>
constexpr auto min_element(Matrix&& mat) {
  using Elem = matrix_element_t<Matrix>;
  return reduce(
      std::numeric_limits<Elem>::max(),
      [](const Elem& a, const Elem& b) noexcept { return min(a, b); },
      std::forward<Matrix>(mat));
}
/// @brief Maximum matrix @p mat element.
template<matrix Matrix>
  requires std::totally_ordered<matrix_element_t<Matrix>>
constexpr auto max_element(Matrix&& mat) {
  using Elem = matrix_element_t<Matrix>;
  return reduce(
      std::numeric_limits<Elem>::lowest(),
      [](const Elem& a, const Elem& b) noexcept { return max(a, b); },
      std::forward<Matrix>(mat));
}

/// @brief Element-wise matrix @p mat \f$ L_{1} \f$-norm.
template<numeric_matrix Matrix>
constexpr auto norm_1(Matrix&& mat) {
  using Result = decltype(abs(std::declval<matrix_element_t<Matrix>>()));
  return reduce(
      Result{0.0}, std::plus{},
      []<class Elem>(Elem&& elem) noexcept {
        return abs(std::forward<Elem>(elem));
      },
      std::forward<Matrix>(mat));
}
/// @brief Element-wise matrix @p mat \f$ L_{2} \f$-norm.
template<numeric_matrix Matrix>
constexpr auto norm_2(Matrix&& mat) {
  using Result = decltype(abs(std::declval<matrix_element_t<Matrix>>()));
  const auto temp = reduce(
      Result{0.0}, std::plus{},
      [](auto elem) noexcept { return real(elem * conj(elem)); },
      std::forward<Matrix>(mat));
  return sqrt(temp);
}
/// @brief Element-wise matrix @p mat \f$ L_{p} \f$-norm.
template<numeric_matrix Matrix>
constexpr auto norm_p(Matrix&& mat, real_t p) {
  STORM_ASSERT_(p > 0.0, "Invalid p-norm parameter!");
  using Result = decltype(abs(std::declval<matrix_element_t<Matrix>>()));
  const auto temp = reduce(
      Result{0.0}, std::plus{},
      [p]<class Elem>(Elem&& elem) noexcept {
        return pow(abs(std::forward<Elem>(elem)), p);
      },
      std::forward<Matrix>(mat));
  return pow(temp, 1.0 / p);
}
/// @brief Element-wise matrix @p mat \f$ L_{\infty} \f$-norm.
template<numeric_matrix Matrix>
constexpr auto norm_inf(Matrix&& mat) {
  using Result = decltype(abs(std::declval<matrix_element_t<Matrix>>()));
  return reduce(
      Result{0.0},
      [](const Result& a, const Result& b) noexcept { return max(a, b); },
      []<class Elem>(Elem&& elem) noexcept {
        return abs(std::forward<Elem>(elem));
      },
      std::forward<Matrix>(mat));
}

/// @brief Element-wise dot product of the matrices @p mat1 and @p mat2.
template<numeric_matrix Matrix1, numeric_matrix Matrix2>
constexpr auto dot_product(Matrix1&& mat1, Matrix2&& mat2) noexcept {
  using Result = decltype(std::declval<matrix_element_t<Matrix1>>() *
                          std::declval<matrix_element_t<Matrix2>>());
  return reduce(
      Result{0.0}, std::plus{},
      []<class Elem1, class Elem2>(Elem1&& elem1, Elem2&& elem2) {
        return std::forward<Elem1>(elem1) * conj(std::forward<Elem2>(elem2));
      },
      std::forward<Matrix1>(mat1), std::forward<Matrix2>(mat2));
}

/// @copydoc norm_2
template<numeric_matrix Matrix>
constexpr auto length(Matrix&& mat) {
  return norm_2(std::forward<Matrix>(mat));
}

} // namespace Storm

#include <Storm/Bittern/MatrixAlgorithms.hpp>
#include <Storm/Bittern/MatrixIo.hpp>
#include <Storm/Bittern/MatrixMath.hpp>
#include <Storm/Bittern/MatrixProduct.hpp>
#include <Storm/Bittern/MatrixSlicing.hpp>
#include <Storm/Bittern/MatrixTranspose.hpp>
