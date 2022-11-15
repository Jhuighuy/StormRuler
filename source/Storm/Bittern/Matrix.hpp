// Copyright © 2020 - 2023 Oleg Butakov
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

#include <Storm/Utils/Math.hpp>

#include <algorithm>
#include <concepts>
#include <limits>
#include <random>
#include <type_traits>
#include <utility>

namespace Storm {

/// @brief Fixed matrix extent.
template<size_t N>
  requires (N != 0)
using fixed_matrix_extent_t = std::integral_constant<size_t, N>;

template<class>
inline constexpr bool is_fixed_matrix_extent_v = false;
template<size_t N>
  requires (N != 0)
inline constexpr bool is_fixed_matrix_extent_v<fixed_matrix_extent_t<N>> = true;

/// @brief Fixed matrix extent concept.
template<class FixedExtent>
concept fixed_matrix_extent = is_fixed_matrix_extent_v<FixedExtent>;

/// @brief Dynamic matrix extent.
using dynamic_matrix_extent_t = size_t;

/// @brief Dynamic matrix extent concept.
template<class DynamicExtent>
concept dynamic_matrix_extent = std::same_as<DynamicExtent, //
                                             dynamic_matrix_extent_t>;

/// @brief Fixed or dynamic matrix extent.
template<class Extent>
concept matrix_extent = fixed_matrix_extent<Extent> || //
                        dynamic_matrix_extent<Extent>;

/// @brief Matrix shape: a pair of the row and column extents.
template<matrix_extent RowsExtent, matrix_extent ColsExtent>
class MatrixShape final {
public:

  /// @brief Number of the matrix rows.
  STORM_NO_UNIQUE_ADDRESS_ RowsExtent num_rows{};

  /// @brief Number of the matrix columns.
  STORM_NO_UNIQUE_ADDRESS_ ColsExtent num_cols{};

  /// @brief Check if the @p row_index and @p col_index are in range.
  [[nodiscard]] constexpr bool in_range(size_t row_index,
                                        size_t col_index) const noexcept {
    return row_index < num_rows && col_index < num_cols;
  }

  /// @brief Comparison operators.
  /// @{
  template<matrix_extent OtherRowsExtent, matrix_extent OtherColsExtent>
  [[nodiscard]] friend constexpr bool
  operator==(MatrixShape shape1,
             MatrixShape<OtherRowsExtent, OtherColsExtent> shape2) noexcept {
    return shape1.num_rows == shape2.num_rows &&
           shape1.num_cols == shape2.num_cols;
  }
  template<matrix_extent OtherRowsExtent, matrix_extent OtherColsExtent>
  [[nodiscard]] friend constexpr bool
  operator!=(MatrixShape shape1,
             MatrixShape<OtherRowsExtent, OtherColsExtent> shape2) noexcept {
    return shape1.num_rows != shape2.num_rows ||
           shape1.num_cols != shape2.num_cols;
  }
  /// @}

}; // class MatrixShape

template<fixed_matrix_extent RowsExtent, std::integral ColsExtent>
MatrixShape(RowsExtent, ColsExtent) -> MatrixShape<RowsExtent, size_t>;
template<std::integral RowsExtent, fixed_matrix_extent ColsExtent>
MatrixShape(RowsExtent, ColsExtent) -> MatrixShape<size_t, ColsExtent>;
template<std::integral RowsExtent, std::integral ColsExtent>
MatrixShape(RowsExtent, ColsExtent) -> MatrixShape<size_t, size_t>;

template<class>
inline constexpr bool is_matrix_shape_v = false;
template<class RowsExtent, class ColsExtent>
inline constexpr bool is_matrix_shape_v<MatrixShape<RowsExtent, //
                                                    ColsExtent>> = true;

/// @brief Matrix shape concept.
template<class MatrixShape>
concept matrix_shape = is_matrix_shape_v<MatrixShape>;

/// @brief Fixed matrix shape.
template<size_t NumRows, size_t NumCols>
  requires (NumRows != 0 && NumCols != 0)
using FixedMatrixShape = MatrixShape<fixed_matrix_extent_t<NumRows>, //
                                     fixed_matrix_extent_t<NumCols>>;

// -----------------------------------------------------------------------------

template<class Matrix>
struct matrix_row_index {
  using type = size_t;
};
/*template<class Matrix>
  requires requires { typename Matrix::row_index_type; }
struct matrix_row_index<Matrix> {
  using type = typename Matrix::row_index_type;
};*/

template<class Matrix>
struct matrix_col_index {
  using type = size_t;
};
/*template<class Matrix>
  requires requires { typename Matrix::col_index_type; }
struct matrix_col_index<Matrix> {
  using type = typename Matrix::col_index_type;
};*/

/// @brief Matrix row index type.
template<class Matrix>
using matrix_row_index_t = typename matrix_row_index<Matrix>::type;
/// @brief Matrix column index type.
template<class Matrix>
using matrix_col_index_t = typename matrix_col_index<Matrix>::type;

/// @brief Matrix: has matrix shape and two subscripts.
// clang-format off
template<class Matrix>
concept matrix = 
    requires(Matrix& mat) {
      { mat.shape() } noexcept -> matrix_shape;
      { mat(std::declval<matrix_row_index_t<Matrix>>(),
            std::declval<matrix_col_index_t<Matrix>>()) } noexcept;
    };
// clang-format on

/// @brief Number of the matrix rows.
template<matrix Matrix>
[[nodiscard]] constexpr size_t num_rows(Matrix&& mat) noexcept {
  return mat.shape().num_rows;
}

/// @brief Number of the matrix columns.
template<matrix Matrix>
[[nodiscard]] constexpr size_t num_cols(Matrix&& mat) noexcept {
  return mat.shape().num_cols;
}

/// @brief Matrix element type, as is.
template<matrix Matrix>
using matrix_element_decltype_t = decltype( //
    std::declval<Matrix>()(std::declval<matrix_row_index_t<Matrix>>(),
                           std::declval<matrix_col_index_t<Matrix>>()));

/// @brief Matrix element type.
template<matrix Matrix>
using matrix_element_t = std::remove_cvref_t<matrix_element_decltype_t<Matrix>>;

/// @brief Matrix element reference type.
template<matrix Matrix>
  requires std::is_lvalue_reference_v<matrix_element_decltype_t<Matrix>>
using matrix_element_ref_t = matrix_element_decltype_t<Matrix>;

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

/// @brief Treat matrices of the specified type as the sparse matrices.
template<class Real>
inline constexpr bool enable_sparse_matrix_v = false;

/// @brief Sparse matrix.
template<class Matrix>
concept sparse_matrix = matrix<Matrix> && enable_sparse_matrix_v<Matrix>;

/// @brief Treat matrices of the specified type as the boolean matrices.
template<class Bool>
inline constexpr bool enable_bool_matrix_for_v = std::is_same_v<Bool, bool>;

/// @brief Matrix with boolean elements.
template<class Matrix>
concept bool_matrix =
    matrix<Matrix> && enable_bool_matrix_for_v<matrix_inner_element_t<Matrix>>;

/// @brief Treat matrices of the specified type as the integer matrices.
template<class Integer>
inline constexpr bool enable_integer_matrix_for_v = std::is_integral_v<Integer>;

/// @brief Matrix with integral elements.
/// We do not have any special integer matrix operations,
/// integer matrices are defined for completeness.
template<class Matrix>
concept integer_matrix =
    matrix<Matrix> &&
    enable_integer_matrix_for_v<matrix_inner_element_t<Matrix>>;

/// @brief Treat matrices of the specified type as the real matrices.
template<class Real>
inline constexpr bool enable_real_matrix_for_v = std::is_floating_point_v<Real>;

/// @brief Matrix with real (floating-point) elements.
template<class Matrix>
concept real_matrix =
    matrix<Matrix> && enable_real_matrix_for_v<matrix_inner_element_t<Matrix>>;

/// @brief Treat matrices of the specified type as the complex matrices.
template<class Complex>
inline constexpr bool enable_complex_matrix_for_v =
    is_complex_floating_point_v<Complex>;

/// @brief Matrix with complex (floating-point) elements.
template<class Matrix>
concept complex_matrix =
    matrix<Matrix> &&
    enable_complex_matrix_for_v<matrix_inner_element_t<Matrix>>;

/// @brief Matrix with real or complex (floating-point) elements.
template<class Matrix>
concept real_or_complex_matrix = real_matrix<Matrix> || complex_matrix<Matrix>;

/// @brief Matrix with numerical elements.
template<class Matrix>
concept numeric_matrix =
    integer_matrix<Matrix> || real_or_complex_matrix<Matrix>;

// -----------------------------------------------------------------------------

/// @brief Assign the matrices.
/// @todo Restrictions!
/// @{
template<matrix OutMatrix, matrix Matrix>
constexpr OutMatrix& assign(OutMatrix& out_mat, Matrix&& mat) noexcept {
  STORM_ASSERT_(out_mat.shape() == mat.shape(), "Matrix shapes doesn't match!");
  for (size_t row_index = 0; row_index < num_rows(out_mat); ++row_index) {
    for (size_t col_index = 0; col_index < num_cols(out_mat); ++col_index) {
      out_mat(row_index, col_index) = mat(row_index, col_index);
    }
  }
  return out_mat;
}
template<matrix OutMatrix, class AssignFunc, matrix... Matrices>
constexpr OutMatrix& assign(OutMatrix& out_mat, AssignFunc assign_func,
                            Matrices&&... mats) noexcept {
  STORM_ASSERT_((out_mat.shape() == mats.shape()) && ...,
                "Matrix shapes doesn't match!");
  for (size_t row_index = 0; row_index < num_rows(out_mat); ++row_index) {
    for (size_t col_index = 0; col_index < num_cols(out_mat); ++col_index) {
      assign_func(out_mat(row_index, col_index), mats(row_index, col_index)...);
    }
  }
  return out_mat;
}
/// @}

/// @brief Assign the matrices.
template<matrix OutMatrix, matrix Matrix>
constexpr OutMatrix& operator<<=(OutMatrix& out_mat, Matrix&& mat) noexcept {
  return assign(out_mat, std::forward<Matrix>(mat));
}

/// @brief Fill the matrix @p out_mat with a scalar @p scal.
template<matrix OutMatrix, std::copyable Scalar>
  requires (!matrix<Scalar>)
constexpr OutMatrix& fill(OutMatrix& out_mat, Scalar scal) {
  return assign(out_mat, [scal = std::move(scal)](auto& out_elem) noexcept {
    out_elem = scal;
  });
}

/// @brief Fill the matrix @p out_mat elements with the random numbers.
/// @warning This is a sequential operation!
template<real_matrix OutMatrix>
constexpr OutMatrix&
fill_randomly(OutMatrix&& out_mat, //
              matrix_element_t<OutMatrix> min = 0,
              matrix_element_t<OutMatrix> max = 1) noexcept {
  static std::mt19937_64 random_engine{};
  std::uniform_real_distribution distribution{min, max};
  for (size_t row_index = 0; row_index < num_rows(out_mat); ++row_index) {
    for (size_t col_index = 0; col_index < num_cols(out_mat); ++col_index) {
      out_mat(row_index, col_index) = distribution(random_engine);
    }
  }
  return out_mat;
}

/// @brief Multiply-assign the matrix @p out_mat by a scalar @p scal.
template<numeric_matrix OutMatrix, std::copyable Scalar>
  requires (!matrix<Scalar>)
constexpr OutMatrix& operator*=(OutMatrix& out_mat, Scalar scal) {
  return assign(out_mat, [scal = std::move(scal)](auto& out_elem) noexcept {
    out_elem *= scal;
  });
}

/// @brief Divide-assign the matrix @p out_mat by a scalar @p scal.
template<numeric_matrix OutMatrix, std::copyable Scalar>
  requires (!matrix<Scalar>)
constexpr OutMatrix& operator/=(OutMatrix& out_mat, Scalar scal) {
  return assign(out_mat, [scal = std::move(scal)](auto& out_elem) noexcept {
    out_elem /= scal;
  });
}

/// @brief Add-assign the matrices @p out_mat and @p mat
template<numeric_matrix OutMatrix, matrix Matrix>
constexpr OutMatrix& operator+=(OutMatrix& out_mat, Matrix&& mat) {
  return assign(
      out_mat,
      []<class Elem>(auto& out_elem, Elem&& elem) noexcept {
        out_elem += std::forward<Elem>(elem);
      },
      std::forward<Matrix>(mat));
}

/// @brief Subtract-assign the matrices @p out_mat and @p mat.
template<numeric_matrix OutMatrix, matrix Matrix>
constexpr OutMatrix& operator-=(OutMatrix& out_mat, Matrix&& mat) {
  return assign(
      out_mat,
      []<class Elem>(auto& out_elem, Elem&& elem) noexcept {
        out_elem -= std::forward<Elem>(elem);
      },
      std::forward<Matrix>(mat));
}

/// @brief Element-wise multiply-assign the matrices @p out_mat and @p mat.
template<numeric_matrix OutMatrix, matrix Matrix>
constexpr OutMatrix& operator*=(OutMatrix& out_mat, Matrix&& mat) {
  return assign(
      out_mat,
      []<class Elem>(auto& out_elem, Elem&& elem) noexcept {
        out_elem *= std::forward<Elem>(elem);
      },
      std::forward<Matrix>(mat));
}

/// @brief Element-wise divide-assign the matrices @p out_mat and @p mat.
template<numeric_matrix OutMatrix, matrix Matrix>
constexpr OutMatrix& operator/=(OutMatrix& out_mat, Matrix&& mat) {
  return assign(
      out_mat,
      []<class Elem>(auto& out_elem, Elem&& elem) noexcept {
        out_elem /= std::forward<Elem>(elem);
      },
      std::forward<Matrix>(mat));
}

// -----------------------------------------------------------------------------

/// @brief Reduce the matrix @p mat coefficients to a single value.
/// @param init Initial reduction value.
/// @param reduce_func Reduction function.
/// @todo Restrictions!
/// @{
template<class Value, class ReduceFunc, matrix Matrix>
[[nodiscard]] constexpr auto reduce(Value init, ReduceFunc reduce_func,
                                    Matrix&& mat) {
  for (size_t row_index = 0; row_index < num_rows(mat); ++row_index) {
    for (size_t col_index = 0; col_index < num_cols(mat); ++col_index) {
      if constexpr (block_matrix<Matrix>) {
        init = reduce(init, reduce_func, mat(row_index, col_index));
      } else {
        init = reduce_func(init, mat(row_index, col_index));
      }
    }
  }
  return init;
}
template<class Value, class ReduceFunc, class Func, //
         matrix Matrix, matrix... RestMatrices>
[[nodiscard]] constexpr auto reduce(Value init, ReduceFunc reduce_func,
                                    Func func, Matrix&& mat,
                                    RestMatrices&&... mats) noexcept {
  STORM_ASSERT_((mat.shape() == mats.shape()) && ...,
                "Matrix shapes doesn't match!");
  for (size_t row_index = 0; row_index < num_rows(mat); ++row_index) {
    for (size_t col_index = 0; col_index < num_cols(mat); ++col_index) {
      if constexpr (block_matrix<Matrix>) {
        init = reduce(init, reduce_func, func,   //
                      mat(row_index, col_index), //
                      mats(row_index, col_index)...);
      } else {
        init = reduce_func(init, func(mat(row_index, col_index),
                                      mats(row_index, col_index)...));
      }
    }
  }
  return init;
}
/// @}

/// @brief Sum the matrix @p mat elements.
template<real_or_complex_matrix Matrix>
[[nodiscard]] constexpr auto sum(Matrix&& mat) {
  return reduce(matrix_inner_element_t<Matrix>{0.0}, std::plus{},
                std::forward<Matrix>(mat));
}

/// @brief Check if all the boolean matrix @p mat elements are true.
template<bool_matrix Matrix>
  requires std::same_as<matrix_element_t<Matrix>, bool>
[[nodiscard]] constexpr auto all(Matrix&& mat) {
  return reduce(true, std::logical_and{}, std::forward<Matrix>(mat));
}
/// @brief Check if any of the boolean matrix @p mat elements is true.
template<bool_matrix Matrix>
  requires std::same_as<matrix_element_t<Matrix>, bool>
[[nodiscard]] constexpr auto any(Matrix&& mat) {
  return reduce(false, std::logical_or{}, std::forward<Matrix>(mat));
}

/// @brief Minimum matrix @p mat element.
template<matrix Matrix>
  requires std::totally_ordered<matrix_element_t<Matrix>>
[[nodiscard]] constexpr auto min_element(Matrix&& mat) {
  using Elem = matrix_element_t<Matrix>;
  return reduce(
      std::numeric_limits<Elem>::max(),
      [](const Elem& a, const Elem& b) noexcept { return min(a, b); },
      std::forward<Matrix>(mat));
}
/// @brief Maximum matrix @p mat element.
template<matrix Matrix>
  requires std::totally_ordered<matrix_element_t<Matrix>>
[[nodiscard]] constexpr auto max_element(Matrix&& mat) {
  using Elem = matrix_element_t<Matrix>;
  return reduce(
      std::numeric_limits<Elem>::lowest(),
      [](const Elem& a, const Elem& b) noexcept { return max(a, b); },
      std::forward<Matrix>(mat));
}

/// @brief Element-wise matrix @p mat \f$ L_{1} \f$-norm.
template<numeric_matrix Matrix>
[[nodiscard]] constexpr auto norm_1(Matrix&& mat) {
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
[[nodiscard]] constexpr auto norm_2(Matrix&& mat) {
  using Result = decltype(abs(std::declval<matrix_element_t<Matrix>>()));
  const auto temp = reduce(
      Result{0.0}, std::plus{},
      [](auto elem) noexcept { return real(elem * conj(elem)); },
      std::forward<Matrix>(mat));
  return sqrt(temp);
}
/// @brief Element-wise matrix @p mat \f$ L_{p} \f$-norm.
template<numeric_matrix Matrix>
[[nodiscard]] constexpr auto norm_p(Matrix&& mat, real_t p) {
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
[[nodiscard]] constexpr auto norm_inf(Matrix&& mat) {
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

} // namespace Storm

// -----------------------------------------------------------------------------

#include "MatrixIo.hpp"
#include <Storm/Bittern/MatrixView_.hpp>
