/// Copyright (C) 2022 Oleg Butakov
///
/// Permission is hereby granted, free of charge, to any person obtaining a copy
/// of this software and associated documentation files (the "Software"), to
/// deal in the Software without restriction, including without limitation the
/// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
/// sell copies of the Software, and to permit persons to whom the Software is
/// furnished to do so, subject to the following conditions:
///
/// The above copyright notice and this permission notice shall be included in
/// all copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
/// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
/// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
/// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
/// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
/// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
/// IN THE SOFTWARE.

// This header should really be <Matrix.hpp>
#pragma once

#include <concepts>
#include <type_traits>

#if _OPENMP
#include <omp.h>
#endif

#include <Storm/Base.hpp>

namespace Storm {

/// @brief Matrix extent type.
template<class MatrixExtent>
concept matrix_extent =
    std::same_as<MatrixExtent, size_t> || is_size_t_constant_v<MatrixExtent>;

/// @brief Matrix shape.
template<matrix_extent Rows, matrix_extent Cols>
class MatrixShape final {
public:

  /// @brief Number of the matrix rows.
  STORM_NO_UNIQUE_ADDRESS_ Rows num_rows{};

  /// @brief Number of the matrix columns.
  STORM_NO_UNIQUE_ADDRESS_ Cols num_cols{};

}; // class MatrixShape

template<size_t NumRows, size_t NumCols>
MatrixShape(size_t_constant<NumRows>, size_t_constant<NumCols>)
    -> MatrixShape<size_t_constant<NumRows>, size_t_constant<NumCols>>;

template<size_t NumRows, std::integral Cols>
MatrixShape(size_t_constant<NumRows>, Cols)
    -> MatrixShape<size_t_constant<NumRows>, size_t>;

template<std::integral Rows, size_t NumCols>
MatrixShape(Rows, size_t_constant<NumCols>)
    -> MatrixShape<size_t, size_t_constant<NumCols>>;

template<std::integral Rows, std::integral Cols>
MatrixShape(Rows, Cols) -> MatrixShape<size_t, size_t>;

/// @brief Compare the matrix shapes @p a and @p b.
/// @{
template<matrix_extent Rows1, matrix_extent Cols1, //
         matrix_extent Rows2, matrix_extent Cols2>
[[nodiscard]] constexpr bool
operator==(const MatrixShape<Rows1, Cols1>& a,
           const MatrixShape<Rows2, Cols2>& b) noexcept {
  return a.num_rows == b.num_rows && a.num_cols == b.num_cols;
}
template<matrix_extent Rows1, matrix_extent Cols1, //
         matrix_extent Rows2, matrix_extent Cols2>
[[nodiscard]] constexpr bool
operator!=(const MatrixShape<Rows1, Cols1>& a,
           const MatrixShape<Rows2, Cols2>& b) noexcept {
  return a.num_rows != b.num_rows || a.num_cols != b.num_cols;
}
/// @}

/// @brief Check if type is a matrix shape.
/// @{
template<class>
inline constexpr bool is_matrix_shape_v = false;
template<matrix_extent Rows, matrix_extent Cols>
inline constexpr bool is_matrix_shape_v<MatrixShape<Rows, Cols>> = true;
/// @}

/// @brief Matrix shape type.
template<class MatrixShape>
concept matrix_shape = is_matrix_shape_v<std::remove_cvref_t<MatrixShape>>;

/// @brief Matrix: has shape and two subscripts.
template<class Matrix>
concept matrix = //
    requires(Matrix& mat) {
      { mat.shape() } noexcept -> matrix_shape;
      { mat(std::declval<size_t>(), std::declval<size_t>()) } noexcept;
    };

/// @brief Number of the matrix rows.
constexpr auto num_rows(const matrix auto& mat) noexcept {
  return mat.shape().num_rows;
}

/// @brief Number of the matrix columns.
constexpr auto num_cols(const matrix auto& mat) noexcept {
  return mat.shape().num_cols;
}

/// @brief Matrix element type, as declared.
template<matrix Matrix>
using matrix_element_decltype_t = decltype( //
    std::declval<Matrix>()(std::declval<size_t>(), std::declval<size_t>()));

/// @brief Matrix element type.
template<matrix Matrix>
using matrix_element_t = std::remove_cvref_t<matrix_element_decltype_t<Matrix>>;

/// @brief Matrix element reference type.
template<matrix Matrix>
  requires std::is_lvalue_reference_v<matrix_element_decltype_t<Matrix>>
using matrix_element_ref_t = matrix_element_decltype_t<Matrix>;
/// @}

// ========================================================================== //
// ========================================================================== //

#if 0
constexpr auto eval_vectorized(auto func, //
                               matrix auto&& mat_lhs,
                               matrix auto&&... mats_rhs) noexcept {
  // When an expression is vectorizable?
  // 1. All the subexpressions are contiguous.

  auto vectorize_func = [](const auto& mat) {
    if constexpr (std::same_as<StormArray<real_t>,
                               std::decay_t<decltype(mat)>>) {
      return make_matrixView(
          num_rows(mat), num_cols(mat),
          [data = mat.data()](size_t row_index,
                              size_t col_index) -> decltype(auto) {
            STORM_ASSERT_(row_index % 4 == 0);
            return reinterpret_cast<SimdBlock<double, 4>&>(data[row_index]);
          });
    } else {
      return mat;
    }
  };

  const size_t num_blocks{mat_lhs.num_rows() / 4};
  [&](auto&& mat_lhs_, auto&&... mats_rhs_) {
#pragma omp parallel for schedule(static)
    for (size_t block_row_index = 0; block_row_index < num_blocks;
         ++block_row_index) {
      const size_t row_index{block_row_index * 4};
      func(mat_lhs_(row_index, 0), mats_rhs_(row_index, 0)...);
    }
  }(forward_as_matrix_view(mat_lhs).transform_tree(vectorize_func),
    forward_as_matrix_view(mats_rhs).transform_tree(vectorize_func)...);
  for (size_t row_index = num_blocks * 4; row_index < mat_lhs.num_rows();
       ++row_index) {
    func(mat_lhs(row_index, 0), mats_rhs(row_index, 0)...);
  }
}
#endif

constexpr auto& eval(auto func, //
                     matrix auto&& mat_lhs,
                     matrix auto&&... mats_rhs) noexcept {
  // constexpr bool verctorizable =
  //     false &&
  //     std::same_as<StormArray<real_t>, std::decay_t<decltype(mat_lhs)>>;
  // if constexpr (verctorizable) {
  // } else {
  //
#pragma omp parallel for schedule(static)
  for (size_t row_index = 0; row_index < num_rows(mat_lhs); ++row_index) {
    for (size_t col_index{0}; col_index < num_cols(mat_lhs); ++col_index) {
      func(mat_lhs(row_index, col_index), mats_rhs(row_index, col_index)...);
    }
  }
  //}

  return mat_lhs;
}

//
// To be carefully reimplemented:
//

constexpr real_t dot_product(real_t v1, real_t v2) {
  return v1 * v2;
}

constexpr real_t dot_product(matrix auto&& mat1, matrix auto&& mat2) {
  real_t d{0.0};
#if _OPENMP
  std::vector<real_t> partial(omp_get_max_threads(), 0.0);
#pragma omp parallel for schedule(static)
  for (size_t row_index = 0; row_index < num_rows(mat1); ++row_index) {
    for (size_t col_index{0}; col_index < num_cols(mat1); ++col_index) {
      partial[omp_get_thread_num()] +=
          dot_product(mat1(row_index, col_index), mat2(row_index, col_index));
    }
  }
  for (auto p : partial) {
    d += p;
  }
#else
  for (size_t row_index = 0; row_index < num_rows(mat1); ++row_index) {
    for (size_t col_index = 0; col_index < num_cols(mat1); ++col_index) {
      d += dot_product(mat1(row_index, col_index), mat2(row_index, col_index));
    }
  }
#endif
  return d;
}

constexpr real_t norm_2(matrix auto&& mat1) {
  return std::sqrt(dot_product(mat1, mat1));
}

} // namespace Storm
