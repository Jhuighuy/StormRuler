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

// This header should really be <Matrix.hpp>
#pragma once

#include <concepts>
#include <type_traits>

#include <omp.h>
#include <ranges>

#include <Storm/Base.hpp>

#include <Storm/Blass/MatrixShape.hpp>
#include <Storm/Utils/SimdBlock.hpp>

namespace Storm {

/// @brief Matrix: has shape and two subscripts.
// clang-format off
template<class Matrix>
concept matrix = 
    requires(Matrix& mat) { { mat.shape() } noexcept -> matrix_shape; } &&
    requires(Matrix& mat, size_t row_index, size_t col_index) {
      { mat(row_index, col_index) } noexcept;
    };
// clang-format on

/// @brief Number of the matrix rows.
constexpr auto num_rows(const matrix auto& mat) noexcept {
  return mat.shape().num_rows();
}

/// @brief Number of the matrix columns.
constexpr auto num_cols(const matrix auto& mat) noexcept {
  return mat.shape().num_cols();
}

template<matrix Matrix>
using matrix_reference_t = decltype( //
    std::declval<Matrix>()(std::declval<size_t>(), std::declval<size_t>()));

template<matrix Matrix>
using const_matrix_reference_t = decltype( //
    std::declval<std::add_const_t<Matrix>>()(std::declval<size_t>(),
                                             std::declval<size_t>()));

/// @brief Contiguous matrix: has pointer to data.
// clang-format off
template<class Matrix>
concept contiguous_matrix =
    matrix<Matrix> && requires(Matrix& mat) {
      { data(mat) } -> 
          std::same_as<std::add_pointer_t<matrix_reference_t<Matrix>>>;
    };
// clang-format on

template<class Matrix>
constexpr inline bool enable_flattenable_matrix_v{false};

template<class Matrix>
concept flattenable_matrix =
    contiguous_matrix<Matrix> || enable_flattenable_matrix_v<Matrix>;

constexpr auto& eval(auto func, matrix auto&& mat_lhs,
                     matrix auto&&... mats_rhs) noexcept {
#if 0
#elif 0

  // When an expression is vectorizable?
  // 1. All the subexpressions are contiguous.

  auto vectorize_func{[](auto& mat) {
    return make_expression(
        mat.num_rows(), mat.num_cols(),
        [data = mat.data()](size_t row_index, [[maybe_unused]] size_t col_index)
            -> decltype(auto) {
          STORM_ASSERT_(row_index % 4 == 0);
          return reinterpret_cast<SimdBlock<double, 4>&>(data[row_index]);
        });
  }};

  [&](auto&& mat_lhs_, auto&&... mats_rhs_) {
    const size_t num_blocks{mat_lhs.num_rows() / 4};
#pragma omp parallel for schedule(static)
    for (size_t block_row_index = 0; block_row_index < num_blocks;
         ++block_row_index) {
      const size_t row_index{block_row_index * 4};
      func(mat_lhs_(row_index, 0), mats_rhs_(row_index, 0)...);
    }
    for (size_t row_index = num_blocks * 4; row_index < mat_lhs.num_rows();
         ++row_index) {
      func(mat_lhs(row_index, 0), mats_rhs(row_index, 0)...);
    }
  }(transform_expression_leafs(mat_lhs, vectorize_func),
    transform_expression_leafs(mats_rhs, vectorize_func)...);

#else

#pragma omp parallel for schedule(static)
  for (size_t row_index = 0; row_index < num_rows(mat_lhs); ++row_index) {
    for (size_t col_index{0}; col_index < num_cols(mat_lhs); ++col_index) {
      func(mat_lhs(row_index, col_index), mats_rhs(row_index, col_index)...);
    }
  }

#endif

  return mat_lhs;
}

//
// To be carefully reimplemented:
//

constexpr real_t dot_product(real_t v1, real_t v2) {
  return v1 * v2;
}

constexpr real_t dot_product(matrix auto&& mat1, matrix auto&& mat2) {
  std::vector<real_t> partial(omp_get_max_threads(), 0.0);
#pragma omp parallel for schedule(static)
  for (size_t row_index = 0; row_index < num_rows(mat1); ++row_index) {
    for (size_t col_index{0}; col_index < num_cols(mat1); ++col_index) {
      partial[omp_get_thread_num()] +=
          dot_product(mat1(row_index, col_index), mat2(row_index, col_index));
    }
  }
  real_t d{0.0};
  for (auto p : partial) {
    d += p;
  }
  return d;
}

constexpr real_t norm_2(matrix auto&& mat1) {
  return std::sqrt(dot_product(mat1, mat1));
}

} // namespace Storm
