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

#include <concepts>
#include <type_traits>

#include <glm/glm.hpp>
#include <omp.h>

#include <stormBase.hxx>
#include <stormUtils/SimdBlock.hxx>

namespace Storm {

/// @brief Check if type is a matrix.
/// @{
template<class>
struct is_matrix_t : std::false_type {};
template<class T>
inline constexpr bool is_matrix_v = is_matrix_t<T>::value;
/// @}

/// @brief Check if type is a matrix view.
/// @{
template<class T>
struct is_matrix_view_t : is_matrix_t<T> {};
template<class T>
inline constexpr bool is_matrix_view_v = is_matrix_view_t<T>::value;
/// @}

template<class T>
using matrix_coeff_t = std::remove_cv_t<decltype(std::declval<T>()(
    std::declval<size_t>(), std::declval<size_t>()))>;

/// @brief Matrix view-like concept.
// clang-format off
template<class MatrixView>
concept is_matrix_view_like = 
    requires(const MatrixView& mat) {
      { mat.num_rows() } -> std::convertible_to<size_t>;
      { mat.num_cols() } -> std::convertible_to<size_t>;
    } && requires(const MatrixView& mat, size_t row_index, size_t col_index) {
      { mat(row_index, col_index) };
    };
// clang-format on

/// @brief Read-write matrix view-like concept.
// clang-format off
template<class Matrix>
concept is_rw_matrix_view_like = 
    is_matrix_view_like<Matrix> && 
    requires(Matrix& mat, size_t row_index, size_t col_index, 
             matrix_coeff_t<Matrix> val) {
      { mat(row_index, col_index) = val };
    };
// clang-format on

/// @brief Matrix concept.
template<class T>
concept is_matrix = is_matrix_v<T> && is_rw_matrix_view_like<T>;

/// @brief Read-write matrix view concept.
template<class T>
concept is_rw_matrix_view = is_matrix_view_v<T> && is_rw_matrix_view_like<T>;

/// @brief Matrix view concept.
template<class T>
concept is_matrix_view = is_matrix_view_v<T> && is_matrix_view_like<T>;

/// @brief Decays to a matrix concept.
template<class T>
concept decays_to_matrix = is_matrix<std::remove_cvref_t<T>>;

/// @brief Decays to a read-write matrix view concept.
template<class T>
concept decays_to_rw_matrix_view =
    is_rw_matrix_view<std::remove_reference_t<T>>;

/// @brief Decays to a matrix view concept.
template<class T>
concept decays_to_matrix_view = is_matrix_view<std::remove_cvref_t<T>>;

/// @brief Matrix view that is not a matrix concept.
template<class T>
concept is_strictly_matrix_view = is_matrix_view<T> && !is_matrix<T>;

constexpr auto& eval(auto func, decays_to_matrix auto&& mat_lhs,
                     const is_matrix_view auto&... mats_rhs) noexcept {
#if 0

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
  for (size_t row_index = 0; row_index < mat_lhs.num_rows(); ++row_index) {
    for (size_t col_index{0}; col_index < mat_lhs.num_cols(); ++col_index) {
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
real_t dot_product(glm::dvec2 v1, glm::dvec2 v2) {
  return glm::dot(v1, v2);
}

constexpr real_t dot_product(const is_matrix_view auto& mat1,
                             const is_matrix_view auto& mat2) {
  std::vector<real_t> partial(omp_get_max_threads(), 0.0);
#pragma omp parallel for schedule(static)
  for (size_t row_index = 0; row_index < mat1.num_rows(); ++row_index) {
    for (size_t col_index{0}; col_index < mat1.num_cols(); ++col_index) {
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

constexpr real_t norm_2(const is_matrix_view auto& mat1) {
  return std::sqrt(dot_product(mat1, mat1));
}

} // namespace Storm
