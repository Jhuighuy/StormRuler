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

#include <stormBase.hxx>

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

// clang-format off

/// @brief Matrix view-like concept.
template<class MatrixView>
concept is_matrix_view_like = 
    requires(const MatrixView& mat) {
      { mat.num_rows() } -> std::convertible_to<size_t>;
      { mat.num_cols() } -> std::convertible_to<size_t>;
    } && requires(const MatrixView& mat, size_t row_index, size_t col_index) {
      { mat(row_index, col_index) };
    };

/// @brief Read-write matrix view-like concept.
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

constexpr auto& eval(auto func, decays_to_rw_matrix_view auto&& mat_lhs,
                     const is_matrix_view auto&... mats_rhs) {
  if (mat_lhs.num_rows() * mat_lhs.num_cols() > 1000) {
#pragma omp parallel for schedule(static)
    for (int row_index = 0; row_index < (int) mat_lhs.num_rows(); ++row_index) {
      for (size_t col_index{0}; col_index < mat_lhs.num_cols(); ++col_index) {
        func(mat_lhs(row_index, col_index), mats_rhs(row_index, col_index)...);
      }
    }
  } else {
    for (size_t row_index{0}; row_index < mat_lhs.num_rows(); ++row_index) {
      for (size_t col_index{0}; col_index < mat_lhs.num_cols(); ++col_index) {
        func(mat_lhs(row_index, col_index), mats_rhs(row_index, col_index)...);
      }
    }
  }
  return mat_lhs;
}

//
// To be carefully reimplemented:
//

constexpr real_t dot_product(const is_matrix_view auto& mat1,
                             const is_matrix_view auto& mat2) {
  real_t d = 0.0;
  if (mat1.num_rows() * mat1.num_cols() > 1000) {
#pragma omp parallel for reduction(+ : d) schedule(static)
    for (int row_index = 0; row_index < (int) mat1.num_rows(); ++row_index) {
      for (size_t col_index{0}; col_index < mat1.num_cols(); ++col_index) {
        d += mat1(row_index, col_index) * mat2(row_index, col_index);
      }
    }
  } else {
    for (size_t row_index{0}; row_index < mat1.num_rows(); ++row_index) {
      for (size_t col_index{0}; col_index < mat1.num_cols(); ++col_index) {
        d += mat1(row_index, col_index) * mat2(row_index, col_index);
      }
    }
  }
  return d;
}

constexpr real_t norm_2(const is_matrix_view auto& mat1) {
  return std::sqrt(dot_product(mat1, mat1));
}

} // namespace Storm
