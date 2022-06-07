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

#include <stdlib.h>

#include <concepts>
#include <type_traits>

#include <stormBase.hxx>

namespace Storm {

/// @brief Check if type is a matrix view.
/// @{
template<class>
struct is_matrix_view_t : std::false_type {};
template<class T>
inline constexpr bool is_matrix_view_v = is_matrix_view_t<T>::value;
/// @}

template<class T>
using matrix_coeff_t = std::remove_cv_t<decltype(std::declval<T>()(
    std::declval<size_t>(), std::declval<size_t>()))>;

// clang-format off

/// @brief Matrix view concept.
template<class MatrixView>
concept is_matrix_view_like = 
    requires(const MatrixView& mat) {
      { mat.num_rows() } -> std::convertible_to<size_t>;
      { mat.num_cols() } -> std::convertible_to<size_t>;
    } && requires(const MatrixView& mat, size_t row_index, size_t col_index) {
      { mat(row_index, col_index) };
    };

/// @brief Matrix concept.
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
concept is_matrix_view = is_matrix_view_v<T> && is_matrix_view_like<T>;

/// @brief Matrix concept.
template<class T>
concept is_rw_matrix_view = is_matrix_view_v<T> && is_rw_matrix_view_like<T>;

constexpr real_t dot_product(const is_matrix_view auto& mat1,
                             const is_matrix_view auto& mat2) {
  real_t d{};
#pragma omp parallel for reduction(+ : d) schedule(static)
  for (int row_index = 0; row_index < (int) mat1.num_rows(); ++row_index) {
    for (size_t col_index{0}; col_index < mat1.num_cols(); ++col_index) {
      d += mat1(row_index, col_index) * mat2(row_index, col_index);
    }
  }
  return d;
}

constexpr real_t norm_2(const is_matrix_view auto& mat1) {
  return std::sqrt(dot_product(mat1, mat1));
}

constexpr auto& operator<<=(is_rw_matrix_view auto& mat1,
                            const is_matrix_view auto& mat2) {
#pragma omp parallel for schedule(static)
  for (int row_index = 0; row_index < (int) mat1.num_rows(); ++row_index) {
    for (size_t col_index{0}; col_index < mat1.num_cols(); ++col_index) {
      mat1(row_index, col_index) = mat2(row_index, col_index);
    }
  }
  return mat1;
}

constexpr auto& operator+=(is_rw_matrix_view auto& mat1,
                           const is_matrix_view auto& mat2) {
  return mat1 <<= mat1 + mat2;
}
constexpr auto& operator-=(is_rw_matrix_view auto& mat1,
                           const is_matrix_view auto& mat2) {
  return mat1 <<= mat1 - mat2;
}

constexpr auto& fill_with(is_rw_matrix_view auto& mat1, const auto& val2) {
#pragma omp parallel for schedule(static)
  for (int row_index = 0; row_index < (int) mat1.num_rows(); ++row_index) {
    for (size_t col_index{0}; col_index < mat1.num_cols(); ++col_index) {
      mat1(row_index, col_index) = val2;
    }
  }
  return mat1;
}

constexpr auto& fill_diag_with(is_rw_matrix_view auto& mat1, const auto& val2) {
#pragma omp parallel for schedule(static)
  for (int row_index = 0; row_index < (int) mat1.num_rows(); ++row_index) {
    for (size_t col_index{0}; col_index < mat1.num_cols(); ++col_index) {
      mat1(row_index, col_index) = (size_t) row_index == col_index ? val2 : 0;
    }
  }
  return mat1;
}

constexpr auto& fill_randomly(is_rw_matrix_view auto& mat1, auto...) {
  srand(1792);
  for (int row_index = 0; row_index < (int) mat1.num_rows(); ++row_index) {
    for (size_t col_index{0}; col_index < mat1.num_cols(); ++col_index) {
      mat1(row_index, col_index) = real_t(rand()) / RAND_MAX;
    }
  }
  return mat1;
}

constexpr auto& operator*=(is_rw_matrix_view auto& mat1, const auto& val2) {
  return mat1 <<= val2 * mat1;
}
constexpr auto& operator/=(is_rw_matrix_view auto& mat1, const auto& val2) {
  return mat1 <<= mat1 / val2;
}

} // namespace Storm
