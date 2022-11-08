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

#pragma once

#include <Storm/Base.hpp>

#include <algorithm>
#include <concepts>
#include <limits>
#include <random>
#include <type_traits>

namespace Storm {

/// @brief Matrix shape.
template<std::convertible_to<size_t> NumRows = size_t,
         std::convertible_to<size_t> NumCols = size_t>
class MatrixShape final {
public:

  /// @brief Number of the matrix rows.
  STORM_NO_UNIQUE_ADDRESS_ NumRows num_rows{};

  /// @brief Number of the matrix columns.
  STORM_NO_UNIQUE_ADDRESS_ NumCols num_cols{};

  /// @brief Compare the matrix shapes.
  constexpr auto operator<=>(const MatrixShape&) const = default;

}; // class MatrixShape

template<class>
inline constexpr bool is_matrix_shape_v = false;
template<class NumRows, class NumCols>
inline constexpr bool is_matrix_shape_v<MatrixShape<NumRows, NumCols>> = true;

template<class MatrixShape>
concept matrix_shape = is_matrix_shape_v<MatrixShape>;

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
template<matrix OutMatrix>
  requires std::floating_point<matrix_element_t<OutMatrix>>
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
template<matrix OutMatrix, std::copyable Scalar>
  requires (!matrix<Scalar>)
constexpr OutMatrix& operator*=(OutMatrix& out_mat, Scalar scal) {
  return assign(out_mat, [scal = std::move(scal)](auto& out_elem) noexcept {
    out_elem *= scal;
  });
}

/// @brief Divide-assign the matrix @p out_mat by a scalar @p scal.
template<matrix OutMatrix, std::copyable Scalar>
  requires (!matrix<Scalar>)
constexpr OutMatrix& operator/=(OutMatrix& out_mat, Scalar scal) {
  return assign(out_mat, [scal = std::move(scal)](auto& out_elem) noexcept {
    out_elem /= scal;
  });
}

/// @brief Add-assign the matrices @p out_mat and @p mat
template<matrix OutMatrix, matrix Matrix>
constexpr OutMatrix& operator+=(OutMatrix& out_mat, Matrix&& mat) {
  return assign(
      out_mat,
      []<class Elem>(auto& out_elem, Elem&& elem) noexcept {
        out_elem += std::forward<Elem>(elem);
      },
      std::forward<Matrix>(mat));
}

/// @brief Subtract-assign the matrices @p out_mat and @p mat.
template<matrix OutMatrix, matrix Matrix>
constexpr OutMatrix& operator-=(OutMatrix& out_mat, Matrix&& mat) {
  return assign(
      out_mat,
      []<class Elem>(auto& out_elem, Elem&& elem) noexcept {
        out_elem -= std::forward<Elem>(elem);
      },
      std::forward<Matrix>(mat));
}

/// @brief Element-wise multiply-assign the matrices @p out_mat and @p mat.
template<matrix OutMatrix, matrix Matrix>
constexpr OutMatrix& operator*=(OutMatrix& out_mat, Matrix&& mat) {
  return assign(
      out_mat,
      []<class Elem>(auto& out_elem, Elem&& elem) noexcept {
        out_elem *= std::forward<Elem>(elem);
      },
      std::forward<Matrix>(mat));
}

/// @brief Element-wise divide-assign the matrices @p out_mat and @p mat.
template<matrix OutMatrix, matrix Matrix>
constexpr OutMatrix& operator/=(OutMatrix& out_mat, Matrix&& mat) {
  return assign(
      out_mat,
      []<class Elem>(auto& out_elem, Elem&& elem) noexcept {
        out_elem /= std::forward<Elem>(elem);
      },
      std::forward<Matrix>(mat));
}

/// @brief Reduce the matrix @p mat coefficients to a single value.
/// @param init Initial reduction value.
/// @param reduce_func Reduction function.
/// @{
template<class Value, class ReduceFunc, matrix Matrix>
[[nodiscard]] constexpr auto reduce(Value init, ReduceFunc reduce_func,
                                    Matrix&& mat) {
  for (size_t row_index = 0; row_index < num_rows(mat); ++row_index) {
    for (size_t col_index = 0; col_index < num_cols(mat); ++col_index) {
      init = reduce_func(init, mat(row_index, col_index));
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
      init = reduce_func(init, func(mat(row_index, col_index), //
                                    mats(row_index, col_index)...));
    }
  }
  return init;
}
/// @}

/// @brief Sum the matrix @p mat elements.
template<matrix Matrix>
[[nodiscard]] constexpr auto sum(Matrix&& mat) {
  return reduce(matrix_element_t<Matrix>{0.0}, std::plus{},
                std::forward<Matrix>(mat));
}

/// @brief Minimum matrix @p mat element.
template<matrix Matrix>
[[nodiscard]] constexpr auto min_element(Matrix&& mat) {
  using Elem = matrix_element_t<Matrix>;
  return reduce(
      std::numeric_limits<Elem>::max(),
      [](const Elem& a, const Elem& b) noexcept { return min(a, b); },
      std::forward<Matrix>(mat));
}
/// @brief Maximum matrix @p mat element.
template<matrix Matrix>
[[nodiscard]] constexpr auto max_element(Matrix&& mat) {
  using Elem = matrix_element_t<Matrix>;
  return reduce(
      std::numeric_limits<Elem>::lowest(),
      [](const Elem& a, const Elem& b) noexcept { return max(a, b); },
      std::forward<Matrix>(mat));
}

/// @brief Element-wise matrix @p mat \f$ L_{1} \f$-norm.
template<matrix Matrix>
[[nodiscard]] constexpr auto norm_1(Matrix&& mat) {
  return reduce(
      matrix_element_t<Matrix>{0.0}, std::plus{},
      []<class Elem>(Elem&& elem) noexcept {
        return abs(std::forward<Elem>(elem));
      },
      std::forward<Matrix>(mat));
}

/// @brief Element-wise matrix @p mat \f$ L_{2} \f$-norm.
template<matrix Matrix>
[[nodiscard]] constexpr auto norm_2(Matrix&& mat) {
  const auto temp = reduce(
      matrix_element_t<Matrix>{0.0}, std::plus{},
      []<class Elem>(Elem&& elem) noexcept {
        /// @todo Possibly a better way here, like `elem * conj(elem)`.
        return abs(pow(std::forward<Elem>(elem), 2));
      },
      std::forward<Matrix>(mat));
  return sqrt(temp);
}

/// @brief Element-wise matrix @p mat \f$ L_{p} \f$-norm.
template<matrix Matrix>
[[nodiscard]] constexpr auto norm_p(Matrix&& mat, real_t p) {
  STORM_ASSERT_(p > 0.0, "Invalid p-norm parameter!");
  const auto temp = reduce(
      matrix_element_t<Matrix>{0.0}, std::plus{},
      [p]<class Elem>(Elem&& elem) noexcept {
        return pow(abs(std::forward<Elem>(elem)), p);
      },
      std::forward<Matrix>(mat));
  return pow(temp, 1.0 / p);
}

/// @brief Element-wise matrix @p mat \f$ L_{\infty} \f$-norm.
template<matrix Matrix>
[[nodiscard]] constexpr auto norm_inf(Matrix&& mat) {
  using Elem = matrix_element_t<Matrix>;
  return reduce(
      Elem{0.0},
      [](const Elem& a, const Elem& b) noexcept { return max(a, b); },
      []<class Elem>(Elem&& elem) noexcept {
        return abs(std::forward<Elem>(elem));
      },
      std::forward<Matrix>(mat));
}

/// @brief Element-wise dot product of the matrices @p mat1 and @p mat2.
template<matrix Matrix1, matrix Matrix2>
constexpr auto dot_product(Matrix1&& mat1, Matrix2&& mat2) noexcept {
  /// @todo Complex inputs!
  using Result = decltype(std::declval<matrix_element_t<Matrix1>>() *
                          std::declval<matrix_element_t<Matrix2>>());
  return reduce(Result{0.0}, std::plus{}, std::multiplies{},
                std::forward<Matrix1>(mat1), std::forward<Matrix2>(mat2));
}

} // namespace Storm
