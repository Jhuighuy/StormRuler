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

#include <array>
#include <memory>
#include <span>
#include <type_traits>

#include <stormBase.hxx>
#include <stormSolvers/MatrixView.hxx>

namespace Storm {

template<class Value, size_t NumRows = std::dynamic_extent,
         size_t NumCols = std::dynamic_extent>
class DenseMatrix {
private:

  std::conditional_t<
      NumRows == std::dynamic_extent || NumCols == std::dynamic_extent,
      std::unique_ptr<Value[]>, std::array<Value, NumRows * NumCols>>
      coeffs_{};
  [[no_unique_address]] std::conditional_t<
      NumRows == std::dynamic_extent, size_t,
      std::integral_constant<size_t, NumRows>>
      num_rows_{};
  [[no_unique_address]] std::conditional_t<
      NumCols == std::dynamic_extent, size_t,
      std::integral_constant<size_t, NumCols>>
      num_cols_{};

public:

  /// @brief Initialize an empty matrix.
  constexpr DenseMatrix() = default;

  constexpr DenseMatrix(size_t num_rows, size_t num_cols) noexcept
      : num_rows_{num_rows}, num_cols_{num_cols} {
    coeffs_ = std::make_unique<Value[]>(size());
  }

  constexpr DenseMatrix(
      std::initializer_list<std::initializer_list<Value>> coeffs) noexcept
      : DenseMatrix(coeffs.size(), coeffs.begin()->size()) {
    *this <<= MatrixExpr(num_rows_, num_cols_, [&](auto i, auto j) {
      return *((coeffs.begin() + i)->begin() + j);
    });
  }

  constexpr auto num_rows() const noexcept {
    return num_rows_;
  }

  constexpr auto num_cols() const noexcept {
    return num_cols_;
  }

  constexpr auto size() const noexcept {
    return num_rows_ * num_cols_;
  }

  /// @brief Get the coefficient at @p row_index and @p col_index.
  /// @{
  constexpr auto& operator()(size_t row_index, size_t col_index) noexcept {
    STORM_ASSERT_(row_index < num_rows_ && col_index < num_cols_ &&
                  "Matrix index out range.");
    return coeffs_[row_index * num_cols_ + col_index];
  }
  constexpr const auto& operator()(size_t row_index,
                                   size_t col_index) const noexcept {
    return const_cast<DenseMatrix&>(*this)(row_index, col_index);
  }
  /// @}

}; // class Matrix

template<class Value, size_t NumRows, size_t NumCols>
struct is_matrix_t<DenseMatrix<Value, NumRows, NumCols>> : std::true_type {};

/// @brief Perform a LU decomposition of a square matrix @p mat.
constexpr void decompose_lu(const is_matrix_view auto& mat,
                            is_rw_matrix_view auto& l_mat,
                            is_rw_matrix_view auto& u_mat) noexcept {
  const auto size{mat.num_rows()};
  fill_diag_with(l_mat, 1.0);
  fill_with(u_mat, 0.0);
  for (size_t ix{0}; ix < size; ++ix) {
    for (size_t iy{0}; iy < ix; ++iy) {
      l_mat(ix, iy) = mat(ix, iy);
      for (size_t iz{0}; iz < iy; ++iz) {
        l_mat(ix, iy) -= l_mat(ix, iz) * u_mat(iz, iy);
      }
      l_mat(ix, iy) /= u_mat(iy, iy);
    }
    for (size_t iy{ix}; iy < size; ++iy) {
      u_mat(ix, iy) = mat(ix, iy);
      for (size_t iz{0}; iz < ix; ++iz) {
        u_mat(ix, iy) -= l_mat(ix, iz) * u_mat(iz, iy);
      }
    }
  }
}

constexpr void inplace_solve_lu(const is_matrix_view auto& l_mat,
                                const is_matrix_view auto& u_mat, auto& vec) {
  const auto size{l_mat.num_rows()};
  for (size_t ix{0}; ix < size; ++ix) {
    for (size_t iy{0}; iy < ix; ++iy) {
      vec(ix) -= l_mat(ix, iy) * vec(iy);
    }
    vec(ix) /= l_mat(ix, ix);
  }
  for (size_t rix{0}; rix < size; ++rix) {
    size_t ix{size - 1 - rix};
    for (size_t iy{ix + 1}; iy < size; ++iy) {
      vec(ix) -= u_mat(ix, iy) * vec(iy);
    }
    vec(ix) /= u_mat(ix, ix);
  }
}

/// @brief Inverse a square matrix @p mat using the LU decomposition.
constexpr void inplace_inverse_lu(const is_matrix_view auto& mat,
                                  is_rw_matrix_view auto& inv_mat) noexcept {
  using Value = std::decay_t<decltype(mat(0, 0))>;
  DenseMatrix<Value> L(mat.num_rows(), mat.num_cols());
  DenseMatrix<Value> U(mat.num_rows(), mat.num_cols());
  decompose_lu(mat, L, U);
  fill_diag_with(inv_mat, 1.0);
  for (size_t iy{0}; iy < mat.num_rows(); ++iy) {
    auto inv_mat_col = [&](size_t ix) -> Value& { return inv_mat(ix, iy); };
    inplace_solve_lu(L, U, inv_mat_col);
  }
}

/// @brief Perform a QR decomposition of a matrix @p mat.
/// @returns A pair of matrices, Q and R factors.
constexpr auto DecomposeQr(auto& mat) noexcept {}

} // namespace Storm
