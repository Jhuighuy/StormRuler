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

#include <Storm/Base.hpp>

#include <Storm/Blass/MatrixView.hpp>

#include <array>
#include <initializer_list>
#include <memory>
#include <type_traits>

namespace Storm {

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Dense matrix.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value, //
         size_t NumRows = std::dynamic_extent,
         size_t NumCols = std::dynamic_extent>
class DenseMatrix {
private:

  STORM_NO_UNIQUE_ADDRESS_
  std::conditional_t<NumRows == std::dynamic_extent, size_t,
                     std::integral_constant<size_t, NumRows>>
      num_rows_{};
  STORM_NO_UNIQUE_ADDRESS_
  std::conditional_t<NumCols == std::dynamic_extent, size_t,
                     std::integral_constant<size_t, NumCols>>
      num_cols_{};
  STORM_NO_UNIQUE_ADDRESS_
  std::conditional_t<NumRows == std::dynamic_extent || //
                         NumCols == std::dynamic_extent,
                     std::unique_ptr<Value[]>,
                     std::array<Value, NumRows * NumCols>>
      coeffs_{};

public:

  /// @brief Construct a dense matrix.
  /// @{
  constexpr DenseMatrix() //
                          /*requires(NumRows != std::dynamic_extent &&
                                   NumCols != std::dynamic_extent)*/
      = default;
  constexpr explicit DenseMatrix(size_t num_rows) //
      requires(NumRows == std::dynamic_extent &&  //
               NumCols != std::dynamic_extent)
      : num_rows_{num_rows}, //
        coeffs_{std::make_unique<Value[]>(num_rows_ * num_cols_)} {}
  constexpr explicit DenseMatrix(size_t num_cols) //
      requires(NumRows != std::dynamic_extent &&  //
               NumCols == std::dynamic_extent)
      : num_cols_{num_cols}, //
        coeffs_{std::make_unique<Value[]>(num_rows_ * num_cols_)} {}
  constexpr DenseMatrix(size_t num_rows, size_t num_cols) //
      requires(NumRows == std::dynamic_extent &&          //
               NumCols == std::dynamic_extent)
      : num_rows_{num_rows}, num_cols_{num_cols}, //
        coeffs_{std::make_unique<Value[]>(num_rows_ * num_cols_)} {}
  /// @}

  /// @brief Construct a dense matrix with coefficients.
  constexpr DenseMatrix(
      std::initializer_list<std::initializer_list<Value>> coeffs) //
      noexcept(NumRows != std::dynamic_extent &&
               NumCols != std::dynamic_extent) {
    if constexpr (NumRows == std::dynamic_extent) {
      num_rows_ = coeffs.size();
    } else {
      STORM_ASSERT_(num_rows_ == coeffs.size() && //
                    "Invalid number of rows.");
    }
    if constexpr (NumCols == std::dynamic_extent) {
      num_cols_ = coeffs.begin()->size();
    } else {
      STORM_ASSERT_(num_rows_ == coeffs.begin()->size() &&
                    "Invalid number of columns.");
    }
    if constexpr (NumRows == std::dynamic_extent ||
                  NumCols == std::dynamic_extent) {
      coeffs_ = std::make_unique<Value[]>(num_rows_ * num_cols_);
    }
    for (size_t i{0}; i < num_rows_; ++i) {
      for (size_t j{0}; j < num_cols_; ++j) {
        (*this)(i, j) = (coeffs.begin()[i]).begin()[j];
      }
    }
  }

  constexpr void assign(size_t num_rows)         //
      requires(NumRows == std::dynamic_extent && //
               NumCols != std::dynamic_extent) {
    num_rows_ = num_rows;
    coeffs_ = std::make_unique<Value[]>(num_rows_ * num_cols_);
  }
  constexpr void assign(size_t num_cols)         //
      requires(NumRows != std::dynamic_extent && //
               NumCols == std::dynamic_extent) {
    num_cols_ = num_cols;
    coeffs_ = std::make_unique<Value[]>(num_rows_ * num_cols_);
  }
  constexpr void assign(size_t num_rows,
                        size_t num_cols)         //
      requires(NumRows == std::dynamic_extent && //
               NumCols == std::dynamic_extent) {
    num_rows_ = num_rows, num_cols_ = num_cols;
    coeffs_ = std::make_unique<Value[]>(num_rows_ * num_cols_);
  }

  /// @brief Number of the matrix rows.
  constexpr auto num_rows() const noexcept {
    return num_rows_;
  }

  /// @brief Number of the matrix columns.
  constexpr auto num_cols() const noexcept {
    return num_cols_;
  }
  auto shape() const noexcept {
    return MatrixShape(num_rows(), num_cols());
  }

  /// @brief Get the matrix coefficient at @p row_index and @p col_index.
  /// @{
  constexpr auto& operator()(size_t row_index, //
                             size_t col_index = 0) noexcept {
    STORM_ASSERT_(row_index < num_rows_ && col_index < num_cols_ &&
                  "Matrix index out range.");
    return coeffs_[row_index * num_cols_ + col_index];
  }
  constexpr const auto& operator()(size_t row_index,
                                   size_t col_index = 0) const noexcept {
    return const_cast<DenseMatrix&>(*this)(row_index, col_index);
  }
  /// @}

}; // class DenseMatrix

template<class Value, size_t NumRows = std::dynamic_extent>
using DenseVector = DenseMatrix<Value, NumRows, 1>;

/// @brief Perform a LU decomposition of a square matrix @p mat.
constexpr void decompose_lu(matrix auto&& mat, //
                            matrix auto&& l_mat, matrix auto&& u_mat) noexcept {
  const auto size{num_rows(mat)};
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

constexpr void inplace_solve_lu(matrix auto&& l_mat, matrix auto&& u_mat,
                                auto& vec) {
  const auto size{num_rows(l_mat)};
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
constexpr void inplace_inverse_lu(matrix auto&& mat,
                                  matrix auto&& inv_mat) noexcept {
  using Value = std::decay_t<decltype(mat(0, 0))>;
  DenseMatrix<Value> L(num_rows(mat), num_cols(mat));
  DenseMatrix<Value> U(num_rows(mat), num_cols(mat));
  decompose_lu(mat, L, U);
  fill_diag_with(inv_mat, 1.0);
  for (size_t iy{0}; iy < num_rows(mat); ++iy) {
    auto inv_mat_col = [&](size_t ix) -> Value& { return inv_mat(ix, iy); };
    inplace_solve_lu(L, U, inv_mat_col);
  }
}

/// @brief Perform a QR decomposition of a matrix @p mat.
/// @returns A pair of matrices, Q and R factors.
constexpr auto DecomposeQr(auto& mat) noexcept {}

} // namespace Storm
