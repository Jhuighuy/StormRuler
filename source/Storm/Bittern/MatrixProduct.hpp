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
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <Storm/Base.hpp>

#include <Storm/Bittern/Math.hpp>
#include <Storm/Bittern/Matrix.hpp>
#include <Storm/Bittern/MatrixView.hpp>

#include <utility>

namespace Storm
{

// -----------------------------------------------------------------------------

/// @brief Matrix product view.
template<matrix_view Matrix1, matrix_view Matrix2>
  requires numeric_matrix<Matrix1> && numeric_matrix<Matrix2>
class MatrixProductView final :
    public MatrixViewInterface<MatrixProductView<Matrix1, Matrix2>>
{
private:

  STORM_NO_UNIQUE_ADDRESS_ Matrix1 mat1_;
  STORM_NO_UNIQUE_ADDRESS_ Matrix2 mat2_;

public:

  /// @brief Construct a matrix product view.
  constexpr explicit MatrixProductView(Matrix1 mat1, Matrix2 mat2)
      : mat1_{std::move(mat1)}, mat2_{std::move(mat2)}
  {
    STORM_ASSERT_(num_cols(mat1_) == num_rows(mat2),
                  "Shapes of the matrix product arguments are invalid!");
  }

  /// @brief Get the matrix shape.
  constexpr auto shape() const noexcept
  {
    return std::array{num_rows(mat1_), num_cols(mat2_)};
  }

  /// @brief Get the matrix element at @p indices.
  constexpr auto operator()(size_t row_index, size_t col_index) const noexcept
  {
    STORM_ASSERT_(in_range(shape(), row_index, col_index),
                  "Indices are out of range!");
    // This is a default very slow generic implementation.
    auto val = mat1_(row_index, 0) * mat2_(0, col_index);
    const size_t cross_size = num_cols(mat1_);
    for (size_t cross_index = 1; cross_index < cross_size; ++cross_index) {
      val += mat1_(row_index, cross_index) * mat2_(cross_index, col_index);
    }
    return val;
  }

}; // class MatrixProductView

template<class Matrix1, class Matrix2>
MatrixProductView(Matrix1&&, Matrix2&&)
    -> MatrixProductView<matrix_view_t<Matrix1>, matrix_view_t<Matrix2>>;

/// @brief Multiply the matrices @p mat1 and @p mat2.
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
  requires numeric_matrix<Matrix1> && numeric_matrix<Matrix2>
constexpr auto matmul(Matrix1&& mat1, Matrix2&& mat2)
{
  return MatrixProductView(std::forward<Matrix1>(mat1),
                           std::forward<Matrix2>(mat2));
}

// -----------------------------------------------------------------------------

/// @brief Cross product of the 3x1 matrices view.
template<matrix_view Matrix1, matrix_view Matrix2>
  requires numeric_matrix<Matrix1> && numeric_matrix<Matrix2>
class CrossProductView : MatrixViewInterface<CrossProductView<Matrix1, Matrix2>>
{
private:

  STORM_NO_UNIQUE_ADDRESS_ Matrix1 mat1_;
  STORM_NO_UNIQUE_ADDRESS_ Matrix2 mat2_;

public:

  /// @brief Construct a matrix product view.
  constexpr explicit CrossProductView(Matrix1 mat1, Matrix2 mat2)
      : mat1_{std::move(mat1)}, mat2_{std::move(mat2)}
  {
    STORM_ASSERT_((mat1_.shape() == shape() && mat2_.shape() == shape()),
                  "Matrices of shape 3x1 are expected!");
  }

  /// @brief Get the matrix shape.
  constexpr static auto shape() noexcept
  {
    return std::array{3_sz, 1_sz};
  }

  /// @brief Get the matrix element at @p indices.
  constexpr auto operator()(size_t row_index, size_t col_index) const noexcept
  {
    STORM_ASSERT_(in_range(shape(), row_index, col_index),
                  "Indices are out of range!");
    switch (row_index) {
      case 0: return mat1_(1, 0) * mat2_(2, 0) - mat1_(2, 0) * mat2_(1, 0);
      case 1: return mat1_(2, 0) * mat2_(0, 0) - mat1_(0, 0) * mat2_(2, 0);
      case 2: return mat1_(0, 0) * mat2_(1, 0) - mat1_(1, 0) * mat2_(0, 0);
      default: STORM_UNREACHABLE_();
    }
  }

}; // class CrossProductView

template<class Matrix1, class Matrix2>
CrossProductView(Matrix1&&, Matrix2&&)
    -> CrossProductView<matrix_view_t<Matrix1>, matrix_view_t<Matrix2>>;

/// @brief Cross product of the 3x1 matrices @p mat1 and @p mat2.
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
  requires numeric_matrix<Matrix1> && numeric_matrix<Matrix2>
constexpr auto cross_product(Matrix1&& mat1, Matrix2&& mat2)
{
  return CrossProductView(std::forward<Matrix1>(mat1),
                          std::forward<Matrix2>(mat2));
}

// -----------------------------------------------------------------------------

} // namespace Storm
