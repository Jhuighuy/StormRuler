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

#include <Storm/Bittern/Matrix.hpp>

#include <concepts>
#include <functional>
#include <tuple>
#include <type_traits>
#include <utility>

namespace Storm {

/// @brief Matrix product view.
template<matrix_view Matrix1, matrix_view Matrix2>
class MatrixProductView final :
    public MatrixViewInterface<MatrixProductView<Matrix1, Matrix2>> {
private:

  STORM_NO_UNIQUE_ADDRESS_ Matrix1 mat1_;
  STORM_NO_UNIQUE_ADDRESS_ Matrix2 mat2_;

public:

  /// @brief Construct a matrix product view.
  constexpr explicit MatrixProductView(Matrix1 mat1, Matrix2 mat2)
      : mat1_{std::move(mat1)}, mat2_{std::move(mat2)} {
    STORM_ASSERT_(num_cols(mat1_) == num_rows(mat2),
                  "Shapes of the matrix product arguments are invalid.");
  }

  /// @copydoc MatrixViewInterface::shape
  [[nodiscard]] constexpr auto shape() const noexcept {
    return MatrixShape{num_rows(mat1_), num_cols(mat2_)};
  }

  /// @copydoc MatrixViewInterface::operator()
  [[nodiscard]] constexpr auto operator()(size_t row_index,
                                          size_t col_index) const noexcept {
    // This is a default very slow generic implementation.
    auto val{mat1_(row_index, 0) * mat2_(0, col_index)};
    const auto cross_size = static_cast<size_t>(num_cols(mat1_));
    for (size_t cross_index{1}; cross_index < cross_size; ++cross_index) {
      val += mat1_(row_index, cross_index) * mat2_(cross_index, col_index);
    }
    return val;
  }

}; // class MatrixProductView

template<class Matrix1, class Matrix2>
MatrixProductView(Matrix1&&, Matrix2&&)
    -> MatrixProductView<forward_as_matrix_view_t<Matrix1>,
                         forward_as_matrix_view_t<Matrix2>>;

/// @brief Multiply the matrices @p mat1 and @p mat2.
template<class Matrix1, class Matrix2>
[[nodiscard]] constexpr auto matmul(Matrix1&& mat1, Matrix2&& mat2) {
  return MatrixProductView(std::forward<Matrix1>(mat1),
                           std::forward<Matrix2>(mat2));
}

} // namespace Storm
