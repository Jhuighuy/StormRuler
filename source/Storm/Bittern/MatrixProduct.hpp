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

#include <Storm/Crow/MathUtils.hpp>

#include <Storm/Bittern/MatrixView.hpp>
#include <Storm/Bittern/Shape.hpp>

#include <concepts>
#include <utility>

namespace Storm {

// -----------------------------------------------------------------------------

#if 0
/// @brief Matrix product view.
template<matrix_view Matrix1, matrix_view Matrix2>
class MatrixProductView final :
    public MatrixViewInterface<MatrixProductView<Matrix1, Matrix2>> {
private:

  STORM_NO_UNIQUE_ADDRESS Matrix1 mat1_;
  STORM_NO_UNIQUE_ADDRESS Matrix2 mat2_;

public:

  /// @brief Construct a matrix product view.
  constexpr explicit MatrixProductView(Matrix1 mat1, Matrix2 mat2)
      : mat1_{std::move(mat1)}, mat2_{std::move(mat2)} {
    STORM_ASSERT(num_cols(mat1_) == num_rows(mat2),
                 "Shapes of the matrix product arguments are invalid!");
  }

  /// @brief Get the matrix shape.
  constexpr auto shape() const noexcept {
    return std::array{num_rows(mat1_), num_cols(mat2_)};
  }

  /// @brief Get the matrix element at @p indices.
  constexpr auto operator()(size_t row_index, size_t col_index) const noexcept {
    STORM_ASSERT(in_range(shape(), row_index, col_index),
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
constexpr auto matmul(Matrix1&& mat1, Matrix2&& mat2) {
  return MatrixProductView(std::forward<Matrix1>(mat1),
                           std::forward<Matrix2>(mat2));
}
#endif

// -----------------------------------------------------------------------------

/// @brief 3D Vector cross product view.
template<matrix_view_r<1> Vector1, matrix_view_r<1> Vector2>
class CrossProductView final :
    public MatrixViewInterface<CrossProductView<Vector1, Vector2>> {
private:

  STORM_NO_UNIQUE_ADDRESS Vector1 vec1_;
  STORM_NO_UNIQUE_ADDRESS Vector2 vec2_;

public:

  /// @brief Construct a cross product view.
  constexpr explicit CrossProductView(Vector1 vec1, Vector2 vec2)
      : vec1_{std::move(vec1)}, vec2_{std::move(vec2)} {
    STORM_ASSERT(vec1_.shape() == shape() && vec2_.shape() == shape(),
                 "3D vectors are expected!");
  }

  /// @brief Get the vector shape.
  constexpr static auto shape() {
    return fixed_shape_t<3>{};
  }

  /// @brief Get the vector element at @p index.
  template<class Index>
    requires compatible_matrix_indices_v<CrossProductView, Index> &&
             std::convertible_to<Index, size_t>
  constexpr auto operator()(Index index) const {
    STORM_ASSERT(in_range(shape(), index), "Index is out of range!");
    switch (static_cast<size_t>(index)) {
      case 0: return vec1_(1_sz) * vec2_(2_sz) - vec1_(2_sz) * vec2_(1_sz);
      case 1: return vec1_(2_sz) * vec2_(0_sz) - vec1_(0_sz) * vec2_(2_sz);
      case 2: return vec1_(0_sz) * vec2_(1_sz) - vec1_(1_sz) * vec2_(0_sz);
      default: STORM_UNREACHABLE();
    }
  }

}; // class CrossProductView

template<class Vector1, class Vector2>
CrossProductView(Vector1&&, Vector2&&)
    -> CrossProductView<matrix_view_t<Vector1>, matrix_view_t<Vector2>>;

/// @brief Cross product of the 3D vectors @p vec1 and @p vec2.
template<viewable_matrix_r<1> Vector1, viewable_matrix_r<1> Vector2>
constexpr auto cross_product(Vector1&& vec1, Vector2&& vec2) noexcept {
  return CrossProductView(std::forward<Vector1>(vec1),
                          std::forward<Vector2>(vec2));
}

// -----------------------------------------------------------------------------

} // namespace Storm
