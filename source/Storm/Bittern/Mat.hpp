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

#include <Storm/Bittern/Matrix.hpp>
#include <Storm/Bittern/MatrixAlgorithms.hpp>

#include <array>
#include <concepts>
#include <utility>

namespace Storm {

/// @brief Statically-sized matrix.
template<class Elem, size_t NumRows, size_t NumCols>
class StaticMatrix final :
    public TargetMatrixInterface<StaticMatrix<Elem, NumRows, NumCols>> {
private:

  std::array<Elem, NumRows * NumCols> elems_{};

public:

  /// @brief Construct a matrix.
  constexpr StaticMatrix(const Elem& init = {}) noexcept {
    fill(init);
  }

  /// @brief Construct a matrix with the elements.
  template<class... RestElems>
    requires (std::convertible_to<RestElems, Elem> && ...) &&
             (sizeof...(RestElems) + 1 == NumRows * NumCols)
  constexpr explicit StaticMatrix(const Elem& first_elem,
                                  const RestElems&... rest_elems) noexcept
      : elems_{first_elem, static_cast<Elem>(rest_elems)...} {}

  /// @brief Construct a matrix with another matrix.
  template<matrix Matrix>
    requires (!std::same_as<StaticMatrix, std::remove_cvref_t<Matrix>>)
  constexpr StaticMatrix(Matrix&& other) noexcept {
    assign(*this, std::forward<Matrix>(other));
  }

  using TargetMatrixInterface<StaticMatrix<Elem, NumRows, NumCols>>::operator=;

  /// @brief Fill the matrix with @p value.
  constexpr void fill(const Elem& value) noexcept {
    for (size_t row_index = 0; row_index < NumRows; ++row_index) {
      for (size_t col_index = 0; col_index < NumCols; ++col_index) {
        (*this)(row_index, col_index) = value;
      }
    }
  }

  /// @brief Matrix shape.
  static constexpr auto shape() noexcept {
    return std::array{NumRows, NumCols};
  }

  /// @brief Get the matrix coefficient at @p row_index and @p col_index.
  /// @{
  constexpr Elem& operator()(size_t row_index, //
                             size_t col_index = 0) noexcept {
    STORM_ASSERT_(in_range(shape(), row_index, col_index),
                  "Indices are out of range!");
    return elems_[row_index * NumCols + col_index];
  }
  constexpr const Elem& operator()(size_t row_index,
                                   size_t col_index = 0) const noexcept {
    STORM_ASSERT_(in_range(shape(), row_index, col_index),
                  "Indices are out of range!");
    return elems_[row_index * NumCols + col_index];
  }
  /// @}

  /// @todo Transition code! Remove me!
  /// @{
  constexpr Elem& operator[](size_t row_index) noexcept {
    return (*this)(row_index);
  }
  constexpr const Elem& operator[](size_t row_index) const noexcept {
    return (*this)(row_index);
  }
  /// @}

}; // class StaticMatrix

/// @brief Statically-sized (small) matrix.
template<class Elem, size_t NumRows, size_t NumCols>
using Mat = StaticMatrix<Elem, NumRows, NumCols>;

/// @brief 2x2 matrix.
template<class Elem>
using Mat2x2 = Mat<Elem, 2, 2>;

/// @brief 3x3 matrix.
template<class Elem>
using Mat3x3 = Mat<Elem, 3, 3>;

/// @brief 4x4 matrix.
template<class Elem>
using Mat4x4 = Mat<Elem, 4, 4>;

/// @brief Statically-sized (small) vector.
template<class Elem, size_t NumRows>
using Vec = Mat<Elem, NumRows, 1>;

/// @brief 2D vector.
template<class Elem>
using Vec2D = Vec<Elem, 2>;

/// @brief 3D vector.
template<class Elem>
using Vec3D = Vec<Elem, 3>;

/// @brief 4D vector.
template<class Elem>
using Vec4D = Vec<Elem, 4>;

} // namespace Storm
