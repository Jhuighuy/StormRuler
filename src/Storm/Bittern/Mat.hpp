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

#include <Storm/Bittern/MatrixView.hpp>

#include <array>
#include <cmath>
#include <iostream>
#include <type_traits>

namespace Storm {

/// @brief Statically-sized matrix.
template<class Value, size_t NumRows, size_t NumCols>
class Mat;

/// @brief 2x2 matrix.
template<class Value>
using Mat2x2 = Mat<Value, 2, 2>;

/// @brief 3x3 matrix.
template<class Value>
using Mat3x3 = Mat<Value, 3, 3>;

/// @brief 4x4 matrix.
template<class Value>
using Mat4x4 = Mat<Value, 4, 4>;

/// @brief Statically-sized vector.
template<class Value, size_t Size>
using Vec = Mat<Value, Size, 1>;

/// @brief 2D vector.
template<class Value>
using Vec2D = Vec<Value, 2>;

/// @brief 3D vector.
template<class Value>
using Vec3D = Vec<Value, 3>;

/// @brief 4D vector.
template<class Value>
using Vec4D = Vec<Value, 4>;

template<class Value, size_t NumRows, size_t NumCols>
class Mat final {
private:

  std::array<std::array<Value, NumCols>, NumRows> coeffs_;

public:

  /// @brief Construct a matrix.
  constexpr Mat(const Value& value = {}) noexcept {
    for (auto& row : coeffs_) {
      row.fill(value);
    }
  }

  /// @brief Construct a matrix with the initializer list.
  constexpr explicit Mat(std::initializer_list<Value> initializer) {
    STORM_ASSERT_(initializer.size() == NumRows * NumCols,
                  "Invalid number of the matrix coefficients!");
    std::copy(initializer.begin(), initializer.end(), data());
  }

  /// @brief Fill the matrix with @p value.
  constexpr void fill(const Value& value) noexcept {
    for (size_t row_index = 0; row_index < NumRows; ++row_index) {
      for (size_t col_index = 0; col_index < NumCols; ++col_index) {
        (*this)(row_index, col_index) = value;
      }
    }
  }

  /// @brief Assign the matrix.
  template<matrix Matrix>
  constexpr auto& operator=(Matrix&& other) noexcept {
    for (size_t row_index = 0; row_index < NumRows; ++row_index) {
      for (size_t col_index = 0; col_index < NumCols; ++col_index) {
        (*this)(row_index, col_index) = other(row_index, col_index);
      }
    }
    return *this;
  }

  /// @brief Number of the matrix rows.
  [[nodiscard]] static constexpr auto num_rows() noexcept {
    return NumRows;
  }
  /// @brief Number of the matrix columns.
  [[nodiscard]] static constexpr auto num_cols() noexcept {
    return NumCols;
  }
  /// @brief Matrix shape.
  [[nodiscard]] static constexpr auto shape() noexcept {
    return MatrixShape{NumRows, NumCols};
  }

  /// @brief Matrix size.
  [[nodiscard]] static constexpr auto size() noexcept {
    return NumRows * NumCols;
  }
  /// @brief Matrix data.
  /// @{
  [[nodiscard]] constexpr Value* data() noexcept {
    return &coeffs_[0][0];
  }
  [[nodiscard]] constexpr const Value* data() const noexcept {
    return &coeffs_[0][0];
  }
  /// @}

  /// @brief Get reference to the component at the index.
  /// @{
  [[nodiscard]] constexpr Value& //
  operator()(size_t row_index, size_t col_index = 0) noexcept {
    STORM_ASSERT_(row_index < NumRows && col_index < NumCols,
                  "Indices are out of range!");
    return (coeffs_[row_index])[col_index];
  }
  [[nodiscard]] constexpr const Value&
  operator()(size_t row_index, size_t col_index = 0) const noexcept {
    STORM_ASSERT_(row_index < NumRows && col_index < NumCols,
                  "Indices are out of range!");
    return (coeffs_[row_index])[col_index];
  }
  /// @}

  /// @todo Transition code! Remove me!
  /// @{
  [[nodiscard]] constexpr Value& //
  operator[](size_t row_index) noexcept {
    return (*this)(row_index);
  }
  [[nodiscard]] constexpr const Value& //
  operator[](size_t row_index) const noexcept {
    return (*this)(row_index);
  }
  /// @}

}; // class Mat

} // namespace Storm
