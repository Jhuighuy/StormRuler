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
#include <cmath>
#include <iostream>
#include <type_traits>

#include <Storm/Base.hpp>

#include <Storm/Blass/MatrixView.hpp>

namespace Storm {

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Statically-sized matrix.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value, size_t SizeX, size_t SizeY>
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

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Statically-sized vector.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
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

template<class Value, size_t SizeX, size_t SizeY>
class Mat final {
private:

  std::array<std::array<Value, SizeY>, SizeX> Coeffs_;

public:

  /// @brief Default constructor.
  constexpr Mat() {
    for (auto& a : Coeffs_)
      a.fill({});
  }

  /// @brief Construct the matrix with the initializer list.
  constexpr explicit Mat(std::initializer_list<Value> initializer) {
    STORM_ASSERT_(initializer.size() == SizeX * SizeY);
    std::copy(initializer.begin(), initializer.end(), Coeffs_[0].data());
  }

  constexpr auto& operator=(const Value& v) noexcept {
    for (size_t ix = 0; ix < SizeX; ++ix) {
      for (size_t iy = 0; iy < SizeY; ++iy) {
        (*this)(ix, iy) = v;
      }
    }
    return *this;
  }

  constexpr auto& operator=(const is_matrix_view auto& other) noexcept {
    for (size_t ix = 0; ix < SizeX; ++ix) {
      for (size_t iy = 0; iy < SizeY; ++iy) {
        (*this)(ix, iy) = other(ix, iy);
      }
    }
    return *this;
  }

  auto num_rows() const noexcept {
    return std::integral_constant<size_t, SizeX>{};
  }
  auto num_cols() const noexcept {
    return std::integral_constant<size_t, SizeY>{};
  }

  /// @brief Get reference to the component at the index.
  /// @{
  constexpr Value& operator()(size_t ix, size_t iy = 0) noexcept {
    STORM_ASSERT_(ix < SizeX && iy < SizeY);
    return (Coeffs_[ix])[iy];
  }
  constexpr const Value& operator()(size_t ix, size_t iy = 0) const noexcept {
    STORM_ASSERT_(ix < SizeX && iy < SizeY);
    return (Coeffs_[ix])[iy];
  }
  /// @}

}; // class Mat

template<class Value, size_t SizeX, size_t SizeY>
struct is_matrix_t<Mat<Value, SizeX, SizeY>> : std::true_type {};

template<class Value, size_t SizeX, size_t SizeY>
constexpr auto& eval(auto func, Mat<Value, SizeX, SizeY>& mat_lhs,
                     const is_matrix_view auto&... mats_rhs) noexcept {
  for (size_t row_index{0}; row_index < mat_lhs.num_rows(); ++row_index) {
    for (size_t col_index{0}; col_index < mat_lhs.num_cols(); ++col_index) {
      func(mat_lhs(row_index, col_index), mats_rhs(row_index, col_index)...);
    }
  }
  return mat_lhs;
}

template<class Value, size_t SizeX, size_t SizeY>
constexpr real_t dot_product(Mat<Value, SizeX, SizeY> m1,
                             Mat<Value, SizeX, SizeY> m2) {
  Value d{};
  for (size_t ix = 0; ix < SizeX; ++ix) {
    for (size_t iy = 0; iy < SizeY; ++iy) {
      d += m1(ix, iy) * m2(ix, iy);
    }
  }
  return d;
}

} // namespace Storm
