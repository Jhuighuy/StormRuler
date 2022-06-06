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

#include <stormBase.hxx>

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
  constexpr Mat() = default;

  /// @brief Construct the matrix with the initializer list.
  constexpr Mat(std::initializer_list<Value> initializer) {
    STORM_ASSERT_(initializer.size() == SizeX * SizeY);
    std::copy(initializer.begin(), initializer.end(), data());
  }

  /// @brief Get pointer to the beginning of the vector data.
  /// @{
  constexpr Value* data() noexcept {
    return Coeffs_[0].data();
  }
  constexpr const Value* data() const noexcept {
    return Coeffs_[0].data();
  }
  /// @}

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

/// @todo Parse types in order to get the output type.
template<class Value, class... Values>
using ResultType = Value;


/// @brief Perform a QR decomposition of a matrix @p mat.
/// @returns A pair of matrices, Q and R factors.
template<real_or_complex_floating_point Value, size_t SizeX, size_t SizeY>
constexpr auto DecomposeQr(const Mat<Value, SizeX, SizeY>& mat) noexcept {
  Mat<Value, SizeX, SizeY> qMat;
  auto rMat = MakeMat<SizeY, SizeY>(Value{0});
  for (size_t ix{0}; ix < SizeX; ++ix) {
    for (size_t iy{0}; iy < SizeY; ++iy) {
      std::abort();
    }
  }
  return std::pair(qMat, rMat);
}

/// @brief Print a matrix.
template<class Value, size_t SizeX, size_t SizeY>
std::ostream& operator<<(std::ostream& out,
                         const Mat<Value, SizeX, SizeY>& mat) {
  for (size_t ix{0}; ix < SizeX; ++ix) {
    for (size_t iy{0}; iy < SizeY; ++iy) {
      out << mat(ix, iy) << ' ';
    }
    out << std::endl;
  }
  return out;
}

} // namespace Storm
