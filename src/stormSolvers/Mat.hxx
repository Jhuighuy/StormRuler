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

#include <cmath>
#include <array>
#include <type_traits>
#include <iostream>

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
  std::array<std::array<Value, SizeY>, SizeX> Data_;

public:

  /// @brief Default constructor.
  constexpr Mat() = default;

  /// @brief Construct the matrix with the initializer list.
  constexpr Mat(std::initializer_list<Value> data) {
    StormAssert(data.size() == SizeX * SizeY);
    std::copy(data.begin(), data.end(), Data());
  }

  /// @brief Get pointer to the beginning of the vector data.
  /// @{
  constexpr Value* Data() noexcept {
    return Data_[0].data();
  }
  constexpr Value const* Data() const noexcept {
    return Data_[0].data();
  }
  /// @}

  /// @brief Get reference to the component at the index. 
  /// @{ 
  constexpr Value& operator()(size_t ix, size_t iy = 0) noexcept {
    StormAssert(ix < SizeX && iy < SizeY);
    return (Data_[ix])[iy];
  }
  constexpr Value const& operator()(size_t ix, size_t iy = 0) const noexcept {
    StormAssert(ix < SizeX && iy < SizeY);
    return (Data_[ix])[iy];
  }
  /// @} 

}; // class Mat<...>

/// @todo Parse types in order to get the output type.
template<class Value, class... Values>
using ResultType = Value; 

/// ----------------------------------------------------------------- ///
/// @brief Make a vector.
/// ----------------------------------------------------------------- ///
template<size_t Size, class Value>
constexpr auto MakeVec(Value const& val) noexcept {

  Vec<Value, Size> vec;

  for (size_t i = 0; i < Size; ++i) {
    vec(i) = val;
  }

  return vec;

} // MakeVec<...>

/// ----------------------------------------------------------------- ///
/// @brief Make an identity matrix. 
/// ----------------------------------------------------------------- ///
template<size_t SizeX, size_t SizeY = SizeX, class Value>
constexpr auto MakeMat(Value const& val) noexcept {

  Mat<Value, SizeX, SizeY> mat;

  for (size_t ix = 0; ix < SizeX; ++ix) {
    for (size_t iy = 0; iy < ix; ++iy) {
      mat(ix, iy) = Value(0.0);
    }
    mat(ix, ix) = val;
    for (size_t iy = ix + 1; iy < SizeY; ++iy) {
      mat(ix, iy) = Value(0.0);
    }
  }

  return mat;

} // MakeMat<...>

/// ----------------------------------------------------------------- ///
/// @brief Transpose a matrix. 
/// ----------------------------------------------------------------- ///
template<class Value, size_t SizeX, size_t SizeY>
constexpr auto Transpose(Mat<Value, SizeX, SizeY> const& mat) noexcept {

  Mat<Value, SizeY, SizeX> out;

  for (size_t ix = 0; ix < SizeX; ++ix) {
    for (size_t iy = 0; iy < SizeY; ++iy) {
      out(iy, ix) = mat(ix, iy);
    }
  }

  return out;

} // Transpose<...>

/// ----------------------------------------------------------------- ///
/// @brief Copy a matrix. 
/// ----------------------------------------------------------------- ///
template<class Value, size_t SizeX, size_t SizeY>
constexpr auto operator+(Mat<Value, SizeX, SizeY> const& mat) noexcept {

  return mat;

} // operator+<...>

/// ----------------------------------------------------------------- ///
/// @brief Negate a matrix. 
/// ----------------------------------------------------------------- ///
template<class Value, size_t SizeX, size_t SizeY>
constexpr auto operator-(Mat<Value, SizeX, SizeY> const& mat) noexcept {

  Mat<ResultType<Value>, SizeX, SizeY> out;

  for (size_t ix = 0; ix < SizeX; ++ix) {
    for (size_t iy = 0; iy < SizeY; ++iy) {
      out(ix, iy) = -mat(ix, iy);
    }
  }

  return out;

} // operator-<...>

/// ----------------------------------------------------------------- ///
/// @brief Add matrices. 
/// ----------------------------------------------------------------- ///
/// @{ 
template<class Value1, class Value2, 
         size_t SizeX, size_t SizeY>
constexpr auto operator+(Mat<Value1, SizeX, SizeY> const& mat1,
                         Mat<Value2, SizeX, SizeY> const& mat2) noexcept {

  Mat<ResultType<Value1, Value2>, SizeX, SizeY> out;

  for (size_t ix = 0; ix < SizeX; ++ix) {
    for (size_t iy = 0; iy < SizeY; ++iy) {
      out(ix, iy) = mat1(ix, iy) + mat2(ix, iy);
    }
  }

  return out;

} // operator+<...>
template<class Value1, class Value2, 
         size_t SizeX, size_t SizeY>
constexpr auto& operator+=(Mat<Value1, SizeX, SizeY>& mat1,
                           Mat<Value2, SizeX, SizeY> const& mat2) noexcept {

  return mat1 = mat1 + mat2;

} // operator+=<...>
/// @} 

/// ----------------------------------------------------------------- ///
/// @brief Subtract matrices. 
/// ----------------------------------------------------------------- ///
/// @{ 
template<class Value1, class Value2, 
         size_t SizeX, size_t SizeY>
constexpr auto operator-(Mat<Value1, SizeX, SizeY> const& mat1,
                         Mat<Value2, SizeX, SizeY> const& mat2) noexcept {

  Mat<ResultType<Value1, Value2>, SizeX, SizeY> out;

  for (size_t ix = 0; ix < SizeX; ++ix) {
    for (size_t iy = 0; iy < SizeY; ++iy) {
      out(ix, iy) = mat1(ix, iy) - mat2(ix, iy);
    }
  }

  return out;

} // operator-<...>
template<class Value1, class Value2, 
         size_t SizeX, size_t SizeY>
constexpr auto& operator-=(Mat<Value1, SizeX, SizeY>& mat1,
                           Mat<Value2, SizeX, SizeY> const& mat2) noexcept {

  return mat1 = mat1 - mat2;

} // operator-=<...>
/// @} 

/// ----------------------------------------------------------------- ///
/// @brief Multiply a matrix by a scalar. 
/// ----------------------------------------------------------------- ///
/// @{ 
template<class Value1, class Value2, 
         size_t SizeX, size_t SizeY>
constexpr auto operator*(Value1 const& val,
                         Mat<Value2, SizeX, SizeY> const& mat) noexcept {

  Mat<ResultType<Value1, Value2>, SizeX, SizeY> out;

  for (size_t ix = 0; ix < SizeX; ++ix) {
    for (size_t iy = 0; iy < SizeY; ++iy) {
      out(ix, iy) = val * mat(ix, iy);
    }
  }

  return out;

} // operator*<...>
template<class Value1, class Value2, 
         size_t SizeX, size_t SizeY>
constexpr auto operator*(Mat<Value1, SizeX, SizeY> const& mat,
                         Value2 const& val) noexcept {

  Mat<ResultType<Value1, Value2>, SizeX, SizeY> out;

  for (size_t ix = 0; ix < SizeX; ++ix) {
    for (size_t iy = 0; iy < SizeY; ++iy) {
      out(ix, iy) = mat(ix, iy) * val;
    }
  }

  return out;

} // operator*<...>
template<class Value1, class Value2, 
         size_t SizeX, size_t SizeY>
constexpr auto& operator*=(Mat<Value1, SizeX, SizeY>& mat,
                           Value2 const& val) noexcept {

  return mat = mat * val;

} // operator*=<...>
/// @} 

/// ----------------------------------------------------------------- ///
/// @brief Multiply matrices (in component-wise sense). 
/// ----------------------------------------------------------------- ///
/// @{ 
template<class Value1, class Value2, 
         size_t SizeX, size_t SizeY>
constexpr auto operator*(Mat<Value1, SizeX, SizeY> const& mat1,
                         Mat<Value2, SizeX, SizeY> const& mat2) noexcept {

  Mat<ResultType<Value1, Value2>, SizeX, SizeY> out;

  for (size_t ix = 0; ix < SizeX; ++ix) {
    for (size_t iy = 0; iy < SizeY; ++iy) {
      out(ix, iy) = mat1(ix, iy) * mat2(ix, iy);
    }
  }

  return out;

} // operator*<...>
template<class Value1, class Value2, 
         size_t SizeX, size_t SizeY>
constexpr auto& operator*=(Mat<Value1, SizeX, SizeY>& mat1,
                           Mat<Value2, SizeX, SizeY> const& mat2) noexcept {

  return mat1 = mat1 * mat2;

} // operator*=<...>
/// @} 

/// ----------------------------------------------------------------- ///
/// @brief Divide a matrix by a scalar. 
/// ----------------------------------------------------------------- ///
/// @{ 
template<class Value1, class Value2, 
         size_t SizeX, size_t SizeY>
constexpr auto operator/(Mat<Value1, SizeX, SizeY> const& mat,
                         Value2 const& val) noexcept {

  Mat<ResultType<Value1, Value2>, SizeX, SizeY> out;

  for (size_t ix = 0; ix < SizeX; ++ix) {
    for (size_t iy = 0; iy < SizeY; ++iy) {
      out(ix, iy) = mat(ix, iy) / val;
    }
  }

  return out;

} // operator/<...>
template<class Value1, class Value2, 
         size_t SizeX, size_t SizeY>
constexpr auto operator/=(Mat<Value1, SizeX, SizeY>& mat,
                          Value2 const& val) noexcept {

  return mat = mat / val;

} // operator/=<...>
/// @} 

/// ----------------------------------------------------------------- ///
/// @brief Divide matrices (in the component-wise sense). 
/// ----------------------------------------------------------------- ///
template<class Value1, class Value2, 
         size_t SizeX, size_t SizeY>
constexpr auto operator/(Mat<Value1, SizeX, SizeY> const& mat1,
                         Mat<Value2, SizeX, SizeY> const& mat2) noexcept {

  Mat<ResultType<Value1, Value2>, SizeX, SizeY> out;

  for (size_t ix = 0; ix < SizeX; ++ix) {
    for (size_t iy = 0; iy < SizeY; ++iy) {
      out(ix, iy) = mat1(ix, iy) / mat2(ix, iy);
    }
  }

  return out;

} // operator/<...>
template<class Value1, class Value2, 
         size_t SizeX, size_t SizeY>
constexpr auto& operator/=(Mat<Value1, SizeX, SizeY>& mat1,
                           Mat<Value2, SizeX, SizeY> const& mat2) noexcept {

  return mat1 = mat1 / mat2;

} // operator/=<...>
/// @} 

/// ----------------------------------------------------------------- ///
/// @brief Dot product of matrices (in the vector sense). 
/// ----------------------------------------------------------------- ///
template<class Value1, class Value2, 
         size_t SizeX, size_t SizeY>
constexpr auto Dot(Mat<Value1, SizeX, SizeY> const& mat1,
                   Mat<Value2, SizeX, SizeY> const& mat2) noexcept {

  auto out = mat1(0, 0) * mat2(0, 0);
  for (size_t iy = 1; iy < SizeY; ++iy) {
    out += mat1(0, iy) * mat2(0, iy);
  }
  for (size_t ix = 1; ix < SizeX; ++ix) {
    for (size_t iy = 0; iy < SizeY; ++iy) {
      out += mat1(ix, iy) * mat2(ix, iy);
    }
  }

  return out;

} // Dot<...>

/// ----------------------------------------------------------------- ///
/// @brief Frobenius norm of a matrix. 
/// ----------------------------------------------------------------- ///
template<class Value, size_t SizeX, size_t SizeY>
constexpr auto Norm(const Mat<Value, SizeX, SizeY>& mat) noexcept {

  return std::sqrt(Dot(mat, mat));

} // Norm<...>

/// ----------------------------------------------------------------- ///
/// @brief Multiply matrices (in normal sense). 
/// ----------------------------------------------------------------- ///
template<class Value1, class Value2, 
         size_t SizeX, size_t SizeY, size_t SizeZ>
constexpr auto MatMul(Mat<Value1, SizeX, SizeY> const& mat1,
                      Mat<Value2, SizeY, SizeZ> const& mat2) noexcept {

  Mat<ResultType<Value1, Value2>, SizeX, SizeZ> out;

  for (size_t ix = 0; ix < SizeX; ++ix) {
    for (size_t iz = 0; iz < SizeZ; ++iz) {
      out(ix, iz) = mat1(ix, 0) * mat2(0, iz);
      for (size_t iy = 1; iy < SizeY; ++iy) {
        out(ix, iz) += mat1(ix, iy) * mat2(iy, iz);
      }
    }
  }

  return out;

} // MatMul<...>

/// ----------------------------------------------------------------- ///
/// @brief Perform a LU decomposition of a square matrix. 
/// ----------------------------------------------------------------- ///
template<class Value, size_t Size>
  requires(std::is_floating_point_v<Value>)
constexpr auto DecomposeLu(Mat<Value, Size, Size> const& mat) noexcept {

  auto lMat = MakeMat<Size>(Value(1.0));
  auto uMat = MakeMat<Size>(Value(0.0));

  for (size_t ix = 0; ix < Size; ++ix) {
    for (size_t iy = 0; iy < ix; ++iy) {
      lMat(ix, iy) = mat(ix, iy);
      for (size_t iz = 0; iz < iy; ++iz) {
        lMat(ix, iy) -= lMat(ix, iz) * uMat(iz, iy);
      }
      lMat(ix, iy) /= uMat(iy, iy);
    }
    for (size_t iy = ix; iy < Size; ++iy) {
      uMat(ix, iy) = mat(ix, iy);
      for (size_t iz = 0; iz < ix; ++iz) {
        uMat(ix, iy) -= lMat(ix, iz) * uMat(iz, iy);
      }
    }
  }

  return std::make_pair(lMat, uMat);

} // DecomposeLu<...>

/// ----------------------------------------------------------------- ///
/// @brief Perform a QR decomposition of a matrix. 
/// ----------------------------------------------------------------- ///
template<class Value, size_t SizeX, size_t SizeY>
  requires(std::is_floating_point_v<Value>)
constexpr auto DecomposeQr(Mat<Value, SizeX, SizeY> const& mat) noexcept {

  Mat<Value, SizeX, SizeY> qMat;
  auto rMat = MakeMat<SizeY, SizeY>(Value(0.0));

  for (size_t ix = 0; ix < SizeX; ++ix) {
    for (size_t iy = 0; iy < SizeY; ++iy) {
      std::abort();
    }
  }

  return std::make_pair(qMat, rMat);

} // DecomposeQr<...>

/// ----------------------------------------------------------------- ///
/// @brief Inverse a square matrix using the LU decomposition. 
/// ----------------------------------------------------------------- ///
template<class Value, size_t Size>
  requires(std::is_floating_point_v<Value>)
constexpr auto InverseLU(Mat<Value, Size, Size> const& mat) noexcept;

/// ----------------------------------------------------------------- ///
/// @brief Print a matrix. 
/// ----------------------------------------------------------------- ///
template<class Value, size_t SizeX, size_t SizeY>
std::ostream& operator<<(std::ostream& out, 
                         Mat<Value, SizeX, SizeY> const& mat) {

  for (size_t ix = 0; ix < SizeX; ++ix) {
    for (size_t iy = 0; iy < SizeY; ++iy) {
      out << mat(ix, iy) << ' ';
    }
    out << std::endl;
  }

  return out;

} // operator<<<...>

} // namespace Storm
