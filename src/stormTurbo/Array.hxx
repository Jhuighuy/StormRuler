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
#include <cassert>

#include <stormBase.hxx>

namespace Storm::Turbo {

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief A simple non-copyable 2D array.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value>
class Array2D final : public NonCopyable {
private:
  size_t SizeX_, SizeY_;
  std::unique_ptr<Value[]> Data_;

public:

  /// @brief Initialize an empty 2D array.
  Array2D() : SizeX_(0), SizeY_(0) {}

  /// @brief Initialize a 2D array with the specified size.
  Array2D(size_t sizeX, size_t sizeY) : Array2D() {
    Assign(sizeX, sizeY);
  }

  /// @brief Assign size to the array.
  void Assign(size_t sizeX, size_t sizeY) {
    SizeX_ = sizeX, SizeY_ = sizeY;
    Data_.reset(new Value[SizeX_ * SizeY_]{});
  }

  /// @brief Size in X direction.
  size_t SizeX() const noexcept {
    return SizeX_;
  }

  /// @brief Size in Y direction.
  size_t SizeY() const noexcept {
    return SizeY_;
  }

  /// @brief Total size.
  size_t Size() const noexcept {
    return SizeX_ * SizeY_;
  }

  /// @brief Check if array is empty.
  bool Empty() const noexcept {
    return SizeX_ | SizeY_ != 0;
  }

  /// @brief Get pointer to the beginning of the array data.
  /// @{
  Value* Data() noexcept {
    return Data_.get();
  }
  Value const* Data() const noexcept {
    return Data_.get();
  }
  /// @}

  /// @brief Get reference to the value at the 2D and block index.
  /// @{
  Value& operator()(size_t ix, size_t iy) noexcept {
    StormAssert(ix < SizeX_ && iy < SizeY_);
    return Data_[ix * SizeY_ + iy];
  }
  Value const& operator()(size_t ix, size_t iy) const noexcept {
    StormAssert(ix < SizeX_ && iy < SizeY_);
    return Data_[ix * SizeY_ + iy];
  }
  /// @}

}; // class Array2D<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief A simple non-copyable 3D array.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value>
class Array3D final : public NonCopyable {
private:
  size_t SizeX_, SizeY_, SizeZ_;
  std::unique_ptr<Value[]> Data_;

public:

  /// @brief Initialize an empty 3D array.
  Array3D() : SizeX_(0), SizeY_(0), SizeZ_(0) {}

  /// @brief Initialize a 3D array with the specified size.
  Array3D(size_t sizeX, size_t sizeY, size_t sizeZ) : Array3D() {
    Assign(sizeX, sizeY, sizeZ);
  }

  /// @brief Assign size to the array.
  void Assign(size_t sizeX, size_t sizeY, size_t sizeZ) {
    SizeX_ = sizeX, SizeY_ = sizeY, SizeZ_ = sizeZ;
    Data_.reset(new Value[SizeX_ * SizeY_ * SizeZ_]{});
  }

  /// @brief Size in X direction.
  size_t SizeX() const noexcept {
    return SizeX_;
  }

  /// @brief Size in Y direction.
  size_t SizeY() const noexcept {
    return SizeY_;
  }

  /// @brief Size in Z direction.
  size_t SizeZ() const noexcept {
    return SizeZ_;
  }

  /// @brief Total size.
  size_t Size() const noexcept {
    return SizeX_ * SizeY_ * SizeZ_;
  }

  /// @brief Check if array is empty.
  bool Empty() const noexcept {
    return SizeX_ | SizeY_ | SizeZ_ != 0;
  }

  /// @brief Get pointer to the beginning of the array data.
  /// @{
  Value* Data() noexcept {
    return Data_.get();
  }
  Value const* Data() const noexcept {
    return Data_.get();
  }
  /// @}

  /// @brief Get reference to the value at the 2D and block index.
  /// @{
  Value& operator()(size_t ix, size_t iy, size_t iz) noexcept {
    StormAssert(ix < SizeX_ && iy < SizeY_ && iz < SizeZ_);
    return Data_[(ix * SizeY_ + iy) * SizeZ_ + iz];
  }
  Value const& operator()(size_t ix, size_t iy, size_t iz) const noexcept {
    StormAssert(ix < SizeX_ && iy < SizeY_ && iz < SizeZ_);
    return Data_[(ix * SizeY_ + iy) * SizeZ_ + iz];
  }
  /// @}

}; // class Array3D<...>

} // namespace Storm::Turbo
