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
#include <span>
#include <tuple>
#include <type_traits>

#include <stormBase.hxx>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Multidimensional shape object, \
///   with a support for both static and dynamic extents.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<size_t... Extents>
class stormShape {};

template<size_t Extent, size_t... RestExtents>
class stormShape<Extent, RestExtents...> {
private:

  template<size_t... OtherExtents>
  friend class stormShape;

  [[no_unique_address]] stormShape<RestExtents...> Base_ = {};
  [[no_unique_address]] std::conditional_t<
      Extent == std::dynamic_extent, size_t,
      std::integral_constant<size_t, Extent>>
      Extent_ = {};

public:

  /// @brief Construct a shape object with zero dynamic extents.
  constexpr stormShape() = default;

  // clang-format off

  /// @brief Construct a shape object with the specified dynamic extents.
  /// @{
  template<class... SizeTypes>
    requires(Extent != std::dynamic_extent) 
  constexpr explicit stormShape(
      SizeTypes... dynamicExtents)
      : Base_(dynamicExtents...), Extent_{} {}
  template<class... SizeTypes>
    requires(Extent == std::dynamic_extent) 
  constexpr explicit stormShape(
      size_t frontDynamicExtent, SizeTypes... restDynamicExtents)
      : Base_(restDynamicExtents...), Extent_{frontDynamicExtent} {}
  /// @}

  // clang-format on

  /// @brief Rank of the shape object.
  static constexpr size_t Rank() noexcept {
    return sizeof...(RestExtents) + 1;
  }

  /// @brief Size of the shape object.
  constexpr size_t Size() const noexcept {
    if constexpr (Rank() == 1) {
      return Extent_;
    } else {
      return Extent_ * Base_.Size();
    }
  }

  /// @brief Offset of the shaped data at the specified index.
  template<class... SizeTypes>
  requires((sizeof...(SizeTypes) + 1) == Rank()) constexpr size_t
  operator()(size_t frontIndex, SizeTypes... restIndices) const noexcept {
    stormAssert(frontIndex < Extent_);
    if constexpr (Rank() == 1) {
      return frontIndex;
    } else {
      return frontIndex * Base_.Extent_ + Base_(restIndices...);
    }
  }

}; // class stormShape

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Dense tensor with parametrized storage.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class ValueType, class ShapeType,
         class StorageType = std::unique_ptr<ValueType[]>>
class stormBaseTensor {
private:

  StorageType Storage_;
  [[no_unique_address]] ShapeType Shape_;

public:

  /// @brief Construct a tensor with \
  ///   unassigned storage and zero dynamic extents.
  stormBaseTensor() = default;

  /// @brief Destroy the tensor.
  ~stormBaseTensor() = default;

  /// @brief Assign a shape to the tensor object \
  ///   and allocate the storage (if the storage class is dynamic).
  /// @{
  void Assign(ShapeType const& shape)
  /*requires stormIsOwnedStorage<StorageType> &&
           stormIsDynamicStorage<StorageType>*/
  {
    Shape_ = shape;
    Storage_ = std::make_unique<ValueType[]>(Shape_.Size());
    // Storage_.Allocate(Shape_.Size());
  }
  template<class... Indices>
  void Assign(Indices... dynamicExtents)
  /*requires stormIsOwnedStorage<StorageType> &&
           stormIsDynamicStorage<StorageType>*/
  {
    Assign(ShapeType(dynamicExtents...));
  }
  /// @}

  /// @brief Get tensor data pointer.
  /// @{
  constexpr ValueType* Data() noexcept {
    return Storage_.get();
    // return Storage_.Data();
  }
  constexpr ValueType const* Data() const noexcept {
    return Storage_.get();
    // return Storage_.Data();
  }
  /// @}

  /// @brief Get tensor shape.
  constexpr ShapeType const& Shape() const noexcept {
    return Shape_;
  }

  /// @brief Size of the tensor.
  constexpr size_t Size() const noexcept {
    return Shape_.Size();
  }

  /// @brief Access tensor at index.
  /// @{
  template<class... IndexTypes>
  requires(sizeof...(IndexTypes) == ShapeType::Rank()) constexpr ValueType&
  operator()(IndexTypes... indices) noexcept {
    return Storage_[Shape_(indices...)];
  }
  template<class... IndexTypes>
  requires(sizeof...(IndexTypes) ==
           ShapeType::Rank()) constexpr ValueType const&
  operator()(IndexTypes... indices) const noexcept {
    return Storage_[Shape_(indices...)];
  }
  /// @}

}; // class stormBaseTensor

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Dense tensor.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value, class Shape>
class stormTensor final : public stormBaseTensor<Value, Shape> {
public:
}; // class stormTensor

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Dense vector with dynamic storage.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value>
using stormVector = stormTensor<Value, stormShape<std::dynamic_extent>>;

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Dense matrix with dynamic storage.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value>
using stormMatrix =
    stormTensor<Value, stormShape<std::dynamic_extent, std::dynamic_extent>>;

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Dense rank 2 tensor with dynamic storage.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value>
using stormTensor2R =
    stormTensor<Value, stormShape<std::dynamic_extent, std::dynamic_extent>>;

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Dense rank 3 tensor with dynamic storage.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value>
using stormTensor3R =
    stormTensor<Value, stormShape<std::dynamic_extent, std::dynamic_extent,
                                  std::dynamic_extent>>;

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Dense rank 4 tensor with dynamic storage.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value>
using stormTensor4R =
    stormTensor<Value, stormShape<std::dynamic_extent, std::dynamic_extent,
                                  std::dynamic_extent, std::dynamic_extent>>;
