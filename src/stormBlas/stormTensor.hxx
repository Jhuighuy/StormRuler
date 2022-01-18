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
#ifndef _STORM_TENSOR_HXX_
#define _STORM_TENSOR_HXX_

#include <type_traits>
#include <tuple>
#include <array>
#include <memory>

#include <stormBlas/stormShape.hxx>

/// ----------------------------------------------------------------- ///
/// @brief Storage class.
/// ----------------------------------------------------------------- ///
enum class stormStorageClass {

  /// @brief The data is stored inside of the object.
  Static, 

  /// @brief The data is stored on the heap, \
  ///   only the pointer to it is stored.
  Dynamic,

}; // enum class stormStorageClass

/// ----------------------------------------------------------------- ///
/// @brief Base class for the dense tensors.
/// ----------------------------------------------------------------- ///
template<class Value, class ShapeType, stormStorageClass StorageClass>
class stormBaseTensor;

/// ----------------------------------------------------------------- ///
/// @brief Base class for the dense tensors with static storage.
/// ----------------------------------------------------------------- ///
template<class Value, class ShapeType>
class stormBaseTensor<Value, ShapeType, stormStorageClass::Static> {
private:
  using DataStorage = 
    std::array<Value, stormShapeMaxSize(std::declval<ShapeType>())>;
  std::conditional_t<std::is_empty_v<ShapeType>,
    std::tuple<DataStorage>, 
    std::tuple<DataStorage, ShapeType>> CompressedStorage_;

protected:
  
  stormBaseTensor() = default;

  ~stormBaseTensor() = default;

  void BaseAssign(ShapeType const& shape) {
    if constexpr (!std::is_empty_v<ShapeType>) {
      std::get<1>(CompressedStorage_) = shape;
    }
  }

public:

  /// @brief Get tensor data pointer.
  /// @{
  Value* Data() noexcept {
    return std::get<0>(CompressedStorage_).data();
  }
  Value const* Data() const noexcept {
    return std::get<0>(CompressedStorage_).data();
  }
  /// @}

  /// @brief Get tensor shape.
  ShapeType const& Shape() const noexcept {
    if constexpr (std::is_empty_v<ShapeType>) {
      return static_cast<ShapeType const&>(*this);
    } else {
      return std::get<1>(CompressedStorage_);
    }
  }

}; // class stormBaseTensor<..., stormStorageClass::Static>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Base class for the dense tensors with dynamic storage.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value, class ShapeType>
class stormBaseTensor<Value, ShapeType, stormStorageClass::Dynamic> {
private:
  using DataStorage = std::unique_ptr<Value[]>; 
  std::conditional_t<std::is_empty_v<ShapeType>,
    std::tuple<DataStorage>, 
    std::tuple<DataStorage, ShapeType>> CompressedStorage_;

protected:
  
  stormBaseTensor() = default;

  ~stormBaseTensor() = default;

  void BaseAssign(ShapeType const& shape) {
    if constexpr (!std::is_empty_v<ShapeType>) {
      std::get<1>(CompressedStorage_) = shape;
    }
    std::get<0>(CompressedStorage_) = 
      std::make_unique<Value[]>(stormShapeSize(Shape()));
  }

public:

  /// @brief Get tensor shape.
  ShapeType const& Shape() const noexcept {
    if constexpr (std::is_empty_v<ShapeType>) {
      return static_cast<ShapeType const&>(*this);
    } else {
      return std::get<1>(CompressedStorage_);
    }
  }

  /// @brief Get tensor data pointer.
  /// @{
  Value* Data() noexcept {
    return std::get<0>(CompressedStorage_).get();
  }
  Value const* Data() const noexcept {
    return std::get<0>(CompressedStorage_).get();
  }
  /// @}

}; // class stormBaseTensor<..., stormStorageClass::Dynamic>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Dense tensor.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value, class Shape, 
  stormStorageClass StorageClass = stormStorageClass::Dynamic>
class stormTensor final : public stormBaseTensor<Value, Shape, StorageClass> {
public:

  /// @brief Assign a shape to the tensor object \
  ///   and allocate the storage (if the storage class is dynamic).
  /// @{
  void Assign(Shape const& shape) {
    this->BaseAssign(shape);
  }
  template<class... Indices>
  void Assign(Indices... dynamicExtents) {
    this->BaseAssign(Shape(dynamicExtents...));
  }
  /// @}

  /// @brief Size of the tensor.
  stormSize_t Size() const noexcept {
    return stormShapeSize(this->Shape());
  }

  /// @brief Access tensor at index.
  /// @{
  template<class... Indices>
  Value& operator()(Indices... indices) noexcept {
    return this->Data()[stormShapeOffset(this->Shape(), indices...)];
  }
  template<class... Indices>
  Value const& operator()(Indices... indices) const noexcept {
    return this->Data()[stormShapeOffset(this->Shape(), indices...)];
  }
  /// @}

}; // class stormTensor

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Dense vector with dynamic storage.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value>
using stormVector = 
  stormTensor<Value, 
    stormShape<stormDynExtent<>>>; 

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Dense matrix with dynamic storage.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value>
using stormMatrix = 
  stormTensor<Value, 
    stormShape<
      stormDynExtent<>, 
      stormDynExtent<>>>;

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Dense rank 2 tensor with dynamic storage.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value>
using stormTensor2R = 
  stormTensor<Value, 
    stormShape<
      stormDynExtent<>, 
      stormDynExtent<>>>;

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Dense rank 3 tensor with dynamic storage.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value>
using stormTensor3R = 
  stormTensor<Value, 
    stormShape<
      stormDynExtent<>, 
      stormDynExtent<>, 
      stormDynExtent<>>>;

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Dense rank 4 tensor with dynamic storage.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value>
using stormTensor4R = 
  stormTensor<Value, 
    stormShape<
      stormDynExtent<>, 
      stormDynExtent<>, 
      stormDynExtent<>, 
      stormDynExtent<>>>;

#endif // ifndef _STORM_TENSOR_HXX_