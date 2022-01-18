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

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Base class for the dense tensors with dynamic storage.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value, class ShapeType>
class stormBaseTensor {
private:
  using DataStorage = std::unique_ptr<Value[]>; 
  std::tuple<DataStorage, ShapeType> CompressedStorage_;

protected:
  
  stormBaseTensor() = default;

  ~stormBaseTensor() = default;

  void BaseAssign(ShapeType const& shape) {
    std::get<1>(CompressedStorage_) = shape;
    std::get<0>(CompressedStorage_) = 
      std::make_unique<Value[]>(Shape().Size());
  }

public:

  /// @brief Get tensor shape.
  ShapeType const& Shape() const noexcept {
    return std::get<1>(CompressedStorage_);
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
template<class Value, class Shape>
class stormTensor final : public stormBaseTensor<Value, Shape> {
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
    return this->Shape().Size();
  }

  /// @brief Access tensor at index.
  /// @{
  template<class... Indices>
  Value& operator()(Indices... indices) noexcept {
    return this->Data()[this->Shape()(indices...)];
  }
  template<class... Indices>
  Value const& operator()(Indices... indices) const noexcept {
    return this->Data()[this->Shape()(indices...)];
  }
  /// @}

}; // class stormTensor

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Dense vector with dynamic storage.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value>
using stormVector = 
  stormTensor<Value, 
    stormShape<stormDynamicExtent>>; 

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Dense matrix with dynamic storage.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value>
using stormMatrix = 
  stormTensor<Value, 
    stormShape<
      stormDynamicExtent, 
      stormDynamicExtent>>;

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Dense rank 2 tensor with dynamic storage.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value>
using stormTensor2R = 
  stormTensor<Value, 
    stormShape<
      stormDynamicExtent, 
      stormDynamicExtent>>;

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Dense rank 3 tensor with dynamic storage.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value>
using stormTensor3R = 
  stormTensor<Value, 
    stormShape<
      stormDynamicExtent, 
      stormDynamicExtent, 
      stormDynamicExtent>>;

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Dense rank 4 tensor with dynamic storage.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value>
using stormTensor4R = 
  stormTensor<Value, 
    stormShape<
      stormDynamicExtent, 
      stormDynamicExtent, 
      stormDynamicExtent, 
      stormDynamicExtent>>;

#endif // ifndef _STORM_TENSOR_HXX_
