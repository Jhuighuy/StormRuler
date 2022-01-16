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

#include <array>
#include <memory>

#include <stormBlas/stormShape.hxx>

enum class stormStorageClass {
  Static,
  Dynamic,
}; // enum class stormStorageClass

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Dense tensor.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value, class Shape, 
  stormStorageClass StorageClass = stormStorageClass::Dynamic>
class stormTensor;

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Dense vector with dynamic storage.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value>
using stormVector = 
  stormTensor<Value, stormShape<stormDynExtent<>>>; 

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Dense matrix with dynamic storage.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value>
using stormMatrix = 
  stormTensor<Value, stormShape<stormDynExtent<>, stormDynExtent<>>>;

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Dense tensor with static storage.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value, class Shape>
class stormTensor<Value, Shape, stormStorageClass::Static> {
private:
  std::array<Value, stormShapeMaxSize(std::declval<Shape>())> Values;

public:

}; // class stormTensor<..., stormStorageClass::Static>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Dense tensor with dynamic storage.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value, class Shape>
class stormTensor<Value, Shape, stormStorageClass::Dynamic> {
private:
  std::unique_ptr<Value[]> Values_{};
  Shape Shape_{};

public:
  
  stormTensor() = default;

  /// @brief Construct an array with the specified shape.
  explicit stormTensor(Shape const& shape) {
    Assign(shape);
  }

  void Assign(Shape const& shape) {
    Shape_ = shape;
    Values_ = std::make_unique<Value[]>(stormShapeSize(Shape_));
  }
  template<class... Indices>
  void Assign(Indices... dynamicExtents) {
    Assign(Shape(dynamicExtents...));
  }

  /// @brief Access tensor at index.
  /// @{
  template<class... Indices>
  Value& operator()(Indices... indices) {
    return Values_[stormShapeOffset(Shape_, indices...)];
  }
  template<class... Indices>
  Value const& operator()(Indices... indices) const {
    return Values_[stormShapeOffset(Shape_, indices...)];
  }
  /// @}

}; // class stormTensor<..., stormStorageClass::Dynamic>

#endif // ifndef _STORM_TENSOR_HXX_
