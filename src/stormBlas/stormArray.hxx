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
#ifndef _STORM_ARRAY_HXX_
#define _STORM_ARRAY_HXX_

#include <array>
#include <numeric>

#include <stormBlas/stormShape.hxx>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Multiextent array view.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value, class... Extents>
class stormArrayView {
private:
  Value* Values;
  stormShape<Extents...> Shape;

public:

  /// @brief Construct an array view.
  template<class... Indices>
  constexpr explicit stormArrayView(Value* values, 
                                    Indices... dynamicExtents) :
    Values(values), Shape(dynamicExtents...) {
  }

  /// @brief Access array at index.
  /// @{
  template<class... Indices>
  constexpr Value& operator()(Indices... indices) {
    return Values[stormShapeOffset(Shape, indices...)];
  }
  template<class... Indices>
  constexpr Value const& operator()(Indices... indices) const {
    return Values[stormShapeOffset(Shape, indices...)];
  }
  /// @}

}; // class stormArrayView<...>

template<class Value>
using stormVectorView = stormArrayView<Value, stormDynExtent>; 
template<class Value>
using stormMatrixView = stormArrayView<Value, stormDynExtent, stormDynExtent>; 

#endif // ifndef _STORM_ARRAY_HXX_
