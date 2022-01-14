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
#ifndef _STORM_SHAPE_HXX_
#define _STORM_SHAPE_HXX_

#include <StormRuler_API.h>

/// ----------------------------------------------------------------- ///
/// @brief Static array extent tag.
/// ----------------------------------------------------------------- ///
template<stormSize_t>
struct stormExtent;

/// ----------------------------------------------------------------- ///
/// @brief Dynamic array extent tag.
/// ----------------------------------------------------------------- ///
struct stormDynExtent;

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Multiextent shape object, \
///   with a support for both static and dynamic storage classes.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class... Extents>
class stormShape;

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Zero rank shape object.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<>
class stormShape<> {

}; // stormShape<>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Multiextent shape object, with a static leading extent.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<stormSize_t StaticExtent, class... RestExtents>
class stormShape<stormExtent<StaticExtent>, RestExtents...> : 
    public stormShape<RestExtents...> {
public:
  static_assert(StaticExtent > 0, 
    "Static extent of a shape object must be greater than zero.");

  /// @brief Construct a shape object.
  template<class... Indices>
  constexpr explicit stormShape(Indices... restDynamicExtents) :
    stormShape<RestExtents...>(restDynamicExtents...) {
  }

  /// @brief Leading extent of the shape object.
  constexpr stormSize_t FrontExtent() const {
    return StaticExtent;
  }

}; // stormShape<stormExtent<...>, ...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Multiextent shape object, with a dynamic leading extent.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class... RestExtents>
class stormShape<stormDynExtent, RestExtents...> : 
    public stormShape<RestExtents...> {
private:
  stormSize_t const DynamicExtent;

public:
  /// @brief Construct a shape object.
  template<class... Indices>
  constexpr explicit stormShape(stormSize_t dynamicExtent,
                                Indices... restDynamicExtents) :
    stormShape<RestExtents...>(restDynamicExtents...), 
    DynamicExtent(dynamicExtent) {
  }

  /// @brief Front extent of the shape object.
  constexpr stormSize_t FrontExtent() const {
    return DynamicExtent;
  }

}; // stormShape<stormDynExtent, ...>

/// ----------------------------------------------------------------- ///
/// @brief Rank of the shape object.
/// ----------------------------------------------------------------- ///
template<class... Extents>
constexpr stormSize_t stormShapeRank(stormShape<Extents...> const&) {
  return sizeof...(Extents);
}

/// ----------------------------------------------------------------- ///
/// @brief Size of the shape object.
/// ----------------------------------------------------------------- ///
/// @{
constexpr stormSize_t stormShapeSize(
    stormShape<> const&) {
  return 1;
}
template<class FrontExtent, class... RestExtents>
constexpr stormSize_t stormShapeSize(
    stormShape<FrontExtent, RestExtents...> const& shape) {
  return shape.FrontExtent() *
    stormShapeSize(static_cast<stormShape<RestExtents...> const&>(shape));
}
/// @}

/// ----------------------------------------------------------------- ///
/// @brief Offset of the shaped data at the specified index.
/// ----------------------------------------------------------------- ///
/// @{
constexpr stormSize_t stormShapeOffset(
    stormShape<> const&) {
  return 0;
}
template<class FrontExtent, class... RestExtents, class... RestIndices>
constexpr stormSize_t stormShapeOffset(
    stormShape<FrontExtent, RestExtents...> const& shape,
    stormSize_t frontIndex,
    RestIndices... restIndices) {
  _STORM_ASSERT_(frontIndex < shape.FrontExtent());
  return frontIndex * 
      stormShapeSize(static_cast<stormShape<RestExtents...> const&>(shape)) +
    stormShapeOffset(
      static_cast<stormShape<RestExtents...> const&>(shape), restIndices...);
}
/// @}

#endif // ifndef _STORM_SHAPE_HXX_
