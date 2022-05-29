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
/// @brief Dynamic extent specifier.
/// ----------------------------------------------------------------- ///
static constexpr stormSize_t stormDynamicExtent = STORM_SIZE_MAX;

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Multidimensional shape object, \
///   with a support for both static and dynamic extents.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<stormSize_t... Extents>
class stormShape {};

template<stormSize_t Extent, stormSize_t... RestExtents>
class stormShape<Extent, RestExtents...> {
private:

  template<stormSize_t... OtherExtents>
  friend class stormShape;

  [[no_unique_address]] stormShape<RestExtents...> Base_ = {};
  [[no_unique_address]] std::conditional_t<
    Extent == stormDynamicExtent, 
    stormSize_t, 
    std::integral_constant<stormSize_t, Extent>> Extent_ = {};

public:

  /// @brief Construct a shape object with \
  ///   zero dynamic extents.
  constexpr stormShape() = default;

  /// @brief Construct a shape object with \
  ///   the specified dynamic extents.
  /// @{
  template<class... SizeTypes>
    requires (Extent != stormDynamicExtent)
  constexpr explicit stormShape(SizeTypes... dynamicExtents) :
      Base_(dynamicExtents...), Extent_{} {
  }
  template<class... SizeTypes>
    requires (Extent == stormDynamicExtent)
  constexpr explicit stormShape(stormSize_t frontDynamicExtent,
                                SizeTypes... restDynamicExtents) :
      Base_(restDynamicExtents...), Extent_{frontDynamicExtent} {
  }
  /// @}

  /// @brief Rank of the shape object.
  static constexpr stormSize_t Rank() noexcept {
    return sizeof...(RestExtents) + 1;
  }

  /// @brief Size of the shape object.
  constexpr stormSize_t Size() const noexcept {
    if constexpr (Rank() == 1) {
      return Extent_;
    } else {
      return Extent_*Base_.Size();
    }
  }

  /// @brief Offset of the shaped data at the specified index.
  template<class... SizeTypes>
    requires ((sizeof...(SizeTypes) + 1) == Rank())
  constexpr stormSize_t operator()(stormSize_t frontIndex,
                                   SizeTypes... restIndices) const noexcept {
    stormAssert(frontIndex < Extent_);
    if constexpr (Rank() == 1) {
      return frontIndex;
    } else {
      return frontIndex*Base_.Extent_ + Base_(restIndices...);
    }
  }

}; // class stormShape

#endif // ifndef _STORM_SHAPE_HXX_
