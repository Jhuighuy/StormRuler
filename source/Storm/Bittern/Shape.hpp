// Copyright (C) 2020-2023 Oleg Butakov
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <Storm/Base.hpp>

#include <Storm/Crow/TupleUtils.hpp>

#include <array>
#include <tuple>
#include <type_traits>
#include <utility>

namespace Storm {

// -----------------------------------------------------------------------------

/// @brief Extent: convertible to size_t object that is not a floating-point.
template<class Extent>
concept extent = std::is_object_v<Extent> && detail::index_like<Extent>;

/// @brief Fixed extent.
template<class Extent>
concept fixed_extent = extent<Extent> && std::is_empty_v<Extent>;

/// @brief Fixed extent type.
template<size_t Arg>
using fixed_extent_t = fixed_size_t<Arg>;

/// @brief Dynamic extent type.
using dynamic_extent_t = size_t;

// -----------------------------------------------------------------------------

namespace detail {
  template<tuple_like Shape>
  inline constexpr bool has_extents_v =
      base_apply<std::tuple_size_v<Shape>>([](auto... axes) {
        return (extent<std::tuple_element_t<axes, Shape>> && ...);
      });
  template<tuple_like Shape>
  inline constexpr bool has_fixed_extents_v =
      base_apply<std::tuple_size_v<Shape>>([](auto... axes) {
        return (... && fixed_extent<std::tuple_element_t<axes, Shape>>);
      });
} // namespace detail

/// @brief Shape: a tuple (or a tuple-like object) of extents.
template<class Shape>
concept shape = tuple_like<Shape> && detail::has_extents_v<Shape>;

/// @brief Shape rank.
template<shape Shape>
inline constexpr size_t shape_rank_v = std::tuple_size_v<Shape>;

/// @brief Shape extent type.
template<shape Shape, size_t Axis>
  requires (Axis < shape_rank_v<Shape>)
using shape_extent_t = std::tuple_element_t<Axis, Shape>;

/// @brief Fixes shape.
template<class Shape>
concept fixed_shape = shape<Shape> && detail::has_fixed_extents_v<Shape>;

/// @brief Fixed shape type.
template<size_t... Args>
using fixed_shape_t = std::tuple<fixed_extent_t<Args>...>;

/// @brief Make a dynamic shape type.
template<size_t Rank>
using dynamic_shape_t = std::array<size_t, Rank>;

// -----------------------------------------------------------------------------

/// @brief Concatenate shapes.
template<shape... Shapes>
constexpr auto cat_shapes(Shapes... shapes) {
  return std::tuple_cat(std::move(shapes)...);
}

/// @brief Concatenate shape types.
template<shape... Shapes>
using cat_shapes_t = decltype(cat_shapes(std::declval<Shapes>()...));

// -----------------------------------------------------------------------------

/// @brief Common shape.
template<shape... Shapes>
constexpr auto common_shape(Shapes... shapes) {
  const auto common_extents = [](auto... extents_tuple) {
    const auto common_extent = [](auto first_extent, auto... rest_extents) {
      STORM_ASSERT(((first_extent == rest_extents) && ...),
                   "Shape extents does not match!");
      using Extent = std::common_type_t<decltype(first_extent), //
                                        decltype(rest_extents)...>;
      return static_cast<Extent>(first_extent);
    };
    return std::tuple{std::apply(common_extent, extents_tuple)...};
  };
  return zip_apply(common_extents, {}, shapes...);
}

/// @brief Common shape type.
template<shape... Shapes>
using common_shape_t = decltype(common_shape(std::declval<Shapes>()...));

// -----------------------------------------------------------------------------

/// @brief Check if all indices @p indices are in range of the
/// shape @p shape extents.
template<class Shape, class... Indices>
[[nodiscard]] constexpr bool in_range(Shape shape, Indices... indices) {
#ifdef __MINGW64__
  return true; // Looks like a bug in the linker.
#else
  const auto check_axes = [](auto... extent_index_pairs) {
    const auto check_axis = [](auto extent, auto index) {
      return static_cast<size_t>(index) < static_cast<size_t>(extent);
    };
    return (std::apply(check_axis, extent_index_pairs) && ...);
  };
  return zip_apply(check_axes, {}, shape, std::tuple{indices...});
#endif
}

// -----------------------------------------------------------------------------

} // namespace Storm
