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

#include <Storm/Crow/Base/Assert.hpp>
#include <Storm/Crow/Base/Config.hpp>
#include <Storm/Crow/Base/Exception.hpp>
#include <Storm/Crow/Base/Log.hpp>
#include <Storm/Crow/Base/Types.hpp>

#include <concepts>
#include <type_traits>

namespace Storm {

// -----------------------------------------------------------------------------

namespace detail {

  // ---------------------------------------------------------------------------

  // To be used inside `static_asserts`.
  template<class>
  inline constexpr bool always_false = false;

  // ---------------------------------------------------------------------------

  // Index-like object.
  template<class Index>
  concept index_like = (!std::is_reference_v<Index>) &&(
      std::convertible_to<Index, size_t> &&
      !std::floating_point<Index>) &&std::copyable<Index>;

  // Multidimensional `for` statement.
  template<class Func, index_like FirstExtent, index_like... RestExtents>
    requires std::regular_invocable<Func, FirstExtent, RestExtents...>
  constexpr void mdfor(Func&& func, FirstExtent first_extent,
                       RestExtents... rest_extents) {
    for (size_t i = 0; i < static_cast<size_t>(first_extent); ++i) {
      if constexpr (sizeof...(RestExtents) == 0) func(i);
      else mdfor([&](auto... rest) { func(i, rest...); }, rest_extents...);
    }
  }

  // ---------------------------------------------------------------------------

} // namespace detail

// -----------------------------------------------------------------------------

namespace _detail {

  constexpr bool _in_range(auto t, auto min, auto max) {
    return min <= t && t <= max;
  }

} // namespace _detail

// -----------------------------------------------------------------------------

} // namespace Storm
