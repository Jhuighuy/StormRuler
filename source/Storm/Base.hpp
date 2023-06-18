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

namespace _detail {

  constexpr bool _in_range(auto t, auto min, auto max) {
    return min <= t && t <= max;
  }

  constexpr bool _one_of(auto x, auto... vals) {
    return ((x == vals) || ...);
  }

  template<class T1, class T2>
  concept _different_from =
      (!std::same_as<std::remove_cvref_t<T1>, std::remove_cvref_t<T2>>);

} // namespace _detail

// -----------------------------------------------------------------------------

} // namespace Storm
