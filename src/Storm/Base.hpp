/// Copyright (C) 2022 Oleg Butakov
///
/// Permission is hereby granted, free of charge, to any person obtaining a copy
/// of this software and associated documentation files (the "Software"), to
/// deal in the Software without restriction, including without limitation the
/// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
/// sell copies of the Software, and to permit persons to whom the Software is
/// furnished to do so, subject to the following conditions:
///
/// The above copyright notice and this permission notice shall be included in
/// all copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
/// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
/// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
/// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
/// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
/// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
/// IN THE SOFTWARE.

#pragma once

#ifndef STORM_HEADER_ONLY
#error #define STORM_HEADER_ONLY before including Storm!
#endif

#ifdef STORM_HEADER_ONLY
#define FMT_HEADER_ONLY 1
#endif

#include <StormRuler_API.h>

#include <cassert>
#include <cmath>
#include <cstddef>

#include <complex>
#include <concepts>
#include <source_location>
#include <tuple>
#include <type_traits>

/// @todo Reimplement me with std::source_location.
#if (!defined(__PRETTY_FUNCTION__) && !defined(__GNUC__))
#define __PRETTY_FUNCTION__ __FUNCSIG__
#endif

#define STORM_FATAL_ERROR_(error_message)                                      \
  do {                                                                         \
    std::fprintf(stderr, "\nAssertion failed:\n%s:%d %s: \"%s\".\n", __FILE__, \
                 __LINE__, __PRETTY_FUNCTION__, error_message);                \
    std::fflush(stderr);                                                       \
    std::abort();                                                              \
  } while (false)

#define STORM_NOT_IMPLEMENTED_() STORM_FATAL_ERROR_("not implemented")

#ifdef _MSC_VER
#define STORM_ASSUME_(x) __assume(x)
#else
#define STORM_ASSUME_(x)                   \
  do {                                     \
    if (!(x)) { __builtin_unreachable(); } \
  } while (false)
#endif

#define STORM_ENSURE_(x)                  \
  do {                                    \
    if (!(x)) { STORM_FATAL_ERROR_(#x); } \
  } while (false)

#ifdef NDEBUG
#define STORM_ASSERT_(x) STORM_ASSUME_(x)
#else
#define STORM_ASSERT_(x) STORM_ENSURE_(x)
#endif

#ifdef _MSC_VER
#define STORM_FORCE_INLINE_ inline __forceinline
#else
#define STORM_FORCE_INLINE_ inline __attribute__((always_inline))
#endif

#ifdef _MSC_VER
#define STORM_NO_UNIQUE_ADDRESS_ [[msvc::no_unique_address]]
#else
#define STORM_NO_UNIQUE_ADDRESS_ [[no_unique_address]]
#endif

namespace Storm {

/// @brief Contains the implementation details.
namespace Detail_ {}

/// @brief Size type.
using size_t = std::size_t;

/// @brief Pointer difference type.
using ptrdiff_t = std::ptrdiff_t;

/// @brief Real floating-point type.
using real_t = double;

namespace detail_ {

  template<class>
  constexpr inline bool always_false_ = false;

  consteval bool in_range_(auto t, auto min, auto max) {
    return min <= t && t <= max;
  }

  template<class T1, class T2>
  concept different_from_ =
      !std::same_as<std::remove_cvref_t<T1>, std::remove_cvref_t<T2>>;

} // namespace detail_

/// @brief @c size_t compile-time constant.
template<size_t N>
using size_t_constant = std::integral_constant<size_t, N>;

namespace detail_ {
  template<class>
  inline constexpr bool is_complex_floating_point_v_ = false;
  template<class T>
  inline constexpr bool is_complex_floating_point_v_<std::complex<T>> =
      std::is_floating_point_v<T>;
} // namespace detail_

/// @brief Real or complex floaing point numbers,
///   e.g. @c real_t or @c std::complex<real_t>.
template<class T>
concept real_or_complex_floating_point =
    std::floating_point<T> || detail_::is_complex_floating_point_v_<T>;

} // namespace Storm
