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

#include <complex>
#include <concepts>
#include <cstdint>
#include <cstdlib>
#include <source_location>
#include <type_traits>

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <StormRuler_API.h>

// Detect the C++ version.
#if __cplusplus < 202002L
#error Storm requires C++20 support!
#endif
#if __cplusplus > 202002L
#define STORM_CPP23_ 1
#else
#define STORM_CPP23_ 0
#endif

// Detect the C++ compiler.
/// @todo Implement me!
#ifdef __GNUC__
#define STORM_COMPILER_GCC_ 1
#else
#define STORM_COMPILER_GCC_ 0
#endif
#ifdef __clang__
#define STORM_COMPILER_CLANG_ 1
#else
#define STORM_COMPILER_CLANG_ 0
#endif
#ifdef _MSC_VER
#define STORM_COMPILER_MSVC_ 1
#else
#define STORM_COMPILER_MSVC_ 0
#endif
#if (STORM_COMPILER_GCC_ + STORM_COMPILER_CLANG_ + STORM_COMPILER_MSVC_) != 1
#warning Storm has detected multiple compilers, something is terribly wrong...
#endif

// Force (kindly ask) the compiler to inline the function.
#if STORM_COMPILER_MSVC_
#define STORM_FORCE_INLINE_ inline __forceinline
#else
#define STORM_FORCE_INLINE_ inline __attribute__((always_inline))
#endif

// Cross-compiler version of `[[no_unique_address]]`.
#if STORM_COMPILER_MSVC_
#define STORM_NO_UNIQUE_ADDRESS_ [[msvc::no_unique_address]]
#else
#define STORM_NO_UNIQUE_ADDRESS_ [[no_unique_address]]
#endif

// Assume the expression is always true.
#if STORM_COMPILER_MSVC_
#define STORM_ASSUME_(expression) __assume(expression)
#else
#define STORM_ASSUME_(expression)                   \
  do {                                              \
    if (!(expression)) { __builtin_unreachable(); } \
  } while (false)
#endif

// Report a trace message.
#define STORM_TRACE_(message, ...) \
  (spdlog::trace(message __VA_OPT__(, __VA_ARGS__)))

// Report a debug message.
#define STORM_DEBUG_(message, ...) \
  (spdlog::debug(message __VA_OPT__(, __VA_ARGS__)))

// Report an info message.
#define STORM_INFO_(message, ...) \
  (spdlog::info(message __VA_OPT__(, __VA_ARGS__)))

// Report a warning.
#define STORM_WARNING_(message, ...) \
  (spdlog::warn(message __VA_OPT__(, __VA_ARGS__)))

// Report an error.
#define STORM_ERROR_(message, ...) \
  (spdlog::error(message __VA_OPT__(, __VA_ARGS__)))

// Report a critical error.
#define STORM_CRITICAL_(message, ...) \
  (spdlog::critical(message __VA_OPT__(, __VA_ARGS__)))

#ifndef __FUNCSIG__
#define __FUNCSIG__ __PRETTY_FUNCTION__
#endif

// Report a fatal error and exit the application.
#define STORM_FATAL_ERROR_(message, ...)                             \
  do {                                                               \
    STORM_CRITICAL_("Fatal error:");                                 \
    STORM_CRITICAL_("{}:{} {}: {}", __FILE__, __LINE__, __FUNCSIG__, \
                    fmt::format(message __VA_OPT__(, __VA_ARGS__))); \
    spdlog::shutdown(), std::abort();                                \
  } while (false)

// Check the expression, exit with fatal error if it fails.
#define STORM_ENSURE_(expression, message, ...)              \
  do {                                                       \
    if (!(expression)) {                                     \
      STORM_CRITICAL_("\"{}\" failed!", #expression);        \
      STORM_FATAL_ERROR_(message __VA_OPT__(, __VA_ARGS__)); \
    }                                                        \
  } while (false)

// Check the expression in the debug mode, exit with fatal error if it fails.
#ifdef NDEBUG
#define STORM_ASSERT_(expression, message, ...) STORM_ASSUME_(expression)
#else
#define STORM_ASSERT_(expression, message, ...) \
  STORM_ENSURE_(expression, message __VA_OPT__(, __VA_ARGS__))
#endif

namespace Storm {

/// Library for the math functions.
namespace math {}

/// A small metaprogramming library.
namespace meta {}

/// @brief Contains the internal implementation details.
namespace detail_ {}

/// @brief Size type.
using size_t = std::size_t;

/// @brief Pointer difference type.
using ptrdiff_t = std::ptrdiff_t;

/// @brief Real floating-point type.
using real_t = double;

namespace detail_ {

  template<class>
  constexpr inline bool always_false_ = false;

  constexpr bool in_range_(auto t, auto min, auto max) {
    return min <= t && t <= max;
  }

  template<class T1, class T2>
  concept different_from_ =
      !std::same_as<std::remove_cvref_t<T1>, std::remove_cvref_t<T2>>;

} // namespace detail_

/// @brief @c size_t compile-time constant.
template<size_t N>
using size_t_constant = std::integral_constant<size_t, N>;

} // namespace Storm
