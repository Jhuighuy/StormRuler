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

#include <concepts>
#include <cstdint>
#include <cstdlib>
#include <type_traits>

#include <fmt/format.h>
#include <spdlog/spdlog.h>

// Detect the C++ version.
#if __cplusplus < 202002L && !defined(_MSC_VER)
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

// C++23 allows `thread_local` inside of the `constexpr` functions.
#if STORM_CPP23_
#define STORM_CPP23_THREAD_LOCAL_ thread_local
#else
#define STORM_CPP23_THREAD_LOCAL_
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
#define STORM_TERMINATE_(message, ...)                               \
  do {                                                               \
    STORM_CRITICAL_("Fatal error:");                                 \
    STORM_CRITICAL_("{}:{} {}: {}", __FILE__, __LINE__, __FUNCSIG__, \
                    fmt::format(message __VA_OPT__(, __VA_ARGS__))); \
    spdlog::shutdown(), std::abort();                                \
  } while (false)

// Check the expression, exit with fatal error if it fails.
#define STORM_ENSURE_(expression, message, ...)            \
  do {                                                     \
    if (!(expression)) {                                   \
      STORM_CRITICAL_("\"{}\" failed!", #expression);      \
      STORM_TERMINATE_(message __VA_OPT__(, __VA_ARGS__)); \
    }                                                      \
  } while (false)

// Check the expression in the debug mode, exit with fatal error if it fails.
#ifdef NDEBUG
#define STORM_ASSERT_(expression, message, ...) STORM_ASSUME_(expression)
#else
#define STORM_ASSERT_(expression, message, ...) \
  STORM_ENSURE_(expression, message __VA_OPT__(, __VA_ARGS__))
#endif

// Include the configuration file.
#if __has_include("Config.hpp")
#include "Config.hpp"
#endif
#include "ConfigDefault.hpp"

#define STORM_THROW_BASE_(error, message, ...)                   \
  do {                                                           \
    throw error(fmt::format(message __VA_OPT__(, __VA_ARGS__))); \
  } while (false)

#define STORM_THROW_(message, ...) \
  STORM_THROW_BASE_(Storm::Error, message __VA_OPT__(, __VA_ARGS__))

#define STORM_THROW_IO_(message, ...) \
  STORM_THROW_BASE_(Storm::IoError, message __VA_OPT__(, __VA_ARGS__))

#define STORM_THROW_GL_(message, ...) \
  STORM_THROW_BASE_(Storm::GlError, message __VA_OPT__(, __VA_ARGS__))

namespace Storm {

/// Library for the math functions.
namespace math {}

/// A small metaprogramming library.
namespace meta {}

/// Shapes collection.
namespace shapes {}

/// @brief Contains the internal implementation details.
namespace detail_ {}

/// @brief Size type.
using size_t = std::size_t;

/// @brief Pointer difference type.
using ptrdiff_t = std::ptrdiff_t;

/// @brief Real floating-point type.
using real_t = double;

namespace detail_ {

  struct noncopyable_ {
    noncopyable_() = default;
    noncopyable_(noncopyable_&&) = default;
    noncopyable_& operator=(noncopyable_&&) = default;
    noncopyable_(const noncopyable_&) = delete;
    noncopyable_& operator=(const noncopyable_&) = delete;
  }; // struct noncopyable_

  [[nodiscard]] constexpr bool in_range_(auto t, auto min, auto max) {
    return min <= t && t <= max;
  }

  [[nodiscard]] constexpr bool one_of_(auto x, auto... vals) {
    return ((x == vals) || ...);
  }

  template<class T1, class T2>
  concept different_from_ =
      (!std::same_as<std::remove_cvref_t<T1>, std::remove_cvref_t<T2>>);

} // namespace detail_

/// @brief @c size_t compile-time constant.
template<size_t N>
using size_t_constant = std::integral_constant<size_t, N>;

/// @brief Check if type is a @c size_t constant.
/// @{
template<class>
inline constexpr bool is_size_t_constant_v = false;
template<size_t N>
inline constexpr bool is_size_t_constant_v<size_t_constant<N>> = true;
/// @}

using Error = std::runtime_error;
using IoError = std::runtime_error;
using GlError = std::runtime_error;

} // namespace Storm
