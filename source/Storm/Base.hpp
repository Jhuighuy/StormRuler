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
// FITNESS FOR Allocator PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT
// SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
// OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

#pragma once

#include <complex>
#include <concepts>
#include <cstdint>
#include <cstdlib>
#include <type_traits>
#include <utility>

#include <fmt/format.h>
#include <spdlog/spdlog.h>

// -----------------------------------------------------------------------------

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
#ifdef __GNUC__
#ifdef __clang__
#define STORM_COMPILER_CLANG_ 1
#define STORM_COMPILER_GCC_ 0
#else
#define STORM_COMPILER_GCC_ 1
#define STORM_COMPILER_CLANG_ 0
#endif
#endif
#ifdef _MSC_VER
#define STORM_COMPILER_MSVC_ 1
#else
#define STORM_COMPILER_MSVC_ 0
#endif
#if (STORM_COMPILER_GCC_ + STORM_COMPILER_CLANG_ + STORM_COMPILER_MSVC_) != 1
#warning Storm has detected multiple compilers, something is terribly wrong!
#endif

// C++23 constexpr.
#if STORM_CPP23_
#define STORM_CPP23_CONSTEXPR_ constexpr
#else
#define STORM_CPP23_CONSTEXPR_
#endif

// C++23 allows `static` inside of the `constexpr` functions.
#if STORM_CPP23_
#define STORM_CPP23_STATIC_ static
#else
#define STORM_CPP23_STATIC_
#endif

// C++23 allows `thread_local` inside of the `constexpr` functions.
#if STORM_CPP23_
#define STORM_CPP23_THREAD_LOCAL_ thread_local
#else
#define STORM_CPP23_THREAD_LOCAL_
#endif

// C++23 `std::unreachable()`.
#if STORM_CPP23_
#define STORM_UNREACHABLE_() std::unreachable()
#else
#if STORM_COMPILER_MSVC_
#define STORM_UNREACHABLE_() __assume(false)
#else
#define STORM_UNREACHABLE_() __builtin_unreachable()
#endif
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

// Include the configuration file.
#if __has_include("./Config.hpp")
#include "./Config.hpp"
#endif
#include "./ConfigDefault.hpp"

// Is the code parsed by Doxygen?
#ifndef STORM_DOXYGEN_
#define STORM_DOXYGEN_ 0
#endif

// -----------------------------------------------------------------------------

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

// Assume the expression is always true.
#if STORM_COMPILER_MSVC_
#define STORM_ASSUME_(expression) __assume(expression)
#else
#define STORM_ASSUME_(expression)               \
  do {                                          \
    if (!(expression)) __builtin_unreachable(); \
  } while (false)
#endif

#ifndef __FUNCSIG__
#define __FUNCSIG__ static_cast<const char*>(__PRETTY_FUNCTION__)
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
#define STORM_ENSURE_(expression, message, ...)                         \
  do {                                                                  \
    if (!(expression)) {                                                \
      STORM_TERMINATE_("Condition \"{}\" failed! {}", #expression,      \
                       fmt::format(message __VA_OPT__(, __VA_ARGS__))); \
    }                                                                   \
  } while (false)

// Enable or disable asserts.
#ifndef STORM_ENABLE_ASSERTS_
#ifdef NDEBUG
#define STORM_ENABLE_ASSERTS_ 0
#else
#define STORM_ENABLE_ASSERTS_ 1
#endif
#endif

// Check the expression, exit with fatal error if it fails.
#ifdef NDEBUG
#define STORM_ASSERT_(expression, message, ...) STORM_ASSUME_(expression)
#else
#define STORM_ASSERT_(expression, message, ...) \
  STORM_ENSURE_(expression, message __VA_OPT__(, __VA_ARGS__))
#endif

// -----------------------------------------------------------------------------

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

// -----------------------------------------------------------------------------

namespace Storm
{

// -----------------------------------------------------------------------------

/// A small metaprogramming library.
namespace meta
{
}

/// Shapes collection.
namespace shapes
{
}

/// @brief Contains the internal implementation details.
namespace detail_
{
}

// -----------------------------------------------------------------------------

/// @brief Size type.
using size_t = std::size_t;

/// @brief Pointer difference type.
using ptrdiff_t = std::ptrdiff_t;

/// @brief Real floating-point type.
using real_t = double;

/// @brief Complex floating-point type.
using complex_t = std::complex<real_t>;

/// @brief size_t literal.
constexpr size_t operator""_sz(unsigned long long arg) noexcept // NOLINT
{
  return static_cast<size_t>(arg);
}

/// @brief size_t constant.
template<size_t Arg>
using size_t_constant = std::integral_constant<size_t, Arg>;

/// @brief real_t literal.
constexpr real_t operator""_dp(long double arg) noexcept
{
  return static_cast<real_t>(arg);
}

// -----------------------------------------------------------------------------

/// @brief Non movable (and non copyable) object interface.
class NonMovable
{
protected:

  /// @brief Non movable object is default initializible.
  NonMovable() = default;

  /// @brief Non copyable object is NOT move constructible.
  NonMovable(NonMovable&&) = default;

  /// @brief Non copyable object is NOT copy constructible.
  NonMovable(const NonMovable&) = delete;

  /// @brief Non copyable object is NOT move assignable.
  NonMovable& operator=(NonMovable&&) = delete;

  /// @brief Non copyable object is NOT copy assignable.
  NonMovable& operator=(const NonMovable&) = delete;

  /// @brief Non movable object is destructible.
  ~NonMovable() = default;

}; // class NonMovable

/// @brief Non copyable object interface.
class NonCopyable
{
protected:

  /// @brief Non copyable object is default initializible.
  NonCopyable() = default;

  /// @brief Non copyable object is move constructible.
  NonCopyable(NonCopyable&&) = default;

  /// @brief Non copyable object is NOT copy constructible.
  NonCopyable(const NonCopyable&) = delete;

  /// @brief Non copyable object is move assignable.
  NonCopyable& operator=(NonCopyable&&) = default;

  /// @brief Non copyable object is NOT copy assignable.
  NonCopyable& operator=(const NonCopyable&) = delete;

  /// @brief Non copyable object is destructible.
  ~NonCopyable() = default;

}; // class NonCopyable

/// @brief Non assignable object interface.
class NonAssignable
{
protected:

  /// @brief Non assignable object is default initializible.
  NonAssignable() = default;

  /// @brief Non assignable object is move constructible.
  NonAssignable(NonAssignable&&) = default;

  /// @brief Non assignable object is copy constructible.
  NonAssignable(const NonAssignable&) = default;

  /// @brief Non assignable object is NOT move assignable.
  NonAssignable& operator=(NonAssignable&&) = delete;

  /// @brief Non assignable object is NOT copy assignable.
  NonAssignable& operator=(const NonAssignable&) = delete;

  /// @brief Non assignable object is destructible.
  ~NonAssignable() = default;

}; // class NonAssignable

// -----------------------------------------------------------------------------

namespace detail_
{

  using noncopyable_ = NonCopyable;

  constexpr bool in_range_(auto t, auto min, auto max)
  {
    return min <= t && t <= max;
  }

  constexpr bool one_of_(auto x, auto... vals)
  {
    return ((x == vals) || ...);
  }

  template<class T1, class T2>
  concept different_from_ =
      (!std::same_as<std::remove_cvref_t<T1>, std::remove_cvref_t<T2>>);

} // namespace detail_

// -----------------------------------------------------------------------------

using Error = std::runtime_error;
using IoError = std::runtime_error;
using GlError = std::runtime_error;

// -----------------------------------------------------------------------------

} // namespace Storm
