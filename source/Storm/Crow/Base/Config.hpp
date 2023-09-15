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

// -----------------------------------------------------------------------------

// Detect the C++ compiler.
#ifdef _MSC_VER
#  define STORM_COMPILER_MSVC 1
#else
#  define STORM_COMPILER_MSVC 0
#endif
#ifdef __clang__
#  define STORM_COMPILER_CLANG 1
#else
#  define STORM_COMPILER_CLANG 0
#endif
#if defined(__GNUC__) && !STORM_COMPILER_CLANG
#  define STORM_COMPILER_GCC 1
#else
#  define STORM_COMPILER_GCC 0
#endif

#if (STORM_COMPILER_GCC + STORM_COMPILER_CLANG + STORM_COMPILER_MSVC) != 1
#  warning Storm has detected multiple compilers, something is terribly wrong!
#endif

// -----------------------------------------------------------------------------

// Detect the CPU architecture.
#if STORM_COMPILER_MSVC
#  ifdef _M_X64
#    define STORM_ARCH_X64 1
#  else
#    define STORM_ARCH_X64 0
#  endif
#  ifdef _M_ARM64
#    define STORM_ARCH_ARM64 1
#  else
#    define STORM_ARCH_ARM64 0
#  endif
#elif STORM_COMPILER_CLANG || STORM_COMPILER_GCC
#  ifdef __amd64__
#    define STORM_ARCH_X64 1
#  else
#    define STORM_ARCH_X64 0
#  endif
#  ifdef __aarch64__
#    define STORM_ARCH_ARM64 1
#  else
#    define STORM_ARCH_ARM64 0
#  endif
#endif

#if (STORM_ARCH_X64 + STORM_ARCH_ARM64) != 1
#  warning Storm has detected multiple CPU architectures, \
           something is terribly wrong!
#endif

// -----------------------------------------------------------------------------

// Detect the C++ version.
#if STORM_COMPILER_MSVC
#  define STORM_CPLUSPLUS _MSVC_LANG
#else
#  define STORM_CPLUSPLUS __cplusplus
#endif

#if STORM_CPLUSPLUS < 202002L
#  error Storm requires C++20 support!
#endif
#if STORM_CPLUSPLUS > 202002L
#  define STORM_CPP23 1
#else
#  define STORM_CPP23 0
#endif

// C++23 constexpr.
#if STORM_CPP23
#  define STORM_CPP23_CONSTEXPR constexpr
#else
#  define STORM_CPP23_CONSTEXPR
#endif

// C++23 allows `static` inside of the `constexpr` functions.
#if STORM_CPP23
#  define STORM_CPP23_STATIC static
#else
#  define STORM_CPP23_STATIC
#endif

// C++23 allows `thread_local` inside of the `constexpr` functions.
#if STORM_CPP23
#  define STORM_CPP23_THREAD_LOCAL thread_local
#else
#  define STORM_CPP23_THREAD_LOCAL
#endif

// C++23 `std::unreachable()`.
#if STORM_CPP23
#  define STORM_UNREACHABLE() std::unreachable()
#else
#  if STORM_COMPILER_MSVC
#    define STORM_UNREACHABLE() __assume(false)
#  elif STORM_COMPILER_CLANG || STORM_COMPILER_GCC
#    define STORM_UNREACHABLE() __builtin_unreachable()
#  endif
#endif

// C++23 `[[assume(..)]]`.
// Double parentheses are needed here for fold expressions to work.
#if STORM_CPP23
#  define STORM_ASSUME(expression) [[assume((expression))]]
#else
#  if STORM_COMPILER_MSVC
#    define STORM_ASSUME(expression) __assume((expression))
#  elif STORM_COMPILER_CLANG || STORM_COMPILER_GCC
#    define STORM_ASSUME(expression) \
      (static_cast<bool>((expression)) ? void(0) : __builtin_unreachable())
#  endif
#endif

// Force (kindly ask) the compiler to inline the function.
#if STORM_COMPILER_MSVC
#  define STORM_FORCE_INLINE inline __forceinline
#elif STORM_COMPILER_CLANG || STORM_COMPILER_GCC
#  define STORM_FORCE_INLINE inline __attribute__((always_inline))
#endif

// Cross-compiler version of `[[no_unique_address]]`.
#if STORM_COMPILER_MSVC
#  define STORM_NO_UNIQUE_ADDRESS [[msvc::no_unique_address]]
#else
#  define STORM_NO_UNIQUE_ADDRESS [[no_unique_address]]
#endif

// Current function name.
#if STORM_COMPILER_MSVC
#  define STORM_CURRENT_FUNCTION __FUNCSIG__
#elif STORM_COMPILER_CLANG || STORM_COMPILER_GCC
#  define STORM_CURRENT_FUNCTION __PRETTY_FUNCTION__
#else
#  define STORM_CURRENT_FUNCTION __FUNCTION__
#endif
#define STORM_CURRENT_FUNCTION_STRING \
  static_cast<const char*>(STORM_CURRENT_FUNCTION)

// Likely/unlikely for ternary operators.
#if STORM_COMPILER_CLANG || STORM_COMPILER_GCC
#  define STORM_LIKELY(expression) __builtin_expect((expression), 1)
#  define STORM_UNLIKELY(expression) __builtin_expect((expression), 0)
#else
#  define STORM_LIKELY(expression) (expression)
#  define STORM_UNLIKELY(expression) (expression)
#endif

// -----------------------------------------------------------------------------

// Include the configuration file.
#if __has_include(<Storm/Config.hpp>)
#  include <Storm/Config.hpp>
#endif
#include <Storm/ConfigDefault.hpp>

// -----------------------------------------------------------------------------

// Is the code parsed by Doxygen?
#ifndef STORM_DOXYGEN
#  define STORM_DOXYGEN 0
#endif

// Is coverage enables?
#ifndef STORM_COVERAGE
#  define STORM_COVERAGE 0
#endif

// Are optimizations enabled?
#ifdef NDEBUG
#  define STORM_NDEBUG 1
#else
#  define STORM_NDEBUG 0
#endif

// -----------------------------------------------------------------------------
