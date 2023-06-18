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

#include <Storm/Crow/Base/Config.hpp>

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <cstdlib>
#include <utility>

// -----------------------------------------------------------------------------

// Report a fatal error and exit the application.
#define STORM_ABORT(message, ...)                                           \
  (STORM_CRITICAL("{}:{} {}: Critical error! " message, __FILE__, __LINE__, \
                  STORM_CURRENT_FUNCTION_STRING __VA_OPT__(, __VA_ARGS__)), \
   std::abort())

// Check the expression, exit with fatal error if it fails.
#define STORM_ENSURE(expression, message, ...)              \
  (STORM_LIKELY(static_cast<bool>((expression))) ?          \
       void(0) :                                            \
       STORM_ABORT("Condition `{}` was violated! " message, \
                   #expression __VA_OPT__(, __VA_ARGS__)))

// Enable or disable asserts.
#ifndef STORM_ENABLE_ASSERTS
#  if STORM_NDEBUG
#    define STORM_ENABLE_ASSERTS 0
#  else
#    define STORM_ENABLE_ASSERTS 1
#  endif
#endif

// If enabled, check the expression, exit with fatal error if it fails.
#if STORM_ENABLE_ASSERTS
#  define STORM_ASSERT(expression, message, ...) \
    STORM_ENSURE(expression, message __VA_OPT__(, __VA_ARGS__))
#else
#  define STORM_ASSERT(expression, message, ...) STORM_ASSUME(expression)
#endif

// -----------------------------------------------------------------------------
