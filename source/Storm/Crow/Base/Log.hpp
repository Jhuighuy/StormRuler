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

#include <fmt/format.h>
#include <spdlog/spdlog.h>

// -----------------------------------------------------------------------------

// Format a string.
#define STORM_FORMAT(message, ...) \
  (fmt::format(message __VA_OPT__(, __VA_ARGS__)))

// Report a trace message.
#define STORM_TRACE(message, ...) \
  (spdlog::trace(message __VA_OPT__(, __VA_ARGS__)))

// Report a debug message.
#define STORM_DEBUG(message, ...) \
  (spdlog::debug(message __VA_OPT__(, __VA_ARGS__)))

// Report an info message.
#define STORM_INFO(message, ...) \
  (spdlog::info(message __VA_OPT__(, __VA_ARGS__)))

// Report a warning.
#define STORM_WARNING(message, ...) \
  (spdlog::warn(message __VA_OPT__(, __VA_ARGS__)))

// Report an error.
#define STORM_ERROR(message, ...) \
  (spdlog::error(message __VA_OPT__(, __VA_ARGS__)))

// Report a critical error.
#define STORM_CRITICAL(message, ...) \
  (spdlog::critical(message __VA_OPT__(, __VA_ARGS__)))

// -----------------------------------------------------------------------------
