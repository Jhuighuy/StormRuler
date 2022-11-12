// Copyright © 2020 - 2023 Oleg Butakov
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

#include <Storm/Base.hpp>

#include <fmt/format.h>

namespace Storm {

/// @brief Print the epic banner.
static void print_banner() {
  // clang-format off
  fmt::print("\n");
  fmt::print(R"(  ╔═══════════════════════════════════════════════════════════════════╗)" "\n");
  fmt::print(R"(  ║ ╔═══════════════════════════════════════════════════════════════╗ ║)" "\n");
  fmt::print(R"(  ║ ║      _____ __                       ____        __            ║ ║)" "\n");
  fmt::print(R"(  ║ ║     / ___// /_____   ____ ___ ___  / __ \__  __/ ___  _____   ║ ║)" "\n");
  fmt::print(R"(  ║ ║     \__ \/ __/`__ \/ ___/`__ `__ \/ /_/ / / / / / _ \/ ___/   ║ ║)" "\n");
  fmt::print(R"(  ║ ║    ___/ / /_/ /_/ / /  / / / / / / _, _/ /_/ / /  __/ /       ║ ║)" "\n");
  fmt::print(R"(  ║ ║   /____/\__/\____/_/  /_/ /_/ /_/_/ |_|\__,_/_/\___/_/        ║ ║)" "\n");
  fmt::print(R"(  ║ ║                                                               ║ ║)" "\n");
  fmt::print(R"(  ║ ╚═══════════════════════════════════════════════════════════════╝ ║)" "\n");
  fmt::print(R"(  ╚═══════════════════════════════════════════════════════════════════╝)" "\n");
  fmt::print("\n");
  STORM_INFO_("Built on {}, {}.", __DATE__, __TIME__);
  // clang-format on
}

} // namespace Storm
