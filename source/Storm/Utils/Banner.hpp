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

#include <Storm/Base.hpp>

#include <cstdio>

namespace Storm
{

// -----------------------------------------------------------------------------

/// @brief Print the epic banner.
static void print_banner()
{
  // clang-format off
  std::puts("\n");
  std::puts(R"(  ╔═══════════════════════════════════════════════════════════════════╗  )");
  std::puts(R"(  ║ ╔═══════════════════════════════════════════════════════════════╗ ║  )");
  std::puts(R"(  ║ ║      _____ __                       ____        __            ║ ║  )");
  std::puts(R"(  ║ ║     / ___// /_____   ____ ___ ___  / __ \__  __/ ___  _____   ║ ║  )");
  std::puts(R"(  ║ ║     \__ \/ __/`__ \/ ___/`__ `__ \/ /_/ / / / / / _ \/ ___/   ║ ║  )");
  std::puts(R"(  ║ ║    ___/ / /_/ /_/ / /  / / / / / / _, _/ /_/ / /  __/ /       ║ ║  )");
  std::puts(R"(  ║ ║   /____/\__/\____/_/  /_/ /_/ /_/_/ |_|\__,_/_/\___/_/        ║ ║  )");
  std::puts(R"(  ║ ║                                                               ║ ║  )");
  std::puts(R"(  ║ ╚═══════════════════════════════════════════════════════════════╝ ║  )");
  std::puts(R"(  ╚═══════════════════════════════════════════════════════════════════╝  )");
  std::puts("\n");
  // clang-format on
  STORM_INFO_("Built on {}, {}.", __DATE__, __TIME__);
}

// -----------------------------------------------------------------------------

} // namespace Storm
