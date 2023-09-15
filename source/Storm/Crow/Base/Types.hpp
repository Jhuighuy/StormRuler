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

#include <complex>
#include <cstddef>
#include <type_traits>

namespace Storm {

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
constexpr size_t operator""_sz(unsigned long long arg) noexcept { // NOLINT
  return static_cast<size_t>(arg);
}

/// @brief size_t constant.
template<size_t Arg>
using fixed_size_t = std::integral_constant<size_t, Arg>;

/// @brief real_t literal.
constexpr real_t operator""_dp(long double arg) noexcept {
  return static_cast<real_t>(arg);
}

// -----------------------------------------------------------------------------

/// @brief Universal placeholder type.
class Underscore final {
public:

  /// @brief Construct the placeholder.
  explicit Underscore() = default;

}; // class Underscore

/// @brief Universal placeholder.
inline constexpr Underscore _{};

// -----------------------------------------------------------------------------

/// @brief Non movable (and non copyable) object interface.
class NonMovableInterface {
protected:

  NonMovableInterface() = default;
  NonMovableInterface(NonMovableInterface&&) = default;
  NonMovableInterface(const NonMovableInterface&) = delete;
  NonMovableInterface& operator=(NonMovableInterface&&) = delete;
  NonMovableInterface& operator=(const NonMovableInterface&) = delete;
  ~NonMovableInterface() = default;

}; // class NonMovableInterface

/// @brief Non copyable object interface.
class NonCopyableInterface {
protected:

  NonCopyableInterface() = default;
  NonCopyableInterface(NonCopyableInterface&&) = default;
  NonCopyableInterface(const NonCopyableInterface&) = delete;
  NonCopyableInterface& operator=(NonCopyableInterface&&) = default;
  NonCopyableInterface& operator=(const NonCopyableInterface&) = delete;
  ~NonCopyableInterface() = default;

}; // class NonCopyableInterface

/// @brief Non assignable object interface.
class NonAssignableInterface {
protected:

  NonAssignableInterface() = default;
  NonAssignableInterface(NonAssignableInterface&&) = default;
  NonAssignableInterface(const NonAssignableInterface&) = default;
  NonAssignableInterface& operator=(NonAssignableInterface&&) = delete;
  NonAssignableInterface& operator=(const NonAssignableInterface&) = delete;
  ~NonAssignableInterface() = default;

}; // class NonAssignableInterface

// -----------------------------------------------------------------------------

} // namespace Storm
