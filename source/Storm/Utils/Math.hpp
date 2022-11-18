// Copyright (C) 2020 - 2023 Oleg Butakov
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

#include <cmath>
#include <complex>
#include <concepts>
#include <type_traits>

namespace Storm {

template<class>
inline constexpr bool is_complex_floating_point_v = false;
template<class T>
inline constexpr bool is_complex_floating_point_v<std::complex<T>> =
    std::is_floating_point_v<T>;

/// @brief Complex floating-point numbers.
template<class T>
concept complex_floating_point = is_complex_floating_point_v<T>;

/// @brief Real or complex floaing point numbers.
template<class T>
concept real_or_complex_floating_point = std::floating_point<T> || //
                                         complex_floating_point<T>;

// -----------------------------------------------------------------------------

static real_t deg2rad(real_t) noexcept;
static real_t rad2deg(real_t) noexcept;

// -----------------------------------------------------------------------------

/// @brief If @p y is zero, return zero,
/// else return value of @p x divided by @p y.
template<real_or_complex_floating_point Value>
auto safe_divide(Value x, Value y) {
  static constexpr Value zero{0.0};
  return y == zero ? zero : (x / y);
}

// -----------------------------------------------------------------------------

/// @brief Conjugate of complex number @p x.
/// @{
[[nodiscard]] constexpr real_t conj(real_t x) noexcept {
  return x;
}
[[nodiscard]] constexpr complex_t conj(const complex_t& x) noexcept {
  return std::conj(x);
}
/// @}

using std::real;

using std::imag;

// -----------------------------------------------------------------------------

using std::abs;

template<class Real>
constexpr Real sign(Real) noexcept;

using std::min;

using std::max;

// -----------------------------------------------------------------------------

using std::exp;

using std::exp2;

using std::log;

using std::log2;

using std::log10;

// -----------------------------------------------------------------------------

using std::pow;

using std::sqrt;

using std::cbrt;

using std::hypot;

// -----------------------------------------------------------------------------

using std::sin;

using std::cos;

using std::tan;

using std::asin;

using std::acos;

using std::atan;

using std::atan2;

// -----------------------------------------------------------------------------

using std::sinh;

using std::cosh;

using std::tanh;

using std::asinh;

using std::acosh;

using std::atanh;

// -----------------------------------------------------------------------------

/// @brief Generate the Givens rotation.
template<real_or_complex_floating_point Value>
auto sym_ortho(Value a, Value b) {
  // Compute:
  // ----------------------
  // 𝑟𝑟 ← (𝑎² + 𝑏²)¹ᐟ²,
  // 𝑐𝑠 ← 𝑎/𝑟𝑟, 𝑠𝑛 ← 𝑏/𝑟𝑟.
  // ----------------------
  Value cs, sn, rr;
  rr = std::hypot(a, b);
  if (rr > Value{0.0}) {
    cs = a / rr;
    sn = b / rr;
  } else {
    cs = 1.0;
    sn = 0.0;
  }

  return std::tuple(cs, sn, rr);

} // sym_ortho

} // namespace Storm
