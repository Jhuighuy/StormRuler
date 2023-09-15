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

#include <algorithm>
#include <cmath>
#include <complex>
#include <concepts>
#include <numbers>
#include <tuple>

namespace Storm {

// -----------------------------------------------------------------------------

/// @brief Floored modulo function.
template<std::integral Dividend, std::integral Divisor>
constexpr auto floored_mod(Dividend dividend, Divisor divisor) noexcept {
  if constexpr (std::unsigned_integral<Dividend> &&
                std::unsigned_integral<Divisor>) {
    return dividend % divisor;
  } else {
    return ((dividend % divisor) + divisor) % divisor;
  }
}

/// @brief If @p y is zero, return zero,
/// otherwise return value of @p x divided by @p y.
template<std::floating_point Arg>
constexpr auto safe_divide(Arg x, Arg y) noexcept {
  return (y == Arg{0}) ? Arg{0} : (x / y);
}

/// @brief If @p y is zero, return zero, otherwise return onverse @p y.
template<std::floating_point Arg>
constexpr auto safe_inverse(Arg y) noexcept {
  return safe_divide(Arg{1}, y);
}

// -----------------------------------------------------------------------------

/// @brief @f$ \pi @f$ constant.
inline constexpr real_t pi = std::numbers::pi_v<real_t>;

/// @brief Convert degrees @p degrees to radians.
template<std::floating_point Arg>
constexpr auto deg2rad(Arg degrees) noexcept {
  constexpr auto pi_over_180 = std::numbers::pi_v<Arg> / Arg{180};
  return pi_over_180 * degrees;
}

/// @brief Convert radians @p radians to degrees.
template<std::floating_point Arg>
constexpr auto rad2deg(Arg radians) noexcept {
  constexpr auto _180_over_pi = Arg{180} / std::numbers::pi_v<Arg>;
  return _180_over_pi * radians;
}

// -----------------------------------------------------------------------------

using std::conj;
using std::imag;
using std::real;

/// @brief @f$ i = \sqrt{-1} @f$.
inline constexpr complex_t i{0.0_dp, 1.0_dp};

/// @brief Dot product of two scalars.
/// @{
template<class Arg1, class Arg2>
  requires (std::integral<Arg1> || std::floating_point<Arg1>) ||
           (std::integral<Arg2> || std::floating_point<Arg2>)
constexpr auto dot_product(Arg1 arg1, Arg2 arg2) noexcept {
  return arg1 * arg2;
}
template<class Arg1, class Arg2>
constexpr auto dot_product(Arg1 arg1, std::complex<Arg2> arg2) noexcept {
  return arg1 * std::conj(arg2);
}
/// @}

// -----------------------------------------------------------------------------

using std::abs;

/// @brief Sign of the number @p arg.
template<class Arg>
  requires std::signed_integral<Arg> || std::floating_point<Arg>
constexpr int sign(Arg arg) noexcept {
  return (Arg{0} < arg) - (arg < Arg{0});
}

// -----------------------------------------------------------------------------

using std::max;
using std::min;

// -----------------------------------------------------------------------------

using std::cbrt;
using std::hypot;
using std::pow;
using std::sqrt;

/// @brief @f$ \sqrt{2} @f$ constant.
inline constexpr real_t sqrt2 = std::numbers::sqrt2_v<real_t>;

/// @brief @f$ \sqrt{3} @f$ constant.
inline constexpr real_t sqrt3 = std::numbers::sqrt3_v<real_t>;

// -----------------------------------------------------------------------------

using std::exp;
using std::exp2;
using std::log;
using std::log10;
using std::log2;

/// @brief @f$ e @f$ constant.
inline constexpr real_t e = std::numbers::e_v<real_t>;

// -----------------------------------------------------------------------------

using std::acos;
using std::asin;
using std::atan;
using std::atan2;
using std::cos;
using std::sin;
using std::tan;

// -----------------------------------------------------------------------------

using std::acosh;
using std::asinh;
using std::atanh;
using std::cosh;
using std::sinh;
using std::tanh;

// -----------------------------------------------------------------------------

/// @brief Generate the Givens rotation.
template<class Value>
constexpr auto sym_ortho(Value a, Value b) noexcept {
  // Compute:
  // ----------------------
  // 𝑟𝑟 ← (𝑎² + 𝑏²)¹ᐟ²,
  // 𝑐𝑠 ← 𝑎/𝑟𝑟, 𝑠𝑛 ← 𝑏/𝑟𝑟.
  // ----------------------
  Value cs, sn, rr;
  rr = hypot(a, b);
  if (rr > Value{0}) {
    cs = a / rr, sn = b / rr;
  } else {
    cs = Value{1}, sn = Value{0};
  }
  return std::tuple{cs, sn, rr};
}

// -----------------------------------------------------------------------------

} // namespace Storm
