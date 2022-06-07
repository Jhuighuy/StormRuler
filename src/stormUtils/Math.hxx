/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// Copyright (C) 2022 Oleg Butakov
///
/// Permission is hereby granted, free of charge, to any person
/// obtaining a copy of this software and associated documentation
/// files (the "Software"), to deal in the Software without
/// restriction, including without limitation the rights  to use,
/// copy, modify, merge, publish, distribute, sublicense, and/or
/// sell copies of the Software, and to permit persons to whom the
/// Software is furnished to do so, subject to the following
/// conditions:
///
/// The above copyright notice and this permission notice shall be
/// included in all copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
/// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
/// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
/// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
/// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
/// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
/// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
/// OTHER DEALINGS IN THE SOFTWARE.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///

#pragma once

#include <stormBase.hxx>

namespace Storm::math {

/// @name Basic operations.
/// @{

constexpr auto abs(auto arg) noexcept {
  return std::abs(arg);
}

/// @brief If @p y is zero, return zero,
///   else return value of @p x divided by @p y.
template<real_or_complex_floating_point Value>
auto safe_divide(Value x, Value y) {
  static constexpr Value zero{0.0};
  return y == zero ? zero : (x / y);
}

/// @}

/// @name Exponential functions.
/// @{

constexpr auto exp(auto x) noexcept {
  return std::exp(x);
}

constexpr auto exp2(auto x) noexcept {
  return std::exp2(x);
}

constexpr auto log(auto x) noexcept {
  return std::log(x);
}

constexpr auto log2(auto x) noexcept {
  return std::log2(x);
}

constexpr auto log10(auto x) noexcept {
  return std::log10(x);
}

/// @}

/// @name Power functions.
/// @{

constexpr auto pow(auto x, auto y) noexcept {
  return std::pow(x, y);
}

constexpr auto sqrt(auto x) noexcept {
  return std::sqrt(x);
}

constexpr auto cbrt(auto x) noexcept {
  return std::cbrt(x);
}

constexpr auto hypot(auto x, auto y) noexcept {
  return std::hypot(x, y);
}

constexpr auto hypot(auto x, auto y, auto z) noexcept {
  return std::hypot(x, y, z);
}

/// @}

/// @name Trigonometric functions.
/// @{

constexpr auto sin(auto x) noexcept {
  return std::sin(x);
}

constexpr auto cos(auto x) noexcept {
  return std::cos(x);
}

constexpr auto tan(auto x) noexcept {
  return std::tan(x);
}

constexpr auto asin(auto x) noexcept {
  return std::asin(x);
}

constexpr auto acos(auto x) noexcept {
  return std::acos(x);
}

constexpr auto atan(auto x) noexcept {
  return std::atan(x);
}

constexpr auto atan2(auto y, auto x) noexcept {
  return std::atan2(y, x);
}

/// @}

/// @name Hyperbolic functions.
/// @{

constexpr auto sinh(auto x) noexcept {
  return std::sinh(x);
}

constexpr auto cosh(auto x) noexcept {
  return std::cosh(x);
}

constexpr auto tanh(auto x) noexcept {
  return std::tanh(x);
}

constexpr auto asinh(auto x) noexcept {
  return std::asin(x);
}

constexpr auto acosh(auto x) noexcept {
  return std::acos(x);
}

constexpr auto atanh(auto x) noexcept {
  return std::atanh(x);
}

/// @}

/// @name Givens rotations.
/// @{

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

/// @}

} // namespace Storm::math
