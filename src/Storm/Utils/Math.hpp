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

#include <Storm/Base.hpp>

#include <cmath>
#include <concepts>
#include <type_traits>

namespace Storm {

namespace detail_ {
  template<class>
  inline constexpr bool is_complex_floating_point_v_ = false;
  template<class T>
  inline constexpr bool is_complex_floating_point_v_<std::complex<T>> =
      std::is_floating_point_v<T>;
} // namespace detail_

/// @brief Real or complex floaing point numbers,
///   e.g. @c real_t or @c std::complex<real_t>.
template<class T>
concept real_or_complex_floating_point =
    std::floating_point<T> || detail_::is_complex_floating_point_v_<T>;

namespace math {

  using std::abs;

  /// @brief If @p y is zero, return zero,
  ///   else return value of @p x divided by @p y.
  template<real_or_complex_floating_point Value>
  auto safe_divide(Value x, Value y) {
    static constexpr Value zero{0.0};
    return y == zero ? zero : (x / y);
  }

  /// @name Exponential functions.
  /// @{

  using std::exp;

  using std::exp2;

  using std::log;

  using std::log2;

  using std::log10;

  /// @} // Exponential functions.

  /// @name Power functions.
  /// @{

  using std::pow;

  using std::sqrt;

  using std::cbrt;

  using std::hypot;

  /// @} // Power functions.

  /// @name Trigonometric functions.
  /// @{

  using std::sin;

  using std::cos;

  using std::tan;

  using std::asin;

  using std::acos;

  using std::atan;

  using std::atan2;

  /// @} // Trigonometric functions.

  /// @name Hyperbolic functions.
  /// @{

  using std::sinh;

  using std::cosh;

  using std::tanh;

  using std::asinh;

  using std::acosh;

  using std::atanh;

  /// @} // Hyperbolic functions.

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

  /// @} // Givens rotations.

} // namespace math
} // namespace Storm
