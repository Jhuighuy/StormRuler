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

#include <cmath>
#include <complex>
#include <concepts>
#include <numbers>
#include <type_traits>

namespace Storm {

// -----------------------------------------------------------------------------

/// @brief Check if the specified class is an instantion of a template.
/// @todo Move me somewhere else!
/// @{
template<class, template<class...> class>
inline constexpr bool is_instantiation_v = false;
template<class... Args, template<class...> class Template>
inline constexpr bool is_instantiation_v<Template<Args...>, Template> = true;
/// @}

// -----------------------------------------------------------------------------

/// @brief Treat the specified type as boolean.
template<class Bool>
inline constexpr bool enable_bool_type_v = std::is_same_v<Bool, bool>;

/// @brief Bool type.
template<class Bool>
concept bool_type = enable_bool_type_v<std::remove_cvref_t<Bool>>;

/// @brief Treat the specified type as the integer.
template<class Integer>
inline constexpr bool enable_integer_type_v = std::is_integral_v<Integer>;

/// @brief Integer type.
template<class Integer>
concept integer_type = enable_integer_type_v<std::remove_cvref_t<Integer>>;

/// @brief Treat the specified type as real (floating-point).
template<class Real>
inline constexpr bool enable_real_type_v = std::is_floating_point_v<Real>;

/// @brief Real (floating-point) type.
template<class Real>
concept real_type = enable_real_type_v<std::remove_cvref_t<Real>>;

/// @brief Treat the specified type as complex (floating-point).
template<class Complex>
inline constexpr bool enable_complex_type_v =
    is_instantiation_v<Complex, std::complex> /*&&
    enable_real_type_v<typename Complex::value_type>*/
    ;

/// @brief Complex (floating-point) type.
template<class Complex>
concept complex_type = enable_complex_type_v<std::remove_cvref_t<Complex>>;

/// @brief Real or complex (floating-point) type.
template<class Type>
concept real_or_complex_type = real_type<Type> || complex_type<Type>;

/// @brief Numberic type.
template<class Type>
concept numeric_type = integer_type<Type> || real_or_complex_type<Type>;

// -----------------------------------------------------------------------------

/// @brief If @p y is zero, return zero,
/// else return value of @p x divided by @p y.
template<real_or_complex_type Value>
constexpr auto safe_divide(Value x, Value y) {
  constexpr Value zero{0.0};
  return y == zero ? zero : (x / y);
}

// -----------------------------------------------------------------------------

/// @brief @f$ \pi @f$ constant.
inline constexpr real_t pi = std::numbers::pi_v<real_t>;

/// @brief Convert degrees to radians.
template<real_type Real>
constexpr auto deg2rad(Real&& degrees) noexcept {
  return (pi / 180.0) * std::forward<Real>(degrees);
}

/// @brief Convert radians to degrees.
template<real_type Real>
constexpr auto rad2deg(Real&& radians) noexcept {
  return (180.0 / pi) * std::forward<Real>(radians);
}

// -----------------------------------------------------------------------------

/// @brief @f$ i = \sqrt{-1} @f$.
inline constexpr complex_t i{0.0_dp, 1.0_dp};

using std::real;

using std::imag;

/// @brief Conjugate of real number argument @p arg (noop).
/// @note Unlike standard function, this function preserves real type.
template<real_type Arg>
constexpr auto conj(Arg&& arg) noexcept {
  return std::forward<Arg>(arg);
}
/// @brief Conjugate of complex number argument @p arg.
template<real_type Type>
constexpr auto conj(const std::complex<Type>& arg) noexcept {
  return std::conj(arg);
}

// -----------------------------------------------------------------------------

template<class Arg>
  requires real_type<Arg> || integer_type<Arg>
constexpr int sign(const Arg& arg) noexcept {
  return (Arg{0} < arg) - (arg < Arg{0});
}

using std::abs;

using std::min;

using std::max;

// -----------------------------------------------------------------------------

using std::pow;

/// @brief @f$ \sqrt{2} @f$ constant.
inline constexpr real_t sqrt2 = std::numbers::sqrt2_v<real_t>;

/// @brief @f$ \sqrt{3} @f$ constant.
inline constexpr real_t sqrt3 = std::numbers::sqrt3_v<real_t>;

using std::sqrt;

using std::cbrt;

using std::hypot;

// -----------------------------------------------------------------------------

/// @brief @f$ e @f$ constant.
inline constexpr real_t e = std::numbers::e_v<real_t>;

using std::exp;

using std::exp2;

using std::log;

using std::log2;

using std::log10;

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
/// @todo Move me somewhere else!
template<real_or_complex_type Value>
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

// -----------------------------------------------------------------------------

} // namespace Storm
