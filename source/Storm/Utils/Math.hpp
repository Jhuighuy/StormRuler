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

/// @brief Check if the specified class is an instantion of a template.
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
concept bool_type = enable_bool_type_v<Bool>;

/// @brief Treat the specified type as the integer.
template<class Integer>
inline constexpr bool enable_integer_type_v = std::is_integral_v<Integer>;

/// @brief Integer type.
template<class Integer>
concept integer_type = enable_integer_type_v<Integer>;

/// @brief Treat the specified type as real (floating-point).
template<class Real>
inline constexpr bool enable_real_type_v = std::is_floating_point_v<Real>;

/// @brief Real (floating-point) type.
template<class Real>
concept real_type = enable_real_type_v<Real>;

/// @brief Treat the specified type as complex (floating-point).
template<class Complex>
inline constexpr bool enable_complex_type_v =
    is_instantiation_v<Complex, std::complex> &&
    enable_real_type_v<typename Complex::value_type>;

/// @brief Complex (floating-point) type.
template<class Complex>
concept complex_type = enable_complex_type_v<Complex>;

/// @brief Real or complex (floating-point) type.
template<class Type>
concept real_or_complex_type = real_type<Type> || complex_type<Type>;

/// @brief Numberic type.
template<class Type>
concept numeric_type = integer_type<Type> || real_or_complex_type<Type>;

// -----------------------------------------------------------------------------

static real_t deg2rad(real_t) noexcept;
static real_t rad2deg(real_t) noexcept;

// -----------------------------------------------------------------------------

/// @brief If @p y is zero, return zero,
/// else return value of @p x divided by @p y.
template<real_or_complex_type Value>
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

} // namespace Storm
