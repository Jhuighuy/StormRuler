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

#include <cassert>
#include <cmath>
#include <cstddef>

#include <complex>
#include <concepts>
#include <source_location>
#include <tuple>
#include <type_traits>

#include <StormRuler_API.h>

#ifdef _MSC_VER
#define StormAssume(x) __assume(x)
#else
#define StormAssume(x)                     \
  do {                                     \
    if (!(x)) { __builtin_unreachable(); } \
  } while (false)
#endif

/// @todo Reimplement me with std::source_location.
#if (!defined(__PRETTY_FUNCTION__) && !defined(__GNUC__))
#define __PRETTY_FUNCTION__ __FUNCSIG__
#endif

#define StormEnsure(x)                                                 \
  do {                                                                 \
    if (!(x)) {                                                        \
      std::fprintf(stderr, "\nAssertion failed:\n%s:%d %s: \"%s\".\n", \
                   __FILE__, __LINE__, __PRETTY_FUNCTION__, #x);       \
      std::fflush(stderr);                                             \
      std::abort();                                                    \
    }                                                                  \
  } while (false)

#ifdef NDEBUG
#define StormAssert(x) StormAssume(x)
#else
#define StormAssert(x) StormEnsure(x)
#endif

namespace Storm {

/// @brief Size type.
using size_t = std::size_t;

/// @brief Pointer difference type.
using ptrdiff_t = std::ptrdiff_t;

/// @brief Real floating-point type.
using real_t = double;

/// @brief Check if type is a complex floating point.
/// @{
template<class>
struct is_complex_floating_point_t : std::false_type {};
template<class T>
struct is_complex_floating_point_t<std::complex<T>> :
    std::bool_constant<std::is_floating_point_v<T>> {};
template<class T>
inline constexpr bool is_complex_floating_point_v =
    is_complex_floating_point_t<T>::value;
/// @}

/// @brief Concept for the real or complex floaing point numbers,
///   e.g. @c real_t or @c std::complex<real_t>.
template<class T>
concept real_or_complex_floating_point =
    std::floating_point<T> || is_complex_floating_point_v<T>;

namespace Utils {

/// @brief If @p y is zero, return zero,
///   else return value of @p x divided by @p y.
template<real_or_complex_floating_point Value>
auto SafeDivide(Value x, Value y) {
  static constexpr Value zero{0.0};
  return y == zero ? zero : (x / y);
}

/// @brief If @p y is zero, assign to @p x and return zero,
///   else assign to @p x and return value of @p x divided by @p y.
auto& SafeDivideEquals(real_or_complex_floating_point auto& x,
                       real_or_complex_floating_point auto y) {
  return x = SafeDivide(x, y);
}

/// @brief Generate the Givens rotation.
template<real_or_complex_floating_point Value>
auto SymOrtho(Value a, Value b) {
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

} // SymOrtho

} // namespace Utils

} // namespace Storm
