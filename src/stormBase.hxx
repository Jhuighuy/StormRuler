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
#ifndef _STORM_BASE_HXX_
#define _STORM_BASE_HXX_

#include <cstddef>

#include <array>
#include <complex>
#include <concepts>
#include <iostream>
#include <type_traits>
#include <vector>

#include <StormRuler_API.h>

#define StormEnabledAssert(x) assert(x)
#define StormDisabledAssert(x) static_cast<void>(x)
#define StormAssert stormAssert

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

/// @brief If @p y is zero, return zero, else
///   return value of @p x divided by @p y.
auto SafeDivide(real_or_complex_floating_point auto x,
                real_or_complex_floating_point auto y) {
  return (y == 0.0) ? 0.0 : (x / y);
}

/// @brief If @p y is zero, assign to @p x and return zero, else
///   assign to @p x and return value of @p x divided by @p y.
auto& SafeDivideEquals(real_or_complex_floating_point auto& x,
                       real_or_complex_floating_point auto y) {
  return x = SafeDivide(x, y);
}

} // namespace Utils

static size_t const DynamicExtent = SIZE_MAX;

template<size_t Extent>
using ExtentSize = std::conditional_t<Extent == DynamicExtent, size_t,
                                      std::integral_constant<size_t, Extent>>;

template<class Value, size_t Extent>
using ExtentArray =
    std::conditional_t<Extent == DynamicExtent, std::vector<Value>,
                       std::array<Value, Extent>>;

/// ----------------------------------------------------------------- ///
/// @brief Constant to mark the the dynamic size in the \
///    SOD (static or dynamic) containers.
/// ----------------------------------------------------------------- ///
static constexpr size_t DynamicSize{SIZE_MAX};

/// ----------------------------------------------------------------- ///
/// @brief SOD @c size_t template.
///
/// Usage:
/// @code
///   [[no_unique_address]] SodSize_t<Size> MySize;
/// @endcode
/// ----------------------------------------------------------------- ///
template<size_t Size>
using SodSize_t = std::conditional_t<Size == DynamicSize, size_t,
                                     std::integral_constant<size_t, Size>>;

/// ----------------------------------------------------------------- ///
/// @brief SOD array template.
///
/// Usage:
/// @code
///   [[no_unique_address]] SodArray<Value, Size> MyArray;
/// @endcode
/// ----------------------------------------------------------------- ///
template<class Value, size_t Size>
using SodArray = std::conditional_t<Size == DynamicSize, std::vector<Value>,
                                    std::array<Value, Size>>;

} // namespace Storm

#endif // ifndef _STORM_BASE_HXX_
