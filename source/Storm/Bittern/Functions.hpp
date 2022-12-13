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

#include <Storm/Bittern/Math.hpp>

#include <algorithm>
#include <concepts>
#include <functional>
#include <limits>
#include <tuple>
#include <utility>

namespace Storm
{

// -----------------------------------------------------------------------------

/// @brief Identity function.
using Identity = std::identity;

/// @brief Bind the first arguments of the function to values.
template<class Func, class... FirstArgs>
  requires std::is_object_v<Func>
class BindFirst final
{
private:

  STORM_NO_UNIQUE_ADDRESS_ Func func_;
  STORM_NO_UNIQUE_ADDRESS_ std::tuple<FirstArgs...> first_args_;

public:

  constexpr explicit BindFirst(Func func, FirstArgs... first_args)
      : func_{std::move(func)}, first_args_{std::move(first_args)...}
  {
  }

  template<class... RestArgs>
    requires std::regular_invocable<Func, FirstArgs..., RestArgs...>
  constexpr decltype(auto) operator()(RestArgs&&... rest_args) const noexcept
  {
    auto call = [this, &rest_args...](const FirstArgs&... first_args) {
      return func_(first_args..., std::forward<RestArgs>(rest_args)...);
    };
    return std::apply(call, first_args_);
  }

}; // class BindFirst

/// @brief Bind the last arguments of the function to values.
template<class Func, class... LastArgs>
  requires std::is_object_v<Func>
class BindLast final
{
private:

  STORM_NO_UNIQUE_ADDRESS_ Func func_;
  STORM_NO_UNIQUE_ADDRESS_ std::tuple<LastArgs...> last_args_;

public:

  constexpr explicit BindLast(Func func, LastArgs... last_args)
      : func_{std::move(func)}, last_args_{std::move(last_args)...}
  {
  }

  template<class... RestArgs>
    requires std::regular_invocable<Func, RestArgs..., LastArgs...>
  constexpr decltype(auto) operator()(RestArgs&&... rest_args) const noexcept
  {
    auto call = [this, &rest_args...](const LastArgs&... last_args) {
      return func_(std::forward<RestArgs>(rest_args)..., last_args...);
    };
    return std::apply(call, last_args_);
  }

}; // class BindLast

/// @brief Compose multiple functions: @f$ f_{1} \circ f_{2} \circ \ldots @f$.
template<class Func, class... RestFuncs>
  requires std::is_object_v<Func> && (... && std::is_object_v<RestFuncs>)
class Compose final
{
private:

  STORM_NO_UNIQUE_ADDRESS_ Func func_;
  STORM_NO_UNIQUE_ADDRESS_ Compose<RestFuncs...> rest_funcs_;

public:

  constexpr Compose(Func func, RestFuncs... rest_funcs)
      : func_{std::move(func)}, rest_funcs_{std::move(rest_funcs)...}
  {
  }

  template<class... Args>
    requires std::regular_invocable<Compose<RestFuncs...>, Args...> &&
             std::regular_invocable<
                 Func, std::invoke_result_t<Compose<RestFuncs...>, Args...>>
  constexpr decltype(auto) operator()(Args&&... args) const noexcept
  {
    return func_(rest_funcs_(std::forward<Args>(args)...));
  }

}; // class Compose

/// @brief Compose two functions: @f$ f_{1} \circ f_{2} @f$.
template<class Func1, class Func2>
  requires std::is_object_v<Func1> && std::is_object_v<Func2>
class Compose<Func1, Func2> final
{
private:

  STORM_NO_UNIQUE_ADDRESS_ Func1 func1_;
  STORM_NO_UNIQUE_ADDRESS_ Func2 func2_;

public:

  constexpr Compose(Func1 func1, Func2 func2)
      : func1_{std::move(func1)}, func2_{std::move(func2)}
  {
  }

  template<class... Args>
    requires std::regular_invocable<Func2, Args...> && //
             std::regular_invocable<Func1, std::invoke_result_t<Func2, Args...>>
  constexpr decltype(auto) operator()(Args&&... args) const noexcept
  {
    return func1_(func2_(std::forward<Args>(args)...));
  }

}; // class Compose

// -----------------------------------------------------------------------------

/// @brief Static cast function.
template<class To>
class Cast final
{
public:

  template<class Arg>
    requires std::convertible_to<std::remove_cvref_t<Arg>, To>
  constexpr To operator()(Arg&& arg) const noexcept
  {
    return static_cast<To>(std::forward<Arg>(arg));
  }

}; // class Cast

// -----------------------------------------------------------------------------

using Not = std::logical_not<>;

/// @brief Logical AND function.
class And final
{
public:

  template<class... Args>
  constexpr bool operator()(Args&&... args) const noexcept
  {
    return (std::forward<Args>(args) && ...);
  }

}; // class And

/// @brief Logical OR function.
class Or final
{
public:

  template<class... Args>
  constexpr bool operator()(Args&&... args) const noexcept
  {
    return (std::forward<Args>(args) || ...);
  }

}; // class Or

/// @brief Ternary operator function.
class Merge final
{
public:

  template<class Cond, class Then, class Else>
    requires std::convertible_to<Cond, bool> && std::common_with<Then, Else>
  constexpr auto operator()(Cond&& cond, //
                            Then&& then_arg, Else&& else_arg) const noexcept
  {
    using Result = std::common_type_t<Then, Else>;
    return static_cast<bool>(std::forward<Cond>(cond)) ?
               static_cast<Result>(std::forward<Then>(then_arg)) :
               static_cast<Result>(std::forward<Else>(else_arg));
  }

}; // class Merge

// -----------------------------------------------------------------------------

using Equal = std::equal_to<>;
using NotEqual = std::not_equal_to<>;
using Less = std::less<>;
using LessEqual = std::less_equal<>;
using Greater = std::greater<>;
using GreaterEqual = std::greater_equal<>;

/// @brief Approximately equals function.
class ApproxEqual final
{
private:

  long double tolerance_;

public:

  constexpr explicit ApproxEqual(
      long double tolerance = sqrt(std::numeric_limits<long double>::epsilon()))
      : tolerance_{tolerance}
  {
    STORM_ASSERT_(tolerance > 0.0l, "Negative tolerance!");
  }

  template<class Arg1, class Arg2>
  constexpr bool operator()(Arg1&& arg1, Arg2&& arg2) const noexcept
  {
    decltype(auto) delta =
        abs(std::forward<Arg1>(arg1) - std::forward<Arg2>(arg2));
    return static_cast<long double>(delta) <= tolerance_;
  }

}; // class ApproxEqual

/// @brief Minimum value function.
class Min final
{
public:

  template<class... Args>
  constexpr auto operator()(Args&&... args) const noexcept
  {
    using Result = std::common_type_t<std::remove_cvref_t<Args>...>;
    return min({static_cast<Result>(std::forward<Args>(args))...});
  }

}; // class Min

/// @brief Maximum value function.
class Max final
{
public:

  template<class... Args>
  constexpr auto operator()(Args&&... args) const noexcept
  {
    using Result = std::common_type_t<std::remove_cvref_t<Args>...>;
    return max({static_cast<Result>(std::forward<Args>(args))...});
  }

}; // class Max

// -----------------------------------------------------------------------------

/// @brief Negate function.
using Negate = std::negate<>;

/// @brief Add function.
class Add final
{
public:

  template<class... Args>
  constexpr auto operator()(Args&&... args) const noexcept
  {
    return (std::forward<Args>(args) + ...);
  }

}; // class Add

/// @brief Subtract function.
using Subtract = std::minus<>;

/// @brief Multiply function.
class Multiply final
{
public:

  template<class... Args>
  constexpr auto operator()(Args&&... args) const noexcept
  {
    return (std::forward<Args>(args) * ...);
  }

}; // class Multiply

/// @brief Divide function.
using Divide = std::divides<>;

/// @brief Assign action function.
class Assign final
{
public:

  template<class Arg1, class Arg2>
  constexpr decltype(auto) operator()(Arg1&& arg1, Arg2&& arg2) const noexcept
  {
    return std::forward<Arg1>(arg1) = std::forward<Arg2>(arg2);
  }

}; // class Assign

/// @brief Add-assign action function.
class AddAssign final
{
public:

  template<class Arg1, class Arg2>
  constexpr decltype(auto) operator()(Arg1&& arg1, Arg2&& arg2) const noexcept
  {
    return std::forward<Arg1>(arg1) += std::forward<Arg2>(arg2);
  }

}; // class AddAssign

/// @brief Subtract-assign action function.
class SubtractAssign final
{
public:

  template<class Arg1, class Arg2>
  constexpr decltype(auto) operator()(Arg1&& arg1, Arg2&& arg2) const noexcept
  {
    return std::forward<Arg1>(arg1) -= std::forward<Arg2>(arg2);
  }

}; // class SubtractAssign

/// @brief Multiply-assign action function.
class MultiplyAssign final
{
public:

  template<class Arg1, class Arg2>
  constexpr decltype(auto) operator()(Arg1&& arg1, Arg2&& arg2) const noexcept
  {
    return std::forward<Arg1>(arg1) *= std::forward<Arg2>(arg2);
  }

}; // class MultiplyAssign

/// @brief Divide-assign action function.
class DivideAssign final
{
public:

  template<class Arg1, class Arg2>
  constexpr decltype(auto) operator()(Arg1&& arg1, Arg2&& arg2) const noexcept
  {
    return std::forward<Arg1>(arg1) /= std::forward<Arg2>(arg2);
  }

}; // class DivideAssign

// -----------------------------------------------------------------------------

/// @brief Real component function.
class Real final
{
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const noexcept
  {
    return real(std::forward<Arg>(arg));
  }

}; // class Real

/// @brief Imaginary component function.
class Imag final
{
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const noexcept
  {
    return imag(std::forward<Arg>(arg));
  }

}; // class Imag

/// @brief Conjugate function.
class Conj final
{
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const noexcept
  {
    return conj(std::forward<Arg>(arg));
  }

}; // class Conj

// -----------------------------------------------------------------------------

/// @brief Absolute value function.
class Abs final
{
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const noexcept
  {
    return abs(std::forward<Arg>(arg));
  }

}; // class Abs

/// @brief Squared absolute value function: @f$ |x|^{2} @f$.
class AbsSquared final
{
public:

  template<class Arg>
  constexpr auto operator()(const Arg& arg) const noexcept
  {
    return real(arg * conj(arg));
  }

}; // class AbsSquared

/// @brief Sign function.
class Sign final
{
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const noexcept
  {
    return sign(std::forward<Arg>(arg));
  }

}; // class Sign

// -----------------------------------------------------------------------------

/// @brief Square root function.
class Sqrt final
{
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const noexcept
  {
    return sqrt(std::forward<Arg>(arg));
  }

}; // class Sqrt

/// @brief Cube root function.
class Cbrt final
{
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const noexcept
  {
    return cbrt(std::forward<Arg>(arg));
  }

}; // class Cbrt

/// @brief Power function.
class Pow final
{
public:

  template<class Arg1, class Arg2>
  constexpr auto operator()(Arg1&& arg1, Arg2&& arg2) const noexcept
  {
    return pow(std::forward<Arg1>(arg1), std::forward<Arg2>(arg2));
  }

}; // class Pow

// -----------------------------------------------------------------------------

/// @brief Exponent function.
class Exp final
{
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const noexcept
  {
    return exp(std::forward<Arg>(arg));
  }

}; // class Exp

/// @brief Exponent (base 2) function.
class Exp2 final
{
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const noexcept
  {
    return exp2(std::forward<Arg>(arg));
  }

}; // class Exp2

/// @brief Logarithm function.
class Log final
{
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const noexcept
  {
    return log(std::forward<Arg>(arg));
  }

}; // class Log

/// @brief Logarithm (base 2) function.
class Log2 final
{
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const noexcept
  {
    return log2(std::forward<Arg>(arg));
  }

}; // class Log2

/// @brief Logarithm (base 10) function.
class Log10 final
{
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const noexcept
  {
    return log10(std::forward<Arg>(arg));
  }

}; // class Log10

// -----------------------------------------------------------------------------

/// @brief Sine function.
class Sin final
{
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const noexcept
  {
    return sin(std::forward<Arg>(arg));
  }

}; // class Sin

/// @brief Cosine function.
class Cos final
{
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const noexcept
  {
    return cos(std::forward<Arg>(arg));
  }

}; // class Cos

/// @brief Tangent function.
class Tan final
{
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const noexcept
  {
    return tan(std::forward<Arg>(arg));
  }

}; // class Tan

/// @brief Inverse sine function.
class Asin final
{
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const noexcept
  {
    return asin(std::forward<Arg>(arg));
  }

}; // class Asin

/// @brief Inverse cosine function.
class Acos final
{
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const noexcept
  {
    return acos(std::forward<Arg>(arg));
  }

}; // class Acos

/// @brief Inverse tangent function.
class Atan final
{
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const noexcept
  {
    return atan(std::forward<Arg>(arg));
  }

}; // class Atan

// -----------------------------------------------------------------------------

/// @brief Hyperbolic sine function.
class Sinh final
{
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const noexcept
  {
    return sinh(std::forward<Arg>(arg));
  }

}; // class Sinh

/// @brief Hyperbolic cosine function.
class Cosh final
{
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const noexcept
  {
    return cosh(std::forward<Arg>(arg));
  }

}; // class Cosh

/// @brief Hyperbolic tangent function.
class Tanh final
{
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const noexcept
  {
    return tanh(std::forward<Arg>(arg));
  }

}; // class Tanh

/// @brief Inverse hyperbolic sine function.
class Asinh final
{
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const noexcept
  {
    return asinh(std::forward<Arg>(arg));
  }

}; // class Asinh

/// @brief Inverse hyperbolic cosine function.
class Acosh final
{
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const noexcept
  {
    return acosh(std::forward<Arg>(arg));
  }

}; // class Acosh

/// @brief Inverse hyperbolic tangent function.
class Atanh final
{
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const noexcept
  {
    return atanh(std::forward<Arg>(arg));
  }

}; // class Atanh

// -----------------------------------------------------------------------------

/// @brief Scalar dot product function: @f$ x \cdot y := x \, \bar{y} @f$.
class DotProduct final
{
public:

  template<class Arg1, class Arg2>
  constexpr auto operator()(Arg1&& arg1, Arg2&& arg2) const noexcept
  {
    return std::forward<Arg1>(arg1) * conj(std::forward<Arg2>(arg2));
  }

}; // class DotProduct

// -----------------------------------------------------------------------------

} // namespace Storm
