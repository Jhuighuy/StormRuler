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

#include <Storm/Crow/MathUtils.hpp>

#include <algorithm>
#include <concepts>
#include <functional>
#include <limits>
#include <tuple>
#include <utility>

namespace Storm {

// -----------------------------------------------------------------------------

/// @brief Identity function.
using Identity = std::identity;

/// @brief Bind the first arguments of the function to values.
template<class Func, class... FirstArgs>
  requires std::is_object_v<Func>
class BindFirst final {
private:

  STORM_NO_UNIQUE_ADDRESS Func _func;
  STORM_NO_UNIQUE_ADDRESS std::tuple<FirstArgs...> _first_args;

public:

  constexpr explicit BindFirst(Func func, FirstArgs... first_args)
      : _func{std::move(func)}, _first_args{std::move(first_args)...} {}

  template<class... RestArgs>
    requires std::regular_invocable<Func, FirstArgs..., RestArgs...>
  constexpr decltype(auto) operator()(RestArgs&&... rest_args) const {
    const auto call = [&](const FirstArgs&... first_args) {
      return _func(first_args..., std::forward<RestArgs>(rest_args)...);
    };
    return std::apply(call, _first_args);
  }

}; // class BindFirst

/// @brief Bind the last arguments of the function to values.
template<class Func, class... LastArgs>
  requires std::is_object_v<Func>
class BindLast final {
private:

  STORM_NO_UNIQUE_ADDRESS Func _func;
  STORM_NO_UNIQUE_ADDRESS std::tuple<LastArgs...> _last_args;

public:

  constexpr explicit BindLast(Func func, LastArgs... last_args)
      : _func{std::move(func)}, _last_args{std::move(last_args)...} {}

  template<class... RestArgs>
    requires std::regular_invocable<Func, RestArgs..., LastArgs...>
  constexpr decltype(auto) operator()(RestArgs&&... rest_args) const {
    const auto call = [&](const LastArgs&... last_args) {
      return _func(std::forward<RestArgs>(rest_args)..., last_args...);
    };
    return std::apply(call, _last_args);
  }

}; // class BindLast

/// @brief Compose multiple functions: @f$ f_{1} \circ f_{2} \circ \ldots @f$.
template<class Func, class... RestFuncs>
  requires std::is_object_v<Func> && (... && std::is_object_v<RestFuncs>)
class Compose final {
private:

  STORM_NO_UNIQUE_ADDRESS Func _func;
  STORM_NO_UNIQUE_ADDRESS Compose<RestFuncs...> _rest_funcs;

public:

  constexpr explicit Compose(Func func, RestFuncs... rest_funcs)
      : _func{std::move(func)}, _rest_funcs{std::move(rest_funcs)...} {}

  template<class... Args>
    requires std::invocable<Compose<RestFuncs...>, Args...> &&
             std::invocable<
                 Func, std::invoke_result_t<Compose<RestFuncs...>, Args...>>
  constexpr decltype(auto) operator()(Args&&... args) const {
    return _func(_rest_funcs(std::forward<Args>(args)...));
  }

}; // class Compose

/// @brief Compose two functions: @f$ f_{1} \circ f_{2} @f$.
template<class Func1, class Func2>
  requires std::is_object_v<Func1> && std::is_object_v<Func2>
class Compose<Func1, Func2> final {
private:

  STORM_NO_UNIQUE_ADDRESS Func1 func1_;
  STORM_NO_UNIQUE_ADDRESS Func2 func2_;

public:

  constexpr Compose(Func1 func1, Func2 func2)
      : func1_{std::move(func1)}, func2_{std::move(func2)} {}

  template<class... Args>
    requires std::invocable<Func2, Args...> && //
             std::invocable<Func1, std::invoke_result_t<Func2, Args...>>
  constexpr decltype(auto) operator()(Args&&... args) const {
    return func1_(func2_(std::forward<Args>(args)...));
  }

}; // class Compose

/// @brief Constant function.
template<class Value>
  requires std::is_object_v<Value>
class Constant final {
private:

  STORM_NO_UNIQUE_ADDRESS Value _value;

public:

  constexpr explicit Constant(Value value) : _value{std::move(value)} {}

  constexpr const Value& operator()(auto&&...) const {
    return _value;
  }

}; // class Constant

/// @brief Eye function.
template<class Value>
  requires std::is_object_v<Value>
class Eye final {
private:

  STORM_NO_UNIQUE_ADDRESS Value _diagonal;
  STORM_NO_UNIQUE_ADDRESS Value _off_diagonal;

public:

  constexpr Eye(Value diagonal, Value off_diagonal)
      : _diagonal{std::move(diagonal)}, _off_diagonal{std::move(off_diagonal)} {
  }

  template<class Arg, class... RestArgs>
  constexpr const Value& operator()(Arg&& arg, RestArgs&&... rest_args) const {
    const bool on_diagonal = [&]() {
      if constexpr (sizeof...(RestArgs) == 0) {
        using Index = std::remove_cvref_t<Arg>;
        static_assert(std::constructible_from<Index, size_t>);
        return std::forward<Arg>(arg) == Index{0};
      } else {
        return ((arg == std::forward<RestArgs>(rest_args)) && ...);
      }
    }();
    return on_diagonal ? _diagonal : _off_diagonal;
  }

}; // class Eye

// -----------------------------------------------------------------------------

/// @brief Move-assign action function.
class MoveAssign final {
public:

  template<class Arg1, class Arg2>
    requires std::assignable_from<Arg1, Arg2>
  constexpr decltype(auto) operator()(Arg1&& arg1, Arg2&& arg2) const {
    return std::forward<Arg1>(arg1) = std::move(std::forward<Arg2>(arg2));
  }

}; // class MoveAssign

/// @brief Copy-assign action function.
class CopyAssign final {
public:

  template<class Arg1, class Arg2>
    requires std::assignable_from<Arg1, Arg2>
  constexpr decltype(auto) operator()(Arg1&& arg1, Arg2&& arg2) const {
    return std::forward<Arg1>(arg1) = std::forward<Arg2>(arg2);
  }

}; // class CopyAssign

// -----------------------------------------------------------------------------

/// @brief Static cast function.
template<class To>
class Cast final {
public:

  template<class Arg>
    requires std::convertible_to<std::remove_cvref_t<Arg>, To>
  constexpr To operator()(Arg&& arg) const {
    return static_cast<To>(std::forward<Arg>(arg));
  }

}; // class Cast

// -----------------------------------------------------------------------------

/// @brief Logical NOT function.
using LogicalNot = std::logical_not<>;

/// @brief Logical AND function.
class LogicalAnd final {
public:

  template<class... Args>
  constexpr bool operator()(Args&&... args) const {
    return (std::forward<Args>(args) && ...);
  }

}; // class LogicalAnd

/// @brief Logical OR function.
class LogicalOr final {
public:

  template<class... Args>
  constexpr bool operator()(Args&&... args) const {
    return (std::forward<Args>(args) || ...);
  }

}; // class LogicalOr

/// @brief Ternary operator function.
class Merge final {
public:

  template<class Cond, class Then, class Else>
    requires std::convertible_to<Cond, bool> && std::common_with<Then, Else>
  constexpr auto operator()(Cond&& cond, //
                            Then&& then_arg, Else&& else_arg) const {
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
class ApproxEqual final {
private:

  long double _tolerance;

public:

  constexpr explicit ApproxEqual(
      long double tolerance = sqrt(std::numeric_limits<long double>::epsilon()))
      : _tolerance{tolerance} {
    STORM_ASSERT(tolerance > 0.0l, "Negative tolerance!");
  }

  template<class Arg1, class Arg2>
  constexpr bool operator()(Arg1&& arg1, Arg2&& arg2) const {
    const auto delta = static_cast<long double>(
        abs(std::forward<Arg1>(arg1) - std::forward<Arg2>(arg2)));
    return delta <= _tolerance;
  }

}; // class ApproxEqual

/// @brief Minimum value function.
class Min final {
public:

  template<class... Args>
  constexpr auto operator()(Args&&... args) const {
    return min(std::forward<Args>(args)...);
  }

}; // class Min

/// @brief Maximum value function.
class Max final {
public:

  template<class... Args>
  constexpr auto operator()(Args&&... args) const {
    return max(std::forward<Args>(args)...);
  }

}; // class Max

// -----------------------------------------------------------------------------

/// @brief Negate function.
using Negate = std::negate<>;

/// @brief Add function.
class Add final {
public:

  template<class... Args>
  constexpr auto operator()(Args&&... args) const {
    return (std::forward<Args>(args) + ...);
  }

}; // class Add

/// @brief Subtract function.
using Subtract = std::minus<>;

/// @brief Multiply function.
class Multiply final {
public:

  template<class... Args>
  constexpr auto operator()(Args&&... args) const {
    return (std::forward<Args>(args) * ...);
  }

}; // class Multiply

/// @brief Divide function.
using Divide = std::divides<>;

/// @brief Assign action function.
class Assign final {
public:

  template<class Arg1, class Arg2>
  constexpr decltype(auto) operator()(Arg1&& arg1, Arg2&& arg2) const noexcept {
    return std::forward<Arg1>(arg1) = std::forward<Arg2>(arg2);
  }

}; // class Assign

/// @brief Add-assign action function.
class AddAssign final {
public:

  template<class Arg1, class Arg2>
  constexpr decltype(auto) operator()(Arg1&& arg1, Arg2&& arg2) const {
    return std::forward<Arg1>(arg1) += std::forward<Arg2>(arg2);
  }

}; // class AddAssign

/// @brief Subtract-assign action function.
class SubtractAssign final {
public:

  template<class Arg1, class Arg2>
  constexpr decltype(auto) operator()(Arg1&& arg1, Arg2&& arg2) const {
    return std::forward<Arg1>(arg1) -= std::forward<Arg2>(arg2);
  }

}; // class SubtractAssign

/// @brief Multiply-assign action function.
class MultiplyAssign final {
public:

  template<class Arg1, class Arg2>
  constexpr decltype(auto) operator()(Arg1&& arg1, Arg2&& arg2) const {
    return std::forward<Arg1>(arg1) *= std::forward<Arg2>(arg2);
  }

}; // class MultiplyAssign

/// @brief Divide-assign action function.
class DivideAssign final {
public:

  template<class Arg1, class Arg2>
  constexpr decltype(auto) operator()(Arg1&& arg1, Arg2&& arg2) const {
    return std::forward<Arg1>(arg1) /= std::forward<Arg2>(arg2);
  }

}; // class DivideAssign

// -----------------------------------------------------------------------------

/// @brief Real component function.
class Real final {
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const {
    return real(std::forward<Arg>(arg));
  }

}; // class Real

/// @brief Imaginary component function.
class Imag final {
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const {
    return imag(std::forward<Arg>(arg));
  }

}; // class Imag

/// @brief Conjugate function.
class Conj final {
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const {
    return conj(std::forward<Arg>(arg));
  }

}; // class Conj

/// @brief Scalar dot product function: @f$ x \cdot y := x \, \bar{y} @f$.
class DotProduct final {
public:

  template<class Arg1, class Arg2>
  constexpr auto operator()(Arg1&& arg1, Arg2&& arg2) const {
    return dot_product(std::forward<Arg1>(arg1), std::forward<Arg2>(arg2));
  }

}; // class DotProduct

// -----------------------------------------------------------------------------

/// @brief Absolute value function.
class Abs final {
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const {
    return abs(std::forward<Arg>(arg));
  }

}; // class Abs

/// @brief Squared absolute value function: @f$ |x|^{2} @f$.
class AbsSquared final {
public:

  template<class Arg>
  constexpr auto operator()(const Arg& arg) const {
    return real(arg * conj(arg));
  }

}; // class AbsSquared

/// @brief Sign function.
class Sign final {
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const {
    return sign(std::forward<Arg>(arg));
  }

}; // class Sign

// -----------------------------------------------------------------------------

/// @brief Square root function.
class Sqrt final {
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const {
    return sqrt(std::forward<Arg>(arg));
  }

}; // class Sqrt

/// @brief Cube root function.
class Cbrt final {
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const {
    return cbrt(std::forward<Arg>(arg));
  }

}; // class Cbrt

/// @brief Power function.
class Pow final {
public:

  template<class Base, class Exponent>
  constexpr auto operator()(Base&& base, Exponent&& exp) const {
    return pow(std::forward<Base>(base), std::forward<Exponent>(exp));
  }

}; // class Pow

// -----------------------------------------------------------------------------

/// @brief Exponent function.
class Exp final {
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const {
    return exp(std::forward<Arg>(arg));
  }

}; // class Exp

/// @brief Exponent (base 2) function.
class Exp2 final {
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const {
    return exp2(std::forward<Arg>(arg));
  }

}; // class Exp2

/// @brief Logarithm function.
class Log final {
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const {
    return log(std::forward<Arg>(arg));
  }

}; // class Log

/// @brief Logarithm (base 2) function.
class Log2 final {
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const {
    return log2(std::forward<Arg>(arg));
  }

}; // class Log2

/// @brief Logarithm (base 10) function.
class Log10 final {
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const {
    return log10(std::forward<Arg>(arg));
  }

}; // class Log10

// -----------------------------------------------------------------------------

/// @brief Sine function.
class Sin final {
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const {
    return sin(std::forward<Arg>(arg));
  }

}; // class Sin

/// @brief Cosine function.
class Cos final {
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const {
    return cos(std::forward<Arg>(arg));
  }

}; // class Cos

/// @brief Tangent function.
class Tan final {
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const {
    return tan(std::forward<Arg>(arg));
  }

}; // class Tan

/// @brief Inverse sine function.
class Asin final {
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const {
    return asin(std::forward<Arg>(arg));
  }

}; // class Asin

/// @brief Inverse cosine function.
class Acos final {
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const {
    return acos(std::forward<Arg>(arg));
  }

}; // class Acos

/// @brief Inverse tangent function.
class Atan final {
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const {
    return atan(std::forward<Arg>(arg));
  }

}; // class Atan

// -----------------------------------------------------------------------------

/// @brief Hyperbolic sine function.
class Sinh final {
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const {
    return sinh(std::forward<Arg>(arg));
  }

}; // class Sinh

/// @brief Hyperbolic cosine function.
class Cosh final {
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const {
    return cosh(std::forward<Arg>(arg));
  }

}; // class Cosh

/// @brief Hyperbolic tangent function.
class Tanh final {
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const {
    return tanh(std::forward<Arg>(arg));
  }

}; // class Tanh

/// @brief Inverse hyperbolic sine function.
class Asinh final {
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const {
    return asinh(std::forward<Arg>(arg));
  }

}; // class Asinh

/// @brief Inverse hyperbolic cosine function.
class Acosh final {
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const {
    return acosh(std::forward<Arg>(arg));
  }

}; // class Acosh

/// @brief Inverse hyperbolic tangent function.
class Atanh final {
public:

  template<class Arg>
  constexpr auto operator()(Arg&& arg) const {
    return atanh(std::forward<Arg>(arg));
  }

}; // class Atanh

// -----------------------------------------------------------------------------

} // namespace Storm
