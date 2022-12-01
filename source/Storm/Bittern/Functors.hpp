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

namespace Storm
{

// -----------------------------------------------------------------------------

/// @brief Function call result type.
template<class Func, class... Args>
using func_result_t = decltype(std::declval<Func>()(std::declval<Args>()...));

/// @brief Bind the first arguments of the function to values.
template<std::copyable Func, std::copyable... FirstArgs>
  requires std::is_object_v<Func>
struct BindFirst {
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
    return std::apply(
        [this, &rest_args...](const FirstArgs&... first_args) {
          return func_(first_args..., std::forward<RestArgs>(rest_args)...);
        },
        first_args_);
  }

}; // struct BindFirst

/// @brief Bind the last arguments of the function to values.
template<std::copyable Func, std::copyable... LastArgs>
  requires std::is_object_v<Func>
struct BindLast {
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
    return std::apply(
        [this, &rest_args...](const LastArgs&... last_args) {
          return func_(std::forward<RestArgs>(rest_args)..., last_args...);
        },
        last_args_);
  }

}; // struct BindLast

/// @brief Compose two functions.
template<std::copyable Func1, std::copyable Func2>
  requires std::is_object_v<Func1> && std::is_object_v<Func2>
struct Compose {
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

}; // struct Compose

// -----------------------------------------------------------------------------

/// @brief Static cast function.
template<class To>
struct Cast {
  template<class Arg>
    requires std::convertible_to<std::remove_cvref_t<Arg>, To>
  constexpr To operator()(Arg&& arg) const noexcept
  {
    return static_cast<To>(std::forward<Arg>(arg));
  }
}; // struct Cast

// -----------------------------------------------------------------------------

using Not = std::logical_not<>;
using And = std::logical_and<>;
using Or = std::logical_or<>;

/// @brief Ternary operator function.
struct Merge {
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
}; // struct Merge

// -----------------------------------------------------------------------------

using Equal = std::equal_to<>;
using NotEqual = std::not_equal_to<>;
using Less = std::less<>;
using LessEqual = std::less_equal<>;
using Greater = std::greater<>;
using GreaterEqual = std::greater_equal<>;

/// @brief Approximately equals function.
struct ApproxEqual {
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

}; // struct ApproxEqual

/// @brief Minimum value function.
struct Min {
  template<class... Args>
  constexpr auto operator()(Args&&... args) const noexcept
  {
    using Result = std::common_type_t<std::remove_cvref_t<Args>...>;
    return min({static_cast<Result>(std::forward<Args>(args))...});
  }
}; // struct Min

/// @brief Maximum value function.
struct Max {
  template<class... Args>
  constexpr auto operator()(Args&&... args) const noexcept
  {
    using Result = std::common_type_t<std::remove_cvref_t<Args>...>;
    return max({static_cast<Result>(std::forward<Args>(args))...});
  }
}; // struct Max

// -----------------------------------------------------------------------------

using Negate = std::negate<>;
using Add = std::plus<>;
using Subtract = std::minus<>;
using Multiply = std::multiplies<>;
using Divide = std::divides<>;

struct AddAssign {
  template<class Arg1, class Arg2>
  constexpr decltype(auto) operator()(Arg1&& arg1, Arg2&& arg2) const noexcept
  {
    return std::forward<Arg1>(arg1) += std::forward<Arg2>(arg2);
  }
}; // struct AddAssign

struct SubtractAssign {
  template<class Arg1, class Arg2>
  constexpr decltype(auto) operator()(Arg1&& arg1, Arg2&& arg2) const noexcept
  {
    return std::forward<Arg1>(arg1) -= std::forward<Arg2>(arg2);
  }
}; // struct SubtractAssign

struct MultiplyAssign {
  template<class Arg1, class Arg2>
  constexpr decltype(auto) operator()(Arg1&& arg1, Arg2&& arg2) const noexcept
  {
    return std::forward<Arg1>(arg1) *= std::forward<Arg2>(arg2);
  }
}; // struct MultiplyAssign

struct DivideAssign {
  template<class Arg1, class Arg2>
  constexpr decltype(auto) operator()(Arg1&& arg1, Arg2&& arg2) const noexcept
  {
    return std::forward<Arg1>(arg1) /= std::forward<Arg2>(arg2);
  }
}; // struct DivideAssign

// -----------------------------------------------------------------------------

/// @brief Absolute value function.
struct Abs {
  template<class Arg>
  constexpr auto operator()(Arg&& arg) const noexcept
  {
    return abs(std::forward<Arg>(arg));
  }
}; // struct Abs

// -----------------------------------------------------------------------------

/// @brief Power function.
struct Pow {
  template<class Arg1, class Arg2>
  constexpr auto operator()(Arg1&& arg1, Arg2&& arg2) const noexcept
  {
    return pow(std::forward<Arg1>(arg1), std::forward<Arg2>(arg2));
  }
}; // struct Pow

// -----------------------------------------------------------------------------

/// @brief Scalar dot product function: @f$ x \cdot y := x \, \bar{y} @f$.
struct DotProduct {
  template<class Arg1, class Arg2>
  constexpr auto operator()(Arg1&& arg1, Arg2&& arg2) const noexcept
  {
    return std::forward<Arg1>(arg1) * conj(std::forward<Arg2>(arg2));
  }
}; // struct DotProduct

// -----------------------------------------------------------------------------

} // namespace Storm
