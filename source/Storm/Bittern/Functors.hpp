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

namespace Storm
{

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

/// @brief Bind the first argument of the function to a value.
template<class FirstArg, class Func>
  requires std::is_object_v<Func>
struct BindFirst {
private:

  STORM_NO_UNIQUE_ADDRESS_ FirstArg first_arg_;
  STORM_NO_UNIQUE_ADDRESS_ Func func_;

public:

  constexpr explicit BindFirst(FirstArg first_arg, Func func)
      : first_arg_{std::move(first_arg)}, func_{std::move(func)}
  {
  }

  template<class... RestArgs>
    requires std::regular_invocable<Func, FirstArg, RestArgs...>
  constexpr decltype(auto) operator()(RestArgs&&... rest_args) const noexcept
  {
    return func_(first_arg_, std::forward<RestArgs>(rest_args)...);
  }

}; // struct BindFirst

/// @brief Bind the last argument of the function to a value.
template<class LastArg, class Func>
  requires std::is_object_v<Func>
struct BindLast {
private:

  STORM_NO_UNIQUE_ADDRESS_ LastArg last_arg_;
  STORM_NO_UNIQUE_ADDRESS_ Func func_;

public:

  constexpr explicit BindLast(LastArg last_arg, Func func)
      : last_arg_{std::move(last_arg)}, func_{std::move(func)}
  {
  }

  template<class... RestArgs>
    requires std::regular_invocable<Func, RestArgs..., LastArg>
  constexpr decltype(auto) operator()(RestArgs&&... rest_args) const noexcept
  {
    return func_(std::forward<RestArgs>(rest_args)..., last_arg_);
  }

}; // struct BindLast

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

} // namespace Storm
