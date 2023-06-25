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

#include <algorithm>
#include <concepts>
#include <tuple>
#include <type_traits>
#include <utility>

namespace Storm {

// -----------------------------------------------------------------------------

// Apply `fixed_size_t<0>{}`, ..., `fixed_size_t<N-1>{}` to function.
template<size_t N, class Func>
constexpr decltype(auto) base_apply(Func&& func) {
  return ([&]<size_t... Is>(std::index_sequence<Is...>) -> decltype(auto) {
    return func(fixed_size_t<Is>{}...);
  }(std::make_index_sequence<N>{}));
}

// Checks if `std::tuple_size` is applicable.
template<class Tuple>
concept has_tuple_size =
    requires { typename std::tuple_size<Tuple>::type; } &&
    std::derived_from<std::tuple_size<Tuple>,
                      std::integral_constant<size_t, std::tuple_size_v<Tuple>>>;

// Checks if `std::tuple_element_t` is applicable.
template<class Tuple, size_t Index>
concept has_tuple_element = requires(Tuple tuple) {
  typename std::tuple_element_t<Index, std::remove_const_t<Tuple>>;
  // clang-format off
  { std::get<Index>(tuple) } -> std::convertible_to<
      const std::tuple_element_t<Index, std::remove_const_t<Tuple>>&>;
  // clang-format on
};

// Checks if `std::tuple_element_t` is applicable.
template<has_tuple_size Tuple>
inline constexpr bool has_tuple_elements_v =
    base_apply<std::tuple_size_v<Tuple>>([](auto... indices) {
      return (has_tuple_element<Tuple, indices> && ...);
    });

// Check if type is tuple-like
// (`std::array`, `std::pair`, `std::tuple`, ...).
template<class Tuple>
concept tuple_like = (!std::is_reference_v<Tuple>) && //
                     has_tuple_size<Tuple> &&
                     has_tuple_elements_v<Tuple>;

// Wrap an argument into a tuple, if it is not a tuple already.
class ToTuple final {
public:

  template<class Tuple>
    requires tuple_like<std::remove_cvref_t<Tuple>>
  constexpr decltype(auto) operator()(Tuple&& tuple) const {
    return std::forward<Tuple>(tuple);
  }
  template<class... Args>
  constexpr auto operator()(Args&&... args) const {
    return std::tuple{std::forward<Args>(args)...};
  }

}; // class ToTuple

// Wrap an argument into a tuple, if it is not a tuple already.
inline constexpr ToTuple to_tuple{};

// -----------------------------------------------------------------------------

// Same as `std::invocable`, but for `std::apply` instead of `std::invoke`.
template<class Func, class Tuple>
concept applicable =
    tuple_like<Tuple> && //
    requires(Func& func, Tuple& tuple) { std::apply(func, tuple); };

// Same as `std::regular_invocable`, but for `std::apply`.
template<class Func, class Tuple>
concept regular_applicable = applicable<Func, Tuple>;

// Same as `std::invoke_result_t`, but for `std::apply`.
template<class Func, tuple_like Tuple>
  requires applicable<Func, Tuple>
using apply_result_t =
    decltype(std::apply(std::declval<Func>(), std::declval<Tuple>()));

// Drop the first N tuple elements and pass the rest to the callback.
template<size_t N, class Func, class Tuple>
  requires tuple_like<std::remove_cvref_t<Tuple>>
constexpr decltype(auto) drop_n_apply(Func&& func, Tuple&& tuple) {
  constexpr size_t Size = std::tuple_size_v<std::remove_cvref_t<Tuple>>;
  static_assert(N <= Size);
  return base_apply<Size - N>([&](auto... indices) -> decltype(auto) {
    return func(std::get<indices + N>(tuple)...);
  });
}
template<size_t N, class Func, class... Args>
constexpr decltype(auto) drop_n(Func&& func, Args&&... args) {
  auto args_tuple = std::forward_as_tuple(std::forward<Args>(args)...);
  return drop_n_apply<N>(func, args_tuple);
}

// Pack the first N tuple elements and pass them and the rest to the callback.
template<size_t N, class Func, class Tuple>
  requires tuple_like<std::remove_cvref_t<Tuple>>
constexpr decltype(auto) pack_n_apply(Func&& func, Tuple&& tuple) {
  constexpr size_t Size = std::tuple_size_v<std::remove_cvref_t<Tuple>>;
  static_assert(N <= Size);
  auto packed = base_apply<N>(
      [&](auto... indices) { return std::tuple{std::get<indices>(tuple)...}; });
  return base_apply<Size - N>([&](auto... indices) -> decltype(auto) {
    return func(std::move(packed), std::get<indices + N>(tuple)...);
  });
}
template<size_t N, class Func, class... Args>
constexpr decltype(auto) pack_n(Func&& func, Args&&... args) {
  auto args_tuple = std::forward_as_tuple(std::forward<Args>(args)...);
  return pack_n_apply<N>(func, args_tuple);
}

// Apply chunks of size N of tuple elements to a function.
template<size_t N, class Func = ToVoid, class Chunk = ToTuple, class Tuple>
  requires tuple_like<std::remove_cvref_t<Tuple>>
constexpr decltype(auto) chunk_n_apply(Func&& func, Chunk&& chunk,
                                       Tuple&& tuple) {
  constexpr size_t Size = std::tuple_size_v<std::remove_cvref_t<Tuple>>;
  constexpr size_t NumChunks = Size / N;
  return base_apply<NumChunks>([&](auto... indices) -> decltype(auto) {
    const auto chunk_at = [&](auto index) -> decltype(auto) {
      return base_apply<N>([&](auto... subindices) -> decltype(auto) {
        return chunk(std::get<index * N + subindices>(tuple)...);
      });
    };
    return func(chunk_at(indices)...);
  });
}
template<size_t N, class Func = ToVoid, class Chunk = ToTuple, class... Args>
constexpr decltype(auto) chunk_n(Func&& func, Chunk&& chunk, Args&&... args) {
  auto args_tuple = std::forward_as_tuple(std::forward<Args>(args)...);
  return chunk_n_apply<N>(func, chunk, args_tuple);
}

// Apply zips of tuple elements to a function.
template<class Func = ToVoid, class Zip = ToTuple, class... Tuples>
  requires (... && tuple_like<std::remove_cvref_t<Tuples>>)
constexpr decltype(auto) zip_apply(Func&& func, Zip&& zip, Tuples&&... tuples) {
  constexpr size_t MinSize =
      std::min({std::tuple_size_v<std::remove_cvref_t<Tuples>>...});
  return base_apply<MinSize>([&](auto... indices) -> decltype(auto) {
    const auto zip_at = [&](auto index) -> decltype(auto) {
      return zip(std::get<index>(tuples)...);
    };
    return func(zip_at(indices)...);
  });
}

// Apply N copies of the argument to a function.
template<size_t N, class Func, class Arg>
constexpr decltype(auto) n_times(Func&& func, Arg&& arg) {
  return base_apply<N>([&](auto... indices) -> decltype(auto) {
    const auto get_at = [&](auto /*index*/) -> decltype(auto) { return arg; };
    return func(get_at(indices)...);
  });
}

// -----------------------------------------------------------------------------

/// @todo Implement `operator<=>`.
#if 0
class SynthThreeWay final {
public:

  template<class Arg1, class Arg2>
  constexpr auto operator()(const Arg1& arg1, const Arg2& arg2) const
    requires requires {
      { arg1 < arg2 } -> boolean_testable;
      { arg2 < arg1 } -> boolean_testable;
    }
  {
    if constexpr (std::three_way_comparable_with<Arg1, Arg2>) {
      return arg1 <=> arg2;
    } else {
      if (arg1 < arg2) return std::weak_ordering::less;
      else if (arg2 < arg1) return std::weak_ordering::greater;
      else return std::weak_ordering::equivalent;
    }
  }

}; // class SynthThreeWay

constexpr inline SynthThreeWay synth_three_way{};

template<class Arg1, class Arg2 = Arg1>
using synth_three_way_t = decltype( //
    synth_three_way(std::declval<Arg1>(), std::declval<Arg2>()));
#endif

template<class Tuple1, class Tuple2>
  requires tuple_like<std::remove_cvref_t<Tuple1>> &&
           tuple_like<std::remove_cvref_t<Tuple2>>
constexpr auto operator==(Tuple1&& tuple1, Tuple2&& tuple2) {
  const auto compare = [](auto&&... pairs) {
    const auto equals = [](auto&& pair) {
      auto&& [first, second] = pair;
      return first == second;
    };
    return (equals(pairs) && ...);
  };
  return zip_apply(compare, {}, tuple1, tuple2);
}

// -----------------------------------------------------------------------------

} // namespace Storm
