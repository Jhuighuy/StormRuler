////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2022 Oleg Butakov                                            //
//                                                                            //
// Permission is hereby granted, free of charge, to any person obtaining a    //
// copy of this software and associated documentation files (the "Software"), //
// to deal in the Software without restriction, including without limitation  //
// the rights to use, copy, modify, merge, publish, distribute, sublicense,   //
// and/or sell copies of the Software, and to permit persons to whom the      //
// Software is furnished to do so, subject to the following conditions:       //
//                                                                            //
// The above copyright notice and this permission notice shall be included in //
// all copies or substantial portions of the Software.                        //
//                                                                            //
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR //
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   //
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    //
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER //
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING    //
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER        //
// DEALINGS IN THE SOFTWARE.                                                  //
////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <Storm/Base.hpp>

#include <stdexcept>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>

namespace Storm::meta {

////////////////////////////////////////////////////////////////////////////////
//                                 Base types                                 //
////////////////////////////////////////////////////////////////////////////////

struct empty_t {};
inline constexpr empty_t empty{};

template<class T>
struct type {};

template<class T>
struct raw {
  using type = T;
}; // struct raw
template<class T>
struct raw<type<T>> {
  using type = T;
}; // struct raw<type<T>>

template<class T>
using raw_t = typename raw<T>::type;

template<class T>
using type_t = type<raw_t<T>>;

template<class T>
inline constexpr type_t<T> type_v{};

template<class... Ts>
struct list {};

template<class... Ts>
using list_t = list<raw_t<Ts>...>;

template<class... Ts>
inline constexpr list_t<Ts...> list_v{};

/// @todo We need `meta::box<X, T = meta::empty>` and
/// `meta::copyable_box<X, T = meta::empty>` instead.
template<class X, class T>
using tagged_std_pair_t = std::pair<raw_t<X>, type_t<T>>;

////////////////////////////////////////////////////////////////////////////////
//                                 Generators                                 //
////////////////////////////////////////////////////////////////////////////////

template<template<size_t> class Type, size_t N>
consteval auto make_list() {
  // clang-format off
  return []<size_t... I>(std::integer_sequence<size_t, I...>) {
    return list_v<Type<I>...>;
  }(std::make_integer_sequence<size_t, N>{});
  // clang-format on
}
template<template<size_t> class Type, size_t N>
using make_list_t = decltype(make_list<Type, N>());

////////////////////////////////////////////////////////////////////////////////
//                                  Queries                                   //
////////////////////////////////////////////////////////////////////////////////

template<class... Ts>
consteval size_t size(list<Ts...>) {
  return sizeof...(Ts);
}
template<class List>
inline constexpr size_t size_v = size(List{});

template<class T, class... Ts>
consteval auto first(list<T, Ts...>) {
  return type_v<T>;
}
template<class List>
using first_t = decltype(first(List{}));

template<class T1, class T2, class... Ts>
consteval auto second(list<T1, T2, Ts...>) {
  return type_v<T2>;
}
template<class List>
using second_t = decltype(second(List{}));

template<class T, class... Ts>
consteval auto rest(list<T, Ts...>) {
  return list_v<Ts...>;
}
template<class List>
using rest_t = decltype(rest(List{}));

template<class X, class... Ts>
consteval bool contains(list<Ts...>) {
  using R = raw_t<X>;
  return (std::is_same_v<R, Ts> || ...);
}
template<class X, class List>
inline constexpr bool contains_v = contains<X>(List{});

consteval auto all_unique(list<>) {
  return true;
}
template<class T, class... Ts>
consteval auto all_unique(list<T, Ts...>) {
  const auto rest = list_v<Ts...>;
  if constexpr (!contains<T>(rest)) {
    return all_unique(rest);
  } else {
    return false;
  }
}
template<class List>
inline constexpr bool all_unique_v = all_unique(List{});

////////////////////////////////////////////////////////////////////////////////
//                                 Algorithms                                 //
////////////////////////////////////////////////////////////////////////////////

template<class Func, class... Ts>
consteval auto transform(Func func, list<Ts...>) {
  return list_v<decltype(func(std::declval<Ts>()))...>;
}
template<class Func, class List>
using transform_t = decltype(transform(Func{}, List{}));

template<class X, class... Ts>
consteval auto append(list<Ts...>) {
  return list_v<Ts..., raw_t<X>>;
}
template<class X, class List>
using append_t = decltype(append<X>(List{}));

template<class X, class... Ts>
consteval auto prepend(list<Ts...>) {
  return list_v<raw_t<X>, Ts...>;
}
template<class X, class List>
using prepend_t = decltype(prepend<X>(List{}));

template<class... Ts>
consteval auto concat(list<Ts...>) {
  return list_v<Ts...>;
}
template<class... Ts, class... Us>
consteval auto concat(list<Ts...>, //
                      list<Us...>, auto... rest) {
  return concat(list_v<Ts..., Us...>, rest...);
}
template<class... Lists>
using concat_t = decltype(concat(Lists{}...));

struct reverse_fn {
  consteval auto operator()(list<>) const {
    return list_v<>;
  }
  template<class T, class... Ts>
  consteval auto operator()(list<T, Ts...>) const {
    return append<T>((*this)(list_v<Ts...>) );
  }
}; // struct reverse_fn
inline constexpr reverse_fn reverse;
template<class List>
using reverse_t = decltype(reverse(List{}));

consteval auto unique(list<>) {
  return list_v<>;
}
template<class T, class... Ts>
consteval auto unique(list<T, Ts...>) {
  const auto rest_unique = unique(list_v<Ts...>);
  if constexpr (contains<T>(rest_unique)) {
    return rest_unique;
  } else {
    return prepend<T>(rest_unique);
  }
}
template<class List>
using unique_t = decltype(unique(List{}));

template<class X, class... Ts>
consteval auto pair_list(list<Ts...>) {
  using R = raw_t<X>;
  return list_v<list_t<R, Ts>...>;
}
template<class X, class List>
using pair_list_t = decltype(pair_list<X>(List{}));

template<class... Ts, class... Us>
consteval auto cartesian_product(list<Ts...>, list<Us...> = list_v<Ts...>) {
  return concat(pair_list<Ts>(list_v<Us...>)...);
}
template<class List1, class List2 = List1>
using cartesian_product_t = decltype(cartesian_product(List1{}, List2{}));

////////////////////////////////////////////////////////////////////////////////
//                                  Casting                                   //
////////////////////////////////////////////////////////////////////////////////

template<template<class, class> class ToPair>
struct pair_cast_fn {
  template<class T, class U>
  consteval auto operator()(list<T, U>) const {
    return type_v<ToPair<T, U>>;
  }
}; // struct pair_cast_fn

template<template<class, class> class ToPair>
inline constexpr pair_cast_fn<ToPair> pair_cast{};

template<template<class, class> class ToPair, class Pair>
using pair_cast_t = decltype(pair_cast<ToPair>(Pair{}));

template<class Pair>
using as_std_pair_t = raw_t<pair_cast_t<std::pair, Pair>>;

template<template<class...> class ToList>
struct list_cast_fn {
  template<class... Ts>
  consteval auto operator()(list<Ts...>) const {
    return type_v<ToList<Ts...>>;
  }
}; // struct list_cast_fn

template<template<class...> class ToList>
inline constexpr list_cast_fn<ToList> list_cast{};

template<template<class...> class ToList, class List>
using list_cast_t = decltype(list_cast<ToList>(List{}));

template<class List>
using as_std_tuple_t = raw_t<list_cast_t<std::tuple, List>>;

////////////////////////////////////////////////////////////////////////////////
//                              Other utilities                               //
////////////////////////////////////////////////////////////////////////////////

/// @brief Get the template parameter name as a string.
/// @tparam T Template parameter.
template<class T>
consteval auto type_name() {
  std::string_view name, prefix, suffix;
#if STORM_COMPILER_GCC_
  name = __PRETTY_FUNCTION__;
  prefix = "consteval auto Storm::meta::type_name() [with T = ";
  suffix = "]";
#elif STORM_COMPILER_CLANG_
  name = __PRETTY_FUNCTION__;
  prefix = "auto Storm::meta::type_name() [T = ";
  suffix = "]";
#elif STORM_COMPILER_MSVC_
  name = __FUNCSIG__;
  prefix = "auto __cdecl Storm::meta::type_name<";
  suffix = ">(void)";
#endif
  if (!name.starts_with(prefix)) { throw std::range_error("Invalid prefix!"); }
  if (!name.ends_with(suffix)) { throw std::range_error("Invalid suffix!"); }
  name.remove_prefix(prefix.size());
  name.remove_suffix(suffix.size());
  return name;
}

template<class T>
consteval auto type_name(T) {
  return type_name<T>();
}

template<class T>
inline constexpr std::string_view type_name_v = type_name<T>();

} // namespace Storm::meta
