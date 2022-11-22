////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2020-2023 Oleg Butakov                                       //
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

///////////////////////////////////
// Empty type

struct empty_t {};

inline constexpr empty_t empty_v{};

//
///////////////////////////////////

///////////////////////////////////
// Type wrapper

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

//
///////////////////////////////////

///////////////////////////////////
// Type list

template<class... Ts>
struct list {};

template<class... Ts>
using list_t = list<raw_t<Ts>...>;

template<class... Ts>
inline constexpr list_t<Ts...> list_v{};

//
///////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                                 Generators                                 //
////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////
// Sequence generator.

template<template<size_t> class Type, size_t First, size_t Last>
consteval auto make_seq() {
  if constexpr (First >= Last) {
    return list_v<>;
  } else {
    // clang-format off
    return []<size_t... I>(std::integer_sequence<size_t, I...>) {
      return list_v<Type<First + I>...>;
    }(std::make_integer_sequence<size_t, Last - First>{});
    // clang-format on
  }
}

template<template<size_t> class Type, size_t First, size_t Last>
using make_seq_t = decltype(make_seq<Type, First, Last>());

//
///////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                                  Queries                                   //
////////////////////////////////////////////////////////////////////////////////

template<class Fn>
struct query_fn {};

///////////////////////////////////
// Size query

template<class... Ts>
consteval size_t size(list<Ts...>) {
  return sizeof...(Ts);
}
template<class List>
inline constexpr size_t size_v = size(List{});

//
///////////////////////////////////

///////////////////////////////////
// List item queries

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

//
///////////////////////////////////

///////////////////////////////////
// Sublist queries

template<class T, class... Ts>
consteval auto drop_first(list<T, Ts...>) {
  return list_v<Ts...>;
}

template<class List>
using drop_first_t = decltype(drop_first(List{}));

//
///////////////////////////////////

///////////////////////////////////
// Contains query.

template<class X, class... Ts>
consteval bool contains(list<Ts...>) {
  using R = raw_t<X>;
  return (std::is_same_v<R, Ts> || ...);
}

template<class X, class List>
inline constexpr bool contains_v = contains<X>(List{});

//
///////////////////////////////////

///////////////////////////////////
// All unique query

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

//
///////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                              Transformations                               //
////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////
// For each item in list

template<class... Ts, class Fn>
constexpr void for_each(list<Ts...>, [[maybe_unused]] Fn fn) {
  (fn(type_v<Ts>), ...);
}
template<class List, class Fn>
constexpr void for_each(Fn fn) {
  for_each(List{}, fn);
}

//
///////////////////////////////////

///////////////////////////////////
// Transform the type list

template<class Fn>
struct transform_fn {};

template<class Fn, class... Ts>
consteval auto transform(transform_fn<Fn>, list<Ts...>) {
  return list_v<decltype(Fn{}(std::declval<Ts>()))...>;
}

template<class Func, class List>
using transform_t = decltype(transform(Func{}, List{}));

//
///////////////////////////////////

///////////////////////////////////
// Append to the type list

template<class X>
struct append_fn : transform_fn<append_fn<X>> {
  template<class... Ts>
  consteval auto operator()(list<Ts...>) const {
    return list_v<Ts..., raw_t<X>>;
  }
}; // struct append_fn

template<class X>
inline constexpr append_fn<X> append{};

template<class X, class List>
using append_t = decltype(append<X>(List{}));

//
///////////////////////////////////

///////////////////////////////////
// Prepend to the type list

template<class X>
struct prepend_fn : transform_fn<prepend_fn<X>> {
  template<class... Ts>
  consteval auto operator()(list<Ts...>) const {
    return list_v<raw_t<X>, Ts...>;
  }
}; // struct prepend_fn

template<class X>
inline constexpr prepend_fn<X> prepend{};

template<class X, class List>
using prepend_t = decltype(prepend<X>(List{}));

//
///////////////////////////////////

///////////////////////////////////
// Reverse the type list

struct reverse_fn : transform_fn<reverse_fn> {
  consteval auto operator()(list<>) const {
    return list_v<>;
  }
  template<class T, class... Ts>
  consteval auto operator()(list<T, Ts...>) const {
    return append<T>((*this)(list_v<Ts...>) );
  }
}; // struct reverse_fn

inline constexpr reverse_fn reverse{};

template<class List>
using reverse_t = decltype(reverse(List{}));

//
///////////////////////////////////

///////////////////////////////////
// Unique types list

struct unique_fn : transform_fn<unique_fn> {
  consteval auto operator()(list<>) const {
    return list_v<>;
  }
  template<class T, class... Ts>
  consteval auto operator()(list<T, Ts...>) const {
    const auto rest_unique = unique(list_v<Ts...>);
    if constexpr (contains<T>(rest_unique)) {
      return rest_unique;
    } else {
      return prepend<T>(rest_unique);
    }
  }
}; // struct unique_fn

inline constexpr unique_fn unique{};

template<class List>
using unique_t = decltype(unique(List{}));

//
///////////////////////////////////

///////////////////////////////////
// Make pair list

template<class X>
struct l_pair_list_fn {
  template<class... Ts>
  consteval auto operator()(list<Ts...>) const {
    using R = raw_t<X>;
    return list_v<list_t<R, Ts>...>;
  }
}; // struct l_pair_list_fn

template<class X>
struct r_pair_list_fn {
  template<class... Ts>
  consteval auto operator()(list<Ts...>) const {
    using R = raw_t<X>;
    return list_v<list_t<Ts, R>...>;
  }
}; // struct r_pair_list_fn

template<class X>
inline constexpr l_pair_list_fn<X> l_pair_list;
template<class X>
inline constexpr r_pair_list_fn<X> r_pair_list;

template<class X, class List>
using l_pair_list_t = decltype(l_pair_list<X>(List{}));
template<class List, class X>
using r_pair_list_t = decltype(r_pair_list<X>(List{}));

template<class, class>
struct pair_list;
template<class X, class... Ts>
struct pair_list<X, list<Ts...>> {
  using type = l_pair_list_t<X, list<Ts...>>;
};
template<class... Ts, class X>
struct pair_list<list<Ts...>, X> {
  using type = r_pair_list_t<list<Ts...>, X>;
};

template<class X, class Y>
using pair_list_t = typename pair_list<X, Y>::type;

template<class X, class Y>
inline constexpr pair_list_t<X, Y> pair_list_v{};

//
///////////////////////////////////

///////////////////////////////////
// Pair cast

template<template<class...> class ToPair>
struct pair_cast_fn : transform_fn<pair_cast_fn<ToPair>> {
  template<class T, class U>
  consteval auto operator()(list<T, U>) const {
    return type_v<ToPair<T, U>>;
  }
}; // struct pair_cast_fn

template<template<class...> class ToPair>
inline constexpr pair_cast_fn<ToPair> pair_cast{};

template<template<class...> class ToPair, class Pair>
using pair_cast_t = decltype(pair_cast<ToPair>(Pair{}));

template<class Pair>
using as_std_pair_t = raw_t<pair_cast_t<std::pair, Pair>>;

//
///////////////////////////////////

///////////////////////////////////
// List cast

template<template<class...> class ToList>
struct list_cast_fn : transform_fn<list_cast_fn<ToList>> {
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

//
///////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                                 Reductions                                 //
////////////////////////////////////////////////////////////////////////////////

template<class Fn>
struct reduce_fn {};

///////////////////////////////////
// Concatenate the type lists.

struct concat_fn : reduce_fn<concat_fn> {
  template<class... Ts>
  consteval auto operator()(list<Ts...>) const {
    return list_v<Ts...>;
  }
  template<class... Ts, class... Us>
  consteval auto operator()(list<Ts...>, list<Us...>, auto... rest) const {
    return (*this)(list_v<Ts..., Us...>, rest...);
  }
}; // struct concat_fn

inline constexpr concat_fn concat{};

template<class... Lists>
using concat_t = decltype(concat(Lists{}...));

//
///////////////////////////////////

///////////////////////////////////
// Cartesian product of the type lists

struct cartesian_product_fn : reduce_fn<cartesian_product_fn> {
  template<class... Ts, class... Us>
  consteval auto operator()(list<Ts...>, list<Us...>) const {
    return concat(l_pair_list<Ts>(list_v<Us...>)...);
  }
}; // struct cartesian_product_fn

inline constexpr cartesian_product_fn cartesian_product;

template<class... Lists>
using cartesian_product_t = decltype(cartesian_product(Lists{}...));

//
///////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                              Other utilities                               //
////////////////////////////////////////////////////////////////////////////////

template<class>
inline constexpr bool always_false = false;

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
  if (!name.starts_with(prefix)) throw std::range_error("Invalid prefix!");
  if (!name.ends_with(suffix)) throw std::range_error("Invalid suffix!");
  name.remove_prefix(prefix.size());
  name.remove_suffix(suffix.size());
  return name;
}

template<class T>
constexpr auto type_name(T&&) {
  return type_name<T>();
}

template<class T>
inline constexpr std::string_view type_name_v = type_name<T>();

} // namespace Storm::meta
