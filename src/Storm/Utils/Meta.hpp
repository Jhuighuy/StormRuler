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
#include <type_traits>
#include <utility>

namespace Storm::meta {

////////////////////////////////////////////////////////////////////////////////
//                                 Base types                                 //
////////////////////////////////////////////////////////////////////////////////

template<class T>
struct type {};

template<class T>
struct raw {
  using type = T;
};
template<class T>
struct raw<type<T>> {
  using type = T;
};

template<class T>
using raw_t = typename raw<T>::type;

template<class T>
using type_t = type<raw_t<T>>;

template<class T>
inline constexpr type_t<T> type_v{};

////////////////////////////////////////////////////////////////////////////////
//                                 Pair types                                 //
////////////////////////////////////////////////////////////////////////////////

template<class T, class U>
struct type_pair {};

template<class T, class U>
using type_pair_t = type_pair<raw_t<T>, raw_t<U>>;

template<class T, class U>
inline constexpr type_pair_t<T, U> type_pair_v{};

template<template<class, class> class ToPair>
struct pair_caster {
  template<class T, class U>
  consteval auto operator()(type_pair<T, U>) const {
    return type_v<ToPair<T, U>>;
  }
}; // struct pair_caster

template<template<class, class> class ToPair>
inline constexpr pair_caster<ToPair> pair_cast{};

template<template<class, class> class ToPair, class Pair>
using pair_cast_t = decltype(pair_cast<ToPair>(Pair{}));

////////////////////////////////////////////////////////////////////////////////
//                                 List types                                 //
////////////////////////////////////////////////////////////////////////////////

template<class... Ts>
struct type_list {};

template<class... Ts>
using type_list_t = type_list<raw_t<Ts>...>;

template<class... Ts>
inline constexpr type_list_t<Ts...> type_list_v{};

template<template<class...> class ToList>
struct list_caster {
  template<class... Ts>
  consteval auto operator()(type_list<Ts...>) const {
    return type_v<ToList<Ts...>>;
  }
}; // struct list_caster

template<template<class...> class ToList>
inline constexpr list_caster<ToList> list_cast{};

template<template<class...> class ToList, class List>
using list_cast_t = decltype(list_cast<ToList>(List{}));

template<class... Ts>
consteval size_t size(type_list<Ts...>) {
  return sizeof...(Ts);
}

template<class List>
inline constexpr size_t size_v = size(List{});

template<class T, class... Ts>
consteval auto head(type_list<T, Ts...>) {
  return type_v<T>;
}

template<class List>
using head_t = decltype(head(List{}));

template<class T, class... Ts>
consteval auto tail(type_list<T, Ts...>) {
  return type_list_v<Ts...>;
}

template<class List>
using tail_t = decltype(tail(List{}));

template<class X, class... Ts>
consteval bool contains(type_list<Ts...>) {
  using R = raw_t<X>;
  return (std::is_same_v<R, Ts> || ...);
}

template<class X, class List>
inline constexpr bool contains_v = contains<X>(List{});

template<class Func, class... Ts>
consteval auto transform(Func func, type_list<Ts...>) {
  return type_list_v<decltype(func(std::declval<Ts>()))...>;
}

template<class Func, class List>
using transform_t = decltype(transform(Func{}, List{}));

template<class X, class... Ts>
consteval auto append(type_list<Ts...>) {
  return type_list_v<Ts..., raw_t<X>>;
}

template<class X, class List>
using append_t = decltype(append<X>(List{}));

template<class X, class... Ts>
consteval auto prepend(type_list<Ts...>) {
  return type_list_v<raw_t<X>, Ts...>;
}

template<class X, class List>
using prepend_t = decltype(prepend<X>(List{}));

consteval auto concat() {
  return type_list_v<>;
}
template<class... Ts>
consteval auto concat(type_list<Ts...>) {
  return type_list_v<Ts...>;
}
template<class... Ts, class... Us>
consteval auto concat(type_list<Ts...>, //
                      type_list<Us...>, auto... rest) {
  return concat(type_list_v<Ts..., Us...>, rest...);
}

template<class... Lists>
using concat_t = decltype(concat(Lists{}...));

consteval auto unique(type_list<>) {
  return type_list_v<>;
}
template<class T, class... Ts>
consteval auto unique(type_list<T, Ts...>) {
  auto rest_unique = unique(type_list_v<Ts...>);
  if constexpr (contains<T>(rest_unique)) {
    return rest_unique;
  } else {
    return prepend<T>(rest_unique);
  }
}

template<class List>
using unique_t = decltype(unique(List{}));

template<class X, class... Ts>
using pair_list_t = type_list_t<type_pair_t<X, Ts>...>;

template<class X, class... Ts>
inline constexpr pair_list_t<X, Ts...> pair_list_v{};

template<class... Ts, class... Us>
consteval auto cartesian_product(type_list<Ts...>, type_list<Us...>) {
  return concat(pair_list_v<Ts, Us...>...);
}

template<class List1, class List2>
using cartesian_product_t = decltype(cartesian_product(List1{}, List2{}));

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
