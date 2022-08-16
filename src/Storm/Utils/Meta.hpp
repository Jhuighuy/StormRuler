/// Copyright (C) 2022 Oleg Butakov
///
/// Permission is hereby granted, free of charge, to any person obtaining a copy
/// of this software and associated documentation files (the "Software"), to
/// deal in the Software without restriction, including without limitation the
/// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
/// sell copies of the Software, and to permit persons to whom the Software is
/// furnished to do so, subject to the following conditions:
///
/// The above copyright notice and this permission notice shall be included in
/// all copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
/// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
/// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
/// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
/// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
/// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
/// IN THE SOFTWARE.

#pragma once

#include <Storm/Base.hpp>

#include <stdexcept>
#include <string_view>

namespace Storm::meta {

/// @brief Type wrapper.
template<class T>
struct type_t {};

/// @brief Type wrapper.
template<class T>
inline constexpr type_t<T> type{};

/// @brief Type list wrapper.
template<class... Ts>
struct type_list_t {};

/// @brief Type list wrapper.
template<class... Ts>
inline constexpr type_list_t<Ts...> type_list{};

/// @brief Concatenate the type lists.
/// @{
consteval auto concat() {
  return type_list<>;
}
template<class... Ts>
consteval auto concat(type_list_t<Ts...>) {
  return type_list<Ts...>;
}
template<class... Ts, class... Us>
consteval auto concat(type_list_t<Ts...>, //
                      type_list_t<Us...>, auto... rest) {
  return concat(type_list<Ts..., Us...>, rest...);
}
/// @}

/// @brief Concatenate the type lists.
template<class... Lists>
using concat_t = decltype(concat(Lists{}...));

/// @brief Type pair.
template<class T, class U>
struct type_pair_t {};

/// @brief Type pair.
template<class T, class U>
inline constexpr type_pair_t<T, U> type_pair{};

/// @brief Generate a `type_list_t<type_pair_t<X, T1>, type_pair_t<X, T2>, ...>`
template<class X, class... Ts>
using pair_list_t = type_list_t<type_pair_t<X, Ts>...>;

/// @brief Generate a `type_list<type_pair_t<X, T1>, type_pair_t<X, T2>, ...>`
template<class X, class... Ts>
inline constexpr pair_list_t<X, Ts...> pair_list{};

/// @brief For given `T1, ..., TN`, `U1, ..., UM` generate:
/// `type_list<
///   type_pair_t<T1, U1>, type_pair_t<T1, U2>, ..., type_pair_t<T1, UM>,
///   type_pair_t<T2, U1>, type_pair_t<T2, U2>, ..., type_pair_t<T2, UM>,
///   ...
///   type_pair_t<TN, U1>, type_pair_t<TN, U2>, ..., type_pair_t<TN, UM>>`.
/// @{
template<class... Ts, class... Us>
consteval auto cartesian_product(type_list_t<Ts...>, type_list_t<Us...>) {
  return concat(pair_list<Ts, Us...>...);
}
/// @}

/// @brief For given `List1 = type_list_t<T1, ..., TN>`,
/// `List2 = type_list_t<U1, ..., UM>` generate:
/// `type_list_t<
///   type_pair_t<T1, U1>, type_pair_t<T1, U2>, ..., type_pair_t<T1, UM>,
///   type_pair_t<T2, U1>, type_pair_t<T2, U2>, ..., type_pair_t<T2, UM>,
///   ...
///   type_pair_t<TN, U1>, type_pair_t<TN, U2>, ..., type_pair_t<TN, UM>>`.
template<class List1, class List2>
using cartesian_product_t = decltype(cartesian_product(List1{}, List2{}));

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
  if (!name.starts_with(prefix)) {
    throw std::runtime_error("Invalid prefix!");
  }
  if (!name.ends_with(suffix)) { //
    throw std::runtime_error("Invalid suffix!");
  }
  name.remove_prefix(prefix.size());
  name.remove_suffix(suffix.size());
  return name;
}

template<class T>
consteval auto type_name(T) {
  return type_name<T>();
}

} // namespace Storm::meta
