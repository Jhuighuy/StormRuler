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

#include <concepts>
#include <type_traits>

namespace Storm {

// -----------------------------------------------------------------------------

/// @brief A semi-valid argument for the CRTP classes.
template<class Derived>
concept crtp_derived = std::is_class_v<Derived> &&
                       std::same_as<Derived, std::remove_cv_t<Derived>>;

// Based on /*is-derived-from-view-interface*/<T> implementation.
namespace detail_ {
  template<template<crtp_derived Derived> class CrtpInterface>
  struct DerivedFromCrtpInterfaceImpl_ {
    template<class T, class U>
    DerivedFromCrtpInterfaceImpl_(const T&,
                                  const CrtpInterface<U>&); // not defined
  };
} // namespace detail_

/// @brief Checks if and only if @c T has exactly one public base class
/// @c CrtpInterface&lt;U&gt; for some type @c U, and @c T has no base classes
/// of type @c CrtpInterface&lt;V&gt; for any other type @c V.
template<class T, template<class Derived> class CrtpInterface>
concept derived_from_crtp_interface = requires(T x) { //
  detail_::DerivedFromCrtpInterfaceImpl_<CrtpInterface>{x, x};
};

// -----------------------------------------------------------------------------

} // namespace Storm
