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

#include <Storm/Crow/ConceptUtils.hpp>
#include <Storm/Crow/TupleUtils.hpp>

#include <array>
#include <utility>

namespace Storm {

// -----------------------------------------------------------------------------

/// @brief Check if type is array-like.
template<class Array>
concept array_like = (!std::is_reference_v<Array>) && //
                     has_tuple_size<Array> &&
                     requires(Array& array, size_t index) {
                       { array[index] } -> can_reference;
                     };

namespace detail {

  // Convert index sequence to an array of indices.
  /// @todo This function should be rethinked.
  template<class Index, Index... Indices>
  consteval auto _to_array(std::integer_sequence<Index, Indices...>) {
    return std::array{Indices...};
  }

  // Convert tuple of indices to an index sequence.
  /// @todo This function should be rethinked.
  template<tuple_like auto IndexTuple>
  consteval auto _to_index_sequence() {
    using Tuple = decltype(IndexTuple);
    return ([]<size_t... Indices>(std::index_sequence<Indices...>) {
      return std::index_sequence< //
          static_cast<size_t>(std::get<Indices>(IndexTuple))...>{};
    }(std::make_index_sequence<std::tuple_size_v<Tuple>>{}));
  }

} // namespace detail

// -----------------------------------------------------------------------------

} // namespace Storm
