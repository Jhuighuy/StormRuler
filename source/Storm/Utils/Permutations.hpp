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

#include <array>
#include <concepts>
#include <iterator>
#include <utility>

/// @brief Permutation utilities.
///
/// Permutation can be used as a mapping from the permuted to the unpermuted
/// indices of a collection. Common usage pattern:
/// @code
/// std::vector<T> values = {...};
/// std::vector<T> permuted_values(values.size());
/// for (size_t i = 0; i < permuted_values.size(); ++i) {
///   permuted_values[i] = values[permutation[i]];
/// }
/// @endcode
///
/// While the common usage pattern for the inverse permutations:
/// @code
/// std::vector<T> values = {...};
/// std::vector<T> permuted_values(values.size());
/// for (size_t j = 0; j < values.size(); ++j) {
///   permuted_values[inverse_permutation[i]] = values[i];
/// }
/// @endcode
namespace Storm::permutations {

// -----------------------------------------------------------------------------

/// @brief Check is specified range is a permutation.
template<class PermutationRange>
constexpr bool is_permutation(PermutationRange&& perm_range) noexcept {
  if (perm_range.size() == 1) {
    return perm_range[0] == 0;
  }
  if (perm_range.size() == 2) {
    return (perm_range[0] == 0 && perm_range[1] == 1) ||
           (perm_range[0] == 1 && perm_range[1] == 0);
  }
  STORM_ABORT("`is_permutation` is not completely implemented yet!");
}

// -----------------------------------------------------------------------------

/// @brief Invert the permutation. Complexity is linear.
/// @todo Add contrains!
/// @{
template<std::input_iterator PermutationIterator,
         std::sentinel_for<PermutationIterator> PermutationSentinel,
         std::random_access_iterator InversePermutationIterator>
  requires std::indirectly_copyable<PermutationIterator,
                                    InversePermutationIterator>
constexpr void invert_permutation(PermutationIterator perm_iterator,
                                  PermutationSentinel perm_sentinel,
                                  InversePermutationIterator iperm_iterator) {
  using Index = std::iter_value_t<PermutationIterator>;
  Index index{};
  for (; perm_iterator != perm_sentinel; ++perm_iterator, ++index) {
    STORM_ASSERT(*perm_iterator != Index{SIZE_MAX},
                 "Invalid permutation iterator value!");
    iperm_iterator[static_cast<size_t>(*perm_iterator)] = index;
  }
}
template<std::ranges::input_range PermutationRange,
         std::random_access_iterator InversePermutationIterator>
  requires std::indirectly_copyable<std::ranges::iterator_t<PermutationRange>,
                                    InversePermutationIterator>
constexpr void invert_permutation(PermutationRange&& perm_range,
                                  InversePermutationIterator iperm_iterator) {
  invert_permutation(std::ranges::begin(perm_range),
                     std::ranges::end(perm_range), std::move(iperm_iterator));
}
template<class Index, size_t Size>
constexpr auto
invert_permutation(const std::array<Index, Size>& perm) noexcept {
  std::array<Index, Size> iperm{};
  invert_permutation(perm, iperm.begin());
  return iperm;
}
/// @}

/// @brief Apply the permutation inplace. Complexity is linear.
/// @{
template<std::random_access_iterator PermutationIterator,
         std::sentinel_for<PermutationIterator> PermutationSentinel,
         std::invocable<std::iter_value_t<PermutationIterator>,
                        std::iter_value_t<PermutationIterator>>
             SwapFunc>
  requires std::permutable<PermutationIterator>
constexpr void permute_inplace(PermutationIterator perm_iterator,
                               PermutationSentinel perm_sentinel,
                               SwapFunc swap_func) {
  // For implementation details see:
  // https://devblogs.microsoft.com/oldnewthing/20170102-00/?p=95095
  using Index = std::iter_value_t<PermutationIterator>;
  Index index{};
  for (; perm_iterator != perm_sentinel; ++perm_iterator, ++index) {
    Index current_index = index;
    PermutationIterator current_perm_iterator = perm_iterator;
    while (*current_perm_iterator != index) {
      STORM_ASSERT(*current_perm_iterator != Index{SIZE_MAX},
                   "Invalid permutation iterator value!");
      const Index new_index = *current_perm_iterator;
      swap_func(current_index, new_index);
      *current_perm_iterator = current_index, current_index = new_index;
      current_perm_iterator = perm_iterator + (current_index - index);
    }
    *current_perm_iterator = current_index;
  }
#if !STORM_NDEBUG
  std::ranges::fill(perm_iterator, perm_sentinel, Index{SIZE_MAX});
#endif
}
template<std::ranges::random_access_range PermutationRange,
         std::invocable<std::ranges::range_value_t<PermutationRange>,
                        std::ranges::range_value_t<PermutationRange>>
             SwapFunc>
  requires std::permutable<std::ranges::iterator_t<PermutationRange>>
constexpr void permute_inplace(PermutationRange&& perm_range,
                               SwapFunc swap_func) {
  permute_inplace(std::ranges::begin(perm_range), //
                  std::ranges::end(perm_range), std::move(swap_func));
}
/// @}

// -----------------------------------------------------------------------------

} // namespace Storm::permutations
