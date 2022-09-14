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
/// FITNESS FOR Allocator PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT
/// SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
/// OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
/// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
/// DEALINGS IN THE SOFTWARE.

#include <Storm/Base.hpp>

#include <concepts>
#include <iterator>

namespace Storm {

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
namespace permutations {}

} // namespace Storm

namespace Storm::permutations {

/// @brief Invert the permutation. Complexity is linear.
/// @todo Add contrains!
/// @{
template<std::input_iterator PermutationIterator,
         std::sentinel_for<PermutationIterator> PermutationSentinel,
         class InversePermutationIterator>
constexpr void invert_permutation(PermutationIterator perm_iterator,
                                  PermutationSentinel perm_sentinel,
                                  InversePermutationIterator iperm_iterator) {
  // clang-format on
  std::iter_value_t<PermutationIterator> index{};
  for (; perm_iterator != perm_sentinel; ++perm_iterator, ++index) {
    STORM_ASSERT_(*perm_iterator != SIZE_MAX,
                  "Invalid permutation iterator value!");
    iperm_iterator[*perm_iterator] = index;
  }
}
template<std::ranges::input_range PermutationRange,
         class InversePermutationIterator>
constexpr void invert_permutation(PermutationRange&& perm_range,
                                  InversePermutationIterator iperm_iterator) {
  invert_permutation(std::ranges::begin(perm_range),
                     std::ranges::end(perm_range.end()), iperm_iterator);
}
/// @}

/// @brief Apply the permutation inplace. Complexity is linear.
/// @todo Add contrains!
/// @{
template<std::random_access_iterator PermutationIterator,
         std::sentinel_for<PermutationIterator> PermutationSentinel,
         std::invocable<std::iter_value_t<PermutationIterator>,
                        std::iter_value_t<PermutationIterator>>
             SwapFunc>
constexpr void permute(PermutationIterator perm_iterator,
                       PermutationSentinel perm_sentinel, SwapFunc swap_func) {
  std::iter_value_t<PermutationIterator> index{};
  for (; perm_iterator != perm_sentinel; ++perm_iterator, ++index) {
    auto current_index = index;
    auto current_perm_iterator = perm_iterator;
    while (*current_perm_iterator != index) {
      STORM_ASSERT_(*current_perm_iterator != SIZE_MAX,
                    "Invalid permutation iterator value!");
      const auto new_index = *current_perm_iterator;
      swap_func(current_index, new_index);
      *current_perm_iterator = current_index;
      current_index = new_index;
      current_perm_iterator = perm_iterator + (current_index - index);
    }
    *current_perm_iterator = current_index;
  }
#ifndef NDEBUG
  std::ranges::fill(perm_iterator, perm_sentinel,
                    std::iter_value_t<PermutationIterator>{SIZE_MAX});
#endif
}
/// @}

} // namespace Storm::permutations
