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

#include <array>
#include <compare>
#include <concepts>
#include <type_traits>
#include <vector>

namespace Storm {

/// @brief Strictly-tagged index.
/// @tparam Tag Any type to distinguish the different index types.
template<class Tag>
class Index final {
private:

  size_t value_;

public:

  /// @brief Construct an undefined index.
  constexpr Index() = default;

  /// @brief Construct an index with value.
  template<std::integral Integral>
  constexpr explicit Index(Integral value) noexcept
      : value_{static_cast<size_t>(value)} {}

  /// @brief Cast to the inegral value.
  template<std::integral Integral>
  constexpr explicit operator Integral() const noexcept {
    return static_cast<Integral>(value_);
  }

  /// @brief Cast into another tagged index.
  template<class OtherTag>
  constexpr explicit operator Index<OtherTag>() const noexcept {
    return Index<OtherTag>{value_};
  }

  /// @brief Comparison operators.
  /// @{
  constexpr auto operator<=>(const Index& other) const = default;
  template<std::integral Integral>
  constexpr auto operator<=>(Integral value) const noexcept {
    return value_ <=> static_cast<size_t>(value);
  }
  /// @}

  /// @brief Increment operators.
  /// @{
  constexpr Index& operator++() noexcept {
    return ++value_, *this;
  }
  constexpr Index operator++(int) noexcept {
    Index copy{*this};
    return ++value_, copy;
  }
  /// @}

  /// @brief Decrement operators.
  /// @{
  constexpr Index& operator--() noexcept {
    return --value_, *this;
  }
  constexpr Index operator--(int) noexcept {
    Index copy{*this};
    return --value_, copy;
  }
  /// @}

  /// @brief Addition operators.
  /// @{
  constexpr Index& operator+=(ptrdiff_t delta) noexcept {
    return value_ += delta, *this;
  }
  friend constexpr Index operator+(Index i, ptrdiff_t d) noexcept {
    return Index{i.value_ + d};
  }
  friend constexpr Index operator+(ptrdiff_t d, Index i) noexcept {
    return Index{d + i.value_};
  }
  /// @}

  /// @brief Subtraction operators.
  /// @{
  constexpr Index& operator-=(ptrdiff_t delta) noexcept {
    return value_ -= delta, *this;
  }
  friend constexpr Index operator-(Index i, ptrdiff_t d) noexcept {
    return Index{i.value_ - d};
  }
  /// @}

  /// @brief Difference operator.
  friend constexpr ptrdiff_t operator-(Index i1, Index i2) noexcept {
    return static_cast<ptrdiff_t>(i1.value_ - i2.value_);
  }

  /// @brief Read @p index from the input @p stream.
  friend std::istream& operator>>(std::istream& stream, Index& i) {
    return stream >> i.value_;
  }

  /// @brief Write @p index to the output @p stream.
  friend std::ostream& operator<<(std::ostream& stream, Index i) {
    return stream << i.value_;
  }

}; // class Index

// -----------------------------------------------------------------------------

template<class>
inline constexpr bool is_index_v = false;
template<class Tag>
inline constexpr bool is_index_v<Index<Tag>> = true;

/// @brief Index type.
template<class Index>
concept index = std::convertible_to<Index, size_t> || is_index_v<Index>;

} // namespace Storm
