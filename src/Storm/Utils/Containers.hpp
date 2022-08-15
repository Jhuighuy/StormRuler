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

#pragma once

#include <Storm/Base.hpp>

#include <array>
#include <compare>
#include <concepts>
#include <type_traits>
#include <vector>

namespace Storm {

/// @brief Strictly-tagged index.
template<class Tag>
class Index final {
private:

  size_t value_;

public:

  /// @brief Construct an undefined index.
  constexpr Index() = default;

  /// @brief Construct an index with value.
  constexpr explicit Index(size_t value) noexcept : value_{value} {}

  /// @brief Cast to the underlying value.
  [[nodiscard]] constexpr explicit operator size_t() const noexcept {
    return value_;
  }

  /// @brief Cast into another tagged index.
  template<class OtherTag>
  [[nodiscard]] constexpr explicit operator Index<OtherTag>() const noexcept {
    return Index<OtherTag>{value_};
  }

  /// @brief Comparison operator.
  /// @{
  [[nodiscard]] constexpr auto operator<=>(const Index& other) const = default;
  [[nodiscard]] constexpr auto operator<=>(size_t value) const noexcept {
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
  [[nodiscard]] friend constexpr Index //
  operator+(Index i, ptrdiff_t d) noexcept {
    return Index{i.value_ + d};
  }
  [[nodiscard]] friend constexpr Index //
  operator+(ptrdiff_t d, Index i) noexcept {
    return Index{d + i.value_};
  }
  /// @}

  /// @brief Subtraction operators.
  /// @{
  constexpr Index& operator-=(ptrdiff_t delta) noexcept {
    return value_ -= delta, *this;
  }
  [[nodiscard]] friend constexpr Index //
  operator-(Index i, ptrdiff_t d) noexcept {
    return Index{i.value_ - d};
  }
  /// @}

  /// @brief Difference operator.
  [[nodiscard]] friend constexpr ptrdiff_t //
  operator-(Index i1, Index i2) noexcept {
    return static_cast<ptrdiff_t>(i1.value_ - i2.value_);
  }

}; // class Index

namespace detail_ {
  template<class>
  inline constexpr bool is_index_v_ = false;
  template<class Tag>
  inline constexpr bool is_index_v_<Index<Tag>> = true;
} // namespace detail_

/// @brief Index type.
template<class Index>
concept index = detail_::is_index_v_<Index>;

/// @brief Wrapper for std::array with strict indexing.
// clang-format off
template<class Value, size_t Size, index Index>
  requires std::is_object_v<Value>
class IndexedArray final : public std::array<Value, Size> {
  // clang-format on
public:

  /// @brief Element at @p index.
  /// @{
  constexpr Value& at(size_t) = delete;
  constexpr const Value& at(size_t) const = delete;
  [[nodiscard]] constexpr Value& at(Index index) {
    return std::array<Value, Size>::at(static_cast<size_t>(index));
  }
  [[nodiscard]] constexpr const Value& at(Index index) const {
    return std::array<Value, Size>::at(static_cast<size_t>(index));
  }
  /// @}

  /// @brief Element at @p index.
  /// @{
  constexpr Value& operator[](size_t) = delete;
  constexpr const Value& operator[](size_t) const = delete;
  [[nodiscard]] constexpr Value& operator[](Index index) noexcept {
    STORM_ASSERT_(index < Size, "Index is out of range!");
    return std::array<Value, Size>::operator[](static_cast<size_t>(index));
  }
  [[nodiscard]] constexpr const Value& operator[](Index index) const noexcept {
    STORM_ASSERT_(index < Size, "Index is out of range!");
    return std::array<Value, Size>::operator[](static_cast<size_t>(index));
  }
  /// @}

}; // class IndexedArray

/// @brief Wrapper for std::vector with strict indexing.
// clang-format off
template<class Value, index Index, class Allocator = std::allocator<Value>>
  requires std::is_object_v<Value> && std::is_object_v<Allocator>
class IndexedVector final : public std::vector<Value, Allocator> {
  // clang-format on
public:

  /// @brief Construct an empty indexed vector.
  constexpr IndexedVector() = default;
  /// @brief Construct an empty indexed vector with a given @p allocator.
  constexpr explicit IndexedVector(const Allocator& allocator) noexcept
      : std::vector<Value, Allocator>(allocator) {}

  /// @brief Construct an indexed vector with @p count copies of @p value.
  constexpr IndexedVector(size_t count, const Value& value,
                          const Allocator& allocator = Allocator{})
      : std::vector<Value, Allocator>(count, value, allocator) {}
  /// @brief Construct an indexed array with @p count default-inserted values.
  constexpr explicit IndexedVector(size_t count,
                                   const Allocator& allocator = Allocator{})
      : std::vector<Value, Allocator>(count, allocator) {}

  /// @brief Construst an indexex vector with the iterator range.
  /// @param allocator Allocator instance.
  template<class Iterator>
  constexpr IndexedVector(Iterator first, Iterator last,
                          const Allocator& allocator = Allocator{})
      : std::vector<Value, Allocator>(first, last, allocator) {}

  /// @brief Copy-construct an indexed vector.
  constexpr IndexedVector(const IndexedVector& other)
      : std::vector<Value, Allocator>(other) {}
  /// @brief Copy-construct an indexed vector.
  /// @param allocator Allocator instance.
  constexpr IndexedVector(const IndexedVector& other,
                          const Allocator& allocator)
      : std::vector<Value, Allocator>(other, allocator) {}

  /// @brief Move-construct an indexed vector.
  constexpr IndexedVector(IndexedVector&& other) noexcept
      : std::vector<Value, Allocator>(
            static_cast<std::vector<Value, Allocator>&&>(other)) {}
  /// @brief Move-construct an indexed vector.
  /// @param allocator Allocator instance.
  constexpr IndexedVector(IndexedVector&& other, const Allocator& allocator)
      : std::vector<Value, Allocator>(
            static_cast<std::vector<Value, Allocator>&&>(other), allocator) {}

  /// @brief Construct an indexed vector with an initializer list @p list.
  /// @param allocator Allocator instance.
  constexpr IndexedVector(std::initializer_list<Value> list,
                          const Allocator& allocator = Allocator{})
      : std::vector<Value, Allocator>(list, allocator) {}

  /// @brief Copy-assignment operator.
  constexpr IndexedVector& operator=(const IndexedVector& other) {
    std::vector<Value, Allocator>::operator=(other);
    return *this;
  }

  /// @brief Move-assignment operator.
  constexpr IndexedVector& operator=(IndexedVector&& other) noexcept(
      noexcept(std::allocator_traits<
                   Allocator>::propagate_on_container_move_assignment::value ||
               std::allocator_traits<Allocator>::is_always_equal::value)) {
    std::vector<Value, Allocator>::operator=(
        static_cast<std::vector<Value, Allocator>&&>(other));
    return *this;
  }

  /// @brief Assign an initializer list @p list.
  constexpr IndexedVector& operator=(std::initializer_list<Value> list) {
    return std::vector<Value, Allocator>::operator=(list), *this;
  }

  /// @todo Document me!
  constexpr void assign(const IndexedVector& other, bool copy = true) {
    if (copy) {
      this->operator=(other);
    } else {
      this->resize(other.size());
    }
  }

  /// @brief Swap the indexed array contents.
  constexpr void swap(IndexedVector& other) noexcept( //
      noexcept(std::allocator_traits<
                   Allocator>::propagate_on_container_swap::value ||
               std::allocator_traits<Allocator>::is_always_equal::value)) {
    std::vector<Value, Allocator>::swap(other);
  }

  /// @brief Element at @p index.
  /// @{
  constexpr Value& at(size_t) = delete;
  constexpr const Value& at(size_t) const = delete;
  [[nodiscard]] constexpr Value& at(Index index) {
    return std::vector<Value, Allocator>::at(static_cast<size_t>(index));
  }
  [[nodiscard]] constexpr const Value& at(Index index) const {
    return std::vector<Value, Allocator>::at(static_cast<size_t>(index));
  }
  /// @}

  /// @brief Element at @p index.
  /// @{
  constexpr Value& operator[](size_t) = delete;
  constexpr const Value& operator[](size_t) const = delete;
  [[nodiscard]] constexpr Value& operator[](Index index) noexcept {
    STORM_ASSERT_(index < this->size(), "Index is out of range!");
    return std::vector<Value, Allocator>::operator[](
        static_cast<size_t>(index));
  }
  [[nodiscard]] constexpr const Value& operator[](Index index) const noexcept {
    STORM_ASSERT_(index < this->size(), "Index is out of range!");
    return std::vector<Value, Allocator>::operator[](
        static_cast<size_t>(index));
  }
  /// @}

}; // class IndexedVector

} // namespace Storm
