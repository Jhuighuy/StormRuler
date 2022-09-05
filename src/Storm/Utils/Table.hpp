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

#include <Storm/Utils/Containers.hpp>

#include <algorithm>
#include <ranges>

namespace Storm {

struct OffsetTag;
using OffsetIndex = Index<OffsetTag>;

/// @brief Compressed sparse row table.
template<class Value, index RowIndex, index ColIndex>
class CsrTable final {
private:

  static_assert(std::is_same_v<Value, void>,
                "Non-void CsrTable is not implemented yet!");

  IndexedVector<OffsetIndex, RowIndex> row_offsets_{OffsetIndex{0}};
  IndexedVector<ColIndex, OffsetIndex> col_indices_;

public:

  /// @brief Number of rows.
  [[nodiscard]] constexpr size_t size() const noexcept {
    return row_offsets_.size() - 1;
  }

  /// @brief Get the column indices range of a row @p row_index.
  /// @{
  [[nodiscard]] constexpr auto operator[](RowIndex row_index) noexcept {
    STORM_ASSERT_(row_index < size(), "Row index is out of range!");
    const auto first = static_cast<size_t>(row_offsets_[row_index]);
    const auto last = static_cast<size_t>(row_offsets_[row_index + 1]);
    return std::ranges::subrange(col_indices_.begin() + first,
                                 col_indices_.begin() + last);
  }
  [[nodiscard]] constexpr auto operator[](RowIndex row_index) const noexcept {
    STORM_ASSERT_(row_index < size(), "Row index is out of range!");
    const auto first = static_cast<size_t>(row_offsets_[row_index]);
    const auto last = static_cast<size_t>(row_offsets_[row_index + 1]);
    return std::ranges::subrange(col_indices_.cbegin() + first,
                                 col_indices_.cbegin() + last);
  }
  /// @}

  /// @brief Push back an empty row.
  constexpr void push_back() {
    row_offsets_.emplace_back(col_indices_.size());
  }
  /// @brief Push back a row with @p row_col_indices.
  // clang-format off
  template<std::ranges::range Range>
    requires std::same_as<std::ranges::range_value_t<Range>, ColIndex>
  constexpr void push_back(Range&& row_col_indices) {
    // clang-format on
    col_indices_.insert(col_indices_.end(), //
                        row_col_indices.begin(), row_col_indices.end());
    row_offsets_.emplace_back(col_indices_.size());
  }

  /// @brief Insert a connection at @p row_index, @p col_index.
  void connect(RowIndex row_index, ColIndex col_index) {
    STORM_ASSERT_(row_index < size(), "Row index is out of range!");
    col_indices_.insert(col_indices_.begin() +
                            static_cast<size_t>(row_offsets_[row_index + 1]),
                        col_index);
    std::for_each(row_offsets_.begin() + static_cast<size_t>(row_index) + 1,
                  row_offsets_.end(), [](auto& offset) { offset += 1; });
  }

}; // class CsrTable

/// @brief Compressed sparse row table without column values.
template<index RowIndex, index ColIndex>
using VoidCsrTable = CsrTable<void, RowIndex, ColIndex>;

/// @brief Modified compressed sparse row table.
template<class Value, index RowIndex, index ColIndex>
class McsrTable final {
private:

  static_assert(std::is_same_v<Value, void>,
                "Non-void McsrTable is not implemented yet!");

  size_t row_capacity_ = 75;
  IndexedVector<OffsetIndex, RowIndex> row_offsets_{OffsetIndex{0}};
  IndexedVector<OffsetIndex, RowIndex> row_end_offsets_;
  IndexedVector<ColIndex, OffsetIndex> col_indices_;

public:

  /// @brief Number of rows.
  [[nodiscard]] constexpr size_t size() const noexcept {
    return row_offsets_.size() - 1;
  }

  /// @brief Get the column indices range of a row @p row_index.
  /// @{
  [[nodiscard]] constexpr auto operator[](RowIndex row_index) noexcept {
    STORM_ASSERT_(row_index < size(), "Row index is out of range!");
    const auto first = static_cast<size_t>(row_offsets_[row_index]);
    const auto last = static_cast<size_t>(row_end_offsets_[row_index]);
    return std::ranges::subrange(col_indices_.begin() + first,
                                 col_indices_.begin() + last);
  }
  [[nodiscard]] constexpr auto operator[](RowIndex row_index) const noexcept {
    STORM_ASSERT_(row_index < size(), "Row index is out of range!");
    const auto first = static_cast<size_t>(row_offsets_[row_index]);
    const auto last = static_cast<size_t>(row_end_offsets_[row_index]);
    return std::ranges::subrange(col_indices_.cbegin() + first,
                                 col_indices_.cbegin() + last);
  }
  /// @}

  /// @brief Push back an empty row.
  constexpr void push_back() {
    row_end_offsets_.emplace_back(col_indices_.size());
    col_indices_.insert(col_indices_.end(), row_capacity_, ColIndex{SIZE_MAX});
    row_offsets_.emplace_back(col_indices_.size());
  }
  /// @brief Push back a row with @p row_col_indices.
  // clang-format off
  template<std::ranges::range Range>
    requires std::same_as<std::ranges::range_value_t<Range>, ColIndex>
  constexpr void push_back(Range&& row_col_indices) {
    // clang-format on
    col_indices_.insert(col_indices_.end(), //
                        row_col_indices.begin(), row_col_indices.end());
    row_end_offsets_.emplace_back(col_indices_.size());
    row_offsets_.emplace_back(col_indices_.size());
  }

  /// @brief Insert a connection at @p row_index, @p col_index.
  void connect(RowIndex row_index, ColIndex col_index) {
    STORM_ASSERT_(row_index < size(), "Row index is out of range!");
    if (row_end_offsets_[row_index] < row_offsets_[row_index + 1]) {
      col_indices_[row_end_offsets_[row_index]++] = col_index;
    } else {
      STORM_WARNING_("MCSR row capacity exhausted, row size is {}, {}!",
                     row_end_offsets_[row_index] - row_offsets_[row_index],
                     meta::type_name_v<decltype(*this)>);
      col_indices_.insert(col_indices_.begin() +
                              static_cast<size_t>(row_offsets_[row_index + 1]),
                          col_index);
      std::for_each(row_offsets_.begin() + static_cast<size_t>(row_index) + 1,
                    row_offsets_.end(), [](auto& offset) { offset += 1; });
      std::for_each(row_end_offsets_.begin() + static_cast<size_t>(row_index),
                    row_end_offsets_.end(), [](auto& offset) { offset += 1; });
    }
  }

}; // class McsrTable

/// @brief Modified compressed sparse row table without column values.
template<index RowIndex, index ColIndex>
using VoidMcsrTable = McsrTable<void, RowIndex, ColIndex>;

} // namespace Storm
