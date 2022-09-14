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

#include <Storm/Utils/Index.hpp>
#include <Storm/Utils/IndexedContainers.hpp>
#include <Storm/Utils/Meta.hpp>

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

  IndexedVector<RowIndex, OffsetIndex> row_offsets_{OffsetIndex{0}};
  IndexedVector<OffsetIndex, ColIndex> col_indices_;

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
  template<std::ranges::input_range Range>
    requires std::same_as<std::ranges::range_value_t<Range>, ColIndex>
  constexpr void push_back(Range&& row_col_indices) {
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

  size_t row_size_hint_ = 25;
  IndexedVector<RowIndex, OffsetIndex> row_offsets_{OffsetIndex{0}};
  IndexedVector<RowIndex, OffsetIndex> row_end_offsets_;
  IndexedVector<OffsetIndex, ColIndex> col_indices_;

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
    col_indices_.insert(col_indices_.end(), row_size_hint_, ColIndex{SIZE_MAX});
    row_offsets_.emplace_back(col_indices_.size());
  }
  /// @brief Push back a row with @p row_col_indices.
  template<std::ranges::input_range Range>
    requires std::same_as<std::ranges::range_value_t<Range>, ColIndex>
  constexpr void push_back(Range&& row_col_indices) {
    col_indices_.insert(col_indices_.end(), //
                        row_col_indices.begin(), row_col_indices.end());
    row_end_offsets_.emplace_back(col_indices_.size());
    row_offsets_.emplace_back(col_indices_.size());
  }

  /// @brief Insert a connection at @p row_index, @p col_index.
  void connect(RowIndex row_index, ColIndex col_index) {
    STORM_ASSERT_(row_index < size(), "Row index is out of range!");
    if (row_end_offsets_[row_index] == row_offsets_[row_index + 1]) {
      const auto row_index_sz = static_cast<size_t>(row_index);
      const size_t row_size =
          row_end_offsets_[row_index] - row_offsets_[row_index];
      STORM_TRACE_("MCSR row {} capacity exhausted, row size is {}, {}!",
                   row_index_sz, row_size + 1,
                   meta::type_name_v<decltype(*this)>);
      if (row_size == row_size_hint_) { row_size_hint_ *= 2; }
      const size_t delta = row_size_hint_ - row_size;
      col_indices_.insert(col_indices_.begin() +
                              static_cast<size_t>(row_end_offsets_[row_index]),
                          delta, ColIndex{SIZE_MAX});
      std::for_each(row_offsets_.begin() + row_index_sz + 1, row_offsets_.end(),
                    [&](auto& offset) { offset += delta; });
      std::for_each(row_end_offsets_.begin() + row_index_sz + 1,
                    row_end_offsets_.end(),
                    [&](auto& offset) { offset += delta; });
    }
    col_indices_[row_end_offsets_[row_index]++] = col_index;
  }

}; // class McsrTable

/// @brief Modified compressed sparse row table without column values.
template<index RowIndex, index ColIndex>
using VoidMcsrTable = McsrTable<void, RowIndex, ColIndex>;

/// @brief Vector-of-vectors sparse row table.
template<class Value, index RowIndex, index ColIndex>
class VovTable final {
private:

  static_assert(std::is_same_v<Value, void>,
                "Non-void VovTable is not implemented yet!");

  IndexedVector<RowIndex, IndexedVector<OffsetIndex, ColIndex>> data_;

public:

  /// @brief Number of rows.
  [[nodiscard]] constexpr size_t size() const noexcept {
    return data_.size();
  }

  /// @brief Get the column indices range of a row @p row_index.
  /// @{
  [[nodiscard]] constexpr auto& operator[](RowIndex row_index) noexcept {
    STORM_ASSERT_(row_index < size(), "Row index is out of range!");
    return data_[row_index];
  }
  [[nodiscard]] constexpr const auto&
  operator[](RowIndex row_index) const noexcept {
    STORM_ASSERT_(row_index < size(), "Row index is out of range!");
    return data_[row_index];
  }
  /// @}

  /// @brief Push back an empty row.
  constexpr void push_back() {
    data_.emplace_back();
  }
  /// @brief Push back a row with @p row_col_indices.
  template<std::ranges::input_range Range>
    requires std::same_as<std::ranges::range_value_t<Range>, ColIndex>
  constexpr void push_back(Range&& row_col_indices) {
    data_.emplace_back();
    std::ranges::copy(row_col_indices, std::back_inserter(data_.back()));
  }

  /// @brief Insert a connection at @p row_index, @p col_index.
  void connect(RowIndex row_index, ColIndex col_index) {
    STORM_ASSERT_(row_index < size(), "Row index is out of range!");
    data_[row_index].push_back(col_index);
  }

}; // class VovTable

/// @brief Modified compressed sparse row table without column values.
template<index RowIndex, index ColIndex>
using VoidVovTable = VovTable<void, RowIndex, ColIndex>;

} // namespace Storm
