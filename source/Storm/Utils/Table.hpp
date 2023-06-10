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

#include <Storm/Utils/Crtp.hpp>
#include <Storm/Utils/Index.hpp>
#include <Storm/Utils/IndexedContainers.hpp>
#include <Storm/Utils/Meta.hpp>

#include <algorithm>
#include <ranges>
#include <span>
#include <utility>

namespace Storm {

// -----------------------------------------------------------------------------

template<crtp_derived Derived>
class TableInterface;

/// @brief Types, enabled to be a table.
template<class Table>
inline constexpr bool enable_table_v =
    derived_from_crtp_interface<Table, TableInterface>;

/// @brief Table concept.
template<class Table>
concept table = std::movable<Table> && enable_table_v<Table>;

struct OffsetTag;
using OffsetIndex = Index<OffsetTag>;

/// @brief Copy the @p in_table to @p out_table.
template<table InTable, table OutTable>
void copy_table(const InTable& in_table, OutTable& out_table) {
  if constexpr (std::assignable_from<OutTable, InTable>) {
    out_table = in_table;
  } else {
    out_table.clear();
    out_table.reserve(in_table.size());
    for (auto row_index : in_table.rows()) {
      out_table.push_back(in_table[row_index]);
    }
  }
}
/// @brief Move the @p in_table to @p out_table.
template<table InTable, table OutTable>
void move_table(InTable&& in_table, OutTable& out_table) {
  if constexpr (std::assignable_from<OutTable, InTable>) {
    out_table = std::move(in_table);
  } else {
    out_table.clear();
    out_table.reserve(in_table.size());
    for (auto row_index : in_table.rows()) {
      out_table.push_back(std::move(in_table[row_index]));
    }
    in_table = {};
  }
}

/// @brief CRTP interface to a table.
template<crtp_derived Derived>
class TableInterface {
private:

  constexpr Derived& self_() noexcept {
    static_assert(std::derived_from<Derived, TableInterface>);
    return static_cast<Derived&>(*this);
  }
  constexpr const Derived& self_() const noexcept {
    return const_cast<TableInterface&>(*this).self_();
  }

public:

  /// @brief Number of the table rows.
  constexpr size_t size() const noexcept {
    return std::ranges::size(self_().rows());
  }

  /// @brief Assign the @p table.
  /// @todo Refactor this in some point of time.
  /// @{
  void assign(Derived&& table) {
    self_() = std::move(table);
  }
  void assign(const Derived& table) {
    self_() = table;
  }
  template<table Table>
  void assign(const Table& table) {
    self_().clear();
    self_().reserve(table.size());
    for (auto row_index : table.rows()) {
      self_().push_back(table[row_index]);
    }
  }
  /// @}

}; // TableInterface

// -----------------------------------------------------------------------------

/// @brief Compressed sparse row table.
template<index RowIndex, class ColValue>
  requires std::is_object_v<ColValue>
class CsrTable final : public TableInterface<CsrTable<RowIndex, ColValue>> {
private:

  IndexedVector<RowIndex, OffsetIndex> row_offsets_{OffsetIndex{0}};
  IndexedVector<OffsetIndex, ColValue> col_values_;

public:

  /// @brief Index range of the rows.
  constexpr auto rows() const noexcept {
    return std::views::iota(RowIndex{0}, RowIndex{row_offsets_.size() - 1});
  }

  /// @brief Get the column indices range of a row @p row_index.
  /// @{
  constexpr auto operator[](RowIndex row_index) noexcept {
    STORM_ASSERT_(row_index < this->size(), "Row index is out of range!");
    // Here we do `{col_values_[begin - 1], &col_values_[end - 1] + 1}`,
    // since `col_values_[end - 1]` may be outsize of the vector.
    return std::span{&col_values_[row_offsets_[row_index]],
                     &col_values_[row_offsets_[row_index + 1] - 1] + 1};
  }
  constexpr auto operator[](RowIndex row_index) const noexcept {
    STORM_ASSERT_(row_index < this->size(), "Row index is out of range!");
    return std::span{&col_values_[row_offsets_[row_index]],
                     &col_values_[row_offsets_[row_index + 1] - 1] + 1};
  }
  /// @}

  /// @brief Clear the table.
  constexpr void clear() {
    row_offsets_.clear(), row_offsets_.emplace_back(0);
    col_values_.clear();
  }

  /// @brief Reserve the row storage.
  constexpr void reserve(size_t capacity) {
    row_offsets_.reserve(capacity + 1);
  }

  /// @brief Push back an empty row.
  constexpr void push_back() {
    row_offsets_.emplace_back(col_values_.size());
  }
  /// @brief Push back a row with @p row_values.
  template<std::ranges::input_range Range>
    requires std::same_as<std::ranges::range_value_t<Range>, ColValue>
  constexpr void push_back(Range&& row_values) {
    col_values_.insert(col_values_.end(), row_values.begin(), row_values.end());
    row_offsets_.emplace_back(col_values_.size());
  }

  /// @brief Insert a connection at @p row_index, @p col_value.
  /// @warning Complexity is linear!
  void insert(RowIndex row_index, ColValue col_value) {
    STORM_ASSERT_(row_index < this->size(), "Row index is out of range!");
    const auto offset = static_cast<size_t>(row_offsets_[row_index + 1]);
    col_values_.insert(col_values_.begin() + offset, col_value);
    std::for_each(row_offsets_.begin() + static_cast<size_t>(row_index) + 1,
                  row_offsets_.end(), [](auto& offset) { offset += 1; });
  }

}; // class CsrTable

// -----------------------------------------------------------------------------

/// @brief Vector-of-vectors sparse row table.
template<index RowIndex, class ColValue = void>
  requires std::is_object_v<ColValue>
class VovTable final : public TableInterface<VovTable<RowIndex, ColValue>> {
private:

  IndexedVector<RowIndex, IndexedVector<OffsetIndex, ColValue>> rows_;

public:

  /// @brief Index range of the rows.
  constexpr auto rows() const noexcept {
    return std::views::iota(RowIndex{0}, RowIndex{rows_.size()});
  }

  /// @brief Get the column indices range of a row @p row_index.
  /// @{
  constexpr auto& operator[](RowIndex row_index) noexcept {
    STORM_ASSERT_(row_index < this->size(), "Row index is out of range!");
    return rows_[row_index];
  }
  constexpr const auto& operator[](RowIndex row_index) const noexcept {
    STORM_ASSERT_(row_index < this->size(), "Row index is out of range!");
    return rows_[row_index];
  }
  /// @}

  /// @brief Clear the table.
  constexpr void clear() {
    rows_.clear();
  }

  /// @brief Reserve the row storage.
  constexpr void reserve(size_t capacity) {
    rows_.reserve(capacity);
  }

  /// @brief Push back an empty row.
  constexpr void push_back() {
    rows_.emplace_back();
  }
  /// @brief Push back a row with @p row_col_indices.
  /// @{
  template<std::ranges::input_range Range>
    requires std::same_as<std::ranges::range_value_t<Range>, ColValue>
  constexpr void push_back(Range&& row_values) {
    rows_.emplace_back();
    std::ranges::copy(row_values, std::back_inserter(rows_.back()));
  }
  constexpr void push_back(IndexedVector<OffsetIndex, ColValue>&& row_values) {
    rows_.emplace_back(std::move(row_values));
  }
  /// @}

  /// @brief Insert a connection at @p row_index, @p col_value.
  void insert(RowIndex row_index, ColValue col_value) {
    STORM_ASSERT_(row_index < this->size(), "Row index is out of range!");
    rows_[row_index].push_back(col_value);
  }

}; // class VovTable

// -----------------------------------------------------------------------------

} // namespace Storm
