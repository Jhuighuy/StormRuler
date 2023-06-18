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

  constexpr Derived& _self() noexcept {
    static_assert(std::derived_from<Derived, TableInterface>);
    return static_cast<Derived&>(*this);
  }
  constexpr const Derived& _self() const noexcept {
    return const_cast<TableInterface&>(*this)._self();
  }

public:

  /// @brief Number of the table rows.
  constexpr size_t size() const noexcept {
    return std::ranges::size(_self().rows());
  }

  /// @brief Assign the @p table.
  /// @todo Refactor this in some point of time.
  /// @{
  void assign(Derived&& table) {
    _self() = std::move(table);
  }
  void assign(const Derived& table) {
    _self() = table;
  }
  template<table Table>
  void assign(const Table& table) {
    _self().clear();
    _self().reserve(table.size());
    for (auto row_index : table.rows()) {
      _self().push_back(table[row_index]);
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

  IndexedVector<RowIndex, OffsetIndex> _row_offsets{OffsetIndex{0}};
  IndexedVector<OffsetIndex, ColValue> _col_values;

public:

  /// @brief Index range of the rows.
  constexpr auto rows() const noexcept {
    return std::views::iota(RowIndex{0}, RowIndex{_row_offsets.size() - 1});
  }

  /// @brief Get the column indices range of a row @p row_index.
  /// @{
  constexpr auto operator[](RowIndex row_index) noexcept {
    STORM_ASSERT(row_index < this->size(), "Row index is out of range!");
    // Here we do `{_col_values[begin - 1], &_col_values[end - 1] + 1}`,
    // since `_col_values[end - 1]` may be outsize of the vector.
    return std::span{&_col_values[_row_offsets[row_index]],
                     &_col_values[_row_offsets[row_index + 1] - 1] + 1};
  }
  constexpr auto operator[](RowIndex row_index) const noexcept {
    STORM_ASSERT(row_index < this->size(), "Row index is out of range!");
    return std::span{&_col_values[_row_offsets[row_index]],
                     &_col_values[_row_offsets[row_index + 1] - 1] + 1};
  }
  /// @}

  /// @brief Clear the table.
  constexpr void clear() {
    _row_offsets.clear(), _row_offsets.emplace_back(0);
    _col_values.clear();
  }

  /// @brief Reserve the row storage.
  constexpr void reserve(size_t capacity) {
    _row_offsets.reserve(capacity + 1);
  }

  /// @brief Push back an empty row.
  constexpr void push_back() {
    _row_offsets.emplace_back(_col_values.size());
  }
  /// @brief Push back a row with @p row_values.
  template<std::ranges::input_range Range>
    requires std::same_as<std::ranges::range_value_t<Range>, ColValue>
  constexpr void push_back(Range&& row_values) {
    _col_values.insert(_col_values.end(), row_values.begin(), row_values.end());
    _row_offsets.emplace_back(_col_values.size());
  }

  /// @brief Insert a connection at @p row_index, @p col_value.
  /// @warning Complexity is linear!
  void insert(RowIndex row_index, ColValue col_value) {
    STORM_ASSERT(row_index < this->size(), "Row index is out of range!");
    const auto offset = static_cast<size_t>(_row_offsets[row_index + 1]);
    _col_values.insert(_col_values.begin() + offset, col_value);
    std::for_each(_row_offsets.begin() + static_cast<size_t>(row_index) + 1,
                  _row_offsets.end(), [](auto& offset) { offset += 1; });
  }

}; // class CsrTable

// -----------------------------------------------------------------------------

/// @brief Vector-of-vectors sparse row table.
template<index RowIndex, class ColValue = void>
  requires std::is_object_v<ColValue>
class VovTable final : public TableInterface<VovTable<RowIndex, ColValue>> {
private:

  IndexedVector<RowIndex, IndexedVector<OffsetIndex, ColValue>> _rows;

public:

  /// @brief Index range of the rows.
  constexpr auto rows() const noexcept {
    return std::views::iota(RowIndex{0}, RowIndex{_rows.size()});
  }

  /// @brief Get the column indices range of a row @p row_index.
  /// @{
  constexpr auto& operator[](RowIndex row_index) noexcept {
    STORM_ASSERT(row_index < this->size(), "Row index is out of range!");
    return _rows[row_index];
  }
  constexpr const auto& operator[](RowIndex row_index) const noexcept {
    STORM_ASSERT(row_index < this->size(), "Row index is out of range!");
    return _rows[row_index];
  }
  /// @}

  /// @brief Clear the table.
  constexpr void clear() {
    _rows.clear();
  }

  /// @brief Reserve the row storage.
  constexpr void reserve(size_t capacity) {
    _rows.reserve(capacity);
  }

  /// @brief Push back an empty row.
  constexpr void push_back() {
    _rows.emplace_back();
  }
  /// @brief Push back a row with @p row_col_indices.
  /// @{
  template<std::ranges::input_range Range>
    requires std::same_as<std::ranges::range_value_t<Range>, ColValue>
  constexpr void push_back(Range&& row_values) {
    _rows.emplace_back();
    std::ranges::copy(row_values, std::back_inserter(_rows.back()));
  }
  constexpr void push_back(IndexedVector<OffsetIndex, ColValue>&& row_values) {
    _rows.emplace_back(std::move(row_values));
  }
  /// @}

  /// @brief Insert a connection at @p row_index, @p col_value.
  void insert(RowIndex row_index, ColValue col_value) {
    STORM_ASSERT(row_index < this->size(), "Row index is out of range!");
    _rows[row_index].push_back(col_value);
  }

}; // class VovTable

// -----------------------------------------------------------------------------

} // namespace Storm
