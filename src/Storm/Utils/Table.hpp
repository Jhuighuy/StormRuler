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

#include <ranges>

namespace Storm {

/// @brief Compressed sparse row table.
template<class Value, index RowIndex, index ColIndex>
class CsrTable final {
private:

  struct RowOffsetTag_;
  using RowOffsetIndex_ = Index<RowOffsetTag_>;

  IndexedVector<RowOffsetIndex_, RowIndex> row_offsets_{RowOffsetIndex_{0}};
  IndexedVector<ColIndex, RowOffsetIndex_> col_indices_;

public:

  [[nodiscard]] constexpr size_t size() const noexcept {
    return row_offsets_.size() - 1;
  }

  /// @{
  [[nodiscard]] constexpr auto operator[](RowIndex row_index) noexcept {
    STORM_ASSERT_(row_index < size(), "Row index is out of range!");
    return std::ranges::subrange( //
        col_indices_.begin() + (size_t) row_offsets_[row_index],
        col_indices_.begin() + (size_t) row_offsets_[row_index + 1]);
  }
  [[nodiscard]] constexpr auto operator[](RowIndex row_index) const noexcept {
    STORM_ASSERT_(row_index < size(), "Row index is out of range!");
    return std::ranges::subrange( //
        col_indices_.cbegin() + (size_t) row_offsets_[row_index],
        col_indices_.cbegin() + (size_t) row_offsets_[row_index + 1]);
  }
  /// @}

  void insert(RowIndex rowIndex, ColIndex columnIndex) {
    STORM_ASSERT_(rowIndex < size(), "");
    col_indices_.insert(col_indices_.begin() +
                            (size_t) row_offsets_[rowIndex + 1],
                        columnIndex);
    std::for_each(row_offsets_.begin() + (size_t) rowIndex + 1,
                  row_offsets_.end(),
                  [](RowOffsetIndex_& offset) { offset += 1; });
  }

  constexpr void emplace_back() {
    row_offsets_.emplace_back(col_indices_.size());
  }
  constexpr void emplace_back(auto&& range) {
    col_indices_.insert(col_indices_.end(), range.begin(), range.end());
    row_offsets_.emplace_back(col_indices_.size());
  }

}; // class CsrTable

template<index RowIndex, index ColIndex>
using VoidCsrTable = CsrTable<void, RowIndex, ColIndex>;

} // namespace Storm
