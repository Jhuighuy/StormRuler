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
/// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
/// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
/// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
/// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
/// IN THE SOFTWARE.

#ifndef STORM_INSIDE_MATRIX_VIEW_HPP_
#error Do not include this header directly, \
       use <Storm/Blass/MatrixView.hpp> instead.
#endif

#include <array>
#include <concepts>
#include <utility>

namespace Storm {

/// @name Slicing views.
/// @{

/// @brief All indices range.
class AllIndices {
public:

  constexpr auto size() const = delete;

  constexpr size_t operator[](size_t index) const noexcept {
    return index;
  }

}; // AllIndices

/// @brief Selected indices range.
template<size_t Size>
class SelectedIndices {
private:

  STORM_NO_UNIQUE_ADDRESS_ std::array<size_t, Size> selected_indices_;

public:

  constexpr explicit SelectedIndices(
      const std::array<size_t, Size>& selected_indices) noexcept
      : selected_indices_(selected_indices) {}

  constexpr auto size() const noexcept {
    return Size;
  }

  constexpr size_t operator[](size_t index) const noexcept {
    STORM_ASSERT_(index < Size && "Index is out of range.");
    return selected_indices_[index];
  }

}; // class SelectedIndices

template<size_t Size>
SelectedIndices(std::array<size_t, Size>) -> SelectedIndices<Size>;

/// @brief Sliced indices range.
class SlicedIndices {
private:

  size_t from_, to_, stride_;

public:

  constexpr SlicedIndices(size_t from, size_t to, size_t stride) noexcept
      : from_{from}, to_{to}, stride_{stride} {}

  constexpr auto size() const noexcept {
    return (from_ - to_) / stride_;
  }

  constexpr size_t operator[](size_t index) const noexcept {
    STORM_ASSERT_(index < size() && "Index is out of range.");
    return from_ + stride_ * index;
  }

}; // class SlicedIndices

/// @brief Submatrix view.
// clang-format off
template<matrix Matrix, class RowIndices, class ColIndices>
  requires std::is_object_v<RowIndices> && std::is_object_v<ColIndices>
class SubmatrixView :
    public MatrixViewInterface<SubmatrixView<Matrix, RowIndices, ColIndices>> {
  // clang-format on
private:

  STORM_NO_UNIQUE_ADDRESS_ Matrix mat_;
  STORM_NO_UNIQUE_ADDRESS_ RowIndices row_indices_;
  STORM_NO_UNIQUE_ADDRESS_ ColIndices col_indices_;

  // clang-format off
  constexpr auto num_rows_() const noexcept 
      requires requires { std::declval<RowIndices>().size(); } {
    // clang-format on
    return row_indices_.size();
  }
  constexpr auto num_rows_() const noexcept {
    return num_rows(mat_);
  }

  // clang-format off
  constexpr auto num_cols_() const noexcept 
      requires requires { std::declval<ColIndices>().size(); } {
    // clang-format on
    return col_indices_.size();
  }
  constexpr auto num_cols_() const noexcept {
    return num_cols(mat_);
  }

public:

  /// @brief Construct a matrix rows view.
  constexpr SubmatrixView(Matrix mat, //
                          RowIndices row_indices,
                          ColIndices col_indices) noexcept
      : mat_{std::move(mat)},                 //
        row_indices_{std::move(row_indices)}, //
        col_indices_{std::move(col_indices)} {}

  /// @copydoc MatrixViewInterface::shape
  constexpr matrix_shape_t shape() const noexcept {
    return {num_rows_(), num_cols_()};
  }

  /// @copydoc MatrixViewInterface::operator()
  /// @{
  constexpr decltype(auto) operator()(size_t row_index,
                                      size_t col_index) noexcept {
    return mat_(row_indices_[row_index], col_indices_[col_index]);
  }
  constexpr decltype(auto) operator()(size_t row_index,
                                      size_t col_index) const noexcept {
    return mat_(row_indices_[row_index], col_indices_[col_index]);
  }
  /// @}

}; // class SubmatrixView

template<class Matrix, class RowIndices, class ColIndices>
SubmatrixView(Matrix&&, RowIndices, ColIndices)
    -> SubmatrixView<forward_as_matrix_view_t<Matrix>, RowIndices, ColIndices>;

/// @brief Select the matrix @p mat rows with @p row_indices view.
constexpr auto
select_rows(viewable_matrix auto&& mat,
            std::convertible_to<size_t> auto... row_indices) noexcept {
  STORM_ASSERT_((static_cast<size_t>(row_indices) < num_rows(mat)) && ... &&
                "Row indices are out of range.");
  return SubmatrixView(
      STORM_FORWARD_(mat),
      SelectedIndices(std::array{static_cast<size_t>(row_indices)...}),
      AllIndices{});
}

/// @brief Select the matrix @p mat columns with @p col_index view.
constexpr auto
select_cols(viewable_matrix auto&& mat,
            std::convertible_to<size_t> auto... col_indices) noexcept {
  STORM_ASSERT_((static_cast<size_t>(col_indices) < mat.num_cols()) && ... &&
                "Columns indices are out of range.");
  return SubmatrixView(
      STORM_FORWARD_(mat), //
      AllIndices{},
      SelectedIndices(std::array{static_cast<size_t>(col_indices)...}));
}

/// @brief Slice the matrix @p mat rows from index @p rows_from
///   to index @p rows_to (not including) with a stride @p row_stride view.
/// @{
constexpr auto slice_rows(viewable_matrix auto&& mat, //
                          size_t rows_from, size_t rows_to,
                          size_t row_stride = 1) noexcept {
  STORM_ASSERT_(rows_from < rows_to &&
                static_cast<size_t>(rows_to) <= num_rows(mat) &&
                "Invalid row range.");
  return SubmatrixView(STORM_FORWARD_(mat),
                       SlicedIndices(rows_from, rows_to, row_stride),
                       AllIndices{});
}
template<size_t RowsFrom, size_t RowsTo, size_t RowStride = 1>
constexpr auto slice_rows(viewable_matrix auto&& mat) {
  return slice_rows(STORM_FORWARD_(mat), size_t_constant<RowsFrom>{},
                    size_t_constant<RowsTo>{}, size_t_constant<RowStride>{});
}
/// @}

/// @brief Slice the matrix @p mat columns from index @p cols_from
///   to index @p cols_to (not including) with a stride @p col_stride view.
/// @{
constexpr auto slice_cols(viewable_matrix auto&& mat, //
                          size_t cols_from, size_t cols_to,
                          size_t col_stride = 1) noexcept {
  STORM_ASSERT_(cols_from < cols_to &&
                static_cast<size_t>(cols_to) <= num_cols(mat) &&
                "Invalid column range.");
  return SubmatrixView(STORM_FORWARD_(mat), //
                       AllIndices{},
                       SlicedIndices(cols_from, cols_to, col_stride));
}
template<size_t ColsFrom, size_t ColsTo, size_t ColStride = 1>
constexpr auto slice_cols(viewable_matrix auto&& mat) {
  return slice_cols(STORM_FORWARD_(mat), size_t_constant<ColsFrom>{},
                    size_t_constant<ColsTo>{}, size_t_constant<ColStride>{});
}
/// @}

/// @todo Implement me!
constexpr auto diag(viewable_matrix auto&& mat) noexcept;

/// @todo Implement me!
constexpr auto lower_triangle(viewable_matrix auto&& mat) noexcept;

/// @todo Implement me!
constexpr auto upper_triangle(viewable_matrix auto&& mat) noexcept;

/// @} // Slicing views.

} // namespace Storm
