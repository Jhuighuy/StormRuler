/// Copyright (C) 2020 - 2023 Oleg Butakov
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
    STORM_ASSERT_(index < Size, "Index is out of range.");
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
    STORM_ASSERT_(index < size(), "Index is out of range.");
    return from_ + stride_ * index;
  }

}; // class SlicedIndices

/// @brief Submatrix view.
template<matrix_view Matrix, class RowIndices, class ColIndices>
  requires std::is_object_v<RowIndices> && std::is_object_v<ColIndices>
class SubmatrixView :
    public MatrixViewInterface<SubmatrixView<Matrix, RowIndices, ColIndices>> {
private:

  STORM_NO_UNIQUE_ADDRESS_ Matrix mat_;
  STORM_NO_UNIQUE_ADDRESS_ RowIndices row_indices_;
  STORM_NO_UNIQUE_ADDRESS_ ColIndices col_indices_;

  // clang-format off
  constexpr auto num_rows_() const noexcept 
      requires requires { row_indices_.size(); } {
    // clang-format on
    return row_indices_.size();
  }
  constexpr auto num_rows_() const noexcept {
    return num_rows(mat_);
  }

  constexpr auto num_cols_() const noexcept
    requires requires { col_indices_.size(); }
  {
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
  constexpr auto shape() const noexcept {
    return MatrixShape{num_rows_(), num_cols_()};
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
template<viewable_matrix Matrix, //
         std::convertible_to<size_t>... RowIndices>
[[nodiscard]] constexpr auto select_rows(Matrix&& mat,
                                         RowIndices&&... row_indices) {
  return SubmatrixView(std::forward<Matrix>(mat),
                       SelectedIndices(std::array{static_cast<size_t>(
                           std::forward<RowIndices>(row_indices))...}),
                       AllIndices{});
}

/// @brief Select the matrix @p mat columns with @p col_index view.
template<viewable_matrix Matrix, //
         std::convertible_to<size_t>... ColIndices>
[[nodiscard]] constexpr auto select_cols(Matrix&& mat,
                                         ColIndices&&... col_indices) {
  return SubmatrixView(std::forward<Matrix>(mat), //
                       AllIndices{},
                       SelectedIndices(std::array{static_cast<size_t>(
                           std::forward<ColIndices>(col_indices))...}));
}

/// @brief Slice the matrix @p mat rows from index @p rows_from
///   to index @p rows_to (not including) with a stride @p row_stride view.
template<viewable_matrix Matrix>
[[nodiscard]] constexpr auto slice_rows(Matrix&& mat, //
                                        size_t rows_from, size_t rows_to,
                                        size_t row_stride = 1) {
  return SubmatrixView(std::forward<Matrix>(mat),
                       SlicedIndices(rows_from, rows_to, row_stride),
                       AllIndices{});
}

/// @brief Slice the matrix @p mat columns from index @p cols_from
///   to index @p cols_to (not including) with a stride @p col_stride view.
template<viewable_matrix Matrix>
[[nodiscard]] constexpr auto slice_cols(Matrix&& mat, //
                                        size_t cols_from, size_t cols_to,
                                        size_t col_stride = 1) {
  return SubmatrixView(std::forward<Matrix>(mat), //
                       AllIndices{},
                       SlicedIndices(cols_from, cols_to, col_stride));
}

/// @todo Implement me!
constexpr auto diag(viewable_matrix auto&& mat) noexcept;

/// @todo Implement me!
constexpr auto lower_triangle(viewable_matrix auto&& mat) noexcept;

/// @todo Implement me!
constexpr auto upper_triangle(viewable_matrix auto&& mat) noexcept;

/// @} // Slicing views.

} // namespace Storm
