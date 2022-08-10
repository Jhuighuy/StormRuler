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

#pragma once

#include <array>
#include <concepts>
#include <functional>
#include <ostream>
#include <tuple>
#include <type_traits>
#include <utility>

#include <Storm/Base.hpp>

#include <Storm/Utils/Math.hpp>

#include <Storm/Blass/MatrixBase.hpp>

namespace Storm {

// clang-format off
template<class Derived>
  requires std::is_class_v<Derived> &&
           std::same_as<Derived, std::remove_cv_t<Derived>>
class MatrixViewInterface;
// clang-format on

namespace detail_ {
  // clang-format off
  template<class T, class U>
    requires(!std::same_as<T, MatrixViewInterface<U>>)
  void is_derived_from_matrix_view_interface_impl_(
      const T&, const MatrixViewInterface<U>&); // not defined
  template<class T>
  concept is_derived_from_matrix_view_interface_ =
      requires(T x) { is_derived_from_matrix_view_interface_impl_(x, x); };
  // clang-format on
} // namespace detail_

/// @brief Types, enabled to be a matrix view.
template<class T>
inline constexpr bool enable_matrix_view_v =
    detail_::is_derived_from_matrix_view_interface_<T>;

/// @brief Matrix view concept.
/// @todo In order to add the `movable` constraint, we
///   need to box the functor inside the `MapMatrixView`.
template<class MatrixView>
concept matrix_view =     //
    matrix<MatrixView> && // std::movable<MatrixView> &&
    enable_matrix_view_v<MatrixView>;

/// @brief Matrix that can be safely casted into a matrix view.
// clang-format off
template<class Matrix>
concept viewable_matrix = matrix<Matrix> &&
    ((matrix_view<std::remove_cvref_t<Matrix>> &&
      std::constructible_from<std::remove_cvref_t<Matrix>, Matrix>) ||
     (!matrix_view<std::remove_cvref_t<Matrix>> &&
      (std::is_lvalue_reference_v<Matrix> ||
       std::movable<std::remove_reference_t<Matrix>>)));
// clang-format on

/// @brief Base class for all matrix views.
// clang-format off
template<class Derived>
  requires std::is_class_v<Derived> &&
           std::same_as<Derived, std::remove_cv_t<Derived>>
class MatrixViewInterface {
  // clang-format on
private:

  constexpr Derived& self_() noexcept {
    static_assert(std::derived_from<Derived, MatrixViewInterface>);
    static_assert(matrix_view<Derived>);
    return static_cast<Derived&>(*this);
  }
  constexpr const Derived& self_() const noexcept {
    return static_cast<MatrixViewInterface&>(*this).self_();
  }

public:

  /// @brief Shape of the matrix.
  constexpr auto shape() const noexcept {
    return self_().shape();
  }

  /// @brief Get the vector coefficient at @p row_index.
  /// @{
  constexpr decltype(auto) operator()(size_t row_index) noexcept {
    return self_()(row_index, 0);
  }
  constexpr decltype(auto) operator()(size_t row_index) const noexcept {
    return self_()(row_index, 0);
  }
  /// @}

  /// @brief Get the matrix coefficient at @p row_index and @p col_index.
  /// @{
  constexpr decltype(auto) operator()(size_t row_index,
                                      size_t col_index) noexcept {
    return self_()(row_index, col_index);
  }
  constexpr decltype(auto) operator()(size_t row_index,
                                      size_t col_index) const noexcept {
    return self_()(row_index, col_index);
  }
  /// @}

}; // class MatrixViewInterface

/// @name Matrix views.
/// @{

/// @name Forwarding views.
/// @{

/// @brief Matrix reference view.
template<matrix Matrix>
class MatrixRefView : public MatrixViewInterface<MatrixRefView<Matrix>> {
private:

  Matrix* p_mat_;

  static void funс_(Matrix&); // not defined
  static void funс_(Matrix&&) = delete;

public:

  /// @brief Construct a matrix reference view.
  // clang-format off
  template<detail_::different_from_<MatrixRefView> OtherMatrix>
    requires std::convertible_to<OtherMatrix, Matrix&> &&
             requires { funс_(std::declval<OtherMatrix>()); }
  constexpr MatrixRefView(OtherMatrix&& mat) noexcept 
      : p_mat_{std::addressof(static_cast<Matrix&>(std::forward<OtherMatrix>(mat)))} {
    // clang-format on
  }

  /// @copydoc MatrixViewInterface::shape
  constexpr auto shape() const noexcept {
    return p_mat_->shape();
  }

  /// @copydoc MatrixViewInterface::operator()
  /// @{
  constexpr decltype(auto) operator()(size_t row_index,
                                      size_t col_index) noexcept {
    return (*p_mat_)(row_index, col_index);
  }
  constexpr decltype(auto) operator()(size_t row_index,
                                      size_t col_index) const noexcept {
    return std::as_const(*p_mat_)(row_index, col_index);
  }
  /// @}

}; // class MatrixRefView

/// @brief Matrix owning view.
// clang-format off
template<matrix Matrix>
  requires std::movable<Matrix>
class MatrixOwningView : public MatrixViewInterface<MatrixOwningView<Matrix>> {
  // clang-format on
private:

  STORM_NO_UNIQUE_ADDRESS_ Matrix mat_{};

public:

  /// @brief Construct an owning view.
  constexpr MatrixOwningView(Matrix&& mat) noexcept : mat_{std::move(mat)} {}

  /// @copydoc MatrixViewInterface::shape
  constexpr auto shape() const noexcept {
    return mat_.shape();
  }

  /// @copydoc MatrixViewInterface::operator()
  /// @{
  constexpr decltype(auto) operator()(size_t row_index,
                                      size_t col_index) noexcept {
    return mat_(row_index, col_index);
  }
  constexpr decltype(auto) operator()(size_t row_index,
                                      size_t col_index) const noexcept {
    return mat_(row_index, col_index);
  }
  /// @}

}; // class MatrixOwningView

template<class Matrix>
MatrixRefView(Matrix&) -> MatrixRefView<Matrix>;

namespace detail_ {
  // clang-format off
  template<class Matrix>
  concept matrix_can_ref_view_ = 
      requires { MatrixRefView{std::declval<Matrix>()}; };
  template<class Matrix>
  concept matrix_can_owning_view_ = 
      requires { MatrixOwningView{std::declval<Matrix>()}; };
  // clang-format on
} // namespace detail_

/// @brief Forward the viewable matrix @p mat as a matrix view.
// clang-format off
template<viewable_matrix Matrix>
  requires matrix_view<std::decay_t<Matrix>> || 
           detail_::matrix_can_ref_view_<Matrix> || 
           detail_::matrix_can_owning_view_<Matrix>
constexpr auto forward_as_matrix_view(Matrix&& mat) noexcept {
  // clang-format on
  if constexpr (matrix_view<std::decay_t<Matrix>>) {
    return std::forward<Matrix>(mat);
  } else if constexpr (detail_::matrix_can_ref_view_<Matrix>) {
    return MatrixRefView{std::forward<Matrix>(mat)};
  } else if constexpr (detail_::matrix_can_owning_view_<Matrix>) {
    return MatrixOwningView{std::forward<Matrix>(mat)};
  }
}

/// @brief Suitable matrix view type for a vieable matrix.
template<viewable_matrix Matrix>
using forward_as_matrix_view_t =
    decltype(forward_as_matrix_view(std::declval<Matrix>()));

/// @} // Forwarding views.

/// @name Generating views.
/// @{

/// @brief Matrix generating view.
/// @todo This is completely overengineered.
// clang-format off
template<class Shape, std::copy_constructible Func>
  requires is_matrix_shape_v<Shape> && std::is_object_v<Func> &&
           std::regular_invocable<Func, size_t, size_t>
class MakeMatrixView : public MatrixViewInterface<MakeMatrixView<Shape, Func>> {
  // clang-format on
private:

  STORM_NO_UNIQUE_ADDRESS_ Shape shape_;
  STORM_NO_UNIQUE_ADDRESS_ Func func_;

public:

  /// @brief Construct a generating view.
  constexpr MakeMatrixView(Shape shape, Func func) noexcept
      : shape_{shape}, func_{std::move(func)} {}

  /// @copydoc MatrixViewInterface::shape
  constexpr auto shape() const noexcept {
    return shape_;
  }

  /// @copydoc MatrixViewInterface::operator()
  constexpr decltype(auto) operator()(size_t row_index,
                                      size_t col_index) const noexcept {
    STORM_ASSERT_(row_index < shape_.num_rows() &&
                  col_index < shape_.num_cols() && "Indices are out of range.");
    return func_(row_index, col_index);
  }

}; // class MakeMatrixView

template<class Shape, class Func>
MakeMatrixView(Shape, Func) -> MakeMatrixView<Shape, Func>;

/// @brief Generate a constant matrix with @p num_rows and @p num_cols.
/// @param value Matrix element value.
template<std::copyable Value>
constexpr auto make_constant_matrix(matrix_extent auto num_rows,
                                    matrix_extent auto num_cols,
                                    Value value) noexcept {
  return MakeMatrixView( //
      MatrixShape(num_rows, num_cols),
      [value = std::move(value)](size_t, size_t) -> const Value& {
        return value;
      });
}

/// @brief Generate a fixed-sized constant matrix.
/// @tparam NumRows Number of the matrix rows.
/// @tparam NumCols Number of the matrix columns.
/// @param value Matrix element value.
template<size_t NumRows, size_t NumCols, std::copyable Value>
constexpr auto make_constant_matrix(Value value) {
  return make_constant_matrix<Value>( //
      size_t_constant<NumRows>{}, size_t_constant<NumCols>{}, std::move(value));
}

/// @brief Generate a constant vector with @p num_rows.
/// @param value Vector element value.
template<std::copyable Value>
constexpr auto make_constant_vector(matrix_extent auto num_rows,
                                    Value value) noexcept {
  return make_constant_matrix<Value>( //
      num_rows, size_t_constant<1>{}, std::move(value));
}

/// @brief Generate a fixed-sized constant vector.
/// @tparam NumRows Number of the vector rows.
/// @param value Vector element value.
template<size_t NumRows, std::copyable Value>
constexpr auto make_constant_matrix(Value value) noexcept {
  return make_constant_vector<Value>( //
      size_t_constant<NumRows>{}, std::move(value));
}

/// @brief Generate a diagonal matrix with @p num_rows, @p num_cols.
/// @param diagonal Matrix diagonal element value.
/// @param off_diagonal Matrix off-diagonal element value.
template<std::copyable Value>
constexpr auto make_diagonal_matrix(matrix_extent auto num_rows,
                                    matrix_extent auto num_cols, //
                                    Value diagonal,
                                    Value off_diagonal = Value{}) noexcept {
  // clang-format on
  return MakeMatrixView(
      MatrixShape(num_rows, num_cols),
      [diagonal = std::move(diagonal), off_diagonal = std::move(off_diagonal)](
          size_t row_index, size_t col_index) -> const Value& {
        return (row_index == col_index) ? diagonal : off_diagonal;
      });
}

/// @brief Generate a fixed-sized diagonal matrix.
/// @tparam NumRows Number of the matrix rows.
/// @tparam NumCols Number of the matrix columns.
/// @param diagonal Matrix diagonal element value.
/// @param off_diagonal Matrix off-diagonal element value.
template<size_t NumRows, size_t NumCols, std::movable Value>
constexpr auto make_diagonal_matrix(Value diagonal,
                                    Value off_diagonal = Value{}) noexcept {
  return make_diagonal_matrix<Value>(
      size_t_constant<NumRows>{}, size_t_constant<NumCols>{}, //
      std::move(diagonal), std::move(off_diagonal));
}

/// @} // Generating views.

#pragma region NOT THAT GOOD
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
    return size_t_constant<Size>{};
  }

  constexpr size_t operator[](size_t index) const noexcept {
    STORM_ASSERT_(index < Size && "Index is out of range.");
    return selected_indices_[index];
  }

}; // class SelectedIndices

template<size_t Size>
SelectedIndices(std::array<size_t, Size>) -> SelectedIndices<Size>;

namespace detail_ {
  constexpr static size_t range_size_(size_t from, //
                                      size_t to,   //
                                      size_t stride) noexcept {
    STORM_ASSERT_(stride != 0 && "Stride cannot be zero.");
    return (to - from) / stride;
  }
  template<size_t From, size_t To, size_t Stride>
  constexpr static auto range_size_(size_t_constant<From>, //
                                    size_t_constant<To>,
                                    size_t_constant<Stride>) noexcept {
    static_assert(Stride != 0 && "Stride cannot be zero.");
    return size_t_constant<(To - From) / Stride>{};
  }
} // namespace detail_

/// @brief Sliced indices range.
template<matrix_extent From, matrix_extent To, matrix_extent Stride>
class SlicedIndices {
private:

  STORM_NO_UNIQUE_ADDRESS_ From from_;
  STORM_NO_UNIQUE_ADDRESS_ To to_;
  STORM_NO_UNIQUE_ADDRESS_ Stride stride_;

public:

  constexpr SlicedIndices(From from, To to, Stride stride) noexcept
      : from_{from}, to_{to}, stride_{stride} {}

  constexpr auto size() const noexcept {
    return detail_::range_size_(from_, to_, stride_);
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
constexpr auto
slice_rows(viewable_matrix auto&& mat, //
           matrix_extent auto rows_from, matrix_extent auto rows_to,
           matrix_extent auto row_stride = size_t_constant<1>{}) noexcept {
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
constexpr auto slice_cols(viewable_matrix auto&& mat,
                          std::convertible_to<size_t> auto cols_from,
                          std::convertible_to<size_t> auto cols_to,
                          std::convertible_to<size_t> auto col_stride =
                              size_t_constant<1>{}) noexcept {
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
#pragma endregion

/// @} // Matrix views.

/// @brief Print a @p mat.
std::ostream& operator<<(std::ostream& out, matrix auto&& mat) {
  for (size_t row_index{0}; row_index < num_rows(mat); ++row_index) {
    out << "( ";
    for (size_t col_index{0}; col_index < num_cols(mat); ++col_index) {
      out << mat(row_index, col_index) << " ";
    }
    out << ")" << std::endl;
  }
  return out;
}

} // namespace Storm

//#include "MatrixViewFunctional.hpp"

#include "MatrixAction.hpp"
