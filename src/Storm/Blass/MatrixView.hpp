/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// Copyright (C) 2022 Oleg Butakov
///
/// Permission is hereby granted, free of charge, to any person
/// obtaining a copy of this software and associated documentation
/// files (the "Software"), to deal in the Software without
/// restriction, including without limitation the rights  to use,
/// copy, modify, merge, publish, distribute, sublicense, and/or
/// sell copies of the Software, and to permit persons to whom the
/// Software is furnished to do so, subject to the following
/// conditions:
///
/// The above copyright notice and this permission notice shall be
/// included in all copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
/// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
/// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
/// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
/// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
/// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
/// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
/// OTHER DEALINGS IN THE SOFTWARE.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///

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

namespace Detail_ {
  // clang-format off
  template<class T, class U>
    requires (!std::same_as<T, MatrixViewInterface<U>>)
  void is_derived_from_matrix_view_interface_impl_(
      const T&, const MatrixViewInterface<U>&); // not defined
  template<class T>
  concept is_derived_from_matrix_view_interface_ =
      requires(T x) { is_derived_from_matrix_view_interface_impl_(x, x); };
  // clang-format on
} // namespace Detail_

/// @brief Types, enabled to be a matrix view.
template<class T>
inline constexpr bool enable_matrix_view_v =
    Detail_::is_derived_from_matrix_view_interface_<T>;

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

/// ----------------------------------------------------------------- ///
/// @brief Base class for all matrix views.
/// ----------------------------------------------------------------- ///
// clang-format off
template<class Derived>
  requires std::is_class_v<Derived> &&
           std::same_as<Derived, std::remove_cv_t<Derived>>
class MatrixViewInterface {
  // clang-format on
private:

  auto& self_() noexcept {
    static_assert(std::derived_from<Derived, MatrixViewInterface<Derived>>);
    // static_assert(matrix_view<Derived>);
    return static_cast<Derived&>(*this);
  }
  const auto& self_() const noexcept {
    static_assert(std::derived_from<Derived, MatrixViewInterface<Derived>>);
    // static_assert(matrix_view<Derived>);
    return static_cast<const Derived&>(*this);
  }

public:

  /// @brief Shape of the matrix.
  constexpr auto shape() const noexcept {
    return Storm::shape(self_());
  }

  /// @brief Number of the matrix rows.
  constexpr auto num_rows() const noexcept {
    return Storm::num_rows(self_());
  }

  /// @brief Number of the matrix columns.
  constexpr auto num_cols() const noexcept {
    return Storm::num_cols(self_());
  }

  /// @brief Number of the matrix coefficients.
  constexpr auto size() const noexcept {
    return num_rows() * num_cols();
  }

  /// @brief Get the vector coefficient at @p row_index.
  /// @{
  constexpr auto operator()(size_t row_index) noexcept -> decltype(auto) {
    return self_()(row_index, 0);
  }
  constexpr auto operator()(size_t row_index) const noexcept -> decltype(auto) {
    return self_()(row_index, 0);
  }
  /// @}

  /// @brief Get the matrix coefficient at @p row_index and @p col_index.
  /// @{
  constexpr auto operator()(size_t row_index, size_t col_index) noexcept
      -> decltype(auto) {
    return self_()(row_index, col_index);
  }
  constexpr auto operator()(size_t row_index, size_t col_index) const noexcept
      -> decltype(auto) {
    return self_()(row_index, col_index);
  }
  /// @}

}; // class MatrixViewInterface

/// @name Matrix views.
/// @{

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @name Forwarding views.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @{

/// ----------------------------------------------------------------- ///
/// @brief Matrix reference view.
/// ----------------------------------------------------------------- ///
template<matrix Matrix>
class MatrixRefView : public MatrixViewInterface<MatrixRefView<Matrix>> {
private:

  Matrix* mat_;

  static void funс_(Matrix&); // not defined
  static void funс_(Matrix&&) = delete;

public:

  /// @brief Construct a matrix reference view.
  // clang-format off
  template<Detail_::different_from_<MatrixRefView> OtherMatrix>
    requires std::convertible_to<OtherMatrix, Matrix&> &&
             requires { funс_(std::declval<OtherMatrix>()); }
  constexpr MatrixRefView(OtherMatrix&& mat) noexcept 
      : mat_{std::addressof(static_cast<Matrix&>(std::forward<OtherMatrix>(mat)))} {}
  // clang-format on

  /// @copydoc MatrixViewInterface::shape
  constexpr auto shape() const noexcept {
    return Storm::shape(*mat_);
  }

  /// @copydoc MatrixViewInterface::operator()
  /// @{
  constexpr auto operator()(size_t row_index, size_t col_index) noexcept
      -> decltype(auto) {
    return (*mat_)(row_index, col_index);
  }
  constexpr auto operator()(size_t row_index, size_t col_index) const noexcept
      -> decltype(auto) {
    return std::as_const(*mat_)(row_index, col_index);
  }
  /// @}

}; // class MatrixRefView

/// ----------------------------------------------------------------- ///
/// @brief Matrix owning view.
/// ----------------------------------------------------------------- ///
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
    return Storm::shape(mat_);
  }

  /// @copydoc MatrixViewInterface::operator()
  /// @{
  constexpr auto operator()(size_t row_index, size_t col_index) noexcept
      -> decltype(auto) {
    return mat_(row_index, col_index);
  }
  constexpr auto operator()(size_t row_index, size_t col_index) const noexcept
      -> decltype(auto) {
    return mat_(row_index, col_index);
  }
  /// @}

}; // class MatrixOwningView

template<class Matrix>
MatrixRefView(Matrix&) -> MatrixRefView<Matrix>;

namespace Detail_ {
  // clang-format off
  template<class Matrix>
  concept matrix_ref_viewable_ = 
      requires { MatrixRefView{std::declval<Matrix>()}; };
  template<class Matrix>
  concept matrix_ownable_ = 
      requires { MatrixOwningView{std::declval<Matrix>()}; };
  // clang-format on
} // namespace Detail_

/// @brief Forward the viewable matrix as a matrix view.
// clang-format off
template<viewable_matrix Matrix>
  requires matrix_view<std::decay_t<Matrix>> || 
           Detail_::matrix_ref_viewable_<Matrix> || 
           Detail_::matrix_ownable_<Matrix>
constexpr auto forward_as_matrix_view(Matrix&& mat) noexcept {
  // clang-format on
  if constexpr (matrix_view<std::decay_t<Matrix>>) {
    return std::forward<Matrix>(mat);
  } else if constexpr (Detail_::matrix_ref_viewable_<Matrix>) {
    return MatrixRefView{std::forward<Matrix>(mat)};
  } else if constexpr (Detail_::matrix_ownable_<Matrix>) {
    return MatrixOwningView{std::forward<Matrix>(mat)};
  }
}

/// @brief Suitable matrix view type for a vieable matrix.
template<viewable_matrix Matrix>
using forward_as_matrix_view_t =
    decltype(forward_as_matrix_view(std::declval<Matrix>()));

/// @} // Forwarding views.

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @name Generating views.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @{

/// ----------------------------------------------------------------- ///
/// @brief Matrix generating view.
/// ----------------------------------------------------------------- ///
// clang-format off
template<matrix_shape Shape, std::regular_invocable<size_t, size_t> Func>
  requires std::is_object_v<Shape> && std::is_object_v<Func> &&
           Detail_::can_reference_<std::invoke_result_t<Func, size_t, size_t>>
class GenerateMatrixView :
  public MatrixViewInterface<GenerateMatrixView<Shape, Func>> {
  // clang-format on
private:

  STORM_NO_UNIQUE_ADDRESS_ Shape shape_;
  STORM_NO_UNIQUE_ADDRESS_ Func func_;

public:

  /// @brief Construct a generating view.
  constexpr GenerateMatrixView(Shape shape, Func func) noexcept
      : shape_{shape}, func_{std::move(func)} {}

  /// @copydoc MatrixViewInterface::shape
  constexpr auto shape() const noexcept {
    return shape_;
  }

  /// @copydoc MatrixViewInterface::operator()
  constexpr auto operator()(size_t row_index, size_t col_index) const noexcept {
    STORM_ASSERT_(row_index < this->num_rows() &&
                  col_index < this->num_cols() && "Indices are out of range.");
    return func_(row_index, col_index);
  }

}; // class GenerateMatrixView

template<class Shape, class Func>
GenerateMatrixView(Shape, Func) -> GenerateMatrixView<Shape, Func>;

/// @brief Generate a constant matrix with @p num_rows and @p num_cols.
/// @{
template<class Value>
constexpr auto make_constant_matrix( //
    std::convertible_to<size_t> auto num_rows,
    std::convertible_to<size_t> auto num_cols, Value scal) {
  return GenerateMatrixView(
      std::pair(num_rows, num_cols),
      [scal](size_t row_index, size_t col_index) { return scal; });
}
template<size_t NumRows, size_t NumCols, class Value>
constexpr auto make_constant_matrix(Value scal) {
  return make_constant_matrix<Value>(size_t_constant<NumRows>{}, //
                                     size_t_constant<NumCols>{}, scal);
}
/// @}

/// @brief Generate a diagonal matrix with @p num_rows and @p num_cols.
/// @{
template<class Value, class Tag = void>
constexpr auto make_diagonal_matrix( //
    std::convertible_to<size_t> auto num_rows,
    std::convertible_to<size_t> auto num_cols, Value scal) {
  return GenerateMatrixView( //
      std::pair(num_rows, num_cols),
      [scal](size_t row_index, size_t col_index) {
        return row_index == col_index ? scal : Value{};
      });
}
template<size_t NumRows, size_t NumCols, class Value>
constexpr auto make_diagonal_matrix(Value scal) {
  return make_diagonal_matrix<Value>(size_t_constant<NumRows>{}, //
                                     size_t_constant<NumCols>{}, scal);
}
/// @}

/// @} // Generating views.

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @name Slicing views.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @{

/// ----------------------------------------------------------------- ///
/// @brief All indices range.
/// ----------------------------------------------------------------- ///
class AllIndices {
public:

  constexpr auto size() const = delete;

  constexpr size_t operator[](size_t index) const noexcept {
    return index;
  }

}; // AllIndices

/// ----------------------------------------------------------------- ///
/// @brief Selected indices range.
/// ----------------------------------------------------------------- ///
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

/// ----------------------------------------------------------------- ///
/// @brief Sliced indices range.
/// ----------------------------------------------------------------- ///
// clang-format off
template<std::convertible_to<size_t> FromType,
         std::convertible_to<size_t> ToType,
         std::convertible_to<size_t> StrideType>
  requires std::is_object_v<FromType> && std::is_object_v<ToType> && 
           std::is_object_v<StrideType> 
class SlicedIndices {
  // clang-format on
private:

  STORM_NO_UNIQUE_ADDRESS_ FromType from_;
  STORM_NO_UNIQUE_ADDRESS_ ToType to_;
  STORM_NO_UNIQUE_ADDRESS_ StrideType stride_;

  constexpr static auto compute_size_(auto from, auto to, auto stride) {
    return static_cast<size_t>((to - from) / stride);
  }
  template<class FromType_, class ToType_, class StrideType_, //
           FromType_ From, ToType_ To, StrideType_ Stride>
  constexpr static auto
  compute_size_(std::integral_constant<FromType_, From>,
                std::integral_constant<ToType_, To>,
                std::integral_constant<StrideType_, Stride>) {
    return size_t_constant<(To - From) / Stride>{};
  }

public:

  constexpr SlicedIndices( //
      FromType from, ToType to, StrideType stride) noexcept
      : from_{from}, to_{to}, stride_{stride} {}

  constexpr auto size() const noexcept {
    return compute_size_(from_, to_, stride_);
  }

  constexpr size_t operator[](size_t index) const noexcept {
    STORM_ASSERT_(index < size() && "Index is out of range.");
    return from_ + stride_ * index;
  }

}; // class SlicedIndices

/// ----------------------------------------------------------------- ///
/// @brief Submatrix view.
/// ----------------------------------------------------------------- ///
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

public:

  /// @brief Construct a matrix rows view.
  constexpr SubmatrixView(Matrix mat, //
                          RowIndices row_indices,
                          ColIndices col_indices) noexcept
      : mat_{std::move(mat)},                 //
        row_indices_{std::move(row_indices)}, //
        col_indices_{std::move(col_indices)} {}

  /// @copydoc MatrixViewInterface::num_rows
  /// @{
  // clang-format off
  constexpr auto num_rows() const noexcept 
      requires requires { std::declval<RowIndices>().size(); } {
    // clang-format on
    return row_indices_.size();
  }
  constexpr auto num_rows() const noexcept {
    return Storm::num_rows(mat_);
  }
  /// @}

  /// @copydoc MatrixViewInterface::num_cols
  /// @{
  // clang-format off
  constexpr auto num_cols() const noexcept 
      requires requires { std::declval<ColIndices>().size(); } {
    // clang-format on
    return col_indices_.size();
  }
  constexpr auto num_cols() const noexcept {
    return Storm::num_cols(mat_);
  }
  /// @}

  /// @copydoc MatrixViewInterface::shape
  constexpr auto shape() const noexcept {
    return std::pair(num_rows(), num_cols());
  }

  /// @copydoc MatrixViewInterface::operator()
  /// @{
  constexpr auto operator()(size_t row_index, size_t col_index) noexcept
      -> decltype(auto) {
    return mat_(row_indices_[row_index], col_indices_[col_index]);
  }
  constexpr auto operator()(size_t row_index, size_t col_index) const noexcept
      -> decltype(auto) {
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
  STORM_ASSERT_((static_cast<size_t>(row_indices) < mat.num_rows()) && ... &&
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
constexpr auto slice_rows(viewable_matrix auto&& mat,
                          std::convertible_to<size_t> auto rows_from,
                          std::convertible_to<size_t> auto rows_to,
                          std::convertible_to<size_t> auto row_stride =
                              size_t_constant<1>{}) noexcept {
  STORM_ASSERT_(rows_from < rows_to &&
                static_cast<size_t>(rows_to) <= mat.num_rows() &&
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
                static_cast<size_t>(cols_to) <= mat.num_cols() &&
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

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @name Functional views.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @{

/// ----------------------------------------------------------------- ///
/// @brief Component-wise product of function to matrices view.
/// ----------------------------------------------------------------- ///
// clang-format off
template<std::copy_constructible Func, matrix... Matrices>
  requires std::is_object_v<Func> &&
           std::regular_invocable<Func, matrix_reference_t<Matrices>...> &&
           Detail_::can_reference_<
              std::invoke_result_t<Func, matrix_reference_t<Matrices>...>>
class MapMatrixView :
    public MatrixViewInterface<MapMatrixView<Func, Matrices...>> {
  // clang-format on
private:

  STORM_NO_UNIQUE_ADDRESS_ Func func_;
  STORM_NO_UNIQUE_ADDRESS_ std::tuple<Matrices...> mats_;

public:

  /// @brief Construct a map view.
  constexpr MapMatrixView(Func func, Matrices... mats) noexcept
      : func_{std::move(func)}, mats_{std::move(mats)...} {}

  /// @copydoc MatrixViewInterface::shape
  constexpr auto shape() const noexcept {
    return Storm::shape(std::get<0>(mats_));
  }

  /// @copydoc MatrixViewInterface::operator()
  constexpr auto operator()(size_t row_index, size_t col_index) const noexcept
      -> decltype(auto) {
    return std::apply(
        [&](const auto&... mats) {
          return func_(mats(row_index, col_index)...);
        },
        mats_);
  }

}; // class MapMatrixView

template<class Func, class... Matrices>
MapMatrixView(Func, Matrices&&...)
    -> MapMatrixView<Func, forward_as_matrix_view_t<Matrices>...>;

/// @brief Make a component-wise product of function @p func
///   to matrices @p mat1, @p mats view.
constexpr auto map(auto func, //
                   viewable_matrix auto&& mat1,
                   viewable_matrix auto&&... mats) noexcept {
  STORM_ASSERT_(((shape(mat1) == shape(mats)) && ...) &&
                "Shapes of the matrix arguments should be the same.");
  return MapMatrixView(func, STORM_FORWARD_(mat1), STORM_FORWARD_(mats)...);
}

/// @brief "+" the matrix @p mat.
constexpr auto operator+(viewable_matrix auto&& mat) noexcept {
  return map([](auto&& m) { return +STORM_FORWARD_(m); }, STORM_FORWARD_(mat));
}

/// @brief Negate the matrix @p mat.
constexpr auto operator-(viewable_matrix auto&& mat) noexcept {
  return map(std::negate{}, STORM_FORWARD_(mat));
}

/// @brief Multiply the matrix @p mat by a scalar @p scal.
/// @{
// clang-format off
template<class Scalar>
  requires (!matrix<Scalar>)
constexpr auto operator*(Scalar scal, viewable_matrix auto&& mat) noexcept {
  // clang-format on
  return map(
      [scal = std::move(scal)](auto&& m) { return scal * STORM_FORWARD_(m); },
      STORM_FORWARD_(mat));
}
// clang-format off
template<class Scalar>
  requires (!matrix<Scalar>)
constexpr auto operator*(viewable_matrix auto&& mat, Scalar scal) noexcept {
  // clang-format on
  return map(
      [scal = std::move(scal)](auto&& m) { return STORM_FORWARD_(m) * scal; },
      STORM_FORWARD_(mat));
}
/// @}

/// @brief Divide the matrix @p mat by a scalar @p scal.
// clang-format off
template<class Scalar>
  requires (!matrix<Scalar>)
constexpr auto operator/(viewable_matrix auto&& mat, Scalar scal) noexcept {
  // clang-format on
  return map(
      [scal = std::move(scal)](auto&& m) { return STORM_FORWARD_(m) / scal; },
      STORM_FORWARD_(mat));
}

/// @brief Add the matrices @p mat1 and @p mat2.
constexpr auto operator+(viewable_matrix auto&& mat1,
                         viewable_matrix auto&& mat2) noexcept {
  return map(std::plus{}, STORM_FORWARD_(mat1), STORM_FORWARD_(mat2));
}

/// @brief Subtract the matrices @p mat1 and @p mat2.
constexpr auto operator-(viewable_matrix auto&& mat1,
                         viewable_matrix auto&& mat2) noexcept {
  return map(std::minus{}, STORM_FORWARD_(mat1), STORM_FORWARD_(mat2));
}

/// @brief Component-wise multiply the matrices @p mat1 and @p mat2.
constexpr auto operator*(viewable_matrix auto&& mat1,
                         viewable_matrix auto&& mat2) noexcept {
  return map(std::multiplies{}, STORM_FORWARD_(mat1), STORM_FORWARD_(mat2));
}

/// @brief Component-wise divide the matrices @p mat1 and @p mat2.
constexpr auto operator/(viewable_matrix auto&& mat1,
                         viewable_matrix auto&& mat2) noexcept {
  return map(std::divides{}, STORM_FORWARD_(mat1), STORM_FORWARD_(mat2));
}

namespace math {

  /// @brief Component-wise @c abs of the matrix @p mat.
  constexpr auto abs(viewable_matrix auto&& mat) noexcept {
    return map([](auto&& m) { return math::abs(STORM_FORWARD_(m)); }, mat);
  }

  /// @name Power functions.
  /// @{

  // clang-format off
  template<class Scalar>
    requires (!matrix<Scalar>)
  constexpr auto pow(viewable_matrix auto&& x_mat, Scalar y) noexcept {
    // clang-format on
    return map( //
        [y = std::move(y)](auto&& x) {
          return math::pow(STORM_FORWARD_(x), y);
        },
        STORM_FORWARD_(x_mat));
  }

  // clang-format off
  template<class Scalar>
    requires (!matrix<Scalar>)
  constexpr auto pow(Scalar x, viewable_matrix auto&& y_mat) noexcept {
    // clang-format on
    return map( //
        [x = std::move(x)](auto&& y) {
          return math::pow(x, STORM_FORWARD_(y));
        },
        STORM_FORWARD_(y_mat));
  }

  /// @brief Component-wise @c atan2 of the matriсes @p x_mat and @p y_mat.
  constexpr auto pow(viewable_matrix auto&& x_mat,
                     viewable_matrix auto&& y_mat) noexcept {
    return map(
        [](auto&& x, auto&& y) {
          return math::pow(STORM_FORWARD_(x), STORM_FORWARD_(y));
        },
        STORM_FORWARD_(x_mat), STORM_FORWARD_(y_mat));
  }

  /// @brief Component-wise @c sqrt of the matrix @p mat.
  constexpr auto sqrt(viewable_matrix auto&& mat) noexcept {
    return map([](auto&& m) { return math::sqrt(STORM_FORWARD_(m)); },
               STORM_FORWARD_(mat));
  }

  /// @brief Component-wise @c cbrt of the matrix @p mat.
  constexpr auto cbrt(viewable_matrix auto&& mat) noexcept {
    return map([](auto&& m) { return math::cbrt(STORM_FORWARD_(m)); },
               STORM_FORWARD_(mat));
  }

  /// @brief Component-wise @c hypot of the matriсes
  ///   @p x_mat and @p y_mat.
  constexpr auto hypot(viewable_matrix auto&& x_mat,
                       viewable_matrix auto&& y_mat) noexcept {
    return map(
        [](auto&& x, auto&& y) {
          return math::hypot(STORM_FORWARD_(x), STORM_FORWARD_(y));
        },
        STORM_FORWARD_(x_mat), STORM_FORWARD_(y_mat));
  }

  /// @brief Component-wise @c hypot of the matriсes
  ///   @p x_mat, @p y_mat and @p y_mat.
  constexpr auto hypot(viewable_matrix auto&& x_mat,
                       viewable_matrix auto&& y_mat,
                       viewable_matrix auto&& z_mat) noexcept {
    return map(
        [](auto&& x, auto&& y, auto&& z) {
          return math::hypot( //
              STORM_FORWARD_(x), STORM_FORWARD_(y), STORM_FORWARD_(z));
        },
        STORM_FORWARD_(x_mat), STORM_FORWARD_(y_mat), STORM_FORWARD_(z_mat));
  }

  /// @} // Power functions.

  /// @name Exponential functions.
  /// @{

  /// @brief Component-wise @c exp of the matrix @p mat.
  constexpr auto exp(viewable_matrix auto&& mat) noexcept {
    return map([](auto&& m) { return math::exp(STORM_FORWARD_(m)); },
               STORM_FORWARD_(mat));
  }

  /// @brief Component-wise @c exp2 of the matrix @p mat.
  constexpr auto exp2(viewable_matrix auto&& mat) noexcept {
    return map([](auto&& m) { return math::exp2(STORM_FORWARD_(m)); },
               STORM_FORWARD_(mat));
  }

  /// @brief Component-wise @c log of the matrix @p mat.
  constexpr auto log(viewable_matrix auto&& mat) noexcept {
    return map([](auto&& m) { return math::log(STORM_FORWARD_(m)); },
               STORM_FORWARD_(mat));
  }

  /// @brief Component-wise @c log2 of the matrix @p mat.
  constexpr auto log2(viewable_matrix auto&& mat) noexcept {
    return map([](auto&& m) { return math::log2(STORM_FORWARD_(m)); },
               STORM_FORWARD_(mat));
  }

  /// @brief Component-wise @c log10 of the matrix @p mat.
  constexpr auto log10(viewable_matrix auto&& mat) noexcept {
    return map([](auto&& m) { return math::log10(STORM_FORWARD_(m)); },
               STORM_FORWARD_(mat));
  }

  /// @} // Exponential functions.

  /// @name Trigonometric functions.
  /// @{

  /// @brief Component-wise @c sin of the matrix @p mat.
  constexpr auto sin(viewable_matrix auto&& mat) noexcept {
    return map([](auto&& m) { return math::sin(STORM_FORWARD_(m)); },
               STORM_FORWARD_(mat));
  }

  /// @brief Component-wise @c cos of the matrix @p mat.
  constexpr auto cos(viewable_matrix auto&& mat) noexcept {
    return map([](auto&& m) { return math::cos(STORM_FORWARD_(m)); },
               STORM_FORWARD_(mat));
  }

  /// @brief Component-wise @c tan of the matrix @p mat.
  constexpr auto tan(viewable_matrix auto&& mat) noexcept {
    return map([](auto&& m) { return math::tan(STORM_FORWARD_(m)); },
               STORM_FORWARD_(mat));
  }

  /// @brief Component-wise @c asin of the matrix @p mat.
  constexpr auto asin(viewable_matrix auto&& mat) noexcept {
    return map([](auto&& m) { return math::asin(STORM_FORWARD_(m)); },
               STORM_FORWARD_(mat));
  }

  /// @brief Component-wise @c acos of the matrix @p mat.
  constexpr auto acos(viewable_matrix auto&& mat) noexcept {
    return map([](auto&& m) { return math::acos(STORM_FORWARD_(m)); },
               STORM_FORWARD_(mat));
  }

  /// @brief Component-wise @c atan of the matrix @p mat.
  constexpr auto atan(viewable_matrix auto&& mat) noexcept {
    return map([](auto&& m) { return math::atan(STORM_FORWARD_(m)); },
               STORM_FORWARD_(mat));
  }

  /// @brief Component-wise @c atan2 of the matriсes @p y_mat and @p x_mat.
  constexpr auto atan2(viewable_matrix auto&& y_mat,
                       viewable_matrix auto&& x_mat) noexcept {
    return map(
        [](auto&& y, auto&& x) {
          return math::atan2(STORM_FORWARD_(y), STORM_FORWARD_(x));
        },
        STORM_FORWARD_(y_mat), STORM_FORWARD_(x_mat));
  }

  /// @} // Trigonometric functions.

  /// @name Hyperbolic functions.
  /// @{

  /// @brief Component-wise @p sinh of the matrix @p mat.
  constexpr auto sinh(viewable_matrix auto&& mat) noexcept {
    return map([](auto&& m) { return math::sinh(STORM_FORWARD_(m)); },
               STORM_FORWARD_(mat));
  }

  /// @brief Component-wise @c cosh of the matrix @p mat.
  constexpr auto cosh(viewable_matrix auto&& mat) noexcept {
    return map([](auto&& m) { return math::cosh(STORM_FORWARD_(m)); },
               STORM_FORWARD_(mat));
  }

  /// @brief Component-wise @c tanh of the matrix @p mat.
  constexpr auto tanh(viewable_matrix auto&& mat) noexcept {
    return map([](auto&& m) { return math::tanh(STORM_FORWARD_(m)); },
               STORM_FORWARD_(mat));
  }

  /// @brief Component-wise @c asinh of the matrix @p mat.
  constexpr auto asinh(viewable_matrix auto&& mat) noexcept {
    return map([](auto&& m) { return math::asinh(STORM_FORWARD_(m)); },
               STORM_FORWARD_(mat));
  }

  /// @brief Component-wise @c acosh of the matrix @p mat.
  constexpr auto acosh(viewable_matrix auto&& mat) noexcept {
    return map([](auto&& m) { return math::acosh(STORM_FORWARD_(m)); },
               STORM_FORWARD_(mat));
  }

  /// @brief Component-wise @c atanh of the matrix @p mat.
  constexpr auto atanh(viewable_matrix auto&& mat) noexcept {
    return map([](auto&& m) { return math::atanh(STORM_FORWARD_(m)); },
               STORM_FORWARD_(mat));
  }

  /// @} // Hyperbolic functions.

} // namespace math

/// ----------------------------------------------------------------- ///
/// @brief Matrix transpose view.
/// ----------------------------------------------------------------- ///
template<matrix Matrix>
class MatrixTransposeView :
    public MatrixViewInterface<MatrixTransposeView<Matrix>> {
private:

  STORM_NO_UNIQUE_ADDRESS_ Matrix mat_;

public:

  /// @brief Construct a matrix transpose view.
  constexpr explicit MatrixTransposeView(Matrix mat) noexcept
      : mat_{std::move(mat)} {}

  /// @copydoc MatrixViewInterface::shape
  constexpr auto shape() const noexcept {
    return std::pair(num_cols(mat_), num_rows(mat_));
  }

  /// @copydoc MatrixViewInterface::operator()
  /// @{
  constexpr auto operator()(size_t row_index, size_t col_index) noexcept
      -> decltype(auto) {
    return mat_(col_index, row_index);
  }
  constexpr auto operator()(size_t row_index, size_t col_index) const noexcept
      -> decltype(auto) {
    return mat_(col_index, row_index);
  }
  /// @}

}; // MatrixTransposeView

template<class Matrix>
MatrixTransposeView(Matrix&&)
    -> MatrixTransposeView<forward_as_matrix_view_t<Matrix>>;

/// @brief Transpose the matrix @p mat.
constexpr auto transpose(viewable_matrix auto&& mat) noexcept {
  return MatrixTransposeView(STORM_FORWARD_(mat));
}

/// ----------------------------------------------------------------- ///
/// @brief Matrix product view.
/// ----------------------------------------------------------------- ///
template<matrix Matrix1, matrix Matrix2>
class MatrixProductView :
    public MatrixViewInterface<MatrixProductView<Matrix1, Matrix2>> {
private:

  STORM_NO_UNIQUE_ADDRESS_ Matrix1 mat1_;
  STORM_NO_UNIQUE_ADDRESS_ Matrix2 mat2_;

public:

  /// @brief Construct a matrix product view.
  constexpr explicit MatrixProductView(Matrix1 mat1, Matrix2 mat2) noexcept
      : mat1_{std::move(mat1)}, mat2_{std::move(mat2)} {}

  /// @copydoc MatrixViewInterface::shape
  constexpr auto shape() const noexcept {
    return std::pair(num_rows(mat1_), num_cols(mat2_));
  }

  /// @copydoc MatrixViewInterface::operator()
  constexpr auto operator()(size_t row_index, size_t col_index) const noexcept
      -> decltype(auto) {
    const auto cross_size{num_cols(mat1_)};
    auto val{mat1_(row_index, 0) * mat2_(0, col_index)};
    for (size_t cross_index{1}; cross_index < cross_size; ++cross_index) {
      val += mat1_(row_index, cross_index) * mat2_(cross_index, col_index);
    }
    return val;
  }

}; // class MatrixProductView

template<class Matrix1, class Matrix2>
MatrixProductView(Matrix1&&, Matrix2&&)
    -> MatrixProductView<forward_as_matrix_view_t<Matrix1>,
                         forward_as_matrix_view_t<Matrix2>>;

/// @brief Multiply the matrices @p mat1 and @p mat2.
constexpr auto matmul(viewable_matrix auto&& mat1,
                      viewable_matrix auto&& mat2) noexcept {
  STORM_ASSERT_(num_cols(mat1) == num_rows(mat2) &&
                "The first matrix should have the same number of columns "
                "as the second matrix has rows.");
  return MatrixProductView(STORM_FORWARD_(mat1), STORM_FORWARD_(mat2));
}

/// @} // Functional views.

/// @} // Matrix views.

/// @brief Print a @p mat.
std::ostream& operator<<(std::ostream& out, matrix auto&& mat) {
  for (size_t row_index{0}; row_index < mat.num_rows(); ++row_index) {
    out << "( ";
    for (size_t col_index{0}; col_index < mat.num_cols(); ++col_index) {
      out << mat(row_index, col_index) << " ";
    }
    out << ")" << std::endl;
  }
  return out;
}

} // namespace Storm

#include "MatrixAction.hpp"
