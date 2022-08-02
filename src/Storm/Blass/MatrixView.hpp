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

template<class>
inline constexpr size_t tensor_rank{1};

// clang-format off
template<class Derived>
  requires std::is_class_v<Derived> &&
           std::same_as<Derived, std::remove_cv_t<Derived>>
class TensorViewInterface;
// clang-format on

namespace Detail_ {
  // clang-format off
  template<class T, class U>
    requires (!std::same_as<T, TensorViewInterface<U>>)
  void is_derived_from_tensor_view_interface_impl_(
      const T&, const TensorViewInterface<U>&); // not defined
  template<class T>
  concept is_derived_from_tensor_view_interface_ =
      requires(T x) { is_derived_from_tensor_view_interface_impl_(x, x); };
  // clang-format on
} // namespace Detail_

/// @brief Types, enabled to be a tensor view.
template<class T>
inline constexpr bool enable_tensor_view_v =
    Detail_::is_derived_from_tensor_view_interface_<T>;

/// @brief Tensor view concept.
/// @todo In order to add the `movable` constraint, we
///   need to box the functor inside the `MapTensorView`.
template<class TensorView>
concept tensor_view =     //
    matrix<TensorView> && // std::movable<TensorView> &&
    enable_tensor_view_v<TensorView>;

/// @brief Tensor that can be safely casted into a tensor view.
// clang-format off
template<class Tensor>
concept viewable_tensor = matrix<Tensor> &&
    ((tensor_view<std::remove_cvref_t<Tensor>> &&
      std::constructible_from<std::remove_cvref_t<Tensor>, Tensor>) ||
     (!tensor_view<std::remove_cvref_t<Tensor>> &&
      (std::is_lvalue_reference_v<Tensor> ||
       std::movable<std::remove_reference_t<Tensor>>)));
// clang-format on

/// ----------------------------------------------------------------- ///
/// @brief Base class for all matrix views.
/// ----------------------------------------------------------------- ///
// clang-format off
template<class Derived>
  requires std::is_class_v<Derived> &&
           std::same_as<Derived, std::remove_cv_t<Derived>>
class TensorViewInterface {
  // clang-format on
private:

  auto& self_() noexcept {
    static_assert(std::derived_from<Derived, TensorViewInterface<Derived>>);
    // static_assert(tensor_view<Derived>);
    return static_cast<Derived&>(*this);
  }
  const auto& self_() const noexcept {
    static_assert(std::derived_from<Derived, TensorViewInterface<Derived>>);
    // static_assert(tensor_view<Derived>);
    return static_cast<const Derived&>(*this);
  }

public:

  /// @brief Shape of the tensor.
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

  /// @brief Get the tensor coefficient at @p row_index.
  /// @{
  constexpr auto
  operator()(std::convertible_to<size_t> auto&&... indices) noexcept
      -> decltype(auto) {
    return self_()(STORM_FORWARD_(indices)...);
  }
  constexpr auto
  operator()(std::convertible_to<size_t> auto&&... indices) const noexcept
      -> decltype(auto) {
    return self_()(STORM_FORWARD_(indices)...);
  }
  /// @}

}; // class TensorViewInterface

/// @name Tensor views.
/// @{

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @name Forwarding views.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @{

/// ----------------------------------------------------------------- ///
/// @brief Tensor reference view.
/// ----------------------------------------------------------------- ///
template<matrix Tensor>
class TensorRefView : public TensorViewInterface<TensorRefView<Tensor>> {
private:

  Tensor* p_ten_;

  static void funс_(Tensor&); // not defined
  static void funс_(Tensor&&) = delete;

public:

  /// @brief Construct a matrix reference view.
  // clang-format off
  template<Detail_::different_from_<TensorRefView> OtherTensor>
    requires std::convertible_to<OtherTensor, Tensor&> &&
        requires { funс_(std::declval<OtherTensor>()); }
  constexpr TensorRefView(OtherTensor&& ten) noexcept 
      : p_ten_{std::addressof(static_cast<Tensor&>(std::forward<OtherTensor>(ten)))} {
    // clang-format on
  }

  /// @copydoc TensorViewInterface::shape
  constexpr auto shape() const noexcept {
    return Storm::shape(*p_ten_);
  }

  /// @copydoc TensorViewInterface::operator()
  /// @{
  constexpr auto
  operator()(std::convertible_to<size_t> auto&&... indices) noexcept
      -> decltype(auto) {
    return (*p_ten_)(STORM_FORWARD_(indices)...);
  }
  constexpr auto
  operator()(std::convertible_to<size_t> auto&&... indices) const noexcept
      -> decltype(auto) {
    return std::as_const(*p_ten_)(STORM_FORWARD_(indices)...);
  }
  /// @}

}; // class TensorRefView

/// ----------------------------------------------------------------- ///
/// @brief Tensor owning view.
/// ----------------------------------------------------------------- ///
// clang-format off
template<matrix Tensor>
  requires std::movable<Tensor>
class TensorOwningView : public TensorViewInterface<TensorOwningView<Tensor>> {
  // clang-format on
private:

  STORM_NO_UNIQUE_ADDRESS_ Tensor ten_{};

public:

  /// @brief Construct an owning view.
  constexpr TensorOwningView(Tensor&& ten) noexcept : ten_{std::move(ten)} {}

  /// @copydoc TensorViewInterface::shape
  constexpr auto shape() const noexcept {
    return Storm::shape(ten_);
  }

  /// @copydoc TensorViewInterface::operator()
  /// @{
  constexpr auto
  operator()(std::convertible_to<size_t> auto&&... indices) noexcept
      -> decltype(auto) {
    return ten_(STORM_FORWARD_(indices)...);
  }
  constexpr auto
  operator()(std::convertible_to<size_t> auto&&... indices) const noexcept
      -> decltype(auto) {
    return ten_(STORM_FORWARD_(indices)...);
  }
  /// @}

}; // class TensorOwningView

template<class Tensor>
TensorRefView(Tensor&) -> TensorRefView<Tensor>;

namespace Detail_ {
  // clang-format off
  template<class Tensor>
  concept tensor_ref_viewable_ = 
      requires { TensorRefView{std::declval<Tensor>()}; };
  template<class Tensor>
  concept tensor_ownable_ = 
      requires { TensorOwningView{std::declval<Tensor>()}; };
  // clang-format on
} // namespace Detail_

/// @brief Forward the viewable tensor as a tensor view.
// clang-format off
template<viewable_tensor Tensor>
  requires tensor_view<std::decay_t<Tensor>> || 
           Detail_::tensor_ref_viewable_<Tensor> || 
           Detail_::tensor_ownable_<Tensor>
constexpr auto forward_as_tensor_view(Tensor&& ten) noexcept {
  // clang-format on
  if constexpr (tensor_view<std::decay_t<Tensor>>) {
    return std::forward<Tensor>(ten);
  } else if constexpr (Detail_::tensor_ref_viewable_<Tensor>) {
    return TensorRefView{std::forward<Tensor>(ten)};
  } else if constexpr (Detail_::tensor_ownable_<Tensor>) {
    return TensorOwningView{std::forward<Tensor>(ten)};
  }
}

/// @brief Suitable tensor view type for a viewable tensor.
template<viewable_tensor Tensor>
using forward_as_tensor_view_t =
    decltype(forward_as_tensor_view(std::declval<Tensor>()));

/// @} // Forwarding views.

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @name Generating views.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @{

/// ----------------------------------------------------------------- ///
/// @brief Tensor generating view.
/// ----------------------------------------------------------------- ///
// clang-format off
template<matrix_shape Shape, std::regular_invocable<size_t, size_t> Func>
  requires std::is_object_v<Shape> && std::is_object_v<Func> &&
           Detail_::can_reference_<std::invoke_result_t<Func, size_t, size_t>>
class GenerateTensorView :
  public TensorViewInterface<GenerateTensorView<Shape, Func>> {
  // clang-format on
private:

  STORM_NO_UNIQUE_ADDRESS_ Shape shape_;
  STORM_NO_UNIQUE_ADDRESS_ Func func_;

public:

  /// @brief Construct a generating view.
  constexpr GenerateTensorView(Shape shape, Func func) noexcept
      : shape_{shape}, func_{std::move(func)} {}

  /// @copydoc TensorViewInterface::shape
  constexpr auto shape() const noexcept {
    return shape_;
  }

  /// @copydoc TensorViewInterface::operator()
  constexpr auto operator()(size_t row_index, size_t col_index) const noexcept {
    STORM_ASSERT_(row_index < this->num_rows() &&
                  col_index < this->num_cols() && "Indices are out of range.");
    return func_(row_index, col_index);
  }

}; // class GenerateTensorView

template<class Shape, class Func>
GenerateTensorView(Shape, Func) -> GenerateTensorView<Shape, Func>;

/// @brief Generate a constant tensor with @p extents.
/// @{
template<class Value>
constexpr auto
make_constant_matrix(Value scal, std::convertible_to<size_t> auto... extents) {
  return GenerateTensorView(
      std::tuple(extents...),
      [scal]([[maybe_unused]] std::same_as<size_t> auto... indices) {
        return scal;
      });
}
template<size_t... Extents, class Value>
constexpr auto make_constant_matrix(Value scal) {
  return make_constant_matrix<Value>(size_t_constant<Extents>{}..., scal);
}
/// @}

/// @brief Generate a diagonal matrix with @p num_rows and @p num_cols.
/// @{
template<class Value>
constexpr auto make_diagonal_matrix( //
    std::convertible_to<size_t> auto num_rows,
    std::convertible_to<size_t> auto num_cols, Value scal) {
  return GenerateTensorView( //
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
    public TensorViewInterface<SubmatrixView<Matrix, RowIndices, ColIndices>> {
  // clang-format on
private:

  STORM_NO_UNIQUE_ADDRESS_ Matrix mat_;
  STORM_NO_UNIQUE_ADDRESS_ RowIndices row_indices_;
  STORM_NO_UNIQUE_ADDRESS_ ColIndices col_indices_;

public:

  /// @brief Construct a matrix rows view.
  constexpr SubmatrixView(Matrix ten, //
                          RowIndices row_indices,
                          ColIndices col_indices) noexcept
      : mat_{std::move(ten)},                 //
        row_indices_{std::move(row_indices)}, //
        col_indices_{std::move(col_indices)} {}

  /// @copydoc TensorViewInterface::num_rows
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

  /// @copydoc TensorViewInterface::num_cols
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

  /// @copydoc TensorViewInterface::shape
  constexpr auto shape() const noexcept {
    return std::pair(num_rows(), num_cols());
  }

  /// @copydoc TensorViewInterface::operator()
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
    -> SubmatrixView<forward_as_tensor_view_t<Matrix>, RowIndices, ColIndices>;

/// @brief Select the matrix @p ten rows with @p row_indices view.
constexpr auto
select_rows(viewable_tensor auto&& ten,
            std::convertible_to<size_t> auto... row_indices) noexcept {
  STORM_ASSERT_((static_cast<size_t>(row_indices) < ten.num_rows()) && ... &&
                "Row indices are out of range.");
  return SubmatrixView(
      STORM_FORWARD_(ten),
      SelectedIndices(std::array{static_cast<size_t>(row_indices)...}),
      AllIndices{});
}

/// @brief Select the matrix @p ten columns with @p col_index view.
constexpr auto
select_cols(viewable_tensor auto&& ten,
            std::convertible_to<size_t> auto... col_indices) noexcept {
  STORM_ASSERT_((static_cast<size_t>(col_indices) < ten.num_cols()) && ... &&
                "Columns indices are out of range.");
  return SubmatrixView(
      STORM_FORWARD_(ten), //
      AllIndices{},
      SelectedIndices(std::array{static_cast<size_t>(col_indices)...}));
}

/// @brief Slice the matrix @p ten rows from index @p rows_from
///   to index @p rows_to (not including) with a stride @p row_stride view.
/// @{
constexpr auto slice_rows(viewable_tensor auto&& ten,
                          std::convertible_to<size_t> auto rows_from,
                          std::convertible_to<size_t> auto rows_to,
                          std::convertible_to<size_t> auto row_stride =
                              size_t_constant<1>{}) noexcept {
  STORM_ASSERT_(rows_from < rows_to &&
                static_cast<size_t>(rows_to) <= ten.num_rows() &&
                "Invalid row range.");
  return SubmatrixView(STORM_FORWARD_(ten),
                       SlicedIndices(rows_from, rows_to, row_stride),
                       AllIndices{});
}
template<size_t RowsFrom, size_t RowsTo, size_t RowStride = 1>
constexpr auto slice_rows(viewable_tensor auto&& ten) {
  return slice_rows(STORM_FORWARD_(ten), size_t_constant<RowsFrom>{},
                    size_t_constant<RowsTo>{}, size_t_constant<RowStride>{});
}
/// @}

/// @brief Slice the matrix @p ten columns from index @p cols_from
///   to index @p cols_to (not including) with a stride @p col_stride view.
/// @{
constexpr auto slice_cols(viewable_tensor auto&& ten,
                          std::convertible_to<size_t> auto cols_from,
                          std::convertible_to<size_t> auto cols_to,
                          std::convertible_to<size_t> auto col_stride =
                              size_t_constant<1>{}) noexcept {
  STORM_ASSERT_(cols_from < cols_to &&
                static_cast<size_t>(cols_to) <= ten.num_cols() &&
                "Invalid column range.");
  return SubmatrixView(STORM_FORWARD_(ten), //
                       AllIndices{},
                       SlicedIndices(cols_from, cols_to, col_stride));
}
template<size_t ColsFrom, size_t ColsTo, size_t ColStride = 1>
constexpr auto slice_cols(viewable_tensor auto&& ten) {
  return slice_cols(STORM_FORWARD_(ten), size_t_constant<ColsFrom>{},
                    size_t_constant<ColsTo>{}, size_t_constant<ColStride>{});
}
/// @}

constexpr auto slice(viewable_tensor auto&& ten, //
                     std::vector<size_t> rows, size_t cols = 0) noexcept {
  return select_cols(
      slice_rows(STORM_FORWARD_(ten), rows[0], rows[1], size_t_constant<1>{}),
      cols);
}

/// @todo Implement me!
constexpr auto diag(viewable_tensor auto&& ten) noexcept;

/// @todo Implement me!
constexpr auto lower_triangle(viewable_tensor auto&& ten) noexcept;

/// @todo Implement me!
constexpr auto upper_triangle(viewable_tensor auto&& ten) noexcept;

/// @} // Slicing views.

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @name Functional views.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @{

/// ----------------------------------------------------------------- ///
/// @brief Coefficient-wise product of function to tensor view.
/// ----------------------------------------------------------------- ///
// clang-format off
template<std::copy_constructible Func, matrix... Tensors>
  requires std::is_object_v<Func> &&
           std::regular_invocable<Func, matrix_reference_t<Tensors>...> &&
           Detail_::can_reference_<
              std::invoke_result_t<Func, matrix_reference_t<Tensors>...>>
class MapTensorView :
    public TensorViewInterface<MapTensorView<Func, Tensors...>> {
  // clang-format on
private:

  STORM_NO_UNIQUE_ADDRESS_ Func func_;
  STORM_NO_UNIQUE_ADDRESS_ std::tuple<Tensors...> tens_;

public:

  /// @brief Construct a map view.
  constexpr MapTensorView(Func func, Tensors... tens) noexcept
      : func_{std::move(func)}, tens_{std::move(tens)...} {}

  /// @copydoc TensorViewInterface::shape
  constexpr auto shape() const noexcept {
    return Storm::shape(std::get<0>(tens_));
  }

  /// @copydoc TensorViewInterface::operator()
  constexpr auto
  operator()(std::convertible_to<size_t> auto&&... indices) const noexcept
      -> decltype(auto) {
    return std::apply(
        [&](const auto&... tens) {
          return func_(tens(STORM_FORWARD_(indices)...)...);
        },
        tens_);
  }

}; // class MapTensorView

template<class Func, class... Tensors>
MapTensorView(Func, Tensors&&...)
    -> MapTensorView<Func, forward_as_tensor_view_t<Tensors>...>;

/// @brief Make a coefficient-wise product of function @p func
///   to tensors @p ten1, @p tens.
constexpr auto map(auto func, //
                   viewable_tensor auto&& ten1,
                   viewable_tensor auto&&... tens) noexcept {
  STORM_ASSERT_((shape(ten1) == shape(tens)) && ... &&
                "Shapes of the tensor arguments should be the same.");
  return MapTensorView(func, STORM_FORWARD_(ten1), STORM_FORWARD_(tens)...);
}

/// @brief "+" the tensor @p ten.
constexpr auto operator+(viewable_tensor auto&& ten) noexcept {
  return map([](auto&& t) { return +STORM_FORWARD_(t); }, STORM_FORWARD_(ten));
}

/// @brief Negate the tensor @p ten.
constexpr auto operator-(viewable_tensor auto&& ten) noexcept {
  return map(std::negate{}, STORM_FORWARD_(ten));
}

/// @brief Multiply the tensor @p ten by a scalar @p scal.
/// @{
// clang-format off
template<class Scalar>
  requires (!matrix<Scalar>)
constexpr auto operator*(Scalar scal, viewable_tensor auto&& ten) noexcept {
  // clang-format on
  return map(
      [scal = std::move(scal)](auto&& t) { return scal * STORM_FORWARD_(t); },
      STORM_FORWARD_(ten));
}
// clang-format off
template<class Scalar>
  requires (!matrix<Scalar>)
constexpr auto operator*(viewable_tensor auto&& ten, Scalar scal) noexcept {
  // clang-format on
  return map(
      [scal = std::move(scal)](auto&& t) { return STORM_FORWARD_(t) * scal; },
      STORM_FORWARD_(ten));
}
/// @}

/// @brief Divide the tensor @p ten by a scalar @p scal.
// clang-format off
template<class Scalar>
  requires (!matrix<Scalar>)
constexpr auto operator/(viewable_tensor auto&& ten, Scalar scal) noexcept {
  // clang-format on
  return map(
      [scal = std::move(scal)](auto&& t) { return STORM_FORWARD_(t) / scal; },
      STORM_FORWARD_(ten));
}

/// @brief Add the tensors @p ten1 and @p ten2.
constexpr auto operator+(viewable_tensor auto&& ten1,
                         viewable_tensor auto&& ten2) noexcept {
  return map(std::plus{}, STORM_FORWARD_(ten1), STORM_FORWARD_(ten2));
}

/// @brief Subtract the tensors @p ten1 and @p ten2.
constexpr auto operator-(viewable_tensor auto&& ten1,
                         viewable_tensor auto&& ten2) noexcept {
  return map(std::minus{}, STORM_FORWARD_(ten1), STORM_FORWARD_(ten2));
}

/// @brief Coefficient-wise multiply the matrices @p ten1 and @p ten2.
constexpr auto operator*(viewable_tensor auto&& ten1,
                         viewable_tensor auto&& ten2) noexcept {
  return map(std::multiplies{}, STORM_FORWARD_(ten1), STORM_FORWARD_(ten2));
}

/// @brief Coefficient-wise divide the matrices @p ten1 and @p ten2.
constexpr auto operator/(viewable_tensor auto&& ten1,
                         viewable_tensor auto&& ten2) noexcept {
  return map(std::divides{}, STORM_FORWARD_(ten1), STORM_FORWARD_(ten2));
}

namespace math {

  constexpr auto abs(viewable_tensor auto&& ten) noexcept {
    return map([](auto&& t) { return math::abs(STORM_FORWARD_(t)); }, ten);
  }

  /// @name Power functions.
  /// @{

  // clang-format off
  template<class Scalar>
    requires (!matrix<Scalar>)
  constexpr auto pow(viewable_tensor auto&& ten_x, Scalar y) noexcept {
    // clang-format on
    return map( //
        [y = std::move(y)](auto&& x) {
          return math::pow(STORM_FORWARD_(x), y);
        },
        STORM_FORWARD_(ten_x));
  }

  // clang-format off
  template<class Scalar>
    requires (!matrix<Scalar>)
  constexpr auto pow(Scalar x, viewable_tensor auto&& ten_y) noexcept {
    // clang-format on
    return map( //
        [x = std::move(x)](auto&& y) {
          return math::pow(x, STORM_FORWARD_(y));
        },
        STORM_FORWARD_(ten_y));
  }

  constexpr auto pow(viewable_tensor auto&& ten_x,
                     viewable_tensor auto&& ten_y) noexcept {
    return map(
        [](auto&& x, auto&& y) {
          return math::pow(STORM_FORWARD_(x), STORM_FORWARD_(y));
        },
        STORM_FORWARD_(ten_x), STORM_FORWARD_(ten_y));
  }

  constexpr auto sqrt(viewable_tensor auto&& ten) noexcept {
    return map([](auto&& t) { return math::sqrt(STORM_FORWARD_(t)); },
               STORM_FORWARD_(ten));
  }

  constexpr auto cbrt(viewable_tensor auto&& ten) noexcept {
    return map([](auto&& t) { return math::cbrt(STORM_FORWARD_(t)); },
               STORM_FORWARD_(ten));
  }

  constexpr auto hypot(viewable_tensor auto&& ten_x,
                       viewable_tensor auto&& ten_y) noexcept {
    return map(
        [](auto&& x, auto&& y) {
          return math::hypot(STORM_FORWARD_(x), STORM_FORWARD_(y));
        },
        STORM_FORWARD_(ten_x), STORM_FORWARD_(ten_y));
  }
  constexpr auto hypot(viewable_tensor auto&& ten_x,
                       viewable_tensor auto&& ten_y,
                       viewable_tensor auto&& ten_z) noexcept {
    return map(
        [](auto&& x, auto&& y, auto&& z) {
          return math::hypot( //
              STORM_FORWARD_(x), STORM_FORWARD_(y), STORM_FORWARD_(z));
        },
        STORM_FORWARD_(ten_x), STORM_FORWARD_(ten_y), STORM_FORWARD_(ten_z));
  }

  /// @} // Power functions.

  /// @name Exponential functions.
  /// @{

  constexpr auto exp(viewable_tensor auto&& ten) noexcept {
    return map([](auto&& t) { return math::exp(STORM_FORWARD_(t)); },
               STORM_FORWARD_(ten));
  }

  constexpr auto exp2(viewable_tensor auto&& ten) noexcept {
    return map([](auto&& t) { return math::exp2(STORM_FORWARD_(t)); },
               STORM_FORWARD_(ten));
  }

  constexpr auto log(viewable_tensor auto&& ten) noexcept {
    return map([](auto&& t) { return math::log(STORM_FORWARD_(t)); },
               STORM_FORWARD_(ten));
  }

  constexpr auto log2(viewable_tensor auto&& ten) noexcept {
    return map([](auto&& t) { return math::log2(STORM_FORWARD_(t)); },
               STORM_FORWARD_(ten));
  }

  constexpr auto log10(viewable_tensor auto&& ten) noexcept {
    return map([](auto&& t) { return math::log10(STORM_FORWARD_(t)); },
               STORM_FORWARD_(ten));
  }

  /// @} // Exponential functions.

  /// @name Trigonometric functions.
  /// @{

  constexpr auto sin(viewable_tensor auto&& ten) noexcept {
    return map([](auto&& t) { return math::sin(STORM_FORWARD_(t)); },
               STORM_FORWARD_(ten));
  }

  constexpr auto cos(viewable_tensor auto&& ten) noexcept {
    return map([](auto&& t) { return math::cos(STORM_FORWARD_(t)); },
               STORM_FORWARD_(ten));
  }

  constexpr auto tan(viewable_tensor auto&& ten) noexcept {
    return map([](auto&& t) { return math::tan(STORM_FORWARD_(t)); },
               STORM_FORWARD_(ten));
  }

  constexpr auto asin(viewable_tensor auto&& ten) noexcept {
    return map([](auto&& t) { return math::asin(STORM_FORWARD_(t)); },
               STORM_FORWARD_(ten));
  }

  constexpr auto acos(viewable_tensor auto&& ten) noexcept {
    return map([](auto&& t) { return math::acos(STORM_FORWARD_(t)); },
               STORM_FORWARD_(ten));
  }

  constexpr auto atan(viewable_tensor auto&& ten) noexcept {
    return map([](auto&& t) { return math::atan(STORM_FORWARD_(t)); },
               STORM_FORWARD_(ten));
  }

  constexpr auto atan2(viewable_tensor auto&& ten_y,
                       viewable_tensor auto&& ten_x) noexcept {
    return map(
        [](auto&& y, auto&& x) {
          return math::atan2(STORM_FORWARD_(y), STORM_FORWARD_(x));
        },
        STORM_FORWARD_(ten_y), STORM_FORWARD_(ten_x));
  }

  /// @} // Trigonometric functions.

  /// @name Hyperbolic functions.
  /// @{

  constexpr auto sinh(viewable_tensor auto&& ten) noexcept {
    return map([](auto&& t) { return math::sinh(STORM_FORWARD_(t)); },
               STORM_FORWARD_(ten));
  }

  constexpr auto cosh(viewable_tensor auto&& ten) noexcept {
    return map([](auto&& t) { return math::cosh(STORM_FORWARD_(t)); },
               STORM_FORWARD_(ten));
  }

  constexpr auto tanh(viewable_tensor auto&& ten) noexcept {
    return map([](auto&& t) { return math::tanh(STORM_FORWARD_(t)); },
               STORM_FORWARD_(ten));
  }

  constexpr auto asinh(viewable_tensor auto&& ten) noexcept {
    return map([](auto&& t) { return math::asinh(STORM_FORWARD_(t)); },
               STORM_FORWARD_(ten));
  }

  constexpr auto acosh(viewable_tensor auto&& ten) noexcept {
    return map([](auto&& t) { return math::acosh(STORM_FORWARD_(t)); },
               STORM_FORWARD_(ten));
  }

  constexpr auto atanh(viewable_tensor auto&& ten) noexcept {
    return map([](auto&& t) { return math::atanh(STORM_FORWARD_(t)); },
               STORM_FORWARD_(ten));
  }

  /// @} // Hyperbolic functions.

} // namespace math

/// ----------------------------------------------------------------- ///
/// @brief Matrix transpose view.
/// ----------------------------------------------------------------- ///
template<matrix Matrix>
class MatrixTransposeView :
    public TensorViewInterface<MatrixTransposeView<Matrix>> {
private:

  STORM_NO_UNIQUE_ADDRESS_ Matrix mat_;

public:

  /// @brief Construct a matrix transpose view.
  constexpr explicit MatrixTransposeView(Matrix ten) noexcept
      : mat_{std::move(ten)} {}

  /// @copydoc TensorViewInterface::shape
  constexpr auto shape() const noexcept {
    return std::pair(num_cols(mat_), num_rows(mat_));
  }

  /// @copydoc TensorViewInterface::operator()
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
    -> MatrixTransposeView<forward_as_tensor_view_t<Matrix>>;

/// @brief Transpose the matrix @p ten.
constexpr auto transpose(viewable_tensor auto&& ten) noexcept {
  return MatrixTransposeView(STORM_FORWARD_(ten));
}

/// ----------------------------------------------------------------- ///
/// @brief Matrix product view.
/// ----------------------------------------------------------------- ///
template<matrix Matrix1, matrix Matrix2>
class MatrixProductView :
    public TensorViewInterface<MatrixProductView<Matrix1, Matrix2>> {
private:

  STORM_NO_UNIQUE_ADDRESS_ Matrix1 mat1_;
  STORM_NO_UNIQUE_ADDRESS_ Matrix2 mat2_;

public:

  /// @brief Construct a matrix product view.
  constexpr explicit MatrixProductView(Matrix1 mat1, Matrix2 mat2) noexcept
      : mat1_{std::move(mat1)}, mat2_{std::move(mat2)} {}

  /// @copydoc TensorViewInterface::shape
  constexpr auto shape() const noexcept {
    return std::pair(num_rows(mat1_), num_cols(mat2_));
  }

  /// @copydoc TensorViewInterface::operator()
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
    -> MatrixProductView<forward_as_tensor_view_t<Matrix1>,
                         forward_as_tensor_view_t<Matrix2>>;

/// @brief Multiply the matrices @p mat1 and @p mat2.
constexpr auto matmul(viewable_tensor auto&& mat1,
                      viewable_tensor auto&& mat2) noexcept {
  STORM_ASSERT_(num_cols(mat1) == num_rows(mat2) &&
                "The first matrix should have the same number of columns "
                "as the second matrix has rows.");
  return MatrixProductView(STORM_FORWARD_(mat1), STORM_FORWARD_(mat2));
}

/// @} // Functional views.

/// @} // Tensor views.

/// @brief Print a @p ten.
std::ostream& operator<<(std::ostream& out, matrix auto&& ten) {
  for (size_t row_index{0}; row_index < ten.num_rows(); ++row_index) {
    out << "( ";
    for (size_t col_index{0}; col_index < ten.num_cols(); ++col_index) {
      out << ten(row_index, col_index) << " ";
    }
    out << ")" << std::endl;
  }
  return out;
}

} // namespace Storm

#include "MatrixAction.hpp"
