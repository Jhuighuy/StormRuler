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

#include <Storm/Base.hpp>

#include <Storm/Utils/Math.hpp>

#include <Storm/Blass/Matrix.hpp>

#include <concepts>
#include <type_traits>
#include <utility>

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

/// @brief CRTP interface to a matrix views.
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

protected:

  /// @brief Construct a matrix view interface.
  constexpr MatrixViewInterface() = default;

public:

  /// @brief Shape of the matrix.
  [[nodiscard]] constexpr matrix_shape_t shape() const noexcept {
    return self_().shape();
  }

  /// @brief Get the vector coefficient at @p row_index.
  /// @{
  [[nodiscard]] constexpr decltype(auto) //
  operator()(size_t row_index) noexcept {
    return self_()(row_index, 0);
  }
  [[nodiscard]] constexpr decltype(auto)
  operator()(size_t row_index) const noexcept {
    return self_()(row_index, 0);
  }
  /// @}

  /// @brief Get the matrix coefficient at @p row_index and @p col_index.
  /// @{
  [[nodiscard]] constexpr decltype(auto) //
  operator()(size_t row_index, size_t col_index) noexcept {
    return self_()(row_index, col_index);
  }
  [[nodiscard]] constexpr decltype(auto)
  operator()(size_t row_index, size_t col_index) const noexcept {
    return self_()(row_index, col_index);
  }
  /// @}

  /// @brief Traverse the matrix expression tree with a @p visitor.
  /// @tparam LeafsToRoot If true, the tree would be traversed from leafs to
  /// root, otherwise from root to leafs.
  template<bool LeafsToRoot = true>
  constexpr void traverse_tree(const auto& visitor) const {
    self_().traverse_tree(visitor);
  }

  /// @brief Transform the matrix expression tree with a @p transformer.
  /// The tree would be traversed and transformed from leafs to root.
  [[nodiscard]] constexpr auto transform_tree(const auto& transformer) const {
    return self_().transform_tree(transformer);
  }

protected:

  template<bool LeafsToRoot>
  constexpr void traverse_self_(const auto& visitor,
                                const matrix auto&... child_nodes) const {
    if constexpr (LeafsToRoot) {
      (visitor(child_nodes), ...), visitor(self_());
    } else {
      visitor(self_()), (visitor(child_nodes), ...);
    }
  }

}; // class MatrixViewInterface

/// @brief Matrix reference view.
template<matrix Matrix>
class MatrixRefView final : public MatrixViewInterface<MatrixRefView<Matrix>> {
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
      : p_mat_{std::addressof(
            static_cast<Matrix&>(std::forward<OtherMatrix>(mat)))} {
    // clang-format on
  }

  /// @copydoc MatrixViewInterface::shape
  [[nodiscard]] constexpr matrix_shape_t shape() const noexcept {
    return p_mat_->shape();
  }

  /// @copydoc MatrixViewInterface::operator()
  /// @{
  [[nodiscard]] constexpr decltype(auto) //
  operator()(size_t row_index, size_t col_index) noexcept {
    return (*p_mat_)(row_index, col_index);
  }
  [[nodiscard]] constexpr decltype(auto)
  operator()(size_t row_index, size_t col_index) const noexcept {
    return std::as_const(*p_mat_)(row_index, col_index);
  }
  /// @}

  /// @copydoc MatrixViewInterface::traverse_tree
  constexpr void traverse_tree(const auto& visitor) const {
    p_mat_->tranverse_tree(visitor);
  }

  /// @copydoc MatrixViewInterface::transform_tree
  [[nodiscard]] constexpr auto transform_tree(const auto& transformer) const {
    return transformer(*p_mat_);
  }

}; // class MatrixRefView

template<class Matrix>
MatrixRefView(Matrix&) -> MatrixRefView<Matrix>;

/// @brief Matrix owning view.
// clang-format off
template<matrix Matrix>
  requires std::movable<Matrix>
class MatrixOwningView final : 
    public MatrixViewInterface<MatrixOwningView<Matrix>> {
  // clang-format on
private:

  STORM_NO_UNIQUE_ADDRESS_ Matrix mat_{};

public:

  /// @brief Construct an owning view.
  constexpr MatrixOwningView(Matrix&& mat) : mat_{std::move(mat)} {}

  /// @copydoc MatrixViewInterface::shape
  [[nodiscard]] constexpr matrix_shape_t shape() const noexcept {
    return mat_.shape();
  }

  /// @copydoc MatrixViewInterface::operator()
  /// @{
  [[nodiscard]] constexpr decltype(auto) //
  operator()(size_t row_index, size_t col_index) noexcept {
    return mat_(row_index, col_index);
  }
  [[nodiscard]] constexpr decltype(auto)
  operator()(size_t row_index, size_t col_index) const noexcept {
    return mat_(row_index, col_index);
  }
  /// @}

  /// @copydoc MatrixViewInterface::traverse_tree
  constexpr void traverse_tree(const auto& visitor) const {
    mat_.tranverse_tree(visitor);
  }

  /// @copydoc MatrixViewInterface::transform_tree
  [[nodiscard]] constexpr auto transform_tree(const auto& transformer) const {
    return transformer(mat_);
  }

}; // class MatrixOwningView

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
[[nodiscard]] constexpr auto forward_as_matrix_view(Matrix&& mat) {
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

} // namespace Storm

#define STORM_INSIDE_MATRIX_VIEW_HPP_
#include <Storm/Blass/MatrixViewFunctional.inl>
#include <Storm/Blass/MatrixViewMaking.inl>
#include <Storm/Blass/MatrixViewSlicing.inl>
#undef STORM_INSIDE_MATRIX_VIEW_HPP_

#include "MatrixAction.hpp"
#include "MatrixIo.hpp"
