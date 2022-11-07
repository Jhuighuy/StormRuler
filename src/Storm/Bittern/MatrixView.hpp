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

#include <Storm/Utils/Crtp.hpp>
#include <Storm/Utils/Math.hpp>

#include <Storm/Bittern/Matrix.hpp>

#include <concepts>
#include <type_traits>
#include <utility>

namespace Storm {

template<crtp_derived Derived>
class MatrixViewInterface;

/// @brief Types, enabled to be a matrix view.
template<class T>
inline constexpr bool enable_matrix_view_v =
    derived_from_crtp_interface<T, MatrixViewInterface>;

/// @brief Matrix view concept.
/// @todo In order to add the `movable` constraint, we
///   need to box the functor inside the `MapMatrixView`.
template<class MatrixView>
concept matrix_view = matrix<MatrixView> && // std::movable<MatrixView> &&
                      enable_matrix_view_v<MatrixView>;

/// @brief Matrix that can be safely casted into a matrix view.
template<class Matrix>
concept viewable_matrix =
    matrix<Matrix> &&
    ((matrix_view<std::remove_cvref_t<Matrix>> &&
      std::constructible_from<std::remove_cvref_t<Matrix>, Matrix>) ||
     (!matrix_view<std::remove_cvref_t<Matrix>> &&
      (std::is_lvalue_reference_v<Matrix> ||
       std::movable<std::remove_reference_t<Matrix>>) ));

/// @brief CRTP interface to a matrix views.
template<crtp_derived Derived>
class MatrixViewInterface {
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
  [[nodiscard]] constexpr auto shape() const noexcept {
    return self_().shape();
  }

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
  template<detail_::different_from_<MatrixRefView> OtherMatrix>
    requires std::convertible_to<OtherMatrix, Matrix&> &&
             requires { funс_(std::declval<OtherMatrix>()); }
  constexpr MatrixRefView(OtherMatrix&& mat) noexcept
      : p_mat_{std::addressof(
            static_cast<Matrix&>(std::forward<OtherMatrix>(mat)))} {}

  /// @copydoc MatrixViewInterface::shape
  [[nodiscard]] constexpr auto shape() const noexcept {
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

}; // class MatrixRefView

template<class Matrix>
MatrixRefView(Matrix&) -> MatrixRefView<Matrix>;

/// @brief Matrix owning view.
template<matrix Matrix>
  requires std::movable<Matrix>
class MatrixOwningView final :
    public MatrixViewInterface<MatrixOwningView<Matrix>> {
private:

  STORM_NO_UNIQUE_ADDRESS_ Matrix mat_{};

public:

  /// @brief Construct an owning view.
  constexpr MatrixOwningView(Matrix&& mat) : mat_{std::move(mat)} {}

  /// @copydoc MatrixViewInterface::shape
  [[nodiscard]] constexpr auto shape() const noexcept {
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

}; // class MatrixOwningView

namespace detail_ {
  template<class Matrix>
  concept matrix_can_ref_view_ =
      requires { MatrixRefView{std::declval<Matrix>()}; };
  template<class Matrix>
  concept matrix_can_owning_view_ =
      requires { MatrixOwningView{std::declval<Matrix>()}; };
} // namespace detail_

/// @brief Forward the viewable matrix @p mat as a matrix view.
template<viewable_matrix Matrix>
  requires matrix_view<std::decay_t<Matrix>> ||
           detail_::matrix_can_ref_view_<Matrix> ||
           detail_::matrix_can_owning_view_<Matrix>
[[nodiscard]] constexpr auto forward_as_matrix_view(Matrix&& mat) {
  if constexpr (matrix_view<std::decay_t<Matrix>>) {
    return std::forward<Matrix>(mat);
  } else if constexpr (detail_::matrix_can_ref_view_<Matrix>) {
    return MatrixRefView{std::forward<Matrix>(mat)};
  } else if constexpr (detail_::matrix_can_owning_view_<Matrix>) {
    return MatrixOwningView{std::forward<Matrix>(mat)};
  }
}

/// @brief Suitable matrix view type for a viewable matrix.
template<viewable_matrix Matrix>
using forward_as_matrix_view_t =
    decltype(forward_as_matrix_view(std::declval<Matrix>()));

} // namespace Storm

#define STORM_INSIDE_MATRIX_VIEW_HPP_
#include <Storm/Bittern/MatrixViewFunctional.inl>
#include <Storm/Bittern/MatrixViewMaking.inl>
#include <Storm/Bittern/MatrixViewSlicing.inl>
#undef STORM_INSIDE_MATRIX_VIEW_HPP_

#include "MatrixAction.hpp"
#include "MatrixIo.hpp"
