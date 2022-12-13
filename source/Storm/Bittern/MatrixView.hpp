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
// FITNESS FOR Allocator PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT
// SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
// OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

#pragma once

#include <Storm/Base.hpp>

#include <Storm/Utils/Crtp.hpp>

#include <Storm/Bittern/Matrix.hpp>

#include <concepts>
#include <type_traits>
#include <utility>

namespace Storm
{

// -----------------------------------------------------------------------------

template<crtp_derived MatrixView>
class MatrixViewInterface;

/// @brief Matrix view: a movable matrix with lazy elements evaluation.
template<class MatrixView>
concept matrix_view =
    matrix<MatrixView> && std::movable<MatrixView> &&
    derived_from_crtp_interface<MatrixView, MatrixViewInterface>;

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
template<crtp_derived MatrixView>
class MatrixViewInterface
{
private:

  constexpr MatrixView& self_() noexcept
  {
    static_assert(std::derived_from<MatrixView, MatrixViewInterface>);
    static_assert(matrix_view<MatrixView>);
    return static_cast<MatrixView&>(*this);
  }
  constexpr const MatrixView& self_() const noexcept
  {
    return const_cast<MatrixViewInterface&>(*this).self_();
  }

public:

  // Something to be added.

}; // class MatrixViewInterface

// -----------------------------------------------------------------------------

/// @brief Matrix reference view.
template<matrix Matrix>
class MatrixRefView final : public MatrixViewInterface<MatrixRefView<Matrix>>
{
private:

  Matrix* p_mat_;

  static consteval void rvalue_shall_not_pass_(Matrix&); // not defined
  static consteval void rvalue_shall_not_pass_(Matrix&&) = delete;

public:

  /// @brief Construct a matrix reference view.
  template<detail_::different_from_<MatrixRefView> OtherMatrix>
    requires (std::convertible_to<OtherMatrix, Matrix&> &&
              requires { rvalue_shall_not_pass_(std::declval<OtherMatrix>()); })
  constexpr MatrixRefView(OtherMatrix&& mat) noexcept // NOLINT
      : p_mat_{std::addressof(
            static_cast<Matrix&>(std::forward<OtherMatrix>(mat)))}
  {
  }

  /// @brief Get the matrix shape.
  constexpr auto shape() const noexcept
  {
    return p_mat_->shape();
  }

  /// @brief Get the matrix element at @p indices.
  /// @{
  template<class... Indices>
    requires compatible_matrix_indices_v<MatrixRefView, Indices...>
  constexpr decltype(auto) operator()(Indices... indices) noexcept
  {
    STORM_ASSERT_(in_range(shape(), indices...), "Indices are out of range!");
    return (*p_mat_)(indices...);
  }
  template<class... Indices>
    requires compatible_matrix_indices_v<MatrixRefView, Indices...>
  constexpr decltype(auto) operator()(Indices... indices) const noexcept
  {
    STORM_ASSERT_(in_range(shape(), indices...), "Indices are out of range!");
    return std::as_const(*p_mat_)(indices...);
  }
  /// @}

}; // class MatrixRefView

template<class Matrix>
MatrixRefView(Matrix&) -> MatrixRefView<Matrix>;

// -----------------------------------------------------------------------------

/// @brief Matrix owning view.
/// Matrix should be noexcept-movable.
template<matrix Matrix>
  requires std::move_constructible<Matrix>
class MatrixOwningView final :
    public MatrixViewInterface<MatrixOwningView<Matrix>>,
    public NonCopyable
{
private:

  STORM_NO_UNIQUE_ADDRESS_ Matrix mat_;

public:

  /// @brief Construct an owning view.
  constexpr MatrixOwningView(Matrix&& mat) noexcept // NOLINT
      : mat_{std::move(mat)}
  {
  }

  /// @brief Get the matrix shape.
  constexpr auto shape() const noexcept
  {
    return mat_.shape();
  }

  /// @brief Get the matrix element at @p indices.
  /// @{
  template<class... Indices>
    requires compatible_matrix_indices_v<MatrixOwningView, Indices...>
  constexpr decltype(auto) operator()(Indices... indices) noexcept
  {
    STORM_ASSERT_(in_range(shape(), indices...), "Indices are out of range!");
    return mat_(indices...);
  }
  template<class... Indices>
    requires compatible_matrix_indices_v<MatrixOwningView, Indices...>
  constexpr decltype(auto) operator()(Indices... indices) const noexcept
  {
    STORM_ASSERT_(in_range(shape(), indices...), "Indices are out of range!");
    return mat_(indices...);
  }
  /// @}

}; // class MatrixOwningView

// -----------------------------------------------------------------------------

namespace detail_
{
  template<class Matrix>
  concept can_matrix_ref_view_ =
      requires { MatrixRefView{std::declval<Matrix>()}; };
  template<class Matrix>
  concept can_matrix_owning_view_ =
      requires { MatrixOwningView{std::declval<Matrix>()}; };
} // namespace detail_

/// @brief Wrap the viewable matrix @p mat into a matrix view.
template<viewable_matrix Matrix>
  requires matrix_view<std::decay_t<Matrix>> ||
           detail_::can_matrix_ref_view_<Matrix> ||
           detail_::can_matrix_owning_view_<Matrix>
constexpr auto make_matrix_view(Matrix&& mat) noexcept
{
  if constexpr (matrix_view<std::decay_t<Matrix>>) {
    return std::forward<Matrix>(mat);
  } else if constexpr (detail_::can_matrix_ref_view_<Matrix>) {
    return MatrixRefView{std::forward<Matrix>(mat)};
  } else if constexpr (detail_::can_matrix_owning_view_<Matrix>) {
    return MatrixOwningView{std::forward<Matrix>(mat)};
  }
}

/// @brief Suitable matrix view type for a viewable matrix.
template<viewable_matrix Matrix>
using matrix_view_t = decltype(make_matrix_view(std::declval<Matrix>()));

// -----------------------------------------------------------------------------

} // namespace Storm
