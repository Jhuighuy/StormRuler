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

#include <Storm/Bittern/Functions.hpp>
#include <Storm/Bittern/Matrix.hpp>
#include <Storm/Bittern/MatrixAlgorithms.hpp>
#include <Storm/Bittern/MatrixView.hpp>

#include <concepts>
#include <utility>

namespace Storm
{

// -----------------------------------------------------------------------------

template<crtp_derived TargetMatrix>
class TargetMatrixInterface;

/// @brief Target matrix: an output matrix that is an assignment target.
/// Additional requirement: target matrix overrides (or uses) the operator=
/// from the TargetMatrixInterface.
template<class TargetMatrix>
concept target_matrix =
    output_matrix<TargetMatrix> &&
    derived_from_crtp_interface<TargetMatrix, TargetMatrixInterface>;

/// @brief CRTP interface to a target matrix.
template<crtp_derived TargetMatrix>
class TargetMatrixInterface
{
private:

  constexpr TargetMatrix& self_() noexcept
  {
    static_assert(std::derived_from<TargetMatrix, TargetMatrixInterface>);
    static_assert(target_matrix<TargetMatrix>);
    return static_cast<TargetMatrix&>(*this);
  }
  constexpr const TargetMatrix& self_() const noexcept
  {
    return const_cast<TargetMatrixInterface&>(*this).self_();
  }

public:

  /// @brief Assign the current matrix coefficients from matrix @p mat.
  template<matrix Matrix>
    requires std::assignable_from< //
        matrix_element_ref_t<TargetMatrix>, matrix_element_t<Matrix>>
  constexpr decltype(auto) operator=(Matrix&& mat) noexcept
  {
    return assign(self_(), std::forward<Matrix>(mat));
  }

  /// @brief Multiply-assign the current matrix by a scalar @p scal.
  template<std::copyable Scalar>
    requires numeric_matrix<TargetMatrix> && numeric_type<Scalar>
  constexpr decltype(auto) operator*=(Scalar scal)
  {
    return assign(BindLast{MultiplyAssign{}, std::move(scal)}, *this);
  }

  /// @brief Divide-assign the current matrix by a scalar @p scal.
  template<std::copyable Scalar>
    requires numeric_matrix<TargetMatrix> && numeric_type<Scalar>
  constexpr decltype(auto) operator/=(Scalar scal)
  {
    return assign(BindLast{DivideAssign{}, std::move(scal)}, *this);
  }

  /// @brief Add-assign the matrix @p mat to the current matrix.
  template<matrix Matrix>
    requires compatible_matrices_v<TargetMatrix, Matrix> && //
             numeric_matrix<TargetMatrix> && numeric_matrix<Matrix>
  constexpr decltype(auto) operator+=(Matrix&& mat)
  {
    return assign(AddAssign{}, self_(), std::forward<Matrix>(mat));
  }

  /// @brief Subtract-assign the matrix @p mat from the current matrix.
  template<matrix Matrix>
    requires compatible_matrices_v<TargetMatrix, Matrix> && //
             numeric_matrix<TargetMatrix> && numeric_matrix<Matrix>
  constexpr decltype(auto) operator-=(Matrix&& mat)
  {
    return assign(SubtractAssign{}, self_(), std::forward<Matrix>(mat));
  }

  /// @brief Element-wise multiply-assign the current matrix by
  /// the matrix @p mat.
  template<matrix Matrix>
    requires compatible_matrices_v<TargetMatrix, Matrix> && //
             numeric_matrix<TargetMatrix> && numeric_matrix<Matrix>
  constexpr decltype(auto) operator*=(Matrix&& mat)
  {
    return assign(MultiplyAssign{}, self_(), std::forward<Matrix>(mat));
  }

  /// @brief Element-wise divide-assign the current matrix by the matrix @p mat.
  template<matrix Matrix>
    requires compatible_matrices_v<TargetMatrix, Matrix> && //
             numeric_matrix<TargetMatrix> && numeric_matrix<Matrix>
  constexpr decltype(auto) operator/=(Matrix&& mat)
  {
    return assign(DivideAssign{}, self_(), std::forward<Matrix>(mat));
  }

}; // class TargetMatrixInterface

// -----------------------------------------------------------------------------

/// @brief Output matrix as target wrapper.
template<matrix_view Matrix>
  requires output_matrix<Matrix>
class TargetMatrixView final :
    public MatrixViewInterface<TargetMatrixView<Matrix>>,
    public TargetMatrixInterface<TargetMatrixView<Matrix>>,
    public NonMovable
{
private:

  STORM_NO_UNIQUE_ADDRESS_ Matrix mat_;

public:

  /// @brief Construct a matrix target view.
  constexpr explicit TargetMatrixView(Matrix mat) : mat_{std::move(mat)} {}

  /// @brief Get the matrix shape.
  constexpr auto shape() const noexcept
  {
    return mat_.shape();
  }

  /// @brief Get the matrix element at @p indices.
  template<class... Indices>
    requires compatible_matrix_indices_v<TargetMatrixView, Indices...>
  constexpr auto operator()(Indices... indices) const noexcept
  {
    STORM_ASSERT_(in_range(shape(), indices...), "Indices are out of range!");
    return mat_(indices...);
  }

  using TargetMatrixInterface<TargetMatrixView<Matrix>>::operator=;

}; // class TargetMatrixView

template<class Matrix>
TargetMatrixView(Matrix) -> TargetMatrixView<matrix_view_t<Matrix>>;

/// @brief Wrap the output matrix @p mat as a target.
template<viewable_matrix Matrix>
  requires output_matrix<Matrix>
constexpr auto as_target(Matrix&& mat)
{
  return TargetMatrixView{std::forward<Matrix>(mat)};
}

// -----------------------------------------------------------------------------

} // namespace Storm
