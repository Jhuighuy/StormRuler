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

// Automatically included with `Matrix.hpp`.
#pragma once

#include <Storm/Base.hpp>

#include <Storm/Bittern/Matrix.hpp>

#include <concepts>
#include <type_traits>
#include <utility>

namespace Storm {

/// @brief Matrix making view.
template<std::copy_constructible Func>
  requires std::is_object_v<Func> &&
           std::regular_invocable<Func, size_t, size_t>
class MakeMatrixView final : public MatrixViewInterface<MakeMatrixView<Func>> {
private:

  MatrixShape _shape;
  STORM_NO_UNIQUE_ADDRESS Func _func;

public:

  /// @brief Construct a making view.
  constexpr MakeMatrixView(MatrixShape shape, Func func)
      : _shape{shape}, _func{std::move(func)} {}

  /// @brief Get the matrix shape.
  constexpr auto shape() const noexcept {
    return _shape;
  }

  /// @brief Get the matrix element at @p indices.
  constexpr auto operator()(size_t row_index, size_t col_index) const noexcept {
    STORM_ASSERT(_shape.in_range(row_index, col_index),
                 "Indices are out of range.");
    return _func(row_index, col_index);
  }

}; // class MakeMatrixView

// -----------------------------------------------------------------------------

/// @brief Make a constant matrix of @p shape.
/// @param value Matrix element value.
template<std::copyable Element>
constexpr auto make_constant_matrix(MatrixShape shape, Element value) {
  return MakeMatrixView( //
      shape,
      [value = std::move(value)]( //
          [[maybe_unused]] size_t row_index,
          [[maybe_unused]] size_t col_index) noexcept -> const Element& {
        return value;
      });
}

/// @brief Make a diagonal matrix of @p shape.
/// @param diagonal Matrix diagonal element value.
/// @param off_diagonal Matrix off-diagonal element value.
template<std::copyable Element>
constexpr auto make_diagonal_matrix(MatrixShape shape, Element diagonal,
                                    Element off_diagonal = Element{}) {
  return MakeMatrixView(
      shape,
      [diagonal = std::move(diagonal), //
       off_diagonal = std::move(off_diagonal)](
          size_t row_index, size_t col_index) noexcept -> const Element& {
        return (row_index == col_index) ? diagonal : off_diagonal;
      });
}

/// @todo Find a better place (and a better name) for me.
constexpr auto& fill_diag_with(matrix auto&& mat, auto scal) noexcept {
  return mat <<= make_diagonal_matrix(mat.shape(), scal);
}

} // namespace Storm
