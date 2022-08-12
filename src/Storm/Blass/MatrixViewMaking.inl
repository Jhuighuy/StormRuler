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

#include <concepts>
#include <type_traits>
#include <utility>

namespace Storm {

/// @brief Matrix generating view.
/// @todo This is completely overengineered.
// clang-format off
template<class Shape, std::copy_constructible Func>
  requires is_matrix_shape_v<Shape> && std::is_object_v<Func> &&
           std::regular_invocable<Func, size_t, size_t>
class MakeMatrixView final : 
    public MatrixViewInterface<MakeMatrixView<Shape, Func>> {
  // clang-format on
private:

  STORM_NO_UNIQUE_ADDRESS_ Shape shape_;
  STORM_NO_UNIQUE_ADDRESS_ Func func_;

public:

  /// @brief Construct a generating view.
  constexpr MakeMatrixView(Shape shape, Func func)
      : shape_{shape}, func_{std::move(func)} {}

  /// @copydoc MatrixViewInterface::shape
  [[nodiscard]] constexpr auto shape() const noexcept {
    return shape_;
  }

  /// @copydoc MatrixViewInterface::operator()
  [[nodiscard]] constexpr decltype(auto)
  operator()(size_t row_index, size_t col_index) const noexcept {
    STORM_ASSERT_(row_index < shape_.num_rows() &&
                  col_index < shape_.num_cols() && "Indices are out of range.");
    return func_(row_index, col_index);
  }

  /// @copydoc MatrixViewInterface::traverse_tree
  constexpr void traverse_tree(const auto& visitor) const {
    visitor(*this);
  }

  /// @copydoc MatrixViewInterface::transform_tree
  [[nodiscard]] constexpr auto transform_tree(const auto& transformer) const {
    return transformer(*this);
  }

}; // class MakeMatrixView

/// @brief Generate a constant matrix with @p num_rows and @p num_cols.
/// @param value Matrix element value.
template<std::copyable Value>
[[nodiscard]] constexpr auto make_constant_matrix(matrix_extent auto num_rows,
                                                  matrix_extent auto num_cols,
                                                  Value value) {
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
[[nodiscard]] constexpr auto make_constant_matrix(Value value) {
  return make_constant_matrix<Value>( //
      size_t_constant<NumRows>{}, size_t_constant<NumCols>{}, std::move(value));
}

/// @brief Generate a constant vector with @p num_rows.
/// @param value Vector element value.
template<std::copyable Value>
[[nodiscard]] constexpr auto make_constant_vector(matrix_extent auto num_rows,
                                                  Value value) {
  return make_constant_matrix<Value>( //
      num_rows, size_t_constant<1>{}, std::move(value));
}

/// @brief Generate a fixed-sized constant vector.
/// @tparam NumRows Number of the vector rows.
/// @param value Vector element value.
template<size_t NumRows, std::copyable Value>
[[nodiscard]] constexpr auto make_constant_matrix(Value value) {
  return make_constant_vector<Value>( //
      size_t_constant<NumRows>{}, std::move(value));
}

/// @brief Generate a diagonal matrix with @p num_rows, @p num_cols.
/// @param diagonal Matrix diagonal element value.
/// @param off_diagonal Matrix off-diagonal element value.
template<std::copyable Value>
[[nodiscard]] constexpr auto
make_diagonal_matrix(matrix_extent auto num_rows,
                     matrix_extent auto num_cols, //
                     Value diagonal, Value off_diagonal = Value{}) {
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
[[nodiscard]] constexpr auto
make_diagonal_matrix(Value diagonal, Value off_diagonal = Value{}) {
  return make_diagonal_matrix<Value>(
      size_t_constant<NumRows>{}, size_t_constant<NumCols>{}, //
      std::move(diagonal), std::move(off_diagonal));
}

} // namespace Storm
