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

#include <concepts>
#include <span>
#include <type_traits>

namespace Storm {

inline constexpr size_t dynamic_value = std::dynamic_extent;

/// @brief Matrix shape.
template<size_t NumRows = dynamic_value, size_t NumCols = dynamic_value>
class MatrixShape {
private:

  using NumRows_ = std::conditional_t<NumRows == dynamic_value, //
                                      size_t, size_t_constant<NumRows>>;
  using NumCols_ = std::conditional_t<NumCols == dynamic_value, //
                                      size_t, size_t_constant<NumCols>>;

  STORM_NO_UNIQUE_ADDRESS_ NumRows_ num_rows_;
  STORM_NO_UNIQUE_ADDRESS_ NumCols_ num_cols_;

public:

  /// @brief Construct a matrix shape with @p num_rows and @p num_cols.
  constexpr explicit MatrixShape(NumRows_ num_rows = NumRows_{},
                                 NumCols_ num_cols = NumCols_{}) noexcept
      : num_rows_{num_rows}, num_cols_{num_cols} {}

  /// @brief Number of rows.
  constexpr NumRows_ num_rows() const noexcept {
    return num_rows_;
  }

  /// @brief Number of columns.
  constexpr NumCols_ num_cols() const noexcept {
    return num_cols_;
  }

  template<size_t OtherNumRows, size_t OtherNumCols>
  constexpr auto operator==(
      const MatrixShape<OtherNumRows, OtherNumCols>& other) const noexcept {
    return num_rows() == other.num_rows() && num_cols() == other.num_cols();
  }

}; // class MatrixShape

template<size_t NumRows>
MatrixShape(size_t_constant<NumRows>, size_t) -> MatrixShape<NumRows>;

template<size_t NumCols>
MatrixShape(size_t, size_t_constant<NumCols>)
    -> MatrixShape<dynamic_value, NumCols>;

template<size_t NumRows, size_t NumCols>
MatrixShape(size_t_constant<NumRows>, size_t_constant<NumCols>)
    -> MatrixShape<NumRows, NumCols>;

/// @brief Vector shape.
template<size_t NumRows = dynamic_value>
using VectorShape = MatrixShape<NumRows, 1>;

template<class>
inline constexpr bool is_matrix_shape_v = false;
template<size_t NumRows, size_t NumCols>
inline constexpr bool is_matrix_shape_v<MatrixShape<NumRows, NumCols>> = true;

/// @brief Matrix shape.
template<class MatrixShape>
concept matrix_shape = is_matrix_shape_v<MatrixShape>;

/// @brief Matrix extent.
template<class MatrixExtent>
concept matrix_extent = std::convertible_to<MatrixExtent, size_t>;

} // namespace Storm
