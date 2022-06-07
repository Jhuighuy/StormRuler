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

#include <concepts>
#include <ostream>
#include <type_traits>
#include <utility>

#include <stormBase.hxx>
#include <stormSolvers/MatrixBase.hxx>
#include <stormUtils/Math.hxx>

namespace Storm {

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Matrix view.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class RowsSize, class ColsSize, class Indexable>
class MatrixView final {
private:

  [[no_unique_address]] RowsSize num_rows_;
  [[no_unique_address]] ColsSize num_cols_;
  [[no_unique_address]] Indexable indexable_;

public:

  /// @brief Construct a matrix view.
  /// @{
  constexpr MatrixView(RowsSize num_rows, ColsSize num_cols,
                       Indexable&& indexable)
      : num_rows_{num_rows}, num_cols_{num_cols},
        indexable_{std::forward<Indexable>(indexable)} {}
  /// @}

  /// @brief Number of the matrix rows.
  constexpr auto num_rows() const noexcept {
    return num_rows_;
  }

  /// @brief Number of the matrix columns.
  constexpr auto num_cols() const noexcept {
    return num_cols_;
  }

  /// @brief Get the coefficient at @p row_index and @p col_index.
  /// @{
  constexpr decltype(auto) operator()(size_t row_index,
                                      size_t col_index) noexcept {
    STORM_ASSERT_(row_index < num_rows_ && col_index < num_cols_ &&
                  "Indices are out of range.");
    return indexable_(row_index, col_index);
  }
  constexpr decltype(auto) operator()(size_t row_index,
                                      size_t col_index) const noexcept {
    static_cast<void>(const_cast<MatrixView&>(*this)(row_index, col_index));
    return indexable_(row_index, col_index);
  }
  /// @}

}; // class MatrixView

template<class RowsSizeRef, class ColsSizeRef, class IndexableRef>
MatrixView(RowsSizeRef, ColsSizeRef, IndexableRef)
    -> MatrixView<std::decay_t<RowsSizeRef>, std::decay_t<ColsSizeRef>,
                  std::decay_t<IndexableRef>>;

template<class RowsSize, class ColsSize, class Indexable>
struct is_matrix_view_t<MatrixView<RowsSize, ColsSize, Indexable>> :
    std::true_type {};

/// @name Matrix views.
/// @{

/// @name General views.
/// @{

/// @brief Wrap the matrix @p mat into a view.
/// @{
constexpr auto as_view(is_matrix_view auto& mat) {
  // Capture the unknown view by reference: it may be a container.
  auto reference_body{
      [&](size_t row_index, size_t col_index) -> decltype(auto) {
        return mat(row_index, col_index);
      }};
  return MatrixView(mat.num_rows(), mat.num_cols(), std::move(reference_body));
}
constexpr auto as_view(const is_matrix_view auto& mat) {
  // Capture the unknown view by reference: it may be a container.
  auto reference_body{
      [&](size_t row_index, size_t col_index) -> decltype(auto) {
        return mat(row_index, col_index);
      }};
  return MatrixView(mat.num_rows(), mat.num_cols(), std::move(reference_body));
}
template<class RowsSize, class ColsSize, class Indexable>
constexpr auto as_view(const MatrixView<RowsSize, ColsSize, Indexable>& mat) {
  // Force copy a temporary view.
  return mat;
}
/// @}

/// @brief Component-wise apply a @p func
///   to the matrix arguments @p mat1, @p mats.
constexpr auto apply(auto func, const is_matrix_view auto& mat1,
                     const is_matrix_view auto&... mats) noexcept {
  STORM_ASSERT_(((mat1.num_rows() == mats.num_rows()) && ...) &&
                ((mat1.num_cols() == mats.num_cols()) && ...) &&
                "Shapes of the matrix arguments should be the same.");
  auto apply_body{[func, mat1 = as_view(mat1), ... mats = as_view(mats)](
                      size_t row_index, size_t col_index) {
    return func(mat1(row_index, col_index), mats(row_index, col_index)...);
  }};
  return MatrixView(mat1.num_rows(), mat1.num_cols(), std::move(apply_body));
}

/// @}

/// @name Arithmetic operations views.
/// @{

/// @brief "+" the matrix @p mat.
constexpr auto operator+(const is_matrix_view auto& mat) noexcept {
  return apply([](const auto& val) { return +val; }, mat);
}

/// @brief Negate the matrix @p mat.
constexpr auto operator-(const is_matrix_view auto& mat) noexcept {
  return apply([](const auto& val) { return -val; }, mat);
}

/// @brief Multiply the matrix @p mat by a scalar @p scal.
/// @{
constexpr auto operator*(auto scal, const is_matrix_view auto& mat) noexcept {
  return apply([scal](const auto& val) { return scal * val; }, mat);
}
constexpr auto operator*(const is_matrix_view auto& mat, auto scal) noexcept {
  return apply([scal](const auto& val) { return val * scal; }, mat);
}
/// @}

/// @brief Divide the matrix @p mat by a scalar @p scal.
constexpr auto operator/(const is_matrix_view auto& mat,
                         const auto& scal) noexcept {
  return apply([scal](const auto& val) { return scal * val; }, mat);
}

/// @brief Add the matrices @p mat1 and @p mat2.
constexpr auto operator+(const is_matrix_view auto& mat1,
                         const is_matrix_view auto& mat2) noexcept {
  return apply([](const auto& val1, const auto& val2) { return val1 + val2; },
               mat1, mat2);
}

/// @brief Subtract the matrices @p mat1 and @p mat2.
constexpr auto operator-(const is_matrix_view auto& mat1,
                         const is_matrix_view auto& mat2) noexcept {
  return apply([](const auto& val1, const auto& val2) { return val1 - val2; },
               mat1, mat2);
}

/// @brief Component-wise multiply the matrices @p mat1 and @p mat2.
constexpr auto operator*(const is_matrix_view auto& mat1,
                         const is_matrix_view auto& mat2) noexcept {
  return apply([](const auto& val1, const auto& val2) { return val1 * val2; },
               mat1, mat2);
}

/// @brief Component-wise divide the matrices @p mat1 and @p mat2.
constexpr auto operator/(const is_matrix_view auto& mat1,
                         const is_matrix_view auto& mat2) noexcept {
  return apply([](const auto& val1, const auto& val2) { return val1 / val2; },
               mat1, mat2);
}

/// @}

/// @name Math operations views.
/// @{

namespace math {

  /// @name Basic operations.
  /// @{

  constexpr auto abs(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::abs(val); }, mat);
  }

  /// @}

  /// @name Power functions.
  /// @{

  constexpr auto pow(const is_matrix_view auto& x_mat, auto y) noexcept {
    return apply([y](const auto& x) { return math::pow(x, y); }, x_mat);
  }

  constexpr auto pow(auto x, const is_matrix_view auto& y_mat) noexcept {
    return apply([x](const auto& y) { return math::pow(x, y); }, y_mat);
  }

  constexpr auto pow(const is_matrix_view auto& x_mat,
                     const is_matrix_view auto& y_mat) noexcept {
    return apply([](const auto& x, const auto& y) { return math::pow(x, y); },
                 x_mat, y_mat);
  }

  constexpr auto sqrt(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::sqrt(val); }, mat);
  }

  constexpr auto cbrt(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::cbrt(val); }, mat);
  }

  constexpr auto hypot(const is_matrix_view auto& x_mat,
                       const is_matrix_view auto& y_mat) noexcept {
    return apply([](const auto& x, const auto& y) { return math::hypot(x, y); },
                 x_mat, y_mat);
  }

  constexpr auto hypot(const is_matrix_view auto& x_mat,
                       const is_matrix_view auto& y_mat,
                       const is_matrix_view auto& z_mat) noexcept {
    return apply([](const auto& x, const auto& y,
                    const auto& z) { return math::hypot(x, y, z); },
                 x_mat, y_mat, z_mat);
  }

  /// @}

  /// @name Exponential functions.
  /// @{

  constexpr auto exp(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::exp(val); }, mat);
  }

  constexpr auto exp2(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::exp2(val); }, mat);
  }

  constexpr auto log(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::log(val); }, mat);
  }

  constexpr auto log2(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::log2(val); }, mat);
  }

  constexpr auto log10(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::log10(val); }, mat);
  }

  /// @}

  /// @name Trigonometric functions.
  /// @{

  constexpr auto sin(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::sin(val); }, mat);
  }

  constexpr auto cos(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::cos(val); }, mat);
  }

  constexpr auto tan(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::tan(val); }, mat);
  }

  constexpr auto asin(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::asin(val); }, mat);
  }

  constexpr auto acos(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::acos(val); }, mat);
  }

  constexpr auto atan(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::atan(val); }, mat);
  }

  constexpr auto atan2(const is_matrix_view auto& y_mat,
                       const is_matrix_view auto& x_mat) noexcept {
    return apply([](const auto& y, const auto& x) { return math::atan2(y, x); },
                 y_mat, x_mat);
  }

  /// @}

  /// @name Hyperbolic functions.
  /// @{

  constexpr auto sinh(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::sinh(val); }, mat);
  }

  constexpr auto cosh(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::cosh(val); }, mat);
  }

  constexpr auto tanh(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::tanh(val); }, mat);
  }

  constexpr auto asinh(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::asinh(val); }, mat);
  }

  constexpr auto acosh(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::acosh(val); }, mat);
  }

  constexpr auto atanh(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::atanh(val); }, mat);
  }

  /// @}

} // namespace math

/// @}

/// @name Matrix-specific views.
/// @{

/// @brief Transpose the matrix @p mat.
constexpr auto transpose(const is_matrix_view auto& mat) noexcept {
  auto transpose_body{[mat = as_view(mat)](size_t row_index, size_t col_index) {
    return mat(col_index, row_index);
  }};
  return MatrixView(mat.num_cols(), mat.num_rows(), std::move(transpose_body));
}

/// @brief Multiply the matrices @p mat1 and @p mat2.
constexpr auto matmul(const is_matrix_view auto& mat1,
                      const is_matrix_view auto& mat2) noexcept {
  STORM_ASSERT_(mat1.num_cols() == mat2.num_rows() &&
                "The first matrix should have the same number of columns "
                "as the second matrix has rows.");
  auto matmul_body{[mat1 = as_view(mat1),
                    mat2 = as_view(mat2)](size_t row_index, size_t col_index) {
    const auto cross_size{mat1.num_cols()};
    auto val = mat1(row_index, 0) * mat2(0, col_index);
    for (size_t cross_index{1}; cross_index < cross_size; ++cross_index) {
      val += mat1(row_index, cross_index) * mat2(cross_index, col_index);
    }
    return val;
  }};
  return MatrixView(mat1.num_rows(), mat2.num_cols(), std::move(matmul_body));
}

/// @}

/// @name Slice views.
/// @{

/// @brief Slice the matrix @p mat rows from index @p from to index @p to
///   (not including) with a stride @p stride.
/// @{
template<is_matrix_view matrix_view>
constexpr auto slice_rows(matrix_view&& mat, size_t from, size_t to,
                          size_t stride = 1) {
  const size_t slice_num_rows{(to - from) / stride};
  STORM_ASSERT_((from < to && to <= mat.num_rows() && slice_num_rows != 0) &&
                "Invalid rows range.");
  auto slice_body{
      [=, mat = std::forward<matrix_view>(mat)](
          size_t slice_row_index, size_t col_index) -> decltype(auto) {
        STORM_ASSERT_(slice_row_index < slice_num_rows &&
                      "Row index is out of range.");
        const size_t row_index{from + slice_row_index * stride};
        return mat(row_index, col_index);
      }};
  return MatrixView(slice_num_rows, mat.num_cols(), std::move(slice_body));
}
constexpr auto slice_rows(is_matrix_view auto& mat, size_t from, size_t to,
                          size_t stride = 1) {
  return slice_rows(make_view(mat), from, to, stride);
}
/// @}

/// @brief Slice the matrix @p mat columns from index @p from to index @p to
///   (not including) with a stride @p stride.
/// @{
template<is_matrix_view matrix_view>
constexpr auto slice_cols(matrix_view&& mat, size_t from, size_t to,
                          size_t stride = 1) {
  const size_t slice_num_cols{(to - from) / stride};
  STORM_ASSERT_((from < to && to <= mat.num_cols() && slice_num_cols != 0) &&
                "Invalid columns range.");
  auto slice_body{
      [=, mat = std::forward<matrix_view>(mat)](
          size_t row_index, size_t slice_col_index) -> decltype(auto) {
        STORM_ASSERT_(slice_col_index < slice_num_cols &&
                      "Column index is out of range.");
        const size_t col_index{from + slice_col_index * stride};
        return mat(row_index, col_index);
      }};
  return MatrixView(slice_num_cols, mat.num_cols(), std::move(slice_body));
}
constexpr auto slice_cols(is_matrix_view auto& mat, size_t from, size_t to,
                          size_t stride = 1) {
  return slice_cols(as_view(mat), from, to, stride);
}
/// @}

/// @brief Select the matrix @p mat rows with @p row_indices.
/// @{
template<is_matrix_view matrix_view>
constexpr auto select_rows(matrix_view&& mat,
                           std::integral auto... row_indices) noexcept {
  STORM_ASSERT_((static_cast<size_t>(row_indices) < mat.num_rows()) && ... &&
                "Row indices are out of range.");
  auto select_body{
      [mat = std::forward<matrix_view>(mat),
       row_indices = std::array{static_cast<size_t>(row_indices)...}](
          size_t slice_row_index, size_t col_index) -> decltype(auto) {
        STORM_ASSERT_(slice_row_index < row_indices.size() &&
                      "Row index is out of range.");
        return mat(row_indices[slice_row_index], col_index);
      }};
  constexpr auto slice_num_rows{
      std::integral_constant<size_t, sizeof...(row_indices)>{}};
  return MatrixView(slice_num_rows, mat.num_cols(), std::move(select_body));
}
constexpr auto select_rows(is_matrix_view auto& mat,
                           std::integral auto... row_indices) noexcept {
  return select_rows(as_view(mat), row_indices...);
}
/// @}

/// @brief Select the matrix @p mat columns with @p col_index.
/// @{
template<is_matrix_view matrix_view>
constexpr auto select_cols(matrix_view&& mat,
                           std::integral auto... col_indices) noexcept {
  STORM_ASSERT_((static_cast<size_t>(col_indices) < mat.num_cols()) && ... &&
                "Columns indices are out of range.");
  auto select_body{
      [mat = std::forward<matrix_view>(mat),
       col_indices = std::array{static_cast<size_t>(col_indices)...}](
          size_t row_index, size_t slice_col_index) -> decltype(auto) {
        STORM_ASSERT_(slice_col_index < col_indices.size() &&
                      "Column index is out of range.");
        return mat(row_index, col_indices[slice_col_index]);
      }};
  constexpr auto slice_num_cols{
      std::integral_constant<size_t, sizeof...(col_indices)>{}};
  return MatrixView(mat.num_rows(), slice_num_cols, std::move(select_body));
}
constexpr auto select_cols(is_matrix_view auto& mat,
                           std::integral auto... col_indices) noexcept {
  return select_cols(as_view(mat), col_indices...);
}
/// @}

/// @}

/// @}

/// @brief Print a @p mat.
std::ostream& operator<<(std::ostream& out, const is_matrix_view auto& mat) {
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
