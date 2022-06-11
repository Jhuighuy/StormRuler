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
#include <tuple>
#include <type_traits>
#include <utility>

#include <stormBase.hxx>
#include <stormSolvers/MatrixBase.hxx>
#include <stormUtils/Math.hxx>

namespace Storm {

/// @name Matrix views.
/// @{

/// @name Cast to view.
/// @{

/// ----------------------------------------------------------------- ///
/// @brief Matrix as a matrix view wrapper.
/// ----------------------------------------------------------------- ///
template<is_maybe_const_matrix Matrix>
class MatrixAsView final {
private:

  Matrix& mat_;

public:

  /// @brief Construct a matrix as view wrapper.
  constexpr MatrixAsView(Matrix& mat) noexcept : mat_{mat} {}

  /// @brief Number of the matrix rows.
  constexpr auto num_rows() const noexcept {
    return mat_.num_rows();
  }

  /// @brief Number of the matrix columns.
  constexpr auto num_cols() const noexcept {
    return mat_.num_cols();
  }

  /// @brief Get the coefficient at @p row_index and @p col_index.
  /// @{
  constexpr auto operator()(size_t row_index, size_t col_index) noexcept
      -> decltype(auto) {
    return mat_(row_index, col_index);
  }
  constexpr auto operator()(size_t row_index, size_t col_index) const noexcept
      -> decltype(auto) {
    return std::as_const(mat_)(row_index, col_index);
  }
  /// @}

}; // class MatrixAsView

template<class Matrix>
MatrixAsView(Matrix&) -> MatrixAsView<Matrix>;

template<class Matrix>
struct is_matrix_view_t<MatrixAsView<Matrix>> : std::true_type {};

/// @brief Wrap the matrix @p mat in view a matrix view.
constexpr auto forward_as_view(is_matrix_ref auto&& mat) noexcept {
  return MatrixAsView{std::forward<decltype(mat)>(mat)};
}

/// @brief Forward the matrix view @p mat as an r-value reference.
constexpr auto&& 
forward_as_view(is_strictly_matrix_view auto&& mat) noexcept {
  return std::forward<decltype(mat)>(mat);
}

/// @brief Copy the matrix view @p mat.
constexpr auto
forward_as_view(const is_strictly_matrix_view auto& mat) noexcept {
  return mat;
}

/// @}

/// @name Matrix expression views.
/// @{

/// ----------------------------------------------------------------- ///
/// @brief Matrix expression view.
/// ----------------------------------------------------------------- ///
// clang-format off
template<std::convertible_to<size_t> RowsSize,
         std::convertible_to<size_t> ColsSize, //
         class ExprFunc, is_matrix_view... ExprArgs>
  requires(std::invocable<ExprFunc, RowsSize, ColsSize, ExprArgs&...> &&
           std::invocable<ExprFunc, RowsSize, ColsSize, 
                          std::add_const_t<ExprArgs>&...>)
class MatrixView {
private:

  [[no_unique_address]] RowsSize num_rows_;
  [[no_unique_address]] ColsSize num_cols_;
  [[no_unique_address]] ExprFunc expr_func_;
  [[no_unique_address]] std::tuple<ExprArgs...> expr_args_;

public:

  /// @brief Construct a matrix extression view.
  constexpr MatrixView(RowsSize num_rows, ColsSize num_cols,
                       ExprFunc&& expr_func, ExprArgs&&... expr_args)
      : num_rows_{num_rows}, num_cols_{num_cols},
        expr_func_{std::forward<ExprFunc>(expr_func)},
        expr_args_{std::forward<ExprArgs>(expr_args)...} {}

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
  constexpr auto operator()(size_t row_index, size_t col_index) noexcept 
      -> decltype(auto) {
    STORM_ASSERT_(row_index < num_rows_ && col_index < num_cols_ &&
                  "Indices are out of range.");
    return std::apply(
        [&expr_func_, row_index, col_index](auto&... expr_args) -> decltype(auto) {
          return expr_func_(row_index, col_index, expr_args...);
        },
        expr_args_);
  }
  constexpr auto operator()(size_t row_index, size_t col_index) const noexcept 
      -> decltype(auto) {
    STORM_ASSERT_(row_index < num_rows_ && col_index < num_cols_ &&
                  "Indices are out of range.");
    return std::apply(
        [&expr_func_, row_index, col_index](const auto&... expr_args) -> decltype(auto) {
          return expr_func_(row_index, col_index, expr_args...);
        },
        expr_args_);
  }
  /// @}

}; // class MatrixView
// clang-format on

template<class RowsSize, class ColsSize, class ExprFunc, class... ExprArgs>
MatrixView(RowsSize, ColsSize, ExprFunc&&, ExprArgs&&...)
    -> MatrixView<RowsSize, ColsSize, ExprFunc, ExprArgs...>;

template<class RowsSize, class ColsSize, class ExprFunc, class... ExprArgs>
struct is_matrix_view_t<MatrixView<RowsSize, ColsSize, ExprFunc, ExprArgs...>> :
    std::true_type {};

/// @name General views.
/// @{

/// @brief Component-wise apply a @p func
///   to the matrix arguments @p mat1, @p mats.
constexpr auto apply(auto func, const is_matrix_view auto& mat1,
                     const is_matrix_view auto&... mats) noexcept {
  STORM_ASSERT_(((mat1.num_rows() == mats.num_rows()) && ...) &&
                ((mat1.num_cols() == mats.num_cols()) && ...) &&
                "Shapes of the matrix arguments should be the same.");
  return MatrixView(
      mat1.num_rows(), mat1.num_cols(),
      [func](size_t row_index, size_t col_index,
             const is_matrix_view auto& mat1,
             const is_matrix_view auto&... mats) {
        return func(mat1(row_index, col_index), mats(row_index, col_index)...);
      },
      forward_as_view(mat1), forward_as_view(mats)...);
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
constexpr auto& operator*=(is_rw_matrix_view_ref auto&& mat1,
                           const auto& val2) {
  return mat1 <<= val2 * mat1;
}
/// @}

/// @brief Divide the matrix @p mat by a scalar @p scal.
/// @{
constexpr auto operator/(const is_matrix_view auto& mat, auto scal) noexcept {
  return apply([scal](const auto& val) { return val / scal; }, mat);
}
constexpr auto& operator/=(is_rw_matrix_view_ref auto&& mat1,
                           const auto& val2) {
  return mat1 <<= mat1 / val2;
}
/// @}

/// @brief Add the matrices @p mat1 and @p mat2.
/// @{
constexpr auto operator+(const is_matrix_view auto& mat1,
                         const is_matrix_view auto& mat2) noexcept {
  return apply([](const auto& val1, const auto& val2) { return val1 + val2; },
               mat1, mat2);
}
constexpr auto& operator+=(is_rw_matrix_view_ref auto&& mat1,
                           const is_matrix_view auto& mat2) {
  return mat1 <<= mat1 + mat2;
}
/// @}

/// @brief Subtract the matrices @p mat1 and @p mat2.
/// @{
constexpr auto operator-(const is_matrix_view auto& mat1,
                         const is_matrix_view auto& mat2) noexcept {
  return apply([](const auto& val1, const auto& val2) { return val1 - val2; },
               mat1, mat2);
}
constexpr auto& operator-=(is_rw_matrix_view_ref auto&& mat1,
                           const is_matrix_view auto& mat2) {
  return mat1 <<= mat1 - mat2;
}
/// @}

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

  /// @brief Component-wise @c abs of the matrix @p mat.
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

  /// @brief Component-wise @c sqrt of the matrix @p mat.
  constexpr auto sqrt(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::sqrt(val); }, mat);
  }

  /// @brief Component-wise @c cbrt of the matrix @p mat.
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

  /// @brief Component-wise @c exp of the matrix @p mat.
  constexpr auto exp(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::exp(val); }, mat);
  }

  /// @brief Component-wise @c exp2 of the matrix @p mat.
  constexpr auto exp2(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::exp2(val); }, mat);
  }

  /// @brief Component-wise @c log of the matrix @p mat.
  constexpr auto log(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::log(val); }, mat);
  }

  /// @brief Component-wise @c log2 of the matrix @p mat.
  constexpr auto log2(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::log2(val); }, mat);
  }

  /// @brief Component-wise @c log10 of the matrix @p mat.
  constexpr auto log10(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::log10(val); }, mat);
  }

  /// @}

  /// @name Trigonometric functions.
  /// @{

  /// @brief Component-wise @c sin of the matrix @p mat.
  constexpr auto sin(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::sin(val); }, mat);
  }

  /// @brief Component-wise @c cos of the matrix @p mat.
  constexpr auto cos(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::cos(val); }, mat);
  }

  /// @brief Component-wise @c tan of the matrix @p mat.
  constexpr auto tan(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::tan(val); }, mat);
  }

  /// @brief Component-wise @c asin of the matrix @p mat.
  constexpr auto asin(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::asin(val); }, mat);
  }

  /// @brief Component-wise @c acos of the matrix @p mat.
  constexpr auto acos(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::acos(val); }, mat);
  }

  /// @brief Component-wise @c atan of the matrix @p mat.
  constexpr auto atan(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::atan(val); }, mat);
  }

  /// @brief Component-wise @c atan2 of the matriсes @p y_mat and @p x_mat.
  constexpr auto atan2(const is_matrix_view auto& y_mat,
                       const is_matrix_view auto& x_mat) noexcept {
    return apply([](const auto& y, const auto& x) { return math::atan2(y, x); },
                 y_mat, x_mat);
  }

  /// @}

  /// @name Hyperbolic functions.
  /// @{

  /// @brief Component-wise @p sinh of the matrix @p mat.
  constexpr auto sinh(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::sinh(val); }, mat);
  }

  /// @brief Component-wise @c cosh of the matrix @p mat.
  constexpr auto cosh(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::cosh(val); }, mat);
  }

  /// @brief Component-wise @c tanh of the matrix @p mat.
  constexpr auto tanh(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::tanh(val); }, mat);
  }

  /// @brief Component-wise @c asinh of the matrix @p mat.
  constexpr auto asinh(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::asinh(val); }, mat);
  }

  /// @brief Component-wise @c acosh of the matrix @p mat.
  constexpr auto acosh(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::acosh(val); }, mat);
  }

  /// @brief Component-wise @c atanh of the matrix @p mat.
  constexpr auto atanh(const is_matrix_view auto& mat) noexcept {
    return apply([](const auto& val) { return math::atanh(val); }, mat);
  }

  /// @}

} // namespace math

/// @}

/// @name Matrix functions views.
/// @{

/// @brief Transpose the matrix @p mat.
constexpr auto transpose(const is_matrix_view auto& mat) noexcept {
  return MatrixView(
      mat.num_cols(), mat.num_rows(),
      [](size_t row_index, size_t col_index, const is_matrix_view auto& mat) {
        return mat(col_index, row_index);
      },
      forward_as_view(mat));
}

/// @brief Multiply the matrices @p mat1 and @p mat2.
constexpr auto matmul(const is_matrix_view auto& mat1,
                      const is_matrix_view auto& mat2) noexcept {
  STORM_ASSERT_(mat1.num_cols() == mat2.num_rows() &&
                "The first matrix should have the same number of columns "
                "as the second matrix has rows.");
  return MatrixView(
      mat1.num_rows(), mat2.num_cols(),
      [](size_t row_index, size_t col_index, //
         const is_matrix_view auto& mat1, const is_matrix_view auto& mat2) {
        const auto cross_size{mat1.num_cols()};
        auto val = mat1(row_index, 0) * mat2(0, col_index);
        for (size_t cross_index{1}; cross_index < cross_size; ++cross_index) {
          val += mat1(row_index, cross_index) * mat2(cross_index, col_index);
        }
        return val;
      },
      forward_as_view(mat1), forward_as_view(mat2));
}

constexpr auto diag(const is_matrix_view auto& mat1) noexcept;

constexpr auto lower_triangle(const is_matrix_view auto& mat1) noexcept;

constexpr auto upper_triangle(const is_matrix_view auto& mat1) noexcept;

/// @}

/// @name Slice views.
/// @{

/// @brief Slice the matrix @p mat rows from index @p from to index @p to
///   (not including) with a stride @p stride.
constexpr auto slice_rows(is_matrix_view_ref auto&& mat, //
                          size_t from, size_t to, size_t stride = 1) noexcept {
  STORM_ASSERT_((from < to && to <= mat.num_rows()) && "Invalid row range.");
  const size_t slice_num_rows{(to - from) / stride};
  return MatrixView(
      slice_num_rows, mat.num_cols(),
      [=](size_t slice_row_index, size_t col_index,
          is_matrix_view_ref auto&& mat) -> decltype(auto) {
        const size_t row_index{from + slice_row_index * stride};
        return mat(row_index, col_index);
      },
      forward_as_view(mat));
}

/// @brief Slice the matrix @p mat columns from index @p from to index @p to
///   (not including) with a stride @p stride.
constexpr auto slice_cols(is_matrix_view_ref auto&& mat, //
                          size_t from, size_t to, size_t stride = 1) noexcept {
  STORM_ASSERT_((from < to && to <= mat.num_cols()) && "Invalid column range.");
  const size_t slice_num_cols{(to - from) / stride};
  return MatrixView(
      slice_num_cols, mat.num_cols(),
      [=](size_t row_index, size_t slice_col_index,
          is_matrix_view_ref auto&& mat) -> decltype(auto) {
        const size_t col_index{from + slice_col_index * stride};
        return mat(row_index, col_index);
      },
      forward_as_view(mat));
}

/// @brief Select the matrix @p mat rows with @p row_indices.
constexpr auto select_rows(is_matrix_view_ref auto&& mat,
                           std::integral auto... row_indices) noexcept {
  STORM_ASSERT_((static_cast<size_t>(row_indices) < mat.num_rows()) && ... &&
                "Row indices are out of range.");
  constexpr size_t_constant<sizeof...(row_indices)> slice_num_rows{};
  return MatrixView(
      slice_num_rows, mat.num_cols(),
      [row_indices = std::array{static_cast<size_t>(row_indices)...}](
          size_t slice_row_index, size_t col_index,
          is_matrix_view_ref auto&& mat) -> decltype(auto) {
        return mat(row_indices[slice_row_index], col_index);
      },
      forward_as_view(mat));
}

/// @brief Select the matrix @p mat columns with @p col_index.
constexpr auto select_cols(is_matrix_view_ref auto&& mat,
                           std::integral auto... col_indices) noexcept {
  STORM_ASSERT_((static_cast<size_t>(col_indices) < mat.num_cols()) && ... &&
                "Columns indices are out of range.");
  constexpr size_t_constant<sizeof...(col_indices)> slice_num_cols{};
  return MatrixView(
      mat.num_rows(), slice_num_cols,
      [col_indices = std::array{static_cast<size_t>(col_indices)...}](
          size_t row_index, size_t slice_col_index,
          is_matrix_view_ref auto&& mat) -> decltype(auto) {
        return mat(row_index, col_indices[slice_col_index]);
      },
      forward_as_view(mat));
}

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
