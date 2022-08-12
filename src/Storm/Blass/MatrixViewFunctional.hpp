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

#include <concepts>
#include <functional>
#include <tuple>
#include <type_traits>

#include <Storm/Base.hpp>

#include <Storm/Utils/Math.hpp>

#include <Storm/Blass/MatrixViewBase.hpp>

namespace Storm {

/// @brief Component-wise product of function to matrices view.
// clang-format off
template<std::copy_constructible Func, matrix_view... Matrices>
  requires std::is_object_v<Func> && (sizeof...(Matrices) >= 1) /*&&
           std::regular_invocable<Func, matrix_element_t<Matrices>...>*/
class MapMatrixView final :
    public MatrixViewInterface<MapMatrixView<Func, Matrices...>> {
  // clang-format on
private:

  STORM_NO_UNIQUE_ADDRESS_ Func func_;
  STORM_NO_UNIQUE_ADDRESS_ std::tuple<Matrices...> mats_;

public:

  /// @brief Construct a map view.
  constexpr MapMatrixView(Func func, Matrices... mats)
      : func_{std::move(func)}, mats_{std::move(mats)...} {
    STORM_ASSERT_( //
        std::apply(
            [](const auto& first_mat, const auto&... rest_mats) {
              return ((first_mat.shape() == rest_mats.shape()) && ...);
            },
            mats_) &&
        "Shapes of the matrix arguments are mismatched.");
  }

  /// @copydoc MatrixViewInterface::shape
  [[nodiscard]] constexpr auto shape() const noexcept {
    return std::get<0>(mats_).shape();
  }

  /// @copydoc MatrixViewInterface::operator()
  [[nodiscard]] constexpr decltype(auto)
  operator()(size_t row_index, size_t col_index) const noexcept {
    return std::apply(
        [&](const Matrices&... mats) {
          return func_(mats(row_index, col_index)...);
        },
        mats_);
  }

  /// @copydoc MatrixViewInterface::traverse_tree
  constexpr void traverse_tree(const auto& visitor) const {
    std::apply(
        [&](const Matrices&... mats) { //
          (mats.tranverse_tree(visitor), ...);
        },
        mats_);
    visitor(*this);
  }

  /// @copydoc MatrixViewInterface::transform_tree
  [[nodiscard]] constexpr auto transform_tree(const auto& transformer) const {
    return std::apply(
        [&](const Matrices&... mats) {
          return transformer(
              Storm::MapMatrixView(func_, mats.transform_tree(transformer)...));
        },
        mats_);
  }

}; // class MapMatrixView

template<class Func, class... Matrices>
MapMatrixView(Func, Matrices&&...)
    -> MapMatrixView<Func, forward_as_matrix_view_t<Matrices>...>;

namespace detail_ {
  // clang-format off
  template<class Func, class... Matrices>
  concept can_map_matrix_view_ = 
      requires { 
        MapMatrixView(std::declval<Func>(), std::declval<Matrices>()...); 
      };
  // clang-format on
} // namespace detail_

/// @brief Make a component-wise product of function @p func
///   to matrices @p mats.
// clang-format off
template<class Func, viewable_matrix... Matrices>
  requires detail_::can_map_matrix_view_<Func, Matrices...>
[[nodiscard]] constexpr auto map(Func&& func, Matrices&&... mats) {
  // clang-format on
  return MapMatrixView(std::forward<Func>(func),
                       std::forward<Matrices>(mats)...);
}

/// @brief Cast the @p mat elements to another type.
/// @tparam To Type, the elements would be cased to.
template<class To, viewable_matrix Matrix>
[[nodiscard]] constexpr auto matrix_cast(Matrix&& mat) {
  return MapMatrixView(
      []<class Elem>(Elem&& elem) noexcept {
        return static_cast<To>(std::forward<Elem>(elem));
      },
      std::forward<Matrix>(mat));
}

/// @brief "+" the matrix @p mat.
template<viewable_matrix Matrix>
[[nodiscard]] constexpr auto operator+(Matrix&& mat) {
  // Since operator "+" does nothing on artimetic types,
  // simply forward matrices with artimetic elements as matrix views.
  if constexpr (std::is_arithmetic_v<
                    std::remove_cvref_t<matrix_element_t<Matrix>>>) {
    return forward_as_matrix_view(std::forward<Matrix>(mat));
  } else {
    return MapMatrixView(
        []<class Elem>(Elem&& elem) noexcept {
          return +std::forward<Elem>(elem);
        },
        std::forward<Matrix>(mat));
  }
}

/// @brief Negate the matrix @p mat.
template<viewable_matrix Matrix>
[[nodiscard]] constexpr auto operator-(Matrix&& mat) {
  return MapMatrixView(std::negate{}, std::forward<Matrix>(mat));
}

/// @brief Multiply the matrix @p mat by a scalar @p scal.
/// @{
// clang-format off
template<viewable_matrix Matrix, std::copyable Scalar>
  requires (!matrix<Scalar>)
[[nodiscard]] constexpr auto operator*(Matrix&& mat, Scalar scal) {
  // clang-format on
  return MapMatrixView(
      [scal = std::move(scal)]<class Elem>(Elem&& elem) noexcept {
        return std::forward<Elem>(elem) * scal;
      },
      std::forward<Matrix>(mat));
}
// clang-format off
template<std::copyable Scalar, viewable_matrix Matrix>
  requires (!matrix<Scalar>)
[[nodiscard]] constexpr auto operator*(Scalar scal, Matrix&& mat) {
  // clang-format on
  return MapMatrixView(
      [scal = std::move(scal)]<class Elem>(Elem&& elem) noexcept {
        return scal * std::forward<Elem>(elem);
      },
      std::forward<Matrix>(mat));
}
/// @}

/// @brief Divide the matrix @p mat by a scalar @p scal.
// clang-format off
template<viewable_matrix Matrix, std::copyable Scalar>
  requires (!matrix<Scalar>)
[[nodiscard]] constexpr auto operator/(Matrix&& mat, Scalar scal) {
  // clang-format on
  return MapMatrixView(
      [scal = std::move(scal)]<class Elem>(Elem&& elem) noexcept {
        return std::forward<Elem>(elem) / scal;
      },
      std::forward<Matrix>(mat));
}

/// @brief Divide the scalar @p scal by a matrix @p mat.
// clang-format off
template<std::copyable Scalar, viewable_matrix Matrix>
  requires (!matrix<Scalar>)
[[nodiscard]] constexpr auto operator/(Scalar scal, Matrix&& mat) {
  // clang-format on
  return MapMatrixView(
      [scal = std::move(scal)]<class Elem>(Elem&& elem) noexcept {
        return scal / std::forward<Elem>(elem);
      },
      std::forward<Matrix>(mat));
}

/// @brief Add the matrices @p mat1 and @p mat2.
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
[[nodiscard]] constexpr auto operator+(Matrix1&& mat1, Matrix2&& mat2) {
  return MapMatrixView(std::plus{}, //
                       std::forward<Matrix1>(mat1),
                       std::forward<Matrix2>(mat2));
}

/// @brief Subtract the matrices @p mat1 and @p mat2.
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
[[nodiscard]] constexpr auto operator-(Matrix1&& mat1, Matrix2&& mat2) {
  return MapMatrixView(std::minus{}, //
                       std::forward<Matrix1>(mat1),
                       std::forward<Matrix2>(mat2));
}

/// @brief Component-wise multiply the matrices @p mat1 and @p mat2.
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
[[nodiscard]] constexpr auto operator*(Matrix1&& mat1, Matrix2&& mat2) {
  return MapMatrixView(std::multiplies{}, //
                       std::forward<Matrix1>(mat1),
                       std::forward<Matrix2>(mat2));
}

/// @brief Component-wise divide the matrices @p mat1 and @p mat2.
template<viewable_matrix Matrix1, viewable_matrix Matrix2>
[[nodiscard]] constexpr auto operator/(Matrix1&& mat1, Matrix2&& mat2) {
  return MapMatrixView(std::divides{}, //
                       std::forward<Matrix1>(mat1),
                       std::forward<Matrix2>(mat2));
}

#define MAKE_UNARY_MATRIX_FUNC_(func)               \
  template<viewable_matrix Matrix>                  \
  [[nodiscard]] constexpr auto func(Matrix&& mat) { \
    return MapMatrixView(                           \
        []<class Elem>(Elem&& elem) noexcept {      \
          return func(std::forward<Elem>(elem));    \
        },                                          \
        std::forward<Matrix>(mat));                 \
  }

// clang-format off
#define MAKE_BINARY_SCALAR_MATRIX_FUNC_(func)                        \
  template<std::copyable Scalar, viewable_matrix Matrix>             \
    requires(!matrix<Scalar>)                                        \
  [[nodiscard]] constexpr auto func(Scalar scal, Matrix&& mat) {     \
    return MapMatrixView(                                            \
        [scal = std::move(scal)]<class Elem>(Elem&& elem) noexcept { \
          return func(scal, std::forward<Elem>(elem));               \
        },                                                           \
        std::forward<Matrix>(mat));                                  \
  }
// clang-format on

// clang-format off
#define MAKE_BINARY_MATRIX_SCALAR_FUNC_(func)                        \
  template<viewable_matrix Matrix, std::copyable Scalar>             \
    requires(!matrix<Scalar>)                                        \
  [[nodiscard]] constexpr auto func(Matrix&& mat, Scalar scal) {     \
    return MapMatrixView(                                            \
        [scal = std::move(scal)]<class Elem>(Elem&& elem) noexcept { \
          return func(std::forward<Elem>(elem), scal);               \
        },                                                           \
        std::forward<Matrix>(mat));                                  \
  }
// clang-format on

#define MAKE_BINARY_MATRIX_MATRIX_FUNC_(func)                                  \
  template<viewable_matrix Matrix1, viewable_matrix Matrix2>                   \
  [[nodiscard]] constexpr auto func(Matrix1&& mat1, Matrix2&& mat2) {          \
    return MapMatrixView(                                                      \
        []<class Elem1, class Elem2>(Elem1&& elem1, Elem2&& elem2) noexcept {  \
          return func(std::forward<Elem1>(elem1), std::forward<Elem2>(elem2)); \
        },                                                                     \
        std::forward<Matrix1>(mat1), std::forward<Matrix2>(mat2));             \
  }

namespace math {

  MAKE_UNARY_MATRIX_FUNC_(abs)

  /// @name Power functions.
  /// @{

  MAKE_BINARY_MATRIX_SCALAR_FUNC_(pow)

  MAKE_BINARY_SCALAR_MATRIX_FUNC_(pow)

  MAKE_BINARY_MATRIX_MATRIX_FUNC_(pow)

  MAKE_UNARY_MATRIX_FUNC_(sqrt)

  MAKE_UNARY_MATRIX_FUNC_(cbrt)

  // clang-format off
  template<viewable_matrix... Matrices>
    requires(detail_::in_range_(sizeof...(Matrices), 2, 3))
  [[nodiscard]] constexpr auto hypot(Matrices&&... mats) {
    // clang-format on
    return MapMatrixView(
        []<class... Elems>(Elems&&... elems) noexcept {
          return hypot(std::forward<Elems>(elems)...);
        },
        std::forward<Matrices>(mats)...);
  }

  /// @} // Power functions.

  /// @name Exponential functions.
  /// @{

  MAKE_UNARY_MATRIX_FUNC_(exp)

  MAKE_UNARY_MATRIX_FUNC_(exp2)

  MAKE_UNARY_MATRIX_FUNC_(log)

  MAKE_UNARY_MATRIX_FUNC_(log2)

  MAKE_UNARY_MATRIX_FUNC_(log10)

  /// @} // Exponential functions.

  /// @name Trigonometric functions.
  /// @{

  MAKE_UNARY_MATRIX_FUNC_(sin)

  MAKE_UNARY_MATRIX_FUNC_(cos)

  MAKE_UNARY_MATRIX_FUNC_(tan)

  MAKE_UNARY_MATRIX_FUNC_(asin)

  MAKE_UNARY_MATRIX_FUNC_(acos)

  MAKE_UNARY_MATRIX_FUNC_(atan)

  MAKE_BINARY_MATRIX_MATRIX_FUNC_(atan2)

  /// @} // Trigonometric functions.

  /// @name Hyperbolic functions.
  /// @{

  MAKE_UNARY_MATRIX_FUNC_(sinh)

  MAKE_UNARY_MATRIX_FUNC_(cosh)

  MAKE_UNARY_MATRIX_FUNC_(tanh)

  MAKE_UNARY_MATRIX_FUNC_(asinh)

  MAKE_UNARY_MATRIX_FUNC_(acosh)

  MAKE_UNARY_MATRIX_FUNC_(atanh)

  /// @} // Hyperbolic functions.

} // namespace math

#undef MAKE_UNARY_MATRIX_FUNC_
#undef MAKE_BINARY_MATRIX_SCALAR_FUNC_
#undef MAKE_BINARY_SCALAR_MATRIX_FUNC_
#undef MAKE_BINARY_MATRIX_MATRIX_FUNC_

/// @brief Matrix transpose view.
template<matrix_view Matrix>
class MatrixTransposeView final :
    public MatrixViewInterface<MatrixTransposeView<Matrix>> {
private:

  STORM_NO_UNIQUE_ADDRESS_ Matrix mat_;

public:

  /// @brief Construct a matrix transpose view.
  constexpr explicit MatrixTransposeView(Matrix mat) : mat_{std::move(mat)} {}

  /// @copydoc MatrixViewInterface::shape
  [[nodiscard]] constexpr auto shape() const noexcept {
    return MatrixShape(num_cols(mat_), num_rows(mat_));
  }

  /// @copydoc MatrixViewInterface::operator()
  /// @{
  [[nodiscard]] constexpr decltype(auto) //
  operator()(size_t row_index, size_t col_index) noexcept {
    return mat_(col_index, row_index);
  }
  [[nodiscard]] constexpr decltype(auto)
  operator()(size_t row_index, size_t col_index) const noexcept {
    return mat_(col_index, row_index);
  }
  /// @}

  /// @copydoc MatrixViewInterface::traverse_tree
  constexpr void traverse_tree(const auto& visitor) const {
    mat_.tranverse_tree(visitor);
    visitor(*this);
  }

  /// @copydoc MatrixViewInterface::transform_tree
  [[nodiscard]] constexpr auto transform_tree(const auto& transformer) const {
    return transformer(
        Storm::MatrixTransposeView(mat_.transform_tree(transformer)));
  }

}; // MatrixTransposeView

template<class Matrix>
MatrixTransposeView(Matrix&&)
    -> MatrixTransposeView<forward_as_matrix_view_t<Matrix>>;

/// @brief Transpose the matrix @p mat.
template<viewable_matrix Matrix>
[[nodiscard]] constexpr auto transpose(Matrix&& mat) {
  return MatrixTransposeView(std::forward<Matrix>(mat));
}

/// @brief Matrix product view.
template<matrix_view Matrix1, matrix_view Matrix2>
class MatrixProductView final :
    public MatrixViewInterface<MatrixProductView<Matrix1, Matrix2>> {
private:

  STORM_NO_UNIQUE_ADDRESS_ Matrix1 mat1_;
  STORM_NO_UNIQUE_ADDRESS_ Matrix2 mat2_;

public:

  /// @brief Construct a matrix product view.
  constexpr explicit MatrixProductView(Matrix1 mat1, Matrix2 mat2)
      : mat1_{std::move(mat1)}, mat2_{std::move(mat2)} {
    STORM_ASSERT_(num_cols(mat1_) == num_rows(mat2) &&
                  "Shapes of the matrix product arguments are invalid.");
  }

  /// @copydoc MatrixViewInterface::shape
  [[nodiscard]] constexpr auto shape() const noexcept {
    return MatrixShape(num_rows(mat1_), num_cols(mat2_));
  }

  /// @copydoc MatrixViewInterface::operator()
  [[nodiscard]] constexpr auto operator()(size_t row_index,
                                          size_t col_index) const noexcept {
    // This is a default very slow generic implementation.
    auto val{mat1_(row_index, 0) * mat2_(0, col_index)};
    const auto cross_size = static_cast<size_t>(num_cols(mat1_));
    for (size_t cross_index{1}; cross_index < cross_size; ++cross_index) {
      val += mat1_(row_index, cross_index) * mat2_(cross_index, col_index);
    }
    return val;
  }

  /// @copydoc MatrixViewInterface::traverse_tree
  constexpr void traverse_tree(const auto& visitor) const {
    mat1_.tranverse_tree(visitor), mat2_.tranverse_tree(visitor);
    visitor(*this);
  }

  /// @copydoc MatrixViewInterface::transform_tree
  [[nodiscard]] constexpr auto transform_tree(const auto& transformer) const {
    return transformer( //
        Storm::MatrixProductView(mat1_.transform_tree(transformer),
                                 mat2_.transform_tree(transformer)));
  }

}; // class MatrixProductView

template<class Matrix1, class Matrix2>
MatrixProductView(Matrix1&&, Matrix2&&)
    -> MatrixProductView<forward_as_matrix_view_t<Matrix1>,
                         forward_as_matrix_view_t<Matrix2>>;

/// @brief Multiply the matrices @p mat1 and @p mat2.
template<class Matrix1, class Matrix2>
[[nodiscard]] constexpr auto matmul(Matrix1&& mat1, Matrix2&& mat2) {
  return MatrixProductView(std::forward<Matrix1>(mat1),
                           std::forward<Matrix2>(mat2));
}

} // namespace Storm
