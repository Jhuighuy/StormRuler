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

#include <type_traits>

#include <stormBase.hxx>

namespace Storm {

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Base class for all matrix views.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Derived>
class BaseMatrixView {
public:

  /// @brief Get the matrix shape.
  constexpr auto shape() const noexcept {
    return static_cast<const Derived&>(*this).shape();
  }

  /// @brief Number of rows.
  constexpr auto num_rows() const noexcept {
    return shape().first;
  }

  /// @brief Number of columns.
  constexpr auto num_cols() const noexcept {
    return shape().second;
  }

  /// @brief Get the matrix coefficient at @p row_index and @p col_index.
  constexpr auto operator()(auto row_index, auto col_index) const noexcept {
    return static_cast<const Derived&>(*this)(row_index, col_index);
  }

}; // class BaseMatrixView

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Matrix view.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Shape, class Indexable>
class MatrixView : public BaseMatrixView<MatrixView<Shape, Indexable>> {
private:

  [[no_unique_address]] Shape shape_;
  [[no_unique_address]] Indexable indexable_;

public:

  /// @brief Construct a matrix view.
  constexpr MatrixView(Shape shape, Indexable indexable)
      : shape_{shape}, indexable_{indexable} {}

  /// @copydoc BaseMatrixView::shape
  constexpr auto shape() const noexcept {
    return shape_;
  }

  /// @copydoc BaseMatrixView::operator()
  constexpr auto operator()(auto row_index, auto col_index) const noexcept {
    return indexable_(row_index, col_index);
  }

}; // class MatrixView

template<class ShapeRef, class IndexableRef>
MatrixView(ShapeRef, IndexableRef)
    -> MatrixView<std::decay_t<ShapeRef>, std::decay_t<IndexableRef>>;

/// @brief Transpose the matrix @p mat.
template<class T>
constexpr auto transpose(const BaseMatrixView<T>& mat) noexcept {
  const std::pair transposed_shape{mat.shape().second, mat.shape().first};
  return MatrixView(transposed_shape, [&](auto row_index, auto col_index) {
    return mat(col_index, row_index);
  });
}

/// @brief Apply a @p func to the matrix arguments @p mat1, @p mats.
template<class T1, class... TN>
constexpr auto apply(auto func, const BaseMatrixView<T1>& mat1,
                     const BaseMatrixView<TN>&... mats) noexcept {
  STORM_ASSERT_((mat1.shape() == mats.shape()) && ... &&
                "Shapes of the matrix arguments should be the same");
  return MatrixView(mat1.shape(), [&, func](auto row_index, auto col_index) {
    return func(mat1(row_index, col_index), mats(row_index, col_index)...);
  });
}

/// @brief Multiply the matrix @p mat by a scalar @p scal.
/// @{
template<class T>
constexpr auto operator*(const auto& scal,
                         const BaseMatrixView<T>& mat) noexcept {
  return apply([scal](const auto& val) { return scal * val; }, mat);
}
template<class T>
constexpr auto operator*(const BaseMatrixView<T>& mat,
                         const auto& scal) noexcept {
  return apply([scal](const auto& val) { return val * scal; }, mat);
}
/// @}

/// @brief Divide the matrix @p mat by a scalar @p scal.
template<class T>
constexpr auto operator/(const BaseMatrixView<T>& mat, auto&& scal) noexcept {
  return apply([scal](const auto& val) { return scal * val; }, mat);
}

/// @brief Add the matrices @p mat1 and @p mat2.
template<class T1, class T2>
constexpr auto operator+(const BaseMatrixView<T1>& mat1,
                         const BaseMatrixView<T2>& mat2) noexcept {
  return apply([](const auto& val1, const auto& val2) { return val1 + val2; },
               mat1, mat2);
}

/// @brief Subtract the matrices @p mat1 and @p mat2.
template<class T1, class T2>
constexpr auto operator-(const BaseMatrixView<T1>& mat1,
                         const BaseMatrixView<T2>& mat2) noexcept {
  return apply([](const auto& val1, const auto& val2) { return val1 - val2; },
               mat1, mat2);
}

/// @brief Component-wise multiply the matrices @p mat1 and @p mat2.
template<class T1, class T2>
constexpr auto operator*(const BaseMatrixView<T1>& mat1,
                         const BaseMatrixView<T2>& mat2) noexcept {
  return apply([](const auto& val1, const auto& val2) { return val1 * val2; },
               mat1, mat2);
}

/// @brief Component-wise divide the matrices @p mat1 and @p mat2.
template<class T1, class T2>
constexpr auto operator/(const BaseMatrixView<T1>& mat1,
                         const BaseMatrixView<T2>& mat2) noexcept {
  return apply([](const auto& val1, const auto& val2) { return val1 / val2; },
               mat1, mat2);
}

/// @brief Multiply the matrices @p mat1 and @p mat2.
template<class T1, class T2>
constexpr auto matmul(const BaseMatrixView<T1>& mat1,
                      const BaseMatrixView<T2>& mat2) noexcept {
  STORM_ASSERT_(mat1.shape().second == mat2.shape().first &&
                "mat1 should have the same number or columns as mat2 has rows");
  const auto reduction_size{mat1.shape().second};
  const std::pair product_shape{mat1.shape().first, mat2.shape().second};
  return MatrixView(product_shape,
                    [&, reduction_size](size_t row_index, size_t col_index) {
                      auto val = mat1(row_index, 0) * mat2(0, col_index);
                      for (size_t iz{1}; iz < reduction_size; ++iz) {
                        val += mat1(row_index, iz) * mat2(iz, col_index);
                      }
                      return val;
                    });
}

/// @brief Print a @p mat.
template<class T>
std::ostream& operator<<(std::ostream& out, const BaseMatrixView<T>& mat) {
  for (size_t row_index{0}; row_index < mat.num_rows(); ++row_index) {
    for (size_t col_index{0}; col_index < mat.num_cols(); ++col_index) {
      out << mat(row_index, col_index) << ' ';
    }
    out << std::endl;
  }
  return out;
}

} // namespace Storm
