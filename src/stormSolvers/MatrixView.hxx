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
    return static_cast<Derived const&>(*this).shape();
  }

  /// @brief Number of rows.
  constexpr size_t NumRows() const noexcept {
    return shape().first;
  }

  /// @brief Number of columns.
  constexpr size_t NumCols() const noexcept {
    return shape().second;
  }

  /// @brief Get the matrix coefficient at @p rowIndex and @p colIndex.
  constexpr auto operator()(size_t rowIndex, size_t colIndex) const noexcept {
    return static_cast<Derived const&>(*this)(rowIndex, colIndex);
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
  constexpr auto operator()(size_t rowIndex, size_t colIndex) const noexcept {
    return indexable_(rowIndex, colIndex);
  }

}; // class MatrixView

template<class Shape, class Indexable>
MatrixView(Shape, Indexable)
    -> MatrixView<std::decay_t<Shape>, std::decay_t<Indexable>>;

/// @{
template<class T1, class... TN>
constexpr auto copy_shape_(BaseMatrixView<T1> const& mat1,
                           BaseMatrixView<TN> const&... mats) {
  StormAssert((mat1.shape() == mats.shape()) && ... &&
              "Shapes of the matrix arguments should be the same");
  return mat1.shape();
}
/// @}

template<class T1, class T2>
constexpr auto merge_shape_(BaseMatrixView<T1> const& mat1,
                            BaseMatrixView<T2> const& mat2) {
  return std::pair(std::pair(mat1.shape().first, mat2.shape().second),
                   mat1.shape().second);
}

/// @brief Add the matrices @p mat1 and @p mat2.
template<class T1, class T2>
constexpr auto operator+(BaseMatrixView<T1> const& mat1,
                         BaseMatrixView<T2> const& mat2) noexcept {
  return MatrixView(
      copy_shape_(mat1, mat2), [&](size_t rowIndex, size_t colIndex) {
        return mat1(rowIndex, colIndex) + mat2(rowIndex, colIndex);
      });
}

/// @brief Subtract the matrices @p mat1 and @p mat2.
template<class T1, class T2>
constexpr auto operator-(BaseMatrixView<T1> const& mat1,
                         BaseMatrixView<T2> const& mat2) noexcept {
  return MatrixView(
      copy_shape_(mat1, mat2), [&](size_t rowIndex, size_t colIndex) {
        return mat1(rowIndex, colIndex) - mat2(rowIndex, colIndex);
      });
}

template<class V1, class T2>
constexpr auto operator*(V1 const& val1,
                         BaseMatrixView<T2> const& mat2) noexcept {
  return MatrixView(copy_shape_(mat2), [&](size_t rowIndex, size_t colIndex) {
    return val1 * mat2(rowIndex, colIndex);
  });
}

template<class T1, class V2>
constexpr auto operator/(BaseMatrixView<T1> const& mat1,
                         V2 const& val2) noexcept {
  return MatrixView(copy_shape_(mat1), [&](size_t rowIndex, size_t colIndex) {
    return mat1(rowIndex, colIndex) / val2;
  });
}

/// @brief Component-wise multiply the matrices @p mat1 and @p mat2.
template<class T1, class T2>
constexpr auto operator*(BaseMatrixView<T1> const& mat1,
                         BaseMatrixView<T2> const& mat2) noexcept {
  return MatrixView(
      copy_shape_(mat1, mat2), [&](size_t rowIndex, size_t colIndex) {
        return mat1(rowIndex, colIndex) * mat2(rowIndex, colIndex);
      });
}

/// @brief Component-wise divide the matrices @p mat1 and @p mat2.
template<class T1, class T2>
constexpr auto operator/(BaseMatrixView<T1> const& mat1,
                         BaseMatrixView<T2> const& mat2) noexcept {
  return MatrixView(
      copy_shape_(mat1, mat2), [&](size_t rowIndex, size_t colIndex) {
        return mat1(rowIndex, colIndex) / mat2(rowIndex, colIndex);
      });
}

/// @brief Transpose the matrix @p mat.
template<class T>
constexpr auto transpose(BaseMatrixView<T> const& mat) noexcept {
  return MatrixView(std::pair(mat.shape().second, mat.shape().first),
                    [&](size_t rowIndex, size_t colIndex) {
                      return mat(colIndex, rowIndex);
                    });
}

/// @brief Multiply the matrices @p mat1 and @p mat2.
template<class T1, class T2>
constexpr auto matmul(BaseMatrixView<T1> const& mat1,
                      BaseMatrixView<T2> const& mat2) noexcept {
  auto const [product_shape, cross_size] = merge_shape_(mat1, mat2);
  return MatrixView(product_shape,
                    [&, cross_size](size_t rowIndex, size_t colIndex) {
                      auto out = mat1(rowIndex, 0) * mat2(0, colIndex);
                      for (size_t iz{1}; iz < cross_size; ++iz) {
                        out += mat1(rowIndex, iz) * mat2(iz, colIndex);
                      }
                      return out;
                    });
}

/// @brief Print a @p mat.
template<class T>
std::ostream& operator<<(std::ostream& out, BaseMatrixView<T> const& mat) {
  for (size_t rowIndex{0}; rowIndex < mat.NumRows(); ++rowIndex) {
    for (size_t colIndex{0}; colIndex < mat.NumCols(); ++colIndex) {
      out << mat(rowIndex, colIndex) << ' ';
    }
    out << std::endl;
  }
  return out;
}

} // namespace Storm
