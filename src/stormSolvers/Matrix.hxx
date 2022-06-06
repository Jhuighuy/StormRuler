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

#include <stormBase.hxx>
#include <stormSolvers/MatrixView.hxx>

namespace Storm {

template<class Derived>
class BaseMatrix : public BaseMatrixView<Derived> {
public:

  /// @copydoc BaseMatrixView::operator()
  constexpr auto& operator()(size_t rowIndex, size_t colIndex) noexcept {
    return static_cast<Derived&>(*this)(rowIndex, colIndex);
  }

  template<class T>
  constexpr auto& assign(BaseMatrixView<T> const& mat) noexcept {
    STORM_ASSERT_(mat.shape() == this->shape());
    for (size_t rowIndex{0}; rowIndex < this->NumRows(); ++rowIndex) {
      for (size_t colIndex{0}; colIndex < this->NumCols(); ++colIndex) {
        (*this)(rowIndex, colIndex) = mat(rowIndex, colIndex);
      }
    }
    return static_cast<Derived&>(*this);
  }

}; // class BaseMatrix

template<class T1, class T2>
constexpr real_t dot_product(BaseMatrixView<T1> const& mat1,
                             BaseMatrixView<T2> const& mat2) {
  real_t d{};
#pragma omp parallel for reduction(+ : d)
  for (int rowIndex = 0; rowIndex < (int) mat1.NumRows(); ++rowIndex) {
    for (size_t colIndex{0}; colIndex < mat1.NumCols(); ++colIndex) {
      d += mat1(rowIndex, colIndex) * mat2(rowIndex, colIndex);
    }
  }
  return d;
}

template<class T1>
constexpr real_t norm_2(BaseMatrixView<T1> const& mat1) {
  return std::sqrt(dot_product(mat1, mat1));
}

template<class T1, class T2>
constexpr auto& operator<<=(BaseMatrix<T1>& mat1,
                            BaseMatrixView<T2> const& mat2) {
#pragma omp parallel for
  for (int rowIndex = 0; rowIndex < (int) mat1.NumRows(); ++rowIndex) {
    for (size_t colIndex{0}; colIndex < mat1.NumCols(); ++colIndex) {
      mat1(rowIndex, colIndex) = mat2(rowIndex, colIndex);
    }
  }
  return mat1;
}

template<class T1, class T2>
constexpr auto& operator+=(BaseMatrix<T1>& mat1,
                           BaseMatrixView<T2> const& mat2) {
  return mat1 <<= mat1 + mat2;
}
template<class T1, class T2>
constexpr auto& operator-=(BaseMatrix<T1>& mat1,
                           BaseMatrixView<T2> const& mat2) {
  return mat1 <<= mat1 - mat2;
}

template<class T1, class V2>
constexpr auto& operator*=(BaseMatrix<T1>& mat1, V2 const& val2) {
  return mat1 <<= val2 * mat1;
}
template<class T1, class V2>
constexpr auto& operator/=(BaseMatrix<T1>& mat1, V2 const& val2) {
  return mat1 <<= mat1 / val2;
}

} // namespace Storm
