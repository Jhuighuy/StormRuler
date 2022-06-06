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
  constexpr auto& operator()(size_t row_index, size_t col_index) noexcept {
    return static_cast<Derived&>(*this)(row_index, col_index);
  }

  template<class T>
  constexpr auto& assign(const BaseMatrixView<T>& mat) noexcept {
    STORM_ASSERT_(mat.shape() == this->shape());
    for (size_t row_index{0}; row_index < this->num_rows(); ++row_index) {
      for (size_t col_index{0}; col_index < this->num_cols(); ++col_index) {
        (*this)(row_index, col_index) = mat(row_index, col_index);
      }
    }
    return static_cast<Derived&>(*this);
  }

}; // class BaseMatrix

template<class T1, class T2>
constexpr real_t dot_product(const BaseMatrixView<T1>& mat1,
                             const BaseMatrixView<T2>& mat2) {
  real_t d{};
  //#pragma omp parallel for reduction(+ : d)
  for (int row_index = 0; row_index < (int) mat1.num_rows(); ++row_index) {
    for (size_t col_index{0}; col_index < mat1.num_cols(); ++col_index) {
      d += mat1(row_index, col_index) * mat2(row_index, col_index);
    }
  }
  return d;
}

template<class T1>
constexpr real_t norm_2(const BaseMatrixView<T1>& mat1) {
  return std::sqrt(dot_product(mat1, mat1));
}

template<class T1, class T2>
constexpr auto& operator<<=(BaseMatrix<T1>& mat1,
                            const BaseMatrixView<T2>& mat2) {
#pragma omp parallel for
  for (int row_index = 0; row_index < (int) mat1.num_rows(); ++row_index) {
    for (size_t col_index{0}; col_index < mat1.num_cols(); ++col_index) {
      mat1(row_index, col_index) = mat2(row_index, col_index);
    }
  }
  return mat1;
}

template<class T1, class T2>
constexpr auto& operator+=(BaseMatrix<T1>& mat1,
                           const BaseMatrixView<T2>& mat2) {
  return mat1 <<= mat1 + mat2;
}
template<class T1, class T2>
constexpr auto& operator-=(BaseMatrix<T1>& mat1,
                           const BaseMatrixView<T2>& mat2) {
  return mat1 <<= mat1 - mat2;
}

template<class T1, class V2>
constexpr auto& fill_with(BaseMatrix<T1>& mat1, const V2& val2) {
#pragma omp parallel for
  for (int row_index = 0; row_index < (int) mat1.num_rows(); ++row_index) {
    for (size_t col_index{0}; col_index < mat1.num_cols(); ++col_index) {
      mat1(row_index, col_index) = val2;
    }
  }
  return mat1;
}

template<class T1, class V2>
constexpr auto& fill_diag_with(BaseMatrix<T1>& mat1, const V2& val2) {
#pragma omp parallel for
  for (int row_index = 0; row_index < (int) mat1.num_rows(); ++row_index) {
    for (size_t col_index{0}; col_index < mat1.num_cols(); ++col_index) {
      mat1(row_index, col_index) = (size_t) row_index == col_index ? val2 : 0;
    }
  }
  return mat1;
}

template<class T1, class... V2>
constexpr auto& fill_randomly(BaseMatrix<T1>& mat1, V2...) {
  STORM_ENSURE_(!"Not implemented");
  return mat1;
}

template<class T1, class V2>
constexpr auto& operator*=(BaseMatrix<T1>& mat1, const V2& val2) {
  return mat1 <<= val2 * mat1;
}
template<class T1, class V2>
constexpr auto& operator/=(BaseMatrix<T1>& mat1, const V2& val2) {
  return mat1 <<= mat1 / val2;
}

} // namespace Storm
