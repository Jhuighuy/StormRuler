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

#include <Storm/Base.hpp>

#include <Storm/Bittern/MatrixView.hpp>

namespace Storm {

/// @name Matrix actions.
/// @{

/// @brief Assign the matrix @p mat2 to @p mat1.
constexpr auto& operator<<=(matrix auto&& mat1, matrix auto&& mat2) noexcept {
  return eval([](auto& val1, const auto& val2) { val1 = val2; }, mat1, mat2);
}

/// @name Functional actions.
/// @{

/// @brief Multiply-assign the matrix @p mat by a scalar @p scal.
/// @todo `scal` should really be auto!
constexpr auto& operator*=(matrix auto&& mat, real_t scal) noexcept {
  return eval([scal](auto& val) { val *= scal; }, mat);
}

/// @brief Divide-assign the matrix @p mat by a scalar @p scal.
/// @todo `scal` should really be auto!
constexpr auto& operator/=(matrix auto&& mat, real_t scal) noexcept {
  return eval([scal](auto& val) { val /= scal; }, mat);
}

/// @brief Add-assign the matrices @p mat1 and @p mat2.
constexpr auto& operator+=(matrix auto&& mat1, matrix auto&& mat2) noexcept {
  return eval([](auto& val1, const auto& val2) { val1 += val2; }, mat1, mat2);
}

/// @brief Subtract-assign the matrices @p mat1 and @p mat2.
constexpr auto& operator-=(matrix auto&& mat1, matrix auto&& mat2) noexcept {
  return eval([](auto& val1, const auto& val2) { val1 -= val2; }, mat1, mat2);
}

/// @brief Component-wise multiply-assign the matrices @p mat1 and @p mat2.
constexpr auto& operator*=(matrix auto&& mat1, matrix auto&& mat2) noexcept {
  return eval([](auto& val1, const auto& val2) { val1 *= val2; }, mat1, mat2);
}

/// @brief Component-wise divide-assign the matrices @p mat1 and @p mat2.
constexpr auto& operator/=(matrix auto&& mat1, matrix auto&& mat2) noexcept {
  return eval([](auto& val1, const auto& val2) { val1 /= val2; }, mat1, mat2);
}

/// @} // Functional actions.

/// @name Filling actions.
/// @{

constexpr auto& fill_diag_with(matrix auto&& mat, auto scal) noexcept {
  return mat <<= make_diagonal_matrix(mat.shape(), scal);
}

/// @brief Fill the matrix @p mat with a scalar @p scal.
constexpr auto& fill_with(matrix auto&& mat, auto scal) noexcept {
  return eval([scal](auto& val) { val = scal; }, mat);
}

/// @brief Fill the matrix @p mat with the random numbers.
constexpr auto& fill_randomly(matrix auto&& mat) noexcept {
  for (size_t row_index{0}; row_index < num_rows(mat); ++row_index) {
    for (size_t col_index{0}; col_index < num_cols(mat); ++col_index) {
      mat(row_index, col_index) = 2.0 * (real_t(rand()) / RAND_MAX) - 1.0;
    }
  }
  return mat;
}

/// @} // Filling actions.

/// @} // Matrix actions.

} // namespace Storm