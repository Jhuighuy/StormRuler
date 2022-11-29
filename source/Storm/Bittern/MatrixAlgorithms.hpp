// Copyright (C) 2020-2023 Oleg Butakov
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR Allocator PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT
// SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
// OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

#pragma once

#include <Storm/Base.hpp>

#include <Storm/Bittern/Matrix.hpp>

#include <algorithm>
#include <ostream>
#include <random>
#include <utility>

namespace Storm
{

/// @brief Assign the matrices.
/// @{
template<output_matrix OutMatrix, matrix Matrix>
constexpr OutMatrix& assign(OutMatrix&& out_mat, Matrix&& mat) noexcept
{
  STORM_ASSERT_(out_mat.shape() == mat.shape(), "Matrix shapes do not match!");
  for (size_t row_index = 0; row_index < num_rows(out_mat); ++row_index) {
    for (size_t col_index = 0; col_index < num_cols(out_mat); ++col_index) {
      out_mat(row_index, col_index) = mat(row_index, col_index);
    }
  }
  return out_mat;
}
template<output_matrix OutMatrix, class AssignFunc, matrix... Matrices>
constexpr OutMatrix& assign(AssignFunc assign_func, //
                            OutMatrix&& out_mat, Matrices&&... mats) noexcept
{
  STORM_ASSERT_((out_mat.shape() == mats.shape()) && ...,
                "Matrix shapes do not match!");
  for (size_t row_index = 0; row_index < num_rows(out_mat); ++row_index) {
    for (size_t col_index = 0; col_index < num_cols(out_mat); ++col_index) {
      assign_func(out_mat(row_index, col_index), mats(row_index, col_index)...);
    }
  }
  return out_mat;
}
/// @}

/// @brief Assign the matrices.
/// @todo To be removed!
template<output_matrix OutMatrix, matrix Matrix>
constexpr OutMatrix& operator<<=(OutMatrix&& out_mat, Matrix&& mat) noexcept
{
  return assign(std::forward<OutMatrix>(out_mat), std::forward<Matrix>(mat));
}

// -----------------------------------------------------------------------------

/// @brief Fill the matrix @p out_mat with a scalar @p scal.
template<matrix OutMatrix, std::copyable Scalar>
  requires (!matrix<Scalar>)
constexpr OutMatrix& fill(OutMatrix& out_mat, Scalar scal)
{
  return assign(out_mat, [scal = std::move(scal)](auto& out_elem) noexcept {
    out_elem = scal;
  });
}

/// @brief Fill the matrix @p out_mat elements with the random numbers.
/// @warning This is a sequential operation!
template<output_matrix OutMatrix>
  requires (real_matrix<OutMatrix>)
constexpr OutMatrix& fill_randomly(OutMatrix&& out_mat, //
                                   matrix_element_t<OutMatrix> min = 0,
                                   matrix_element_t<OutMatrix> max = 1) noexcept
{
  STORM_CPP23_STATIC_ std::mt19937_64 random_engine{};
  std::uniform_real_distribution distribution{min, max};
  for (size_t row_index = 0; row_index < num_rows(out_mat); ++row_index) {
    for (size_t col_index = 0; col_index < num_cols(out_mat); ++col_index) {
      out_mat(row_index, col_index) = distribution(random_engine);
    }
  }
  return out_mat;
}

// -----------------------------------------------------------------------------

/// @brief Print a @p mat.
/// @todo This is too trivial implementation, we need something fancier :)
std::ostream& operator<<(std::ostream& out, const matrix auto& mat)
{
  for (size_t row_index = 0; row_index < num_rows(mat); ++row_index) {
    out << "( ";
    for (size_t col_index = 0; col_index < num_cols(mat); ++col_index) {
      out << mat(row_index, col_index) << " ";
    }
    out << ")" << std::endl;
  }
  return out;
}

} // namespace Storm
