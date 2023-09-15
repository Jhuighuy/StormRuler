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
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <Storm/Base.hpp>

#include <Storm/Crow/ConceptUtils.hpp>
#include <Storm/Utils/Index.hpp>

#include <Storm/Bittern/Shape.hpp>

#include <concepts>
#include <tuple>
#include <type_traits>
#include <utility>

namespace Storm {

// -----------------------------------------------------------------------------

/// @brief Matrix: has shape and subscripts.
template<class Matrix>
concept matrix = //
    requires(Matrix& mat) {
      { mat.shape() } -> shape;
      { std::apply(mat, mat.shape()) } -> can_reference;
    };

/// @brief Scalar: copyable object that is not a matrix.
template<class Scalar>
concept scalar =
    std::is_object_v<Scalar> && std::copyable<Scalar> && (!matrix<Scalar>);

// -----------------------------------------------------------------------------

/// @brief Matrix shape type.
template<class Matrix>
using matrix_shape_t = decltype(std::declval<Matrix>().shape());

/// @brief Matrix rank value.
template<class Matrix>
inline constexpr size_t matrix_rank_v = shape_rank_v<matrix_shape_t<Matrix>>;

/// @brief Matrix of a specified rank.
template<class Matrix, size_t Rank>
concept matrix_r = matrix<Matrix> && (matrix_rank_v<Matrix> == Rank);

/// @brief Check for matrices to be of the same ranks.
template<matrix Matrix, matrix... RestMatrices>
inline constexpr bool compatible_matrices_v =
    ((matrix_rank_v<Matrix> == matrix_rank_v<RestMatrices>) &&...);

template<class Matrix, index... Indices>
inline constexpr bool compatible_matrix_indices_v =
    (matrix_rank_v<Matrix> == sizeof...(Indices)) &&
    (detail::index_like<Indices> && ...);

template<matrix Matrix, size_t... Axes>
inline constexpr bool compatible_matrix_axes_v =
    ((matrix_rank_v<Matrix> > Axes) && ...);

/// @brief Number of the matrix rows.
template<matrix Matrix>
  requires (matrix_rank_v<Matrix> >= 1)
constexpr size_t num_rows(Matrix&& mat) {
  return std::get<0>(mat.shape());
}

/// @brief Number of the matrix columns.
template<matrix Matrix>
  requires (matrix_rank_v<Matrix> >= 2)
constexpr size_t num_cols(Matrix&& mat) {
  return std::get<1>(mat.shape());
}

// -----------------------------------------------------------------------------

/// @brief Matrix element reference type.
template<matrix Matrix>
using matrix_element_ref_t = decltype( //
    std::apply(std::declval<Matrix&>(), std::declval<Matrix>().shape()));

/// @brief Matrix element type.
template<matrix Matrix>
using matrix_element_t = std::remove_cvref_t<matrix_element_ref_t<Matrix>>;

/// @brief Matrix with assignable elements.
template<class Matrix>
concept output_matrix =
    matrix<Matrix> &&
    std::same_as<matrix_element_ref_t<Matrix>,
                 std::add_lvalue_reference_t<matrix_element_t<Matrix>>>;

/// @brief Matrix elements are assignable from the elements of another matrix.
template<class OutMatrix, class Matrix>
concept assignable_matrix =
    output_matrix<OutMatrix> && matrix<Matrix> &&
    compatible_matrices_v<OutMatrix, Matrix> &&
    std::assignable_from<matrix_element_ref_t<OutMatrix>,
                         matrix_element_t<Matrix>>;

// -----------------------------------------------------------------------------

} // namespace Storm

#include <Storm/Bittern/MatrixAlgorithms.hpp>
#include <Storm/Bittern/MatrixTarget.hpp>

#include <Storm/Bittern/MatrixBuilder.hpp>
#include <Storm/Bittern/MatrixDense.hpp>
#include <Storm/Bittern/MatrixManipulation.hpp>
#include <Storm/Bittern/MatrixMath.hpp>
#include <Storm/Bittern/MatrixPart.hpp>
#include <Storm/Bittern/MatrixProduct.hpp>
