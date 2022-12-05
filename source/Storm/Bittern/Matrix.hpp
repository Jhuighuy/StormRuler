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

#include <Storm/Utils/Index.hpp>

#include <Storm/Bittern/Math.hpp>

#include <concepts>
#include <tuple>
#include <type_traits>
#include <utility>

namespace Storm
{

// -----------------------------------------------------------------------------

/// @brief Shape.
template<size_t Rank>
using shape_t = std::array<size_t, Rank>;

/// @brief Matrix shape: instantiation of shape_t.
template<class MatrixShape>
concept matrix_shape =
    requires { std::tuple_size_v<MatrixShape>; } &&
    (std::tuple_size_v<MatrixShape> != 0) &&
    std::same_as<MatrixShape, shape_t<std::tuple_size_v<MatrixShape>>>;

/// @brief Matrix: has shape and subscripts.
template<class Matrix>
concept matrix = //
    requires(Matrix& mat) {
      // clang-format off
      { mat.shape() } -> matrix_shape;
      { std::apply(mat, mat.shape()) } /*-> referenceable */;
      // clang-format on
    };

// -----------------------------------------------------------------------------

/// @brief Matrix shape type.
template<class Matrix>
using matrix_shape_t = decltype(std::declval<Matrix>().shape());

/// @brief Matrix rank value.
template<class Matrix>
inline constexpr size_t matrix_rank_v =
    std::tuple_size_v<matrix_shape_t<Matrix>>;

/// @brief Matrix of a specified rank.
template<class Matrix, size_t Rank>
concept matrix_r = matrix<Matrix> && (matrix_rank_v<Matrix> == Rank);

/// @brief Number of the matrix rows.
template<matrix Matrix>
constexpr size_t num_rows(Matrix&& mat) noexcept
{
  return mat.shape()[0];
}

/// @brief Number of the matrix columns.
template<matrix Matrix>
  requires (matrix_rank_v<Matrix> >= 2)
constexpr size_t num_cols(Matrix&& mat) noexcept
{
  return mat.shape()[1];
}

// -----------------------------------------------------------------------------

template<matrix Matrix, matrix... RestMatrices>
inline constexpr bool compatible_matrices_v =
    ((matrix_rank_v<Matrix> == matrix_rank_v<RestMatrices>) &&...);

template<class Matrix, class... Indices>
inline constexpr bool compatible_matrix_indices_v =
    matrix_rank_v<Matrix> == sizeof...(Indices);

/// @brief Check if all indices are in range.
template<size_t Rank, index... Indices>
  requires (Rank == sizeof...(Indices))
constexpr bool in_range(const shape_t<Rank>&, Indices...)
{
  return true; // to be implemented.
}

// -----------------------------------------------------------------------------

/// @brief Matrix element type, as is.
template<matrix Matrix>
using matrix_element_decltype_t =
    decltype(std::declval<Matrix>()(size_t{}, size_t{}));

/// @brief Matrix element type.
template<matrix Matrix>
using matrix_element_t = std::remove_cvref_t<matrix_element_decltype_t<Matrix>>;

/// @brief Matrix with assignable elements.
template<class Matrix>
concept output_matrix =
    matrix<Matrix> && /// @todo Check if reference is not const!
    std::is_lvalue_reference_v<matrix_element_decltype_t<Matrix>>;

/// @brief Matrix element reference type.
template<output_matrix Matrix>
using matrix_element_ref_t = matrix_element_decltype_t<Matrix>;

/// @brief Matrix with boolean elements.
template<class Matrix>
concept bool_matrix = matrix<Matrix> && bool_type<matrix_element_t<Matrix>>;

/// @brief Matrix with integral elements.
/// We do not have any special integer matrix operations,
/// integer matrices are defined for completeness.
template<class Matrix>
concept integer_matrix =
    matrix<Matrix> && integer_type<matrix_element_t<Matrix>>;

/// @brief Matrix with real (floating-point) elements.
template<class Matrix>
concept real_matrix = matrix<Matrix> && real_type<matrix_element_t<Matrix>>;

/// @brief Matrix with complex (floating-point) elements.
template<class Matrix>
concept complex_matrix =
    matrix<Matrix> && complex_type<matrix_element_t<Matrix>>;

/// @brief Matrix with real or complex (floating-point) elements.
template<class Matrix>
concept real_or_complex_matrix = real_matrix<Matrix> || complex_matrix<Matrix>;

/// @brief Matrix with numerical elements.
template<class Matrix>
concept numeric_matrix =
    integer_matrix<Matrix> || real_or_complex_matrix<Matrix>;

// -----------------------------------------------------------------------------

} // namespace Storm

#include <Storm/Bittern/MatrixAlgorithms.hpp>
#include <Storm/Bittern/MatrixTarget.hpp>

// #include <Storm/Bittern/MatrixGenerator.hpp>
#include <Storm/Bittern/MatrixMath.hpp>
#include <Storm/Bittern/MatrixProduct.hpp>
#include <Storm/Bittern/MatrixTranspose.hpp>
// #include <Storm/Bittern/MatrixSlicing.hpp>
