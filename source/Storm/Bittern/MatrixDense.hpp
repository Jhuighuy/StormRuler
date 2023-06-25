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
#include <Storm/Crow/TupleUtils.hpp>

#include <Storm/Bittern/Matrix.hpp>
#include <Storm/Bittern/MatrixTarget.hpp>
#include <Storm/Bittern/Shape.hpp>

#include <algorithm>
#include <array>
#include <concepts>
#include <initializer_list>
#include <ranges>
#include <type_traits>
#include <utility>

namespace Storm {

// -----------------------------------------------------------------------------

/// @brief Dense matrix.
template<scalar Elem, shape Shape>
class DenseMatrix;

/// @brief Dense matrix with fixed shape.
template<scalar Elem, size_t... Extents>
using FixedMatrix = DenseMatrix<Elem, fixed_shape_t<Extents...>>;

/// @brief Dense vector with fixed shape.
template<scalar Elem, size_t Extent>
using FixedVector = FixedMatrix<Elem, Extent>;

#if 0
/// @brief Dense matrix with dynamic shape.
template<scalar Elem, size_t Rank = 2>
using DynamicMatrix = DenseMatrix<Elem, dynamic_shape_t<Rank>>;

/// @brief Dense vector with dynamic shape.
template<scalar Elem>
using DynamicVector = DynamicMatrix<Elem, 1>;
#endif

// -----------------------------------------------------------------------------

/// @brief Construct a dense matrix from the arguments @p args.
template<class... Args>
constexpr auto to_matrix(Args&&... args) noexcept {
  return DenseMatrix{std::forward<Args>(args)...};
}

/// @brief Dense matrix type for the arguments.
template<class... Args>
using matrix_t = decltype(to_matrix(std::declval<Args>()...));

template<matrix Matrix>
DenseMatrix(Matrix&&)
    -> DenseMatrix<matrix_element_t<Matrix>, matrix_shape_t<Matrix>>;

// Construct a vector with elements.
// E.g. `DenseMatrix{1, 2, 3}`.
template<scalar... Elems>
DenseMatrix(Elems...) -> DenseMatrix<std::common_type_t<Elems...>,
                                     fixed_shape_t<sizeof...(Elems)>>;
// E.g. `DenseMatrix{std::array{1, 2, 3}}`.
template<std::ranges::range ElemsRange>
  requires scalar<std::ranges::range_value_t<ElemsRange>>
DenseMatrix(ElemsRange&&)
    -> DenseMatrix<std::ranges::range_value_t<ElemsRange>, dynamic_shape_t<1>>;

namespace detail {
  // Can a non-1x1 matrix be constructed from this argument?
  template<class Slice>
  concept _can_slice = requires { typename matrix_t<Slice>; } &&
                       (!std::same_as<matrix_t<Slice>, FixedMatrix<Slice, 1>>);
} // namespace detail

// Construct a matrix with slices.
// E.g. `DenseMatrix{std::array{1, 2, 3}, std::array{4, 5, 6}}`.
// E.g. `DenseMatrix{DenseMatrix{1, 2, 3}, DenseMatrix{4, 5, 6}}`.
template<class... Slices>
  requires ((sizeof...(Slices) > 1) && (detail::_can_slice<Slices> && ...))
DenseMatrix(Slices&&...) -> DenseMatrix<
    std::common_type_t<matrix_element_t<matrix_t<Slices>>...>,
    cat_shapes_t<fixed_shape_t<sizeof...(Slices)>,
                 common_shape_t<matrix_shape_t<matrix_t<Slices>>...>>>;
// E.g. `DenseMatrix{std::array{std::array{1, 2, 3}, std::array{4, 5, 6}}}`.
// E.g. `DenseMatrix{std::array{DenseMatrix{1, 2, 3}, DenseMatrix{4, 5, 6}}}`.
template<std::ranges::range SlicesRange>
  requires detail::_can_slice<std::ranges::range_value_t<SlicesRange>>
DenseMatrix(SlicesRange&&) -> DenseMatrix<
    matrix_element_t<matrix_t<std::ranges::range_value_t<SlicesRange>>>,
    cat_shapes_t<
        dynamic_shape_t<1>,
        matrix_shape_t<matrix_t<std::ranges::range_value_t<SlicesRange>>>>>;

namespace detail {
  template<class Type, size_t Size>
  using array_cref_t_ = const Type (&)[Size]; // NOSONAR
} // namespace detail

template<scalar... Elems, size_t Size>
DenseMatrix(detail::array_cref_t_<Elems, Size>...)
    -> DenseMatrix<std::common_type_t<Elems...>,
                   fixed_shape_t<sizeof...(Elems), Size>>;
template<class... Slices, size_t Size>
  requires ((sizeof...(Slices) > 1) && (detail::_can_slice<Slices> && ...))
DenseMatrix(detail::array_cref_t_<Slices, Size>...) -> DenseMatrix<
    std::common_type_t<matrix_element_t<matrix_t<Slices>>...>,
    cat_shapes_t<fixed_shape_t<sizeof...(Slices), Size>,
                 common_shape_t<matrix_shape_t<matrix_t<Slices>>...>>>;

// -----------------------------------------------------------------------------

/// @brief Dense matrix (fixed shape specialization).
template<scalar Elem, fixed_shape Shape>
class DenseMatrix<Elem, Shape> final :
    public TargetMatrixInterface<DenseMatrix<Elem, Shape>> {
private:

  static constexpr bool _Multirank = shape_rank_v<Shape> > 1;

  static constexpr auto _FirstExtent =
      static_cast<size_t>(shape_extent_t<Shape, 0>{});

  using _SliceShape = decltype( //
      drop_n_apply<1>(to_tuple, std::declval<Shape>()));

  using _Slice =
      std::conditional_t<_Multirank, DenseMatrix<Elem, _SliceShape>, Elem>;

  STORM_NO_UNIQUE_ADDRESS std::array<_Slice, _FirstExtent> _slices{};

public:

  // ---------------------------------------------------------------------------

  /// @brief Construct a matrix.
  DenseMatrix() = default;

  // ---------------------------------------------------------------------------

  /// @brief Construct a matrix with elements of the matrix @p mat.
  template<matrix Matrix>
    requires assignable_matrix<DenseMatrix, Matrix> &&
             different_from<DenseMatrix, Matrix>
  constexpr DenseMatrix(Matrix&& mat) { // NOSONAR
    this->assign(std::forward<Matrix>(mat));
  }

  /// @brief Assign the current matrix elements from matrix @p mat.
  template<matrix Matrix>
    requires assignable_matrix<DenseMatrix, Matrix> &&
             different_from<DenseMatrix, Matrix>
  constexpr DenseMatrix& operator=(Matrix&& mat) {
    return this->assign(std::forward<Matrix>(mat));
  }

  // ---------------------------------------------------------------------------

  /// @brief Construct a matrix with the slice range @p slices.
  /// @{
  constexpr DenseMatrix(std::initializer_list<_Slice> slices) {
    std::ranges::copy(slices, _slices.begin());
  }
  template<std::ranges::input_range SlicesRange>
  constexpr explicit DenseMatrix(SlicesRange&& slices) {
    std::ranges::copy(slices, _slices.begin());
  }
  /// @}

  /// @brief Assign the current matrix to from slices range @p slices.
  /// @{
  constexpr DenseMatrix& operator=(std::initializer_list<_Slice> slices) {
    std::ranges::copy(slices, _slices.begin());
    return *this;
  }
  template<std::ranges::input_range SlicesRange>
  constexpr DenseMatrix& operator=(SlicesRange&& slices) {
    std::ranges::copy(slices, _slices.begin());
    return *this;
  }
  /// @}

  // ---------------------------------------------------------------------------

  /// @brief Get the matrix shape.
  constexpr auto shape() const noexcept {
    return Shape{};
  }

  /// @brief Get the slice at @p index.
  /// @{
  constexpr _Slice& operator[](size_t slice_index) noexcept {
    STORM_ASSERT(slice_index < _FirstExtent, "Slice index is out of range!");
    return _slices[slice_index];
  }
  constexpr const _Slice& operator[](size_t slice_index) const noexcept {
    STORM_ASSERT(slice_index < _FirstExtent, "Slice index is out of range!");
    return _slices[slice_index];
  }
  /// @}

  /// @brief Get the matrix element at @p indices.
  /// @{
  template<class... Indices>
    requires compatible_matrix_indices_v<DenseMatrix, Indices...>
  constexpr Elem& operator()(Indices... indices) noexcept {
    STORM_ASSERT(in_range(shape(), indices...), "Indices are out of range!");
    return [&](auto first_index, auto... rest_indices) -> Elem& {
      _Slice& slice = _slices[first_index];
      if constexpr (_Multirank) return slice(rest_indices...);
      else return slice;
    }(indices...);
  }
  template<class... Indices>
    requires compatible_matrix_indices_v<DenseMatrix, Indices...>
  constexpr const Elem& operator()(Indices... indices) const noexcept {
    STORM_ASSERT(in_range(shape(), indices...), "Indices are out of range!");
    return [&](auto first_index, auto... rest_indices) -> const Elem& {
      const _Slice& slice = _slices[first_index];
      if constexpr (_Multirank) return slice(rest_indices...);
      else return slice;
    }(indices...);
  }
  /// @}

}; // class DenseMatrix

// -----------------------------------------------------------------------------

/// @brief Dense matrix.
template<scalar Elem, shape Shape>
class DenseMatrix final {
  static_assert(detail::always_false<Shape>,
                "Dynamic matrices are not implemented yet!");
}; // class DenseMatrix

// -----------------------------------------------------------------------------

} // namespace Storm
