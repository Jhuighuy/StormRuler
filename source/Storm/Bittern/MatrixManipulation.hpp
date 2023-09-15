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

#include <Storm/Crow/ArrayUtils.hpp>
#include <Storm/Crow/ConceptUtils.hpp>
#include <Storm/Crow/TupleUtils.hpp>
#include <Storm/Utils/Permutations.hpp>

#include <Storm/Bittern/MatrixMath.hpp>
#include <Storm/Bittern/MatrixView.hpp>
#include <Storm/Bittern/Shape.hpp>

#include <array>
#include <concepts>
#include <tuple>
#include <utility>

namespace Storm {

// -----------------------------------------------------------------------------

/// @brief Matrix roll view.
/// @todo This should be re-thinked.
template<array_like auto AxesArray, //
         matrix_view Matrix, tuple_like Shifts>
  requires boxable<Shifts>
class MatrixRollView final :
    public MatrixViewInterface<MatrixRollView<AxesArray, Matrix, Shifts>> {
private:

  static_assert(std::copyable<Shifts>, "Boxing is not implemented yet!");

  STORM_NO_UNIQUE_ADDRESS Matrix _mat;
  STORM_NO_UNIQUE_ADDRESS Shifts _shifts;

public:

  /// @brief Construct a matrix roll view.
  constexpr explicit MatrixRollView(Matrix mat, Shifts shifts)
      : _mat{std::move(mat)}, _shifts{std::move(shifts)} {}

  /// @brief Get the matrix shape.
  constexpr auto shape() const {
    return _mat.shape();
  }

  /// @brief Get the matrix element at @p indices.
  template<class... Indices>
    requires compatible_matrix_indices_v<MatrixRollView, Indices...>
  constexpr auto operator()(Indices... indices) const {
    STORM_ASSERT(in_range(shape(), indices...), "Indices are out of range!");
    const auto extents_tuple = _mat.shape();
    std::array indices_array{static_cast<size_t>(indices)...};
    ([&]<size_t... Axes>(std::index_sequence<Axes...>) {
      const auto roll_axes = [&](auto... axis_shift_pairs) {
        const auto roll_axis = [&](auto axis, ptrdiff_t shift) {
          auto& index = indices_array[axis];
          const auto extent = std::get<axis>(extents_tuple);
          index = (index + extent - shift) % extent;
        };
        (std::apply(roll_axis, std::move(axis_shift_pairs)), ...);
      };
      const std::tuple axes{fixed_size_t<Axes>{}...};
      zip_apply(roll_axes, {}, axes, _shifts);
    }(detail::_to_index_sequence<AxesArray>()));
    return std::apply(_mat, indices_array);
  }

}; // class MatrixRollView

/// @brief Create a matrix roll view.
template<array_like auto AxesArray, //
         viewable_matrix Matrix, tuple_like Shifts>
constexpr auto make_roll_view(Matrix&& mat, Shifts shifts) {
  using View = MatrixRollView<AxesArray, matrix_view_t<Matrix>, Shifts>;
  return View{std::forward<Matrix>(mat), std::move(shifts)};
}

/// @brief Roll the matrix @p mat elements along given axes.
template<size_t... Axes, viewable_matrix Matrix,
         std::convertible_to<ptrdiff_t>... Shifts>
  requires compatible_matrix_axes_v<Matrix, Axes...>
constexpr auto roll(Matrix&& mat, Shifts... shifts) {
  constexpr auto AxesArray = [] {
    constexpr size_t NumAxes = sizeof...(Axes);
    constexpr size_t NumShifts = sizeof...(Shifts);
#if STORM_COMPILER_MSVC
    // Bug in MSVC.
    static_cast<void>(NumAxes);
    static_cast<void>(NumShifts);
#endif
    if constexpr (NumAxes == 0) {
      return detail::_to_array(std::make_index_sequence<NumShifts>{});
    } else if constexpr (NumAxes == NumShifts) {
      // We have to use the `static_cast` due to the bug in MSVC.
      return std::array{static_cast<size_t>(Axes)...};
    } else {
      static_assert(detail::always_false<Matrix>,
                    "Incompatible amount of the specified axes and shifts!");
    }
  }();
  return make_roll_view<AxesArray>( //
      std::forward<Matrix>(mat), std::tuple{shifts...});
}

// -----------------------------------------------------------------------------

/// @brief Vector transpose view.
/// @todo This should be re-thinked.
template<matrix_view_r<1> Vector>
class VectorTransposeView final :
    public MatrixViewInterface<VectorTransposeView<Vector>> {
private:

  STORM_NO_UNIQUE_ADDRESS Vector _vec;

public:

  /// @brief Construct a vector transpose view.
  constexpr explicit VectorTransposeView(Vector vec) : _vec{std::move(vec)} {}

  /// @brief Get the matrix shape.
  constexpr auto shape() const {
    return cat_shapes(fixed_shape_t<1>{}, _vec.shape());
  }

  /// @brief Get the matrix element at @p row_index, @p col_index.
  template<class Index1, class Index2>
    requires compatible_matrix_indices_v<VectorTransposeView, Index1, Index2>
  constexpr auto operator()(Index1 index1, Index2 index2) const {
    STORM_ASSERT(in_range(shape(), index1, index2),
                 "Indices are out of range!");
    return _vec(index2);
  }

}; // class VectorTransposeView

template<class Vector>
VectorTransposeView(Vector&&) -> VectorTransposeView<matrix_view_t<Vector>>;

/// @brief Transpose the vector @p vec.
template<viewable_matrix_r<1> Vector>
constexpr auto transpose(Vector&& vec) {
  return VectorTransposeView{std::forward<Vector>(vec)};
}

// -----------------------------------------------------------------------------

/// @brief Matrix transpose view.
/// @tparam AxesArray Array-like permutation of the extents.
/// @todo This should be re-thinked.
template<array_like auto AxesArray, matrix_view Matrix>
class MatrixTransposeView final :
    public MatrixViewInterface<MatrixTransposeView<AxesArray, Matrix>> {
private:

  STORM_NO_UNIQUE_ADDRESS Matrix _mat;

  static constexpr auto _AxesInverseArray = [] {
    static_assert(permutations::is_permutation(AxesArray),
                  "Axes form an invalid permutation!");
    return permutations::invert_permutation(AxesArray);
  }();

public:

  /// @brief Construct a matrix transpose view.
  constexpr explicit MatrixTransposeView(Matrix mat) : _mat{std::move(mat)} {}

  /// @brief Get the matrix shape.
  constexpr auto shape() const {
    return ([&]<size_t... Axes>(std::index_sequence<Axes...>) {
      const auto extents = _mat.shape();
      return std::tuple{std::get<Axes>(extents)...};
    }(detail::_to_index_sequence<AxesArray>()));
  }

  /// @brief Get the matrix element at @p indices.
  template<class... Indices>
    requires compatible_matrix_indices_v<MatrixTransposeView, Indices...>
  constexpr auto operator()(Indices... indices) const {
    STORM_ASSERT(in_range(shape(), indices...), "Indices are out of range!");
    return ([&]<size_t... Axes>(std::index_sequence<Axes...>) {
      const std::tuple indices_tuple{std::move(indices)...};
      return _mat(std::get<Axes>(indices_tuple)...);
    }(detail::_to_index_sequence<_AxesInverseArray>()));
  }

}; // MatrixTransposeView

/// @brief Create a matrix transpose view.
template<array_like auto AxesArray, viewable_matrix Matrix>
constexpr auto make_matrix_transpose(Matrix&& mat) {
  using View = MatrixTransposeView<AxesArray, matrix_view_t<Matrix>>;
  return View{std::forward<Matrix>(mat)};
}

/// @brief Permute the matrix's @p mat extents.
/// @tparam Axes Permutation of the extents, if present.
template<size_t... Axes, viewable_matrix Matrix>
  requires compatible_matrix_axes_v<Matrix, Axes...> &&
           (sizeof...(Axes) != 0 || matrix_rank_v<Matrix> != 1)
constexpr auto transpose(Matrix&& mat) {
  constexpr auto AxesArray = ([] {
    constexpr size_t NumAxes = sizeof...(Axes);
    constexpr size_t Rank = matrix_rank_v<Matrix>;
#if STORM_COMPILER_MSVC
    // Bug in MSVC.
    static_cast<void>(NumAxes);
    static_cast<void>(Rank);
#endif
    if constexpr (NumAxes == 0) {
      return ([]<size_t... Extents>(std::index_sequence<Extents...>) {
        return std::array{(Rank - Extents - 1_sz)...};
      }(std::make_index_sequence<Rank>{}));
    } else if constexpr (NumAxes == Rank) {
      // We have to use the `static_cast` due to the bug in MSVC.
      return std::array{static_cast<size_t>(Axes)...}; // NOSONAR
    } else {
      static_assert(detail::always_false<Matrix>,
                    "Invalid numer of axes specified!");
    }
  }());
  return make_matrix_transpose<AxesArray>(std::forward<Matrix>(mat));
}

// -----------------------------------------------------------------------------

/// @brief Conjugate-transpose the matrix @p mat.
/// @tparam Axes Permutation of the extents.
template<size_t... Axes, viewable_matrix Matrix>
  requires compatible_matrix_axes_v<Matrix, Axes...>
constexpr auto conj_transpose(Matrix&& mat) {
  return conj(transpose<Axes...>(std::forward<Matrix>(mat)));
}

// -----------------------------------------------------------------------------

} // namespace Storm
