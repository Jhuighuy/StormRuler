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
#include <Storm/Crow/FunctionalUtils.hpp>
#include <Storm/Crow/TupleUtils.hpp>

#include <Storm/Bittern/Matrix.hpp>
#include <Storm/Bittern/MatrixView.hpp>
#include <Storm/Bittern/Shape.hpp>

#include <concepts>
#include <tuple>
#include <utility>

namespace Storm {

// -----------------------------------------------------------------------------

/// @brief Matrix builder view.
template<shape Shape, boxable Func>
  requires regular_applicable<Func, Shape> &&
           can_reference<apply_result_t<Func, Shape>>
class MatrixBuilderView final :
    public MatrixViewInterface<MatrixBuilderView<Shape, Func>> {
private:

  // static_assert(std::copyable<Func>, "Boxing is not implemented yet!");
  static_assert(!matrix<apply_result_t<Func, Shape>>,
                "Block matrix builders are not supported yet!");

  STORM_NO_UNIQUE_ADDRESS Shape _shape;
  STORM_NO_UNIQUE_ADDRESS Func _func;

public:

  /// @brief Construct a builder view.
  constexpr MatrixBuilderView(Shape shape, Func func)
      : _shape{shape}, _func{std::move(func)} {}

  /// @brief Get the matrix shape.
  constexpr auto shape() const {
    return _shape;
  }

  /// @brief Get the matrix element at @p indices.
  template<class... Indices>
    requires compatible_matrix_indices_v<MatrixBuilderView, Indices...>
  constexpr auto operator()(Indices... indices) const {
    STORM_ASSERT(in_range(shape(), indices...), "Indices are out of range!");
    return _func(indices...);
  }

}; // class MatrixBuilderView

namespace detail {
  template<class Shape, class Func>
  concept can_matrix_builder_view = requires {
    MatrixBuilderView{std::declval<Shape>(), std::declval<Func>()};
  };
} // namespace detail

/// @brief Build a matrix of shape @p shape with function @p func.
template<shape Shape, class Func>
  requires detail::can_matrix_builder_view<Shape, Func>
constexpr auto build(Shape shape, Func func) {
  return MatrixBuilderView{shape, std::move(func)};
}

// -----------------------------------------------------------------------------

/// @brief Build a constant matrix with shape @p shape
/// filled with elements equal to @p elem.
template<scalar Elem, shape Shape>
constexpr auto values(Shape shape, Elem elem) noexcept {
  return MatrixBuilderView{shape, Constant{std::move(elem)}};
}

/// @brief Build a matrix with shape @p shape filled with zeroes.
template<scalar Elem = real_t, shape Shape>
  requires std::constructible_from<Elem, int>
constexpr auto zeroes(Shape shape) noexcept {
  return values(shape, Elem{0});
}
/// @brief Build a matrix with extents @p extents filled with zeroes.
template<scalar Elem = real_t, extent... Extents>
  requires std::constructible_from<Elem, int>
constexpr auto zeroes(Extents... extents) noexcept {
  return zeroes<Elem>(std::tuple{extents...});
}

/// @brief Build a matrix of shape @p shape filled with ones.
template<scalar Elem = real_t, shape Shape>
  requires std::constructible_from<Elem, int>
constexpr auto ones(Shape shape) noexcept {
  return values(shape, Elem{1});
}
/// @brief Build a matrix with extents @p extents filled with ones.
template<scalar Elem = real_t, extent... Extents>
  requires std::constructible_from<Elem, int>
constexpr auto ones(Extents... extents) noexcept {
  return ones<Elem>(std::tuple{extents...});
}

// -----------------------------------------------------------------------------

/// @brief Build a diagonal matrix with @p shape with a diagonal element
/// value @p diagonal and off-diagonal element value @p off_diagonal.
template<scalar Elem = real_t, shape Shape>
constexpr auto eye(Shape shape, //
                   Elem diagonal = Elem{1},
                   Elem off_diagonal = Elem{0}) noexcept {
  Eye eye{std::move(diagonal), std::move(off_diagonal)};
  return MatrixBuilderView{shape, std::move(eye)};
}
/// @brief Build an identity matrix with extents @p extents.
template<scalar Elem = real_t, extent... Extents>
constexpr auto eye(Extents... extents) noexcept {
  return eye<Elem>(std::tuple{extents...});
}

// -----------------------------------------------------------------------------

/// @brief Build a matrix with @p extent evenly spaced elements
/// over the specified interval @p start, @p stop.
template<scalar Elem = real_t, extent Extent>
constexpr auto linspace(Elem start, Elem stop, Extent extent) noexcept {
  const auto step = (stop - start) / (extent + 1);
  const auto func = [start, step](size_t i) { return start + i * step; };
  return MatrixBuilderView{std::tuple{extent}, std::move(func)};
}

/// @brief Build a matrix with @p extent evenly spaced elements
/// (on a log scale) over the specified interval @p start, @p stop.
template<scalar Elem = real_t, extent Extent>
constexpr auto logspace(Elem start, Elem stop, Extent extent) noexcept {
  static_assert(detail::always_false<Elem>,
                "`linspace` is not implemented yet!");
}

// -----------------------------------------------------------------------------

} // namespace Storm
