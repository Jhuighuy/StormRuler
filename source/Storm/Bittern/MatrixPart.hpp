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
#include <Storm/Bittern/MatrixView.hpp>
#include <Storm/Bittern/Shape.hpp>

#include <algorithm>
#include <array>
#include <concepts>
#include <ranges>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace Storm {

// -----------------------------------------------------------------------------

/// @brief New axis placeholder type.
class NewAxis final {
public:

  /// @brief Construct the new axis placeholder.
  explicit NewAxis() = default;

}; // class NewAxis

/// @brief New axis placeholder.
inline constexpr NewAxis new_axis{};

// -----------------------------------------------------------------------------

/// @brief Mask matrix view.
template<matrix_view Matrix>
  requires boolean_testable<matrix_element_t<Matrix>>
class MaskMatrixView final :
    public MatrixViewInterface<MaskMatrixView<Matrix>> {
private:

  using _Indices = std::array<size_t, matrix_rank_v<Matrix>>;

  STORM_NO_UNIQUE_ADDRESS Matrix _mat;
  std::vector<_Indices> _true_indices;

public:

  /// @brief Construct a mask matrix view.
  template<matrix OtherMatrix>
    requires different_from<MaskMatrixView, OtherMatrix> &&
             std::constructible_from<Matrix, OtherMatrix>
  constexpr MaskMatrixView(OtherMatrix&& mat) // NOSONAR
      : _mat{std::forward<OtherMatrix>(mat)} {
    const auto cache_extents = [&](auto... extents) {
      const auto cache_indices = [&](auto... indices) {
        if (_mat(indices...)) _true_indices.push_back(_Indices{indices...});
      };
      detail::mdfor(cache_indices, extents...);
    };
    std::apply(cache_extents, shape());
  }

  /// @brief Get the matrix shape.
  constexpr auto shape() const {
    return _mat.shape();
  }

  /// @brief Get the matrix element at @p indices.
  template<class... Indices>
    requires compatible_matrix_indices_v<MaskMatrixView, Indices...>
  constexpr auto operator()(Indices... indices) const {
    STORM_ASSERT(in_range(shape(), indices...), "Indices are out of range!");
    return _mat(indices...);
  }

  /// @brief Get the range of indices,
  /// where the matrix contains a trueish elements.
  constexpr const auto& true_indices() const noexcept {
    return _true_indices;
  }

}; // class MaskMatrixView

template<class Matrix>
MaskMatrixView(Matrix&&) -> MaskMatrixView<matrix_view_t<Matrix>>;

// -----------------------------------------------------------------------------

/// @brief Index-like matrix coordate.
template<class Coord>
concept index_coord = detail::index_like<std::remove_cvref_t<Coord>>;

/// @brief Placeholder matrix coordinate.
template<class Coord>
concept placeholder_coord =
    std::same_as<std::remove_cvref_t<Coord>, Underscore>;

/// @brief New axis matrix coordinate.
template<class Coord>
concept new_axis_coord = std::same_as<std::remove_cvref_t<Coord>, NewAxis>;

/// @brief Range matrix coordinate:
/// sized random-access range with index-like values.
template<class Coord>
concept range_coord = std::ranges::sized_range<Coord> &&
                      std::ranges::random_access_range<Coord> &&
                      detail::index_like<std::ranges::range_value_t<Coord>>;

/// @brief Mask matrix coordinate: matrix with boolean-testable elements.
template<class Coord>
concept mask_coord = matrix<Coord> && //
                     boolean_testable<matrix_element_t<Coord>>;

/// @brief Simple matrix coordinate: index like or a placeholder.
/// This is a coordinate type that would not be wrapped into a view.
template<class Coord>
concept simple_coord =
    index_coord<Coord> || placeholder_coord<Coord> || new_axis_coord<Coord>;

/// @brief Matrix coordinate: either is a simple coordinate,
/// either a random-access range with index-like values,
/// either a matrix with boolean-testable elements.
template<class Coord>
concept coord = simple_coord<Coord> || range_coord<Coord> || mask_coord<Coord>;

/// @brief Matrix coordinate rank (how much indices it eats).
/// @{
template<coord Coord>
inline constexpr size_t coord_rank_v = [] {
  if constexpr (new_axis_coord<Coord>) return 0_sz;
  else if constexpr (mask_coord<Coord>) return matrix_rank_v<Coord>;
  else return 1_sz;
}();
template<coord... Coords>
inline constexpr size_t coords_rank_v = (0 + ... + coord_rank_v<Coords>);
/// @}

/// @brief Matrix coordinate that can be safely casted
/// into a matrix coordinate view.
template<class Coord>
concept viewable_coord =
    simple_coord<Coord> ||
    (range_coord<Coord> && std::ranges::viewable_range<Coord>) ||
    (mask_coord<Coord> && viewable_matrix<Coord>);

/// @brief Matrix coordinate view.
template<class Coord>
concept coord_view = coord<Coord> && std::movable<Coord> &&
                     ((simple_coord<Coord> && std::is_object_v<Coord>) ||
                      (range_coord<Coord> && std::ranges::view<Coord>) ||
                      (mask_coord<Coord> && matrix_view<Coord>) );

/// @brief Wrap the viewable coordinate @p coord into a coordinate view.
template<viewable_coord Coord>
constexpr auto to_coord_view(Coord&& coord) {
  if constexpr (coord_view<std::remove_cvref_t<Coord>>) {
    return std::forward<Coord>(coord);
  } else if constexpr (std::ranges::viewable_range<Coord>) {
    return std::views::all(std::forward<Coord>(coord));
  } else if constexpr (viewable_matrix<Coord>) {
    return MaskMatrixView{std::forward<Coord>(coord)};
  } else {
    static_assert(detail::always_false<Coord>, "Unexpected coordinate type!");
  }
}

/// @brief Suitable coordinate view type for a viewable coordinate.
template<viewable_coord Coord>
using coord_view_t = decltype(to_coord_view(std::declval<Coord>()));

// -----------------------------------------------------------------------------

/// @brief Matrix part view.
template<matrix_view Matrix, coord_view... Coords>
  requires (matrix_rank_v<Matrix> == (coord_rank_v<Coords> + ...))
class MatrixPartView final :
    public MatrixViewInterface<MatrixPartView<Matrix, Coords...>> {
private:

  STORM_NO_UNIQUE_ADDRESS Matrix _mat;
  STORM_NO_UNIQUE_ADDRESS std::tuple<Coords...> _coords;

public:

  /// @brief Construct a matrix part view.
  constexpr MatrixPartView(Matrix mat, Coords... coords)
      : _mat{std::move(mat)}, _coords{std::move(coords)...} {
    const auto check_coords = [&](auto... extents) {
      _check_coords(extents...);
    };
    std::apply(check_coords, _mat.shape());
  }

  /// @brief Get the matrix shape.
  constexpr auto shape() const {
    const auto map_extents = [&](auto... extents) {
      return _map_extents(extents...);
    };
    return std::apply(map_extents, _mat.shape());
  }

  /// @brief Get the matrix element at @p indices.
  /// @{
  template<class... Indices>
    requires compatible_matrix_indices_v<MatrixPartView, Indices...>
  constexpr decltype(auto) operator()(Indices... indices) {
    STORM_ASSERT(in_range(shape(), indices...), "Indices are out of range!");
    return std::apply(_mat, _map_indices(indices...));
  }
  template<class... Indices>
    requires compatible_matrix_indices_v<MatrixPartView, Indices...>
  constexpr decltype(auto) operator()(Indices... indices) const {
    STORM_ASSERT(in_range(shape(), indices...), "Indices are out of range!");
    return std::apply(_mat, _map_indices(indices...));
  }
  /// @}

private:

  static constexpr size_t _LastCoordIndex = sizeof...(Coords) - 1;

  // Check that coords are in range.
  template<size_t CoordIndex = 0>
  constexpr void _check_coords(auto... extents) const {
    const auto recure = [&](auto... remaining) {
      if constexpr (CoordIndex < _LastCoordIndex) {
        _check_coords<CoordIndex + 1>(remaining...);
      }
    };

    const auto& coord = std::get<CoordIndex>(_coords);
    using Coord = std::remove_cvref_t<decltype(coord)>;
    if constexpr (index_coord<Coord>) {
      [&](auto current_extent, auto... rest_extents) {
        const auto extent_sz = static_cast<size_t>(current_extent);
        STORM_ASSERT(static_cast<size_t>(coord) < extent_sz,
                     "Index-like coordinate is out of range!");
        recure(rest_extents...);
      }(extents...);
    } else if constexpr (new_axis_coord<Coord>) {
      recure(extents...);
    } else if constexpr (range_coord<Coord>) {
      [&](auto current_extent, auto... rest_extents) {
        const auto max_coord_element_sz = *std::ranges::max_element(
            coord, {}, [](auto value) { return static_cast<size_t>(value); });
        const auto extent_sz = static_cast<size_t>(current_extent);
        STORM_ASSERT(max_coord_element_sz < extent_sz,
                     "Range-like coordinate maximum element is out of range!");
        recure(rest_extents...);
      }(extents...);
    } else if constexpr (mask_coord<Coord>) {
      const auto check_mask = [&](auto current_extents, auto... rest_extents) {
        STORM_ASSERT(coord.shape() == current_extents,
                     "Mask-like coordinate shape is mismatched!");
        recure(rest_extents...);
      };
      constexpr size_t CoordRank = matrix_rank_v<Coord>;
      pack_n<CoordRank>(check_mask, extents...);
    }
  }

  // Map the original matrix extents to the view extents.
  template<size_t CoordIndex = 0>
  constexpr auto _map_extents(auto... extents) const {
    const auto recure = [&](auto done, auto... remaining) {
      const auto _done = to_tuple(done);
      if constexpr (CoordIndex == _LastCoordIndex) {
        return _done;
      } else {
        const auto _remaining = _map_extents<CoordIndex + 1>(remaining...);
        return std::tuple_cat(_done, _remaining);
      }
    };

    const auto& coord = std::get<CoordIndex>(_coords);
    using Coord = std::remove_cvref_t<decltype(coord)>;
    if constexpr (index_coord<Coord>) {
      return [&](auto /*current_extent*/, auto... rest_extents) {
        return recure(std::tuple{}, rest_extents...);
      }(extents...);
    } else if constexpr (placeholder_coord<Coord>) {
      return recure(extents...);
    } else if constexpr (new_axis_coord<Coord>) {
      return recure(fixed_extent_t<1>{}, extents...);
    } else if constexpr (range_coord<Coord>) {
      return [&](auto /*current_extent*/, auto... rest_extents) {
        return recure(std::ranges::size(coord), rest_extents...);
      }(extents...);
    } else if constexpr (mask_coord<Coord>) {
      const auto map_mask = [&](auto... rest_extents) {
        return recure(std::ranges::size(coord.true_indices()), rest_extents...);
      };
      constexpr size_t CoordRank = matrix_rank_v<Coord>;
      return drop_n<CoordRank>(map_mask, extents...);
    } else {
      static_assert(detail::always_false<Coord>, "Unexpected coordinate type!");
    }
  }

  // Map the view indices to the original matrix indices.
  template<size_t CoordIndex = 0>
  constexpr auto _map_indices(auto... indices) const {
    const auto recure = [&](auto done, auto... remaining) {
      const auto _done = to_tuple(done);
      if constexpr (CoordIndex == _LastCoordIndex) {
        return _done;
      } else {
        const auto _remaining = _map_indices<CoordIndex + 1>(remaining...);
        return std::tuple_cat(_done, _remaining);
      }
    };

    const auto& coord = std::get<CoordIndex>(_coords);
    using Coord = std::remove_cvref_t<decltype(coord)>;
    if constexpr (index_coord<Coord>) {
      return recure(coord, indices...);
    } else if constexpr (placeholder_coord<Coord>) {
      return recure(indices...);
    } else if constexpr (new_axis_coord<Coord>) {
      return [&](auto /*current_index*/, auto... rest_indices) {
        return recure(std::tuple{}, rest_indices...);
      }(indices...);
    } else if constexpr (range_coord<Coord>) {
      return [&](auto current_index, auto... rest_indices) {
        return recure(coord[current_index], rest_indices...);
      }(indices...);
    } else if constexpr (mask_coord<Coord>) {
      return [&](auto current_index, auto... rest_indices) {
        return recure(coord.true_indices()[current_index], rest_indices...);
      }(indices...);
    } else {
      static_assert(detail::always_false<Coord>, "Unexpected coordinate type!");
    }
  }

}; // class MatrixPartView

template<class Matrix, class... Coords>
MatrixPartView(Matrix&&, Coords&&...)
    -> MatrixPartView<matrix_view_t<Matrix>, coord_view_t<Coords>...>;

// -----------------------------------------------------------------------------

/// @brief Take a part of the matrix @p mat at coordinates @p coords.
template<viewable_matrix Matrix, viewable_coord... Coords>
  requires (matrix_rank_v<Matrix> >= coords_rank_v<Coords...>)
constexpr auto at(Matrix&& mat, Coords&&... coords) {
  constexpr size_t Delta = matrix_rank_v<Matrix> - coords_rank_v<Coords...>;
  if constexpr (Delta == 0) {
    if constexpr ((index_coord<Coords> && ...)) {
      return mat(std::forward<Coords>(coords)...);
    } else {
      return MatrixPartView{std::forward<Matrix>(mat),
                            std::forward<Coords>(coords)...};
    }
  } else {
    const auto make_view = [&](auto... _s) {
      return MatrixPartView{std::forward<Matrix>(mat),
                            std::forward<Coords>(coords)..., _s...};
    };
    return n_times<Delta>(make_view, _);
  }
}

// -----------------------------------------------------------------------------

} // namespace Storm
