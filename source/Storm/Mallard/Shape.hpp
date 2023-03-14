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

#include <Storm/Mallard/Fwd.hpp>

#include <algorithm>
#include <array>
#include <cstdint>
#include <initializer_list>
#include <ranges>
#include <tuple>
#include <vector>

namespace Storm::shapes
{

// -----------------------------------------------------------------------------

/// @brief Shape type.
enum class Type : std::uint8_t {
  segment,        ///< Segment shape type.
  triangle,       ///< Triangle shape type.
  triangle_strip, ///< Triangle strip shape type.
  quadrangle,     ///< Quadrangle shape type.
  polygon,        ///< Polygon shape type.
  tetrahedron,    ///< Tetrahedron shape type.
  pyramid,        ///< Pyramid shape type.
  pentahedron,    ///< Pentahedron shape type.
  hexahedron,     ///< Hexahedron shape type.
  count_,
}; // enum class shape_type

/// @brief Shape concept.
template<class Shape>
concept shape =
    requires(Shape shape) {
      // clang-format off
      { shape.type() } -> std::same_as<Type>;
      { shape.nodes() } -> std::ranges::forward_range;
      // clang-format off
    } &&
    std::same_as<
        std::ranges::range_value_t<decltype(std::declval<Shape>().nodes())>,
        NodeIndex>;

/// @brief 1D Shape concept.
template<class Shape>
concept shape1D = shape<Shape> && //
                  !requires(Shape shape) { shape.edges(); } &&
                  !requires(Shape shape) { shape.faces(); };

/// @brief 2D Shape concept.
template<class Shape>
concept shape2D = shape<Shape> && //
                  requires(Shape shape) {
                    std::apply([](shape1D auto&&...) {}, shape.edges());
                  } && !requires(Shape shape) { shape.faces(); };

/// @brief 3D Shape concept.
template<class Shape>
concept shape3D = shape<Shape> && //
                  requires(Shape shape) {
                    std::apply([](shape1D auto&&...) {}, shape.edges());
                  } && //
                  requires(Shape shape) {
                    std::apply([](shape2D auto&&...) {}, shape.faces());
                  };

/// @brief Shape topological dimensionality.
template<shape Shape>
  requires (shape1D<Shape> || shape2D<Shape> || shape3D<Shape>)
inline constexpr size_t shape_dim_v = []() {
  if constexpr (shape1D<Shape>) return 1;
  if constexpr (shape2D<Shape>) return 2;
  if constexpr (shape3D<Shape>) return 3;
}();

/// @brief Get the shape parts.
template<size_t Index, shape Shape>
  requires ((Index == 0) ||
            (Index == 1 && (shape2D<Shape> || shape3D<Shape>) ) ||
            (Index == 2 && shape3D<Shape>) )
constexpr auto parts(const Shape& shape) noexcept
{
  if constexpr (Index == 0) {
    return shape.nodes();
  } else if constexpr (Index == 1) {
    return shape.edges();
  } else if constexpr (Index == 2) {
    return shape.faces();
  }
}

// -----------------------------------------------------------------------------

/// @brief Complex shape concept.
template<class Shape, class Mesh>
concept complex_shape =
    shape<Shape> && mesh<Mesh> &&
    requires(Shape shape, const Mesh& mesh) {
      // clang-format off
      { shape.pieces(mesh) } -> std::ranges::forward_range;
      // clang-format on
    } &&
    shape<std::ranges::range_value_t<decltype( //
        std::declval<Shape>().pieces(std::declval<Mesh>()))>>;

/// @brief Piece type of a complex.
template<class Shape, class Mesh>
  requires complex_shape<Shape, Mesh>
using piece_t = std::ranges::range_value_t<decltype( //
    std::declval<Shape>().pieces(std::declval<Mesh>()))>;

namespace detail_
{
  template<class Shape, class Mesh>
  concept can_volume_ =
      requires(Shape shape, const Mesh& mesh) { volume(shape, mesh); };
  template<class Shape, class Mesh>
  concept can_barycenter_ =
      requires(Shape shape, const Mesh& mesh) { barycenter(shape, mesh); };
  template<class Shape, class Mesh>
  concept can_normal_ =
      requires(Shape shape, const Mesh& mesh) { normal(shape, mesh); };
} // namespace detail_

/// @brief Compute the complex @p shape "volume"
/// (length in 1D, area in 2D, volume in 3D).
template<shape Shape, mesh Mesh>
  requires (complex_shape<Shape, Mesh> &&
            detail_::can_volume_<piece_t<Shape, Mesh>, Mesh>)
constexpr real_t volume(const Shape& shape, const Mesh& mesh) noexcept
{
  const auto pieces = shape.pieces(mesh);
  auto vol = volume(pieces.front(), mesh);
  for (const auto& piece : pieces | std::views::drop(1)) {
    vol += volume(piece, mesh);
  }
  return vol;
}

/// @brief Compute the simplex @p shape barycenter.
template<shape Shape, mesh Mesh>
  requires (!complex_shape<Shape, Mesh> &&
            (mesh_dim_v<Mesh> >= shape_dim_v<Shape>) )
constexpr mesh_vec_t<Mesh> barycenter(const Shape& shape,
                                      const Mesh& mesh) noexcept
{
  const auto nodes = shape.nodes();
  mesh_vec_t<Mesh> sum_position{mesh.position(nodes.front())};
  for (const auto& node : nodes | std::views::drop(1)) {
    sum_position += mesh.position(node);
  }
  return sum_position / nodes.size();
}

/// @brief Compute the complex @p shape barycenter.
template<shape Shape, mesh Mesh>
  requires (complex_shape<Shape, Mesh> &&
            detail_::can_volume_<piece_t<Shape, Mesh>, Mesh> &&
            detail_::can_barycenter_<piece_t<Shape, Mesh>, Mesh>)
constexpr mesh_vec_t<Mesh> barycenter(const Shape& shape,
                                      const Mesh& mesh) noexcept
{
  const auto pieces = shape.pieces(mesh);
  auto vol = volume(pieces.front(), mesh);
  auto vol_center = vol * barycenter(pieces.front(), mesh);
  for (const auto& piece : pieces | std::views::drop(1)) {
    const auto dv = volume(piece, mesh);
    vol += dv, vol_center += dv * barycenter(piece, mesh);
  }
  return vol_center / vol;
}

/// @brief Compute the complex @p shape normal.
template<shape Shape, mesh Mesh>
  requires (complex_shape<Shape, Mesh> &&
            detail_::can_volume_<piece_t<Shape, Mesh>, Mesh> &&
            detail_::can_normal_<piece_t<Shape, Mesh>, Mesh>)
constexpr mesh_vec_t<Mesh> normal(const Shape& shape, const Mesh& mesh) noexcept
{
  const auto pieces = shape.pieces(mesh);
  auto vol_normal = volume(pieces.front(), mesh) * normal(pieces.front(), mesh);
  for (const auto& piece : pieces | std::views::drop(1)) {
    vol_normal += volume(piece, mesh) * normal(piece, mesh);
  }
  return normalize(vol_normal);
}

// -----------------------------------------------------------------------------

/// @brief Segment shape.
/// @verbatim
/// +==========================================================================+
/// |                                           | Edges:                       |
/// |               n1 @ f1                     |------------------------------|
/// |                   \                       | e1 = (n1,n2)                 |
/// |                    \                      |                              |
/// |                     v e1                  |------------------------------|
/// |                      \                    | Faces:                       |
/// |                       \                   |------------------------------|
/// |                     n2 @ f2               | f1 = (n1)                    |
/// |                                           | f2 = (n2)                    |
/// +==========================================================================+
/// @endverbatim
class Seg final
{
public:

  /// @brief Segment node.
  /// @{
  NodeIndex n1{}, n2{};
  /// @}

  /// @brief Construct a segment.
  constexpr Seg() = default;
  /// @brief Construct a segment with nodes.
  constexpr Seg(NodeIndex i1, NodeIndex i2) noexcept : n1{i1}, n2{i2} {}

  /// @brief Segment type.
  constexpr static auto type() noexcept
  {
    return Type::segment;
  }

  /// @brief Segment nodes.
  constexpr auto nodes() const noexcept
  {
    return std::array{n1, n2};
  }

}; // class Seg

/// @brief Segment @p seg length ("volume").
template<mesh Mesh>
constexpr real_t volume(const Seg& seg, const Mesh& mesh) noexcept
{
  const auto v1{mesh.position(seg.n1)}, v2{mesh.position(seg.n2)};
  return length(v2 - v1);
}

/// @brief Segment @p seg normal.
template<mesh Mesh>
  requires (mesh_dim_v<Mesh> == 2)
constexpr mesh_vec_t<Mesh> normal(const Seg& seg, const Mesh& mesh) noexcept
{
  const auto v1{mesh.position(seg.n1)}, v2{mesh.position(seg.n2)};
  const auto d = normalize(v2 - v1);
  /// @todo Is sign here correct?
  /// @todo This looks incomplete:
  return -mesh_vec_t<Mesh>{-d(1, 0), d(0, 0)};
}

// -----------------------------------------------------------------------------

/// @brief Triangular shape.
/// @verbatim
/// +==========================================================================+
/// |                                         | Edges /                        |
/// |                    n3                   |      / Faces:                  |
/// |                    @                    |--------------------------------|
/// |                   / \                   | e1 = f1 = (n1,n2)              |
/// |                  /   \                  | e2 = f2 = (n2,n3)              |
/// |           e3/f3 v     ^ e2/f2           | e3 = f3 = (n3,n1)              |
/// |                /       \                |                                |
/// |               /         \               |                                |
/// |           n1 @----->-----@ n2           |                                |
/// |                  e1/f1                  |                                |
/// |                                         |                                |
/// +==========================================================================+
/// @endverbatim
class Triangle final
{
public:

  /// @brief Triangle node.
  /// @{
  NodeIndex n1{}, n2{}, n3{};
  /// @}

  /// @brief Construct a triangle.
  constexpr Triangle() = default;
  /// @brief Construct a triangle with nodes.
  constexpr Triangle(NodeIndex i1, NodeIndex i2, NodeIndex i3) noexcept
      : n1{i1}, n2{i2}, n3{i3}
  {
  }

  /// @brief Triangle type.
  constexpr static auto type() noexcept
  {
    return Type::triangle;
  }

  /// @brief Triangle nodes.
  constexpr auto nodes() const noexcept
  {
    return std::array{n1, n2, n3};
  }

  /// @brief Triangle edges.
  constexpr auto edges() const noexcept
  {
    return std::array{Seg{n1, n2}, Seg{n2, n3}, Seg{n3, n1}};
  }

}; // class Triangle

/// @brief Triangle @p tri "volume" (area).
template<mesh Mesh>
  requires (mesh_dim_v<Mesh> >= 2)
constexpr real_t volume(const Triangle& tri, const Mesh& mesh) noexcept
{
  const auto v1{mesh.position(tri.n1)}, v2{mesh.position(tri.n2)};
  const auto v3{mesh.position(tri.n3)};
  if constexpr (mesh_dim_v<Mesh> == 2) {
    // const glm::dmat2 d{v2 - v1, v3 - v1};
    // return glm::abs(glm::determinant(d)) / 2.0;
    /// @todo This looks incomplete:
    auto d0 = v2 - v1, d1 = v3 - v1;
    return 0.5 * abs(d0(0, 0) * d1(1, 0) - d0(1, 0) * d1(0, 0));
  } else {
    return length(cross_product(v2 - v1, v3 - v1)) / 2.0;
  }
}

/// @brief Triangle @p tri normal.
template<mesh Mesh>
  requires (mesh_dim_v<Mesh> == 3)
constexpr mesh_vec_t<Mesh> normal(const Triangle& tri,
                                  const Mesh& mesh) noexcept
{
  const auto v1{mesh.position(tri.n1)}, v2{mesh.position(tri.n2)};
  const auto v3{mesh.position(tri.n3)};
  return normalize(cross(v2 - v1, v3 - v1));
}

// -----------------------------------------------------------------------------

/// @brief Quadrangular shape.
/// @verbatim
/// +==========================================================================+
/// |                                                  | Edges /               |
/// |                         e4/f4                    |      / Faces:         |
/// |                 n1 @------<------@ n4            |-----------------------|
/// |                   /             /                | e1 = f1 = (n1,n2)     |
/// |            e1/f1 v             /                 | e2 = f2 = (n2,n3)     |
/// |                 /             ^ e3/f3            | e3 = f3 = (n3,n4)     |
/// |                /             /                   | e4 = f4 = (n4,n1)     |
/// |            n2 @------>------@ n3                 |                       |
/// |                   e2/f2                          |                       |
/// |                                                  |                       |
/// +==========================================================================+
/// @endverbatim
class Quadrangle final
{
public:

  /// @brief Quadrangle node.
  /// @{
  NodeIndex n1{}, n2{}, n3{}, n4{};
  /// @}

  /// @brief Construct a quadrangle.
  constexpr Quadrangle() = default;
  /// @brief Construct a quadrangle with nodes.
  constexpr Quadrangle(NodeIndex i1, NodeIndex i2, //
                       NodeIndex i3, NodeIndex i4) noexcept
      : n1{i1}, n2{i2}, n3{i3}, n4{i4}
  {
  }

  /// @brief Quadrangle type.
  constexpr static auto type() noexcept
  {
    return Type::quadrangle;
  }

  /// @brief Quadrangle nodes.
  constexpr auto nodes() const noexcept
  {
    return std::array{n1, n2, n3, n4};
  }

  /// @brief Quadrangle edges.
  constexpr auto edges() const noexcept
  {
    return std::array{Seg{n1, n2}, Seg{n2, n3}, Seg{n3, n4}, Seg{n4, n1}};
  }

  /// @brief Quadrangle pieces.
  template<mesh Mesh>
    requires (mesh_dim_v<Mesh> >= 2)
  constexpr auto pieces(const Mesh& mesh) const noexcept
  {
    static_cast<void>(mesh); /// @todo Convexity check!
    return std::array{Triangle{n1, n2, n3}, Triangle{n3, n4, n1}};
    // return std::array{Triangle{n1, n2, n4}, Triangle{n2, n3, n4}};
  }

}; // class Quadrangle

// -----------------------------------------------------------------------------

/// @brief Triangle strip shape.
/// @verbatim
/// +==========================================================================+
/// |                                        |                                 |
/// |         n2      n4      n6      n8     |                                 |
/// |         @---<---@---<---@---<---@      |                                 |
/// |        / \     / \     / \     /       |                                 |
/// |       v   \   /   \   /   \   ^        |          INCOMPLETE..           |
/// |      /     \ /     \ /     \ /         |                                 |
/// |     @--->---@--->---@--->---@          |                                 |
/// |     n1      n3      n5      n9         |                                 |
/// |                                        |                                 |
/// +==========================================================================+
/// @endverbatim
class TriangleStrip final
{
public:

  /// @brief Triangle strip nodes.
  std::vector<NodeIndex> n{};

  /// @brief Construct a triangle strip.
  constexpr TriangleStrip() = default;
  /// @brief Construct a triangle strip with nodes.
  template<class... NodeIndices>
  constexpr TriangleStrip(NodeIndices... i) : n{i...}
  {
    STORM_ASSERT_(n.size() >= 3, "Triangle strip must have least 3 nodes!");
  }

  /// @brief Triangle strip type.
  /// May fallback to triangle or quadrangle types.
  constexpr auto type() const noexcept
  {
    STORM_ASSERT_(n.size() >= 3, "Triangle strip must have least 3 nodes!");
    switch (n.size()) {
      case 3: return Type::triangle;
      case 4: return Type::quadrangle;
      default: return Type::triangle_strip;
    }
  }

  /// @brief Triangle strip nodes.
  constexpr const auto& nodes() const noexcept
  {
    return n;
  }

  /// @brief Triangle strip edges.
  template<class T = void>
  constexpr auto edges() const noexcept
  {
    static_assert(meta::always_false<T>,
                  "TriangleStrip::edges is not implemented yet!");
    return std::vector<Seg>{};
  }

  /// @brief Triangle strip pieces.
  template<mesh Mesh>
  constexpr auto pieces(const Mesh& mesh) const noexcept
  {
    static_cast<void>(mesh); /// @todo Convexity check!
    return std::views::iota(size_t{2}, n.size()) |
           std::views::transform([&](size_t i) {
             return Triangle{n[i - 2], n[i - 1], n[i]};
           });
  }

}; // class TriangleStrip

// -----------------------------------------------------------------------------

/// @brief Arbitrary polygon shape.
/// @verbatim
/// +==========================================================================+
/// |                                      |                                   |
/// |                  n5  e4/f4           |                                   |
/// |                 .@---<--@ n4         |                                   |
/// |               ./`      /             |                                   |
/// |        e5/f5 v        ^ e3/f3        |                                   |
/// |           ./`        /               |                                   |
/// |        n6 @         @ n3             |           INCOMPLETE..            |
/// |           |          \               |                                   |
/// |     e6/f6 v           ^ e2/f2        |                                   |
/// |           |            \             |                                   |
/// |           @------>------@            |                                   |
/// |           n1     |     n2            |                                   |
/// |                e1/f1                 |                                   |
/// |                                      |                                   |
/// +==========================================================================+
/// @endverbatim
class Polygon final
{
public:

  /// @brief Polygon nodes.
  std::vector<NodeIndex> n{};

  /// @brief Construct a polygon.
  constexpr Polygon() = default;
  /// @brief Construct a polygon with nodes.
  template<class... NodeIndices>
  constexpr Polygon(NodeIndices... i) : n{i...}
  {
    STORM_ASSERT_(n.size() >= 3, "Polygon must have least 3 nodes!");
  }

  /// @brief Polygon type.
  /// May fallback to triangle or quadrangle types.
  constexpr auto type() const noexcept
  {
    STORM_ASSERT_(n.size() >= 3, "Polygon must have least 3 nodes!");
    switch (n.size()) {
      case 3: return Type::triangle;
      case 4: return Type::quadrangle;
      default: return Type::polygon;
    }
  }

  /// @brief Polygon nodes.
  constexpr const auto& nodes() const noexcept
  {
    return n;
  }

  /// @brief Polygon edges.
  constexpr auto edges() const noexcept
  {
    return std::views::iota(size_t{1}, n.size()) |
           std::views::transform([&](size_t i) {
             return Seg{n[i - 1], n[i]};
           });
  }

  /// @brief Polygon pieces.
  template<mesh Mesh>
  constexpr auto pieces(const Mesh& mesh) const
  {
    static_assert(meta::always_false<Mesh>,
                  "Polygon::pieces is not implemented yet!");
    static_cast<void>(mesh); /// @todo Convexity check!
    return std::vector<Triangle>{};
  }

}; // class Polygon

// -----------------------------------------------------------------------------

/// @brief Tetrahedral shape.
/// @verbatim
/// +==========================================================================+
/// |                                                 | Edges:                 |
/// |                         f4                      |------------------------|
/// |                    n4   ^                       | e1 = (n1,n2)           |
/// |                     @   |                       | e2 = (n2,n3)           |
/// |              f2    /|\. |     f3                | e3 = (n3,n1)           |
/// |              ^    / | `\.     ^                 | e4 = (n1,n4)           |
/// |               \  /  |   `\.  /                  | e5 = (n2,n4)           |
/// |                \`   |   | `\/                   | e6 = (n3,n4)           |
/// |                ,\   |   o  /`\                  |                        |
/// |            e4 ^  *  |     *   `^.e6             |                        |
/// |              /   e5 ^           `\              |------------------------|
/// |          n1 @-------|--<----------@ n3          | Faces:                 |
/// |              \      |  e3       ,/              |------------------------|
/// |               \     |     o   ,/`               | f1 = (n1,n3,n2)        |
/// |                \    ^ e5  | ,/`                 | f2 = (n1,n2,n4)        |
/// |              e1 v   |     ,^ e2                 | f3 = (n2,n3,n4)        |
/// |                  \  |   ,/`                     | f4 = (n3,n1,n4)        |
/// |                   \ | ,/` |                     |                        |
/// |                    \|/`   |                     |                        |
/// |                     @ n2  v                     |                        |
/// |                           f1                    |                        |
/// |                                                 |                        |
/// +==========================================================================+
/// @endverbatim
class Tetrahedron final
{
public:

  /// @brief Tetrahedron node.
  /// @{
  NodeIndex n1{}, n2{}, n3{}, n4{};
  /// @}

  /// @brief Construct a tetrahedron.
  constexpr Tetrahedron() = default;
  /// @brief Construct a tetrahedron with nodes.
  constexpr Tetrahedron(NodeIndex i1, NodeIndex i2, //
                        NodeIndex i3, NodeIndex i4) noexcept
      : n1{i1}, n2{i2}, n3{i3}, n4{i4}
  {
  }

  /// @brief Tetrahedron type.
  constexpr static auto type() noexcept
  {
    return Type::tetrahedron;
  }

  /// @brief Tetrahedron nodes.
  constexpr auto nodes() const noexcept
  {
    return std::array{n1, n2, n3, n4};
  }

  /// @brief Tetrahedron edges.
  constexpr auto edges() const noexcept
  {
    return std::array{Seg{n1, n2}, Seg{n2, n3}, Seg{n3, n1},
                      Seg{n1, n4}, Seg{n2, n4}, Seg{n3, n4}};
  }

  /// @brief Tetrahedron faces.
  constexpr auto faces() const noexcept
  {
    return std::array{Triangle{n1, n3, n2}, Triangle{n1, n2, n4},
                      Triangle{n2, n3, n4}, Triangle{n3, n1, n4}};
  }

}; // class Tetrahedron

/// @brief Tetrahedron @p tet barycenter.
template<mesh Mesh>
  requires (mesh_dim_v<Mesh> >= 3)
constexpr mesh_vec_t<Mesh> barycenter(const Tetrahedron& tet,
                                      const Mesh& mesh) noexcept
{
  const auto v1{mesh.position(tet.n1)}, v2{mesh.position(tet.n2)};
  const auto v3{mesh.position(tet.n3)}, v4{mesh.position(tet.n4)};
  return (v1 + v2 + v3 + v4) / 4.0;
}

// -----------------------------------------------------------------------------

/// @brief Pyramidal shape.
/// @verbatim
/// +==========================================================================+
/// |                                                     | Edges:             |
/// |                                                     |--------------------|
/// |                               n5                    | e1 = (n1,n2)       |
/// |                 f4           ,@                     | e2 = (n2,n3)       |
/// |                  ^        ,/`/|\     f2             | e3 = (n3,n4)       |
/// |                   \    ,/`  / | \    ^              | e4 = (n4,n1)       |
/// |                    \,/`    /  |  \  /               | e5 = (n1,n5)       |
/// |               e8 ,/`\     /   |   \/                | e6 = (n2,n5)       |
/// |               ,^`    o   /    |   /\                | e7 = (n3,n5)       |
/// |            ,/`          /     |  *  \               | e8 = (n4,n5)       |
/// | f5 <------------*      /   e7 ^   o--\---------> f3 |                    |
/// |      ,/`              /       |       \             |--------------------|
/// |  n4 @-----<----------/--------@  n3    ^ e6         | Faces:             |
/// |      `\.  e3        /          `\.      \           |--------------------|
/// |         `>.        ^ e5           `\. e2 \          | f1 = (n1,n4,n3,n2) |
/// |         e4 `\.    /       o          `<.  \         | f2 = (n1,n2,n5)    |
/// |               `\./        |             `\.\        | f3 = (n2,n3,n5)    |
/// |              n1 @-------------------->-----@ n2     | f4 = (n3,n4,n5)    |
/// |                           |          e1             | f5 = (n4,n1,n5)    |
/// |                           |                         |                    |
/// |                           v                         |                    |
/// |                           f1                        |                    |
/// |                                                     |                    |
/// +==========================================================================+
/// @endverbatim
class Pyramid final
{
public:

  /// @brief Pyramid node.
  /// @{
  NodeIndex n1{}, n2{}, n3{}, n4{}, n5{};
  /// @}

  /// @brief Construct a tetrahedron.
  constexpr Pyramid() = default;
  /// @brief Construct a tetrahedron with nodes.
  constexpr Pyramid(NodeIndex i1, NodeIndex i2, //
                    NodeIndex i3, NodeIndex i4, NodeIndex i5) noexcept
      : n1{i1}, n2{i2}, n3{i3}, n4{i4}, n5{i5}
  {
  }

  /// @brief Pyramid type.
  constexpr static auto type() noexcept
  {
    return Type::pyramid;
  }

  /// @brief Pyramid nodes.
  constexpr auto nodes() const noexcept
  {
    return std::array{n1, n2, n3, n4, n5};
  }

  /// @brief Pyramid edges.
  constexpr auto edges() const noexcept
  {
    return std::array{Seg{n1, n2}, Seg{n2, n3}, Seg{n3, n4}, Seg{n4, n1},
                      Seg{n1, n5}, Seg{n2, n5}, Seg{n3, n5}, Seg{n4, n5}};
  }

  /// @brief Pyramid faces.
  constexpr auto faces() const noexcept
  {
    return std::tuple{Quadrangle{n1, n4, n3, n2}, Triangle{n1, n2, n5},
                      Triangle{n2, n3, n5}, Triangle{n3, n4, n5},
                      Triangle{n4, n1, n5}};
  }

  /// @brief Pyramid pieces.
  template<mesh Mesh>
    requires (mesh_dim_v<Mesh> >= 3)
  constexpr auto pieces(const Mesh& mesh) const noexcept
  {
    static_cast<void>(mesh); /// @todo Convexity check!
    return std::array{Tetrahedron{n1, n2, n3, n5}, Tetrahedron{n3, n4, n1, n5}};
  }

}; // class Pyramid

// -----------------------------------------------------------------------------

/// @brief Pentahedral shape (triangular prism).
/// @verbatim
/// +==========================================================================+
/// |                                            | Edges:                      |
/// |                  f5                        |-----------------------------|
/// |                  ^  f3                     | e1 = (n1,n2)                |
/// |                  |  ^                      | e2 = (n2,n3)                |
/// |              e9  |  |                      | e3 = (n3,n1)                |
/// |       n4 @---<---|-------------@ n6        | e4 = (n1,n4)                |
/// |          |\      *  |        ,/|           | e5 = (n2,n5)                |
/// |          | \        o      ,/` |           | e6 = (n3,n6)                |
/// |          |  \         e8 ,^`   |           | e7 = (n4,n5)                |
/// |          |   v e7      ,/`     |           | e8 = (n5,n6)                |
/// |       e4 ^    \      ,/`       ^ e6        | e9 = (n6,n4)                |
/// |          |     \   ,/`         |           |                             |
/// |          |      \ /`        *-------> f2   |-----------------------------|
/// |   f1 <-------*   @ n5          |           | Faces:                      |
/// |          |       |             |           |-----------------------------|
/// |       n1 @-------|---------<---@ n3        | f1 = (n1,n2,n5,n4)          |
/// |           \      |        e3 ,/            | f2 = (n2,n3,n6,n5)          |
/// |            \     |     o   ,/`             | f3 = (n3,n1,n4,n6)          |
/// |             \    ^ e5  | ,/`               | f4 = (n1,n3,n2)             |
/// |           e1 v   |     ,^ e2               | f5 = (n4,n5,n6)             |
/// |               \  |   ,/|                   |                             |
/// |                \ | ,/` |                   |                             |
/// |                 \|/`   |                   |                             |
/// |                  @ n2  v                   |                             |
/// |                        f4                  |                             |
/// |                                            |                             |
/// +==========================================================================+
/// @endverbatim
class Pentahedron final
{
public:

  /// @brief Pentahedron node.
  /// @{
  NodeIndex n1{}, n2{}, n3{}, n4{}, n5{}, n6{};
  /// @}

  /// @brief Construct a pentahedron.
  constexpr Pentahedron() = default;
  /// @brief Construct a pentahedron with nodes.
  constexpr Pentahedron(NodeIndex i1, NodeIndex i2, //
                        NodeIndex i3, NodeIndex i4, //
                        NodeIndex i5, NodeIndex i6) noexcept
      : n1{i1}, n2{i2}, n3{i3}, n4{i4}, n5{i5}, n6{i6}
  {
  }

  /// @brief Pentahedron type.
  constexpr static auto type() noexcept
  {
    return Type::pentahedron;
  }

  /// @brief Pentahedron nodes.
  constexpr auto nodes() const noexcept
  {
    return std::array{n1, n2, n3, n4, n5, n6};
  }

  /// @brief Pentahedron edges.
  constexpr auto edges() const noexcept
  {
    return std::array{Seg{n1, n2}, Seg{n2, n3}, Seg{n3, n1},
                      Seg{n1, n4}, Seg{n2, n5}, Seg{n3, n6},
                      Seg{n4, n5}, Seg{n5, n6}, Seg{n6, n4}};
  }

  /// @brief Pentahedron faces.
  constexpr auto faces() const noexcept
  {
    return std::tuple{Quadrangle{n1, n2, n5, n4}, Quadrangle{n2, n3, n6, n5},
                      Quadrangle{n3, n1, n4, n6}, Triangle{n1, n3, n2},
                      Triangle{n4, n5, n6}};
  }

  /// @brief Pentahedron pieces.
  template<mesh Mesh>
    requires (mesh_dim_v<Mesh> >= 3)
  constexpr auto pieces(const Mesh& mesh) const noexcept
  {
    static_cast<void>(mesh); /// @todo Convexity check!
    return std::array{Tetrahedron{n1, n2, n3, n5}, Tetrahedron{n3, n1, n4, n5},
                      Tetrahedron{n4, n6, n3, n5}};
  }

}; // class Pentahedron

// -----------------------------------------------------------------------------

/// @brief Hexahedral shape.
/// @verbatim
/// +==========================================================================+
/// |                                               | Edges:                   |
/// |                       f6                      |--------------------------|
/// |                       ^       f3              |  e1 = (n1,n2)            |
/// |                       |       ^               |  e2 = (n2,n3)            |
/// |                   e10 |      /                |  e3 = (n3,n4)            |
/// |             n7 @---<--|----------@ n6         |  e4 = (n4,n1)            |
/// |               /|      |    /    /|            |  e5 = (n1,n5)            |
/// |              / |      |   o    / |            |  e6 = (n2,n6)            |
/// |         e11 v  |      *    e9 ^  ^ e6         |  e7 = (n3,n7)            |
/// |            /   ^ e7          /   |            |  e8 = (n4,n8)            |
/// |           /    |      e12   /  *-------> f2   |  e9 = (n5,n6)            |
/// |       n8 @------------->---@ n5  |            | e10 = (n6,n7)            |
/// |   f4 <---|--o  |           |     |            | e11 = (n7,n8)            |
/// |          |  n3 @---<-------|-----@ n2         | e12 = (n8,n5)            |
/// |          |    /    e2      |    /             |                          |
/// |       e8 ^   /          e5 ^   /              |--------------------------|
/// |          |  v e3  *        |  ^ e1            | Faces:                   |
/// |          | /     /    o    | /                |--------------------------|
/// |          |/     /     |    |/                 | f1 = (n1,n4,n3,n2)       |
/// |       n4 @-----/-->--------@ n1               | f2 = (n1,n2,n6,n5)       |
/// |               /   e4  |                       | f3 = (n2,n3,n7,n6)       |
/// |              /        v                       | f4 = (n3,n4,n8,n7)       |
/// |             v         f1                      | f5 = (n1,n5,n8,n4)       |
/// |             f5                                | f6 = (n5,n6,n7,n8)       |
/// |                                               |                          |
/// |                                               |                          |
/// +==========================================================================+
/// @endverbatim
class Hexahedron final
{
public:

  /// @brief Hexahedron node.
  /// @{
  NodeIndex n1{}, n2{}, n3{}, n4{}, n5{}, n6{}, n7{}, n8{};
  /// @}

  /// @brief Construct a hexahedron.
  constexpr Hexahedron() = default;
  /// @brief Construct a hexahedron with nodes.
  constexpr Hexahedron(NodeIndex i1, NodeIndex i2, //
                       NodeIndex i3, NodeIndex i4, //
                       NodeIndex i5, NodeIndex i6, //
                       NodeIndex i7, NodeIndex i8) noexcept
      : n1{i1}, n2{i2}, n3{i3}, n4{i4}, n5{i5}, n6{i6}, n7{i7}, n8{i8}
  {
  }

  /// @brief Hexahedron type.
  constexpr static auto type() noexcept
  {
    return Type::hexahedron;
  }

  /// @brief Hexahedron nodes.
  constexpr auto nodes() const noexcept
  {
    return std::array{n1, n2, n3, n4, n5, n6, n7, n8};
  }

  /// @brief Hexahedron edges.
  constexpr auto edges() const noexcept
  {
    return std::array{Seg{n1, n2}, Seg{n2, n3}, Seg{n3, n4}, Seg{n4, n1},
                      Seg{n1, n5}, Seg{n2, n6}, Seg{n3, n7}, Seg{n4, n8},
                      Seg{n5, n6}, Seg{n6, n7}, Seg{n7, n8}, Seg{n8, n5}};
  }

  /// @brief Hexahedron faces.
  constexpr auto faces() const noexcept
  {
    return std::array{Quadrangle{n1, n4, n3, n2}, Quadrangle{n1, n2, n6, n5},
                      Quadrangle{n2, n3, n7, n6}, Quadrangle{n3, n4, n8, n7},
                      Quadrangle{n1, n5, n8, n4}, Quadrangle{n5, n6, n7, n8}};
  }

  /// @brief Hexahedron pieces.
  template<mesh Mesh>
    requires (mesh_dim_v<Mesh> >= 3)
  constexpr auto pieces(const Mesh& mesh) const noexcept
  {
    static_cast<void>(mesh); /// @todo Convexity check!
    return std::array{Tetrahedron{n1, n4, n2, n5}, Tetrahedron{n4, n3, n2, n7},
                      Tetrahedron{n5, n6, n7, n2}, Tetrahedron{n5, n7, n8, n4},
                      Tetrahedron{n5, n4, n2, n7}};
  }

}; // class Hexahedron

// -----------------------------------------------------------------------------

} // namespace Storm::shapes
