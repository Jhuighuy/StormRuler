/// Copyright (C) 2022 Oleg Butakov
///
/// Permission is hereby granted, free of charge, to any person obtaining a copy
/// of this software and associated documentation files (the "Software"), to
/// deal in the Software without restriction, including without limitation the
/// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
/// sell copies of the Software, and to permit persons to whom the Software is
/// furnished to do so, subject to the following conditions:
///
/// The above copyright notice and this permission notice shall be included in
/// all copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
/// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
/// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
/// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
/// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
/// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
/// IN THE SOFTWARE.

#pragma once

#include <Storm/Base.hpp>

#include <Storm/Mallard/Mesh.hpp>

#include <array>
#include <cstdint>
#include <initializer_list>
#include <ranges>
#include <tuple>

namespace Storm::shapes {

/// @brief Axis-aligned bounding box (AABB).
template<class Vec>
class AABB {
private:

  Vec min_{}, max_{};

public:

  /// @brief Construct an AABB.
  constexpr AABB() = default;

  /// @brief Construct an AABB with a vector @p vec.
  constexpr AABB(const Vec& vec) noexcept : min_{vec}, max_{vec} {}

  /// @brief AABB min vector.
  [[nodiscard]] constexpr const Vec& min() const noexcept {
    return min_;
  }
  /// @brief AABB max vector.
  [[nodiscard]] constexpr const Vec& max() const noexcept {
    return max_;
  }

  /// @brief AABB center.
  [[nodiscard]] constexpr Vec center() const noexcept {
    return (max_ + min_) / 2;
  }
  /// @brief AABB extents.
  [[nodiscard]] constexpr Vec extents() const noexcept {
    return max_ - min_;
  }

  /// @brief Extend the AABB.
  void extend(const Vec& vec) noexcept {
    min_ = glm::min(min_, vec), max_ = glm::max(max_, vec);
  }

}; // class AABB

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
// clang-format off
template<class Shape>
concept shape = 
    requires(Shape shape) {
      { shape.type() } noexcept -> std::same_as<Type>;
      { shape.nodes() } -> std::ranges::range;
    } && 
    std::same_as<std::ranges::range_value_t<
        decltype(std::declval<Shape>().nodes())>, NodeIndex>;
// clang-format on

/// @brief 1D Shape concept.
// clang-format off
template<class Shape>
concept shape1D = shape<Shape> &&
    !requires(Shape shape) { shape.edges(); } &&
    !requires(Shape shape) { shape.faces(); };
// clang-format on

/// @brief 2D Shape concept.
// clang-format off
template<class Shape>
concept shape2D = shape<Shape> &&
    requires(Shape shape) { 
      std::apply([](shape1D auto&&...) {}, shape.edges());
    } && !requires(Shape shape) { shape.faces(); };
// clang-format on

/// @brief 3D Shape concept.
// clang-format off
template<class Shape>
concept shape3D = shape<Shape> &&
    requires(Shape shape) { 
      std::apply([](shape1D auto&&...) {}, shape.edges());
    } && requires(Shape shape) { 
      std::apply([](shape2D auto&&...) {}, shape.faces());
    };
// clang-format on

/// @brief Shape topological dimensionality.
// clang-format off
template<shape Shape>
  requires (shape1D<Shape> || shape2D<Shape> || shape3D<Shape>)
inline constexpr size_t shape_dim_v = []() {
  // clang-format on
  if constexpr (shape1D<Shape>) return 1;
  if constexpr (shape2D<Shape>) return 2;
  if constexpr (shape3D<Shape>) return 3;
}();

/// @brief Get the shape parts.
// clang-format off
template<size_t Index, shape Shape>
  requires ((Index == 0) ||
            (Index == 1 && (shape2D<Shape> || shape3D<Shape>)) ||
            (Index == 2 && shape3D<Shape>))
[[nodiscard]] constexpr auto parts(const Shape& shape) noexcept {
  // clang-format on
  if constexpr (Index == 0) {
    return shape.nodes();
  } else if constexpr (Index == 1) {
    return shape.edges();
  } else if constexpr (Index == 2) {
    return shape.faces();
  }
}

/// @brief Complex shape concept.
// clang-format off
template<class Shape, class Mesh>
concept complex_shape = shape<Shape> && mesh<Mesh> && 
    requires(Shape shape, const Mesh& mesh) { 
      { shape.pieces(mesh) } -> std::ranges::range; 
    } && 
    shape<std::ranges::range_value_t<decltype(
        std::declval<Shape>().pieces(std::declval<Mesh>()))>>;
// clang-format on

/// @brief Piece type of a complex.
// clang-format off
template<class Shape, class Mesh>
  requires complex_shape<Shape, Mesh>
using piece_t =
    std::ranges::range_value_t<decltype(
        std::declval<Shape>().pieces(std::declval<Mesh>()))>;
// clang-format on

namespace detail_ {
  // clang-format off
  template<class Shape, class Mesh>
  concept has_volume_ =
      requires(Shape shape, const Mesh& mesh) { shape.volume(mesh); };
  template<class Shape, class Mesh>
  concept has_barycenter_ =
      requires(Shape shape, const Mesh& mesh) { shape.barycenter(mesh); };
  template<class Shape, class Mesh>
  concept has_normal_ =
      requires(Shape shape, const Mesh& mesh) { shape.normal(mesh); };
  // clang-format on

  // clang-format off
  template<class Shape, class Mesh>
  concept can_volume_ =
      requires(Shape shape, const Mesh& mesh) { volume(shape, mesh); };
  template<class Shape, class Mesh>
  concept can_barycenter_ =
      requires(Shape shape, const Mesh& mesh) { barycenter(shape, mesh); };
  template<class Shape, class Mesh>
  concept can_normal_ =
      requires(Shape shape, const Mesh& mesh) { normal(shape, mesh); };
  // clang-format on
} // namespace detail_

/// @brief Compute the shape "volume" (length in 1D, area in 2D, volume in 3D).
// clang-format off
template<shape Shape, mesh Mesh>
  requires (detail_::has_volume_<Shape, Mesh> || 
            (complex_shape<Shape, Mesh> && 
             detail_::can_volume_<piece_t<Shape, Mesh>, Mesh>))
[[nodiscard]] constexpr auto volume(const Shape& shape,
                                    const Mesh& mesh) {
  // clang-format on
  if constexpr (detail_::has_volume_<Shape, Mesh>) {
    return shape.volume(mesh);
  } else if constexpr (complex_shape<Shape, Mesh> &&
                       detail_::can_volume_<piece_t<Shape, Mesh>, Mesh>) {
    // Compute the complex volume.
    const auto pieces = shape.pieces(mesh);
    auto vol = volume(pieces.front(), mesh);
    for (const auto& piece : pieces | std::views::drop(1)) {
      vol += volume(piece, mesh);
    }
    return vol;
  }
}

/// @brief Compute the shape barycenter.
// clang-format off
template<shape Shape, mesh Mesh>
  requires (detail_::has_barycenter_<Shape, Mesh> || 
            (complex_shape<Shape, Mesh> &&
             detail_::can_volume_<piece_t<Shape, Mesh>, Mesh> &&
             detail_::can_barycenter_<piece_t<Shape, Mesh>, Mesh>))
[[nodiscard]] constexpr auto barycenter(const Shape& shape,
                                        const Mesh& mesh) {
  // clang-format on
  if constexpr (detail_::has_barycenter_<Shape, Mesh>) {
    return shape.barycenter(mesh);
  } else if constexpr (complex_shape<Shape, Mesh> &&
                       detail_::can_volume_<piece_t<Shape, Mesh>, Mesh> &&
                       detail_::can_barycenter_<piece_t<Shape, Mesh>, Mesh>) {
    // Compute the complex barycenter.
    const auto pieces = shape.pieces(mesh);
    auto vol = volume(pieces.front(), mesh);
    auto vol_center = vol * barycenter(pieces.front(), mesh);
    for (const auto& piece : pieces | std::views::drop(1)) {
      const auto dv = volume(piece, mesh);
      vol += dv, vol_center += dv * barycenter(piece, mesh);
    }
    return vol_center / vol;
  }
}

/// @brief Compute the shape normal.
// clang-format off
template<shape Shape, mesh Mesh>
  requires (detail_::has_normal_<Shape, Mesh> || 
            (complex_shape<Shape, Mesh> &&
             detail_::can_volume_<piece_t<Shape, Mesh>, Mesh> &&
             detail_::can_normal_<piece_t<Shape, Mesh>, Mesh>))
[[nodiscard]] constexpr auto normal(const Shape& shape,
                                    const Mesh& mesh) {
  // clang-format on
  if constexpr (detail_::has_normal_<Shape, Mesh>) {
    return shape.normal(mesh);
  } else if constexpr (complex_shape<Shape, Mesh> &&
                       detail_::can_normal_<piece_t<Shape, Mesh>, Mesh>) {
    // Compute the complex barycenter.
    const auto pieces = shape.pieces(mesh);
    auto vol_normal =
        volume(pieces.front(), mesh) * normal(pieces.front(), mesh);
    for (const auto& piece : pieces | std::views::drop(1)) {
      vol_normal += volume(piece, mesh) * normal(piece, mesh);
    }
    return glm::normalize(vol_normal);
  }
}

/// @brief Segmental shape.
/// @verbatim
///
///  n1 @ f1
///      \
///       \         e1 = (n1,n2)
///        v e1     f1 = (n1)
///         \       f2 = (n2)
///          \
///        n2 @ f2
///
/// @endverbatim
class Seg final {
public:

  /// @brief Segment nodes.
  NodeIndex n1, n2;

  /// @brief Construct a segment.
  /// @{
  constexpr Seg() = default;
  constexpr Seg(NodeIndex i1, NodeIndex i2) noexcept : n1{i1}, n2{i2} {}
  /// @}

  /// @brief Segment type.
  [[nodiscard]] constexpr static auto type() noexcept {
    return Type::segment;
  }

  /// @brief Segment nodes.
  [[nodiscard]] constexpr auto nodes() const noexcept {
    return std::array{n1, n2};
  }

  /// @brief Segment "volume" (length).
  template<mesh Mesh>
  [[nodiscard]] constexpr real_t volume(const Mesh& mesh) const noexcept {
    const auto v1{mesh.position(n1)}, v2{mesh.position(n2)};
    return glm::length(v2 - v1);
  }

  /// @brief Segment barycenter.
  template<mesh Mesh>
  [[nodiscard]] constexpr auto barycenter(const Mesh& mesh) const noexcept {
    const auto v1{mesh.position(n1)}, v2{mesh.position(n2)};
    return (v1 + v2) / 2.0;
  }

  /// @brief Segment normal.
  // clang-format off
  template<mesh Mesh>
    requires (mesh_dim_v<Mesh> == 2)
  [[nodiscard]] constexpr auto normal(const Mesh& mesh) const noexcept {
    // clang-format on
    const auto v1{mesh.position(n1)}, v2{mesh.position(n2)};
    const glm::dvec2 d = glm::normalize(v2 - v1);
    return glm::vec2(-d.y, d.x);
  }

}; // class Seg

/// Triangular shape.
/// @verbatim
///
///           n3
///           @           e1 = f1 = (n1,n2)
///          / \          e2 = f2 = (n2,n3)
///         /   \         e3 = f3 = (n3,n1)
///  e3/f3 v     ^ e2/f2
///       /       \
///      /         \
///  n1 @----->-----@ n2
///         e1/f1
///
/// @endverbatim
class Triangle final {
public:

  /// @brief Triangle nodes.
  NodeIndex n1, n2, n3;

  /// @brief Construct a triangle.
  /// @{
  constexpr Triangle() = default;
  constexpr Triangle(NodeIndex i1, NodeIndex i2, NodeIndex i3) noexcept
      : n1{i1}, n2{i2}, n3{i3} {}
  /// @}

  /// @brief Triangle type.
  [[nodiscard]] constexpr static auto type() noexcept {
    return Type::triangle;
  }

  /// @brief Triangle nodes.
  [[nodiscard]] constexpr auto nodes() const noexcept {
    return std::array{n1, n2, n3};
  }

  /// @brief Triangle edges.
  [[nodiscard]] constexpr auto edges() const noexcept {
    return std::array{Seg{n1, n2}, Seg{n2, n3}, Seg{n3, n1}};
  }

  /// @brief Triangle "volume" (area).
  // clang-format off
  template<mesh Mesh>
    requires (mesh_dim_v<Mesh> >= 2)
  [[nodiscard]] constexpr real_t volume(const Mesh& mesh) const noexcept {
    // clang-format on
    const auto v1{mesh.position(n1)}, v2{mesh.position(n2)};
    const auto v3{mesh.position(n3)};
    if constexpr (mesh_dim_v<Mesh> == 2) {
      const glm::dmat2 d{v2 - v1, v3 - v1};
      return glm::abs(glm::determinant(d)) / 2.0;
    } else {
      return glm::length(glm::cross(v2 - v1, v3 - v1)) / 2.0;
    }
  }

  /// @brief Triangle barycenter.
  // clang-format off
  template<mesh Mesh>
    requires (mesh_dim_v<Mesh> >= 2)
  [[nodiscard]] constexpr auto barycenter(const Mesh& mesh) const noexcept {
    // clang-format on
    const auto v1{mesh.position(n1)}, v2{mesh.position(n2)};
    const auto v3{mesh.position(n3)};
    return (v1 + v2 + v3) / 3.0;
  }

  /// @brief Triangle normal.
  // clang-format off
  template<mesh Mesh>
    requires (mesh_dim_v<Mesh> == 3)
  [[nodiscard]] constexpr auto normal(const Mesh& mesh) const noexcept {
    // clang-format on
    const auto v1{mesh.position(n1)}, v2{mesh.position(n2)};
    const auto v3{mesh.position(n3)};
    return glm::normalize(glm::cross(v2 - v1, v3 - v1));
  }

}; // class Triangle

/// @brief Quadrangular shape.
/// @verbatim
///
///               e4/f4
///       n1 @-----<-----@ n4    e1 = f1 = (n1,n2)
///         /           /        e2 = f2 = (n2,n3)
///  e1/f1 v           ^ e3/f3   e3 = f3 = (n3,n4)
///       /           /          e4 = f4 = (n4,n1)
///   n2 @----->-----@ n3     split = ((n1,n2,n3),(n3,n4,n1))
///          e2/f2
///
/// @endverbatim
class Quadrangle final {
public:

  /// @brief Quadrangle nodes.
  NodeIndex n1, n2, n3, n4;

  /// @brief Construct a quadrangle.
  /// @{
  constexpr Quadrangle() = default;
  constexpr Quadrangle(NodeIndex i1, NodeIndex i2, //
                       NodeIndex i3, NodeIndex i4) noexcept
      : n1{i1}, n2{i2}, n3{i3}, n4{i4} {}
  /// @}

  /// @brief Quadrangle type.
  [[nodiscard]] constexpr static auto type() noexcept {
    return Type::quadrangle;
  }

  /// @brief Quadrangle nodes.
  [[nodiscard]] constexpr auto nodes() const noexcept {
    return std::array{n1, n2, n3, n4};
  }

  /// @brief Quadrangle edges.
  [[nodiscard]] constexpr auto edges() const noexcept {
    return std::array{Seg{n1, n2}, Seg{n2, n3}, Seg{n3, n4}, Seg{n4, n1}};
  }

  /// @brief Quadrangle pieces.
  // clang-format off
  template<mesh Mesh>
    requires (mesh_dim_v<Mesh> >= 2)
  [[nodiscard]] constexpr auto pieces(const Mesh& mesh) const noexcept {
    // clang-format on
    static_cast<void>(mesh); /// @todo Convexity check!
    return std::array{Triangle{n1, n2, n3}, Triangle{n3, n4, n1}};
    // return std::array{Triangle{n1, n2, n4}, Triangle{n2, n3, n4}};
  }

}; // class Quadrangle

/// @brief Triangle strip shape.
/// @verbatim
///
///           n2      n4      n6      n8
///           @---<---@---<---@---<---@
///          / \     / \     / \     /
///         v   \   /   \   /   \   ^
///        /     \ /     \ /     \ /
///       @--->---@--->---@--->---@
///       n1      n3      n5      n9
///
/// @endverbatim
class TriangleStrip final {
public:

  /// @brief Triangle strip nodes.
  std::vector<NodeIndex> n{};

  /// @brief Construct a triangle strip.
  /// @{
  constexpr TriangleStrip() = default;
  constexpr TriangleStrip(std::initializer_list<NodeIndex> i) : n{i} {}
  /// @}

  /// @brief Triangle strip type.
  /// May fallback to triangle or quadrangle types.
  [[nodiscard]] constexpr auto type() const noexcept {
    STORM_ASSERT_(n.size() >= 3, "Triangle strip must have least 3 nodes!");
    switch (n.size()) {
      case 3: return Type::triangle;
      case 4: return Type::quadrangle;
      default: return Type::triangle_strip;
    }
  }

  /// @brief Triangle strip nodes.
  [[nodiscard]] constexpr const auto& nodes() const noexcept {
    return n;
  }

  /// @brief Triangle strip edges.
  template<class T = void>
  [[nodiscard]] constexpr auto edges() const noexcept {
    static_assert(meta::always_false<T>,
                  "TriangleStrip::edges is not implemented yet!");
    return std::vector<Seg>{};
  }

  /// @brief Triangle strip pieces.
  template<mesh Mesh>
  [[nodiscard]] constexpr auto pieces(const Mesh& mesh) const noexcept {
    static_cast<void>(mesh); /// @todo Convexity check!
    return std::views::iota(size_t{2}, n.size()) |
           std::views::transform([&](size_t i) {
             return Triangle{n[i - 2], n[i - 1], n[i]};
           });
  }

}; // class TriangleStrip

/// @brief Arbitrary polygon shape.
/// @verbatim
///
///            n5
///           .@---<--@ n4
///         ./`      /
///        v        ^
///     ./`        /
///  n6 @         @ n3
///     |          \
///     v           ^ e2/f2
///     |            \
///     @------>------@
///     n1     |     n2
///          e1/f1
///
/// @endverbatim
class Polygon final {
public:

  /// @brief Polygon nodes.
  std::vector<NodeIndex> n{};

  /// @brief Construct a polygon.
  /// @{
  constexpr Polygon() = default;
  constexpr Polygon(std::initializer_list<NodeIndex> i) : n{i} {}
  /// @}

  /// @brief Polygon type.
  /// May fallback to triangle or quadrangle types.
  [[nodiscard]] constexpr auto type() const noexcept {
    STORM_ASSERT_(n.size() >= 3, "Polygon must have least 3 nodes!");
    switch (n.size()) {
      case 3: return Type::triangle;
      case 4: return Type::quadrangle;
      default: return Type::polygon;
    }
  }

  /// @brief Polygon nodes.
  [[nodiscard]] constexpr const auto& nodes() const noexcept {
    return n;
  }

  /// @brief Polygon edges.
  [[nodiscard]] constexpr auto edges() const noexcept {
    return std::views::iota(size_t{1}, n.size()) |
           std::views::transform([&](size_t i) {
             return Seg{n[i - 1], n[i]};
           });
  }

  /// @brief Polygon pieces.
  template<mesh Mesh>
  [[nodiscard]] constexpr auto pieces(const Mesh& mesh) const {
    static_assert(meta::always_false<Mesh>,
                  "Polygon::pieces is not implemented yet!");
    static_cast<void>(mesh); /// @todo Convexity check!
    return std::vector<Triangle>{};
  }

}; // class Polygon

/// @brief Tetrahedral shape.
/// @verbatim
///                 f4
///            n4   ^
///             @   |
///      f2    /|\. |     f3         e1 = (n1,n2)
///      ^    / | `\.     ^          e2 = (n2,n3)
///       \  /  |   `\.  /           e3 = (n3,n1)
///        \`   |   | `\/            e4 = (n1,n4)
///        ,\   |   o  /`\           e5 = (n2,n4)
///    e4 ^  *  |     *   `^.e6      e6 = (n3,n4)
///      /   e5 ^           `\       f1 = (n1,n3,n2)
///  n1 @-------|--<----------@ n3   f2 = (n1,n2,n4)
///      \      |  e3       ,/       f3 = (n2,n3,n4)
///       \     |     o   ,/`        f4 = (n3,n1,n4)
///        \    ^ e5  | ,/`
///      e1 v   |     ,^ e2
///          \  |   ,/`
///           \ | ,/` |
///            \|/`   |
///             @ n2  v
///                   f1
/// @endverbatim
class Tetrahedron final {
public:

  /// @brief Tetrahedron nodes.
  NodeIndex n1, n2, n3, n4;

  /// @brief Construct a tetrahedron.
  /// @{
  constexpr Tetrahedron() = default;
  constexpr Tetrahedron(NodeIndex i1, NodeIndex i2, //
                        NodeIndex i3, NodeIndex i4) noexcept
      : n1{i1}, n2{i2}, n3{i3}, n4{i4} {}
  /// @}

  /// @brief Tetrahedron type.
  [[nodiscard]] constexpr static auto type() noexcept {
    return Type::tetrahedron;
  }

  /// @brief Tetrahedron nodes.
  [[nodiscard]] constexpr auto nodes() const noexcept {
    return std::array{n1, n2, n3, n4};
  }

  /// @brief Tetrahedron edges.
  [[nodiscard]] constexpr auto edges() const noexcept {
    return std::array{Seg{n1, n2}, Seg{n2, n3}, Seg{n3, n1},
                      Seg{n1, n4}, Seg{n2, n4}, Seg{n3, n4}};
  }

  /// @brief Tetrahedron faces.
  [[nodiscard]] constexpr auto faces() const noexcept {
    return std::tuple{Triangle{n1, n3, n2}, Triangle{n1, n2, n4},
                      Triangle{n2, n3, n4}, Triangle{n3, n1, n4}};
  }

  /// @brief Tetrahedron volume.
  // clang-format off
  template<mesh Mesh>
    requires (mesh_dim_v<Mesh> >= 3)
  [[nodiscard]] constexpr real_t volume(const Mesh& mesh) const noexcept {
    // clang-format on
    const auto v1{mesh.position(n1)}, v2{mesh.position(n2)};
    const auto v3{mesh.position(n3)}, v4{mesh.position(n4)};
    return glm::abs(glm::dot(v2 - v1, glm::cross(v3 - v1, v4 - v1))) / 6.0;
  }

  /// @brief Tetrahedron barycenter.
  // clang-format off
  template<mesh Mesh>
    requires (mesh_dim_v<Mesh> >= 3)
  [[nodiscard]] constexpr auto barycenter(const Mesh& mesh) const noexcept {
    // clang-format on
    const auto v1{mesh.position(n1)}, v2{mesh.position(n2)};
    const auto v3{mesh.position(n3)}, v4{mesh.position(n4)};
    return (v1 + v2 + v3 + v4) / 4.0;
  }

}; // class Tetrahedron

/// @brief Pyramidal shape.
/// @verbatim
///                                n5                      e1 = (n1,n2)
///                  f4           ,@                       e2 = (n2,n3)
///                   ^        ,/`/|\     f2               e3 = (n3,n4)
///                    \    ,/`  / | \    ^                e4 = (n4,n1)
///                     \,/`    /  |  \  /                 e5 = (n1,n5)
///                e8 ,/`\     /   |   \/                  e6 = (n2,n5)
///                ,^`    o   /    |   /\                  e7 = (n3,n5)
///             ,/`          /     |  *  \                 e8 = (n4,n5)
///  f5 <------------*      /   e7 ^   o--\---------> f3   f1 = (n1,n4,n3,n2)
///       ,/`              /       |       \               f2 = (n1,n2,n5)
///   n4 @-----<----------/--------@  n3    ^ e6           f3 = (n2,n3,n5)
///       `\.  e3        /          `\.      \             f4 = (n3,n4,n5)
///          `>.        ^ e5           `\. e2 \            f5 = (n4,n1,n5)
///          e4 `\.    /       o          `<.  \       split = ((n1,n2,n3,n5),
///                `\./        |             `\.\               (n3,n4,n1,n5))
///               n1 @-------------------->-----@ n2
///                            |          e1
///                            |
///                            v
///                            f1
/// @endverbatim
class Pyramid final {
public:

  /// @brief Pyramid nodes.
  NodeIndex n1, n2, n3, n4, n5;

  /// @brief Construct a tetrahedron.
  /// @{
  constexpr Pyramid() = default;
  constexpr Pyramid(NodeIndex i1, NodeIndex i2, //
                    NodeIndex i3, NodeIndex i4, //
                    NodeIndex i5) noexcept
      : n1{i1}, n2{i2}, n3{i3}, n4{i4}, n5{i5} {}
  /// @}

  /// @brief Pyramid type.
  [[nodiscard]] constexpr static auto type() noexcept {
    return Type::pyramid;
  }

  /// @brief Pyramid nodes.
  [[nodiscard]] constexpr auto nodes() const noexcept {
    return std::array{n1, n2, n3, n4, n5};
  }

  /// @brief Pyramid edges.
  [[nodiscard]] constexpr auto edges() const noexcept {
    return std::array{Seg{n1, n2}, Seg{n2, n3}, Seg{n3, n4}, Seg{n4, n1},
                      Seg{n1, n5}, Seg{n2, n5}, Seg{n3, n5}, Seg{n4, n5}};
  }

  /// @brief Pyramid faces.
  [[nodiscard]] constexpr auto faces() const noexcept {
    return std::tuple{Quadrangle{n1, n4, n3, n2}, Triangle{n1, n2, n5},
                      Triangle{n2, n3, n5}, Triangle{n3, n4, n5},
                      Triangle{n4, n1, n5}};
  }

  /// @brief Pyramid pieces.
  // clang-format off
  template<mesh Mesh>
    requires (mesh_dim_v<Mesh> >= 3)
  [[nodiscard]] constexpr auto pieces(const Mesh& mesh) const noexcept {
    // clang-format on
    static_cast<void>(mesh); /// @todo Convexity check!
    return std::array{Tetrahedron{n1, n2, n3, n5}, Tetrahedron{n3, n4, n1, n5}};
  }

}; // class Pyramid

/// @brief Pentahedral shape (triangular prism).
/// @verbatim
///                 f5
///                 ^  f3
///                 |  ^
///             e9  |  |
///      n4 @---<---|-------------@ n6        e1 = (n1,n2)
///         |\      *  |        ,/|           e2 = (n2,n3)
///         | \        o      ,/` |           e3 = (n3,n1)
///         |  \         e8 ,^`   |           e4 = (n1,n4)
///         |   v e7      ,/`     |           e5 = (n2,n5)
///      e4 ^    \      ,/`       ^ e6        e6 = (n3,n6)
///         |     \   ,/`         |           e7 = (n4,n5)
///         |      \ /`        *-------> f2   e8 = (n5,n6)
///  f1 <-------*   @ n5          |           e9 = (n6,n4)
///         |       |             |           f1 = (n1,n2,n5,n4)
///      n1 @-------|---------<---@ n3        f2 = (n2,n3,n6,n5)
///          \      |        e3 ,/            f3 = (n3,n1,n4,n6)
///           \     |     o   ,/`             f4 = (n1,n3,n2)
///            \    ^ e5  | ,/`               f5 = (n4,n5,n6)
///          e1 v   |     ,^ e2            split = ((n1,n2,n3,n5),
///              \  |   ,/|                         (n3,n1,n4,n5),
///               \ | ,/` |                         (n4,n6,n3,n5))
///                \|/`   |
///                 @ n2  v
///                       f4
/// @endverbatim
class Pentahedron final {
public:

  /// @brief Pentahedron nodes.
  NodeIndex n1, n2, n3, n4, n5, n6;

  /// @brief Construct a pentahedron.
  /// @{
  constexpr Pentahedron() = default;
  constexpr Pentahedron(NodeIndex i1, NodeIndex i2, //
                        NodeIndex i3, NodeIndex i4, //
                        NodeIndex i5, NodeIndex i6) noexcept
      : n1{i1}, n2{i2}, n3{i3}, n4{i4}, n5{i5}, n6{i6} {}
  /// @}

  /// @brief Pentahedron type.
  [[nodiscard]] constexpr static auto type() noexcept {
    return Type::pentahedron;
  }

  /// @brief Pentahedron nodes.
  [[nodiscard]] constexpr auto nodes() const noexcept {
    return std::array{n1, n2, n3, n4, n5, n6};
  }

  /// @brief Pentahedron edges.
  [[nodiscard]] constexpr auto edges() const noexcept {
    return std::array{Seg{n1, n2}, Seg{n2, n3}, Seg{n3, n1},
                      Seg{n1, n4}, Seg{n2, n5}, Seg{n3, n6},
                      Seg{n4, n5}, Seg{n5, n6}, Seg{n6, n4}};
  }

  /// @brief Pentahedron faces.
  [[nodiscard]] constexpr auto faces() const noexcept {
    return std::tuple{Quadrangle{n1, n2, n5, n4}, Quadrangle{n2, n3, n6, n5},
                      Quadrangle{n3, n1, n4, n6}, Triangle{n1, n3, n2},
                      Triangle{n4, n5, n6}};
  }

  /// @brief Pentahedron pieces.
  // clang-format off
  template<mesh Mesh>
    requires (mesh_dim_v<Mesh> >= 3)
  [[nodiscard]] constexpr auto pieces(const Mesh& mesh) const noexcept {
    // clang-format on
    static_cast<void>(mesh); /// @todo Convexity check!
    return std::array{Tetrahedron{n1, n2, n3, n5}, Tetrahedron{n3, n1, n4, n5},
                      Tetrahedron{n4, n6, n3, n5}};
  }

}; // class Pentahedron

/// @brief Hexahedral shape.
/// @verbatim
///                      f6
///                      ^       f3
///                      |       ^
///                  e10 |      /
///            n7 @---<--|----------@ n6         e1 = (n1,n2)
///              /|      |    /    /|            e2 = (n2,n3)
///             / |      |   o    / |            e3 = (n3,n4)
///        e11 v  |      *    e9 ^  ^ e6         e4 = (n4,n1)
///           /   ^ e7          /   |            e5 = (n1,n5)
///          /    |      e12   /  *-------> f2   e6 = (n2,n6)
///      n8 @------------->---@ n5  |            e7 = (n3,n7)
///  f4 <---|--o  |           |     |            e8 = (n4,n8)
///         |  n3 @---<-------|-----@ n2         e9 = (n5,n6)
///         |    /    e2      |    /            e10 = (n6,n7)
///      e8 ^   /          e5 ^   /             e11 = (n7,n8)
///         |  v e3  *        |  ^ e1           e12 = (n8,n5)
///         | /     /    o    | /                f1 = (n1,n4,n3,n2)
///         |/     /     |    |/                 f2 = (n1,n2,n6,n5)
///      n4 @-----/-->--------@ n1               f3 = (n2,n3,n7,n6)
///              /   e4  |                       f4 = (n3,n4,n8,n7)
///             /        v                       f5 = (n1,n5,n8,n4)
///            v         f1                      f6 = (n5,n6,n7,n8)
///            f5
/// @endverbatim
class Hexahedron final {
public:

  /// @brief Hexahedron nodes.
  NodeIndex n1, n2, n3, n4, n5, n6, n7, n8;

  /// @brief Construct a hexahedron.
  /// @{
  constexpr Hexahedron() = default;
  constexpr Hexahedron(NodeIndex i1, NodeIndex i2, //
                       NodeIndex i3, NodeIndex i4, //
                       NodeIndex i5, NodeIndex i6, //
                       NodeIndex i7, NodeIndex i8) noexcept
      : n1{i1}, n2{i2}, n3{i3}, n4{i4}, n5{i5}, n6{i6}, n7{i7}, n8{i8} {}
  /// @}

  /// @brief Hexahedron type.
  [[nodiscard]] constexpr static auto type() noexcept {
    return Type::hexahedron;
  }

  /// @brief Hexahedron nodes.
  [[nodiscard]] constexpr auto nodes() const noexcept {
    return std::array{n1, n2, n3, n4, n5, n6, n7, n8};
  }

  /// @brief Hexahedron edges.
  [[nodiscard]] constexpr auto edges() const noexcept {
    return std::array{Seg{n1, n2}, Seg{n2, n3}, Seg{n3, n4}, Seg{n4, n1},
                      Seg{n1, n5}, Seg{n2, n6}, Seg{n3, n7}, Seg{n4, n8},
                      Seg{n5, n6}, Seg{n6, n7}, Seg{n7, n8}, Seg{n8, n5}};
  }

  /// @brief Hexahedron faces.
  [[nodiscard]] constexpr auto faces() const noexcept {
    return std::array{Quadrangle{n1, n4, n3, n2}, Quadrangle{n1, n2, n6, n5},
                      Quadrangle{n2, n3, n7, n6}, Quadrangle{n3, n4, n8, n7},
                      Quadrangle{n1, n5, n8, n4}, Quadrangle{n5, n6, n7, n8}};
  }

  /// @brief Hexahedron pieces.
  // clang-format off
  template<mesh Mesh>
    requires (mesh_dim_v<Mesh> >= 3)
  [[nodiscard]] constexpr auto pieces(const Mesh& mesh) const noexcept {
    // clang-format on
    static_cast<void>(mesh); /// @todo Convexity check!
    return std::array{Tetrahedron{n1, n4, n2, n5}, Tetrahedron{n4, n3, n2, n7},
                      Tetrahedron{n5, n6, n7, n2}, Tetrahedron{n5, n7, n8, n4},
                      Tetrahedron{n5, n4, n2, n7}};
  }

}; // class Hexahedron

} // namespace Storm::shapes
