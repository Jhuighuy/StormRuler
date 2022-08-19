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
#include <ranges>
#include <tuple>

namespace Storm {
#if 1 // this should not be here..
namespace detail_ {
  struct LabelTag_;
  template<size_t I>
  struct TopologicalIndexTag_;
} // namespace detail_
using Label = Index<detail_::LabelTag_>;
template<size_t I>
using EntityIndex = Index<detail_::TopologicalIndexTag_<I>>;
using NodeIndex = EntityIndex<0>;
using EdgeIndex = EntityIndex<1>;
#endif
} // namespace Storm

namespace Storm::shapes {

/// @brief Shape concept.
// clang-format off
template<class Shape>
concept shape = 
    requires(Shape& shape) {
      { shape.nodes() } -> std::ranges::range;
    };
// clang-format on

namespace detail_ {
  // clang-format off
  template<class Shape>
  concept has_edges_ =
      requires { std::declval<Shape>().edges(); };
  template<class Shape>
  concept has_faces_ =
      requires { std::declval<Shape>().faces(); };
  // clang-format on
} // namespace detail_

/// @brief Get the shape parts.
// clang-format off
template<size_t TopologicalIndex, shape Shape>
  requires ((TopologicalIndex == 0) ||
            (TopologicalIndex == 1 && detail_::has_edges_<Shape>) ||
            (TopologicalIndex == 2 && detail_::has_faces_<Shape>))
[[nodiscard]] constexpr auto parts(const Shape& shape) noexcept {
  // clang-format on
  if constexpr (TopologicalIndex == 0) {
    return shape.nodes();
  } else if constexpr (TopologicalIndex == 1) {
    return shape.edges();
  } else if constexpr (TopologicalIndex == 2) {
    return shape.faces();
  }
}

/// @brief Complex shape concept.
// clang-format off
template<class Shape>
concept complex_shape =
    shape<Shape> && requires { std::declval<Shape>().simplices(); };
// clang-format on

/// @brief Simplex shape concept.
template<class Shape>
concept simplex_shape = shape<Shape> && !complex_shape<Shape>;

namespace detail_ {
  // clang-format off
    template<class Shape, class Mesh>
    concept has_volume_ =
        requires { std::declval<Shape>().volume(std::declval<Mesh>()); };
    template<class Shape, class Mesh>
    concept has_barycenter_ =
        requires { std::declval<Shape>().barycenter(std::declval<Mesh>()); };
    template<class Shape, class Mesh>
    concept has_normal_ =
        requires { std::declval<Shape>().normal(std::declval<Mesh>()); };
  // clang-format on
} // namespace detail_

real_t volume(auto&&, auto&&) {
  return {};
}

glm::dvec3 barycenter_position(auto&&, auto&&) {
  return {};
}

glm::dvec3 normal(auto&&, auto&&) {
  return {};
}

/// @brief Segmental shape.
/// @verbatim
///
///  n1 O f0
///      \
///       \         e0 = (n1,n2)
///        v e0     f0 = (n1)
///         \       f1 = (n2)
///          \
///        n2 O f1
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

  /// @brief Segment nodes.
  [[nodiscard]] constexpr auto nodes() const noexcept {
    return std::array{n1, n2};
  }

}; // class Seg

/// Triangular shape.
/// @verbatim
///           n3
///           O           e0 = f0 = (n1,n2)
///          / \          e1 = f1 = (n2,n3)
///         /   \         e2 = f2 = (n3,n1)
///  e2/f2 v     ^ e1/f1
///       /       \
///      /         \
///  n1 O----->-----O n2
///         e0/f0
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

  /// @brief Triangle nodes.
  [[nodiscard]] constexpr auto nodes() const noexcept {
    return std::array{n1, n2, n3};
  }

  /// @brief Triangle edges.
  [[nodiscard]] constexpr auto edges() const noexcept {
    return std::array{Seg{n1, n2}, Seg{n2, n3}, Seg{n3, n1}};
  }

}; // class Triangle

/// @brief Quadrangular shape.
/// @verbatim
///               e2/f2
///       n4 O-----<-----O n3    e0 = f0 = (n1,n2)
///         /           /        e1 = f2 = (n2,n3)
///  e3/f3 v           ^ e1/f1   e2 = f2 = (n3,n4)
///       /           /          e3 = f3 = (n4,n1)
///   n1 O----->-----O n2     split = ((n1,n2,n3),(n3,n4,n1))
///          e0/f0
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

  /// @brief Quadrangle nodes.
  [[nodiscard]] constexpr auto nodes() const noexcept {
    return std::array{n1, n2, n3, n4};
  }

  /// @brief Quadrangle edges.
  [[nodiscard]] constexpr auto edges() const noexcept {
    return std::array{Seg{n1, n2}, Seg{n2, n3}, Seg{n3, n4}, Seg{n4, n1}};
  }

  /// @brief Quadrangle simplices.
  [[nodiscard]] constexpr auto simplices() const noexcept {
    return std::array{Triangle{n1, n2, n3}, Triangle{n3, n4, n1}};
  }

}; // class Quadrangle

/// @brief Tetrahedral shape.
/// @verbatim
///                    f3
///               n4   ^
///                O   |
///         f1    /|\. |     f2         e0 = (n1,n2)
///         ^    / | `\.     ^          e1 = (n2,n3)
///          \  /  |   `\.  /           e2 = (n3,n1)
///           \`   |   | `\/            e3 = (n1,n4)
///           ,\   |   o  /`\           e4 = (n2,n4)
///       e3 ^  *  |     *   `^.e5      e5 = (n3,n4)
///         /   e4 ^           `\       f0 = (n1,n3,n2)
///     n1 O-------|--<----------O n3   f1 = (n1,n2,n4)
///         \      |  e2       ,/       f2 = (n2,n3,n4)
///          \     |     o   ,/`        f3 = (n3,n1,n4)
///           \    ^ e4  | ,/`
///         e0 v   |     ,^ e1
///             \  |   ,/`
///              \ | ,/` |
///               \|/`   |
///                O n2  v
///                      f0
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

}; // class Tetrahedron

/// @brief Pyramidal shape.
/// @verbatim
///                                n5                      e0 = (n1,n2)
///                  f3           ,O                       e1 = (n2,n3)
///                   ^        ,/`/|\     f1               e2 = (n3,n4)
///                    \    ,/`  / | \    ^                e3 = (n4,n1)
///                     \,/`    /  |  \  /                 e4 = (n1,n5)
///                e7 ,/`\     /   |   \/                  e5 = (n2,n5)
///                ,^`    o   /    |   /\                  e6 = (n3,n5)
///             ,/`          /     |  *  \                 e7 = (n4,n5)
///  f4 <------------*      /   e6 ^   o--\---------> f2   f0 = (n1,n4,n3,n2)
///       ,/`              /       |       \               f1 = (n1,n2,n5)
///   n4 O-----<----------/--------O  n3    ^ e5           f2 = (n2,n3,n5)
///       `\.  e2        /          `\.      \             f3 = (n3,n4,n5)
///          `>.        ^ e4           `\. e1 \            f4 = (n4,n1,n5)
///          e3 `\.    /       o          `<.  \       split = ((n1,n2,n3,n5),
///                `\./        |             `\.\               (n3,n4,n1,n5))
///               n1 O-------------------->-----O n2
///                            |          e0
///                            |
///                            v
///                            f0
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

  /// @brief Pyramid simplices.
  [[nodiscard]] constexpr auto simplices() const noexcept {
    return std::array{Tetrahedron{n1, n2, n3, n5}, Tetrahedron{n3, n4, n1, n5}};
  }

}; // class Pyramid

/// @brief Pentahedral shape (triangular prism).
/// @verbatim
///                 f4
///                 ^  f2
///                 |  ^
///             e8  |  |
///      n4 O---<---|-------------O n6        e0 = (n1,n2)
///         |\      *  |        ,/|           e1 = (n2,n3)
///         | \        o      ,/` |           e2 = (n3,n1)
///         |  \         e7 ,^`   |           e3 = (n1,n4)
///         |   v e6      ,/`     |           e4 = (n2,n5)
///      e3 ^    \      ,/`       ^ e5        e5 = (n3,n6)
///         |     \   ,/`         |           e6 = (n4,n5)
///         |      \ /`        *-------> f1   e7 = (n5,n6)
///  f0 <-------*   @ n5          |           e8 = (n6,n4)
///         |       |             |           f0 = (n1,n2,n5,n4)
///      n1 O-------|---------<---O n3        f1 = (n2,n3,n6,n5)
///          \      |        e2 ,/            f2 = (n3,n1,n4,n6)
///           \     |     o   ,/`             f3 = (n1,n3,n2)
///            \    ^ e4  | ,/`               f4 = (n4,n5,n6)
///          e0 v   |     ,^ e1            split = ((n1,n2,n3,n5),
///              \  |   ,/|                         (n3,n1,n4,n5),
///               \ | ,/` |                         (n4,n6,n3,n5))
///                \|/`   |
///                 O n2  v
///                       f3
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

  /// @brief Pentahedron simplices.
  [[nodiscard]] constexpr auto simplices() const noexcept {
    return std::array{Tetrahedron{n1, n2, n3, n5}, Tetrahedron{n3, n1, n4, n5},
                      Tetrahedron{n4, n6, n3, n5}};
  }

}; // class Pentahedron

/// @brief Hexahedral shape.
/// @verbatim
///                      f5
///                      ^       f2
///                      |       ^
///                   e9 |      /
///            n7 O---<--|----------O n6         e0 = (n1,n2)
///              /|      |    /    /|            e1 = (n2,n3)
///             / |      |   o    / |            e2 = (n3,n4)
///        e10 v  |      *    e8 ^  ^ e5         e3 = (n4,n1)
///           /   ^ e6          /   |            e4 = (n1,n5)
///          /    |      e11   /  *-------> f1   e5 = (n2,n6)
///      n8 O------------->---O n5  |            e6 = (n3,n7)
///  f3 <---|--o  |           |     |            e7 = (n4,n8)
///         |  n3 O---<-------|-----O n2         e8 = (n5,n6)
///         |    /    e1      |    /             e9 = (n6,n7)
///      e7 ^   /          e4 ^   /             e10 = (n7,n8)
///         |  v e2  *        |  ^ e0           e11 = (n8,n5)
///         | /     /    o    | /                f0 = (n1,n4,n3,n2)
///         |/     /     |    |/                 f1 = (n1,n2,n6,n5)
///      n4 O-----/-->--------O n1               f2 = (n2,n3,n7,n6)
///              /   e3  |                       f3 = (n3,n4,n8,n7)
///             /        v                       f4 = (n1,n5,n8,n4)
///            v         f0                      f5 = (n5,n6,n7,n8)
///            f4
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

  /// @brief Hexahedron simplices.
  [[nodiscard]] constexpr auto simplices() const noexcept {
    return std::array{Tetrahedron{n1, n4, n2, n5}, Tetrahedron{n4, n3, n2, n7},
                      Tetrahedron{n5, n6, n7, n2}, Tetrahedron{n5, n7, n8, n4},
                      Tetrahedron{n5, n4, n2, n7}};
  }

}; // class Hexahedron

} // namespace Storm::shapes
