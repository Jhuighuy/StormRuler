/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// Copyright (C) 2022 Oleg Butakov
///
/// Permission is hereby granted, free of charge, to any person
/// obtaining a copy of this software and associated documentation
/// files (the "Software"), to deal in the Software without
/// restriction, including without limitation the rights  to use,
/// copy, modify, merge, publish, distribute, sublicense, and/or
/// sell copies of the Software, and to permit persons to whom the
/// Software is furnished to do so, subject to the following
/// conditions:
///
/// The above copyright notice and this permission notice shall be
/// included in all copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
/// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
/// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
/// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
/// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
/// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
/// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
/// OTHER DEALINGS IN THE SOFTWARE.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///

#pragma once

#include <span>
#include <vector>
#include <memory>
#include <stdexcept>

#include <stormBase.hxx>
#include <stormMesh/Forward.hxx>

namespace Storm {

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract mesh element.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<size_t Dim>
class Element {
protected:
  std::span<GVec<Dim> const> NodeCoords_;
  std::vector<size_t> NodeIndices_;

  template<typename... SizeTypes>
  ElementDescription Part_(ShapeType shape, SizeTypes... indices) const {
    return { shape, std::vector<size_t>{ NodeIndices_[indices]... } };
  }

public:

  /// @brief Make a new element with description.
  static std::unique_ptr<Element> Make(ElementDescription&& desc,
                                       std::span<GVec<Dim> const> nodeCoords);

  virtual ~Element() = default;

  /// @brief Element node indices.
  std::vector<size_t> const& NodeIndices() const {
    return NodeIndices_;
  }

  /// @brief Node position.
  GVec<Dim> const& NodeCoords(size_t index) const {
    StormAssert(index < NodeIndices_.size());
    return NodeCoords_[NodeIndices_[index]];
  }

  /// @brief Element diameter.
  virtual real_t Diameter() const {
    throw std::runtime_error("Not overridden.");
  }

  /// @brief Element length/area/volume.
  virtual real_t Volume() const {
    throw std::runtime_error("Not overridden.");
  }

  /// @brief Normal to element.
  virtual GVec<Dim> Normal() const {
    throw std::runtime_error("Not overridden.");
  }

  /// @brief Element direction.
  virtual GVec<Dim> Direction() const {
    throw std::runtime_error("Not overridden.");
  }

  /// @brief Element barycenter.
  virtual GVec<Dim> CenterCoords() const {
    throw std::runtime_error("Not overridden.");
  }

  /// @brief Element shape.
  virtual ShapeType Shape() const = 0;

  /// @brief Number of nodes in the element.
  virtual size_t NumNodes() const = 0;

  /// @brief Number of edges in the element.
  size_t NumEdges() const {
    return EdgesDesc().size();
  }

  /// @brief Element edges description.
  virtual ElementDescriptionArray EdgesDesc() const = 0;

  /// @brief  Number of faces in the element.
  size_t NumFaces() const {
    return FacesDesc().size();
  }

  /// @brief Element faces description.
  virtual ElementDescriptionArray FacesDesc() const = 0;

}; // class Element<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract simplex element.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<size_t Dim>
class SimplexElement : public Element<Dim> {
}; // class SimplexElement<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract complex (not simplex) element. 
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<size_t Dim>
class ComplexElement : public Element<Dim> {
public:
  real_t Diameter() const final;
  real_t Volume() const final;
  GVec<Dim> Normal() const final;
  GVec<Dim> CenterCoords() const final;

  /// @brief Get splitting into the simplex parts.
  virtual ElementDescriptionArray SimplicialParts() const = 0;

private:

  template<typename Func>
  void ForEachSimplex_(Func func) const;

}; // class ComplexElement<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Dummy nodal element.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<size_t Dim>
class Node final : public SimplexElement<Dim> {
public:

    real_t Diameter() const override;

    real_t Volume() const override;

    GVec<Dim> Normal() const override;

    GVec<Dim> CenterCoords() const override;

    ShapeType Shape() const override;

    size_t NumNodes() const override;

    ElementDescriptionArray EdgesDesc() const override;

    ElementDescriptionArray FacesDesc() const override;

}; // class Node<Dim>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief 1D/2D/3D Segmental element.
/// @verbatim
///
///  n0 O f0
///      \
///       \         e0 = (n0,n1)
///        v e0     f0 = (n0)
///         \       f1 = (n1)
///          \
///        n1 O f1
///
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<size_t Dim>
  requires (1 <= Dim && Dim <= 3)
class Segment final : public SimplexElement<Dim> {
public:

  real_t Diameter() const override;

  real_t Volume() const override;

  GVec<Dim> Normal() const override;

  GVec<Dim> Direction() const override;

  GVec<Dim> CenterCoords() const override;

  ShapeType Shape() const override { 
    return ShapeType::Segment; 
  }

  size_t NumNodes() const override { 
    return 2; 
  }

  ElementDescriptionArray EdgesDesc() const override;

  ElementDescriptionArray FacesDesc() const override;

}; // class Segment<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief 2D/3D Triangular element.
/// @verbatim
///
///           n2
///           O           e0 = f0 = (n0,n1)
///          / \          e1 = f1 = (n1,n2)
///         /   \         e2 = f2 = (n2,n0)
///  e2/f2 v     ^ e1/f1
///       /       \
///      /         \
///  n0 O----->-----O n1
///        e0/f0
///
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<size_t Dim>
  requires (2 <= Dim && Dim <= 3)
class Triangle final : public SimplexElement<Dim> {
public:

  real_t Diameter() const override;

  real_t Volume() const override;

  GVec<Dim> Normal() const override;

  GVec<Dim> CenterCoords() const override;

  ShapeType Shape() const override {
    return ShapeType::Triangle;
  }

  size_t NumNodes() const override {
    return 3;
  }

  ElementDescriptionArray EdgesDesc() const override;

  ElementDescriptionArray FacesDesc() const override;

}; // class Triangle<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief 2D/3D Quadrangular element.
/// @verbatim
///
///               e2/f2
///       n3 O-----<-----O n2    e0 = f0 = (n0,n1)
///         /           /        e1 = f2 = (n1,n2)
///  e3/f3 v           ^ e1/f1   e2 = f2 = (n2,n3)
///       /           /          e3 = f3 = (n3,n0)
///   n0 O----->-----O n1     split = ((n0,n1,n2),(n2,n3,n0))
///          e0/f0
///
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<size_t Dim>
  requires (2 <= Dim && Dim <= 3)
class Quadrangle final : public ComplexElement<Dim> {
public:

  ShapeType Shape() const override {
    return ShapeType::Quadrangle;
  }

  size_t NumNodes() const override {
    return 4;
  }

  ElementDescriptionArray EdgesDesc() const override;

  ElementDescriptionArray FacesDesc() const override;

  ElementDescriptionArray SimplicialParts() const override;

}; // class Quadrangle<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief 3D Tetrahedral element.
/// @verbatim
///                    f3
///               n3   ^
///                O   |
///         f1    /|\. |     f2         e0 = (n0,n1)
///         ^    / | `\.     ^          e1 = (n1,n2)
///          \  /  |   `\.  /           e2 = (n2,n0)
///           \`   |   | `\/            e3 = (n0,n3)
///           ,\   |   o  /`\           e4 = (n1,n3)
///       e3 ^  *  |     *   `^.e5      e5 = (n2,n3)
///         /   e4 ^           `\       f0 = (n0,n2,n1)
///     n0 O-------|--<----------O n2   f1 = (n0,n1,n3)
///         \      |  e2       ,/       f2 = (n1,n2,n3)
///          \     |     o   ,/`        f3 = (n2,n0,n3)
///           \    ^ e4  | ,/`
///         e0 v   |     ,^ e1
///             \  |   ,/`
///              \ | ,/` |
///               \|/`   |
///                O n1  v
///                      f0
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<size_t Dim>
  requires (Dim == 3)
class Tetrahedron final : public SimplexElement<Dim> {
public:

  real_t Diameter() const override;

  real_t Volume() const override;

  GVec<Dim> CenterCoords() const override;

  ShapeType Shape() const override {
    return ShapeType::Tetrahedron;
  }

  size_t NumNodes() const override {
    return 4;
  }

  ElementDescriptionArray EdgesDesc() const override;

  ElementDescriptionArray FacesDesc() const override;

}; // class Tetrahedron<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief 3D Pyramidal element.
/// @verbatim
///                                n4                      e0 = (n0,n1)
///                  f3           ,O                       e1 = (n1,n2)
///                   ^        ,/`/|\     f1               e2 = (n2,n3)
///                    \    ,/`  / | \    ^                e3 = (n3,n0)
///                     \,/`    /  |  \  /                 e4 = (n0,n4)
///                e7 ,/`\     /   |   \/                  e5 = (n1,n4)
///                ,^`    o   /    |   /\                  e6 = (n2,n4)
///             ,/`          /     |  *  \                 e7 = (n3,n4)
///  f4 <------------*      /   e6 ^   o--\---------> f2   f0 = (n0,n3,n2,n1)
///       ,/`              /       |       \               f1 = (n0,n1,n4)
///   n3 O-----<----------/--------O  n2    ^ e5           f2 = (n1,n2,n4)
///       `\.  e2        /          `\.      \             f3 = (n2,n3,n4)
///          `>.        ^ e4           `\. e1 \            f4 = (n3,n0,n4)
///          e3 `\.    /       o          `<.  \        split = ((n0,n1,n2,n4),
///                `\./        |             `\.\                (n2,n3,n0,n4))
///               n0 O-------------------->-----O n1
///                            |          e0
///                            |
///                            v
///                            f0
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<size_t Dim>
  requires (Dim == 3)
class Pyramid final : public ComplexElement<Dim> {
public:

  ShapeType Shape() const override {
    return ShapeType::Pyramid;
  }

  size_t NumNodes() const override {
    return 5;
  }

  ElementDescriptionArray EdgesDesc() const override;

  ElementDescriptionArray FacesDesc() const override;

  ElementDescriptionArray SimplicialParts() const override;

}; // class Pyramid<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief 3D Pentahedral element (triangular prism).
/// @verbatim
///                 f4
///                 ^  f2
///                 |  ^
///             e8  |  |
///      n3 O---<---|-------------O n5        e0 = (n0,n1)
///         |\      *  |        ,/|           e1 = (n1,n2)
///         | \        o      ,/` |           e2 = (n2,n0)
///         |  \         e7 ,^`   |           e3 = (n0,n3)
///         |   v e6      ,/`     |           e4 = (n1,n4)
///      e3 ^    \      ,/`       ^ e5        e5 = (n2,n5)
///         |     \   ,/`         |           e6 = (n3,n4)
///         |      \ /`        *-------> f1   e7 = (n4,n5)
///  f0 <-------*   O n4          |           e8 = (n5,n3)
///         |       |             |           f0 = (n0,n1,n4,n3)
///      n0 O-------|---------<---O n2        f1 = (n1,n2,n5,n4)
///          \      |        e2 ,/            f2 = (n2,n0,n3,n5)
///           \     |     o   ,/`             f3 = (n0,n2,n1)
///            \    ^ e4  | ,/`               f4 = (n3,n4,n5)
///          e0 v   |     ,^ e1
///              \  |   ,/|
///               \ | ,/` |
///                \|/`   |
///                 O n1  v
///                       f3
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<size_t Dim>
  requires (Dim == 3)
class Prism final : public ComplexElement<Dim> {
public:

  ShapeType Shape() const override {
    return ShapeType::Prism;
  }

  size_t NumNodes() const override {
    return 6;
  }

  ElementDescriptionArray EdgesDesc() const override;

  ElementDescriptionArray FacesDesc() const override;

  ElementDescriptionArray SimplicialParts() const override;

}; // class Prism<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief 3D Hexahedral element.
/// @verbatim
///                      f5
///                      ^   f2
///                      |   ^
///                   e9 |   |
///            n6 O---<--|----------O n5         e0 = (n0,n1)
///              /|      |   |     /|            e1 = (n1,n2)
///             / |      |   o    / |            e2 = (n2,n3)
///        e10 v  |      *    e8 ^  ^ e5         e3 = (n3,n0)
///           /   ^ e6          /   |            e4 = (n0,n4)
///          /    |      e11   /  *-------> f1   e5 = (n1,n5)
///      n7 O------------->---O n4  |            e6 = (n2,n6)
///  f3 <---|--o  |           |     |            e7 = (n3,n7)
///         |  n2 O---<-------|-----O n1         e8 = (n4,n5)
///         |    /    e1      |    /             e9 = (n5,n6)
///      e7 ^   /          e4 ^   /             e10 = (n6,n7)
///         |  v e2  *        |  ^ e0           e11 = (n7,n4)
///         | /      |   o    | /                f0 = (n0,n3,n2,n1)
///         |/       |   |    |/                 f1 = (n0,n1,n5,n4)
///      n3 O--->----|--------O n0               f2 = (n1,n2,n6,n5)
///             e3   |   |                       f3 = (n2,n3,n7,n6)
///                  |   v                       f4 = (n0,n4,n7,n3)
///                  v   f0                      f5 = (n4,n5,n6,n7)
///                  f4
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<size_t Dim>
  requires (Dim == 3)
class Hexahedron final : public ComplexElement<Dim> {
public:

  ShapeType Shape() const override {
    return ShapeType::Hexahedron;
  }

  size_t NumNodes() const override {
    return 8;
  }

  ElementDescriptionArray EdgesDesc() const override;

  ElementDescriptionArray FacesDesc() const override;

  ElementDescriptionArray SimplicialParts() const override;

}; // class Hexahedron<...>

} // namespace Storm
