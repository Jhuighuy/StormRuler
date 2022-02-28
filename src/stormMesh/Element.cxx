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

#include <stormMesh/Element.hxx>

namespace Storm {

template<size_t Dim>
std::unique_ptr<Element<Dim>> MakeElement_(ShapeType shape) {

  if (shape == ShapeType::Node) {
    return std::make_unique<Node<Dim>>();
  }

  if (shape == ShapeType::Segment) {
    return std::make_unique<Segment<Dim>>();
  }

  if constexpr (Dim >= 2) {

    if (shape == ShapeType::Triangle) {
      return std::make_unique<Triangle<Dim>>();
    }

    if (shape == ShapeType::Quadrangle) {
      return std::make_unique<Quadrangle<Dim>>();
    }

  }

  if constexpr (Dim == 3) {

    if (shape == ShapeType::Tetrahedron) {
      return std::make_unique<Triangle<Dim>>();
    }

    if (shape == ShapeType::Pyramid) {
      return std::make_unique<Pyramid<Dim>>();
    }

    if (shape == ShapeType::Prism) {
      return std::make_unique<Prism<Dim>>();
    }

    if (shape == ShapeType::Hexahedron) {
      return std::make_unique<Hexahedron<Dim>>();
    }

  }

  throw std::invalid_argument("Invalid element shape.");

} // MakeElement_<...>

/**
 * Construct a new element object.
 */
template<size_t Dim>
std::unique_ptr<Element<Dim>> Element<Dim>::Make(ElementDescription&& desc,
                                                 std::span<GVec<Dim> const> nodeCoords) {

  auto element = MakeElement_<Dim>(desc.Shape);
  element->NodeCoords_ = nodeCoords;
  element->NodeIndices_ = std::move(desc.NodeIndices);

  return element;

} // Element<...>::Make

template class Element<1>;
template class Element<2>;
template class Element<3>;

template<size_t Dim>
template<typename Func>
void ComplexElement<Dim>::ForEachSimplex_(Func func) const {

  ElementDescriptionArray simplices(SimplicialParts());
  for (ElementDescription& part : simplices) {

    auto const simplex =
      Element<Dim>::Make(std::move(part), this->NodeCoords_);

    func(*simplex);

  }

} // ComplexElement<...>::ForEachSimplex_

template<size_t Dim>
real_t ComplexElement<Dim>::Diameter() const {

  real_t diameter = 0.0;
  ForEachSimplex_(
    [&](Element<Dim> const& simplex) {
      diameter = std::max(diameter, simplex.Diameter());
    });

  return diameter;

} // ComplexElement<...>::get_diameter

template<size_t Dim>
real_t ComplexElement<Dim>::Volume() const {

  real_t volume = 0.0;
  ForEachSimplex_(
    [&](Element<Dim> const& simplex) {
      volume += simplex.Volume();
    });

  return volume;

} // ComplexElement<...>::Volume

template<size_t Dim>
GVec<Dim> ComplexElement<Dim>::Normal() const {

  if constexpr (Dim == 2) {

    throw std::logic_error("Normal to the 2D complex.");

  } else {

    GVec<Dim> normals(0.0);
    ForEachSimplex_(
      [&](Element<Dim> const& simplex) {
        normals += simplex.Volume()*simplex.Normal();
      });

    return glm::normalize(normals);

  }

} // ComplexElement<...>::Normal

template<size_t Dim>
GVec<Dim> ComplexElement<Dim>::CenterCoords() const {

  real_t volume = 0.0;
  GVec<Dim> centerCoords(0.0);
  ForEachSimplex_(
    [&](Element<Dim> const& simplex) {
      real_t const deltaVolume = simplex.Volume();
      volume += deltaVolume;
      centerCoords += deltaVolume*simplex.CenterCoords();
    });

  return centerCoords/volume;

} // ComplexElement<Dim>::CenterCoords

template<size_t Dim>
real_t Segment<Dim>::Diameter() const {

  GVec<Dim> const delta =
    this->NodeCoords(1) - this->NodeCoords(0);

  return glm::length(delta);

} // Segment<...>::Diameter

template<size_t Dim>
real_t Segment<Dim>::Volume() const {

  return Diameter();

} // Segment<...>::Volume

template<size_t Dim>
GVec<Dim> Segment<Dim>::Normal() const {

  _STORM_NOT_IMPLEMENTED_();
//GVec<Dim> const delta =
//  this->NodeCoords(1) - this->NodeCoords(0);
//
//static GVec<3> const up(0.0, 0.0, 1.0);
//return glm::normalize(glm::cross(delta, up));

} // Segment<...>::Normal

template<size_t Dim>
GVec<Dim> Segment<Dim>::Direction() const {

  GVec<Dim> const delta =
    this->NodeCoords(1) - this->NodeCoords(0);

  return glm::normalize(delta);

} // Segment<...>::Direction

template<size_t Dim>
GVec<Dim> Segment<Dim>::CenterCoords() const {

  return (this->NodeCoords(0) + this->NodeCoords(1))/2.0;

} // Segment<...>::CenterCoords

template<size_t Dim>
ElementDescriptionArray Segment<Dim>::EdgesDesc() const {

  return { this->Part_(ShapeType::Segment, 0, 1) };

} // Segment<...>::EdgesDesc

template<size_t Dim>
ElementDescriptionArray Segment<Dim>::FacesDesc() const {

  return { this->Part_(ShapeType::Node, 0), this->Part_(ShapeType::Node, 1) };

} // Segment<...>::FacesDesc

template class Segment<1>;
template class Segment<2>;
template class Segment<3>;

template<size_t Dim>
real_t Triangle<Dim>::Diameter() const {

  GVec<Dim> const delta1 =
      this->NodeCoords(1) - this->NodeCoords(0);

  GVec<Dim> const delta2 =
      this->NodeCoords(2) - this->NodeCoords(1);

  GVec<Dim> const delta3 =
      this->NodeCoords(0) - this->NodeCoords(2);

  return std::max(
    glm::length(delta1),
    std::max(glm::length(delta2), glm::length(delta3)));

} // Triangle<...>::Diameter

template<size_t Dim>
real_t Triangle<Dim>::Volume() const {

  GVec<Dim> const delta1 =
      this->NodeCoords(1) - this->NodeCoords(0);

  GVec<Dim> const delta2 =
      this->NodeCoords(2) - this->NodeCoords(0);

  if constexpr (Dim == 2) {
    return std::abs(glm::determinant(GMat<Dim>(delta1, delta2)))/2.0;
  } else {
    return glm::length(glm::cross(delta1, delta2))/2.0;
  }

} // Triangle<...>::Volume

template<size_t Dim>
GVec<Dim> Triangle<Dim>::Normal() const {

  if constexpr (Dim == 2) {

    throw std::logic_error("Normal to the 2D triangle.");

  } else {

    GVec<Dim> const delta1 =
        this->NodeCoords(1) - this->NodeCoords(0);

    GVec<Dim> const delta2 =
        this->NodeCoords(2) - this->NodeCoords(0);

    return glm::normalize(glm::cross(delta1, delta2));

  }

} // Triangle<...>::Normal

template<size_t Dim>
GVec<Dim> Triangle<Dim>::CenterCoords() const {

  return (
    this->NodeCoords(0) +
    this->NodeCoords(1) + 
    this->NodeCoords(2))/3.0;

} // Triangle<...>::CenterCoords

template<size_t Dim>
ElementDescriptionArray Triangle<Dim>::EdgesDesc() const {

  return { 
    this->Part_(ShapeType::Segment, 0, 1),
    this->Part_(ShapeType::Segment, 1, 2),
    this->Part_(ShapeType::Segment, 2, 0) };

} // Triangle<...>::EdgesDesc

template<size_t Dim>
ElementDescriptionArray Triangle<Dim>::FacesDesc() const {

  return EdgesDesc();

} // Triangle<...>::FacesDesc

template class Triangle<2>;
template class Triangle<3>;

template<size_t Dim>
ElementDescriptionArray Quadrangle<Dim>::EdgesDesc() const {

  return { 
    this->Part_(ShapeType::Segment, 0, 1),
    this->Part_(ShapeType::Segment, 1, 2),
    this->Part_(ShapeType::Segment, 2, 3),
    this->Part_(ShapeType::Segment, 3, 0) };

} // Quadrangle<...>::EdgesDesc

template<size_t Dim>
ElementDescriptionArray Quadrangle<Dim>::FacesDesc() const {

  return EdgesDesc();

} // Quadrangle<...>::FacesDesc

template<size_t Dim>
ElementDescriptionArray Quadrangle<Dim>::SimplicialParts() const {

  return { 
    this->Part_(ShapeType::Triangle, 0, 1, 2),
    this->Part_(ShapeType::Triangle, 2, 3, 0) };

} // Quadrangle<...>::SimplicialParts

template class Quadrangle<2>;
template class Quadrangle<3>;

template<size_t Dim>
real_t Tetrahedron<Dim>::Diameter() const {

  GVec<Dim> const delta1 =
    this->NodeCoords(1) - this->NodeCoords(0);

  GVec<Dim> const delta2 =
    this->NodeCoords(2) - this->NodeCoords(1);

  GVec<Dim> const delta3 =
    this->NodeCoords(0) - this->NodeCoords(2);

  GVec<Dim> const delta4 =
    this->NodeCoords(3) - this->NodeCoords(0);

  GVec<Dim> const delta5 =
    this->NodeCoords(3) - this->NodeCoords(1);

  GVec<Dim> const delta6 =
    this->NodeCoords(3) - this->NodeCoords(2);

  return std::max({
    glm::length(delta1), glm::length(delta2),
    glm::length(delta3), glm::length(delta4),
    glm::length(delta5), glm::length(delta6) });

} // Tetrahedron<...>::Diameter

template<size_t Dim>
real_t Tetrahedron<Dim>::Volume() const {

  GVec<Dim> const delta1 =
    this->NodeCoords(1) - this->NodeCoords(0);

  GVec<Dim> const delta2 =
    this->NodeCoords(2) - this->NodeCoords(0);

  GVec<Dim> const delta3 =
    this->NodeCoords(3) - this->NodeCoords(0);

  return std::abs(glm::dot(delta1, glm::cross(delta2, delta3)))/6.0;

} // Tetrahedron<...>::Volume

template<size_t Dim>
GVec<Dim> Tetrahedron<Dim>::CenterCoords() const {

  return (
    this->NodeCoords(0) + this->NodeCoords(1) +
    this->NodeCoords(2) + this->NodeCoords(3))/4.0;

} // Tetrahedron<...>::CenterCoords

template<size_t Dim>
ElementDescriptionArray Tetrahedron<Dim>::EdgesDesc() const {

  return {
    this->Part_(ShapeType::Segment, 0, 1),
    this->Part_(ShapeType::Segment, 1, 2),
    this->Part_(ShapeType::Segment, 2, 0),
    this->Part_(ShapeType::Segment, 0, 3),
    this->Part_(ShapeType::Segment, 1, 3),
    this->Part_(ShapeType::Segment, 2, 3) };

} // Tetrahedron<...>::EdgesDesc

template<size_t Dim>
ElementDescriptionArray Tetrahedron<Dim>::FacesDesc() const {

  return {
    this->Part_(ShapeType::Triangle, 0, 2, 1),
    this->Part_(ShapeType::Triangle, 0, 1, 3),
    this->Part_(ShapeType::Triangle, 1, 2, 3),
    this->Part_(ShapeType::Triangle, 2, 0, 3) };

} // Tetrahedron<...>::FacesDesc

template class Tetrahedron<3>;

template<size_t Dim>
ElementDescriptionArray Pyramid<Dim>::EdgesDesc() const {

  return {
    this->Part_(ShapeType::Segment, 0, 1),
    this->Part_(ShapeType::Segment, 1, 2),
    this->Part_(ShapeType::Segment, 2, 3),
    this->Part_(ShapeType::Segment, 3, 0),
    this->Part_(ShapeType::Segment, 0, 4),
    this->Part_(ShapeType::Segment, 1, 4),
    this->Part_(ShapeType::Segment, 2, 4),
    this->Part_(ShapeType::Segment, 3, 4) };

} // Pyramid<...>::EdgesDesc

template<size_t Dim>
ElementDescriptionArray Pyramid<Dim>::FacesDesc() const {

  return {
    this->Part_(ShapeType::Quadrangle, 0, 3, 2, 1),
    this->Part_(ShapeType::Triangle,   0, 1, 4   ),
    this->Part_(ShapeType::Triangle,   1, 2, 4   ),
    this->Part_(ShapeType::Triangle,   2, 3, 4   ),
    this->Part_(ShapeType::Triangle,   3, 0, 4   ) };

} // Pyramid<...>::FacesDesc

template<size_t Dim>
ElementDescriptionArray Pyramid<Dim>::SimplicialParts() const {

  return {
    this->Part_(ShapeType::Tetrahedron, 0, 1, 2, 4),
    this->Part_(ShapeType::Tetrahedron, 2, 3, 0, 4) };

} // Pyramid<...>::SimplicialParts

template class Pyramid<3>;

template<size_t Dim>
ElementDescriptionArray Prism<Dim>::EdgesDesc() const {

  return {
    this->Part_(ShapeType::Segment, 0, 1),
    this->Part_(ShapeType::Segment, 1, 2),
    this->Part_(ShapeType::Segment, 2, 0),
    this->Part_(ShapeType::Segment, 0, 3),
    this->Part_(ShapeType::Segment, 1, 4),
    this->Part_(ShapeType::Segment, 2, 5),
    this->Part_(ShapeType::Segment, 3, 4),
    this->Part_(ShapeType::Segment, 4, 5),
    this->Part_(ShapeType::Segment, 5, 3) };

} // Prism<...>::EdgesDesc

template<size_t Dim>
ElementDescriptionArray Prism<Dim>::FacesDesc() const {

  return {
    this->Part_(ShapeType::Quadrangle, 0, 1, 4, 3),
    this->Part_(ShapeType::Quadrangle, 1, 2, 5, 4),
    this->Part_(ShapeType::Quadrangle, 2, 0, 3, 5),
    this->Part_(ShapeType::Triangle,   0, 2, 1   ),
    this->Part_(ShapeType::Triangle,   3, 4, 5   ) };

} // Prism<...>::FacesDesc

template<size_t Dim>
ElementDescriptionArray Prism<Dim>::SimplicialParts() const {

  _STORM_NOT_IMPLEMENTED_();

} // Prism<...>::SimplicialParts

template class Prism<3>;

template<size_t Dim>
ElementDescriptionArray Hexahedron<Dim>::EdgesDesc() const {

  return {
    this->Part_(ShapeType::Segment, 0, 1),
    this->Part_(ShapeType::Segment, 1, 2),
    this->Part_(ShapeType::Segment, 2, 3),
    this->Part_(ShapeType::Segment, 3, 0),
    this->Part_(ShapeType::Segment, 0, 4),
    this->Part_(ShapeType::Segment, 1, 5),
    this->Part_(ShapeType::Segment, 2, 6),
    this->Part_(ShapeType::Segment, 3, 7),
    this->Part_(ShapeType::Segment, 4, 5),
    this->Part_(ShapeType::Segment, 5, 6),
    this->Part_(ShapeType::Segment, 6, 7),
    this->Part_(ShapeType::Segment, 7, 4) };

} // Hexahedron<...>::EdgesDesc

template<size_t Dim>
ElementDescriptionArray Hexahedron<Dim>::FacesDesc() const {

  return {
    this->Part_(ShapeType::Quadrangle, 0, 3, 2, 1),
    this->Part_(ShapeType::Quadrangle, 0, 1, 5, 4),
    this->Part_(ShapeType::Quadrangle, 1, 2, 6, 5),
    this->Part_(ShapeType::Quadrangle, 2, 3, 7, 6),
    this->Part_(ShapeType::Quadrangle, 0, 4, 7, 3),
    this->Part_(ShapeType::Quadrangle, 4, 5, 6, 7) };

} // Hexahedron<...>::FacesDesc

template<size_t Dim>
ElementDescriptionArray Hexahedron<Dim>::SimplicialParts() const {

  return {
    this->Part_(ShapeType::Tetrahedron, 0, 3, 1, 4),
    this->Part_(ShapeType::Tetrahedron, 3, 2, 1, 6),
    this->Part_(ShapeType::Tetrahedron, 4, 5, 6, 1),
    this->Part_(ShapeType::Tetrahedron, 4, 6, 7, 3),
    this->Part_(ShapeType::Tetrahedron, 4, 3, 1, 6) };

} // Hexahedron<...>::SimplicialParts

template class Hexahedron<3>;

#if !STORM_DOXYGEN
/// Alternative pyramid sketches:
/// @verbatim
///                            n4
///                            O    f1
///                           /|\   ^
///                          /`|`\./
///                        ,// | `/
///                       /,/  | /\.
///                     ,/ |`  |/ `\
///          f3     e7 ^` ,|   /   ^ e5
///           ^      ,/   /`  /|   \.
///            \    /`   /   * |   `\
///             \ ,/   ,/      |    |
///              /`    |`   e6 ^    \.
///            ,/ \    ^ e4    |  o-`|---------> f2
///           /`   o  ,|       |     |
///   <----------*    /`       |     \.
///       ,/`        /         |     `\
///   n3 O----<----,/----------O n2   |
///       \.  e2   |`           \.    \.
///        `\.    ,|          e1 `^.  `\
///      e3 `v.   /`     o         `\. |
///           `\ /       |            `\|
///         n0 O------------------->---O n1
///                      |        e0
///                      |
///                      v
///                      f0
/// @endverbatim
///
/// @verbatim
///                           f2
///                     n4    ^
///                     O.   /
///                   ,/|`\./
///              e7 ,^/ |  /\.
///               ,/ /  | *  ^ e6
///             ,/`,/   ^ e5 `\.
///           ,/` /`    |      `\
///       n3 O---/------|---<---O n2
///         /   ^ e4    |  e2 ,/`
///        /  ,/ *      |   ,/`
///    e3 v /   / o     | ,^` e1
///      /,/   /  |     |/`
///  n0 O-----/----->---O n1
///          /    | e0
///         /     |
///        v      v
///       f1      f0
/// @endverbatim
#endif

} // namespace Storm
