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

#include <stormBase.hxx>
#include <stormUtils/Glm.hxx>
#include <stormUtils/TaggedIndex.hxx>

// Use a different assert macro to debug the
// mesh algorithms separate from the others.
#ifdef _STORM_MESH_DEBUG_
#define StormMeshAssert StormEnabledAssert
#else
#define StormMeshAssert StormDisabledAssert
#endif

namespace Storm {

template<size_t Dim>
class NodeView;
template<size_t Dim>
class EdgeView;
template<size_t Dim>
class FaceView;
template<size_t Dim>
class CellView;

/// ----------------------------------------------------------------- ///
/// @brief Shape of a mesh element.
/// ----------------------------------------------------------------- ///
enum class ShapeType {

  /// @brief Undefined shape.
  Undefined,

  /// @brief Node shape.
  Node,

  /// @brief Segment shape.
  Segment,

  /// @brief Triangle shape.
  Triangle,

  /// @brief Quadrangle shape.
  Quadrangle,

  /// @brief Tetrahedron shape.
  Tetrahedron,

  /// @brief Pyramid shape.
  Pyramid,

  /// @brief Prism (or wedge, penthedron) shape.
  Prism,

  /// @brief Hexahedron shape.
  Hexahedron,

}; // enum class ShapeType

/// ----------------------------------------------------------------- ///
/// @brief Element description.
/// ----------------------------------------------------------------- ///
struct ElementDescription {
  ShapeType Shape;
  std::vector<size_t> NodeIndices;
}; // struct ElementDescription

/// ----------------------------------------------------------------- ///
/// @brief Array of element descriptions.
/// ----------------------------------------------------------------- ///
using ElementDescriptionArray = std::vector<ElementDescription>;

} // namespace Storm
