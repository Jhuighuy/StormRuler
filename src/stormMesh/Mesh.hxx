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

#include <tuple>
#include <utility>
#include <ranges>

#include <stormBase.hxx>
#include <stormUtils/TaggedIndex.hxx>
#include <stormUtils/Array.hxx>

// Use a different assert macro to debug the
// mesh algorithms separate from the others.
#ifdef _STORM_MESH_DEBUG_
#define StormMeshAssert StormEnabledAssert
#else
#define StormMeshAssert StormDisabledAssert
#endif

namespace Storm {

class NodeTag_;
class EdgeTag_;
class FaceTag_;
class CellTag_;
template<class Tag>
class MarkTag_;

template<size_t Dim>
class Mesh;

template<size_t Dim>
class NodeView;
template<size_t Dim>
class EdgeView;
template<size_t Dim>
class FaceView;
template<size_t Dim>
class CellView;

/// @brief Index of a node in the mesh.
using NodeIndex = TaggedIndex<NodeTag_>;

/// @brief Index of an edge in the mesh.
using EdgeIndex = TaggedIndex<EdgeTag_>;

/// @brief Index of a face in the mesh.
using FaceIndex = TaggedIndex<FaceTag_>;

/// @brief Index of a cell in the mesh.
using CellIndex = TaggedIndex<CellTag_>;

/// @brief Index of a node mark in the mesh.
using NodeMarkIndex = TaggedIndex<MarkTag_<NodeTag_>>;

/// @brief Index of an edge mark in the mesh.
using EdgeMarkIndex = TaggedIndex<MarkTag_<EdgeTag_>>;

/// @brief Index of a face mark in the mesh.
using FaceMarkIndex = TaggedIndex<MarkTag_<FaceTag_>>;

/// @brief Index of a cell mark in the mesh.
using CellMarkIndex = TaggedIndex<MarkTag_<CellTag_>>;

/// ----------------------------------------------------------------- ///
/// @brief Shape of a mesh element.
/// ----------------------------------------------------------------- ///
enum class Shape {

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

}; // enum class Shape

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Immutable hybrid unstructured multidimensional mesh.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<size_t Dim>
class Mesh : public BaseObject {
private:

  template<size_t>
  friend class NodeView;
  template<size_t>
  friend class EdgeView;
  template<size_t>
  friend class FaceView;
  template<size_t>
  friend class CellView;

protected:

  NodeIndex NumAllNodes_;
  EdgeIndex NumAllEdges_;
  FaceIndex NumAllFaces_;
  CellIndex NumAllCells_;
  Array<NodeIndex, NodeMarkIndex> NodeRanges_;
  Array<EdgeIndex, EdgeMarkIndex> EdgeRanges_;
  Array<FaceIndex, FaceMarkIndex> FaceRanges_;
  Array<CellIndex, CellMarkIndex> CellRanges_;
  Table<void, NodeIndex, NodeIndex> NodeNodes_;
  Table<void, NodeIndex, EdgeIndex> NodeEdges_;
  Table<void, NodeIndex, FaceIndex> NodeFaces_;
  Table<void, NodeIndex, CellIndex> NodeCells_;
  Array<std::pair<NodeIndex, NodeIndex>, EdgeIndex> EdgeNodes_;
  Table<void, EdgeIndex, EdgeIndex> EdgeEdges_;
  Table<void, EdgeIndex, FaceIndex> EdgeFaces_;
  Table<void, EdgeIndex, CellIndex> EdgeCells_;
  Table<void, FaceIndex, NodeIndex> FaceNodes_;
  Table<void, FaceIndex, EdgeIndex> FaceEdges_;
  Table<void, FaceIndex, FaceIndex> FaceFaces_;
  Array<std::pair<CellIndex, CellIndex>, FaceIndex> FaceCells_;
  Table<void, CellIndex, NodeIndex> CellNodes_;
  Table<void, CellIndex, EdgeIndex> CellEdges_;
  Table<void, CellIndex, FaceIndex> CellFaces_;
  Table<void, CellIndex, CellIndex> CellCells_;
  Array<Vec<Dim>, NodeIndex> NodeCoords_;
  Array<std::tuple<Vec<Dim>, real_t>, EdgeIndex> 
                                EdgeDirsAndLengths_;
  Array<std::tuple<Vec<Dim>, Vec<Dim>, real_t>, FaceIndex> 
                          FaceCenterCoordsAndNormalsAndAreas_;
  Array<std::tuple<Vec<Dim>, real_t>, CellIndex> 
                        CellCenterCoordsAndVolumes_;
  Array<Shape, FaceIndex> FaceShapes_;
  Array<Shape, CellIndex> CellShapes_;

  template<class View, class Range>
  auto ViewRange_(Range range) const noexcept {
    return std::views::transform(range,
      [this](auto index){ return View(this, index); });
  }
  template<class View, class Index>
  auto IotaViewRange_(Index first, Index last) const noexcept {
    return ViewRange_<View>(std::views::iota(first, last));
  }

public:

  /// @brief Total number of nodes in the mesh.
  NodeIndex NumAllNodes() const noexcept {
    return NumAllNodes_;
  }

  /// @brief Total number of edges in the mesh.
  EdgeIndex NumAllEdges() const noexcept {
    return NumAllEdges_;
  }

  /// @brief Total number of faces in the mesh.
  FaceIndex NumAllFaces() const noexcept {
    return NumAllFaces_;
  }

  /// @brief Total number of cells in the mesh.
  CellIndex NumAllCells() const noexcept {
    return NumAllCells_;
  }

  /// @brief Get node at the index.
  auto Node(NodeIndex nodeIndex) const noexcept {
    StormMeshAssert(nodeIndex < NumAllNodes_);
    return NodeView<Dim>(this, nodeIndex);
  }

  /// @brief Get edge at the index.
  auto Edge(EdgeIndex edgeIndex) const noexcept {
    StormMeshAssert(edgeIndex < NumAllEdges_);
    return EdgeView<Dim>(this, edgeIndex);
  }

  /// @brief Get face at the index.
  auto Face(FaceIndex faceIndex) const noexcept {
    StormMeshAssert(faceIndex < NumAllFaces_);
    return FaceView<Dim>(this, faceIndex);
  }

  /// @brief Get cell at the index.
  auto Cell(CellIndex cellIndex) const noexcept {
    StormMeshAssert(cellIndex < NumAllCells_);
    return CellView<Dim>(this, cellIndex);
  }

  /// @brief Number of node marks in the mesh.
  NodeMarkIndex NumNodeMarks() const noexcept {
    return NodeRanges_.Size() - 1;
  }

  /// @brief Number of edge marks in the mesh.
  EdgeMarkIndex NumEdgeMarks() const noexcept {
    return EdgeRanges_.Size() - 1;
  }

  /// @brief Number of face marks in the mesh.
  FaceMarkIndex NumFaceMarks() const noexcept {
    return FaceRanges_.Size() - 1;
  }

  /// @brief Number of cell marks in the mesh.
  CellMarkIndex NumCellMarks() const noexcept {
    return CellRanges_.Size() - 1;
  }

  /// @brief Range of nodes with the specified mark.
  auto Nodes(NodeMarkIndex nodeMark = 0) const noexcept {
    StormMeshAssert(nodeMark < NumNodeMarks());
    return IotaViewRange_<NodeView<Dim>>(
      NodeRanges_[nodeMark], NodeRanges_[nodeMark + 1]);
  }

  /// @brief Range of edges with the specified mark.
  auto Edges(EdgeMarkIndex edgeMark = 0) const noexcept {
    StormMeshAssert(edgeMark < NumEdgeMarks());
    return IotaViewRange_<EdgeView<Dim>>(
      EdgeRanges_[edgeMark], EdgeRanges_[edgeMark + 1]);
  }

  /// @brief Range of faces with the specified mark.
  auto Faces(FaceMarkIndex faceMark = 0) const noexcept {
    StormMeshAssert(faceMark < NumFaceMarks());
    return IotaViewRange_<FaceView<Dim>>(
      FaceRanges_[faceMark], FaceRanges_[faceMark + 1]);
  }

  /// @brief Range of faces with the specified mark.
  auto Cells(CellMarkIndex cellMark = 0) const noexcept {
    StormMeshAssert(cellMark < NumCellMarks());
    return IotaViewRange_<CellView<Dim>>(
      CellRanges_[cellMark], CellRanges_[cellMark + 1]);
  }

}; // class Mesh<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief View to a node of the unstructured mesh.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<size_t Dim>
class NodeView final {
private:
  Mesh<Dim> const* Mesh_;
  NodeIndex NodeIndex_;

  void SelfCheck_() const noexcept {
    StormMeshAssert(
      Mesh_ != nullptr && NodeIndex_ < Mesh_->NumAllNodes());
  }

public:

  /// @brief Construct pointer to a node.
  explicit NodeView(Mesh<Dim> const* mesh,
                    NodeIndex nodeIndex) noexcept :
      Mesh_{mesh}, NodeIndex_{nodeIndex} {
    SelfCheck_();
  }

  /// @brief Coordinates of the current node.
  Vec<Dim> const& Coords() const noexcept {
    SelfCheck_();
    return Mesh_->NodeCoords_[NodeIndex_];
  }

  /// @brief Range of nodes that share an edge with the current node.
  auto Nodes() const noexcept {
    SelfCheck_();
    return Mesh_->template ViewRange_<NodeView<Dim>>(Mesh_->NodeNodes_[NodeIndex_]);
  }

  /// @brief Range of edges that contain the current node.
  auto Edges() const noexcept {
    SelfCheck_();
    return Mesh_->template ViewRange_<EdgeView<Dim>>(Mesh_->NodeEdges_[NodeIndex_]);
  }

  /// @brief Range of faces that contain the current node.
  auto Faces() const noexcept {
    SelfCheck_();
    return Mesh_->template ViewRange_<FaceView<Dim>>(Mesh_->NodeFaces_[NodeIndex_]);
  }

  /// @brief Range of cells that contain the current node.
  auto Cells() const noexcept {
    SelfCheck_();
    return Mesh_->template ViewRange_<CellView<Dim>>(Mesh_->NodeCells_[NodeIndex_]);
  }

}; // class NodeView<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief View to an edge of the unstructured mesh.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<size_t Dim>
class EdgeView final {
private:
  Mesh<Dim> const* Mesh_;
  EdgeIndex EdgeIndex_;

  void SelfCheck_() const noexcept {
    StormMeshAssert(
      Mesh_ != nullptr && EdgeIndex_ < Mesh_->NumAllEdges());
  }

public:

  /// @brief Construct pointer to an edge.
  explicit EdgeView(Mesh<Dim> const* mesh,
                    EdgeIndex edgeIndex) noexcept :
      Mesh_{mesh}, EdgeIndex_{edgeIndex} {
    SelfCheck_();
  }

  /// @brief Unit direction to the current edge.
  Vec<Dim> const& Direction() const noexcept {
    SelfCheck_();
    return std::get<0>(Mesh_->EdgeDirsAndLengths_[EdgeIndex_]);
  }

  /// @brief Length of the current edge.
  real_t Length() const noexcept {
    SelfCheck_();
    return std::get<1>(Mesh_->EdgeDirsAndLengths_[EdgeIndex_]);
  }

  /// @brief Pair of nodes of the edge.
  auto NodesPair() const noexcept {
    SelfCheck_();
    auto const& nodeIndices = Mesh_->EdgeNodes_[EdgeIndex_];
    return std::make_pair(
      Mesh_->Node(nodeIndices.first), Mesh_->Node(nodeIndices.second));
  }

  /// @brief Ranges of edges that share a node with the current edge.
  auto Edges() const noexcept {
    SelfCheck_();
    return Mesh_->template ViewRange_<EdgeView<Dim>>(Mesh_->EdgeEdges_[EdgeIndex_]);
  }

  /// @brief Range of faces that contain the current edge.
  auto Faces() const noexcept {
    SelfCheck_();
    return Mesh_->template ViewRange_<FaceView<Dim>>(Mesh_->EdgeFaces_[EdgeIndex_]);
  }

  /// @brief Range of cells that contain the current cell.
  auto Cells() const noexcept {
    SelfCheck_();
    return Mesh_->template ViewRange_<CellView<Dim>>(Mesh_->EdgeCells_[EdgeIndex_]);
  }

}; // class EdgeView<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Pointer to a face of the unstructured mesh.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<size_t Dim>
class FaceView final {
private:
  Mesh<Dim> const* Mesh_;
  FaceIndex FaceIndex_;

  void SelfCheck_() const noexcept {
    StormMeshAssert(
      Mesh_ != nullptr && FaceIndex_ < Mesh_->NumAllFaces());
  }

public:

  /// @brief Construct pointer to a face.
  explicit FaceView(Mesh<Dim> const* mesh,
                    FaceIndex faceIndex) noexcept :
      Mesh_{mesh}, FaceIndex_{faceIndex} {
    SelfCheck_();
  }

  /// @brief Center coordinates of the current face.
  Vec<Dim> const& CenterCoords() const noexcept {
    SelfCheck_();
    return std::get<0>(
      Mesh_->FaceCenterCoordsAndNormalsAndAreas_[FaceIndex_]);
  }

  /// @brief Unit normal to the current face.
  Vec<Dim> const& Normal() const noexcept {
    SelfCheck_();
    return std::get<1>(
      Mesh_->FaceCenterCoordsAndNormalsAndAreas_[FaceIndex_]);
  }

  /// @brief Area of the current face.
  real_t Area() const noexcept {
    SelfCheck_();
    return std::get<2>(
      Mesh_->FaceCenterCoordsAndNormalsAndAreas_[FaceIndex_]);
  }

  /// @brief Range of all nodes of the face.
  auto Nodes() const noexcept {
    SelfCheck_();
    return Mesh_->template ViewRange_<NodeView<Dim>>(Mesh_->FaceNodes_[FaceIndex_]);
  }

  /// @brief Range of all edges of the face.
  auto Edges() const noexcept {
    SelfCheck_();
    return Mesh_->template ViewRange_<EdgeView<Dim>>(Mesh_->FaceEdges_[FaceIndex_]);
  }

  /// @brief Range of all faces that share an edge with the current face.
  auto Faces() const noexcept {
    SelfCheck_();
    return Mesh_->template ViewRange_<FaceView<Dim>>(Mesh_->FaceFaces_[FaceIndex_]);
  }

  /// @brief Get the inner and outer cells of the face.
  auto CellsPair() const noexcept {
    SelfCheck_();
    auto const& cellIndices = Mesh_->FaceCells_[FaceIndex_];
    return std::make_pair(
      Mesh_->Cell(cellIndices.first), Mesh_->Cell(cellIndices.second));
  }

}; // class FaceView<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Pointer to a cell of the unstructured mesh.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<size_t Dim>
class CellView final {
private:
  Mesh<Dim> const* Mesh_;
  CellIndex CellIndex_;

  void SelfCheck_() const noexcept {
    StormMeshAssert(
      Mesh_ != nullptr && CellIndex_ < Mesh_->NumAllCells());
  }

public:

  /// @brief Construct pointer to a cell.
  explicit CellView(Mesh<Dim> const* mesh,
                    CellIndex cellIndex) noexcept :
      Mesh_{mesh}, CellIndex_{cellIndex} {
    SelfCheck_();
  }

  /// @brief Center coordinates of the current cell.
  Vec<Dim> const& CenterCoords() const noexcept {
    SelfCheck_();
    return std::get<0>(Mesh_->CellCenterCoordsAndVolumes_[CellIndex_]);
  }

  /// @brief Volume of the current cell.
  real_t Volume() const noexcept {
    SelfCheck_();
    return std::get<1>(Mesh_->CellCenterCoordsAndVolumes_[CellIndex_]);
  }

  /// @brief Range of all nodes of the cell.
  auto Nodes() const noexcept {
    SelfCheck_();
    return Mesh_->template ViewRange_<NodeView<Dim>>(Mesh_->CellNodes_[CellIndex_]);
  }

  /// @brief Range of all edges of the cell.
  auto Edges() const noexcept {
    SelfCheck_();
    return Mesh_->template ViewRange_<EdgeView<Dim>>(Mesh_->CellEdges_[CellIndex_]);
  }

  /// @brief Range of all faces of the cell.
  auto Faces() const noexcept {
    SelfCheck_();
    return Mesh_->template ViewRange_<FaceView<Dim>>(Mesh_->CellFaces_[CellIndex_]);
  }

  /// @brief Range of all cells that share a face with the current cell.
  auto Cells() const noexcept {
    SelfCheck_();
    return Mesh_->template ViewRange_<CellView<Dim>>(Mesh_->CellCells_[CellIndex_]);
  }

}; // class CellView<...>

} // namespace Storm
