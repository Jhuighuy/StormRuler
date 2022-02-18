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
#include <algorithm>

#include <stormBase.hxx>
#include <stormUtils/TaggedIndex.hxx>

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

template<class Mesh, class View, class PointerOrIndex>
class ElementIterator;

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

/// @brief Node iterator.
template<size_t Dim>
using NodeIterator = 
  ElementIterator<Mesh<Dim>, NodeView<Dim>, NodeIndex>;

/// @brief Edge iterator.
template<size_t Dim>
using EdgeIterator = 
  ElementIterator<Mesh<Dim>, EdgeView<Dim>, EdgeIndex>;

/// @brief Face iterator.
template<size_t Dim>
using FaceIterator = 
  ElementIterator<Mesh<Dim>, FaceView<Dim>, FaceIndex>;

/// @brief Cell iterator.
template<size_t Dim>
using CellIterator = 
  ElementIterator<Mesh<Dim>, CellView<Dim>, CellIndex>;

/// @brief Node local iterator.
template<size_t Dim>
using NodeLocalIterator = 
  ElementIterator<Mesh<Dim>, NodeView<Dim>, NodeIndex const*>;

/// @brief Edge local iterator.
template<size_t Dim>
using EdgeLocalIterator = 
  ElementIterator<Mesh<Dim>, EdgeView<Dim>, EdgeIndex const*>;

/// @brief Face local iterator.
template<size_t Dim>
using FaceLocalIterator = 
  ElementIterator<Mesh<Dim>, FaceView<Dim>, FaceIndex const*>;

/// @brief Cell local iterator.
template<size_t Dim>
using CellLocalIterator = 
  ElementIterator<Mesh<Dim>, CellView<Dim>, CellIndex const*>;

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

  /// @brief Iterator pointing to \
  ///   the first node with the specified mark.
  auto BeginNode(NodeMarkIndex nodeMark = 0) const noexcept {
    StormMeshAssert(nodeMark < NumNodeMarks());
    return NodeIterator<Dim>(this, NodeRanges_[nodeMark]);
  }

  /// @brief Iterator pointing to \
  ///   the node following the last node with the specified mark.
  auto EndNode(NodeMarkIndex nodeMark = 0) const noexcept {
    StormMeshAssert(nodeMark < NumNodeMarks());
    return NodeIterator<Dim>(this, NodeRanges_[nodeMark + 1]);
  }

  /// @brief Iterator pointing to \
  ///   the first edge with the specified mark.
  auto BeginEdge(EdgeMarkIndex edgeMark = 0) const noexcept {
    StormMeshAssert(edgeMark < NumEdgeMarks());
    return EdgeIterator<Dim>(this, EdgeRanges_[edgeMark]);
  }

  /// @brief Iterator pointing to \
  ///   the edge following the last edge with the specified mark.
  auto EndEdge(EdgeMarkIndex edgeMark = 0) const noexcept {
    StormMeshAssert(edgeMark < NumEdgeMarks());
    return EdgeIterator<Dim>(this, EdgeRanges_[edgeMark + 1]);
  }

  /// @brief Iterator pointing to \
  ///   the first face with the specified mark.
  auto BeginFace(FaceMarkIndex faceMark = 0) const noexcept {
    StormMeshAssert(faceMark < NumFaceMarks());
    return FaceIterator<Dim>(this, FaceRanges_[faceMark]);
  }

  /// @brief Iterator pointing to \
  ///   the face following the last face with the specified mark.
  auto EndFace(FaceMarkIndex faceMark = 0) const noexcept {
    StormMeshAssert(faceMark < NumFaceMarks());
    return FaceIterator<Dim>(this, FaceRanges_[faceMark + 1]);
  }

  /// @brief Iterator pointing to \
  ///   the first cell with the specified mark.
  auto BeginCell(CellMarkIndex cellMark = 0) const noexcept {
    StormMeshAssert(cellMark < NumCellMarks());
    return CellIterator<Dim>(this, CellRanges_[cellMark]);
  }

  /// @brief Iterator pointing to \
  ///   the cell following after the last face with a given mark.
  auto EndCell(CellMarkIndex cellMark = 0) const noexcept {
    StormMeshAssert(cellMark < NumCellMarks());
    return CellIterator<Dim>(this, CellRanges_[cellMark + 1]);
  }

  /// @brief Number of nodes with the specified mark. 
  NodeIndex NumNodes(NodeMarkIndex nodeMark = 0) const noexcept {
    return EndNode(nodeMark) - BeginNode(nodeMark);
  }

  /// @brief Number of edges with the specified mark. 
  EdgeIndex NumEdges(EdgeMarkIndex edgeMark = 0) const noexcept {
    return EndEdge(edgeMark) - BeginEdge(edgeMark);
  }

  /// @brief Number of faces with the specified mark. 
  FaceIndex NumFaces(FaceMarkIndex faceMark = 0) const noexcept {
    return EndFace(faceMark) - BeginFace(faceMark);
  }

  /// @brief Number of cells with the specified mark. 
  CellIndex NumCells(CellMarkIndex cellMark = 0) const noexcept {
    return EndCell(cellMark) - BeginCell(cellMark);
  }

}; // class Mesh<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief View to a node of the unstructured mesh.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<size_t Dim>
class NodeView {
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

  /// @brief Iterator pointing to \
  ///   the first node that shares an edge with the current node.
  auto BeginNode() const noexcept {
    SelfCheck_();
    return NodeLocalIterator<Dim>(
      Mesh_, Mesh_->NodeNodes_[NodeIndex_].Begin());
  }

  /// @brief Iterator pointing to the node following \
  ///   the last node that shares an edge with the current node.
  auto EndNode() const noexcept {
    SelfCheck_();
    return NodeLocalIterator<Dim>(
      Mesh_, Mesh_->NodeNodes_[NodeIndex_].End());
  }

  /// @brief Iterator pointing to \
  ///   the first edge that contains the current node.
  auto BeginEdge() const noexcept {
    SelfCheck_();
    return EdgeLocalIterator<Dim>(
      Mesh_, Mesh_->NodeEdges_[NodeIndex_].Begin());
  }

  /// @brief Iterator pointing to the edge following \
  ///   the last edge that contains the current node.
  auto EndEdge() const noexcept {
    SelfCheck_();
    return EdgeLocalIterator<Dim>(
      Mesh_, Mesh_->NodeEdges_[NodeIndex_].End());
  }

  /// @brief Iterator pointing to \
  ///   the first face that contains the current node.
  auto BeginFace() const noexcept {
    SelfCheck_();
    return FaceLocalIterator<Dim>(
      Mesh_, Mesh_->NodeFaces_[NodeIndex_].Begin());
  }

  /// @brief Iterator pointing to the face following \
  ///   the last face that contains the current node.
  auto EndFace() const noexcept {
    SelfCheck_();
    return FaceLocalIterator<Dim>(
      Mesh_, Mesh_->NodeFaces_[NodeIndex_].End());
  }

  /// @brief Iterator pointing to \
  ///   the first cell that contains the current node.
  auto BeginCell() const noexcept {
    SelfCheck_();
    return CellLocalIterator<Dim>(
      Mesh_, Mesh_->NodeCells_[NodeIndex_].Begin());
  }

  /// @brief Iterator pointing to the cell following \
  ///   the last cell that contains the current node.
  auto EndCell() const noexcept {
    SelfCheck_();
    return CellLocalIterator<Dim>(
      Mesh_, Mesh_->NodeCells_[NodeIndex_].End());
  }

  /// @brief Iterate through all nodes \
  ///   that shares an edge with the current node.
  template<class Func>
  void ForEachNode(Func const& nodeFunc) const noexcept {
    std::for_each(BeginNode(), EndNode(), nodeFunc);
  }

  /// @brief Iterate through all edges \
  ///   that contain the current node.
  template<class Func>
  void ForEachEdge(Func const& edgeFunc) const noexcept {
    std::for_each(BeginEdge(), EndEdge(), edgeFunc);
  }

  /// @brief Iterate through all faces \
  ///   that contain the current node.
  template<class Func>
  void ForEachFace(Func const& faceFunc) const noexcept {
    std::for_each(BeginFace(), EndFace(), faceFunc);
  }

  /// @brief Iterate through all cells \
  ///   that contain the current node.
  template<class Func>
  void ForEachCell(Func const& cellFunc) const noexcept {
    std::for_each(BeginCell(), EndCell(), cellFunc);
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

  /// @brief Get the nodes of the edge.
  auto Nodes() const noexcept {
    SelfCheck_();
    auto const& nodeIndices = Mesh_->EdgeNodes_[EdgeIndex_];
    return std::make_pair(
      Mesh_->Node(nodeIndices.first), Mesh_->Node(nodeIndices.second));
  }

  /// @brief Iterator pointing to \
  ///   the first edge that shares a node with the current edge.
  auto BeginEdge() const noexcept {
    SelfCheck_();
    return EdgeLocalIterator<Dim>(
      Mesh_, Mesh_->EdgeEdges_[EdgeIndex_].Begin());
  }

  /// @brief Iterator pointing to the edge following \
  ///   the last edge that shares a node with the current edge.
  auto EndEdge() const noexcept {
    SelfCheck_();
    return EdgeLocalIterator<Dim>(
      Mesh_, Mesh_->EdgeEdges_[EdgeIndex_].End());
  }

  /// @brief Iterator pointing to \
  ///   the first face that contains the current edge.
  auto BeginFace() const noexcept {
    SelfCheck_();
    return FaceLocalIterator<Dim>(
      Mesh_, Mesh_->EdgeFaces_[EdgeIndex_].Begin());
  }

  /// @brief Iterator pointing to the face following \
  ///   the last face that contains the current edge.
  auto EndFace() const noexcept {
    SelfCheck_();
    return FaceLocalIterator<Dim>(
      Mesh_, Mesh_->EdgeFaces_[EdgeIndex_].End());
  }

  /// @brief Iterator pointing to \
  ///   the first cell that contains the current edge.
  auto BeginCell() const noexcept {
    SelfCheck_();
    return CellLocalIterator<Dim>(
      Mesh_, Mesh_->EdgeCells_[EdgeIndex_].Begin());
  }

  /// @brief Iterator pointing to the cell following \
  ///   the last cell that contains the current edge.
  auto EndCell() const noexcept {
    SelfCheck_();
    return CellLocalIterator<Dim>(
      Mesh_, Mesh_->EdgeCells_[EdgeIndex_].End());
  }

  /// @brief Iterate through all edges \
  ///   that share a node with the current edge.
  template<class Func>
  void ForEachEdge(Func const& edgeFunc) const noexcept {
    std::for_each(BeginEdge(), EndEdge(), edgeFunc);
  }

  /// @brief Iterate through all faces \
  ///   that contain the current edge.
  template<class Func>
  void ForEachFace(Func const& faceFunc) const noexcept {
    std::for_each(BeginFace(), EndFace(), faceFunc);
  }

  /// @brief Iterate through all cells \
  ///   that contain the current cell.
  template<class Func>
  void ForEachCell(Func const& cellFunc) const noexcept {
    std::for_each(BeginCell(), EndCell(), cellFunc);
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

  /// @brief Iterator pointing to \
  ///   the first node of the current face.
  auto BeginNode() const noexcept {
    SelfCheck_();
    return NodeLocalIterator<Dim>(
      Mesh_, Mesh_->FaceNodes_[FaceIndex_].Begin());
  }

  /// @brief Iterator pointing to the node following \
  ///   the last node of the current face.
  auto EndNode() const noexcept {
    SelfCheck_();
    return NodeLocalIterator<Dim>(
      Mesh_, Mesh_->FaceNodes_[FaceIndex_].End());
  }

  /// @brief Iterator pointing to \
  ///   the first edge of the current face.
  auto BeginEdge() const noexcept {
    SelfCheck_();
    return EdgeLocalIterator<Dim>(
      Mesh_, Mesh_->FaceEdges_[FaceIndex_].Begin());
  }

  /// @brief Iterator pointing to the edge following \
  ///   the last edge of the current face.
  auto EndEdge() const noexcept {
    SelfCheck_();
    return EdgeLocalIterator<Dim>(
      Mesh_, Mesh_->FaceEdges_[FaceIndex_].End());
  }

  /// @brief Iterator pointing to \
  ///   the first face that shares an edge with the current face.
  auto BeginFace() const noexcept {
    SelfCheck_();
    return FaceLocalIterator<Dim>(
      Mesh_, Mesh_->FaceFaces_[FaceIndex_].Begin());
  }

  /// @brief Iterator pointing to the face following \
  ///   the last face that share an edge with the current face.
  auto EndFace() const noexcept {
    SelfCheck_();
    return FaceLocalIterator<Dim>(
      Mesh_, Mesh_->FaceFaces_[FaceIndex_].End());
  }

  /// @brief Get the inner and outer cells of the face.
  auto Cells() const noexcept {
    SelfCheck_();
    auto const& cellIndices = Mesh_->FaceCells_[FaceIndex_];
    return std::make_pair(
      Mesh_->Cell(cellIndices.first), Mesh_->Cell(cellIndices.second));
  }

  /// @brief Iterate through all nodes of the face.
  template<class Func>
  void ForEachNode(Func const& nodeFunc) const noexcept {
    std::for_each(BeginNode(), EndNode(), nodeFunc);
  }

  /// @brief Iterate through all edges of the face.
  template<class Func>
  void ForEachEdge(Func const& edgeFunc) const noexcept {
    std::for_each(BeginEdge(), EndEdge(), edgeFunc);
  }

  /// @brief Iterate through all faces \
  ///   that share an edge with the current face.
  template<class Func>
  void ForEachFace(Func const& faceFunc) const noexcept {
    std::for_each(BeginFace(), EndFace(), faceFunc);
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

  /// @brief Iterator pointing to \
  ///   the first node of the current cell.
  auto BeginNode() const noexcept {
    SelfCheck_();
    return NodeLocalIterator<Dim>(
      Mesh_, Mesh_->CellNodes_[CellIndex_].Begin());
  }

  /// @brief Iterator pointing to the node following \
  ///   the last node of the current cell.
  auto EndNode() const noexcept {
    SelfCheck_();
    return NodeLocalIterator<Dim>(
      Mesh_, Mesh_->CellNodes_[CellIndex_].End());
  }

  /// @brief Iterator pointing to \
  ///   the first edge of the current cell.
  auto BeginEdge() const noexcept {
    SelfCheck_();
    return EdgeLocalIterator<Dim>(
      Mesh_, Mesh_->CellEdges_[CellIndex_].Begin());
  }

  /// @brief Iterator pointing to the edge following \
  ///   the last edge of the current cell.
  auto EndEdge() const noexcept {
    SelfCheck_();
    return EdgeLocalIterator<Dim>(
      Mesh_, Mesh_->CellEdges_[CellIndex_].End());
  }

  /// @brief Iterator pointing to \
  ///   the first face of the current cell.
  auto BeginFace() const noexcept {
    SelfCheck_();
    return FaceLocalIterator<Dim>(
      Mesh_, Mesh_->CellFaces_[CellIndex_].Begin());
  }

  /// @brief Iterator pointing to the face following \
  ///   the last face of the current cell.
  auto EndFace() const noexcept {
    SelfCheck_();
    return FaceLocalIterator<Dim>(
      Mesh_, Mesh_->CellFaces_[CellIndex_].End());
  }

  /// @brief Iterator pointing to \
  ///   the first cell that shares a face with the current cell.
  auto BeginCell() const noexcept {
    SelfCheck_();
    return CellLocalIterator<Dim>(
      Mesh_, Mesh_->CellCells_[CellIndex_].Begin());
  }

  /// @brief Iterator pointing to the face following \
  ///   the last cell that share a face with the current cell.
  auto EndCell() const noexcept {
    SelfCheck_();
    return CellLocalIterator<Dim>(
      Mesh_, Mesh_->CellCells_[CellIndex_].End());
  }

  /// @brief Iterate through all nodes of the cell.
  template<class Func>
  void ForEachNode(Func const& nodeFunc) const {
    std::for_each(BeginNode(), EndNode(), nodeFunc);
  }

  /// @brief Iterate through all edges of the cell.
  template<class Func>
  void ForEachEdge(Func const& edgeFunc) const {
    std::for_each(BeginEdge(), EndEdge(), edgeFunc);
  }

  /// @brief Iterate through all faces of the face.
  template<class Func>
  void ForEachFace(Func const& faceFunc) const {
    std::for_each(BeginFace(), EndFace(), faceFunc);
  }

  /// @brief Iterate through all cells \
  ///   that share a face with the current cell.
  template<class Func>
  void ForEachCell(Func const& cellFunc) const {
    std::for_each(BeginCell(), EndCell(), cellFunc);
  }

}; // class CellView<...>

/// ----------------------------------------------------------------- ///
/// @brief Immutable element random access iterator, \
///   can be used for the both mesh elements iteration and \
//    local connectivity iterations.
/// ----------------------------------------------------------------- ///
template<class Mesh, class View, class PointerOrIndex>
class ElementIterator final : 
  public BaseIterator<ElementIterator<Mesh, View, PointerOrIndex>> {
private:
  Mesh const* Mesh_;
  PointerOrIndex PointerOrIndex_;

public:

  /// @brief Construct an iterator.
  explicit ElementIterator(Mesh const* mesh,
                           PointerOrIndex pointerOrIndex) noexcept :
      Mesh_{mesh}, PointerOrIndex_{pointerOrIndex} {
    StormMeshAssert(Mesh_ != nullptr);
  }

  /// @brief Dereference an iterator to a view.
  View operator*() const noexcept {
    if constexpr (std::is_pointer_v<PointerOrIndex>) {
      return View(Mesh_, *PointerOrIndex_);
    } else {
      return View(Mesh_, PointerOrIndex_);
    }
  }

}; // class ElementIterator<...>

} // namespace Storm
