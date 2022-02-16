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

#include <vector>

#include <stormBase.hxx>

namespace Storm {

template<class Tag, class Index = size_t>
using TaggedInteger = Index;
//class TaggedInteger;

template<class Index = size_t>
class CsrTable {
public:

  template<class Func>
  void ForEachRowColumns(Index index, 
                         Func const& func) const noexcept;

}; // class CsrTable<...>

class NodeTag;
class EdgeTag;
class FaceTag;
class CellTag;

using NodeIndex = TaggedInteger<NodeTag>;
using EdgeIndex = TaggedInteger<EdgeTag>;
using FaceIndex = TaggedInteger<FaceTag>;
using CellIndex = TaggedInteger<CellTag>;

template<size_t Dim>
class Vec;

template<size_t Dim>
class NodePtr;
template<size_t Dim>
class EdgePtr;
template<size_t Dim>
class FacePtr;
template<size_t Dim>
class CellPtr;

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Hybrid unstructured multidimensional mesh.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<size_t Dim>
class Mesh : public BaseObject {
private:

  template<size_t>
  friend class NodePtr;
  template<size_t>
  friend class EdgePtr;
  template<size_t>
  friend class FacePtr;
  template<size_t>
  friend class CellPtr;

  size_t NumAllNodes_, NumAllEdges_;
  size_t NumAllFaces_, NumAllCells_;
  std::vector<size_t> NodeMarkRanges_;
  std::vector<size_t> EdgeMarkRanges_;
  std::vector<size_t> FaceMarkRanges_;
  std::vector<size_t> CellMarkRanges_;

  CsrTable<> NodeNodes_, NodeEdges_, NodeFaces_, NodeCells_;
  CsrTable<> EdgeEdges_, EdgeFaces_, EdgeCells_;
  CsrTable<> FaceNodes_, FaceEdges_, FaceFaces_;
  CsrTable<> CellNodes_, CellEdges_, CellFaces_, CellCells_;
  std::vector<std::pair<size_t, size_t>> EdgeNodes_;
  std::vector<std::pair<size_t, size_t>> FaceCells_;

public:

  /// @brief Number of nodes in the mesh.
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
    StormAssert(nodeIndex < NumAllNodes_);
    return NodePtr<Dim>(this, nodeIndex);
  }

  /// @brief Get edge at the index.
  auto Edge(EdgeIndex edgeIndex) const noexcept {
    StormAssert(edgeIndex < NumAllEdges_);
    return EdgePtr<Dim>(this, edgeIndex);
  }

  /// @brief Get face at the index.
  auto Face(FaceIndex faceIndex) const noexcept {
    StormAssert(faceIndex < NumAllFaces_);
    return FacePtr<Dim>(this, faceIndex);
  }

  /// @brief Get cell at the index.
  auto Cell(CellIndex cellIndex) const noexcept {
    StormAssert(cellIndex < NumAllCells_);
    return CellPtr<Dim>(this, cellIndex);
  }

  /// @brief Number of node marks in the mesh.
  size_t NumNodeMarks() const noexcept {
    return NodeMarkRanges_.size() - 1;
  }

  /// @brief Number of edge marks in the mesh.
  size_t NumEdgeMarks() const noexcept {
    return EdgeMarkRanges_.size() - 1;
  }

  /// @brief Number of face marks in the mesh.
  size_t NumFaceMarks() const noexcept {
    return FaceMarkRanges_.size() - 1;
  }

  /// @brief Number of cell marks in the mesh.
  size_t NumCellMarks() const noexcept {
    return CellMarkRanges_.size() - 1;
  }

  /// @brief Pointer to the first node with a given mark.
  auto BeginNode(size_t mark = 0) const noexcept {
    StormAssert(mark < NumNodeMarks());
    return Node(NodeMarkRanges_[mark]);
  }

  /// @brief Pointer to the first edge with a given mark.
  auto BeginEdge(size_t mark = 0) const noexcept {
    StormAssert(mark < NumEdgeMarks());
    return Edge(EdgeMarkRanges_[mark]);
  }

  /// @brief Pointer to the first face with a given mark.
  auto BeginFace(size_t mark = 0) const noexcept {
    StormAssert(mark < NumFaceMarks());
    return Face(FaceMarkRanges_[mark]);
  }

  /// @brief Pointer to the first cell with a given mark.
  auto BeginCell(size_t mark = 0) const noexcept {
    StormAssert(mark < NumCellMarks());
    return Cell(CellMarkRanges_[mark]);
  }

  /// @brief Pointer to a node after the last node with a given mark.
  auto EndNode(size_t mark = 0) const noexcept {
    StormAssert(mark < NumNodeMarks());
    return Node(NodeMarkRanges_[mark + 1]);
  }

  /// @brief Pointer to a edge after the last edge with a given mark.
  auto EndEdge(size_t mark = 0) const noexcept {
    StormAssert(mark < NumEdgeMarks());
    return Edge(EdgeMarkRanges_[mark + 1]);
  }

  /// @brief Pointer to a face after the last face with a given mark.
  auto EndFace(size_t mark = 0) const noexcept {
    StormAssert(mark < NumFaceMarks());
    return Face(FaceMarkRanges_[mark + 1]);
  }

  /// @brief Pointer to a face after the last face with a given mark.
  auto EndCell(size_t mark = 0) const noexcept {
    StormAssert(mark < NumCellMarks());
    return Cell(CellMarkRanges_[mark + 1]);
  }

}; // class Mesh<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Pointer to a node of the unstructured mesh.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<size_t Dim>
class NodePtr {
private:
  Mesh<Dim> const* Mesh_;
  NodeIndex NodeIndex_;

public:

  /// @brief Construct pointer to a node. 
  explicit NodePtr(Mesh<Dim> const* mesh, 
                   NodeIndex nodeIndex) noexcept :
    Mesh_{mesh}, NodeIndex_{nodeIndex} {
  }

  /// @brief Iterate through all nodes \
  ///   that share an edge with the current node.
  template<class Func>
  void ForEachNode(Func const& nodeFunc) const noexcept {
    Mesh_->NodeNodes_.ForEachRowColumns(NodeIndex_,
      [&](size_t nodeIndex) { nodeFunc( Mesh_->Node(nodeIndex) ); });
  }

  /// @brief Iterate through all edges \
  ///   that ...
  template<class Func>
  void ForEachEdge(Func const& edgeFunc) const noexcept {
    Mesh_->NodeEdges_.ForEachRowColumns(NodeIndex_,
      [&](size_t edgeIndex) { edgeFunc( Mesh_->Edge(edgeIndex) ); });
  }

  /// @brief Iterate through all faces \
  ///   that ...
  template<class Func>
  void ForEachFace(Func const& faceFunc) const noexcept {
    Mesh_->NodeFaces_.ForEachRowColumns(NodeIndex_,
      [&](size_t faceIndex) { faceFunc( Mesh_->Face(faceIndex) ); });
  }

  /// @brief Iterate through all cells \
  ///   that ...
  template<class Func>
  void ForEachCell(Func const& cellFunc) const noexcept {
    Mesh_->NodeCells_.ForEachRowColumns(NodeIndex_,
      [&](size_t cellIndex) { cellFunc( Mesh_->Cell(cellIndex) ); });
  }

}; // class NodePtr<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Pointer to an edge of the unstructured mesh.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<size_t Dim>
class EdgePtr final {
private:
  Mesh<Dim> const* Mesh_;
  EdgeIndex EdgeIndex_;

public:

  /// @brief Construct pointer to an edge. 
  explicit EdgePtr(Mesh<Dim> const* mesh, 
                   EdgeIndex edgeIndex) noexcept :
    Mesh_{mesh}, EdgeIndex_{edgeIndex} {
  }

  /// @brief Get the nodes of the edge.
  auto Nodes() const noexcept {
    return std::make_pair(
      Mesh_->Node(Mesh_->EdgeNodes_[EdgeIndex_].first),
      Mesh_->Node(Mesh_->EdgeNodes_[EdgeIndex_].second));
  }

  /// @brief Iterate through all edges \
  ///   that share a node with the current edge.
  template<class Func>
  void ForEachEdge(Func const& edgeFunc) const noexcept {
    Mesh_->EdgeEdges_.ForEachRowColumns(EdgeIndex_,
      [&](size_t edgeIndex) { edgeFunc( Mesh_->Edge(edgeIndex) ); });
  }

  /// @brief Iterate through all faces \
  ///   that ...
  template<class Func>
  void ForEachFace(Func const& faceFunc) const noexcept {
    Mesh_->EdgeFaces_.ForEachRowColumns(EdgeIndex_,
      [&](size_t faceIndex) { faceFunc( Mesh_->Face(faceIndex) ); });
  }

  /// @brief Iterate through all cells \
  ///   that ...
  template<class Func>
  void ForEachCell(Func const& cellFunc) const noexcept {
    Mesh_->EdgeCells_.ForEachRowColumns(EdgeIndex_,
      [&](size_t cellIndex) { cellFunc( Mesh_->Cell(cellIndex) ); });
  }

}; // class EdgePtr<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Pointer to a face of the unstructured mesh.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<size_t Dim>
class FacePtr final {
private:
  Mesh<Dim> const* Mesh_;
  FaceIndex FaceIndex_;

public:

  /// @brief Construct pointer to a face. 
  explicit FacePtr(Mesh<Dim> const* mesh, 
                   FaceIndex faceIndex) noexcept :
    Mesh_{mesh}, FaceIndex_{faceIndex} {
  }

  /// @brief Iterate through all nodes of the face.
  template<class Func>
  void ForEachNode(Func const& nodeFunc) const noexcept {
    Mesh_->FaceNodes_.ForEachRowColumns(FaceIndex_,
      [&](size_t nodeIndex) { nodeFunc( Mesh_->Node(nodeIndex) ); });
  }

  /// @brief Iterate through all edges of the face.
  template<class Func>
  void ForEachEdge(Func const& edgeFunc) const noexcept {
    Mesh_->FaceEdges_.ForEachRowColumns(FaceIndex_,
      [&](size_t edgeIndex) { edgeFunc( Mesh_->Edge(edgeIndex) ); });
  }

  /// @brief Iterate through all faces \
  ///   that ...
  template<class Func>
  void ForEachFace(Func const& faceFunc) const noexcept {
    Mesh_->FaceFaces_.ForEachRowColumns(FaceIndex_,
      [&](size_t faceIndex) { faceFunc( Mesh_->Face(faceIndex) ); });
  }

  /// @brief Get the inner and outer cells of the face.
  auto Cells() const noexcept {
    return std::make_pair(
      Mesh_->Cell(Mesh_->FaceCells_[FaceIndex_].first),
      Mesh_->Cell(Mesh_->FaceCells_[FaceIndex_].second));
  }

}; // class FacePtr<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Pointer to a cell of the unstructured mesh.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<size_t Dim>
class CellPtr final {
private:
  Mesh<Dim> const* Mesh_;
  CellIndex CellIndex_;

public:

  /// @brief Construct pointer to a cell. 
  explicit CellPtr(Mesh<Dim> const* mesh, 
                   CellIndex cellIndex) noexcept :
    Mesh_{mesh}, CellIndex_{cellIndex} {
  }

  /// @brief Iterate through all nodes of the cell.
  template<class Func>
  void ForEachNode(Func const& nodeFunc) const {
    Mesh_->CellNodes_.ForEachRowColumns(CellIndex_,
      [&](size_t nodeIndex) { nodeFunc( Mesh_->Node(nodeIndex) ); });
  }

  /// @brief Iterate through all edges of the cell.
  template<class Func>
  void ForEachEdge(Func const& edgeFunc) const {
    Mesh_->CellEdges_.ForEachRowColumns(CellIndex_,
      [&](size_t edgeIndex) { edgeFunc( Mesh_->Edge(edgeIndex) ); });
  }

  /// @brief Iterate through all faces of the face.
  template<class Func>
  void ForEachFace(Func const& faceFunc) const {
    Mesh_->CellFaces_.ForEachRowColumns(CellIndex_,
      [&](size_t faceIndex) { faceFunc( Mesh_->Face(faceIndex) ); });
  }

  /// @brief Iterate through all cells \
  ///   that share a face with the current cell.
  template<class Func>
  void ForEachCell(Func const& cellFunc) const {
    Mesh_->CellCells_.ForEachRowColumns(CellIndex_,
      [&](size_t cellIndex) { cellFunc( Mesh_->Cell(cellIndex) ); });
  }

}; // class CellPtr<...>

} // namespace Storm
