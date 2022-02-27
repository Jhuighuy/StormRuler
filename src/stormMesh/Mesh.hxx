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

#include <ranges>
#include <tuple>
#include <utility>

#include <stormBase.hxx>
#include <stormMesh/Forward.hxx>

namespace Storm {

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

  size_t NumAllNodes_;
  size_t NumAllEdges_;
  size_t NumAllFaces_;
  size_t NumAllCells_;
  std::vector<size_t> NodeRanges_;
  std::vector<size_t> EdgeRanges_;
  std::vector<size_t> FaceRanges_;
  std::vector<size_t> CellRanges_;
  Table<void> NodeNodes_; // NodeEdges_*EdgeNodes_
  Table<void> NodeEdges_; // (EdgeNodes_)^T
  Table<void> NodeFaces_; // (FaceNodes_)^T
  Table<void> NodeCells_; // (CellNodes_)^T
  std::vector<std::pair<size_t, size_t>> EdgeNodes_; // T
  Table<void> EdgeEdges_; // EdgeNodes_*EdgeNodes_
  Table<void> EdgeFaces_; // (FaceEdges_)^T, Pass 2.
  Table<void> EdgeCells_; // (CellNodes_)^T
  Table<void> FaceNodes_; // T
  Table<void> FaceEdges_; // T, Pass 2.
  Table<void> FaceFaces_; // FaceEdges_*EdgeFaces_, Pass 1.
  std::vector<std::pair<size_t, size_t>> FaceCells_;
  Table<void> CellNodes_; // T, IN
  Table<void> CellEdges_; // T
  Table<void> CellFaces_; // T, Pass 1.
  Table<void> CellCells_; // CellFaces_*FaceCells_
  std::vector<GVec<Dim>> NodeCoords_;
  std::vector<std::tuple<GVec<Dim>, real_t>> EdgeDirsAndLengths_;
  std::vector<std::tuple<GVec<Dim>, GVec<Dim>, real_t>> 
              FaceCenterCoordsAndNormalsAndAreas_;
  std::vector<std::tuple<GVec<Dim>, real_t>> CellCenterCoordsAndVolumes_;
  std::vector<ShapeType> FaceShapes_;
  std::vector<ShapeType> CellShapes_;

private:

  /// @brief Transform the specified \
  ///   index range to the element view range.
  template<class View, class IndexRange>
  auto ViewRange_(IndexRange indexRange) const noexcept {
    return indexRange | std::views::transform(
      [this](auto index){ return View(this, index); });
  }

  /// @brief Transform the iota \
  ///   index range to the element view range.
  template<class View, class Index>
  auto IotaViewRange_(Index first, Index last) const noexcept {
    return ViewRange_<View>(std::views::iota(first, last));
  }

public:

  /// @brief Total number of nodes in the mesh.
  size_t NumAllNodes() const noexcept {
    return NumAllNodes_;
  }

  /// @brief Total number of edges in the mesh.
  size_t NumAllEdges() const noexcept {
    return NumAllEdges_;
  }

  /// @brief Total number of faces in the mesh.
  size_t NumAllFaces() const noexcept {
    return NumAllFaces_;
  }

  /// @brief Total number of cells in the mesh.
  size_t NumAllCells() const noexcept {
    return NumAllCells_;
  }

  /// @brief Get node at the index.
  auto Node(size_t nodeIndex) const noexcept {
    StormMeshAssert(nodeIndex < NumAllNodes_);
    return NodeView<Dim>(this, nodeIndex);
  }

  /// @brief Get edge at the index.
  auto Edge(size_t edgeIndex) const noexcept {
    StormMeshAssert(edgeIndex < NumAllEdges_);
    return EdgeView<Dim>(this, edgeIndex);
  }

  /// @brief Get face at the index.
  auto Face(size_t faceIndex) const noexcept {
    StormMeshAssert(faceIndex < NumAllFaces_);
    return FaceView<Dim>(this, faceIndex);
  }

  /// @brief Get cell at the index.
  auto Cell(size_t cellIndex) const noexcept {
    StormMeshAssert(cellIndex < NumAllCells_);
    return CellView<Dim>(this, cellIndex);
  }

  /// @brief Number of node marks in the mesh.
  size_t NumNodeMarks() const noexcept {
    return NodeRanges_.size() - 1;
  }

  /// @brief Number of edge marks in the mesh.
  size_t NumEdgeMarks() const noexcept {
    return EdgeRanges_.size() - 1;
  }

  /// @brief Number of face marks in the mesh.
  size_t NumFaceMarks() const noexcept {
    return FaceRanges_.size() - 1;
  }

  /// @brief Number of cell marks in the mesh.
  size_t NumCellMarks() const noexcept {
    return CellRanges_.size() - 1;
  }

  /// @brief Range of nodes with the specified mark.
  auto Nodes(size_t nodeMark = 0) const noexcept {
    StormMeshAssert(nodeMark < NumNodeMarks());
    return IotaViewRange_<NodeView<Dim>>(
      NodeRanges_[nodeMark], NodeRanges_[nodeMark + 1]);
  }

  /// @brief Range of edges with the specified mark.
  auto Edges(size_t edgeMark = 0) const noexcept {
    StormMeshAssert(edgeMark < NumEdgeMarks());
    return IotaViewRange_<EdgeView<Dim>>(
      EdgeRanges_[edgeMark], EdgeRanges_[edgeMark + 1]);
  }

  /// @brief Range of faces with the specified mark.
  auto Faces(size_t faceMark = 0) const noexcept {
    StormMeshAssert(faceMark < NumFaceMarks());
    return IotaViewRange_<FaceView<Dim>>(
      FaceRanges_[faceMark], FaceRanges_[faceMark + 1]);
  }

  /// @brief Range of faces with the specified mark.
  auto Cells(size_t cellMark = 0) const noexcept {
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
  size_t NodeIndex_;

  void SelfCheck_() const noexcept {
    StormMeshAssert(
      Mesh_ != nullptr && NodeIndex_ < Mesh_->NumAllNodes());
  }

public:

  /// @brief Construct pointer to a node.
  explicit NodeView(Mesh<Dim> const* mesh,
                    size_t nodeIndex) noexcept :
      Mesh_{mesh}, NodeIndex_{nodeIndex} {
    SelfCheck_();
  }

  /// @brief Implicit cast to index.
  operator size_t() const noexcept {
    SelfCheck_();
    return NodeIndex_;
  }

  /// @brief Coordinates of the current node.
  GVec<Dim> const& Coords() const noexcept {
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
  size_t EdgeIndex_;

  void SelfCheck_() const noexcept {
    StormMeshAssert(
      Mesh_ != nullptr && EdgeIndex_ < Mesh_->NumAllEdges());
  }

public:

  /// @brief Construct pointer to an edge.
  explicit EdgeView(Mesh<Dim> const* mesh,
                    size_t edgeIndex) noexcept :
      Mesh_{mesh}, EdgeIndex_{edgeIndex} {
    SelfCheck_();
  }

  /// @brief Implicit cast to index.
  operator size_t() const noexcept {
    SelfCheck_();
    return EdgeIndex_;
  }

  /// @brief Unit direction to the current edge.
  GVec<Dim> const& Direction() const noexcept {
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
  size_t FaceIndex_;

  void SelfCheck_() const noexcept {
    StormMeshAssert(
      Mesh_ != nullptr && FaceIndex_ < Mesh_->NumAllFaces());
  }

public:

  /// @brief Construct pointer to a face.
  explicit FaceView(Mesh<Dim> const* mesh,
                    size_t faceIndex) noexcept :
      Mesh_{mesh}, FaceIndex_{faceIndex} {
    SelfCheck_();
  }

  /// @brief Implicit cast to index.
  operator size_t() const noexcept {
    SelfCheck_();
    return FaceIndex_;
  }

  /// @brief Center coordinates of the current face.
  GVec<Dim> const& CenterCoords() const noexcept {
    SelfCheck_();
    return std::get<0>(
      Mesh_->FaceCenterCoordsAndNormalsAndAreas_[FaceIndex_]);
  }

  /// @brief Unit normal to the current face.
  GVec<Dim> const& Normal() const noexcept {
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
  size_t CellIndex_;

  void SelfCheck_() const noexcept {
    StormMeshAssert(
      Mesh_ != nullptr && CellIndex_ < Mesh_->NumAllCells());
  }

public:

  /// @brief Construct pointer to a cell.
  explicit CellView(Mesh<Dim> const* mesh,
                    size_t cellIndex) noexcept :
      Mesh_{mesh}, CellIndex_{cellIndex} {
    SelfCheck_();
  }

  /// @brief Implicit cast to index.
  operator size_t() const noexcept {
    SelfCheck_();
    return CellIndex_;
  }

  /// @brief Center coordinates of the current cell.
  GVec<Dim> const& CenterCoords() const noexcept {
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
