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

#include <Storm/Utils/Meta.hpp>

#include <Storm/Mallard/Fwd.hpp>
#include <Storm/Mallard/Shape.hpp>

#include <compare>
#include <concepts>
#include <ranges>

namespace Storm {

struct LabelTag;

/// @brief Label index type.
using Label = Index<LabelTag>;

template<mesh Mesh>
class NodeView;
template<mesh Mesh>
class EdgeView;
template<mesh Mesh>
class FaceView;
template<mesh Mesh>
class CellView;

template<class Mesh>
NodeView(const Mesh&, NodeIndex) -> NodeView<Mesh>;
template<class Mesh>
EdgeView(const Mesh&, EdgeIndex) -> EdgeView<Mesh>;
template<class Mesh>
FaceView(const Mesh&, FaceIndex<Mesh>) -> FaceView<Mesh>;
template<class Mesh>
CellView(const Mesh&, CellIndex<Mesh>) -> CellView<Mesh>;

namespace detail_ {
  constexpr auto to_node_view_(const auto& mesh) noexcept {
    return std::views::transform(
        [&](NodeIndex node_index) { return NodeView(mesh, node_index); });
  }
  constexpr auto to_edge_view_(const auto& mesh) noexcept {
    return std::views::transform(
        [&](EdgeIndex edge_index) { return EdgeView(mesh, edge_index); });
  }
  template<class Mesh>
  constexpr auto to_face_view_(const Mesh& mesh) noexcept {
    return std::views::transform(
        [&](FaceIndex<Mesh> face_index) { return FaceView(mesh, face_index); });
  }
  template<class Mesh>
  constexpr auto to_cell_view_(const Mesh& mesh) noexcept {
    return std::views::transform(
        [&](CellIndex<Mesh> cell_index) { return CellView(mesh, cell_index); });
  }
} // namespace detail_

/// @brief Mesh entity view.
template<mesh Mesh, index Index>
class EntityView;

template<mesh Mesh, size_t I>
  requires std::is_object_v<Mesh>
class EntityView<Mesh, EntityIndex<I>> {
private:

  template<mesh, index>
  friend class EntityView;

  const Mesh* p_mesh_;
  EntityIndex<I> index_;

protected:

  /// @brief Construct an entity view with a @p mesh and @p index.
  constexpr EntityView(const Mesh& mesh, EntityIndex<I> index) noexcept
      : p_mesh_{&mesh}, index_{index} {}

  /// @brief Destroy the entity view.
  constexpr ~EntityView() = default;

public:

  /// @brief Entity mesh.
  /// @{
  [[nodiscard]] constexpr Mesh& mesh() noexcept {
    return *p_mesh_;
  }
  [[nodiscard]] constexpr const Mesh& mesh() const noexcept {
    return *p_mesh_;
  }
  /// @}

  /// @brief Entity index.
  /// @{
  [[nodiscard]] constexpr EntityIndex<I> index() const noexcept {
    return index_;
  }
  [[nodiscard]] constexpr size_t index_sz() const noexcept {
    return static_cast<size_t>(index_);
  }
  /// @}

  /// @brief Cast to entity index operator.
  [[nodiscard]] constexpr operator EntityIndex<I>() const noexcept {
    return index_;
  }

  /// @brief Comparison operator.
  [[nodiscard]] constexpr auto
  operator<=>(const EntityView& other) const noexcept {
    STORM_ASSERT_(p_mesh_ == other.p_mesh_,
                  "Entities correspond to the different meshes!");
    return index_ <=> other.index_;
  }

  /// @brief Entity label.
  [[nodiscard]] constexpr auto label() const noexcept {
    return mesh().label(index_);
  }

  /// @brief Range of the adjacent nodes.
  [[nodiscard]] constexpr auto nodes() const noexcept {
    return mesh().adjacent(index_, meta::type_v<NodeIndex>) |
           detail_::to_node_view_(*p_mesh_);
  }

  /// @brief Range of the adjacent edges.
  [[nodiscard]] constexpr auto edges() const noexcept {
    return mesh().adjacent(index_, meta::type_v<EdgeIndex>) |
           detail_::to_edge_view_(*p_mesh_);
  }

  /// @brief Range of the adjacent faces.
  [[nodiscard]] constexpr auto faces() const noexcept {
    return mesh().adjacent(index_, meta::type_v<FaceIndex<Mesh>>) |
           detail_::to_face_view_(*p_mesh_);
  }

  /// @brief Range of the adjacent cells.
  [[nodiscard]] constexpr auto cells() const noexcept {
    return mesh().adjacent(index_, meta::type_v<CellIndex<Mesh>>) |
           detail_::to_cell_view_(*p_mesh_);
  }

  /// @brief Iterate all the adjacent nodes.
  template<std::invocable<NodeView<Mesh>> Func>
  constexpr void for_each_node(Func func) const noexcept {
    std::ranges::for_each(nodes(), func);
  }

  /// @brief Iterate all the adjacent edges.
  template<std::invocable<EdgeView<Mesh>> Func>
  constexpr void for_each_edge(Func func) const noexcept {
    std::ranges::for_each(edges(), func);
  }

  /// @brief Iterate all the adjacent faces.
  /// @{
  template<std::invocable<FaceView<Mesh>> Func>
  constexpr void for_each_face(Func func) const noexcept {
    std::ranges::for_each(faces(), func);
  }
  template<std::invocable<CellView<Mesh>, CellView<Mesh>> Func>
  constexpr void for_each_face_cells(Func func) const noexcept {
    std::ranges::for_each(faces(), [&](FaceView<Mesh> face) {
      func(face.inner_cell(), face.outer_cell());
    });
  }
  /// @}

  /// @brief Iterate all the adjacent cells.
  template<std::invocable<CellView<Mesh>> Func>
  constexpr void for_each_cell(Func func) const noexcept {
    std::ranges::for_each(cells(), func);
  }

}; // class EntityView

/// @brief Mesh node view.
template<mesh Mesh>
class NodeView final : public EntityView<Mesh, NodeIndex> {
public:

  /// @brief Construct a node view.
  constexpr NodeView(const Mesh& mesh, NodeIndex index) noexcept
      : EntityView<Mesh, NodeIndex>(mesh, index) {}

  /// @brief Node coordinates.
  [[nodiscard]] constexpr auto position() const noexcept {
    return this->mesh().position(this->index());
  }

}; // class NodeView

/// @brief Mesh edge view.
template<mesh Mesh>
class EdgeView final : public EntityView<Mesh, EdgeIndex> {
public:

  /// @brief Construct an edge view.
  constexpr EdgeView(const Mesh& mesh, EdgeIndex index) noexcept
      : EntityView<Mesh, EdgeIndex>(mesh, index) {}

  /// @brief Edge shape type.
  [[nodiscard]] constexpr shapes::Type shape_type() const noexcept {
    return this->mesh().shape_type(this->index());
  }

  /// @brief Edge length.
  [[nodiscard]] constexpr real_t length() const noexcept {
    return this->mesh().volume(this->index());
  }

  /// @brief Edge barycenter position.
  [[nodiscard]] constexpr auto center() const noexcept {
    return this->mesh().position(this->index());
  }

}; // class EdgeView

/// @brief Face view.
template<mesh Mesh>
class FaceView final : public EntityView<Mesh, FaceIndex<Mesh>> {
public:

  /// @brief Construct a face view.
  constexpr FaceView(const Mesh& mesh, FaceIndex<Mesh> index) noexcept
      : EntityView<Mesh, FaceIndex<Mesh>>(mesh, index) {}

  /// @brief Face shape type.
  [[nodiscard]] constexpr shapes::Type shape_type() const noexcept {
    return this->mesh().shape_type(this->index());
  }

  /// @brief Face area (or length in 2D).
  [[nodiscard]] constexpr real_t area() const noexcept {
    return this->mesh().volume(this->index());
  }

  /// @brief Face barycenter position.
  [[nodiscard]] constexpr auto center() const noexcept {
    return this->mesh().position(this->index());
  }

  /// @brief Face normal.
  [[nodiscard]] constexpr auto normal() const noexcept {
    return this->mesh().normal(this->index());
  }

  /// @brief Get the adjacent inner cell.
  [[nodiscard]] constexpr CellView<Mesh> inner_cell() const noexcept {
    STORM_ASSERT_(!this->cells().empty(),
                  "The face does not have an adjacent inner cell!");
    return this->cells().front();
  }

  /// @brief Get the adjacent outer cell.
  [[nodiscard]] constexpr CellView<Mesh> outer_cell() const noexcept {
    STORM_ASSERT_(this->cells().size() == 2,
                  "The face does not have an adjacent outer cell!");
    return this->cells().back();
  }

}; // class FaceView

/// @brief Mesh cell view.
template<mesh Mesh>
class CellView final : public EntityView<Mesh, CellIndex<Mesh>> {
public:

  /// @brief Construct a cell view.
  constexpr CellView(const Mesh& mesh, CellIndex<Mesh> index) noexcept
      : EntityView<Mesh, CellIndex<Mesh>>(mesh, index) {}

  /// @brief Cell shape type.
  [[nodiscard]] constexpr shapes::Type shape_type() const noexcept {
    return this->mesh().shape_type(this->index());
  }

  /// @brief Cell volume (or area in 2D).
  [[nodiscard]] constexpr real_t volume() const noexcept {
    return this->mesh().volume(this->index());
  }

  /// @brief Cell barycenter position.
  [[nodiscard]] constexpr auto center() const noexcept {
    return this->mesh().position(this->index());
  }

  /// @todo REMOVE ME!
  template<std::invocable<CellView<Mesh>> Func>
  constexpr void for_each_cell(Func func) const noexcept {
    this->for_each_face_cells(
        [&](CellView<Mesh> cell_inner, CellView<Mesh> cell_outer) {
          if (cell_inner.index() == this->index()) func(cell_outer);
          if (cell_outer.index() == this->index()) func(cell_inner);
        });
  }

}; // class CellView

/// @brief CRTP interface to a mesh.
template<crtp_derived Derived>
class MeshInterface {
private:

  constexpr Derived& self_() noexcept {
    static_assert(std::derived_from<Derived, MeshInterface>);
    return static_cast<Derived&>(*this);
  }
  constexpr const Derived& self_() const noexcept {
    return const_cast<MeshInterface&>(*this).self_();
  }

public:

  /// @brief Number of the node labels.
  [[nodiscard]] constexpr size_t num_node_labels() const noexcept {
    return self_().num_labels(meta::type_v<NodeIndex>);
  }
  /// @brief Number of the edge labels.
  [[nodiscard]] constexpr size_t num_edge_labels() const noexcept {
    return self_().num_labels(meta::type_v<EdgeIndex>);
  }
  /// @brief Number of the face labels.
  [[nodiscard]] constexpr size_t num_face_labels() const noexcept {
    return self_().num_labels(meta::type_v<FaceIndex<Derived>>);
  }
  /// @brief Number of the cell labels.
  [[nodiscard]] constexpr size_t num_cell_labels() const noexcept {
    return self_().num_labels(meta::type_v<CellIndex<Derived>>);
  }

  /// @brief Number of the nodes.
  [[nodiscard]] constexpr size_t num_nodes() const noexcept {
    return self_().num_entities(meta::type_v<NodeIndex>);
  }
  /// @brief Number of the nodes with a @p label.
  [[nodiscard]] constexpr size_t num_nodes(Label label) const noexcept {
    return self_().num_entities(label, meta::type_v<NodeIndex>);
  }

  /// @brief Number of the edges.
  [[nodiscard]] constexpr size_t num_edges() const noexcept {
    return self_().num_entities(meta::type_v<EdgeIndex>);
  }
  /// @brief Number of the edges with a @p label.
  [[nodiscard]] constexpr size_t num_edges(Label label) const noexcept {
    return self_().num_entities(label, meta::type_v<EdgeIndex>);
  }

  /// @brief Number of the faces.
  [[nodiscard]] constexpr size_t num_faces() const noexcept {
    return self_().num_entities(meta::type_v<FaceIndex<Derived>>);
  }
  /// @brief Number of the faces with a @p label.
  [[nodiscard]] constexpr size_t num_faces(Label label) const noexcept {
    return self_().num_entities(label, meta::type_v<FaceIndex<Derived>>);
  }

  /// @brief Number of the cells.
  [[nodiscard]] constexpr size_t num_cells() const noexcept {
    return self_().num_entities(meta::type_v<CellIndex<Derived>>);
  }
  /// @brief Number of the cells with a @p label.
  [[nodiscard]] constexpr size_t num_cells(Label label) const noexcept {
    return self_().num_entities(label, meta::type_v<CellIndex<Derived>>);
  }

  /// @brief Range of the nodes.
  [[nodiscard]] constexpr auto nodes() const noexcept {
    return self_().entities(meta::type_v<NodeIndex>) |
           detail_::to_node_view_(self_());
  }
  /// @brief Range of the nodes with a @p label.
  [[nodiscard]] constexpr auto nodes(Label label) const noexcept {
    return self_().entities(label, meta::type_v<NodeIndex>) |
           detail_::to_node_view_(self_());
  }

  /// @brief Range of the edges.
  [[nodiscard]] constexpr auto edges() const noexcept {
    return self_().entities(meta::type_v<EdgeIndex>) |
           detail_::to_edge_view_(self_());
  }
  /// @brief Range of the edges with @p label.
  [[nodiscard]] constexpr auto edges(Label label) const noexcept {
    return self_().entities(label, meta::type_v<EdgeIndex>) |
           detail_::to_edge_view_(self_());
  }

  /// @brief Range of the faces.
  [[nodiscard]] constexpr auto faces() const noexcept {
    return self_().entities(meta::type_v<FaceIndex<Derived>>) |
           detail_::to_face_view_(self_());
  }
  /// @brief Range of the faces with @p label.
  [[nodiscard]] constexpr auto faces(Label label) const noexcept {
    return self_().entities(label, meta::type_v<FaceIndex<Derived>>) |
           detail_::to_face_view_(self_());
  }

  /// @brief Range of the cells.
  [[nodiscard]] constexpr auto cells() const noexcept {
    return self_().entities(meta::type_v<CellIndex<Derived>>) |
           detail_::to_cell_view_(self_());
  }
  /// @brief Range of the cells with @p label.
  [[nodiscard]] constexpr auto cells(Label label) const noexcept {
    return self_().entities(label, meta::type_v<CellIndex<Derived>>) |
           detail_::to_cell_view_(self_());
  }

  /// @brief Range of the interior nodes.
  [[nodiscard]] constexpr auto interior_nodes() const noexcept {
    return nodes(Label{0});
  }

  /// @brief Range of the interior edges.
  [[nodiscard]] constexpr auto interior_edges() const noexcept {
    return edges(Label{0});
  }

  /// @brief Range of the interior faces.
  [[nodiscard]] constexpr auto interior_faces() const noexcept {
    return faces(Label{0});
  }

  /// @brief Range of the interior cells.
  [[nodiscard]] constexpr auto interior_cells() const noexcept {
    return cells(Label{0});
  }

  /// @brief Range of the boundary nodes.
  [[nodiscard]] constexpr auto boundary_nodes() const noexcept {
    return nodes() | std::views::drop(num_nodes(Label{0}));
  }

  /// @brief Range of the boundary edges.
  [[nodiscard]] constexpr auto boundary_edges() const noexcept {
    return edges() | std::views::drop(num_edges(Label{0}));
  }

  /// @brief Range of the boundary faces.
  [[nodiscard]] constexpr auto boundary_faces() const noexcept {
    return faces() | std::views::drop(num_faces(Label{0}));
  }

  /// @brief Range of the boundary cells.
  [[nodiscard]] constexpr auto boundary_cells() const noexcept {
    return cells() | std::views::drop(num_cells(Label{0}));
  }

}; // class MeshInterface

} // namespace Storm
