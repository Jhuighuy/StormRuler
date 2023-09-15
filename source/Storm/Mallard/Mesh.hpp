// Copyright (C) 2020-2023 Oleg Butakov
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <Storm/Base.hpp>

#include <Storm/Utils/Meta.hpp>

#include <Storm/Mallard/Fwd.hpp>
#include <Storm/Mallard/Shape.hpp>

#include <compare>
#include <concepts>
#include <ranges>

namespace Storm {

// -----------------------------------------------------------------------------

template<mesh Mesh>
class NodeView;
template<mesh Mesh>
class EdgeView;
template<mesh Mesh>
class FaceView;
template<mesh Mesh>
class CellView;

namespace _detail {
  constexpr auto _to_node_view(const auto& mesh) noexcept {
    return std::views::transform(
        [&](NodeIndex node_index) { return NodeView(mesh, node_index); });
  }
  constexpr auto _to_edge_view(const auto& mesh) noexcept {
    return std::views::transform(
        [&](EdgeIndex edge_index) { return EdgeView(mesh, edge_index); });
  }
  template<class Mesh>
  constexpr auto _to_face_view(const Mesh& mesh) noexcept {
    return std::views::transform(
        [&](FaceIndex<Mesh> face_index) { return FaceView(mesh, face_index); });
  }
  template<class Mesh>
  constexpr auto _to_cell_view(const Mesh& mesh) noexcept {
    return std::views::transform(
        [&](CellIndex<Mesh> cell_index) { return CellView(mesh, cell_index); });
  }
} // namespace _detail

/// @brief Mesh entity view.
template<mesh Mesh, index Index>
class EntityView {
private:

  template<mesh, index>
  friend class EntityView;

  const Mesh* _p_mesh;
  Index _index;

protected:

  /// @brief Construct an entity view with a @p mesh and @p index.
  constexpr EntityView(const Mesh& mesh, Index index) noexcept
      : _p_mesh{&mesh}, _index{index} {}

  /// @brief Destroy the entity view.
  constexpr ~EntityView() = default;

public:

  /// @brief Entity mesh.
  /// @{
  constexpr Mesh& mesh() noexcept {
    return *_p_mesh;
  }
  constexpr const Mesh& mesh() const noexcept {
    return *_p_mesh;
  }
  /// @}

  /// @brief Entity index.
  /// @{
  constexpr Index index() const noexcept {
    return _index;
  }
  constexpr size_t index_sz() const noexcept {
    return static_cast<size_t>(_index);
  }
  /// @}

  /// @brief Cast to entity index operator.
  constexpr operator Index() const noexcept {
    return _index;
  }

  /// @brief Comparison operator.
  constexpr auto operator<=>(const EntityView& other) const noexcept {
    STORM_ASSERT(_p_mesh == other._p_mesh,
                 "Entities correspond to the different meshes!");
    return _index <=> other._index;
  }

  /// @brief Entity label.
  constexpr auto label() const noexcept {
    return mesh().label(_index);
  }

  /// @brief Range of the adjacent nodes.
  constexpr auto nodes() const noexcept {
    return mesh().adjacent(_index, meta::type_v<NodeIndex>) |
           _detail::_to_node_view(*_p_mesh);
  }

  /// @brief Range of the adjacent edges.
  constexpr auto edges() const noexcept {
    return mesh().adjacent(_index, meta::type_v<EdgeIndex>) |
           _detail::_to_edge_view(*_p_mesh);
  }

  /// @brief Range of the adjacent faces.
  constexpr auto faces() const noexcept {
    return mesh().adjacent(_index, meta::type_v<FaceIndex<Mesh>>) |
           _detail::_to_face_view(*_p_mesh);
  }

  /// @brief Range of the adjacent cells.
  constexpr auto cells() const noexcept {
    return mesh().adjacent(_index, meta::type_v<CellIndex<Mesh>>) |
           _detail::_to_cell_view(*_p_mesh);
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

// -----------------------------------------------------------------------------

/// @brief Mesh node view.
template<mesh Mesh>
class NodeView final : public EntityView<Mesh, NodeIndex> {
public:

  /// @brief Construct a node view.
  constexpr NodeView(const Mesh& mesh, NodeIndex index) noexcept
      : EntityView<Mesh, NodeIndex>(mesh, index) {}

  /// @brief Node coordinates.
  constexpr auto position() const noexcept {
    return this->mesh().position(this->index());
  }

}; // class NodeView

template<class Mesh>
NodeView(const Mesh&, NodeIndex) -> NodeView<Mesh>;

// -----------------------------------------------------------------------------

/// @brief Mesh edge view.
template<mesh Mesh>
class EdgeView final : public EntityView<Mesh, EdgeIndex> {
public:

  /// @brief Construct an edge view.
  constexpr EdgeView(const Mesh& mesh, EdgeIndex index) noexcept
      : EntityView<Mesh, EdgeIndex>(mesh, index) {}

  /// @brief Edge shape type.
  constexpr shapes::Type shape_type() const noexcept {
    return this->mesh().shape_type(this->index());
  }

  /// @brief Edge length.
  constexpr real_t length() const noexcept {
    return this->mesh().volume(this->index());
  }

  /// @brief Edge barycenter position.
  constexpr auto center() const noexcept {
    return this->mesh().position(this->index());
  }

}; // class EdgeView

template<class Mesh>
EdgeView(const Mesh&, EdgeIndex) -> EdgeView<Mesh>;

// -----------------------------------------------------------------------------

/// @brief Face view.
template<mesh Mesh>
class FaceView final : public EntityView<Mesh, FaceIndex<Mesh>> {
public:

  /// @brief Construct a face view.
  constexpr FaceView(const Mesh& mesh, FaceIndex<Mesh> index) noexcept
      : EntityView<Mesh, FaceIndex<Mesh>>(mesh, index) {}

  /// @brief Face shape type.
  constexpr shapes::Type shape_type() const noexcept {
    return this->mesh().shape_type(this->index());
  }

  /// @brief Face area (or length in 2D).
  constexpr real_t area() const noexcept {
    return this->mesh().volume(this->index());
  }

  /// @brief Face barycenter position.
  constexpr auto center() const noexcept {
    return this->mesh().position(this->index());
  }

  /// @brief Face normal.
  constexpr auto normal() const noexcept {
    return this->mesh().normal(this->index());
  }

  /// @brief Get the adjacent inner cell.
  constexpr CellView<Mesh> inner_cell() const noexcept {
    STORM_ASSERT(!this->cells().empty(),
                 "The face does not have an adjacent inner cell!");
    return this->cells().front();
  }

  /// @brief Get the adjacent outer cell.
  constexpr CellView<Mesh> outer_cell() const noexcept {
    STORM_ASSERT(this->cells().size() == 2,
                 "The face does not have an adjacent outer cell!");
    return this->cells().back();
  }

}; // class FaceView

template<class Mesh>
FaceView(const Mesh&, FaceIndex<Mesh>) -> FaceView<Mesh>;

// -----------------------------------------------------------------------------

/// @brief Mesh cell view.
template<mesh Mesh>
class CellView final : public EntityView<Mesh, CellIndex<Mesh>> {
public:

  /// @brief Construct a cell view.
  constexpr CellView(const Mesh& mesh, CellIndex<Mesh> index) noexcept
      : EntityView<Mesh, CellIndex<Mesh>>(mesh, index) {}

  /// @brief Cell shape type.
  constexpr shapes::Type shape_type() const noexcept {
    return this->mesh().shape_type(this->index());
  }

  /// @brief Cell volume (or area in 2D).
  constexpr real_t volume() const noexcept {
    return this->mesh().volume(this->index());
  }

  /// @brief Cell barycenter position.
  constexpr auto center() const noexcept {
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

template<class Mesh>
CellView(const Mesh&, CellIndex<Mesh>) -> CellView<Mesh>;

// -----------------------------------------------------------------------------

/// @brief CRTP interface to a mesh.
template<crtp_derived Derived>
class MeshInterface {
private:

  constexpr Derived& _self() noexcept {
    static_assert(std::derived_from<Derived, MeshInterface>);
    return static_cast<Derived&>(*this);
  }
  constexpr const Derived& _self() const noexcept {
    return const_cast<MeshInterface&>(*this)._self();
  }

public:

  /// @brief Number of the node labels.
  constexpr size_t num_node_labels() const noexcept {
    return _self().num_labels(meta::type_v<NodeIndex>);
  }
  /// @brief Number of the edge labels.
  constexpr size_t num_edge_labels() const noexcept {
    return _self().num_labels(meta::type_v<EdgeIndex>);
  }
  /// @brief Number of the face labels.
  constexpr size_t num_face_labels() const noexcept {
    return _self().num_labels(meta::type_v<FaceIndex<Derived>>);
  }
  /// @brief Number of the cell labels.
  constexpr size_t num_cell_labels() const noexcept {
    return _self().num_labels(meta::type_v<CellIndex<Derived>>);
  }

  /// @brief Number of the nodes.
  constexpr size_t num_nodes() const noexcept {
    return _self().num_entities(meta::type_v<NodeIndex>);
  }
  /// @brief Number of the nodes with a @p label.
  constexpr size_t num_nodes(Label label) const noexcept {
    return _self().num_entities(label, meta::type_v<NodeIndex>);
  }

  /// @brief Number of the edges.
  constexpr size_t num_edges() const noexcept {
    return _self().num_entities(meta::type_v<EdgeIndex>);
  }
  /// @brief Number of the edges with a @p label.
  constexpr size_t num_edges(Label label) const noexcept {
    return _self().num_entities(label, meta::type_v<EdgeIndex>);
  }

  /// @brief Number of the faces.
  constexpr size_t num_faces() const noexcept {
    return _self().num_entities(meta::type_v<FaceIndex<Derived>>);
  }
  /// @brief Number of the faces with a @p label.
  constexpr size_t num_faces(Label label) const noexcept {
    return _self().num_entities(label, meta::type_v<FaceIndex<Derived>>);
  }

  /// @brief Number of the cells.
  constexpr size_t num_cells() const noexcept {
    return _self().num_entities(meta::type_v<CellIndex<Derived>>);
  }
  /// @brief Number of the cells with a @p label.
  constexpr size_t num_cells(Label label) const noexcept {
    return _self().num_entities(label, meta::type_v<CellIndex<Derived>>);
  }

  /// @brief Range of the nodes.
  constexpr auto nodes() const noexcept {
    return _self().entities(meta::type_v<NodeIndex>) |
           _detail::_to_node_view(_self());
  }
  /// @brief Range of the nodes with a @p label.
  constexpr auto nodes(Label label) const noexcept {
    return _self().entities(label, meta::type_v<NodeIndex>) |
           _detail::_to_node_view(_self());
  }

  /// @brief Range of the edges.
  constexpr auto edges() const noexcept {
    return _self().entities(meta::type_v<EdgeIndex>) |
           _detail::_to_edge_view(_self());
  }
  /// @brief Range of the edges with @p label.
  constexpr auto edges(Label label) const noexcept {
    return _self().entities(label, meta::type_v<EdgeIndex>) |
           _detail::_to_edge_view(_self());
  }

  /// @brief Range of the faces.
  constexpr auto faces() const noexcept {
    return _self().entities(meta::type_v<FaceIndex<Derived>>) |
           _detail::_to_face_view(_self());
  }
  /// @brief Range of the faces with @p label.
  constexpr auto faces(Label label) const noexcept {
    return _self().entities(label, meta::type_v<FaceIndex<Derived>>) |
           _detail::_to_face_view(_self());
  }

  /// @brief Range of the cells.
  constexpr auto cells() const noexcept {
    return _self().entities(meta::type_v<CellIndex<Derived>>) |
           _detail::_to_cell_view(_self());
  }
  /// @brief Range of the cells with @p label.
  constexpr auto cells(Label label) const noexcept {
    return _self().entities(label, meta::type_v<CellIndex<Derived>>) |
           _detail::_to_cell_view(_self());
  }

  /// @brief Range of the interior nodes.
  constexpr auto interior_nodes() const noexcept {
    return nodes(Label{0});
  }

  /// @brief Range of the interior edges.
  constexpr auto interior_edges() const noexcept {
    return edges(Label{0});
  }

  /// @brief Range of the interior faces.
  constexpr auto interior_faces() const noexcept {
    return faces(Label{0});
  }

  /// @brief Range of the interior cells.
  constexpr auto interior_cells() const noexcept {
    return cells(Label{0});
  }

  /// @brief Range of the boundary nodes.
  constexpr auto boundary_nodes() const noexcept {
    return nodes() | std::views::drop(num_nodes(Label{0}));
  }

  /// @brief Range of the boundary edges.
  constexpr auto boundary_edges() const noexcept {
    return edges() | std::views::drop(num_edges(Label{0}));
  }

  /// @brief Range of the boundary faces.
  constexpr auto boundary_faces() const noexcept {
    return faces() | std::views::drop(num_faces(Label{0}));
  }

  /// @brief Range of the boundary cells.
  constexpr auto boundary_cells() const noexcept {
    return cells() | std::views::drop(num_cells(Label{0}));
  }

}; // class MeshInterface

// -----------------------------------------------------------------------------

} // namespace Storm
