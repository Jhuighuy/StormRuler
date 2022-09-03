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

#include <Storm/Mallard/Mesh.hpp>

#include <compare>
#include <concepts>
#include <ranges>

namespace Storm {

template<mesh Mesh>
class NodeView;
template<class Mesh>
NodeView(Mesh&, NodeIndex) -> NodeView<Mesh>;

template<mesh Mesh>
class EdgeView;
template<class Mesh>
EdgeView(Mesh&, EdgeIndex) -> EdgeView<Mesh>;

template<mesh Mesh>
class FaceView;
template<class Mesh>
FaceView(Mesh&, FaceIndex<Mesh>) -> FaceView<Mesh>;

template<mesh Mesh>
class CellView;
template<class Mesh>
CellView(Mesh&, CellIndex<Mesh>) -> CellView<Mesh>;

namespace detail_ {
  constexpr auto to_node_view_(auto& mesh) noexcept {
    return std::views::transform(
        [&](NodeIndex node_index) { return NodeView(mesh, node_index); });
  }
  constexpr auto to_edge_view_(auto& mesh) noexcept {
    return std::views::transform(
        [&](EdgeIndex edge_index) { return EdgeView(mesh, edge_index); });
  }
  template<class Mesh>
  constexpr auto to_face_view_(Mesh& mesh) noexcept {
    return std::views::transform(
        [&](FaceIndex<Mesh> face_index) { return FaceView(mesh, face_index); });
  }
  template<class Mesh>
  constexpr auto to_cell_view_(Mesh& mesh) noexcept {
    return std::views::transform(
        [&](CellIndex<Mesh> cell_index) { return CellView(mesh, cell_index); });
  }
} // namespace detail_

/// @brief Mesh entity view.
template<mesh Mesh, index Index>
class EntityView;

// clang-format off
template<mesh Mesh, size_t I>
  requires std::is_object_v<Mesh>
class EntityView<Mesh, EntityIndex<I>> {
  // clang-format on
private:

  template<mesh, index>
  friend class EntityView;

  Mesh* p_mesh_;
  EntityIndex<I> index_;

protected:

  /// @brief Construct an entity view with a @p mesh and @p index.
  constexpr EntityView(Mesh& mesh, EntityIndex<I> index) noexcept
      : p_mesh_{&mesh}, index_{index} {}

  /// @brief Copy-construct an entity view.
  /// @todo Do we need this?
  template<class OtherMesh>
  constexpr EntityView(
      const EntityView<OtherMesh, EntityIndex<I>>& other) noexcept
      : p_mesh_{other.p_mesh_}, index_{other.index_} {}

  /// @brief Destroy an entity view.
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
  [[nodiscard]] constexpr EntityIndex<I> index() const noexcept {
    return index_;
  }

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
  constexpr NodeView(Mesh& mesh, NodeIndex index) noexcept
      : EntityView<Mesh, NodeIndex>(mesh, index) {}

  /// @brief Copy-construct a node view.
  constexpr NodeView(const NodeView<std::remove_const_t<Mesh>>& other) noexcept
      : EntityView<Mesh, NodeIndex>(other) {}

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
  constexpr EdgeView(Mesh& mesh, EdgeIndex index) noexcept
      : EntityView<Mesh, EdgeIndex>(mesh, index) {}

  /// @brief Copy-construct an edge view.
  constexpr EdgeView(const EdgeView<std::remove_const_t<Mesh>>& other) noexcept
      : EntityView<Mesh, EdgeIndex>(other) {}

  /// @brief Edge length.
  [[nodiscard]] constexpr real_t length() const noexcept {
    return this->mesh().volume(this->index());
  }

  /// @brief Edge barycenter position.
  [[nodiscard]] constexpr auto barycenter_position() const noexcept {
    return this->mesh().position(this->index());
  }

}; // class EdgeView

/// @brief Face view.
template<mesh Mesh>
class FaceView final : public EntityView<Mesh, FaceIndex<Mesh>> {
public:

  /// @brief Construct a face view.
  constexpr FaceView(Mesh& mesh, FaceIndex<Mesh> index) noexcept
      : EntityView<Mesh, FaceIndex<Mesh>>(mesh, index) {}

  /// @brief Copy-construct a face view.
  constexpr FaceView(const FaceView<std::remove_const_t<Mesh>>& other) noexcept
      : EntityView<Mesh, FaceIndex<Mesh>>(other) {}

  /// @brief Face area (or length in 2D).
  [[nodiscard]] constexpr real_t area() const noexcept {
    return this->mesh().volume(this->index());
  }

  /// @brief Face barycenter position.
  [[nodiscard]] constexpr auto barycenter_position() const noexcept {
    return this->mesh().position(this->index());
  }

  /// @brief Face normal.
  [[nodiscard]] constexpr auto normal() const noexcept {
    return this->mesh().normal(this->index());
  }

  /// @brief Get the adjacent inner cell.
  [[nodiscard]] constexpr CellView<Mesh> inner_cell() const noexcept {
    STORM_ASSERT_(detail_::face_inner_cell_ < this->cells().size(),
                  "The face does not have an adjacent inner cell!");
    return this->cells()[detail_::face_inner_cell_];
  }

  /// @brief Get the adjacent outer cell.
  [[nodiscard]] constexpr CellView<Mesh> outer_cell() const noexcept {
    STORM_ASSERT_(detail_::face_outer_cell_ < this->cells().size(),
                  "The face does not have an adjacent outer cell!");
    return this->cells()[detail_::face_outer_cell_];
  }

}; // class FaceView

/// @brief Mesh cell view.
template<mesh Mesh>
class CellView final : public EntityView<Mesh, CellIndex<Mesh>> {
public:

  /// @brief Construct a cell view.
  constexpr CellView(Mesh& mesh, CellIndex<Mesh> index) noexcept
      : EntityView<Mesh, CellIndex<Mesh>>(mesh, index) {}

  /// @brief Copy-construct a cell view.
  constexpr CellView(const CellView<std::remove_const_t<Mesh>>& other) noexcept
      : EntityView<Mesh, CellIndex<Mesh>>(other) {}

  /// @brief Cell volume (or area in 2D).
  [[nodiscard]] constexpr real_t volume() const noexcept {
    return this->mesh().volume(this->index());
  }

  /// @brief Cell barycenter position.
  [[nodiscard]] constexpr auto barycenter_position() const noexcept {
    return this->mesh().position(this->index());
  }

}; // class CellView

/// @brief Range of the @p mesh nodes.
[[nodiscard]] constexpr auto nodes(auto& mesh) noexcept {
  return mesh.entities(meta::type_v<NodeIndex>) | //
         detail_::to_node_view_(mesh);
}
/// @brief Range of the @p mesh nodes with a @p label.
[[nodiscard]] constexpr auto nodes(auto& mesh, Label label) noexcept {
  return mesh.entities(meta::type_v<NodeIndex>, label) |
         detail_::to_node_view_(mesh);
}

/// @brief Range of the @p mesh edges.
[[nodiscard]] constexpr auto edges(auto& mesh) noexcept {
  return mesh.entities(meta::type_v<EdgeIndex>) | //
         detail_::to_edge_view_(mesh);
}
/// @brief Range of the @p mesh edges with @p label.
[[nodiscard]] constexpr auto edges(auto& mesh, Label label) noexcept {
  return mesh.entities(meta::type_v<EdgeIndex>, label) |
         detail_::to_edge_view_(mesh);
}

/// @brief Range of the @p mesh faces.
template<mesh Mesh>
[[nodiscard]] constexpr auto faces(Mesh& mesh) noexcept {
  return mesh.entities(meta::type_v<FaceIndex<Mesh>>) |
         detail_::to_face_view_(mesh);
}
/// @brief Range of the @p mesh faces with @p label.
template<mesh Mesh>
[[nodiscard]] constexpr auto faces(Mesh& mesh, Label label) noexcept {
  return mesh.entities(meta::type_v<FaceIndex<Mesh>>, label) |
         detail_::to_face_view_(mesh);
}

/// @brief Range of the @p mesh cells.
template<mesh Mesh>
[[nodiscard]] constexpr auto cells(auto& mesh) noexcept {
  return mesh.entities(meta::type_v<CellIndex<Mesh>>) |
         detail_::to_cell_view_(mesh);
}
/// @brief Range of the @p mesh cells with @p label.
template<mesh Mesh>
[[nodiscard]] constexpr auto cells(auto& mesh, Label label) noexcept {
  return mesh.entities(meta::type_v<CellIndex<Mesh>>, label) |
         detail_::to_cell_view_(mesh);
}

/// @brief Range of the interior @p mesh nodes.
[[nodiscard]] constexpr auto int_nodes(auto& mesh) noexcept {
  return nodes(mesh, Label{0});
}

/// @brief Range of the interior @p mesh edges.
[[nodiscard]] constexpr auto int_edges(auto& mesh) noexcept {
  return edges(mesh, Label{0});
}

/// @brief Range of the interior @p mesh faces.
[[nodiscard]] constexpr auto int_faces(auto& mesh) noexcept {
  return faces(mesh, Label{0});
}

/// @brief Range of the interior @p mesh cells.
[[nodiscard]] constexpr auto int_cells(auto& mesh) noexcept {
  return cells(mesh, Label{0});
}

/// @brief Range of the boundary @p mesh nodes.
[[nodiscard]] constexpr auto bnd_nodes(auto& mesh) noexcept {
  return nodes(mesh) | std::views::drop(mesh.num_entities( //
                           meta::type_v<NodeIndex>, Label{0}));
}

/// @brief Range of the boundary @p mesh edges.
[[nodiscard]] constexpr auto bnd_edges(auto& mesh) noexcept {
  return edges(mesh) | std::views::drop(mesh.num_entities( //
                           meta::type_v<EdgeIndex>, Label{0}));
}

/// @brief Range of the boundary @p mesh faces.
template<mesh Mesh>
[[nodiscard]] constexpr auto bnd_faces(Mesh& mesh) noexcept {
  return faces(mesh) | std::views::drop(mesh.num_entities(
                           meta::type_v<FaceIndex<Mesh>>, Label{0}));
}

/// @brief Range of the boundary @p mesh cells.
template<mesh Mesh>
[[nodiscard]] constexpr auto bnd_cells(Mesh& mesh) noexcept {
  return cells(mesh) | std::views::drop(mesh.num_entities(
                           meta::type_v<CellIndex<Mesh>>, Label{0}));
}

} // namespace Storm
