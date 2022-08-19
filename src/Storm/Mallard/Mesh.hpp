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

#include <Storm/Utils/Containers.hpp>
#include <Storm/Utils/Meta.hpp>
#include <Storm/Utils/Table.hpp>

#include <Storm/Mallard/Shape.hpp>

#include <algorithm>
#include <compare>
#include <fstream>
#include <optional>
#include <ranges>
#include <tuple>
#include <vector>

namespace Storm {

#if 0
namespace detail_ {
  struct LabelTag_;
  template<size_t I>
  struct TopologicalIndexTag_;
} // namespace detail_

/// @brief Label index type.
using Label = Index<detail_::LabelTag_>;

/// @brief Topological index type.
template<size_t I>
using EntityIndex = Index<detail_::TopologicalIndexTag_<I>>;

/// @brief Node index type.
using NodeIndex = EntityIndex<0>;

/// @brief Edge index type.
using EdgeIndex = EntityIndex<1>;
#endif

/// @brief Face index type.
template<class Mesh>
using FaceIndex = typename Mesh::FaceIndex;

/// @brief Cell index type.
template<class Mesh>
using CellIndex = typename Mesh::CellIndex;

namespace detail_ {
  inline constexpr size_t face_inner_cell_ = 0;
  inline constexpr size_t face_outer_cell_ = 1;
} // namespace detail_

/// @brief Hybrid unstructured multidimensional mesh.
/// @tparam Dim Spatial dimensionality.
/// @tparam TopologicalDim Topological dimensionality.
/// @tparam Table Connectivity table class.
template<size_t Dim, size_t TopologicalDim = Dim,
         template<class, class> class Table = VoidCsrTable>
class UnstructuredMesh {
public:

  /// @brief Spatial vector type.
  using Vec = FastVector<real_t, Dim>;

  /// @brief Face index type.
  using FaceIndex = EntityIndex<TopologicalDim - 1>;

  /// @brief Cell index type.
  using CellIndex = EntityIndex<TopologicalDim>;

private:

  using EntityIndices_ = meta::make_seq_t<EntityIndex, 0, TopologicalDim + 1>;

  // clang-format off

  meta::as_std_tuple_t<
    meta::transform_t<
      meta::pair_cast_fn<IndexedVector>,
      meta::transform_t<
        meta::reverse_fn,
        meta::pair_list_t<Label, EntityIndices_>
      >
    >
  > label_ranges_tuple_{};

  meta::as_std_tuple_t<
    meta::transform_t<
      meta::pair_cast_fn<IndexedVector>,
      meta::pair_list_t<Label, EntityIndices_>
    >
  > entity_labels_tuple_{};

  meta::as_std_tuple_t<
    meta::prepend_t<
      meta::empty_t,
      meta::transform_t<
        meta::pair_cast_fn<IndexedVector>,
        meta::pair_list_t<real_t, meta::drop_first_t<EntityIndices_>>
      >
    >
  > entity_volumes_tuple_{};

  meta::as_std_tuple_t<
    meta::transform_t<
      meta::pair_cast_fn<IndexedVector>,
      meta::pair_list_t<Vec, EntityIndices_>
    >
  > entity_positions_tuple_{};

  // You don't belong here, huh?
  IndexedVector<Vec, FaceIndex> face_normals_{};

  meta::as_std_tuple_t<
    meta::transform_t<
      meta::pair_cast_fn<Table>,
      meta::cartesian_product_t<EntityIndices_, EntityIndices_>
    >
  > connectivity_tuple_{};

  // clang-format on

public:

  /// @brief Initialize the mesh.
  constexpr UnstructuredMesh() {
    // Initialize the default label (0).
    std::apply([](auto&... ranges) { (ranges.emplace_back(0), ...); },
               label_ranges_tuple_);
  }

  /// @brief Number of entity labels.
  template<size_t I>
  [[nodiscard]] constexpr size_t
  num_labels(meta::type<EntityIndex<I>> = {}) const noexcept {
    return std::get<I>(label_ranges_tuple_).size() - 1;
  }

  /// @brief Number of entities.
  template<size_t I>
  [[nodiscard]] constexpr size_t
  num_entities(meta::type<EntityIndex<I>> = {}) const noexcept {
    return static_cast<size_t>(std::get<I>(label_ranges_tuple_).back());
  }
  /// @brief Number of entities with label @p label.
  template<size_t I>
  [[nodiscard]] constexpr size_t
  num_entities(Label label, meta::type<EntityIndex<I>> = {}) const noexcept {
    const auto& label_ranges = std::get<I>(label_ranges_tuple_);
    return label_ranges[label + 1] - label_ranges[label];
  }

  /// @brief Index range of the entitites.
  template<size_t I>
  [[nodiscard]] constexpr auto //
  entities(meta::type<EntityIndex<I>> = {}) const noexcept {
    const auto& label_ranges = std::get<I>(label_ranges_tuple_);
    return std::views::iota(label_ranges.front(), label_ranges.back());
  }
  /// @brief Index range of the entitites with label @p label.
  template<size_t I>
  [[nodiscard]] constexpr auto
  entities(Label label, meta::type<EntityIndex<I>> = {}) const noexcept {
    const auto& label_ranges = std::get<I>(label_ranges_tuple_);
    return std::views::iota(label_ranges[label], label_ranges[label + 1]);
  }

  /// @brief Label of the entitity at @p index.
  template<size_t I>
  [[nodiscard]] constexpr Label //
  label(EntityIndex<I> index) const noexcept {
    return std::get<I>(entity_labels_tuple_)[index];
  }

  /// @brief "Volume" of the entitity at @p index.
  template<size_t I>
  [[nodiscard]] constexpr real_t //
  volume(EntityIndex<I> index) const noexcept {
    return std::get<I>(entity_volumes_tuple_)[index];
  }

  /// @brief Position of the entitity at @p index.
  template<size_t I>
  [[nodiscard]] constexpr Vec //
  position(EntityIndex<I> index) const noexcept {
    return std::get<I>(entity_positions_tuple_)[index];
  }

  /// @brief Normal to the face at @p face_index.
  [[nodiscard]] constexpr Vec //
  normal(FaceIndex face_index) const noexcept {
    return face_normals_[face_index];
  }

  /// @brief Range of adjacent entity indices of dim J of an entity at @p index.
  template<size_t J, size_t I>
  [[nodiscard]] constexpr auto //
  adjacent(EntityIndex<I> index,
           meta::type<EntityIndex<J>> = {}) const noexcept {
    using T = Table<EntityIndex<I>, EntityIndex<J>>;
    return std::get<T>(connectivity_tuple_)[index];
  }

  /// @brief Find an entity by it's @p node_indices.
  template<size_t I, class Range>
  [[nodiscard]] constexpr std::optional<EntityIndex<I>>
  find(const Range& node_indices, meta::type<EntityIndex<I>> = {}) const {
    // Select the entities that are adjacent to the first node in the list.
    auto adj = adjacent<I>(node_indices.front());
    std::vector<EntityIndex<I>> found(adj.begin(), adj.end());
    std::ranges::sort(found);

    // For the other node indices, select the adjacent enitites,
    // and intersect with the recently found.
    for (std::vector<EntityIndex<I>> temp, update;
         NodeIndex node_index : node_indices | std::views::drop(1)) {
      if (found.empty()) { break; }
      adj = adjacent<I>(node_index);
      temp.assign(adj.begin(), adj.end());
      std::ranges::sort(temp);
      std::ranges::set_intersection(found, temp, std::back_inserter(update));
      found.swap(update), update.clear();
    }

    // Return the result.
    switch (found.size()) {
      case 0: return std::nullopt;
      case 1: return found.front();
      default:
        throw std::range_error(
            "For the specified node list more than one entity found!");
    }
  }

  /// @brief Insert a new or find an existing entity of shape @p shape.
  /// @returns Index of the entity.
  // clang-format off
  template<size_t I, class Shape>
    requires ((I == 0 && std::constructible_from<Vec, const Shape&>) ||
              (I != 0 && shapes::shape<Shape>))
  constexpr EntityIndex<I>
  insert(const Shape& shape, meta::type<EntityIndex<I>> = {}) {
    // clang-format on

    // If an edge or or face is inserted, check if it exists first.
    if constexpr (!(std::is_same_v<EntityIndex<I>, NodeIndex> ||
                    std::is_same_v<EntityIndex<I>, CellIndex>) ) {
      if (const auto found_entity_index = find<I>(shape.nodes());
          found_entity_index.has_value()) {
        return *found_entity_index;
      }
    }
    const EntityIndex<I> entity_index{num_entities<I>()};

    // Assign the default label.
    /// @todo Label!
    std::get<I>(label_ranges_tuple_).back() += 1;
    std::get<I>(entity_labels_tuple_).emplace_back(num_labels<I>() - 1);

    // Assign the geometrical properties.
    if constexpr (std::is_same_v<EntityIndex<I>, NodeIndex>) {
      // Assign the node position.
      std::get<I>(entity_positions_tuple_).emplace_back(shape);
    } else {
      // Assign the entity volume, center position (and normal).
      std::get<I>(entity_volumes_tuple_)
          .emplace_back(shapes::volume(shape, *this));
      std::get<I>(entity_positions_tuple_)
          .emplace_back(shapes::barycenter(shape, *this));
      if constexpr (std::is_same_v<EntityIndex<I>, FaceIndex>) {
        face_normals_.emplace_back(shapes::normal(shape, *this));
      }
    }

    // Assign the topological properties.
    meta::for_each<EntityIndices_>([&]<size_t J>(meta::type<EntityIndex<J>>) {
      // Allocate the empty rows.
      using T = Table<EntityIndex<I>, EntityIndex<J>>;
      std::get<T>(connectivity_tuple_).emplace_back();
    });
    if constexpr (!std::is_same_v<EntityIndex<I>, NodeIndex>) {
      // Connnect the entity with it's nodes.
      using T = Table<EntityIndex<I>, NodeIndex>;
      using U = Table<NodeIndex, EntityIndex<I>>;
      for (NodeIndex node_index : shape.nodes()) {
        std::get<T>(connectivity_tuple_).insert(entity_index, node_index);
        std::get<U>(connectivity_tuple_).insert(node_index, entity_index);
      }
    }
    meta::for_each<meta::make_seq_t<EntityIndex, 1, I>>(
        [&]<size_t J>(meta::type<EntityIndex<J>>) {
          // Connect the entity with it's parts (possibly creating them).
          using T = Table<EntityIndex<I>, EntityIndex<J>>;
          using U = Table<EntityIndex<J>, EntityIndex<I>>;
          const auto process_part = [&](const auto& part) {
            const EntityIndex<J> part_index = insert<J>(part);
            std::get<T>(connectivity_tuple_).insert(entity_index, part_index);
            std::get<U>(connectivity_tuple_).insert(part_index, entity_index);
          };
          std::apply([&](const auto&... parts) { (process_part(parts), ...); },
                     shapes::parts<J>(shape));
        });

    // Connect the entity with it's siblings.
    /// @todo Implement me!
    using T = Table<EntityIndex<I>, EntityIndex<I>>;
    std::get<T>(connectivity_tuple_).insert(entity_index, entity_index);

    STORM_INFO_("inserting ({}, {})", I, (size_t) entity_index);
    return entity_index;
  }

#if 1
  void read_from_triangle(const auto& path) {
    std::string line;

    std::ifstream nodeStream(path + std::string("node"));
    STORM_ENSURE_(nodeStream.is_open(), "");
    size_t numNodes{0}, dim{0};
    nodeStream >> numNodes >> dim;
    std::getline(nodeStream, line);
    for (size_t i{0}; i < numNodes; ++i) {
      NodeIndex nodeIndex{0};
      glm::dvec3 nodeCoords(0.0);
      nodeStream >> nodeIndex >> nodeCoords.x >> nodeCoords.y;
      std::getline(nodeStream, line);
      STORM_ENSURE_(nodeIndex == insert<0>(nodeCoords), "");
    }

#if 0
    std::ifstream faceStream(path + std::string("edge"));
    STORM_ENSURE_(faceStream.is_open(), "");
    size_t numFaces{0};
    faceStream >> numFaces;
    std::getline(faceStream, line);
    for (size_t i{0}; i < numFaces; ++i) {
      FaceIndex faceIndex{0};
      std::vector<NodeIndex> faceNodes(2);
      size_t faceMark{0};
      faceStream >> faceIndex >> faceNodes[0] >> faceNodes[1] >> faceMark;
      STORM_ENSURE_(faceIndex == insert_face({ShapeType::Segment, faceNodes},
                                             FaceMark(faceMark)),
                    "");
      std::getline(faceStream, line);
    }
#endif

#if 1
    std::ifstream cellStream(path + std::string("ele"));
    STORM_ENSURE_(cellStream.is_open(), "");
    size_t numCells{0};
    cellStream >> numCells;
    std::getline(cellStream, line);
    for (size_t i{0}; i < numCells; ++i) {
      size_t cellIndex{0};
      shapes::Triangle cell;
      cellStream >> cellIndex >> cell.n1 >> cell.n2 >> cell.n3;
      insert<2>(cell);
      // STORM_ENSURE_(cellIndex == insert<2>(Triangle{cellNodes}), "");
      std::getline(cellStream, line);
    }
#endif
  }
#endif

}; // class UnstructuredMesh

template<class Mesh>
class NodeView;
template<class Mesh>
NodeView(Mesh&, NodeIndex) -> NodeView<Mesh>;

template<class Mesh>
class EdgeView;
template<class Mesh>
EdgeView(Mesh&, EdgeIndex) -> EdgeView<Mesh>;

template<class Mesh>
class FaceView;
template<class Mesh>
FaceView(Mesh&, FaceIndex<Mesh>) -> FaceView<Mesh>;

template<class Mesh>
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
template<class Mesh, class Index>
class EntityView;

// clang-format off
template<class Mesh, size_t I>
  requires std::is_object_v<Mesh>
class EntityView<Mesh, EntityIndex<I>> {
  // clang-format on
private:

  template<class, class>
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
    return p_mesh_->label(index_);
  }

  /// @brief Range of the adjacent nodes.
  [[nodiscard]] constexpr auto nodes() const noexcept {
    return p_mesh_->adjacent(index_, meta::type_v<NodeIndex>) |
           detail_::to_node_view_(*p_mesh_);
  }

  /// @brief Range of the adjacent edges.
  [[nodiscard]] constexpr auto edges() const noexcept {
    return p_mesh_->adjacent(index_, meta::type_v<EdgeIndex>) |
           detail_::to_edge_view_(*p_mesh_);
  }

  /// @brief Range of the adjacent faces.
  [[nodiscard]] constexpr auto faces() const noexcept {
    return p_mesh_->adjacent(index_, meta::type_v<FaceIndex<Mesh>>) |
           detail_::to_face_view_(*p_mesh_);
  }

  /// @brief Range of the adjacent cells.
  [[nodiscard]] constexpr auto cells() const noexcept {
    return p_mesh_->adjacent(index_, meta::type_v<CellIndex<Mesh>>) |
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
template<class Mesh>
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
template<class Mesh>
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
template<class Mesh>
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
template<class Mesh>
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
template<class Mesh>
[[nodiscard]] constexpr auto faces(Mesh& mesh) noexcept {
  return mesh.entities(meta::type_v<FaceIndex<Mesh>>) |
         detail_::to_face_view_(mesh);
}
/// @brief Range of the @p mesh faces with @p label.
template<class Mesh>
[[nodiscard]] constexpr auto faces(Mesh& mesh, Label label) noexcept {
  return mesh.entities(meta::type_v<FaceIndex<Mesh>>, label) |
         detail_::to_face_view_(mesh);
}

/// @brief Range of the @p mesh cells.
template<class Mesh>
[[nodiscard]] constexpr auto cells(auto& mesh) noexcept {
  return mesh.entities(meta::type_v<CellIndex<Mesh>>) |
         detail_::to_cell_view_(mesh);
}
/// @brief Range of the @p mesh cells with @p label.
template<class Mesh>
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
template<class Mesh>
[[nodiscard]] constexpr auto bnd_faces(Mesh& mesh) noexcept {
  return faces(mesh) | std::views::drop(mesh.num_entities(
                           meta::type_v<FaceIndex<Mesh>>, Label{0}));
}

/// @brief Range of the boundary @p mesh cells.
template<class Mesh>
[[nodiscard]] constexpr auto bnd_cells(Mesh& mesh) noexcept {
  return cells(mesh) | std::views::drop(mesh.num_entities(
                           meta::type_v<CellIndex<Mesh>>, Label{0}));
}

} // namespace Storm
