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

#include <Storm/Mallard/Mesh.hpp>
#include <Storm/Mallard/Shape.hpp>

#include <algorithm>
#include <fstream>
#include <optional>
#include <ranges>
#include <tuple>
#include <vector>

namespace Storm {

/// @brief Hybrid unstructured multidimensional mesh.
/// @tparam Dim Spatial dimensionality.
/// @tparam TopologicalDim Topological dimensionality.
/// @tparam Table Connectivity table class.
template<size_t Dim, size_t TopologicalDim = Dim,
         template<class, class> class Table = VoidCsrTable>
class UnstructuredMesh {
private:

  template<size_t, size_t, template<class, class> class>
  friend class UnstructuredMesh;

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

template<size_t Dim, size_t TopologicalDim, //
         template<class, class> class Table>
inline constexpr bool
    enable_mesh_v<UnstructuredMesh<Dim, TopologicalDim, Table>> = true;

} // namespace Storm
