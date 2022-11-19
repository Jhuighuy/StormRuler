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
// FITNESS FOR Allocator PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT
// SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
// OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

#pragma once

#include <Storm/Base.hpp>

#include <Storm/Utils/AABB.hpp>
#include <Storm/Utils/Index.hpp>
#include <Storm/Utils/IndexedContainers.hpp>
#include <Storm/Utils/Meta.hpp>
#include <Storm/Utils/Permutations.hpp>
#include <Storm/Utils/Table.hpp>

#include <Storm/Mallard/Mesh.hpp>
#include <Storm/Mallard/Shape.hpp>

#include <glm/glm.hpp>

#include <algorithm>
#include <numeric>
#include <optional>
#include <ranges>
#include <tuple>
#include <utility>
#include <vector>

namespace Storm {

/// @brief Hybrid unstructured multidimensional mesh.
/// @tparam Dim Spatial dimensionality.
/// @tparam TopologicalDim Topological dimensionality.
/// @tparam Table Connectivity table class.
template<size_t Dim, size_t TopologicalDim = Dim,
         template<class, class> class Table = CsrTable>
  requires (2 <= TopologicalDim && TopologicalDim <= Dim && Dim <= 3)
class UnstructuredMesh final :
    public MeshInterface<UnstructuredMesh<Dim, TopologicalDim, Table>> {
private:

  template<size_t, size_t, template<class, class> class>
  friend class UnstructuredMesh;

public:

  /// @brief Spatial vector type.
  using Vec = glm::vec<Dim, real_t>;

  /// @brief Spatial vector type.
  using Mat = glm::mat<Dim, Dim, real_t>;

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
      meta::pair_list_t<Label, EntityIndices_>
    >
  > entity_ranges_tuple_{};

  meta::as_std_tuple_t<
    meta::prepend_t<
      meta::empty_t,
      meta::transform_t<
        meta::pair_cast_fn<IndexedVector>,
        meta::pair_list_t<meta::drop_first_t<EntityIndices_>, shapes::Type>
      >
    >
  > entity_shape_types_tuple_{};

  meta::as_std_tuple_t<
    meta::prepend_t<
      meta::empty_t,
      meta::transform_t<
        meta::pair_cast_fn<IndexedVector>,
        meta::pair_list_t<meta::drop_first_t<EntityIndices_>, real_t>
      >
    >
  > entity_volumes_tuple_{};

  meta::as_std_tuple_t<
    meta::transform_t<
      meta::pair_cast_fn<IndexedVector>,
      meta::pair_list_t<EntityIndices_, Vec>
    >
  > entity_positions_tuple_{};

  // You don't belong here, huh?
  IndexedVector<FaceIndex, Vec> face_normals_{};

  AABB<Vec> aabb_{};

  meta::as_std_tuple_t<
    meta::transform_t<
      meta::pair_cast_fn<Table>,
      meta::cartesian_product_t<EntityIndices_, EntityIndices_>
    >
  > connectivity_tuple_{};

  // clang-format on

public:

  /// @brief Construct the mesh.
  constexpr UnstructuredMesh() {
    // Initialize the default label.
    meta::for_each<EntityIndices_>(
        [&]<size_t I>(meta::type<EntityIndex<I>>) { insert_label<I>(); });
  }

  /// @brief Number of entity labels.
  template<size_t I>
  [[nodiscard]] constexpr size_t
  num_labels(meta::type<EntityIndex<I>> = {}) const noexcept {
    return std::get<I>(entity_ranges_tuple_).size() - 1;
  }

  /// @brief Label range.
  template<size_t I>
  [[nodiscard]] constexpr auto
  labels(meta::type<EntityIndex<I>> = {}) const noexcept {
    const auto& entity_ranges = std::get<I>(entity_ranges_tuple_);
    return std::views::iota(Label{0}, Label{entity_ranges.size() - 1});
  }

  /// @brief Number of entities.
  template<size_t I>
  [[nodiscard]] constexpr size_t
  num_entities(meta::type<EntityIndex<I>> = {}) const noexcept {
    return static_cast<size_t>(std::get<I>(entity_ranges_tuple_).back());
  }
  /// @brief Number of entities with label @p label.
  template<size_t I>
  [[nodiscard]] constexpr size_t
  num_entities(Label label, meta::type<EntityIndex<I>> = {}) const noexcept {
    STORM_ASSERT_(label < num_labels<I>(), "Label is out of range!");
    const auto& entity_ranges = std::get<I>(entity_ranges_tuple_);
    return entity_ranges[label + 1] - entity_ranges[label];
  }

  /// @brief Index range of the entitites.
  template<size_t I>
  [[nodiscard]] constexpr auto
  entities(meta::type<EntityIndex<I>> = {}) const noexcept {
    const auto& entity_ranges = std::get<I>(entity_ranges_tuple_);
    return std::views::iota(entity_ranges.front(), entity_ranges.back());
  }
  /// @brief Index range of the entitites with label @p label.
  template<size_t I>
  [[nodiscard]] constexpr auto
  entities(Label label, meta::type<EntityIndex<I>> = {}) const noexcept {
    STORM_ASSERT_(label < num_labels<I>(), "Label is out of range!");
    const auto& entity_ranges = std::get<I>(entity_ranges_tuple_);
    return std::views::iota(entity_ranges[label], entity_ranges[label + 1]);
  }

  /// @brief Label of the entitity at @p index.
  template<size_t I>
  [[nodiscard]] constexpr Label label(EntityIndex<I> index) const noexcept {
    STORM_ASSERT_(index < num_entities<I>(), "Entity index is out of range!");
    // Binary search for entity in the label ranges.
    const auto& entity_ranges = std::get<I>(entity_ranges_tuple_);
    const auto lower_bound = std::ranges::lower_bound(entity_ranges, index);
    return Label{lower_bound - entity_ranges.begin() - 1};
  }

  /// @brief Shape type of the entitity at @p index.
  template<size_t I>
  [[nodiscard]] constexpr shapes::Type
  shape_type(EntityIndex<I> index) const noexcept {
    STORM_ASSERT_(index < num_entities<I>(), "Entity index is out of range!");
    return std::get<I>(entity_shape_types_tuple_)[index];
  }

  /// @brief "Volume" of the entitity at @p index.
  template<size_t I>
  [[nodiscard]] constexpr real_t volume(EntityIndex<I> index) const noexcept {
    STORM_ASSERT_(index < num_entities<I>(), "Entity index is out of range!");
    return std::get<I>(entity_volumes_tuple_)[index];
  }

  /// @brief Position of the entitity at @p index.
  template<size_t I>
  [[nodiscard]] constexpr Vec position(EntityIndex<I> index) const noexcept {
    STORM_ASSERT_(index < num_entities<I>(), "Entity index is out of range!");
    return std::get<I>(entity_positions_tuple_)[index];
  }

  /// @brief Normal to the face at @p face_index.
  [[nodiscard]] constexpr Vec normal(FaceIndex face_index) const noexcept {
    STORM_ASSERT_(face_index < num_entities(meta::type_v<FaceIndex>),
                  "Face index is out of range!");
    return face_normals_[face_index];
  }

  /// @brief Mesh AABB.
  [[nodiscard]] constexpr const auto& aabb() const noexcept {
    return aabb_;
  }

  /// @brief Range of adjacent entity indices of dim J of an entity at @p index.
  template<size_t J, size_t I>
  [[nodiscard]] constexpr auto
  adjacent(EntityIndex<I> index,
           meta::type<EntityIndex<J>> = {}) const noexcept {
    STORM_ASSERT_(index < num_entities<I>(), "Entity index is out of range!");
    using T = Table<EntityIndex<I>, EntityIndex<J>>;
    return std::get<T>(connectivity_tuple_)[index];
  }

  /// @brief Find an entity by it's @p node_indices. Complexity is constant.
  template<size_t I, std::ranges::input_range Range>
    requires (I != 0) &&
             std::same_as<std::ranges::range_value_t<Range>, NodeIndex>
  [[nodiscard]] constexpr std::optional<EntityIndex<I>> find(
      Range&& node_indices, meta::type<EntityIndex<I>> = {}) const {
    // Select the entities that are adjacent to the first node in the list.
    auto adj = adjacent<I>(node_indices.front());
    STORM_CPP23_THREAD_LOCAL_ std::vector<EntityIndex<I>> found{};
    found.assign(adj.begin(), adj.end());
    std::ranges::sort(found);

    // For the other node indices, select the adjacent enitites,
    // and intersect with the recently found.
    STORM_CPP23_THREAD_LOCAL_ std::vector<EntityIndex<I>> temp{}, update{};
    for (NodeIndex node_index : node_indices | std::views::drop(1)) {
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
    }
    STORM_THROW_("For the specified node list more than one entity found!");
  }

  /// @brief Copy-assign the unstructured @p mesh with different table type.
  template<template<class, class> class OtherTable>
  constexpr void
  assign(const UnstructuredMesh<Dim, TopologicalDim, OtherTable>& mesh) {
    // Copy the ranges and shape properties.
    entity_ranges_tuple_ = mesh.entity_ranges_tuple_;
    entity_shape_types_tuple_ = mesh.entity_shape_types_tuple_;
    entity_volumes_tuple_ = mesh.entity_volumes_tuple_;
    entity_positions_tuple_ = mesh.entity_positions_tuple_;
    face_normals_ = mesh.face_normals_;
    aabb_ = mesh.aabb_;

    // Copy the connectivity tables.
    meta::for_each<EntityIndices_>([&]<size_t I>(meta::type<EntityIndex<I>>) {
      meta::for_each<EntityIndices_>([&]<size_t J>(meta::type<EntityIndex<J>>) {
        using T = Table<EntityIndex<I>, EntityIndex<J>>;
        using U = OtherTable<EntityIndex<I>, EntityIndex<J>>;
        std::get<T>(connectivity_tuple_)
            .assign(std::get<U>(mesh.connectivity_tuple_));
      });
    });
  }

  /// @brief Move-assign the unstructured @p mesh with different table type.
  template<template<class, class> class OtherTable>
  constexpr void
  assign(UnstructuredMesh<Dim, TopologicalDim, OtherTable>&& mesh) {
    // Move the ranges and shape properties.
    /// @todo Ranges of the source mesh would be invalid.
    entity_ranges_tuple_ = std::move(mesh.entity_ranges_tuple_);
    entity_shape_types_tuple_ = std::move(mesh.entity_shape_types_tuple_);
    entity_volumes_tuple_ = std::move(mesh.entity_volumes_tuple_);
    entity_positions_tuple_ = std::move(mesh.entity_positions_tuple_);
    face_normals_ = std::move(mesh.face_normals_);
    aabb_ = std::move(mesh.aabb_);

    // Move the connectivity tables.
    meta::for_each<EntityIndices_>([&]<size_t I>(meta::type<EntityIndex<I>>) {
      meta::for_each<EntityIndices_>([&]<size_t J>(meta::type<EntityIndex<J>>) {
        using T = Table<EntityIndex<I>, EntityIndex<J>>;
        using U = OtherTable<EntityIndex<I>, EntityIndex<J>>;
        std::get<T>(connectivity_tuple_)
            .assign(std::move(std::get<U>(mesh.connectivity_tuple_)));
      });
    });
  }

  /// @brief Assign the arbitrary @p mesh.
  template<mesh Mesh>
  constexpr void assign(const Mesh& mesh); // not implemented yet!

  /// @brief Reverve memory for the entities.
  template<size_t I>
  constexpr void reserve(size_t capacity, meta::type<EntityIndex<I>> = {}) {
    std::get<I>(entity_positions_tuple_).reserve(capacity);
    if constexpr (!std::is_same_v<EntityIndex<I>, NodeIndex>) {
      std::get<I>(entity_shape_types_tuple_).reserve(capacity);
      std::get<I>(entity_volumes_tuple_).reserve(capacity);
    }
    if constexpr (std::is_same_v<EntityIndex<I>, FaceIndex>) {
      face_normals_.reserve(capacity);
    }
    meta::for_each<EntityIndices_>([&]<size_t J>(meta::type<EntityIndex<J>>) {
      using T = Table<EntityIndex<I>, EntityIndex<J>>;
      std::get<T>(connectivity_tuple_).reserve(capacity);
    });
  }

  /// @brief Insert a new entity label.
  template<size_t I>
  constexpr void insert_label(meta::type<EntityIndex<I>> = {}) {
    auto& entity_ranges = std::get<I>(entity_ranges_tuple_);
    if (entity_ranges.empty()) {
      entity_ranges.reserve(2);
      entity_ranges.emplace_back(0), entity_ranges.emplace_back(0);
    } else {
      entity_ranges.push_back(entity_ranges.back());
    }
  }

  /// @brief Insert a new entity of shape @p shape.
  /// The last existing entity label will be assigned to it.
  /// @returns Index of the entity.
  template<size_t I, class Shape>
    requires ((I == 0 && std::constructible_from<Vec, Shape>) ||
              (shapes::shape<Shape> && I == shapes::shape_dim_v<Shape>) )
  constexpr EntityIndex<I> insert(Shape&& shape,
                                  meta::type<EntityIndex<I>> = {}) {
    // Allocate the entity index,
    // implicitly assigning the last existing label label.
    const EntityIndex<I> index{std::get<I>(entity_ranges_tuple_).back()++};

    // Assign the geometrical properties.
    if constexpr (std::is_same_v<EntityIndex<I>, NodeIndex>) {
      // Assign the node position and update the AABB.
      auto& node_positions = std::get<I>(entity_positions_tuple_);
      const Vec& position = node_positions.emplace_back(shape);
      if (node_positions.size() == 1) {
        aabb_ = AABB{position};
      } else {
        aabb_.extend(position);
      }
    } else {
      // Assign the entity shape type, volume, center position (and normal).
      std::get<I>(entity_shape_types_tuple_).emplace_back(shape.type());
      std::get<I>(entity_volumes_tuple_)
          .emplace_back(shapes::volume(shape, *this));
      std::get<I>(entity_positions_tuple_)
          .emplace_back(shapes::barycenter(shape, *this));
      if constexpr (std::is_same_v<EntityIndex<I>, FaceIndex>) {
        if constexpr (Dim == TopologicalDim) {
          face_normals_.emplace_back(shapes::normal(shape, *this));
        } else {
          /// @todo What are the face normals for surface meshes?
          face_normals_.emplace_back({});
        }
      }
    }

    // Assign the connectivity properties.
    meta::for_each<EntityIndices_>([&]<size_t J>(meta::type<EntityIndex<J>>) {
      // Allocate the empty rows.
      using T = Table<EntityIndex<I>, EntityIndex<J>>;
      std::get<T>(connectivity_tuple_).push_back();
    });
    if constexpr (!std::is_same_v<EntityIndex<I>, NodeIndex>) {
      // Connnect the entity with it's nodes.
      using T = Table<EntityIndex<I>, NodeIndex>;
      using U = Table<NodeIndex, EntityIndex<I>>;
      for (NodeIndex node_index : shape.nodes()) {
        std::get<T>(connectivity_tuple_).insert(index, node_index);
        std::get<U>(connectivity_tuple_).insert(node_index, index);
      }
    }
    meta::for_each<meta::make_seq_t<EntityIndex, 1, I>>(
        [&]<size_t J>(meta::type<EntityIndex<J>>) {
          // Connect the entity with it's parts (possibly creating them).
          using T = Table<EntityIndex<I>, EntityIndex<J>>;
          using U = Table<EntityIndex<J>, EntityIndex<I>>;
          const auto process_part = [&](const auto& part) {
            const EntityIndex<J> part_index = find_or_insert<J>(part);
            std::get<T>(connectivity_tuple_).insert(index, part_index);
            std::get<U>(connectivity_tuple_).insert(part_index, index);
            // Fix the face orientation (if needed).
            if constexpr (std::is_same_v<EntityIndex<J>, FaceIndex>) {
              update_face_orientation_(part_index, part.nodes());
            }
          };
          std::apply([&](const auto&... parts) { (process_part(parts), ...); },
                     shapes::parts<J>(shape));
        });

    // Connect the entity with it's siblings.
    /// @todo Implement me!
    using T = Table<EntityIndex<I>, EntityIndex<I>>;
    std::get<T>(connectivity_tuple_).insert(index, index);

    return index;
  }

  /// @brief Find an existing entity of shape @p shape or insert a new one.
  /// @returns Index of the entity.
  template<size_t I, shapes::shape Shape>
    requires (I == shapes::shape_dim_v<Shape>)
  constexpr EntityIndex<I> find_or_insert(Shape&& shape,
                                          meta::type<EntityIndex<I>> = {}) {
    if (const auto index = find<I>(shape.nodes()); index.has_value()) {
      return *index;
    }
    return insert<I>(std::forward<Shape>(shape));
  }

  /// @brief Permute the entities. Complexity is log²-linear.
  /// @param perm Entity permutation range, it may be modified in order to
  ///             preserve the label ranges correctness.
  /// @warning Permutation range values will be spoiled after the call!
  template<size_t I, std::ranges::sized_range Range>
    requires std::ranges::random_access_range<Range> &&
             std::permutable<std::ranges::iterator_t<Range>> &&
             std::same_as<std::ranges::range_value_t<Range>, EntityIndex<I>>
  constexpr void permute(Range&& perm, meta::type<EntityIndex<I>> = {}) {
    STORM_ASSERT_(std::ranges::size(perm) == num_entities<I>(),
                  "Invalid permutation size!");

    // Stable-sort the permutation in order to keep the label ranges correct.
    std::ranges::stable_sort(
        perm, [&](EntityIndex<I> index_1, EntityIndex<I> index_2) {
          return label(index_1) < label(index_2);
        });

    // Permute the entities.
    permute_base_(perm);
  }

  /// @brief Assign the entity labels. Complexity is log²-linear.
  /// @param label Labels range to assign. May be smaller,
  ///              than the amount of entities.
  template<size_t I, std::ranges::sized_range Range>
    requires std::ranges::random_access_range<Range> &&
             std::same_as<std::ranges::range_value_t<Range>, Label>
  constexpr void assign_labels(Range&& labels,
                               meta::type<EntityIndex<I>> = {}) {
    STORM_ASSERT_(std::ranges::size(labels) <= num_entities<I>(),
                  "Invalid labels range size!");

    // Permute the entities according to labels.
    std::vector<EntityIndex<I>> perm(num_entities<I>());
    std::ranges::copy(entities<I>(), perm.begin());
    std::ranges::stable_sort(
        perm, [&](EntityIndex<I> index_1, EntityIndex<I> index_2) {
          const auto new_label = [&](EntityIndex<I> index) {
            const auto index_sz = static_cast<size_t>(index);
            return index >= labels.size() ? Label{0} : labels[index_sz];
          };
          return new_label(index_1) < new_label(index_2);
        });
    permute_base_<I>(perm);

    // Generate the new label ranges.
    const size_t delta = this->num_entities<I>() - std::ranges::size(labels);
    auto& entity_ranges = std::get<I>(entity_ranges_tuple_);
    entity_ranges.clear();
    for (Label label : labels) {
      entity_ranges.resize(std::max(entity_ranges.size(), //
                                    static_cast<size_t>(label) + 2));
      entity_ranges[label + 1] += 1;
    };
    entity_ranges[Label{0} + 1] += delta;
    std::partial_sum(entity_ranges.begin(), entity_ranges.end(),
                     entity_ranges.begin(),
                     [](EntityIndex<I> index_1, EntityIndex<I> index_2) {
                       return index_1 + static_cast<size_t>(index_2);
                     });
  }

private:

  // Update the face orientation such that if it has a single adjacent cell,
  // than the cell is inner.
  template<std::ranges::forward_range Range>
    requires std::ranges::sized_range<Range> &&
             std::same_as<std::ranges::range_value_t<Range>, NodeIndex>
  constexpr void update_face_orientation_(FaceIndex face_index,
                                          Range&& cell_face_nodes) {
    // Detect the face orientation.
    const auto face_nodes = adjacent<0>(face_index);
    bool cell_is_inner, cell_is_outer;
    const size_t num_cell_faces = std::ranges::size(cell_face_nodes);
    STORM_ASSERT_(num_cell_faces >= 2, "This function should be called in 1D!");
    if (num_cell_faces == 2) {
      cell_is_inner = std::ranges::equal(cell_face_nodes, face_nodes);
      cell_is_outer =
          std::ranges::equal(cell_face_nodes, face_nodes | std::views::reverse);
    } else {
      // Need to take into the account possible rotations.
      STORM_CPP23_THREAD_LOCAL_ std::vector<NodeIndex> temp{};
      temp.assign(cell_face_nodes.begin(), cell_face_nodes.end());
      std::ranges::rotate(temp, std::ranges::find(temp, face_nodes.front()));
      cell_is_inner = std::ranges::equal(temp, face_nodes);
      cell_is_outer = !cell_is_inner && [&] {
        std::ranges::rotate(temp, std::ranges::find(temp, face_nodes.back()));
        return std::ranges::equal(temp, face_nodes | std::views::reverse);
      }();
    }
    // Check and update the orientation (if needed).
    const auto cell_faces = adjacent<TopologicalDim>(face_index);
    if (cell_faces.size() == 1) {
      STORM_ENSURE_(cell_is_inner || cell_is_outer,
                    "Face has a single adjacent cell, "
                    "but it can be neither inner nor outer!");
      // Flip the face if the cell is outer.
      if (cell_is_outer) {
        face_normals_[face_index] = -face_normals_[face_index];
        meta::for_each<meta::make_seq_t<EntityIndex, 0, TopologicalDim - 1>>(
            [&]<size_t I>(meta::type<EntityIndex<I>>) {
              using T = Table<FaceIndex, EntityIndex<I>>;
              std::ranges::reverse(
                  std::get<T>(connectivity_tuple_)[face_index]);
            });
      }
    } else if (cell_faces.size() == 2) {
      STORM_ENSURE_(cell_is_outer,
                    "Face has two adjacent cells, "
                    "but the second cell cannot be the outer one!");
    } else {
      STORM_TERMINATE_("Invalid number of the face cells!");
    }
  }

  // Actually permute the entities.
  template<size_t I, std::ranges::random_access_range Range>
    requires std::permutable<std::ranges::iterator_t<Range>> &&
             std::same_as<std::ranges::range_value_t<Range>, EntityIndex<I>>
  constexpr void permute_base_(Range&& perm, meta::type<EntityIndex<I>> = {}) {
    { // Update the connectivity with inverse permutations.
      IndexedVector<EntityIndex<I>, EntityIndex<I>> iperm(num_entities<I>());
      permutations::invert_permutation(perm, iperm.begin());
      meta::for_each<EntityIndices_>([&]<size_t J>(meta::type<EntityIndex<J>>) {
        /// @todo Use parallel loop here!
        std::ranges::for_each(entities<J>(), [&](EntityIndex<J> index) {
          using T = Table<EntityIndex<J>, EntityIndex<I>>;
          decltype(auto) adj = std::get<T>(connectivity_tuple_)[index];
          for (EntityIndex<I>& adjacent_index : adj) {
            adjacent_index = iperm[adjacent_index];
          }
        });
      });
    }

    // Permute the entity-* connectivity.
    meta::for_each<EntityIndices_>([&]<size_t J>(meta::type<EntityIndex<J>>) {
      using T = Table<EntityIndex<I>, EntityIndex<J>>;
      auto& table = std::get<T>(connectivity_tuple_);
      T permuted_table{};
      permuted_table.reserve(num_entities<I>());
      for (EntityIndex<I> index : entities<I>()) {
        const auto index_sz = static_cast<size_t>(index);
        const EntityIndex<I> permuted_index = perm[index_sz];
        permuted_table.push_back(std::move(table[permuted_index]));
      }
      table = std::move(permuted_table);
    });

    // Permute the entity shape properties.
    permutations::permute_inplace(
        perm, [&](EntityIndex<I> index_1, EntityIndex<I> index_2) {
          const auto gather_shape_properties = [&](EntityIndex<I> index) {
            auto& positions = std::get<I>(entity_positions_tuple_);
            if constexpr (std::is_same_v<EntityIndex<I>, NodeIndex>) {
              return std::tie(positions[index]);
            } else {
              auto& shape_types = std::get<I>(entity_shape_types_tuple_);
              auto& volumes = std::get<I>(entity_volumes_tuple_);
              if constexpr (std::is_same_v<EntityIndex<I>, FaceIndex>) {
                return std::tie(shape_types[index], volumes[index],
                                positions[index], face_normals_[index]);
              } else {
                return std::tie(shape_types[index], volumes[index],
                                positions[index]);
              }
            }
          };
          auto properties_1 = gather_shape_properties(index_1);
          auto properties_2 = gather_shape_properties(index_2);
          std::swap(properties_1, properties_2);
        });
  }

}; // class UnstructuredMesh

} // namespace Storm
