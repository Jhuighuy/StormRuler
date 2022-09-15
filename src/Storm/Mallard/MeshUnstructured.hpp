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

#include <Storm/Utils/Index.hpp>
#include <Storm/Utils/IndexedContainers.hpp>
#include <Storm/Utils/Meta.hpp>
#include <Storm/Utils/Permutations.hpp>
#include <Storm/Utils/Table.hpp>

#include <Storm/Mallard/Mesh.hpp>
#include <Storm/Mallard/Shape.hpp>

#include <algorithm>
#include <numeric>
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
         template<class, class> class Table = VoidVovTable>
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

  shapes::AABB<Vec> aabb_{};

  meta::as_std_tuple_t<
    meta::transform_t<
      meta::pair_cast_fn<IndexedVector>,
      meta::pair_list_t<EntityIndices_, Vec>
    >
  > entity_positions_tuple_{};

  // You don't belong here, huh?
  IndexedVector<FaceIndex, Vec> face_normals_{};

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
    meta::for_each<EntityIndices_>([&]<size_t I>(meta::type<EntityIndex<I>>) {
      insert_label<I>(), insert_label<I>();
    });
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

  /// @brief Mesh AABB.
  [[nodiscard]] constexpr const auto& aabb() const noexcept {
    return aabb_;
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
    std::get<I>(entity_ranges_tuple_).emplace_back(0);
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
    const EntityIndex<I> entity_index{
        std::get<I>(entity_ranges_tuple_).back()++};

    // Assign the geometrical properties.
    if constexpr (std::is_same_v<EntityIndex<I>, NodeIndex>) {
      // Assign the node position and update the AABB.
      auto& positions = std::get<I>(entity_positions_tuple_);
      const Vec& pos = positions.emplace_back(shape);
      if (positions.size() == 1) {
        aabb_ = shapes::AABB{pos};
      } else {
        aabb_.extend(pos);
      }
    } else {
      // Assign the entity shape type, volume, center position (and normal).
      std::get<I>(entity_shape_types_tuple_).emplace_back(shape.type());
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
      std::get<T>(connectivity_tuple_).push_back();
    });
    if constexpr (!std::is_same_v<EntityIndex<I>, NodeIndex>) {
      // Connnect the entity with it's nodes.
      using T = Table<EntityIndex<I>, NodeIndex>;
      using U = Table<NodeIndex, EntityIndex<I>>;
      for (NodeIndex node_index : shape.nodes()) {
        std::get<T>(connectivity_tuple_).connect(entity_index, node_index);
        std::get<U>(connectivity_tuple_).connect(node_index, entity_index);
      }
    }
    meta::for_each<meta::make_seq_t<EntityIndex, 1, I>>(
        [&]<size_t J>(meta::type<EntityIndex<J>>) {
          // Connect the entity with it's parts (possibly creating them).
          using T = Table<EntityIndex<I>, EntityIndex<J>>;
          using U = Table<EntityIndex<J>, EntityIndex<I>>;
          const auto process_part = [&](const auto& part) {
            const EntityIndex<J> part_index = find_or_insert<J>(part);
            std::get<T>(connectivity_tuple_).connect(entity_index, part_index);
            std::get<U>(connectivity_tuple_).connect(part_index, entity_index);
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
    std::get<T>(connectivity_tuple_).connect(entity_index, entity_index);

    return entity_index;
  }

  /// @brief Find an existing entity of shape @p shape or insert a new one.
  /// @returns Index of the entity.
  template<size_t I, class Shape>
    requires ((I == 0 && std::constructible_from<Vec, Shape>) ||
              (shapes::shape<Shape> && I == shapes::shape_dim_v<Shape>) )
  constexpr EntityIndex<I> find_or_insert(Shape&& shape,
                                          meta::type<EntityIndex<I>> = {}) {
    if (const auto entity_index = find<I>(shape.nodes());
        entity_index.has_value()) {
      return *entity_index;
    }
    return insert<I>(std::forward<Shape>(shape));
  }

#if 0
  /// @brief Insert a ghost cell.
  /// @returns Index of the ghost cell.
  constexpr CellIndex insert_ghost(FaceIndex face_index) {
    constexpr size_t I = TopologicalDim;
    const CellIndex cell_index{num_entities<I>()};

    // Assign the last existing label label.
    std::get<I>(entity_ranges_tuple_).back() += 1;
    std::get<I>(entity_labels_tuple_).emplace_back(num_labels<I>() - 1);

    // Assign the entity shape type, volume, center position (and normal).
    const CellIndex mirror_cell_index = adjacent<I>(face_index).front();
    std::get<I>(entity_shape_types_tuple_).emplace_back(shapes::Type::ghost);
    std::get<I>(entity_volumes_tuple_).emplace_back(volume(mirror_cell_index));
    std::get<I>(entity_positions_tuple_)
        .emplace_back(2.0 * position(face_index) - position(mirror_cell_index));

    // Allocate the empty connectivity rows, connect the face with ghost.
    meta::for_each<EntityIndices_>([&]<size_t J>(meta::type<EntityIndex<J>>) {
      using T = Table<CellIndex, EntityIndex<J>>;
      std::get<T>(connectivity_tuple_).push_back();
    });
    using T = Table<FaceIndex, CellIndex>;
    std::get<T>(connectivity_tuple_).connect(face_index, cell_index);

    return cell_index;
  }
#endif

  /// @brief Permute the entities.
  /// @param perm Entity permutation range, it may be modified in order to
  ///             preserve the label ranges correctness.
  /// @warning This operation is very slow!
  template<size_t I, std::ranges::random_access_range Range>
    requires std::permutable<std::ranges::iterator_t<Range>> &&
             std::same_as<std::ranges::range_value_t<Range>, EntityIndex<I>>
  constexpr void permute(Range&& perm, meta::type<EntityIndex<I>> = {}) {
    if constexpr (std::ranges::sized_range<Range>) {
      STORM_ASSERT_(std::ranges::size(perm) == num_entities<I>(),
                    "Invalid permutation size!");
    }

    // Stable-sort the permutation in order to keep the label ranges correct.
    std::ranges::stable_sort(perm, [&](EntityIndex<I> entity_index_1,
                                       EntityIndex<I> entity_index_2) {
      return label(entity_index_1) < label(entity_index_2);
    });

    // Permute the entities.
    permute_base_(perm);
  }

  /// @brief Assign the entity labels.
  /// @param label Labels range to assign. May be smaller,
  ///              than the amount of entities.
  /// @warning This operation is very slow!
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
    std::ranges::stable_sort(perm, [&](EntityIndex<I> entity_index_1,
                                       EntityIndex<I> entity_index_2) {
      const auto new_label = [&](EntityIndex<I> entity_index) {
        const auto entity_index_sz = static_cast<size_t>(entity_index);
        return entity_index >= labels.size() ? Label{0} :
                                               labels[entity_index_sz];
      };
      return new_label(entity_index_1) < new_label(entity_index_2);
    });
    permute_base_<I>(perm);

    // Regenerate label ranges.
    const size_t delta = this->num_entities<I>() - std::ranges::size(labels);
    auto& entity_ranges = std::get<I>(entity_ranges_tuple_);
    entity_ranges.clear();
    std::ranges::for_each(labels, [&](Label label) {
      entity_ranges.resize(
          std::max(entity_ranges.size(), static_cast<size_t>(label) + 2));
      entity_ranges[label + 1] += 1;
    });
    entity_ranges[Label{0} + 1] += delta;
    std::partial_sum(
        entity_ranges.begin(), entity_ranges.end(), entity_ranges.begin(),
        [](EntityIndex<I> entity_index_1, EntityIndex<I> entity_index_2) {
          return entity_index_1 + static_cast<size_t>(entity_index_2);
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
    if (std::ranges::size(cell_face_nodes) == 2) {
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

  template<size_t I, std::ranges::random_access_range Range>
    requires std::permutable<std::ranges::iterator_t<Range>> &&
             std::same_as<std::ranges::range_value_t<Range>, EntityIndex<I>>
  constexpr void permute_base_(Range&& perm, meta::type<EntityIndex<I>> = {}) {
    { // Generate the inverse permutations and fix the connectivity.
      IndexedVector<EntityIndex<I>, EntityIndex<I>> iperm(num_entities<I>());
      permutations::invert_permutation(perm, iperm.begin());
      meta::for_each<EntityIndices_>([&]<size_t J>(meta::type<EntityIndex<J>>) {
        /// @todo Use parallel loop here!
        std::ranges::for_each(entities<J>(), [&](EntityIndex<J> entity_index) {
          using T = Table<EntityIndex<J>, EntityIndex<I>>;
          std::ranges::for_each(std::get<T>(connectivity_tuple_)[entity_index],
                                [&](EntityIndex<I>& adjacent_index) {
                                  adjacent_index = iperm[adjacent_index];
                                });
        });
      });
    }

    // Permute the entity-* connectivity.
    /// @todo This is actually the slowest part of the entire function.
    /// @todo For some table types, when the rows are swappable,
    ///       so the inplace permutation may be used.
    meta::for_each<EntityIndices_>([&]<size_t J>(meta::type<EntityIndex<J>>) {
      using T = Table<EntityIndex<I>, EntityIndex<J>>;
      auto& table = std::get<T>(connectivity_tuple_);
      T permuted_table{};
      permuted_table.reserve(num_entities<I>());
      std::ranges::for_each(entities<I>(), [&](EntityIndex<I> entity_index) {
        const auto entity_index_sz = static_cast<size_t>(entity_index);
        const EntityIndex<I> permuted_entity_index = perm[entity_index_sz];
        permuted_table.push_back(std::move(table[permuted_entity_index]));
      });
      table = std::move(permuted_table);
    });

    // Permute the entity shape properties.
    permutations::permute_inplace(perm, [&](EntityIndex<I> entity_index_1,
                                            EntityIndex<I> entity_index_2) {
      const auto gather_shape_properties_ = [&](EntityIndex<I> entity_index) {
        auto& positions = std::get<I>(entity_positions_tuple_);
        if constexpr (std::is_same_v<EntityIndex<I>, NodeIndex>) {
          return std::tie(labels[entity_index], positions[entity_index]);
        } else {
          auto& shape_types = std::get<I>(entity_shape_types_tuple_);
          auto& volumes = std::get<I>(entity_volumes_tuple_);
          if constexpr (std::is_same_v<EntityIndex<I>, FaceIndex>) {
            return std::tie(shape_types[entity_index], volumes[entity_index],
                            positions[entity_index],
                            face_normals_[entity_index]);
          } else {
            return std::tie(shape_types[entity_index], volumes[entity_index],
                            positions[entity_index]);
          }
        }
      };
      auto properties_tuple_1 = gather_shape_properties_(entity_index_1);
      auto properties_tuple_2 = gather_shape_properties_(entity_index_2);
      std::swap(properties_tuple_1, properties_tuple_2);
    });
  }

}; // class UnstructuredMesh

template<size_t Dim, size_t TopologicalDim, template<class, class> class Table>
inline constexpr bool
    enable_mesh_v<UnstructuredMesh<Dim, TopologicalDim, Table>> = true;

} // namespace Storm
