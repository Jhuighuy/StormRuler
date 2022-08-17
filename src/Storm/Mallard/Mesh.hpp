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

#include <fstream>
#include <ranges>
#include <tuple>

namespace Storm {

enum class ShapeType {};

namespace detail_ {
  struct label_tag_t_;
  template<size_t I>
  struct topological_index_tag_t_;
} // namespace detail_

/// @brief Label index type.
using Label = Index<detail_::label_tag_t_>;

/// @brief Topological index type.
template<size_t I>
using EntityIndex = Index<detail_::topological_index_tag_t_<I>>;

/// @brief Node index type.
using NodeIndex = EntityIndex<0>;

/// @brief Edge index type.
using EdgeIndex = EntityIndex<1>;

struct Triangle {
  std::vector<NodeIndex> Indices;

  const auto& node_indices() const {
    return Indices;
  }
  real_t compute_volume(auto&&) const {
    return {};
  }
  glm::dvec3 compute_center_position(auto&&) {
    return {};
  }
};

/// @brief Hybrid unstructured multidimensional mesh.
/// @tparam SpDim Spatial dimensionality.
/// @tparam TpDim Topological dimensionality.
/// @tparam Table ??
template<size_t SpDim, size_t TpDim, //
         template<class, class> class Table = VoidCsrTable>
class BaseBaseBaseMesh {
public:

  /// @brief Spatial vector type.
  using SpVec = FastVector<real_t, SpDim>;

  /// @brief Node index type.
  using NodeIndex = EntityIndex<0>;

  /// @brief Edge index type.
  using EdgeIndex = EntityIndex<1>;

  /// @brief Face index type.
  using FaceIndex = EntityIndex<TpDim - 1>;

  /// @brief Cell index type.
  using CellIndex = EntityIndex<TpDim>;

private:

  using EntityIndices_ = meta::make_list_t<EntityIndex, TpDim + 1>;

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
        meta::pair_list_t<real_t, meta::rest_t<EntityIndices_>>
      >
    >
  > entity_volumes_tuple_{};

  meta::as_std_tuple_t<
    meta::transform_t<
      meta::pair_cast_fn<IndexedVector>,
      meta::pair_list_t<SpVec, EntityIndices_>
    >
  > entity_positions_tuple_{};

  // You don't belong here, huh?
  IndexedVector<SpVec, FaceIndex> face_normals_{};

  meta::as_std_tuple_t<
    meta::transform_t<
      meta::pair_cast_fn<Table>,
      meta::cartesian_product_t<EntityIndices_>
    >
  > entity_connectivity_hypertuple_{};

  // clang-format on

public:

  /// @brief Initialize the mesh.
  constexpr BaseBaseBaseMesh() {
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
  [[nodiscard]] constexpr SpVec //
  position(EntityIndex<I> index) const noexcept {
    return std::get<I>(entity_positions_tuple_)[index];
  }

  /// @brief Normal to the face at @p face_index.
  [[nodiscard]] constexpr SpVec //
  normal(FaceIndex face_index) const noexcept {
    return face_normals_[face_index];
  }

  /// @brief Range of adjacent entity indices of dim J of an entity at @p index.
  template<size_t J, size_t I>
  [[nodiscard]] constexpr auto //
  adjacent(EntityIndex<I> index, meta::type<EntityIndex<J>> = {}) noexcept {
    using T = Table<EntityIndex<I>, EntityIndex<J>>;
    return std::get<T>(entity_connectivity_hypertuple_)[index];
  }

  template<size_t I, class Data>
  constexpr EntityIndex<I> insert(Data&& data,
                                  meta::type<EntityIndex<I>> = {}) {
    const EntityIndex<I> index{num_entities<I>()};

    // Assign the label.
    std::get<I>(label_ranges_tuple_).back() += 1;
    std::get<I>(entity_labels_tuple_).emplace_back(num_labels<I>() - 1);

    // Assign the geometrical properties.
    if constexpr (std::is_same_v<EntityIndex<I>, NodeIndex>) {
      static_assert(std::constructible_from<SpVec, Data&&>,
                    "Invalid node position type.");
      std::get<I>(entity_positions_tuple_)
          .emplace_back(std::forward<Data>(data));
    } else {
      static_assert(true, //
                    "Here should be some check...");
      std::get<I>(entity_volumes_tuple_)
          .emplace_back(data.compute_volume(*this));
      std::get<I>(entity_positions_tuple_)
          .emplace_back(data.compute_center_position(*this));
      if constexpr (std::is_same_v<EntityIndex<I>, FaceIndex>) {
        face_normals_.emplace_back(data.compute_normal(*this));
      }
    }

    // Assign the topological properties.
    /// @todo This should be implemented correctly.
    if constexpr (!std::is_same_v<EntityIndex<I>, NodeIndex>) {
      using T = Table<EntityIndex<I>, EntityIndex<0>>;
      std::get<T>(entity_connectivity_hypertuple_).emplace(data.node_indices());
    }

    return index;
  }

  void read_from_triangle(std::string const& path) {
    std::string line;

    std::ifstream nodeStream(path + std::string("node"));
    STORM_ENSURE_(nodeStream.is_open(), "");
    size_t numNodes{0}, dim{0};
    nodeStream >> numNodes >> dim;
    std::getline(nodeStream, line);
    for (size_t i{0}; i < numNodes; ++i) {
      NodeIndex nodeIndex{0};
      SpVec nodeCoords(0.0);
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
      CellIndex cellIndex{0};
      std::vector<NodeIndex> cellNodes(3);
      cellStream >> cellIndex >> cellNodes[0] >> cellNodes[1] >> cellNodes[2];
      insert<2>(Triangle{cellNodes});
      // STORM_ENSURE_(cellIndex == insert<2>(Triangle{cellNodes}), "");
      std::getline(cellStream, line);
    }
#endif
  }

}; // class BaseBaseBaseMesh

} // namespace Storm
