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

namespace Storm {

struct LabelTag;

/// @brief Label index type.
using Label = Index<LabelTag>;

template<size_t I>
struct TopologicalIndexTag;

/// @brief Topological index type.
template<size_t I>
using EntityIndex = Index<TopologicalIndexTag<I>>;

/// @brief Node index type.
using NodeIndex = EntityIndex<0>;

/// @brief Edge index type.
using EdgeIndex = EntityIndex<1>;

/// @brief Face index type.
template<class Mesh>
using FaceIndex = typename Mesh::FaceIndex;

/// @brief Cell index type.
template<class Mesh>
using CellIndex = typename Mesh::CellIndex;

/// @brief Types, enabled to be a mesh.
template<class>
inline constexpr bool enable_mesh_v = false;

/// @brief Mesh concept.
// clang-format off
template<class Mesh>
concept mesh = enable_mesh_v<std::remove_const_t<Mesh>> &&
    requires {
      typename FaceIndex<Mesh>;
      typename CellIndex<Mesh>;
    };
// clang-format on

/// @brief Mesh spatial vector type.
template<mesh Mesh>
using mesh_vec_t = std::remove_cvref_t< //
    decltype(std::declval<Mesh>().position(std::declval<NodeIndex>()))>;

/// @brief Mesh spatial dimensionality.
template<mesh Mesh>
inline constexpr size_t mesh_dim_v = fast_vector_size_v<mesh_vec_t<Mesh>>;

namespace detail_ {
  inline constexpr size_t face_inner_cell_ = 0;
  inline constexpr size_t face_outer_cell_ = 1;
} // namespace detail_

} // namespace Storm
