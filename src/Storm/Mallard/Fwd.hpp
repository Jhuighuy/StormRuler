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

#include <Storm/Utils/Crtp.hpp>
#include <Storm/Utils/Index.hpp>
#include <Storm/Utils/Meta.hpp>

#include <compare>
#include <concepts>
#include <ranges>

namespace Storm {

struct LabelTag;

/// @brief Label index type.
using Label = Index<LabelTag>;

template<crtp_derived Derived>
class MeshInterface;

/// @brief Types, enabled to be a mesh.
template<class Mesh>
inline constexpr bool enable_mesh_v =
    derived_from_crtp_interface<Mesh, MeshInterface>;

/// @brief Mesh concept.
template<class Mesh>
concept mesh = enable_mesh_v<std::remove_const_t<Mesh>>;

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
template<mesh Mesh>
using FaceIndex = typename Mesh::FaceIndex;

/// @brief Cell index type.
template<mesh Mesh>
using CellIndex = typename Mesh::CellIndex;

/// @brief Mesh spatial vector type.
template<mesh Mesh>
using mesh_vec_t = std::remove_cvref_t< //
    decltype(std::declval<Mesh>().position(std::declval<NodeIndex>()))>;

/// @brief Mesh spatial dimensionality.
template<mesh Mesh>
inline constexpr size_t mesh_dim_v = fast_vector_size_v<mesh_vec_t<Mesh>>;

} // namespace Storm
