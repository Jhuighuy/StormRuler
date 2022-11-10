///  ______  ______   ______   ______  __  __   ______   ______   ______
/// /\  ___\/\  ___\ /\  __ \ /\__  _\/\ \_\ \ /\  ___\ /\  __ \ /\  ___\ 
/// \ \  __\\ \  _\  \ \  __ \\/_/\ \/\ \  __ \\ \  __\ \ \  __/ \ \___  \ 
///  \ \_\   \ \_____\\ \_\ \_\  \ \_\ \ \_\ \_\\ \_____\\ \_\ \_\\/\_____\ 
///   \/_/    \/_____/ \/_/\/_/   \/_/  \/_/\/_/ \/_____/ \/_/ /_/ \/_____/
///
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

#include <Storm/Bittern/Mat.hpp>

#include <Storm/Mallard/MeshUnstructured.hpp>

#include <glm/glm.hpp>

#include <memory>
#include <ranges>

namespace Storm::Feathers {

using vec2_t = glm::dvec2;
using vec3_t = glm::dvec3;
using mat2_t = glm::dmat2;
using mat3_t = glm::dmat3;
template<typename type_t>
using tObject = std::enable_shared_from_this<type_t>;
using Mesh = UnstructuredMesh<2, 2, CsrTable>;

/// @brief Mesh field entry.
template<class Value, size_t NumVars = 1>
using Subfield = std::conditional_t<NumVars == 1, //
                                    Value, StaticVector<Value, NumVars>>;

/// @brief Generic mesh field.
template<mesh Mesh, index Index, class Value = real_t, size_t NumVars = 1>
class Field final {
private:

  std::string name_;
  const Mesh* p_mesh_ = nullptr;
  IndexedVector<Index, Subfield<Value, NumVars>> data_;

public:

  constexpr Field() = default;

  constexpr Field(const Mesh& mesh)
      : p_mesh_{&mesh}, data_(p_mesh_->num_entities(meta::type_v<Index>)) {}

  /// @brief Field shape.
  [[nodiscard]] constexpr auto shape() const noexcept {
    return MatrixShape{p_mesh_->num_entities(meta::type_v<Index>), NumVars};
  }

  /// @todo Document me!
  constexpr void assign(const Field& other, bool copy = true) {
    *this = Field{*other.p_mesh_};
  }

  /// @brief Element at @p index.
  /// @{
  [[nodiscard]] constexpr auto& operator[](Index index) noexcept {
    STORM_ASSERT_(index < p_mesh_->num_entities(meta::type_v<Index>),
                  "Index is out of range!");
    return data_[index];
  }
  [[nodiscard]] constexpr const auto& operator[](Index index) const noexcept {
    STORM_ASSERT_(index < p_mesh_->num_entities(meta::type_v<Index>),
                  "Index is out of range!");
    return data_[index];
  }
  /// @}

  /// @brief Element at @p row_index and @p col_index.
  /// @{
  [[nodiscard]] constexpr Value& //
  operator()(size_t row_index, size_t col_index = 0) noexcept {
    return data_[Index{row_index}]; //[col_index];
  }
  [[nodiscard]] constexpr const Value&
  operator()(size_t row_index, size_t col_index = 0) const noexcept {
    return data_[Index{row_index}]; //[col_index];
  }
  /// @}

}; // class Field

/// @brief Node-centered field.
template<mesh Mesh, class Value, size_t NumVars = 1>
using NodeField = Field<Mesh, NodeIndex, Value, NumVars>;
/// @brief Node-centered vector field.
template<mesh Mesh, class Real, size_t NumVars = 1>
using NodeVectorField = NodeField<Mesh, mesh_vec_t<Mesh, Real>, NumVars>;
/// @brief Node-centered matrix field.
template<mesh Mesh, class Real, size_t NumVars = 1>
using NodeMatrixField = NodeField<Mesh, mesh_mat_t<Mesh, Real>, NumVars>;

/// @brief Edge-centered field.
template<mesh Mesh, class Value, size_t NumVars = 1>
using EdgeField = Field<Mesh, EdgeIndex, Value, NumVars>;
/// @brief Edge-centered vector field.
template<mesh Mesh, class Real, size_t NumVars = 1>
using EdgeVector = EdgeField<Mesh, mesh_vec_t<Mesh, Real>, NumVars>;
/// @brief Edge-centered matrix field.
template<mesh Mesh, class Real, size_t NumVars = 1>
using EdgeMatrixField = EdgeField<Mesh, mesh_mat_t<Mesh, Real>, NumVars>;

/// @brief Face-centered field.
template<mesh Mesh, class Value, size_t NumVars = 1>
using FaceField = Field<Mesh, CellIndex<Mesh>, Value, NumVars>;
/// @brief Face-centered vector field.
template<mesh Mesh, class Real, size_t NumVars = 1>
using FaceVectorField = FaceField<Mesh, mesh_vec_t<Mesh, Real>, NumVars>;
/// @brief Face-centered matrix field.
template<mesh Mesh, class Real, size_t NumVars = 1>
using FaceMatrixField = FaceField<Mesh, mesh_mat_t<Mesh, Real>, NumVars>;

/// @brief Cell-centered field.
template<mesh Mesh, class Value, size_t NumVars = 1>
using CellField = Field<Mesh, CellIndex<Mesh>, Value, NumVars>;
/// @brief Cell-centered vector field.
template<mesh Mesh, class Real, size_t NumVars = 1>
using CellVectorField = CellField<Mesh, mesh_vec_t<Mesh, Real>, NumVars>;
/// @brief Cell-centered matrix field.
template<mesh Mesh, class Real, size_t NumVars = 1>
using CellMatrixField = CellField<Mesh, mesh_mat_t<Mesh, Real>, NumVars>;

template<size_t N = 5>
struct sFieldDesc {
  const char* name;
  size_t var_index;
  CellField<Mesh, real_t, N>* scalar;
}; // struct sFieldDesc

} // namespace Storm::Feathers
