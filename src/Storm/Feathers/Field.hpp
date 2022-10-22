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

template<class T, size_t N = 1>
using Subfield = std::array<T, N>;

template<class Mesh, class T, size_t N = 1>
using CellField = IndexedVector<CellIndex<Mesh>, Subfield<T, N>>;
template<class Mesh, class T, size_t N = 1>
using CellVecField = IndexedVector<CellIndex<Mesh>, Subfield<vec2_t, N>>;
template<class Mesh, class T, size_t N = 1>
using CellMatField = IndexedVector<CellIndex<Mesh>, Subfield<mat2_t, N>>;

struct sFieldDesc {
  const char* name;
  size_t var_index;
  CellField<Mesh, real_t, 5>* scalar;
}; // struct sFieldDesc

} // namespace Storm::Feathers
