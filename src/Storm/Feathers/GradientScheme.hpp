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

#include "Field.hpp"

namespace Storm::Feathers {

/// @brief Weighted Least-Squares gradient estimation scheme, cell-based:
/// computes cell-centered gradients based on the cell-centered values.
template<mesh Mesh>
class LeastSquaresGradientScheme final {
private:

  const Mesh* p_mesh_;
  CellMatField<Mesh, real_t> m_inverse_matrices;

public:

  /// @brief Construct the gradient scheme.
  constexpr explicit LeastSquaresGradientScheme(const Mesh& mesh)
      : p_mesh_{&mesh}, m_inverse_matrices(p_mesh_->num_cells()) {
    // Compute the least-squares problem matrices.
    std::ranges::for_each(p_mesh_->interior_cells(), [&](CellView<Mesh> cell) {
      auto& mat = m_inverse_matrices[cell][0];
      mat = {};
      cell.for_each_face_cells(
          [&](CellView<Mesh> cell_inner, CellView<Mesh> cell_outer) {
            const vec3_t dr = cell_outer.center3D() - cell_inner.center3D();
            mat += glm::outerProduct(dr, dr);
          });
    });

    // Compute the inverse of the least squares problem matrices.
    // (Matrix is stabilized by a small number, added to the diagonal.)
    std::ranges::for_each(p_mesh_->interior_cells(), [&](CellView<Mesh> cell) {
      static const mat3_t eps(1e-14);
      auto& mat = m_inverse_matrices[cell][0];
      mat = glm::inverse(mat + eps);
    });
  }

  /// @brief Compute cell-centered gradients.
  template<class Real, size_t NumVars>
  void operator()(CellVecField<Mesh, Real, NumVars>& grad_u,
                  const CellField<Mesh, Real, NumVars>& u) const noexcept {
    // Compute the least-squares problem right-hand statements.
    std::ranges::for_each(p_mesh_->interior_cells(), [&](CellView<Mesh> cell) {
      grad_u[cell].fill(vec3_t(0.0));
      cell.for_each_face_cells(
          [&](CellView<Mesh> cell_inner, CellView<Mesh> cell_outer) {
            const vec3_t dr = cell_outer.center3D() - cell_inner.center3D();
            for (size_t i = 0; i < NumVars; ++i) {
              grad_u[cell][i] += (u[cell_outer][i] - u[cell_inner][i]) * dr;
            }
          });
    });

    // Solve the least-squares problem.
    std::ranges::for_each(p_mesh_->cells(), [&](CellView<Mesh> cell) {
      for (size_t i = 0; i < NumVars; ++i) {
        const mat3_t& mat = m_inverse_matrices[cell][0];
        grad_u[cell][i] = mat * grad_u[cell][i];
      }
    });
  }

}; // class LeastSquaresGradientScheme

} // namespace Storm::Feathers
