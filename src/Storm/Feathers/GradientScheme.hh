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

#include "Field.hh"

namespace Storm::Feathers {

/// @brief Weighted Least-Squares gradient estimation scheme, cell-based:
/// computes cell-centered gradients based on the cell-centered values.
template<mesh Mesh>
class LeastSquaresGradientScheme final {
private:

  const Mesh* p_mesh_;
  bool initialized_ = false;
  tMatrixField m_inverse_matrices;

public:

  /// @brief Construct the gradient scheme.
  constexpr explicit LeastSquaresGradientScheme(const Mesh& mesh)
      : p_mesh_{&mesh}, m_inverse_matrices(1, p_mesh_->num_cells()) {}

private:

  void init_gradients_() {
    /* Compute the least-squares
     * problem matrices for the interior Cells. */
    std::ranges::for_each(p_mesh_->interior_cells(), [&](CellView<Mesh> cell) {
      mat3_t& mat = (m_inverse_matrices[cell][0] = mat3_t(0.0));
      cell.for_each_face_cells(
          [&](CellView<Mesh> cell_inner, CellView<Mesh> cell_outer) {
            const vec3_t dr = cell_outer.center3D() - cell_inner.center3D();
            mat += glm::outerProduct(dr, dr);
          });
    });

    /* Compute the least squares problem right-hand statements for the boundary
     * cells. Use the same stencil as for the interior cell, but centered to a
     * boundary cell. */
    std::ranges::for_each(p_mesh_->boundary_faces(), [&](FaceView<Mesh> face) {
      const CellView<Mesh> cell_inner = face.inner_cell();
      const CellView<Mesh> cell_outer = face.inner_cell();
      mat3_t& mat = (m_inverse_matrices[cell_outer][0] = mat3_t(0.0));
      const vec3_t dr = cell_outer.center3D() - cell_inner.center3D();
      mat += glm::outerProduct(dr, dr);
      cell_inner.for_each_face_cells([&](CellView<Mesh> cell_inner_inner,
                                         CellView<Mesh> cell_inner_outer) {
        if (cell_inner_outer == cell_inner) {
          std::swap(cell_inner_inner, cell_inner_outer);
        }
        const vec3_t dr_inner =
            cell_inner_outer.center3D() - cell_inner.center3D();
        mat += glm::outerProduct(dr_inner, dr_inner);
      });
    });

    /* Compute the inverse of the least squares problem matrices.
     * ( Matrix is stabilized by a small number, added to the diagonal. ) */
    std::ranges::for_each(p_mesh_->cells(), [&](CellView<Mesh> cell) {
      static const mat3_t eps(1e-14);
      mat3_t& mat = m_inverse_matrices[cell][0];
      mat = glm::inverse(mat + eps);
    });

    initialized_ = true;
  }

public:

  /** Compute cell-centered gradients. */
  void get_gradients(size_t num_vars, tVectorField& grad_u,
                     const tScalarField& u) const {
    if (!initialized_)
      const_cast<LeastSquaresGradientScheme*>(this)->init_gradients_();

    /* Compute the least-squares
     * problem right-hand statements for the interior Cells. */
    std::ranges::for_each(p_mesh_->interior_cells(), [&](CellView<Mesh> cell) {
      grad_u[cell].fill(vec3_t(0.0));
      cell.for_each_face_cells(
          [&](CellView<Mesh> cell_inner, CellView<Mesh> cell_outer) {
            const vec3_t dr = cell_outer.center3D() - cell_inner.center3D();
            for (size_t i = 0; i < num_vars; ++i) {
              grad_u[cell][i] += (u[cell_outer][i] - u[cell_inner][i]) * dr;
            }
          });
    });

    /* Compute the least squares problem right-hand statements for the boundary
     * cells. Use the same stencil as for the interior cell, but centered to a
     * boundary cell. */
    std::ranges::for_each(p_mesh_->boundary_faces(), [&](FaceView<Mesh> face) {
      const CellView<Mesh> cell_inner = face.inner_cell();
      const CellView<Mesh> cell_outer = face.inner_cell();
      grad_u[cell_outer].fill(vec3_t(0.0));
      const vec3_t dr = cell_outer.center3D() - cell_inner.center3D();
      for (size_t i = 0; i < num_vars; ++i) {
        grad_u[cell_outer][i] += (u[cell_outer][i] - u[cell_inner][i]) * dr;
      }
      cell_inner.for_each_face_cells([&](CellView<Mesh> cell_inner_inner,
                                         CellView<Mesh> cell_inner_outer) {
        if (cell_inner_outer == cell_inner) {
          std::swap(cell_inner_inner, cell_inner_outer);
        }
        const vec3_t dr_inner =
            cell_inner_outer.center3D() - cell_inner.center3D();
        for (size_t i = 0; i < num_vars; ++i) {
          grad_u[cell_outer][i] +=
              (u[cell_inner_outer][i] - u[cell_inner][i]) * dr_inner;
        }
      });
    });

    /* Solve the least-squares problem. */
    std::ranges::for_each(p_mesh_->cells(), [&](CellView<Mesh> cell) {
      for (size_t i = 0; i < num_vars; ++i) {
        const mat3_t& mat = m_inverse_matrices[cell][0];
        grad_u[cell][i] = mat * grad_u[cell][i];
      }
    });
  }

}; // class LeastSquaresGradientScheme

} // namespace Storm::Feathers
