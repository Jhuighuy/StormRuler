/*
 *  ______  ______   ______   ______  __  __   ______   ______   ______
 * /\  ___\/\  ___\ /\  __ \ /\__  _\/\ \_\ \ /\  ___\ /\  __ \ /\  ___\
 * \ \  __\\ \  _\  \ \  __ \\/_/\ \/\ \  __ \\ \  __\ \ \  __/ \ \___  \
 *  \ \_\   \ \_____\\ \_\ \_\  \ \_\ \ \_\ \_\\ \_____\\ \_\ \_\\/\_____\
 *   \/_/    \/_____/ \/_/\/_/   \/_/  \/_/\/_/ \/_____/ \/_/ /_/ \/_____/
 *
 * Copyright (c) 2021 Oleg Butakov
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "GradientScheme.hh"

namespace Storm {

/**
 * Init the gradient scheme.
 */
void cLeastSquaresGradientScheme::init_gradients_() {
  /* Compute the least-squares
   * problem matrices for the interior Cells. */
  ForEach(int_cell_views(*m_mesh), [&](CellView cell) {
    mat3_t& mat = (m_inverse_matrices[cell][0] = mat3_t(0.0));
    cell.for_each_face_cells([&](CellView cell_inner, CellView cell_outer) {
      const vec3_t dr = cell_outer.center() - cell_inner.center();
      mat += glm::outerProduct(dr, dr);
    });
  });

  /* Compute the least squares problem
   * right-hand statements for the boundary Cells.
   * Use the same stencil as for the interior cell, but centered to a boundary
   * cell. */
  for_each_bnd_face_cells(*m_mesh, [&](CellView cell_inner,
                                       CellView cell_outer) {
    mat3_t& mat = (m_inverse_matrices[cell_outer][0] = mat3_t(0.0));
    const vec3_t dr = cell_outer.center() - cell_inner.center();
    mat += glm::outerProduct(dr, dr);
    cell_inner.for_each_face_cells([&](CellView cell_inner_inner,
                                       CellView cell_inner_outer) {
      if (cell_inner_outer == cell_inner) {
        std::swap(cell_inner_inner, cell_inner_outer);
      }
      const vec3_t dr_inner = cell_inner_outer.center() - cell_inner.center();
      mat += glm::outerProduct(dr_inner, dr_inner);
    });
  });

  /* Compute the inverse of the least squares problem matrices.
   * ( Matrix is stabilized by a small number, added to the diagonal. ) */
  ForEach(cell_views(*m_mesh), [&](CellView cell) {
    static const mat3_t eps(1e-14);
    mat3_t& mat = m_inverse_matrices[cell][0];
    mat = glm::inverse(mat + eps);
  });

} // cLeastSquaresGradientScheme::init_gradients_

/**
 * Compute cell-centered gradients using the Weighted Least-Squares, cell-based
 * version.
 */
void cLeastSquaresGradientScheme::get_gradients(size_t num_vars,
                                                tVectorField& grad_u,
                                                const tScalarField& u) const {
  /* Compute the least-squares
   * problem right-hand statements for the interior Cells. */
  ForEach(int_cell_views(*m_mesh), [&](CellView cell) {
    grad_u[cell].fill(vec3_t(0.0));
    cell.for_each_face_cells([&](CellView cell_inner, CellView cell_outer) {
      const vec3_t dr = cell_outer.center() - cell_inner.center();
      for (size_t i = 0; i < num_vars; ++i) {
        grad_u[cell][i] += (u[cell_outer][i] - u[cell_inner][i]) * dr;
      }
    });
  });

  /* Compute the least squares problem
   * right-hand statements for the boundary Cells.
   * Use the same stencil as for the interior cell, but centered to a boundary
   * cell. */
  for_each_bnd_face_cells(*m_mesh, [&](CellView cell_inner,
                                       CellView cell_outer) {
    grad_u[cell_outer].fill(vec3_t(0.0));
    const vec3_t dr = cell_outer.center() - cell_inner.center();
    for (ptrdiff_t i = 0; i < num_vars; ++i) {
      grad_u[cell_outer][i] += (u[cell_outer][i] - u[cell_inner][i]) * dr;
    }
    cell_inner.for_each_face_cells([&](CellView cell_inner_inner,
                                       CellView cell_inner_outer) {
      if (cell_inner_outer == cell_inner) {
        std::swap(cell_inner_inner, cell_inner_outer);
      }
      const vec3_t dr_inner = cell_inner_outer.center() - cell_inner.center();
      for (size_t i = 0; i < num_vars; ++i) {
        grad_u[cell_outer][i] +=
            (u[cell_inner_outer][i] - u[cell_inner][i]) * dr_inner;
      }
    });
  });

  /* Solve the least-squares problem. */
  ForEach(cell_views(*m_mesh), [&](CellView cell) {
    for (size_t i = 0; i < num_vars; ++i) {
      const mat3_t& mat = m_inverse_matrices[cell][0];
      grad_u[cell][i] = mat * grad_u[cell][i];
    }
  });
} // cLeastSquaresGradientScheme::get_gradients

} // namespace Storm
