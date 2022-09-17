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

#include "ConvectionScheme.hh"

namespace Storm {

/**
 * Compute the first-order upwind convection.
 */
void cUpwindConvectionScheme::get_cell_convection(size_t num_vars,
                                                  tScalarField& div_f,
                                                  const tScalarField& u) const {
  /* Compute the first order numerical fluxes. */
  tScalarField flux_u(num_vars, m_mesh->faces().size());
  ForEach(face_views(*m_mesh), [&](FaceView face) {
    const CellView cell_outer = face.outer_cell();
    const CellView cell_inner = face.inner_cell();

    tScalarSubField flux = (flux_u[face] = {});
    m_flux->get_numerical_flux(num_vars, face.normal(), u[cell_outer],
                               u[cell_inner], flux);
  });

  /* Compute the first order convection. */
  ForEach(int_cell_views(*m_mesh), [&](CellView cell) {
    div_f[cell] = {};
    cell.for_each_face([&](FaceView face) {
      const CellView cell_outer = face.outer_cell();
      const CellView cell_inner = face.inner_cell();
      const real_t ds = face.area();
      if (cell_outer == cell) {
        for (size_t i = 0; i < num_vars; ++i) {
          div_f[cell][i] -= flux_u[face][i] * ds;
        }
      } else if (cell_inner == cell) {
        for (size_t i = 0; i < num_vars; ++i) {
          div_f[cell][i] += flux_u[face][i] * ds;
        }
      }
    });
    const real_t inv_dv = 1.0 / cell.volume();
    for (size_t i = 0; i < num_vars; ++i) {
      div_f[cell][i] *= inv_dv;
    }
  });

} // cUpwindConvectionScheme::get_cell_convection

/**
 * Compute the second-order upwind convection.
 */
void cUpwind2ConvectionScheme::get_cell_convection(
    size_t num_vars, tScalarField& div_f, const tScalarField& u) const {
  /* Compute the second order limited gradients. */
  tVectorField grad_u(num_vars, m_mesh->cells().size());
  m_gradient_scheme->get_gradients(num_vars, grad_u, u);

  tScalarField lim_u(num_vars, m_mesh->cells().size());
  m_gradient_limiter_scheme->get_cell_limiter(num_vars, lim_u, u, grad_u);

  /* Compute the second order numerical fluxes:
   * integrate the numerical flux over the face Nodes. */
  tScalarField flux_f(num_vars, m_mesh->faces().size());
  ForEach(face_views(*m_mesh), [&](FaceView face) {
    const CellView cell_outer = face.outer_cell();
    const CellView cell_inner = face.inner_cell();
    const vec3_t dr_outer = face.center() - cell_outer.center();
    const vec3_t dr_inner = face.center() - cell_inner.center();
    FEATHERS_TMP_SCALAR_FIELD(u_outer, num_vars);
    FEATHERS_TMP_SCALAR_FIELD(u_inner, num_vars);
    for (size_t i = 0; i < num_vars; ++i) {
      u_outer[i] =
          u[cell_outer][i] +
          lim_u[cell_outer][i] * glm::dot(grad_u[cell_outer][i], dr_outer);
      u_inner[i] =
          u[cell_inner][i] +
          lim_u[cell_inner][i] * glm::dot(grad_u[cell_inner][i], dr_inner);
    }

    tScalarSubField flux = (flux_f[face] = {});
    m_flux->get_numerical_flux(num_vars, face.normal(), u_outer, u_inner, flux);
  });

  /* Compute the second order convection. */
  ForEach(int_cell_views(*m_mesh), [&](CellView cell) {
    div_f[cell] = {};
    cell.for_each_face([&](FaceView face) {
      const CellView cell_outer = face.outer_cell();
      const CellView cell_inner = face.inner_cell();
      const real_t ds = face.area();
      if (cell_outer == cell) {
        for (size_t i = 0; i < num_vars; ++i) {
          div_f[cell][i] -= flux_f[face][i] * ds;
        }
      } else if (cell_inner == cell) {
        for (size_t i = 0; i < num_vars; ++i) {
          div_f[cell][i] += flux_f[face][i] * ds;
        }
      }
    });
    const real_t inv_dv = 1.0 / cell.volume();
    for (size_t i = 0; i < num_vars; ++i) {
      div_f[cell][i] *= inv_dv;
    }
  });
} // cUpwind2ConvectionScheme::get_cell_convection

} // namespace Storm
