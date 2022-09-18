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

#include "ConvectionScheme.hh"
#include "Field.hh"
#include "FluxScheme.hh"
#include "GradientLimiterScheme.hh"
#include "SkunkFvBC.hh"

#include <map>

namespace Storm::Feathers {

/**
 * @brief A Finite volume solver.
 */
template<typename MhdPhysicsT>
class MhdFvSolverT :
    public std::enable_shared_from_this<MhdFvSolverT<MhdPhysicsT>> {
public:

  using MhdFluidStateT = typename MhdPhysicsT::MhdFluidStateT;
  static constexpr size_t num_vars = MhdPhysicsT::num_vars;

private:

  std::shared_ptr<Mesh> m_mesh;
  std::shared_ptr<iConvectionScheme> m_conv;
  std::map<Label, std::shared_ptr<MhdFvBcPT<MhdPhysicsT>>> m_bcs;

public:

  explicit MhdFvSolverT(std::shared_ptr<Mesh> mesh)
      : m_mesh(mesh),
        m_conv(new LinearUpwindConvectionScheme(
            *mesh, tLaxFriedrichsFluxScheme<MhdPhysicsT>{},
            LeastSquaresGradientScheme{*mesh},
            GradientLimiterScheme{*mesh, //
                                  CubicSlopeLimiter{}, CubicSecondLimiter{}})) {
    m_bcs[Label{1}] = std::make_shared<MhdFvBcFarFieldT<MhdPhysicsT>>();
    m_bcs[Label{2}] = std::make_shared<MhdFvBcSlipT<MhdPhysicsT>>();
  }

  /**
   * @brief Compute spacial discretization.
   */
  void calc_func(tScalarField& u, tScalarField& u_out) const {
    /*
     * Clear fields and apply boundary conditions.
     */
    ForEach(m_mesh->cells(),
            [&](CellView<Mesh> cell) { u_out[cell].fill(0.0); });
    for (size_t mark = 1; mark < m_mesh->num_face_labels(); ++mark) {
      const Label label{mark};
      const auto& bc = m_bcs.at(label);
      ForEach(m_mesh->faces(label), [&](FaceView<Mesh> face) {
        bc->get_ghost_state(
            face.normal3D(), //
            face.inner_cell().center3D(), face.outer_cell().center3D(),
            u[face.inner_cell()].data(), u[face.outer_cell()].data());
      });
    }

    m_conv->get_cell_convection(5, u_out, u);
  }

  /*
   * Compute time step.
   */
  void calc_step(real_t& dt, tScalarField& u, tScalarField& u_hat) const {
    calc_func(u, u_hat);
    ForEach(m_mesh->interior_cells(), [&](CellView<Mesh> cell) {
      for (size_t i = 0; i < num_vars; ++i) {
        u_hat[cell][i] = u[cell][i] - dt * u_hat[cell][i];
      }
    });
  }

}; // class MhdFvSolverT

} // namespace Storm::Feathers
