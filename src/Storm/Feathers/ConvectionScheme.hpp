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

#include "FluxScheme.hpp"
#include "GradientLimiterScheme.hpp"
#include "GradientScheme.hpp"
#include "SkunkFvBC.hpp"

namespace Storm::Feathers {

/**
 * Abstract convection scheme.
 */
class iConvectionScheme : public tObject<iConvectionScheme> {
public:

  std::map<Label, std::shared_ptr<MhdFvBcPT<tGasPhysics>>>* bcs_;

  /** Compute the nonlinear convection. */
  virtual void
  get_cell_convection(CellField<Mesh, real_t, 5>& div_f,
                      const CellField<Mesh, real_t, 5>& u) const = 0;

}; // class iConvectionScheme

/// @brief Piecewise-constant upwind convection scheme.
template<mesh Mesh, class FluxScheme>
  requires std::is_object_v<FluxScheme>
class UpwindConvectionScheme final : public iConvectionScheme {
private:

  const Mesh* p_mesh_;
  STORM_NO_UNIQUE_ADDRESS_ FluxScheme flux_scheme_;

public:

  /// @brief Construct the first order upwind convection scheme.
  constexpr explicit UpwindConvectionScheme(const Mesh& mesh,
                                            FluxScheme flux_scheme) noexcept
      : p_mesh_{&mesh}, flux_scheme_{std::move(flux_scheme)} {}

  void get_cell_convection(CellField<Mesh, real_t, 5>& div_f,
                           const CellField<Mesh, real_t, 5>& u) const final {
    (*this)(div_f, u);
  }

  /// @brief Compute the first order upwind nonlinear convection.
  template<class Real, size_t NumVars>
  void operator()(CellField<Mesh, Real, NumVars>& div_f,
                  const CellField<Mesh, Real, NumVars>& u) const noexcept {
    // Compute the fluxes for the interior faces.
    std::ranges::for_each(p_mesh_->interior_faces(), [&](FaceView<Mesh> face) {
      const CellView<Mesh> cell_inner = face.inner_cell();
      const CellView<Mesh> cell_outer = face.outer_cell();

      // Compute the flux.
      const auto flux =
          flux_scheme_(face.normal(), u[cell_outer], u[cell_inner]);
      div_f[cell_inner] += flux * face.area() / cell_inner.volume();
      div_f[cell_outer] -= flux * face.area() / cell_outer.volume();
    });

    // Compute the fluxes for the boundary faces.
    for (const auto& [label, bc] : *bcs_) {
      std::ranges::for_each(p_mesh_->faces(label), [&](FaceView<Mesh> face) {
        const CellView<Mesh> cell_inner = face.inner_cell();
        Subfield<Real, NumVars> u_outer{};
        bc->get_ghost_state(face.normal(), face.center(),
                            u[face.inner_cell()].data(), u_outer.data());

        // Compute the flux.
        const auto flux = flux_scheme_(face.normal(), u_outer, u[cell_inner]);
        div_f[cell_inner] += flux * face.area() / cell_inner.volume();
      });
    }
  }

}; // class UpwindConvectionScheme

/// @brief Piecewise-linear upwind convection scheme.
template<mesh Mesh, class FluxScheme, //
         class GradientScheme, class GradientLimiterScheme>
  requires std::is_object_v<FluxScheme> && std::is_object_v<GradientScheme> &&
           std::is_object_v<GradientLimiterScheme>
class LinearUpwindConvectionScheme final : public iConvectionScheme {
public:

  const Mesh* p_mesh_;
  STORM_NO_UNIQUE_ADDRESS_ FluxScheme flux_scheme_;
  STORM_NO_UNIQUE_ADDRESS_ GradientScheme gradient_scheme_;
  STORM_NO_UNIQUE_ADDRESS_ GradientLimiterScheme gradient_limiter_scheme_;

public:

  /// @brief Construct the second order upwind convection scheme.
  constexpr explicit LinearUpwindConvectionScheme(
      const Mesh& mesh, FluxScheme flux_scheme, GradientScheme gradient_scheme,
      GradientLimiterScheme gradient_limiter_scheme) noexcept
      : p_mesh_{&mesh}, flux_scheme_{flux_scheme}, //
        gradient_scheme_{std::move(gradient_scheme)},
        gradient_limiter_scheme_{std::move(gradient_limiter_scheme)} {}

  void get_cell_convection(CellField<Mesh, real_t, 5>& div_f,
                           const CellField<Mesh, real_t, 5>& u) const final {
    (*this)(div_f, u);
  }

  /// @brief Compute the second-order upwind nonlinear convection.
  template<class Real, size_t NumVars>
  void operator()(CellField<Mesh, Real, NumVars>& div_f,
                  const CellField<Mesh, Real, NumVars>& u) const noexcept {
    // Compute the gradients.
    CellVecField<Mesh, Real, NumVars> grad_u(p_mesh_->num_cells());
    CellField<Mesh, Real, NumVars> lim_u(p_mesh_->num_cells());
    gradient_scheme_(grad_u, u);
    gradient_limiter_scheme_(lim_u, u, grad_u);
    std::ranges::for_each(p_mesh_->interior_cells(), [&](CellView<Mesh> cell) {
      grad_u[cell] *= lim_u[cell];
    });

    // Compute the fluxes for the interior faces.
    std::ranges::for_each(p_mesh_->interior_faces(), [&](FaceView<Mesh> face) {
      const CellView<Mesh> cell_inner = face.inner_cell();
      const CellView<Mesh> cell_outer = face.outer_cell();

      // Reconstruct the face values.
      const auto dr_outer = face.center() - cell_outer.center();
      const auto dr_inner = face.center() - cell_inner.center();
      Subfield<Real, NumVars> u_outer{}, u_inner{};
      for (size_t i = 0; i < NumVars; ++i) {
        u_inner[i] =
            u[cell_inner][i] + glm::dot(grad_u[cell_inner][i], dr_inner);
        u_outer[i] =
            u[cell_outer][i] + glm::dot(grad_u[cell_outer][i], dr_outer);
      }

      // Compute the flux.
      const auto flux = flux_scheme_(face.normal(), u_outer, u_inner);
      div_f[cell_inner] += flux * face.area() / cell_inner.volume();
      div_f[cell_outer] -= flux * face.area() / cell_outer.volume();
    });

    // Compute the fluxes for the boundary faces.
    for (const auto& [label, bc] : *bcs_) {
      std::ranges::for_each(p_mesh_->faces(label), [&](FaceView<Mesh> face) {
        const CellView<Mesh> cell_inner = face.inner_cell();

        // Reconstruct the face values.
        const auto dr_inner = face.center() - cell_inner.center();
        Subfield<Real, NumVars> u_outer{}, u_inner{};
        for (size_t i = 0; i < NumVars; ++i) {
          u_inner[i] = u[cell_inner][i] +
                       glm::dot(0.0 * grad_u[cell_inner][i], dr_inner);
        }
        bc->get_ghost_state(face.normal(), face.center(), //
                            u_inner.data(), u_outer.data());

        // Compute the flux.
        const auto flux = flux_scheme_(face.normal(), u_outer, u_inner);
        div_f[cell_inner] += flux * face.area() / cell_inner.volume();
      });
    }
  }

}; // class UpwindConvectionScheme

} // namespace Storm::Feathers
