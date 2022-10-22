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
  virtual void get_cell_convection(Field<real_t, 5>& div_f,
                                   const Field<real_t, 5>& u) const = 0;

}; // class iConvectionScheme

/// @brief Piecewise-constant upwind convection scheme.
/// This is a first order scheme.
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

  /// @brief Compute the first order upwind nonlinear convection.
  template<class Real, size_t NumVars>
  void get_cell_convection(Field<Real, NumVars>& div_f,
                           const Field<Real, NumVars>& u) const {
    // Compute the fluxes for the interior faces.
    std::ranges::for_each(p_mesh_->interior_faces(), [&](FaceView<Mesh> face) {
      const CellView<Mesh> cell_inner = face.inner_cell();
      const CellView<Mesh> cell_outer = face.outer_cell();

      // Compute the flux.
      Subfield<Real, NumVars> flux{};
      flux_scheme_.get_numerical_flux(face.normal3D(), //
                                      u[cell_outer], u[cell_inner], flux);
      const real_t ds = face.area();
      for (size_t i = 0; i < NumVars; ++i) {
        div_f[cell_inner][i] += flux[i] * ds / cell_inner.volume();
        div_f[cell_outer][i] -= flux[i] * ds / cell_outer.volume();
      }
    });

    // Compute the fluxes for the boundary faces.
    for (const auto& [label, bc] : *bcs_) {
      std::ranges::for_each(p_mesh_->faces(label), [&](FaceView<Mesh> face) {
        const CellView<Mesh> cell_inner = face.inner_cell();
        Subfield<Real, NumVars> u_outer{};
        bc->get_ghost_state(face.normal3D(), face.center3D(),
                            u[face.inner_cell()].data(), u_outer.data());

        // Compute the flux.
        Subfield<Real, NumVars> flux{};
        flux_scheme_.get_numerical_flux(face.normal3D(), //
                                        u_outer, u[cell_inner], flux);
        const real_t ds = face.area();
        for (size_t i = 0; i < NumVars; ++i) {
          div_f[cell_inner][i] += flux[i] * ds / cell_inner.volume();
        }
      });
    }
  }

}; // class UpwindConvectionScheme

/// @brief Piecewise-linear upwind convection scheme.
/// This is a second order scheme.
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

  void get_cell_convection(Field<real_t, 5>& div_f,
                           const Field<real_t, 5>& u) const final {
    get_cell_convection1(div_f, u);
  }

  /// @brief Compute the second-order upwind nonlinear convection.
  template<class Real, size_t NumVars>
  void get_cell_convection1(Field<Real, NumVars>& div_f,
                            const Field<Real, NumVars>& u) const {
    // Compute the gradients.
    Field<Vec<Real>, NumVars> grad_u(p_mesh_->num_cells());
    gradient_scheme_.get_gradients(grad_u, u);
    Field<Real, NumVars> phi_u(p_mesh_->num_cells());
    gradient_limiter_scheme_.get_cell_limiter(phi_u, u, grad_u);

    // Compute the fluxes for the interior faces.
    std::ranges::for_each(p_mesh_->interior_faces(), [&](FaceView<Mesh> face) {
      const CellView<Mesh> cell_inner = face.inner_cell();
      const CellView<Mesh> cell_outer = face.outer_cell();

      // Reconstruct the face values.
      const vec3_t dr_outer = face.center3D() - cell_outer.center3D();
      const vec3_t dr_inner = face.center3D() - cell_inner.center3D();
      Subfield<Real, NumVars> u_outer{}, u_inner{};
      for (size_t i = 0; i < NumVars; ++i) {
        u_inner[i] =
            u[cell_inner][i] +
            phi_u[cell_inner][i] * glm::dot(grad_u[cell_inner][i], dr_inner);
        u_outer[i] =
            u[cell_outer][i] +
            phi_u[cell_outer][i] * glm::dot(grad_u[cell_outer][i], dr_outer);
      }

      // Compute the flux.
      Subfield<Real, NumVars> flux{};
      flux_scheme_.get_numerical_flux(face.normal3D(), //
                                      u_outer, u_inner, flux);
      const real_t ds = face.area();
      for (size_t i = 0; i < NumVars; ++i) {
        div_f[cell_inner][i] += flux[i] * ds / cell_inner.volume();
        div_f[cell_outer][i] -= flux[i] * ds / cell_outer.volume();
      }
    });

    // Compute the fluxes for the boundary faces.
    for (const auto& [label, bc] : *bcs_) {
      std::ranges::for_each(p_mesh_->faces(label), [&](FaceView<Mesh> face) {
        const CellView<Mesh> cell_inner = face.inner_cell();
        const vec3_t dr_inner = face.center3D() - cell_inner.center3D();

        // Reconstruct the face values.
        Subfield<Real, NumVars> u_outer{}, u_inner{};
        for (size_t i = 0; i < NumVars; ++i) {
          u_inner[i] =
              u[cell_inner][i] +
              phi_u[cell_inner][i] * glm::dot(grad_u[cell_inner][i], dr_inner);
        }
        bc->get_ghost_state(face.normal3D(), face.center3D(), //
                            u_inner.data(), u_outer.data());

        // Compute the flux.
        Subfield<Real, NumVars> flux{};
        flux_scheme_.get_numerical_flux(face.normal3D(), //
                                        u_outer, u_inner, flux);
        const real_t ds = face.area();
        for (size_t i = 0; i < NumVars; ++i) {
          div_f[cell_inner][i] += flux[i] * ds / cell_inner.volume();
        }
      });
    }
  }

}; // class UpwindConvectionScheme

} // namespace Storm::Feathers
