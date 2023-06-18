// ╔═══════════════════════════════════════════════════════════════════════════╗
// ║   ______  ______   ______   ______  __  __   ______   ______   ______     ║
// ║  /\  ___\/\  ___\ /\  __ \ /\__  _\/\ \_\ \ /\  ___\ /\  __ \ /\  ___\    ║
// ║  \ \  __\\ \  _\  \ \  __ \\/_/\ \/\ \  __ \\ \  __\ \ \  __/ \ \___  \   ║
// ║   \ \_\   \ \_____\\ \_\ \_\  \ \_\ \ \_\ \_\\ \_____\\ \_\ \_\\/\_____\  ║
// ║    \/_/    \/_____/ \/_/\/_/   \/_/  \/_/\/_/ \/_____/ \/_/ /_/ \/_____/  ║
// ╚═══════════════════════════════════════════════════════════════════════════╝
//
// Copyright (C) 2020-2023 Oleg Butakov
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include "./FluxScheme.hpp"
#include "./GradientLimiterScheme.hpp"
#include "./GradientScheme.hpp"
#include "./SkunkFvBC.hpp"

#include <map>

namespace Storm::Feathers {

/**
 * Abstract convection scheme.
 */
class iConvectionScheme : public tObject<iConvectionScheme> {
public:

  std::map<Label, std::shared_ptr<MhdFvBcPT<tGasPhysics>>>* _bcs;

  virtual ~iConvectionScheme() = default;

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

  const Mesh* _p_mesh;
  STORM_NO_UNIQUE_ADDRESS FluxScheme _flux_scheme;

public:

  /// @brief Construct the first order upwind convection scheme.
  constexpr explicit UpwindConvectionScheme(const Mesh& mesh,
                                            FluxScheme flux_scheme) noexcept
      : _p_mesh{&mesh}, _flux_scheme{std::move(flux_scheme)} {}

  void get_cell_convection(CellField<Mesh, real_t, 5>& div_f,
                           const CellField<Mesh, real_t, 5>& u) const final {
    (*this)(div_f, u);
  }

  /// @brief Compute the first order upwind nonlinear convection.
  template<class Real, size_t NumVars>
  void operator()(CellField<Mesh, Real, NumVars>& div_f,
                  const CellField<Mesh, Real, NumVars>& u) const noexcept {
    // Compute the fluxes for the interior faces.
    std::ranges::for_each(_p_mesh->interior_faces(), [&](FaceView<Mesh> face) {
      const CellView<Mesh> cell_inner = face.inner_cell();
      const CellView<Mesh> cell_outer = face.outer_cell();

      // Compute the flux.
      const auto flux =
          _flux_scheme(face.normal(), u[cell_outer], u[cell_inner]);
      div_f[cell_inner] += (face.area() / cell_inner.volume()) * flux;
      div_f[cell_outer] -= (face.area() / cell_outer.volume()) * flux;
    });

    // Compute the fluxes for the boundary faces.
    for (const auto& [label, bc] : *_bcs) {
      std::ranges::for_each(_p_mesh->faces(label), [&](FaceView<Mesh> face) {
        const CellView<Mesh> cell_inner = face.inner_cell();
        Subfield<Real, NumVars> u_outer{};
        bc->get_ghost_state(face.normal(), face.center(),
                            u[face.inner_cell()].data(), u_outer.data());

        // Compute the flux.
        const auto flux = _flux_scheme(face.normal(), u_outer, u[cell_inner]);
        div_f[cell_inner] += (face.area() / cell_inner.volume()) * flux;
      });
    }
  }

}; // class UpwindConvectionScheme

// -----------------------------------------------------------------------------

/// @brief Piecewise-linear upwind convection scheme.
template<mesh Mesh, class FluxScheme, //
         class GradientScheme, class GradientLimiterScheme>
  requires std::is_object_v<FluxScheme> && std::is_object_v<GradientScheme> &&
           std::is_object_v<GradientLimiterScheme>
class LinearUpwindConvectionScheme final : public iConvectionScheme {
public:

  const Mesh* _p_mesh;
  STORM_NO_UNIQUE_ADDRESS FluxScheme _flux_scheme;
  STORM_NO_UNIQUE_ADDRESS GradientScheme _gradient_scheme;
  STORM_NO_UNIQUE_ADDRESS GradientLimiterScheme _gradient_limiter_scheme;

public:

  /// @brief Construct the second order upwind convection scheme.
  constexpr explicit LinearUpwindConvectionScheme(
      const Mesh& mesh, FluxScheme flux_scheme, GradientScheme gradient_scheme,
      GradientLimiterScheme gradient_limiter_scheme) noexcept
      : _p_mesh{&mesh}, _flux_scheme{flux_scheme}, //
        _gradient_scheme{std::move(gradient_scheme)},
        _gradient_limiter_scheme{std::move(gradient_limiter_scheme)} {}

  void get_cell_convection(CellField<Mesh, real_t, 5>& div_f,
                           const CellField<Mesh, real_t, 5>& u) const final {
    (*this)(div_f, u);
  }

  /// @brief Compute the second-order upwind nonlinear convection.
  template<class Real, size_t NumVars>
  void operator()(CellField<Mesh, Real, NumVars>& div_f,
                  const CellField<Mesh, Real, NumVars>& u) const noexcept {
    // Compute the gradients.
    CellVectorField<Mesh, Real, NumVars> grad_u(_p_mesh->num_cells());
    CellField<Mesh, Real, NumVars> lim_u(_p_mesh->num_cells());
    _gradient_scheme(grad_u, u);
    _gradient_limiter_scheme(lim_u, u, grad_u);
    std::ranges::for_each(_p_mesh->interior_cells(), [&](CellView<Mesh> cell) {
      grad_u[cell] *= lim_u[cell];
    });

    // Compute the fluxes for the interior faces.
    std::ranges::for_each(_p_mesh->interior_faces(), [&](FaceView<Mesh> face) {
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
      const auto flux = _flux_scheme(face.normal(), u_outer, u_inner);
      div_f[cell_inner] += (face.area() / cell_inner.volume()) * flux;
      div_f[cell_outer] -= (face.area() / cell_outer.volume()) * flux;
    });

    // Compute the fluxes for the boundary faces.
    for (const auto& [label, bc] : *_bcs) {
      std::ranges::for_each(_p_mesh->faces(label), [&](FaceView<Mesh> face) {
        const CellView<Mesh> cell_inner = face.inner_cell();

        // Reconstruct the face values.
        const auto dr_inner = face.center() - cell_inner.center();
        Subfield<Real, NumVars> u_outer{}, u_inner{};
        for (size_t i = 0; i < NumVars; ++i) {
          u_inner[i] =
              u[cell_inner][i] + glm::dot(grad_u[cell_inner][i], dr_inner);
        }
        bc->get_ghost_state(face.normal(), face.center(), //
                            u_inner.data(), u_outer.data());

        // Compute the flux.
        const auto flux = _flux_scheme(face.normal(), u_outer, u_inner);
        div_f[cell_inner] += (face.area() / cell_inner.volume()) * flux;
      });
    }
  }

}; // class UpwindConvectionScheme

} // namespace Storm::Feathers
