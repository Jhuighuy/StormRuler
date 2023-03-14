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

#include "./Field.hpp"

namespace Storm::Feathers
{

/// @brief Least Squares gradient estimation scheme, cell-based:
/// computes cell-centered gradients based on the cell-centered values.
template<mesh Mesh>
class LeastSquaresGradientScheme final
{
private:

  const Mesh* p_mesh_;
  CellMatrixField<Mesh, real_t> g_mats_;

public:

  /// @brief Construct the gradient scheme.
  constexpr explicit LeastSquaresGradientScheme(const Mesh& mesh)
      : p_mesh_{&mesh}, g_mats_(p_mesh_->num_cells())
  {
    // Compute the inverse least squares matrices.
    std::ranges::for_each(p_mesh_->interior_cells(), [&](CellView<Mesh> cell) {
      // Compute the direct least squares matrices.
      g_mats_[cell][0] = {};
      cell.for_each_cell([&](CellView<Mesh> adj_cell) {
        const auto dr = adj_cell.center() - cell.center();
        g_mats_[cell][0] += glm::outerProduct(dr, dr);
      });

      // Compute the inverse.
      static const mat2_t eps(1e-14);
      g_mats_[cell][0] = glm::inverse(eps + g_mats_[cell][0]);
    });
  }

  /// @brief Compute cell-centered gradients.
  template<class Real, size_t NumVars>
  void operator()(CellVectorField<Mesh, Real, NumVars>& grad_u,
                  const CellField<Mesh, Real, NumVars>& u) const noexcept
  {
    std::ranges::for_each(p_mesh_->interior_cells(), [&](CellView<Mesh> cell) {
      // Compute the least squares RHS.
      grad_u[cell].fill({});
      cell.for_each_cell([&](CellView<Mesh> adj_cell) {
        const auto dr = adj_cell.center() - cell.center();
        for (size_t i = 0; i < NumVars; ++i) {
          grad_u[cell][i] += (u[adj_cell][i] - u[cell][i]) * dr;
        }
      });

      // Compute the least squares gradients.
      for (size_t i = 0; i < NumVars; ++i) {
        grad_u[cell][i] = g_mats_[cell][0] * grad_u[cell][i];
      }
    });
  }

}; // class LeastSquaresGradientScheme

} // namespace Storm::Feathers
