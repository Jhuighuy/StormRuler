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

#include <Storm/Utils/Math.hpp>

#include "Field.hpp"

#include <type_traits>

namespace Storm::Feathers {

/// @brief Barth-Jespersen (minmod) slope limiter.
/// This is a non-differentiable limiter, so it may
/// affect convergence properties of the implicit scheme.
class MinmodSlopeLimiter {
public:

  /// @brief Compute limiter value.
  template<class Real>
  [[nodiscard]] Real operator()(Real du_min,  //
                                Real du_max,  //
                                Real du_face, //
                                [[maybe_unused]] Real eps_sqr) const noexcept {
    // Compute deltas: [1], page 4.
    const Real delta_neg = du_face;
    Real delta_pos;
    if (delta_neg < 0.0) {
      delta_pos = du_min;
    } else if (delta_neg > 0.0) {
      delta_pos = du_max;
    } else {
      return 1.0;
    }
    // Compute limiter: [1], page 3.
    const Real y_cur = delta_pos / delta_neg;
    const Real limiter = std::min(1.0, y_cur);
    return limiter;
  }

}; // class MinmodSlopeLimiter

/// @brief Venkatakrishnan slope limiter.
/// This is a differentiable limiter.
class VenkatakrishnanSlopeLimiter {
public:

  /// @brief Compute limiter value.
  template<class Real>
  [[nodiscard]] Real operator()(Real du_min,  //
                                Real du_max,  //
                                Real du_face, //
                                Real eps_sqr) const noexcept {
    // Compute deltas: [1], page 4.
    const Real delta_neg = du_face;
    Real delta_pos;
    if (delta_neg < 0.0) {
      delta_pos = du_min;
    } else if (delta_neg > 0.0) {
      delta_pos = du_max;
    } else {
      return 1.0;
    }
    // Compute limiter: [1], page 4.
    const Real delta_pos_sqr = pow(delta_pos, 2);
    const Real delta_neg_sqr = pow(delta_neg, 2);
    const Real delta_pos_neg = delta_pos * delta_neg;
    const Real limiter =
        (delta_pos_sqr + 2.0 * delta_pos_neg + eps_sqr) /
        (delta_pos_sqr + 2.0 * delta_neg_sqr + delta_pos_neg + eps_sqr);
    return limiter;
  }

}; // class VenkatakrishnanSlopeLimiter

/// @brief Michalak Ollivier-Gooch (cubic) slope limiter.
/// This is a differentiable limiter.
class CubicSlopeLimiter {
public:

  /// @brief Compute limiter value.
  template<class Real>
  [[nodiscard]] Real operator()(Real du_min,  //
                                Real du_max,  //
                                Real du_face, //
                                [[maybe_unused]] Real eps_sqr) const noexcept {
    // Compute deltas: [1], page 4.
    const Real delta_neg = du_face;
    Real delta_pos;
    if (delta_neg < 0.0) {
      delta_pos = du_min;
    } else if (delta_neg > 0.0) {
      delta_pos = du_max;
    } else {
      return 1.0;
    }
    // Compute limiter: [1], page 5.
    const Real y_cur = delta_pos / delta_neg;
    const Real y_thr = 1.75;
    auto limiter = 1.0;
    if (y_cur < y_thr) {
      const Real y_div = y_cur / y_thr;
      limiter =
          y_cur + pow(y_div, 2) * (3.0 - 2.0 * y_thr + (y_thr - 2.0) * y_div);
    }
    return limiter;
  }

}; // class CubicSlopeLimiter

/// @brief Dummy second slope limiter.
/// This is a differentiable limiter.
class DummySecondLimiter {
public:

  /// @brief Compute limiter value.
  template<class Real>
  [[nodiscard]] Real operator()(Real limiter, //
                                [[maybe_unused]] Real du_min,
                                [[maybe_unused]] Real du_max,
                                [[maybe_unused]] Real eps_sqr) const noexcept {
    const Real second_limiter = limiter;
    return second_limiter;
  }

}; // class DummySecondLimiter

/// @brief Michalak Ollivier-Gooch cubic second slope limiter.
/// This is a differentiable limiter.
class CubicSecondLimiter {
public:

  /// @brief Compute limiter value.
  template<class Real>
  [[nodiscard]] Real operator()(Real limiter, //
                                Real du_min,  //
                                Real du_max,  //
                                Real eps_sqr) const noexcept {
    const Real du_sqr = std::pow(du_max - du_min, 2);
    if (du_sqr <= eps_sqr) { return 1.0; }
    if (du_sqr >= 2.0 * eps_sqr) { return limiter; }
    // Compute weight: [1], page 5.
    const Real dy = (du_sqr - eps_sqr) / eps_sqr;
    const Real dy_sqr = pow(dy, 2);
    const Real weight = (2.0 * dy - 3.0) * dy_sqr + 1.0;
    const Real second_limiter = weight * 1.0 + (1.0 - weight) * limiter;
    return second_limiter;
  }

}; // class CubicSecondLimiter

/// @brief Gradient cell limiter estimation scheme.
///
/// References:
/// @verbatim
/// [1] Krzysztof Michalak, Carl Ollivier-Gooch,
///     "Limiters for Unstructured Higher-Order Accurate
///      Solutions of the Euler Equations" (2008).
/// @endverbatim
template<mesh Mesh, class SlopeLimiter, class SecondLimiter>
  requires std::is_object_v<SlopeLimiter> &&
           std::is_object_v<DummySecondLimiter>
class GradientLimiterScheme final {
private:

  const Mesh* p_mesh_;
  STORM_NO_UNIQUE_ADDRESS_ SlopeLimiter slope_limiter_;
  STORM_NO_UNIQUE_ADDRESS_ SecondLimiter second_limiter_;

public:

  /// @brief Construct the gradient limiter.
  constexpr explicit GradientLimiterScheme(
      const Mesh& mesh, SlopeLimiter slope_limiter = {},
      SecondLimiter second_limiter = {}) noexcept
      : p_mesh_{&mesh}, slope_limiter_{std::move(slope_limiter)},
        second_limiter_{std::move(second_limiter)} {}

  /// @brief Compute cell-centered gradient limiter coefficients.
  template<class Real, size_t NumVars>
  void operator()( //
      CellField<Mesh, Real, NumVars>& lim_u,
      const CellField<Mesh, Real, NumVars>& u,
      const CellVectorField<Mesh, Real, NumVars>& grad_u) const noexcept {
    std::ranges::for_each(p_mesh_->interior_cells(), [&](CellView<Mesh> cell) {
      // Find the largest negative and positive differences
      // between values of the current cell and the adjacent cells.
      auto du_min = u[cell], du_max = u[cell];
      cell.for_each_cell([&](CellView<Mesh> adj_cell) {
        du_min = min(du_min, u[adj_cell]);
        du_max = max(du_max, u[adj_cell]);
      });
      du_min = min(0.0, du_min - u[cell]);
      du_max = max(0.0, du_max - u[cell]);

      // Compute slope limiting coefficients:
      // clamp the face delta with computed local delta extrema.
      static const real_t k = 0.1;
      const real_t eps_sqr = std::pow(k * cell.volume(), 3);
      lim_u[cell].fill(1.0);
      cell.for_each_face([&](FaceView<Mesh> face) {
        const auto dr = face.center() - cell.center();
        for (size_t i = 0; i < NumVars; ++i) {
          const Real du_face = glm::dot(grad_u[cell][i], dr);
          const Real limiter =
              slope_limiter_(du_min[i], du_max[i], du_face, eps_sqr);
          lim_u[cell][i] = std::min(lim_u[cell][i], limiter);
        }
      });

      // Compute secondary limiting coefficients by disabling
      // limiting near the smooth regions.
      for (size_t i = 0; i < NumVars; ++i) {
        const Real limiter =
            second_limiter_(lim_u[cell][i], du_min[i], du_max[i], eps_sqr);
        lim_u[cell][i] = limiter;
      }
    });
  }

}; // class GradientLimiterScheme

} // namespace Storm::Feathers
