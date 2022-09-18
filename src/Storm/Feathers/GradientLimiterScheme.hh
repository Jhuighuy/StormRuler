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

#include "Field.hh"

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
                                [[maybe_unused]] Real eps_sqr) const {
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
                                Real eps_sqr) const {
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
    const Real delta_pos_sqr = math::pow(delta_pos, 2);
    const Real delta_neg_sqr = math::pow(delta_neg, 2);
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
                                [[maybe_unused]] Real eps_sqr) const {
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
      limiter = y_cur + math::pow(y_div, 2) *
                            (3.0 - 2.0 * y_thr + (y_thr - 2.0) * y_div);
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
                                [[maybe_unused]] Real eps_sqr) const {
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
                                Real eps_sqr) const {
    /* Compute weight:
     * [1], page 5. */
    const Real du_sqr = std::pow(du_max - du_min, 2);
    if (du_sqr <= eps_sqr) { return 1.0; }
    if (eps_sqr < du_sqr && du_sqr < 2.0 * eps_sqr) {
      const Real dy = (du_sqr - eps_sqr) / eps_sqr;
      const Real dy_sqr = std::pow(dy, 2);
      const Real weight = (2.0 * dy - 3.0) * dy_sqr + 1.0;
      const Real second_limiter = weight + (1.0 - weight) * limiter;
      return second_limiter;
    }
    if (du_sqr >= 2.0 * eps_sqr) { return limiter; }
    STORM_THROW_( //
        "Broken cubic second limiter, du^2 = {}, eps^2 = {}!", du_sqr, eps_sqr);
  }

}; // class CubicSecondLimiter

/**
 * Gradient cell limiter estimation scheme:
 * computes cell-centered limiters and averages based on the cell-centered
 * expansion.
 */
class iGradientLimiterScheme : public tObject<iGradientLimiterScheme> {
public:

  /** Compute cell-centered gradient limiter coefficients. */
  virtual void get_cell_limiter(size_t num_vars, tScalarField& lim_u,
                                const tScalarField& u,
                                const tVectorField& grad_u) const = 0;

}; // class iGradientLimiterScheme

/**
 * Gradient cell limiter estimation scheme:
 * computes cell-centered limiters and averages based on the cell-centered
 * expansion.
 */
template<class tSlopeLimiter, class tSecondLimiter = DummySecondLimiter>
class tGradientLimiterScheme final : public iGradientLimiterScheme {
public:

  std::shared_ptr<const Mesh> m_mesh;
  tSlopeLimiter m_slope_limiter;
  tSecondLimiter m_second_limiter;

public:

  /** Initialize the limiting scheme. */
  explicit tGradientLimiterScheme(std::shared_ptr<const Mesh> mesh,
                                  const tSlopeLimiter& slope_limiter = {},
                                  const tSecondLimiter& second_limiter = {})
      : m_mesh(std::move(mesh)), m_slope_limiter(slope_limiter),
        m_second_limiter(second_limiter) {}

  /** Compute cell-centered gradient limiter coefficients. */
  void get_cell_limiter(size_t num_vars, tScalarField& lim_u,
                        const tScalarField& u,
                        const tVectorField& grad_u) const final {
    /* Compute the cell-centered
     * limiting coefficients and averages. */
    ForEach(m_mesh->interior_cells(), [&](CellView<Mesh> cell) {
      static const real_t k = 0.1;
      const real_t eps_sqr = std::pow(k * cell.volume(), 3);
      /* Find the largest negative and positive differences
       * between values of and neighbor Cells and the current cell. */
      FEATHERS_TMP_SCALAR_FIELD(du_min, num_vars);
      du_min = u[cell];
      FEATHERS_TMP_SCALAR_FIELD(du_max, num_vars);
      du_max = u[cell];
      cell.for_each_face_cells([&](CellView<Mesh> cell_inner,
                                   CellView<Mesh> cell_outer) {
        for (size_t i = 0; i < num_vars; ++i) {
          du_min[i] =
              std::min(du_min[i], std::min(u[cell_outer][i], u[cell_inner][i]));
          du_max[i] =
              std::max(du_max[i], std::max(u[cell_outer][i], u[cell_inner][i]));
        }
      });
      for (size_t i = 0; i < num_vars; ++i) {
        du_min[i] = std::min(0.0, du_min[i] - u[cell][i]);
        du_max[i] = std::max(0.0, du_max[i] - u[cell][i]);
      }

      /* Compute slope limiting coefficients:
       * clamp the Node delta with computed local delta extrema. */
      lim_u[cell].fill(1.0);
      cell.for_each_face([&](FaceView<Mesh> face) {
        const vec3_t dr = face.center3D() - cell.center3D();
        for (size_t i = 0; i < num_vars; ++i) {
          const real_t du_face = glm::dot(grad_u[cell][i], dr);
          const real_t limiter =
              m_slope_limiter(du_min[i], du_max[i], du_face, eps_sqr);
          lim_u[cell][i] = std::min(lim_u[cell][i], limiter);
        }
      });

      /* Compute secondary limiting coefficients:
       * disable limiting near smooth regions. */
      for (size_t i = 0; i < num_vars; ++i) {
        const real_t limiter =
            m_second_limiter(lim_u[cell][i], du_min[i], du_max[i], eps_sqr);
        lim_u[cell][i] = limiter;
      }
    });
  }

}; // class tGradientLimiterScheme

using cMinmodGradientLimiterScheme = tGradientLimiterScheme<MinmodSlopeLimiter>;

using cVenkatakrishnanGradientLimiterScheme =
    tGradientLimiterScheme<VenkatakrishnanSlopeLimiter>;

using cVenkatakrishnan2GradientLimiterScheme =
    tGradientLimiterScheme<VenkatakrishnanSlopeLimiter, CubicSecondLimiter>;

using cCubicGradientLimiterScheme = tGradientLimiterScheme<CubicSlopeLimiter>;

using cCubic2GradientLimiterScheme =
    tGradientLimiterScheme<CubicSlopeLimiter, CubicSecondLimiter>;

} // namespace Storm::Feathers
