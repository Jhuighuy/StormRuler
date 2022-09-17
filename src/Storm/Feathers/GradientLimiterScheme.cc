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

#include "GradientLimiterScheme.hh"

namespace Storm {

/**
 * Compute local slope coefficient.
 * @verbatim
 * [1] Krzysztof Michalak, Carl Ollivier-Gooch,
 *     "Limiters for Unstructured Higher-Order Accurate
 *      Solutions of the Euler Equations" (2008).
 * @endverbatim
 */
inline real_t
tMinmodSlopeLimiter::operator()(real_t du_min, real_t du_max, real_t du_face,
                                real_t FEATHERS_NOT_USED(eps_sqr)) const {
  /* Compute deltas:
   * [1], page 4. */
  const real_t delta_neg = du_face;
  real_t delta_pos;
  if (delta_neg < 0.0) {
    delta_pos = du_min;
  } else if (delta_neg > 0.0) {
    delta_pos = du_max;
  } else {
    return 1.0;
  }
  /* Compute limiter:
   * [1], page 3. */
  const real_t y_cur = delta_pos / delta_neg;
  const real_t limiter = std::min(1.0, y_cur);
  return limiter;
} // tMinmodSlopeLimiter::get_limiter_coefficient

/**
 * Compute local slope coefficient.
 * @verbatim
 * [1] Krzysztof Michalak, Carl Ollivier-Gooch,
 *     "Limiters for Unstructured Higher-Order Accurate
 *      Solutions of the Euler Equations" (2008).
 * @endverbatim
 */
inline real_t tVenkatakrishnanSlopeLimiter::operator()(real_t du_min,
                                                       real_t du_max,
                                                       real_t du_face,
                                                       real_t eps_sqr) const {
  /* Compute deltas:
   * [1], page 4. */
  const real_t delta_neg = du_face;
  real_t delta_pos;
  if (delta_neg < 0.0) {
    delta_pos = du_min;
  } else if (delta_neg > 0.0) {
    delta_pos = du_max;
  } else {
    return 1.0;
  }
  /* Compute limiter:
   * [1], page 4. */
  const real_t delta_pos_sqr = std::pow(delta_pos, 2);
  const real_t delta_neg_sqr = std::pow(delta_neg, 2);
  const real_t delta_pos_neg = delta_pos * delta_neg;
  const real_t limiter =
      (delta_pos_sqr + 2.0 * delta_pos_neg + eps_sqr) /
      (delta_pos_sqr + 2.0 * delta_neg_sqr + delta_pos_neg + eps_sqr);
  return limiter;
} // tVenkatakrishnanSlopeLimiter::operator()

/**
 * Compute local slope coefficient.
 * @verbatim
 * [1] Krzysztof Michalak, Carl Ollivier-Gooch,
 *     "Limiters for Unstructured Higher-Order Accurate
 *      Solutions of the Euler Equations" (2008).
 * @endverbatim
 */
inline real_t
tCubicSlopeLimiter::operator()(real_t du_min, real_t du_max, real_t du_face,
                               real_t FEATHERS_NOT_USED(eps_sqr)) const {
  /* Calculate deltas:
   * [1], page 4. */
  const real_t delta_neg = du_face;
  real_t delta_pos;
  if (delta_neg < 0.0) {
    delta_pos = du_min;
  } else if (delta_neg > 0.0) {
    delta_pos = du_max;
  } else {
    return 1.0;
  }
  /* Compute limiter:
   * [1], page 5. */
  const real_t y_cur = delta_pos / delta_neg;
  const real_t y_thr = 1.75;
  auto limiter = 1.0;
  if (y_cur < y_thr) {
    const real_t y_div = y_cur / y_thr;
    limiter = y_cur +
              std::pow(y_div, 2) * (3.0 - 2.0 * y_thr + (y_thr - 2.0) * y_div);
  }
  return limiter;
} // tCubicSlopeLimiter::operator()

// ------------------------------------------------------------------------------------
// //
// ------------------------------------------------------------------------------------
// //

/**
 * Compute second slope coefficient.
 */
inline real_t tDummySecondLimiter::operator()(
    real_t limiter, real_t FEATHERS_NOT_USED(du_min),
    real_t FEATHERS_NOT_USED(du_max), real_t FEATHERS_NOT_USED(eps_sqr)) const {
  const real_t second_limiter = limiter;
  return second_limiter;
} // tDummySecondLimiter::operator()

/**
 * Compute second slope coefficient.
 * @verbatim
 * [1] Krzysztof Michalak, Carl Ollivier-Gooch,
 *     "Limiters for Unstructured Higher-Order Accurate
 *      Solutions of the Euler Equations" (2008).
 * @endverbatim
 */
inline real_t tCubicSecondLimiter::operator()(real_t limiter, real_t du_min,
                                              real_t du_max,
                                              real_t eps_sqr) const {
  /* Compute weight:
   * [1], page 5. */
  const real_t du_sqr = std::pow(du_max - du_min, 2);
  if (du_sqr <= eps_sqr) { return 1.0; }
  if (eps_sqr < du_sqr && du_sqr < 2.0 * eps_sqr) {
    const real_t dy = (du_sqr - eps_sqr) / eps_sqr;
    const real_t dy_sqr = std::pow(dy, 2);
    const real_t weight = (2.0 * dy - 3.0) * dy_sqr + 1.0;
    const real_t second_limiter = weight + (1.0 - weight) * limiter;
    return second_limiter;
  }
  if (du_sqr >= 2.0 * eps_sqr) { return limiter; }
  FEATHERS_ENSURE(!"Broken cubic second limiter.");
} // tCubicSecondLimiter::get_slope_coefficient

// ------------------------------------------------------------------------------------
// //
// ------------------------------------------------------------------------------------
// //

/**
 * Compute cell-centered gradient limiter coefficients.
 */
template<class tSlopeLimiter, class tSecondLimiter>
void tGradientLimiterScheme<tSlopeLimiter, tSecondLimiter>::get_cell_limiter(
    size_t num_vars, tScalarField& lim_u, const tScalarField& u,
    const tVectorField& grad_u) const {
  /* Compute the cell-centered
   * limiting coefficients and averages. */
  ForEach(int_cell_views(*m_mesh), [&](CellView cell) {
    static const real_t k = 0.1;
    const real_t eps_sqr = std::pow(k * cell.volume(), 3);
    /* Find the largest negative and positive differences
     * between values of and neighbor Cells and the current cell. */
    FEATHERS_TMP_SCALAR_FIELD(du_min, num_vars);
    du_min = u[cell];
    FEATHERS_TMP_SCALAR_FIELD(du_max, num_vars);
    du_max = u[cell];
    cell.for_each_face_cells([&](CellView cell_inner, CellView cell_outer) {
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
    cell.for_each_face([&](FaceView face) {
      const vec3_t dr = face.center() - cell.center();
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
} // tGradientLimiterScheme::get_cell_limiter

template class tGradientLimiterScheme<tMinmodSlopeLimiter>;
template class tGradientLimiterScheme<tVenkatakrishnanSlopeLimiter>;
template class tGradientLimiterScheme<tVenkatakrishnanSlopeLimiter,
                                      tCubicSecondLimiter>;
template class tGradientLimiterScheme<tCubicSlopeLimiter>;
template class tGradientLimiterScheme<tCubicSlopeLimiter, tCubicSecondLimiter>;

} // namespace Storm
