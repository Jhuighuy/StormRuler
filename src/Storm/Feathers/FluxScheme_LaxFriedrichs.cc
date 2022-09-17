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

#include "FluxScheme.hh"

namespace Storm {

/**
 * Calculate the Local Lax-Friedrichs (Rusanov) numerical flux.
 *
 * @verbatim
 * [1] Eleuterio F. Toro,
 *     "Riemann Solvers and Numerical Methods
 *      for Fluid Dynamics" (Third Edition, 2009).
 * @endverbatim
 */
void tLaxFriedrichsFluxScheme<tGasPhysics>::get_numerical_flux(
    size_t num_vars, const vec3_t& n, tScalarConstSubField cons_r,
    tScalarConstSubField cons_l, tScalarSubField flux) const {
  const tGasPhysics::tFluidState ur(n, cons_r.data());
  const tGasPhysics::tFluidState ul(n, cons_l.data());

  /* Approximate |J| with its maximum eigenvalue.
   * [1] Eq. (10.55-10.56). */
  const real_t ss =
      std::max(std::abs(ur.vel_n) + ur.c_snd, std::abs(ul.vel_n) + ul.c_snd);
  FEATHERS_TMP_SCALAR_FIELD(flux_r, num_vars);
  FEATHERS_TMP_SCALAR_FIELD(flux_l, num_vars);
  ur.make_flux(num_vars, n, flux_r.data());
  ul.make_flux(num_vars, n, flux_l.data());
  for (size_t i = 0; i < num_vars; ++i) {
    flux[i] = 0.5 * ((flux_r[i] + flux_l[i]) - ss * (cons_r[i] - cons_l[i]));
  }
} // tLaxFriedrichsFluxScheme::get_numerical_flux

} // namespace Storm
