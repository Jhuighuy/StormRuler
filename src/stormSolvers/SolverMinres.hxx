/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// Copyright (C) 2022 Oleg Butakov
///
/// Permission is hereby granted, free of charge, to any person
/// obtaining a copy of this software and associated documentation
/// files (the "Software"), to deal in the Software without
/// restriction, including without limitation the rights  to use,
/// copy, modify, merge, publish, distribute, sublicense, and/or
/// sell copies of the Software, and to permit persons to whom the
/// Software is furnished to do so, subject to the following
/// conditions:
///
/// The above copyright notice and this permission notice shall be
/// included in all copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
/// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
/// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
/// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
/// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
/// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
/// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
/// OTHER DEALINGS IN THE SOFTWARE.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
#if 0
#pragma once

#include <cmath>

#include <stormBase.hxx>
#include <stormSolvers/Solver.hxx>

namespace Storm {

namespace Blas {

/// @brief Generate Givens rotation.
template<class Real>
inline auto sym_ortho(Real a, Real b) {
  // ----------------------
  // 𝑟𝑟 ← (𝑎² + 𝑏²)¹ᐟ²,
  // 𝑐𝑠 ← 𝑎/𝑟𝑟, 𝑠𝑛 ← 𝑏/𝑟𝑟.
  // ----------------------
  Real cs, sn, rr;
  rr = std::hypot(a, b);
  if (rr > 0.0) {
    cs = a / rr;
    sn = b / rr;
  } else {
    cs = 1.0;
    sn = 0.0;
  }

  return std::make_tuple(cs, sn, rr);

} // sym_ortho

} // namespace Blas

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief The @c MINRES (Minimal Residual) linear self-adjoint 
///   (indefinite) operator equation solver.
///
/// @c MINRES can be applied to the singular problems, and the 
///   self-adjoint least squares problems, although convergeance to 
///   minimum norm solution is not guaranteed.
///
/// @note Despite 𝓐 may be indefinite, a positive-definite 
///   preconditioner 𝓟 is explicitly required.
///
/// References:
/// @verbatim
/// [1] Paige, C. and M. Saunders.
///     “Solution of Sparse Indefinite Systems of Linear Equations.”
///     SIAM Journal on Numerical Analysis 12 (1975): 617-629.
/// [2] Choi, S.-C. T.
///     “Iterative Methods for Singular Linear Equations and
///     Least-Squares Problems” PhD thesis, ICME, Stanford University.
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Vector>
class MinresSolver final : public IterativeSolver<Vector> {
private:

  real_t alpha, beta, betaBar, gamma, delta, deltaBar, epsilon, epsilonBar, tau,
      phi, phiTilde, cs, sn;
  Vector p_vec, q_vec, qBarVec, w_vec, wBarVec, wBarBarVec, z_vec, zBarVec,
      zBarBarVec;

  real_t init(const Vector& x_vec, const Vector& b_vec,
              const Operator<Vector>& lin_op,
              const Preconditioner<Vector>* pre_op) override;

  real_t iterate(Vector& x_vec, const Vector& b_vec,
                 const Operator<Vector>& lin_op,
                 const Preconditioner<Vector>* pre_op) override;

}; // class MinresSolver

template<class Vector>
real_t MinresSolver<Vector>::init(const Vector& x_vec, const Vector& b_vec,
                                  const Operator<Vector>& lin_op,
                                  const Preconditioner<Vector>* pre_op) {
  assert(pre_op != nullptr && "MINRES requires preconditioning for now.");

  p_vec.assign(x_vec, false);
  w_vec.assign(x_vec, false);
  wBarVec.assign(x_vec, false);
  wBarBarVec.assign(x_vec, false);
  z_vec.assign(x_vec, false);
  zBarVec.assign(x_vec, false);
  zBarBarVec.assign(x_vec, false);
  if (pre_op != nullptr) {
    q_vec.assign(x_vec, false);
    qBarVec.assign(x_vec, false);
  }

  // ----------------------
  // Initialize:
  // 𝒘̅ ← {𝟢}ᵀ,
  // 𝒘̿ ← {𝟢}ᵀ,
  // 𝒛̅ ← 𝓐𝒙,     // Modification in order to
  // 𝒛̅ ← 𝒃 - 𝒛̅,  // utilize the initial guess.
  // 𝒛̿ ← {𝟢}ᵀ,
  // 𝒒 ← [𝓟]𝒛̅,
  // 𝛽̅ ← 𝟣, 𝛽 ← √<𝒒⋅𝒛̅>,
  // 𝜑 ← 𝛽, 𝛿 ← 𝟢, 𝜀 ← 𝟢,
  // 𝑐𝑠 ← -𝟣, 𝑠𝑛 ← 𝟢.
  // ----------------------
  fill_with(wBarVec, 0.0);
  fill_with(wBarBarVec, 0.0);
  lin_op.mul(zBarVec, x_vec);
  Blas::Sub(zBarVec, b_vec, zBarVec);
  fill_with(zBarBarVec, 0.0);
  if (pre_op != nullptr) {
    pre_op->mul(q_vec, zBarVec);
  } else {
    // q_vec = zBarVec
  }
  betaBar = 1.0;
  beta = std::sqrt(dot_product(q_vec, zBarVec));
  phi = beta;
  delta = 0.0;
  epsilon = 0.0;
  cs = -1.0;
  sn = 0.0;

  return phi;

} // MinresSolver::init

template<class Vector>
real_t MinresSolver<Vector>::iterate(Vector& x_vec, const Vector& b_vec,
                                     const Operator<Vector>& lin_op,
                                     const Preconditioner<Vector>* pre_op) {
  assert(pre_op != nullptr && "MINRES requires preconditioning for now.");

  // ----------------------
  // Continue the Lanczos process:
  // 𝒑 ← 𝓐𝒒,
  // 𝛼 ← <𝒒⋅𝒑>/𝛽²,
  // 𝒛 ← (𝟣/𝛽)⋅𝒑 - (𝛼/𝛽)⋅𝒛̅,
  // 𝒛 ← 𝒛 - (𝛽/𝛽̅)⋅𝒛̿,
  // 𝒒̅ ← 𝒒, 𝒒 ← [𝓟]𝒛,
  // 𝛽̅ ← 𝛽, 𝛽 ← <𝒒⋅𝒛>¹ᐟ²,
  // 𝒛̿ ← 𝒛̅, 𝒛̅ ← 𝒛.
  // ----------------------
  lin_op.mul(p_vec, q_vec);
  alpha = dot_product(q_vec, p_vec) * std::pow(beta, -2);
  Blas::Sub(z_vec, p_vec, 1.0 / beta, zBarVec, alpha / beta);
  Blas::Sub(z_vec, z_vec, zBarBarVec, beta / betaBar);
  if (pre_op != nullptr) {
    std::swap(qBarVec, q_vec);
    pre_op->mul(q_vec, z_vec);
  } else {
    // qBarVec = q_vec; q_vec = z_vec
  }
  betaBar = beta, beta = std::sqrt(dot_product(q_vec, z_vec));
  std::swap(zBarBarVec, zBarVec), std::swap(zBarVec, z_vec);

  // ----------------------
  // Construct and apply rotations:
  // 𝛿̅ ← 𝑐𝑠⋅𝛿 + 𝑠𝑛⋅𝛼, 𝛾 ← 𝑠𝑛⋅𝛿 - 𝑐𝑠⋅𝛼,
  // 𝜀̅ ← 𝜀, 𝜀 ← 𝑠𝑛⋅𝛽, 𝛿 ← -𝑐𝑠⋅𝛽,
  // 𝑐𝑠, 𝑠𝑛, 𝛾 ← 𝘚𝘺𝘮𝘖𝘳𝘵𝘩𝘰(𝛾, 𝛽),
  // 𝜏 ← 𝑐𝑠⋅𝜑, 𝜑 ← 𝑠𝑛⋅𝜑.
  // ----------------------
  deltaBar = cs * delta + sn * alpha, gamma = sn * delta - cs * alpha;
  epsilonBar = epsilon, epsilon = sn * beta, delta = -cs * beta;
  std::tie(cs, sn, gamma) =  math::sym_ortho(gamma, beta);
  tau = cs * phi, phi = sn * phi;

  // ----------------------
  // Update the solution:
  // 𝒘 ← (𝟣/(𝛽̅𝛾))𝒒̅ - (𝛿̅/𝛾)𝒘̅,
  // 𝒘 ← 𝒘 - (𝜀̅/𝛾)𝒘̿,
  // 𝒙 ← 𝒙 + 𝜏𝒘,
  // 𝒘̿ ← 𝒘̅, 𝒘̅ ← 𝒘.
  // ----------------------
  Blas::Sub(w_vec, qBarVec, 1.0 / (betaBar * gamma), wBarVec, deltaBar / gamma);
  Blas::Sub(w_vec, w_vec, wBarBarVec, epsilonBar / gamma);
  Blas::Add(x_vec, x_vec, w_vec, tau);
  std::swap(wBarBarVec, wBarVec), std::swap(wBarVec, w_vec);

  return phi;

} // MinresSolver::iterate

} // namespace Storm
#endif