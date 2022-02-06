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
#ifndef _STORM_SOLVER_MINRES_HXX_
#define _STORM_SOLVER_MINRES_HXX_

#include <cmath>

#include <stormSolvers/stormSolver.hxx>

namespace stormBlas {

  /// @brief Generate Givens rotation.
  inline auto SymOrtho(stormReal_t a, stormReal_t b) {

    // ----------------------
    // 𝑟𝑟 ← (𝑎² + 𝑏²)¹ᐟ²,
    // 𝑐𝑠 ← 𝑎/𝑟𝑟, 𝑠𝑛 ← 𝑏/𝑟𝑟.
    // ----------------------
    stormReal_t cs, sn, rr;
    rr = std::hypot(a, b);
    if (rr > 0.0) {
      cs = a/rr; sn = b/rr;
    } else {
      cs = 1.0; sn = 0.0;
    }

    return std::make_tuple(cs, sn, rr);

  } // SymOrtho

} // namespace stormBlas

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Solve a linear self-adjoint indefinite operator equation \
///   with the @c MINRES method.
///
/// @c MINRES can be applied to the singular problems, and the \
///   self-adjoint least squares problems, although convergeance to \
///   minimum norm solution is not guaranteed.
///
/// @note Despite 𝓐 may be indefinite, a positive-definite \
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
class stormMinresSolver final : public stormIterativeSolver<Vector> {
private:
  stormReal_t alpha, beta, betaBar, gamma, delta, deltaBar,
    epsilon, epsilonBar, tau, phi, phiTilde, cs, sn;
  Vector pVec, qVec, qBarVec,
    wVec, wBarVec, wBarBarVec, zVec, zBarVec, zBarBarVec;

  stormReal_t Init(Vector const& xVec,
                   Vector const& bVec,
                   stormOperator<Vector> const& linOp,
                   stormPreconditioner<Vector> const* preOp) override;

  stormReal_t Iterate(Vector& xVec,
                      Vector const& bVec,
                      stormOperator<Vector> const& linOp,
                      stormPreconditioner<Vector> const* preOp) override;

}; // class stormMinresSolver<...>

template<class Vector>
stormReal_t stormMinresSolver<Vector>::Init(Vector const& xVec,
                                            Vector const& bVec,
                                            stormOperator<Vector> const& linOp,
                                            stormPreconditioner<Vector> const* preOp) {

  assert(preOp != nullptr && "MINRES requires preconditioning for now.");

  pVec.Assign(xVec, false);
  wVec.Assign(xVec, false); 
  wBarVec.Assign(xVec, false); 
  wBarBarVec.Assign(xVec, false); 
  zVec.Assign(xVec, false); 
  zBarVec.Assign(xVec, false); 
  zBarBarVec.Assign(xVec, false); 
  if (preOp != nullptr) {
    qVec.Assign(xVec, false); 
    qBarVec.Assign(xVec, false); 
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
  wBarVec.Fill(0.0);
  wBarBarVec.Fill(0.0);
  linOp.MatVec(zBarVec, xVec);
  zBarVec.Sub(bVec, zBarVec);
  zBarBarVec.Fill(0.0);
  if (preOp != nullptr) {
    preOp->MatVec(qVec, zBarVec);
  } else {
    //qVec = zBarVec
  }
  betaBar = 1.0; beta = std::sqrt(stormBlas::Dot(qVec, zBarVec));
  phi = beta; delta = 0.0; epsilon = 0.0;
  cs = -1.0; sn = 0.0;

  return phi;

} // stormMinresSolver<...>::Init

template<class Vector>
stormReal_t stormMinresSolver<Vector>::Iterate(Vector& xVec,
                                               Vector const& bVec,
                                               stormOperator<Vector> const& linOp,
                                               stormPreconditioner<Vector> const* preOp) {

  assert(preOp != nullptr && "MINRES requires preconditioning for now.");

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
  linOp.MatVec(pVec, qVec);
  alpha = stormBlas::Dot(qVec, pVec)*std::pow(beta, -2);
  zVec.Sub(pVec, 1.0/beta, zBarVec, alpha/beta);
  zVec.Sub(zVec, zBarBarVec, beta/betaBar);
  if (preOp != nullptr) {
    std::swap(qBarVec, qVec);
    preOp->MatVec(qVec, zVec);
  } else {
    //qBarVec = qVec; qVec = zVec
  }
  betaBar = beta, beta = std::sqrt(stormBlas::Dot(qVec, zVec));
  std::swap(zBarBarVec, zBarVec), std::swap(zBarVec, zVec);

  // ----------------------
  // Construct and apply rotations:
  // 𝛿̅ ← 𝑐𝑠⋅𝛿 + 𝑠𝑛⋅𝛼, 𝛾 ← 𝑠𝑛⋅𝛿 - 𝑐𝑠⋅𝛼,
  // 𝜀̅ ← 𝜀, 𝜀 ← 𝑠𝑛⋅𝛽, 𝛿 ← -𝑐𝑠⋅𝛽,
  // 𝑐𝑠, 𝑠𝑛, 𝛾 ← 𝘚𝘺𝘮𝘖𝘳𝘵𝘩𝘰(𝛾, 𝛽),
  // 𝜏 ← 𝑐𝑠⋅𝜑, 𝜑 ← 𝑠𝑛⋅𝜑.
  // ----------------------
  deltaBar = cs*delta + sn*alpha, gamma = sn*delta - cs*alpha;
  epsilonBar = epsilon, epsilon = sn*beta, delta = -cs*beta;
  std::tie(cs, sn, gamma) = stormBlas::SymOrtho(gamma, beta);
  tau = cs*phi, phi = sn*phi;

  // ----------------------
  // Update the solution:
  // 𝒘 ← (𝟣/(𝛽̅𝛾))𝒒̅ - (𝛿̅/𝛾)𝒘̅,
  // 𝒘 ← 𝒘 - (𝜀̅/𝛾)𝒘̿,
  // 𝒙 ← 𝒙 + 𝜏𝒘,
  // 𝒘̿ ← 𝒘̅, 𝒘̅ ← 𝒘.
  // ----------------------
  wVec.Sub(qBarVec, 1.0/(betaBar*gamma), wBarVec, deltaBar/gamma);
  wVec.Sub(wVec, wBarBarVec, epsilonBar/gamma);
  xVec.Add(wVec, tau);
  std::swap(wBarBarVec, wBarVec), std::swap(wBarVec, wVec);

  return phi;

} // stormMinresSolver<...>::Iterate

#endif // ifndef _STORM_SOLVER_MINRES_HXX_
