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
/// @c MINRES can be applied to the singular problems, and the self-adjoint
/// least squares problems: ‖[𝓜](𝓐[𝓜ᵀ]𝒚 - 𝒃)‖₂ → 𝘮𝘪𝘯, 𝒙 = [𝓜ᵀ]𝒚,
/// although convergeance to minimum norm solution is not guaranteed.
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
template<class tArray>
class stormMinresSolver final : public stormIterativeSolver<tArray> {
private:
  stormReal_t alpha, beta, betaBar, gamma, delta, deltaBar,
    epsilon, epsilonBar, tau, phi, phiTilde, cs, sn;
  tArray pArr, qArr, qBarArr,
    wArr, wBarArr, wBarBarArr, zArr, zBarArr, zBarBarArr;

  stormReal_t Init(tArray& xArr,
                   tArray const& bArr,
                   stormOperator<tArray> const& linOp,
                   stormPreconditioner<tArray> const* preOp) override;

  stormReal_t Iterate(tArray& xArr,
                      tArray const& bArr,
                      stormOperator<tArray> const& linOp,
                      stormPreconditioner<tArray> const* preOp) override;

}; // class stormMinresSolver<...>

template<class tArray>
stormReal_t stormMinresSolver<tArray>::Init(tArray& xArr,
                                            tArray const& bArr,
                                            stormOperator<tArray> const& linOp,
                                            stormPreconditioner<tArray> const* preOp) {
  assert(preOp != nullptr && "MINRES requires preconditioning for now.");

  // ----------------------
  // Allocate the intermediate arrays:
  // ----------------------
  stormUtils::AllocLike(xArr, pArr,
    wArr, wBarArr, wBarBarArr, zArr, zBarArr, zBarBarArr);
  if (preOp != nullptr) {
    stormUtils::AllocLike(xArr, qArr, qBarArr);
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
  stormBlas::Fill(wBarArr, 0.0);
  stormBlas::Fill(wBarBarArr, 0.0);
  linOp.MatVec(zBarArr, xArr);
  stormBlas::Sub(zBarArr, bArr, zBarArr);
  stormBlas::Fill(zBarBarArr, 0.0);
  if (preOp != nullptr) {
    preOp->MatVec(qArr, zBarArr);
  } else {
    //qArr = zBarArr
  }
  betaBar = 1.0; beta = std::sqrt(stormBlas::Dot(qArr, zBarArr));
  phi = beta; delta = 0.0; epsilon = 0.0;
  cs = -1.0; sn = 0.0;

  return phi;

} // stormMinresSolver<...>::Init

template<class tArray>
stormReal_t stormMinresSolver<tArray>::Iterate(tArray& xArr,
                                               tArray const& bArr,
                                               stormOperator<tArray> const& linOp,
                                               stormPreconditioner<tArray> const* preOp) {
  assert(preOp != nullptr && "MINRES requires preconditioning for now.");

  // ----------------------
  // Continue the Lanczos process:
  // 𝒑 ← 𝓐𝒒,
  // 𝛼 ← <𝒒⋅𝒑>/𝛽²,
  // 𝒛 ← (𝟣/𝛽)𝒑 - (𝛼/𝛽)𝒛̅,
  // 𝒛 ← 𝒛 - (𝛽/𝛽̅)𝒛̿,
  // 𝒒̅ ← 𝒒, 𝒒 ← [𝓟]𝒛,
  // 𝛽̅ ← 𝛽, 𝛽 ← √<𝒒⋅𝒛>,
  // 𝒛̿ ← 𝒛̅, 𝒛̅ ← 𝒛.
  // ----------------------
  linOp.MatVec(pArr, qArr);
  alpha = stormBlas::Dot(qArr, pArr)*std::pow(beta, -2);
  stormBlas::Sub(zArr, pArr, zBarArr, alpha/beta, 1.0/beta);
  stormBlas::Sub(zArr, zArr, zBarBarArr, beta/betaBar);
  if (preOp != nullptr) {
    std::swap(qBarArr, qArr);
    preOp->MatVec(qArr, zArr);
  } else {
    //qBarArr = qArr; qArr = zArr
  }
  betaBar = beta, beta = std::sqrt(stormBlas::Dot(qArr, zArr));
  std::swap(zBarBarArr, zBarArr), std::swap(zBarArr, zArr);

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
  // Update solution:
  // 𝒘 ← (𝟣/(𝛽̅𝛾))𝒒̅ - (𝛿̅/𝛾)𝒘̅,
  // 𝒘 ← 𝒘 - (𝜀̅/𝛾)𝒘̿,
  // 𝒙 ← 𝒙 + 𝜏𝒘,
  // 𝒘̿ ← 𝒘̅, 𝒘̅ ← 𝒘.
  // ----------------------
  stormBlas::Sub(wArr, qBarArr, wBarArr, deltaBar/gamma, 1.0/(betaBar*gamma));
  stormBlas::Sub(wArr, wArr, wBarBarArr, epsilonBar/gamma);
  stormBlas::Add(xArr, xArr, wArr, tau);
  std::swap(wBarBarArr, wBarArr), std::swap(wBarArr, wArr);

  return phi;

} // stormMinresSolver<...>::Iterate

#endif // ifndef _STORM_SOLVER_MINRES_HXX_
