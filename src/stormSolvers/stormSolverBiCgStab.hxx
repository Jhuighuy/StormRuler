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
#ifndef _STORM_SOLVER_BICGSTAB_
#define _STORM_SOLVER_BICGSTAB_

#include <stormSolvers/stormSolver.hxx>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Solve a linear operator equation [𝓟]𝓐𝒙 = [𝓟]𝒃, using \
///   the good old @c BiCGStab (Biconjugate Gradients Stabilized) method.
///
/// Residual norm, ‖𝒓‖, where 𝒓 = 𝒃 - 𝓐𝒙, is reported.
///
/// @c BiCGStab may be applied to the consistent singular problems,
/// it converges towards..
///
/// References:
/// @verbatim
/// [1] van der Vorst, Henk A.
///     “Bi-CGSTAB: A Fast and Smoothly Converging Variant of Bi-CG
///      for the Solution of Nonsymmetric Linear Systems.”
///     SIAM J. Sci. Comput. 13 (1992): 631-644.
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class tArray>
class stormBiCgStabSolver final : public stormIterativeSolver<tArray> {
private:
  stormReal_t alpha, beta, rho, omega;
  tArray pArr, rArr, rTildeArr, sArr, tArr, vArr, wArr, yArr, zArr;

  stormReal_t Init(tArray& xArr,
                   const tArray& bArr,
                   const stormOperator<tArray>& linOp,
                   const stormPreconditioner<tArray>* preOp) override;

  stormReal_t Iterate(tArray& xArr,
                      const tArray& bArr,
                      const stormOperator<tArray>& linOp,
                      const stormPreconditioner<tArray>* preOp) override;

}; // class stormBiCgStabSolver<...>

template<class tArray>
stormReal_t stormBiCgStabSolver<tArray>::Init(tArray& xArr,
                                              const tArray& bArr,
                                              const stormOperator<tArray>& linOp,
                                              const stormPreconditioner<tArray>* preOp) {

  // ----------------------
  // Allocate the intermediate arrays:
  // ----------------------
  stormUtils::AllocLike(xArr, pArr, rArr, rTildeArr, sArr, tArr, vArr);
  if (preOp != nullptr) {
    stormUtils::AllocLike(xArr, wArr, yArr, zArr);
  }

  // ----------------------
  // 𝒓 ← 𝓐𝒙,
  // 𝒓 ← 𝒃 - 𝒓,
  // 𝛿 ← ‖𝒓‖,
  // ----------------------
  linOp.MatVec(rArr, xArr);
  stormUtils::Sub(rArr, bArr, rArr);
  const stormReal_t delta = stormUtils::Norm2(rArr);

  // ----------------------
  // 𝒓̃ ← 𝒓,
  // 𝒑 ← {𝟢}ᵀ, 𝒗 ← {𝟢}ᵀ,
  // 𝜌 ← 𝟣, 𝛼 ← 𝟣, 𝜔 ← 𝟣.
  // ----------------------
  stormUtils::Set(rTildeArr, rArr);
  stormUtils::Fill(pArr, 0.0);
  stormUtils::Fill(vArr, 0.0);
  rho = 1.0, alpha = 1.0, omega = 1.0;

  return delta;

} // stormBiCgStabSolver<...>::Init

template<class tArray>
stormReal_t stormBiCgStabSolver<tArray>::Iterate(tArray& xArr,
                                                 const tArray& bArr,
                                                 const stormOperator<tArray>& linOp,
                                                 const stormPreconditioner<tArray>* preOp) {

  // ----------------------
  // Iterate:
  // 𝜌̂ ← 𝜌,
  // 𝜌 ← <𝒓̃⋅𝒓>,
  // 𝛽 ← (𝜌/𝜌̂)⋅(𝛼/𝜔),
  // ----------------------
  stormReal_t rhoHat = rho; 
  rho = stormUtils::Dot(rTildeArr, rArr);
  beta = stormUtils::SafeDivide(rho, rhoHat)*stormUtils::SafeDivide(alpha, omega);

  // ----------------------
  // 𝒑 ← 𝒑 - 𝜔𝒗,
  // 𝒑 ← 𝒓 + 𝛽𝒑,
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝒚 ← 𝓟𝒑,
  //   𝒗 ← 𝓐𝒚.
  // 𝗲𝗹𝘀𝗲:
  //   𝒗 ← 𝓐𝒑.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  stormUtils::Sub(pArr, pArr, vArr, omega);
  stormUtils::Add(pArr, rArr, pArr, beta);
  if (preOp != nullptr) {
    preOp->MatVec(yArr, pArr);
    linOp.MatVec(vArr, yArr);
  } else {
    linOp.MatVec(vArr, pArr);
  }

  // ----------------------
  // 𝛼 ← 𝜌/<𝒓̃⋅𝒗>,
  // 𝒔 ← 𝒓 - 𝛼𝒗,
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝒛 ← 𝓟𝒔,
  //   𝒕 ← 𝓐𝒛.
  // 𝗲𝗹𝘀𝗲:
  //   𝒕 ← 𝓐𝒔.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  alpha = stormUtils::SafeDivide(rho, stormUtils::Dot(rTildeArr, vArr));
  stormUtils::Sub(sArr, rArr, vArr, alpha);
  if (preOp != nullptr) {
    preOp->MatVec(zArr, sArr);
    linOp.MatVec(tArr, zArr);
  } else {
    linOp.MatVec(tArr, sArr);
  }

  // ----------------------
  // Update the solution:
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝒘 ← 𝓟𝒕,
  //   𝜔 ← <𝒘⋅𝒛>/<𝒘⋅𝒘>,
  //   𝒙 ← 𝒙 + 𝜔𝒛,
  //   𝒙 ← 𝒙 + 𝛼𝒚,
  // 𝗲𝗹𝘀𝗲:
  //   𝜔 ← <𝒕⋅𝒔>/<𝒕⋅𝒕>,
  //   𝒙 ← 𝒙 + 𝜔𝒔,
  //   𝒙 ← 𝒙 + 𝛼𝒑,
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  if (preOp != nullptr) {
    preOp->MatVec(wArr, tArr);
    omega = stormUtils::SafeDivide(
      stormUtils::Dot(wArr, zArr), stormUtils::Dot(wArr, wArr));
    stormUtils::Add(xArr, xArr, zArr, omega);
    stormUtils::Add(xArr, xArr, yArr, alpha);
  } else {
    omega = stormUtils::SafeDivide(
      stormUtils::Dot(tArr, sArr), stormUtils::Dot(tArr, tArr));
    stormUtils::Add(xArr, xArr, sArr, omega);
    stormUtils::Add(xArr, xArr, pArr, alpha);
  }

  // ----------------------
  // Update residual:
  // 𝒓 ← 𝒔 - 𝜔𝒕,
  // 𝛿 ← ‖𝒓‖.
  // ----------------------
  stormUtils::Sub(rArr, sArr, tArr, omega);
  const stormReal_t delta = stormUtils::Norm2(rArr);

  return delta;

} // stormBiCgStabSolver<...>::Iterate

#endif // ifndef _STORM_SOLVER_BICGSTAB_
