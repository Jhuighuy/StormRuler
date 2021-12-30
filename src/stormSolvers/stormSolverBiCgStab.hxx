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
/// @brief Solve a linear operator equation with the good old \
///   @c BiCGStab (Biconjugate Gradients Stabilized) method.
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
  tArray pArr, rArr, rTildeArr, sArr, tArr, vArr, yArr, zArr;

  stormReal_t Init(tArray& xArr,
                   tArray const& bArr,
                   stormOperator<tArray> const& linOp,
                   stormPreconditioner<tArray> const* preOp) override;

  stormReal_t Iterate(tArray& xArr,
                      tArray const& bArr,
                      stormOperator<tArray> const& linOp,
                      stormPreconditioner<tArray> const* preOp) override;

}; // class stormBiCgStabSolver<...>

template<class tArray>
stormReal_t stormBiCgStabSolver<tArray>::Init(tArray& xArr,
                                              tArray const& bArr,
                                              stormOperator<tArray> const& linOp,
                                              stormPreconditioner<tArray> const* preOp) {

  // ----------------------
  // Allocate the intermediate arrays:
  // ----------------------
  stormUtils::AllocLike(xArr, pArr, rArr, rTildeArr, sArr, tArr, vArr);
  if (preOp != nullptr) {
    stormUtils::AllocLike(xArr, yArr, zArr);
  }

  // ----------------------
  // 𝒓 ← 𝓐𝒙,
  // 𝒓 ← 𝒃 - 𝒓,
  // 𝜑 ← ‖𝒓‖,
  // ----------------------
  linOp.MatVec(rArr, xArr);
  stormBlas::Sub(rArr, bArr, rArr);
  const stormReal_t phi = stormBlas::Norm2(rArr);

  // ----------------------
  // 𝒓̃ ← 𝒓,
  // 𝒑 ← {𝟢}ᵀ, 𝒗 ← {𝟢}ᵀ,
  // 𝜌 ← 𝟣, 𝛼 ← 𝟣, 𝜔 ← 𝟣.
  // ----------------------
  stormBlas::Set(rTildeArr, rArr);
  stormBlas::Fill(pArr, 0.0);
  stormBlas::Fill(vArr, 0.0);
  rho = 1.0, alpha = 1.0, omega = 1.0;

  return phi;

} // stormBiCgStabSolver<...>::Init

template<class tArray>
stormReal_t stormBiCgStabSolver<tArray>::Iterate(tArray& xArr,
                                                 tArray const& bArr,
                                                 stormOperator<tArray> const& linOp,
                                                 stormPreconditioner<tArray> const* preOp) {

  // ----------------------
  // Iterate:
  // 𝜌̂ ← 𝜌,
  // 𝜌 ← <𝒓̃⋅𝒓>,
  // 𝛽 ← (𝜌/𝜌̂)⋅(𝛼/𝜔),
  // ----------------------
  stormReal_t rhoHat = rho; 
  rho = stormBlas::Dot(rTildeArr, rArr);
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
  stormBlas::Sub(pArr, pArr, vArr, omega);
  stormBlas::Add(pArr, rArr, pArr, beta);
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
  alpha = stormUtils::SafeDivide(rho, stormBlas::Dot(rTildeArr, vArr));
  stormBlas::Sub(sArr, rArr, vArr, alpha);
  if (preOp != nullptr) {
    preOp->MatVec(zArr, sArr);
    linOp.MatVec(tArr, zArr);
  } else {
    linOp.MatVec(tArr, sArr);
  }

  // ----------------------
  // Update the solution:
  // 𝜔 ← <𝒕⋅𝒔>/<𝒕⋅𝒕>,
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝒙 ← 𝒙 + 𝜔𝒛,
  //   𝒙 ← 𝒙 + 𝛼𝒚,
  // 𝗲𝗹𝘀𝗲:
  //   𝒙 ← 𝒙 + 𝜔𝒔,
  //   𝒙 ← 𝒙 + 𝛼𝒑,
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  omega = stormUtils::SafeDivide(
    stormBlas::Dot(tArr, sArr), stormBlas::Dot(tArr, tArr));
  if (preOp != nullptr) {
    stormBlas::Add(xArr, xArr, zArr, omega);
    stormBlas::Add(xArr, xArr, yArr, alpha);
  } else {
    stormBlas::Add(xArr, xArr, sArr, omega);
    stormBlas::Add(xArr, xArr, pArr, alpha);
  }

  // ----------------------
  // Update residual:
  // 𝒓 ← 𝒔 - 𝜔𝒕,
  // 𝜑 ← ‖𝒓‖.
  // ----------------------
  stormBlas::Sub(rArr, sArr, tArr, omega);
  stormReal_t const phi = stormBlas::Norm2(rArr);

  return phi;

} // stormBiCgStabSolver<...>::Iterate

#endif // ifndef _STORM_SOLVER_BICGSTAB_
