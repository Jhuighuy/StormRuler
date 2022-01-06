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
/// Both right and left preconditioning is supported, left
/// preconditioning has slightly higher memory requirements and 
/// uses an additional preconditioning operator application per 
/// iteration.
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
public:
  stormPreconditionerSide PreSide = stormPreconditionerSide::Left;

private:
  stormReal_t alpha, rho, omega;
  tArray pArr, rArr, rTildeArr, tArr, vArr, zArr, sArr;

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
  stormUtils::AllocLike(xArr, pArr, rArr, rTildeArr, tArr, vArr);
  if (preOp != nullptr) {
    stormUtils::AllocLike(xArr, zArr);
    if (PreSide == stormPreconditionerSide::Left) {
      stormUtils::AllocLike(xArr, sArr);
    }
  }

  // ----------------------
  // 𝒓 ← 𝓐𝒙,
  // 𝒓 ← 𝒃 - 𝒓,
  // 𝒓̃ ← 𝒓.
  // ----------------------
  linOp.MatVec(rArr, xArr);
  stormBlas::Sub(rArr, bArr, rArr);
  stormBlas::Set(rTildeArr, rArr);

  return stormBlas::Norm2(rArr);

} // stormBiCgStabSolver<...>::Init

template<class tArray>
stormReal_t stormBiCgStabSolver<tArray>::Iterate(tArray& xArr,
                                                 tArray const& bArr,
                                                 stormOperator<tArray> const& linOp,
                                                 stormPreconditioner<tArray> const* preOp) {

  // ----------------------
  // Continue the iterations:
  // 𝜌̅ ← 𝜌,
  // 𝜌 ← <𝒓̃⋅𝒓>.
  // ----------------------
  stormReal_t rhoBar = rho; 
  rho = stormBlas::Dot(rTildeArr, rArr);

  // ----------------------
  // 𝗶𝗳 𝘍𝘪𝘳𝘴𝘵𝘐𝘵𝘦𝘳𝘢𝘵𝘪𝘰𝘯:
  //   𝒑 ← 𝒓.
  // 𝗲𝗹𝘀𝗲:
  //   𝛽 ← (𝜌/𝜌̅)⋅(𝛼/𝜔),
  //   𝒑 ← 𝒑 - 𝜔⋅𝒗,
  //   𝒑 ← 𝒓 + 𝛽⋅𝒑.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  bool const firstIteration = this->Iteration == 0;
  if (firstIteration) {
    stormBlas::Set(pArr, rArr);
  } else {
    stormReal_t const beta = 
      stormUtils::SafeDivide(rho, rhoBar)*stormUtils::SafeDivide(alpha, omega);
    stormBlas::Sub(pArr, pArr, vArr, omega);
    stormBlas::Add(pArr, rArr, pArr, beta);
  }

  // ----------------------
  // Update the solution and the residual:
  // 𝒗, 𝒛 ← 𝓐[𝓟]𝒑, [𝓟𝒑],
  // 𝛼 ← 𝜌/<𝒓̃⋅𝒗>,
  // 𝒙 ← 𝒙 + 𝛼⋅(𝓟 ≠ 𝗻𝗼𝗻𝗲 ? 𝒛 : 𝒑),
  // 𝒓 ← 𝒓 - 𝛼⋅𝒗.
  // ----------------------
  stormUtils::MatVecRightPre(vArr, zArr, pArr, linOp, preOp);
  alpha = stormUtils::SafeDivide(rho, stormBlas::Dot(rTildeArr, vArr));
  stormBlas::Add(xArr, xArr, (preOp != nullptr) ? zArr : pArr, alpha);
  stormBlas::Sub(rArr, rArr, vArr, alpha);

  /// @todo Check the residual norm here!
  //return stormBlas::Norm2(rArr);

  // ----------------------
  // Update the solution and the residual again:
  // 𝒕, 𝒛 ← 𝓐[𝓟]𝒓, [𝓟𝒓],
  // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
  //   𝒔 ← 𝓟𝒕,
  //   𝜔 ← <𝒔⋅𝒛>/<𝒔⋅𝒔>,
  // 𝗲𝗹𝘀𝗲:
  //   𝜔 ← <𝒕⋅𝒓>/<𝒕⋅𝒕>,
  // 𝗲𝗻𝗱 𝗶𝗳
  // 𝒙 ← 𝒙 + 𝜔⋅(𝓟 ≠ 𝗻𝗼𝗻𝗲 ? 𝒛 : 𝒓),
  // 𝒓 ← 𝒓 - 𝜔⋅𝒕.
  // ----------------------
  stormUtils::MatVecRightPre(tArr, zArr, rArr, linOp, preOp);
  bool const leftPre = (preOp != nullptr) && 
    (PreSide == stormPreconditionerSide::Left);
  if (leftPre) {
    preOp->MatVec(sArr, tArr);
    omega = stormUtils::SafeDivide(
      stormBlas::Dot(sArr, zArr), stormBlas::Dot(sArr, sArr));
  } else {
    omega = stormUtils::SafeDivide(
      stormBlas::Dot(tArr, rArr), stormBlas::Dot(tArr, tArr));
  }
  stormBlas::Add(xArr, xArr, (preOp != nullptr) ? zArr : rArr, omega);
  stormBlas::Sub(rArr, rArr, tArr, omega);

  return stormBlas::Norm2(rArr);

} // stormBiCgStabSolver<...>::Iterate

#endif // ifndef _STORM_SOLVER_BICGSTAB_
