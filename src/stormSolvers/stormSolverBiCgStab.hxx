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
  stormReal_t alpha, rho, omega;
  tArray pArr, rArr, rTildeArr, tArr, vArr, yArr, zArr;

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
  // ----------------------
  stormBlas::Set(rTildeArr, rArr);

  return phi;

} // stormBiCgStabSolver<...>::Init

template<class tArray>
stormReal_t stormBiCgStabSolver<tArray>::Iterate(tArray& xArr,
                                                 tArray const& bArr,
                                                 stormOperator<tArray> const& linOp,
                                                 stormPreconditioner<tArray> const* preOp) {

  // ----------------------
  // Iterate:
  // 𝜌̅ ← 𝜌,
  // 𝜌 ← <𝒓̃⋅𝒓>.
  // ----------------------
  stormReal_t rhoBar = rho; 
  rho = stormBlas::Dot(rTildeArr, rArr);

  // ----------------------
  // 𝗶𝗳 𝘍𝘪𝘳𝘴𝘵𝘐𝘵𝘦𝘳𝘢𝘵𝘪𝘰𝘯:
  //   𝒑 ← 𝒓,
  // 𝗲𝗹𝘀𝗲:
  //   𝛽 ← (𝜌/𝜌̅)⋅(𝛼/𝜔),
  //   𝒑 ← 𝒑 - 𝜔⋅𝒗,
  //   𝒑 ← 𝒓 + 𝛽⋅𝒑,
  // 𝗲𝗻𝗱 𝗶𝗳
  // 𝒗, 𝒚 ← 𝓐[𝓟]𝒑, [𝓟𝒑].
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
  stormUtils::MatVecRightPre(vArr, yArr, pArr, linOp, preOp);

  // ----------------------
  // 𝛼 ← 𝜌/<𝒓̃⋅𝒗>,
  // 𝒓 ← 𝒓 - 𝛼⋅𝒗,
  // 𝒕, 𝒛 ← 𝓐[𝓟]𝒓, [𝓟𝒓].
  // ----------------------
  alpha = stormUtils::SafeDivide(rho, stormBlas::Dot(rTildeArr, vArr));
  stormBlas::Sub(rArr, rArr, vArr, alpha);
  stormUtils::MatVecRightPre(tArr, zArr, rArr, linOp, preOp);

  // ----------------------
  // Update the solution:
  // 𝜔 ← <𝒕⋅𝒓>/<𝒕⋅𝒕>,
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝒙 ← 𝒙 + 𝜔⋅𝒛,
  //   𝒙 ← 𝒙 + 𝛼⋅𝒚,
  // 𝗲𝗹𝘀𝗲:
  //   𝒙 ← 𝒙 + 𝜔⋅𝒓,
  //   𝒙 ← 𝒙 + 𝛼⋅𝒑,
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  omega = stormUtils::SafeDivide(
    stormBlas::Dot(tArr, rArr), stormBlas::Dot(tArr, tArr));
  if (preOp != nullptr) {
    stormBlas::Add(xArr, xArr, zArr, omega);
    stormBlas::Add(xArr, xArr, yArr, alpha);
  } else {
    stormBlas::Add(xArr, xArr, rArr, omega);
    stormBlas::Add(xArr, xArr, pArr, alpha);
  }

  // ----------------------
  // Update residual:
  // 𝒓 ← 𝒓 - 𝜔⋅𝒕,
  // 𝜑 ← ‖𝒓‖.
  // ----------------------
  stormBlas::Sub(rArr, rArr, tArr, omega);
  stormReal_t const phi = stormBlas::Norm2(rArr);

  return phi;

} // stormBiCgStabSolver<...>::Iterate

#endif // ifndef _STORM_SOLVER_BICGSTAB_
