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
#ifndef _STORM_SOLVER_CGS_HXX_
#define _STORM_SOLVER_CGS_HXX_

#include <cmath>

#include <stormSolvers/stormSolver.hxx>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Solve a non-singular operator equation \
///   equation with the @c CGS (Conjugate Gradients Squared) method.
///
/// References:
/// @verbatim
/// [1] ???
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class tArray>
class stormCgsSolver final : public stormIterativeSolver<tArray> {
private:
  stormReal_t rho;
  tArray pArr, qArr, rArr, rTildeArr, uArr, vArr;

  stormReal_t Init(tArray& xArr,
                   tArray const& bArr,
                   stormOperator<tArray> const& linOp,
                   stormPreconditioner<tArray> const* preOp) override;

  stormReal_t Iterate(tArray& xArr,
                      tArray const& bArr,
                      stormOperator<tArray> const& linOp,
                      stormPreconditioner<tArray> const* preOp) override;

}; // class stormCgsSolver<...>

template<class tArray>
stormReal_t stormCgsSolver<tArray>::Init(tArray& xArr,
                                         tArray const& bArr,
                                         stormOperator<tArray> const& linOp,
                                         stormPreconditioner<tArray> const* preOp) {

  // ----------------------
  // Allocate the intermediate arrays:
  // ----------------------
  stormUtils::AllocLike(xArr, pArr, qArr, rArr, rTildeArr, uArr, vArr);

  // ----------------------
  // 𝒓 ← 𝓐𝒙,
  // 𝒓 ← 𝒃 - 𝒓,
  // 𝒓̃ ← 𝒓.
  // ----------------------
  linOp.MatVec(rArr, xArr);
  stormBlas::Sub(rArr, bArr, rArr);
  stormBlas::Set(rTildeArr, rArr);

  return stormBlas::Norm2(rArr);

} // stormCgsSolver<...>::Init

template<class tArray>
stormReal_t stormCgsSolver<tArray>::Iterate(tArray& xArr,
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
  //   𝒖 ← 𝒓,
  //   𝒑 ← 𝒖.
  // 𝗲𝗹𝘀𝗲:
  //   𝛽 ← 𝜌/𝜌̅,
  //   𝒖 ← 𝒓 + 𝛽⋅𝒒,
  //   𝒑 ← 𝒒 + 𝛽⋅𝒑,
  //   𝒑 ← 𝒖 + 𝛽⋅𝒑.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  bool const firstIteration = this->Iteration == 0;
  if (firstIteration) {
    stormBlas::Set(uArr, rArr);
    stormBlas::Set(pArr, uArr);
  } else {
    stormReal_t const beta = 
      stormUtils::SafeDivide(rho, rhoBar);
    stormBlas::Add(uArr, rArr, qArr, beta);
    stormBlas::Add(pArr, qArr, pArr, beta);
    stormBlas::Add(pArr, uArr, pArr, beta);
  }

  // ----------------------
  // Update the solution and the residual:
  // 𝒗, 𝒒 ← 𝓐[𝓟]𝒑, [𝓟𝒑].
  // 𝛼 ← 𝜌/<𝒓̃⋅𝒗>,
  // 𝒒 ← 𝒖 - 𝛼⋅𝒗,
  // 𝒗 ← 𝒖 + 𝒒,
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝒗, 𝒖 ← 𝓐𝓟𝒗, 𝓟𝒗,
  //   𝒙 ← 𝒙 + 𝛼⋅𝒖,
  //   𝒓 ← 𝒓 - 𝛼⋅𝒗.
  // 𝗲𝗹𝘀𝗲:
  //   𝒖 ← 𝓐𝒗,
  //   𝒙 ← 𝒙 + 𝛼⋅𝒗,
  //   𝒓 ← 𝒓 - 𝛼⋅𝒖.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  stormUtils::MatVecRightPre(vArr, qArr, pArr, linOp, preOp);
  stormReal_t const alpha = 
    stormUtils::SafeDivide(rho, stormBlas::Dot(rTildeArr, vArr));
  stormBlas::Sub(qArr, uArr, vArr, alpha);
  stormBlas::Add(vArr, uArr, qArr);
  if (preOp != nullptr) {
    stormUtils::MatVecRightPre(vArr, uArr, vArr, linOp, preOp);
    stormBlas::Add(xArr, xArr, uArr, alpha);
    stormBlas::Sub(rArr, rArr, vArr, alpha);
  } else {
    linOp.MatVec(uArr, vArr);
    stormBlas::Add(xArr, xArr, vArr, alpha);
    stormBlas::Sub(rArr, rArr, uArr, alpha);
  }

  return stormBlas::Norm2(rArr);

} // stormCgsSolver<...>::Iterate

#endif // ifndef _STORM_SOLVER_CGS_HXX_
