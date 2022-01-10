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

#include <stormSolvers/stormSolver.hxx>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Solve a non-singular operator equation \
///   equation with the @c CGS (Conjugate Gradients Squared) method.
///
/// References:
/// @verbatim
/// [1] Sonneveld, Peter. 
///     “CGS, A Fast Lanczos-Type Solver for Nonsymmetric Linear systems.” 
///     SIAM J. Sci. Stat. Comput., 10:36-52, 1989.
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class tArray>
class stormCgsSolver final : public stormIterativeSolver<tArray> {
public:
  stormPreconditionerSide PreSide = stormPreconditionerSide::Right;

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

  bool const leftPre = 
    (preOp != nullptr) && (PreSide == stormPreconditionerSide::Left);

  stormUtils::AllocLike(xArr, pArr, qArr, rArr, rTildeArr, uArr, vArr);

  // ----------------------
  // 𝒓 ← 𝓐𝒙,
  // 𝒓 ← 𝒃 - 𝒓,
  // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
  //   𝒖 ← 𝒓,
  //   𝒓 ← 𝓟𝒖,
  // 𝗲𝗻𝗱 𝗶𝗳
  // 𝒓̃ ← 𝒓.
  // ----------------------
  linOp.MatVec(rArr, xArr);
  stormBlas::Sub(rArr, bArr, rArr);
  if (leftPre) {
    std::swap(uArr, rArr);
    preOp->MatVec(rArr, uArr);
  }
  stormBlas::Set(rTildeArr, rArr);

  return stormBlas::Norm2(rArr);

} // stormCgsSolver<...>::Init

template<class tArray>
stormReal_t stormCgsSolver<tArray>::Iterate(tArray& xArr,
                                            tArray const& bArr,
                                            stormOperator<tArray> const& linOp,
                                            stormPreconditioner<tArray> const* preOp) {

  bool const leftPre = 
    (preOp != nullptr) && (PreSide == stormPreconditionerSide::Left);
  bool const rightPre = 
    (preOp != nullptr) && (PreSide == stormPreconditionerSide::Right);

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
  // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
  //   𝒗 ← 𝓟(𝒒 ← 𝓐𝒑),
  // 𝗲𝗹𝘀𝗲 𝗶𝗳 𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦:
  //   𝒗 ← 𝓐(𝒒 ← 𝓟𝒑),
  // 𝗲𝗹𝘀𝗲:
  //   𝒗 ← 𝓐𝒑,
  // 𝗲𝗻𝗱 𝗶𝗳
  // 𝛼 ← 𝜌/<𝒓̃⋅𝒗>,
  // 𝒒 ← 𝒖 - 𝛼⋅𝒗,
  // 𝒗 ← 𝒖 + 𝒒.
  // ----------------------
  if (leftPre) {
    stormBlas::MatVec(vArr, *preOp, qArr, linOp, pArr);
  } else if (rightPre) {
    stormBlas::MatVec(vArr, linOp, qArr, *preOp, pArr);
  } else {
    linOp.MatVec(vArr, pArr);
  }
  stormReal_t const alpha = 
    stormUtils::SafeDivide(rho, stormBlas::Dot(rTildeArr, vArr));
  stormBlas::Sub(qArr, uArr, vArr, alpha);
  stormBlas::Add(vArr, uArr, qArr);

  // ----------------------
  // Update the solution and the residual:
  // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
  //   𝒙 ← 𝒙 + 𝛼⋅𝒗,
  //   𝒗 ← 𝓟(𝒖 ← 𝓐𝒗),
  //   𝒓 ← 𝒓 - 𝛼⋅𝒗.
  // 𝗲𝗹𝘀𝗲 𝗶𝗳 𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦:
  //   𝒗 ← 𝓐(𝒖 ← 𝓟𝒗),
  //   𝒙 ← 𝒙 + 𝛼⋅𝒖,
  //   𝒓 ← 𝒓 - 𝛼⋅𝒗.
  // 𝗲𝗹𝘀𝗲:
  //   𝒖 ← 𝓐𝒗,
  //   𝒙 ← 𝒙 + 𝛼⋅𝒗,
  //   𝒓 ← 𝒓 - 𝛼⋅𝒖.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  if (leftPre) {
    stormBlas::Add(xArr, xArr, vArr, alpha);
    stormBlas::MatVec(vArr, *preOp, uArr, linOp, vArr);
    stormBlas::Sub(rArr, rArr, vArr, alpha);
  } else if (rightPre) {
    stormBlas::MatVec(vArr, linOp, uArr, *preOp, vArr);
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
