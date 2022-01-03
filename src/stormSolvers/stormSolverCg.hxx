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
#ifndef _STORM_SOLVER_CG_HXX_
#define _STORM_SOLVER_CG_HXX_

#include <cmath>

#include <stormSolvers/stormSolver.hxx>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Solve a linear self-adjoint definite operator equation \
///   equation with the @c CG (Conjugate Gradients) method.
///
/// @c CG may be applied to the consistent singular problems,
/// it converges towards..
///
/// References:
/// @verbatim
/// [1] Hestenes, Magnus R. and Eduard Stiefel.
///     “Methods of conjugate gradients for solving linear systems.”
///     Journal of research of the National Bureau of Standards 49 (1952): 409-435.
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class tArray>
class stormCgSolver final : public stormIterativeSolver<tArray> {
private:
  stormReal_t alpha;
  tArray pArr, rArr, zArr;

  stormReal_t Init(tArray& xArr,
                   tArray const& bArr,
                   stormOperator<tArray> const& linOp,
                   stormPreconditioner<tArray> const* preOp) override;

  stormReal_t Iterate(tArray& xArr,
                      tArray const& bArr,
                      stormOperator<tArray> const& linOp,
                      stormPreconditioner<tArray> const* preOp) override;

}; // class stormCgSolver<...>

template<class tArray>
stormReal_t stormCgSolver<tArray>::Init(tArray& xArr,
                                        tArray const& bArr,
                                        stormOperator<tArray> const& linOp,
                                        stormPreconditioner<tArray> const* preOp) {

  // ----------------------
  // Allocate the intermediate arrays:
  // ----------------------
  stormUtils::AllocLike(xArr, pArr, rArr, zArr);

  // ----------------------
  // Initialize:
  // 𝒓 ← 𝓐𝒙,
  // 𝒓 ← 𝒃 - 𝒓.
  // ----------------------
  linOp.MatVec(rArr, xArr);
  stormBlas::Sub(rArr, bArr, rArr);

  // ----------------------
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝒛 ← 𝓟𝒓,
  //   𝒑 ← 𝒛,
  //   𝛼 ← <𝒓⋅𝒛>,
  // 𝗲𝗹𝘀𝗲:
  //   𝒑 ← 𝒓,
  //   𝛼 ← <𝒓⋅𝒓>.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  if (preOp != nullptr) {
    preOp->MatVec(zArr, rArr);
    stormBlas::Set(pArr, zArr);
    alpha = stormBlas::Dot(rArr, zArr);
  } else {
    stormBlas::Set(pArr, rArr);
    alpha = stormBlas::Dot(rArr, rArr);
  }

  return (preOp != nullptr) ? stormBlas::Norm2(rArr) : std::sqrt(alpha);

} // stormCgSolver<...>::Init

template<class tArray>
stormReal_t stormCgSolver<tArray>::Iterate(tArray& xArr,
                                           tArray const& bArr,
                                           stormOperator<tArray> const& linOp,
                                           stormPreconditioner<tArray> const* preOp) {

  // ----------------------
  // Iterate:
  // 𝒛 ← 𝓐𝒑,
  // 𝛼̅ ← 𝛼,
  // 𝛼 ← 𝛼/<𝒑⋅𝒛>,
  // 𝒙 ← 𝒙 + 𝛼⋅𝒑,
  // 𝒓 ← 𝒓 - 𝛼⋅𝒛,
  // ----------------------
  linOp.MatVec(zArr, pArr);
  stormReal_t const alphaBar = alpha;
  stormUtils::SafeDivideEquals(alpha, stormBlas::Dot(pArr, zArr));
  stormBlas::Add(xArr, xArr, pArr, alpha);
  stormBlas::Sub(rArr, rArr, zArr, alpha);

  // ----------------------
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝒛 ← 𝓟𝒓,
  //   𝛼 ← <𝒓⋅𝒛>,
  // 𝗲𝗹𝘀𝗲:
  //   𝛼 ← <𝒓⋅𝒓>.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  if (preOp != nullptr) {
    preOp->MatVec(zArr, rArr);
    alpha = stormBlas::Dot(rArr, zArr);
  } else {
    alpha = stormBlas::Dot(rArr, rArr);
  }

  // ----------------------
  // 𝛽 ← 𝛼/𝛼̅,
  // 𝒑 ← (𝓟 ≠ 𝗻𝗼𝗻𝗲 ? 𝒛 : 𝒓) + 𝛽⋅𝒑.
  // ----------------------
  stormReal_t const beta = stormUtils::SafeDivide(alpha, alphaBar);
  stormBlas::Add(pArr, (preOp != nullptr ? zArr : rArr), pArr, beta);

  return (preOp != nullptr) ? stormBlas::Norm2(rArr) : std::sqrt(alpha);

} // stormCgSolver<...>::Iterate

#endif // ifndef _STORM_SOLVER_CG_HXX_
