/// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ///
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
/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ///
#ifndef STORM_RULER_SOLVER_CG_HXX_
#define STORM_RULER_SOLVER_CG_HXX_

#include <StormRuler_Solver.hxx>

#include <cmath>

/// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ///
/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ///

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Solve a linear self-adjoint definite operator equation \
/// [𝓜]𝓐[𝓜ᵀ]𝒚 = [𝓜]𝒃, 𝒙 = [𝓜ᵀ]𝒚, [𝓜𝓜ᵀ = 𝓟], using the \
/// Conjugate Gradients (CG) method.
///
/// CG may be applied to the consistent singular problems, 
/// it converges towards..
///
/// References:
/// @verbatim
/// [1] Hestenes, Magnus R. and Eduard Stiefel. 
///     “Methods of conjugate gradients for solving linear systems.” 
///     Journal of research of the National Bureau of Standards 49 (1952): 409-435.
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class tArray, class tOperator>
class stormCgSolver final : public stormIterativeSolver<tArray, tOperator> {
private:
  stormReal_t alpha, beta, gamma, delta;
  tArray pArr, rArr, tArr, zArr;

protected:

  /// @brief Initialize the CG iterative solver.
  ///
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param linOp Self-adjoint linear operator, 𝓐(𝒙).
  /// @param preOp Self-adjoint linear preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Preconditioned residual norm, \
  /// square root of <𝒓⋅𝒛>, where 𝒓 = 𝓐𝒙 - 𝒃 and 𝒛 = [𝓟]𝒓.
  stormReal_t Init(tArray& xArr,
                   const tArray& bArr,
                   const tOperator& linOp,
                   const tOperator* preOp) override final;

  /// @brief Iterate the CG solver.
  ///
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param linOp Self-adjoint linear operator, 𝓐(𝒙).
  /// @param preOp Self-adjoint linear preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Preconditioned residual norm, \
  /// square root of <𝒓⋅𝒛>, where 𝒓 = 𝓐𝒙 - 𝒃 and 𝒛 = [𝓟]𝒓.
  stormReal_t Iterate(tArray& xArr,
                      const tArray& bArr,
                      const tOperator& linOp,
                      const tOperator* preOp) override final;

}; // class stormCgSolver<...>

template<class tArray, class tOperator>
stormReal_t stormCgSolver<tArray, tOperator>::Init(tArray& xArr,
                                                   const tArray& bArr,
                                                   const tOperator& linOp,
                                                   const tOperator* preOp) {
  // ----------------------
  // Allocate the intermediate arrays. 
  // ----------------------
  stormUtils::AllocLike(xArr, pArr, rArr, tArr);
  if (preOp != nullptr) {
    stormUtils::AllocLike(rArr, zArr);
  }

  // ----------------------
  // Initialize:
  // 𝒓 ← 𝓐𝒙,
  // 𝒓 ← 𝒃 - 𝒕.
  // ----------------------
  linOp.MatVec(rArr, xArr);
  stormUtils::Sub(rArr, bArr, rArr);

  // ----------------------
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝒛 ← 𝓟𝒓,
  //   𝒑 ← 𝒛,
  //   𝛾 ← <𝒓⋅𝒛>,
  // 𝗲𝗹𝘀𝗲:
  //   𝒑 ← 𝒓, 
  //   𝛾 ← <𝒓⋅𝒓>.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  if (preOp != nullptr) {
    preOp->MatVec(zArr, rArr);
    stormUtils::Set(pArr, zArr);
    gamma = stormUtils::Dot(rArr, zArr);
  } else {
    stormUtils::Set(pArr, rArr);
    gamma = stormUtils::Dot(rArr, rArr);
  }

  return std::sqrt(gamma);

} // stormCgSolver<...>::Init

template<class tArray, class tOperator>
stormReal_t stormCgSolver<tArray, tOperator>::Iterate(tArray& xArr,
                                                      const tArray& bArr,
                                                      const tOperator& linOp,
                                                      const tOperator* preOp) {
  // ----------------------
  // 𝒕 ← 𝓐𝒑,
  // 𝛼 ← 𝛾/<𝒑⋅𝒕>,
  // 𝒙 ← 𝒙 + 𝛼𝒑,
  // 𝒓 ← 𝒓 - 𝛼𝒕,
  // ----------------------
  linOp.MatVec(tArr, pArr);
  alpha = stormUtils::SafeDivide(gamma, stormUtils::Dot(pArr, tArr));
  stormUtils::Add(xArr, xArr, pArr, alpha);
  stormUtils::Sub(rArr, rArr, tArr, alpha);

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
    alpha = stormUtils::Dot(rArr, zArr);
  } else {
    alpha = stormUtils::Dot(rArr, rArr);
  }

  // ----------------------
  // 𝛽 ← 𝛼/𝛾,
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝒑 ← 𝒛 + 𝛽𝒑,
  // 𝗲𝗹𝘀𝗲:
  //   𝒑 ← r + 𝛽𝒑,
  // 𝗲𝗻𝗱 𝗶𝗳  
  // 𝛾 ← 𝛼.
  // ----------------------
  beta = stormUtils::SafeDivide(alpha, gamma);
  stormUtils::Add(pArr, (preOp != nullptr ? zArr : rArr), pArr, beta);
  gamma = alpha;

  return std::sqrt(gamma);

} // stormCgSolver<...>::Iterate

/// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ///
/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ///

template<typename stormMatVecFuncT_t>
STORM_INL void stormLinSolve2(stormMesh_t mesh,
                              stormString_t method,
                              stormString_t preMethod,
                              stormArray_t x,
                              stormArray_t b,
                              stormMatVecFuncT_t matVec) {
  stormArray xx = {mesh, x}, bb = {mesh, b};
  stormLinearOperator<stormArray> op {
    [&](stormArray& yy, const stormArray& xx) {
      matVec(yy.Mesh, yy.Array, xx.Array);
    }
  };

  stormCgSolver<stormArray, stormLinearOperator<stormArray>> cgSolver;
  cgSolver.Solve(xx, bb, op);

} // stormLinSolve

#endif // ifndef STORM_RULER_SOLVER_CG_HXX_
