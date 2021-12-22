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
#ifndef STORM_RULER_SOLVER_LSQR_HXX_
#define STORM_RULER_SOLVER_LSQR_HXX_

#include <StormRuler_Solver.hxx>

#include <cmath>

/// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ///
/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ///

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Solve a right preconditioned linear least squares problem \
///   ‖𝓐[𝓟]𝒚 - 𝒃‖₂ → 𝘮𝘪𝘯, 𝒙 = [𝓟]𝒚, using the @c LSQR method.
///
/// @c LSQR is algebraically equivalent to applying @c CG
/// to the normal equations: (𝓐[𝓟])*𝓐[𝓟]𝒚 = (𝓐[𝓟])*𝒃, 𝒙 = [𝓟]𝒚,
/// (or, equivalently, [𝓟*]𝓐*𝓐[𝓟]𝒚 = [𝓟*]𝓐*𝒃, 𝒙 = [𝓟]𝒚),
/// but has better numerical properties.
///
/// The residual norm ‖𝓐[𝓟]𝒚 - 𝒃‖₂ decreases monotonically, 
/// while the normal equation's residual norm ‖(𝓐[𝓟])*(𝓐[𝓟]𝒚 - 𝒃)‖ 
/// is not guaranteed to decrease.
///
/// @c LSQR is not recommended in the self-adjoint case,
/// please consider @c MINRES instead.
///
/// References:
/// @verbatim
/// [1] Paige, C. and M. Saunders. 
///     “LSQR: An Algorithm for Sparse Linear Equations and 
///     Sparse Least Squares.” ACM Trans. Math. Softw. 8 (1982): 43-71.
/// [2] Karimi, S., D. K. Salkuyeh and F. Toutounian. 
///     “A preconditioner for the LSQR algorithm.” 
///     Journal of applied mathematics & informatics 26 (2008): 213-222.
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class tArray, class tOperator>
class stormLsqrSolver final : public stormIterativeSolver<tArray, tOperator> {
private:
  stormReal_t alpha, beta, rho, rhoBar, theta, phi, phiBar, phiTilde, cs, sn;
  tArray sArr, tArr, rArr, uArr, vArr, wArr, zArr;

protected:

  /// @brief Initialize the @c LSQR solver.
  ///
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param linOp Linear operator, 𝓐(𝒙).
  /// @param preOp Linear preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Residual norm, ‖𝒓‖, where 𝒓 = 𝒃 - 𝓐𝒙.
  stormReal_t Init(tArray& xArr,
                   const tArray& bArr,
                   const tOperator& linOp,
                   const tOperator* preOp) override final;

  /// @brief Iterate the @c LSQR solver.
  ///
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param linOp Linear operator, 𝓐(𝒙).
  /// @param preOp Linear preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Residual norm, ‖𝒓‖, where 𝒓 = 𝒃 - 𝓐𝒙.
  stormReal_t Iterate(tArray& xArr,
                      const tArray& bArr,
                      const tOperator& linOp,
                      const tOperator* preOp) override final;

  /// @brief Finalize the @c LSQR iterations.
  ///
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  /// @param preOp Preconditioner operator, 𝓟(𝒙).
  void Finalize(tArray& xArr,
                const tArray& bArr,
                const tOperator& anyOp,
                const tOperator* preOp = nullptr) override final;

}; // class stormLsqrSolver<...>

template<class tArray, class tOperator>
stormReal_t stormLsqrSolver<tArray, tOperator>::Init(tArray& xArr,
                                                     const tArray& bArr,
                                                     const tOperator& linOp,
                                                     const tOperator* preOp) {
  // ----------------------
  // Allocate the intermediate arrays:
  // ----------------------
  stormUtils::AllocLike(xArr, tArr, rArr, uArr, vArr, wArr, zArr);
  if (preOp != nullptr) {
    stormUtils::AllocLike(xArr, sArr);
  }

  // ----------------------
  // Utilize the initial guess.
  // Consider the decomposition:
  // 𝒙 = 𝒙₀ + 𝒛. (*)
  // Substituting (*) into the equation, we get:
  // 𝓐[𝓟]𝒚 = 𝒓, where: 𝒛 = [𝓟]𝒚, 𝒓 = 𝒃 - 𝓐𝒙₀.
  // The last equations can be solved with 𝒚₀ = {𝟢}ᵀ.
  // ----------------------

  // ----------------------
  // Initialize:
  // 𝒓 ← 𝓐𝒙,
  // 𝒓 ← 𝒃 - 𝒓,
  // 𝛽 ← ‖𝒓‖, 𝒖 ← 𝒓/𝛽,
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲: 
  //   𝒔 ← 𝓐*𝒖, 𝒕 ← 𝓟*𝒔, 
  // 𝗲𝗹𝘀𝗲: 𝒕 ← 𝓐*𝒖, 𝗲𝗻𝗱 𝗶𝗳
  // 𝛼 ← ‖𝒕‖, 𝒗 ← 𝒕/𝛼.
  // ----------------------
  linOp.MatVec(rArr, xArr);
  stormUtils::Sub(rArr, bArr, rArr);
  beta = stormUtils::Norm2(rArr); stormUtils::Scale(uArr, rArr, 1.0/beta);
  if (preOp != nullptr) {
    linOp.ConjMatVec(sArr, uArr);
    preOp->ConjMatVec(tArr, sArr);
  } else {
    linOp.ConjMatVec(tArr, uArr);
  }
  alpha = stormUtils::Norm2(tArr); stormUtils::Scale(vArr, tArr, 1.0/alpha);

  // ----------------------
  // 𝜑̅ ← 𝛽, 𝜌̅ ← 𝛼.
  // 𝒛 ← {𝟢}ᵀ,
  // 𝒘 ← 𝒗,
  // ----------------------
  phiBar = beta; rhoBar = alpha;
  stormUtils::Fill(zArr, 0.0);
  stormUtils::Set(wArr, vArr);

  return phiBar;

} // stormLsqrSolver<...>::Init

template<class tArray, class tOperator>
stormReal_t stormLsqrSolver<tArray, tOperator>::Iterate(tArray& xArr,
                                                        const tArray& bArr,
                                                        const tOperator& linOp,
                                                        const tOperator* preOp) {
  // ----------------------
  // Continue the bidiagonalization:
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲: 
  //   𝒔 ← 𝓟𝒗, 𝒕 ← 𝓐𝒔,
  // 𝗲𝗹𝘀𝗲: 𝒕 ← 𝓐𝒗, 𝗲𝗻𝗱 𝗶𝗳
  // 𝒕 ← 𝒕 - 𝛼𝒖,
  // 𝛽 ← ‖𝒕‖, 𝒖 ← 𝒕/𝛽,
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝒔 ← 𝓐*𝒖, 𝒕 ← 𝓟*𝒔, 
  // 𝗲𝗹𝘀𝗲: 𝒕 ← 𝓐*𝒖, 𝗲𝗻𝗱 𝗶𝗳
  // 𝒕 ← 𝒕 - 𝛽𝒗,
  // 𝛼 ← ‖𝒕‖, 𝒗 ← 𝒕/𝛼.
  // ----------------------
  if (preOp != nullptr) {
    preOp->MatVec(sArr, vArr);
    linOp.MatVec(tArr, sArr);
  } else {
    linOp.MatVec(tArr, vArr);
  }
  stormUtils::Sub(tArr, tArr, uArr, alpha);
  beta = stormUtils::Norm2(tArr); stormUtils::Scale(uArr, tArr, 1.0/beta);
  if (preOp != nullptr) {
    linOp.ConjMatVec(sArr, uArr);
    preOp->ConjMatVec(tArr, sArr);
  } else {
    linOp.ConjMatVec(tArr, uArr);
  }
  stormUtils::Sub(tArr, tArr, vArr, beta);
  alpha = stormUtils::Norm2(tArr); stormUtils::Scale(vArr, tArr, 1.0/alpha);

  // ----------------------
  // Construct and apply rotation:
  // 𝜌 ← (𝜌̅² + 𝛽²)¹ᐟ²,
  // 𝑐𝑠 ← 𝜌̅/𝜌, 𝑠𝑛 ← 𝛽/𝜌,
  // 𝜃 ← 𝑠𝑛⋅𝛼, 𝜌̅ ← -𝑐𝑠⋅𝛼,
  // 𝜑 ← 𝑐𝑠⋅𝜑, 𝜑̅ ← 𝑠𝑛⋅𝜑̅.
  // ----------------------
  rho = std::hypot(rhoBar, beta);
  cs = rhoBar/rho; sn = beta/rho;
  theta = sn*alpha; rhoBar = -cs*alpha;
  phi = cs*phiBar; phiBar = sn*phiBar;

  // ----------------------
  // Update 𝒛-solution:
  // 𝒛 ← 𝒛 + (𝜑/𝜌)𝒘,
  // 𝒘 ← 𝒗 - (𝜃/𝜌)𝒘.
  // ----------------------
  stormUtils::Add(zArr, zArr, wArr, phi/rho);
  stormUtils::Sub(wArr, vArr, wArr, theta/rho);

  return phiBar;

} // stormLsqrSolver<...>::Iterate

template<class tArray, class tOperator>
void stormLsqrSolver<tArray, tOperator>::Finalize(tArray& xArr,
                                                  const tArray& bArr,
                                                  const tOperator& linOp,
                                                  const tOperator* preOp) {
  // ----------------------
  // Compute 𝒙-solution:
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝒕 ← 𝓟𝒛, 𝒙 ← 𝒙 + 𝒕.
  // 𝗲𝗹𝘀𝗲: 𝒙 ← 𝒙 + 𝒛. 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  if (preOp != nullptr) {
    preOp->MatVec(tArr, zArr);
    stormUtils::Add(xArr, xArr, tArr);
  } else {
    stormUtils::Add(xArr, xArr, zArr);
  }

} // stormLsqrSolver<...>::Finalize

/// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ///
/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ///

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Solve a right preconditioned linear least squares problem \
///   ‖𝓐[𝓟]𝒚 - 𝒃‖₂ → 𝘮𝘪𝘯, 𝒙 = [𝓟]𝒚, using the @c LSMR method.
///
/// @c LSMR is algebraically equivalent to applying @c MINRES 
/// to the normal equations: (𝓐[𝓟])*𝓐[𝓟]𝒚 = (𝓐[𝓟])*𝒃, 𝒙 = [𝓟]𝒚, 
/// (or, equivalently, [𝓟*]𝓐*𝓐[𝓟]𝒚 = [𝓟*]𝓐*𝒃, 𝒙 = [𝓟]𝒚),
/// but has better numerical properties.
/// 
/// The normal equation's residual norm ‖(𝓐[𝓟])*(𝓐[𝓟]𝒚 - 𝒃)‖ 
/// decreases monotonically, while the residual norm ‖𝓐[𝓟]𝒚 - 𝒃‖
/// is not guaranteed to decrease (but decreases on practice).
///
/// Using @c LSMR is not recommended in the self-adjoint case,
/// please consider @c MINRES instead.
///
/// References:
/// @verbatim
/// [1] Fong, D. C. and M. Saunders. 
///     “LSMR: An Iterative Algorithm for Sparse Least-Squares Problems.” 
///     SIAM J. Sci. Comput. 33 (2011): 2950-2971.
/// [2] Karimi, S., D. K. Salkuyeh and F. Toutounian. 
///     “A preconditioner for the LSQR algorithm.” 
///     Journal of applied mathematics & informatics 26 (2008): 213-222.
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class tArray, class tOperator>
class stormLsmrSolver final : public stormIterativeSolver<tArray, tOperator> {
private:
  stormReal_t alpha, alphaBar, beta, rho, rhoBar, cs, sn, \
    theta, thetaBar, psi, psiBar, psiTilde, zeta, csBar, snBar;
  tArray rArr, sArr, tArr, wArr, hArr, uArr, vArr, zArr;

protected:

  /// @brief Initialize the @c LSMR solver.
  ///
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param linOp Linear operator, 𝓐(𝒙).
  /// @param preOp Linear preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Normal equation residual norm, ‖(𝓐[𝓟])*𝒓‖, where 𝒓 = 𝒃 - 𝓐𝒙.
  stormReal_t Init(tArray& xArr,
                   const tArray& bArr,
                   const tOperator& linOp,
                   const tOperator* preOp) override final;

  /// @brief Iterate the @c LSMR solver.
  ///
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param linOp Linear operator, 𝓐(𝒙).
  /// @param preOp Linear preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Normal equation residual norm, ‖(𝓐[𝓟])*𝒓‖, where 𝒓 = 𝒃 - 𝓐𝒙.
  stormReal_t Iterate(tArray& xArr,
                      const tArray& bArr,
                      const tOperator& linOp,
                      const tOperator* preOp) override final;

  /// @brief Finalize the @c LSMR iterations.
  ///
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  /// @param preOp Preconditioner operator, 𝓟(𝒙).
  void Finalize(tArray& xArr,
                const tArray& bArr,
                const tOperator& anyOp,
                const tOperator* preOp = nullptr) override final;

}; // class stormLsmrSolver<...>

template<class tArray, class tOperator>
stormReal_t stormLsmrSolver<tArray, tOperator>::Init(tArray& xArr,
                                                     const tArray& bArr,
                                                     const tOperator& linOp,
                                                     const tOperator* preOp) {
  // ----------------------
  // Allocate the intermediate arrays:
  // ----------------------
  stormUtils::AllocLike(xArr, tArr, rArr, uArr, vArr, wArr, hArr, zArr);
  if (preOp != nullptr) {
    stormUtils::AllocLike(xArr, sArr);
  }

  // ----------------------
  // Utilize the initial guess.
  // Consider the decomposition:
  // 𝒙 = 𝒙₀ + 𝒛. (*)
  // Substituting (*) into the equation, we get:
  // 𝓐[𝓟]𝒚 = 𝒓, where: 𝒛 = [𝓟]𝒚, 𝒓 = 𝒃 - 𝓐𝒙₀.
  // The last equations can be solved with 𝒚₀ = {𝟢}ᵀ.
  // ----------------------

  // ----------------------
  // Initialize:
  // 𝒓 ← 𝓐𝒙,
  // 𝒓 ← 𝒃 - 𝒓,
  // 𝛽 ← ‖𝒓‖, 𝒖 ← 𝒓/𝛽,
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲: 
  //   𝒔 ← 𝓐*𝒖, 𝒕 ← 𝓟*𝒔, 
  // 𝗲𝗹𝘀𝗲: 𝒕 ← 𝓐*𝒖, 𝗲𝗻𝗱 𝗶𝗳
  // 𝛼 ← ‖𝒕‖, 𝒗 ← 𝒕/𝛼.
  // ----------------------
  linOp.MatVec(rArr, xArr);
  stormUtils::Sub(rArr, bArr, rArr);
  beta = stormUtils::Norm2(rArr); stormUtils::Scale(uArr, rArr, 1.0/beta);
  if (preOp != nullptr) {
    linOp.ConjMatVec(sArr, uArr);
    preOp->ConjMatVec(tArr, sArr);
  } else {
    linOp.ConjMatVec(tArr, uArr);
  }
  alpha = stormUtils::Norm2(tArr); stormUtils::Scale(vArr, tArr, 1.0/alpha);

  // ----------------------
  // 𝛼̅ ← 𝛼, 𝜓̅ ← 𝛼𝛽,
  // 𝜁 ← 𝟣, 𝑐̅𝑠̅ ← 𝟣, 𝑠̅𝑛̅ ← 𝟢,
  // 𝒛 ← {𝟢}ᵀ,
  // 𝒘 ← 𝒗, 𝒉 ← {𝟢}ᵀ.
  // ----------------------
  alphaBar = alpha; psiBar = alpha*beta;
  zeta = 1.0; csBar = 1.0; snBar = 0.0;
  stormUtils::Fill(zArr, 0.0);
  stormUtils::Set(wArr, vArr); stormUtils::Fill(hArr, 0.0);

  return std::abs(psiBar);

} // stormLsmrSolver<...>::Init

template<class tArray, class tOperator>
stormReal_t stormLsmrSolver<tArray, tOperator>::Iterate(tArray& xArr,
                                                        const tArray& bArr,
                                                        const tOperator& linOp,
                                                        const tOperator* preOp) {
  // ----------------------
  // Continue the bidiagonalization:
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲: 
  //   𝒔 ← 𝓟𝒗, 𝒕 ← 𝓐𝒔,
  // 𝗲𝗹𝘀𝗲: 𝒕 ← 𝓐𝒗, 𝗲𝗻𝗱 𝗶𝗳
  // 𝒕 ← 𝒕 - 𝛼𝒖,
  // 𝛽 ← ‖𝒕‖, 𝒖 ← 𝒕/𝛽,
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝒔 ← 𝓐*𝒖, 𝒕 ← 𝓟*𝒔, 
  // 𝗲𝗹𝘀𝗲: 𝒕 ← 𝓐*𝒖, 𝗲𝗻𝗱 𝗶𝗳
  // 𝒕 ← 𝒕 - 𝛽𝒗,
  // 𝛼 ← ‖𝒕‖, 𝒗 ← 𝒕/𝛼.
  // ----------------------
  if (preOp != nullptr) {
    preOp->MatVec(sArr, vArr);
    linOp.MatVec(tArr, sArr);
  } else {
    linOp.MatVec(tArr, vArr);
  }
  stormUtils::Sub(tArr, tArr, uArr, alpha);
  beta = stormUtils::Norm2(tArr); stormUtils::Scale(uArr, tArr, 1.0/beta);
  if (preOp != nullptr) {
    linOp.ConjMatVec(sArr, uArr);
    preOp->ConjMatVec(tArr, sArr);
  } else {
    linOp.ConjMatVec(tArr, uArr);
  }
  stormUtils::Sub(tArr, tArr, vArr, beta);
  alpha = stormUtils::Norm2(tArr); stormUtils::Scale(vArr, tArr, 1.0/alpha);
  
  // ----------------------
  // Construct and apply rotations:
  // 𝜌 ← (𝛼̅² + 𝛽²)¹ᐟ²,
  // 𝑐𝑠 ← 𝛼̅/𝜌, 𝑠𝑛 ← 𝛽/𝜌,
  // 𝜃 ← 𝑠𝑛⋅𝛼, 𝛼̅ ← 𝑐𝑠⋅𝛼,
  // 𝜃̅ ← 𝑠̅𝑛̅⋅𝜌, 𝜌̅ ← [(𝑐̅𝑠̅⋅𝜌)² + 𝜃²]¹ᐟ²,
  // 𝑐̅𝑠̅ ← 𝑐̅𝑠̅⋅𝜌/𝜌̅, 𝑠̅𝑛̅ ← 𝜃/𝜌̅,
  // 𝜓 ← 𝑐̅𝑠̅⋅𝜓̅, 𝜓̅ ← -𝑠̅𝑛̅⋅𝜓̅.
  // ----------------------
  rho = std::hypot(alphaBar, beta);
  cs = alphaBar/rho; sn = beta/rho;
  theta = sn*alpha; alphaBar = cs*alpha;
  thetaBar = snBar*rho; rhoBar = std::hypot(csBar*rho, theta);
  csBar = csBar*rho/rhoBar; snBar = theta/rhoBar;
  psi = csBar*psiBar; psiBar = -snBar*psiBar;

  // ----------------------
  // Update 𝒛-solution:
  // 𝒉 ← 𝒘 - (𝜃𝜌/𝜁)𝒉, 𝜁 ← 𝜌𝜌̅,
  // 𝒛 ← 𝒛 + (𝜓/𝜁)𝒉,
  // 𝒘 ← 𝒗 - (𝜃/𝜌)𝒘.
  // ----------------------
  stormUtils::Sub(hArr, wArr, hArr, thetaBar*rho/zeta); zeta = rho*rhoBar;
  stormUtils::Add(zArr, zArr, hArr, psi/zeta);
  stormUtils::Sub(wArr, vArr, wArr, theta/rho);

  return std::abs(psiBar);

} // stormLsmrSolver<...>::Iterate

template<class tArray, class tOperator>
void stormLsmrSolver<tArray, tOperator>::Finalize(tArray& xArr,
                                                  const tArray& bArr,
                                                  const tOperator& linOp,
                                                  const tOperator* preOp) {
  // ----------------------
  // Compute 𝒙-solution:
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝒕 ← 𝓟𝒛, 𝒙 ← 𝒙 + 𝒕.
  // 𝗲𝗹𝘀𝗲: 𝒙 ← 𝒙 + 𝒛. 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  if (preOp != nullptr) {
    preOp->MatVec(tArr, zArr);
    stormUtils::Add(xArr, xArr, tArr);
  } else {
    stormUtils::Add(xArr, xArr, zArr);
  }

} // stormLsmrSolver<...>::Finalize

/// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ///
/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ///

#endif // ifndef STORM_RULER_SOLVER_LSQR_HXX_
