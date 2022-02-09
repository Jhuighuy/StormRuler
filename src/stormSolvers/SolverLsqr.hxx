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
#ifndef _STORM_SOLVER_LSQR_HXX_
#define _STORM_SOLVER_LSQR_HXX_

#include <cmath>

#include <stormSolvers/Solver.hxx>

/// ----------------------------------------------------------------- ///
/// @brief Base class for solvers, based on \
///   the Golub-Kahan-Lanczos bidiagonalization procedure.
///
/// @see @c LSQR, @c LSMR.
/// ----------------------------------------------------------------- ///
template<class tInArray, class tOutArray = tInArray>
class stormGolubKahanSolver : public stormIterativeSolver<tInArray, tOutArray> {
  static_assert(std::is_same_v<tInArray, tOutArray>,
            "Non-square case is not implemented yet.");

protected:

  /// @brief Initialize the bidiagonalization procedure.
  static void InitBidiagonalization(tInArray& sArr,
                                    tInArray& tArr,
                                    tInArray& uArr,
                                    tInArray& vArr,
                                    stormReal_t& alpha,
                                    stormReal_t& beta,
                                    tInArray const& bArr,
                                    stormOperator<tInArray, tOutArray> const& linOp,
                                    stormPreconditioner<tInArray> const* preOp);

  /// @brief Continue the bidiagonalization procedure.
  static void ContinueBidiagonalization(tInArray& sArr,
                                        tInArray& tArr,
                                        tInArray& uArr,
                                        tInArray& vArr,
                                        stormReal_t& alpha,
                                        stormReal_t& beta,
                                        stormOperator<tInArray, tOutArray> const& linOp,
                                        stormPreconditioner<tInArray> const* preOp);

}; // class stormGolubKahanSolver<...>

template<class tInArray, class tOutArray>
void stormGolubKahanSolver<tInArray, tOutArray>::
              InitBidiagonalization(tInArray& sArr,
                                    tInArray& tArr,
                                    tInArray& uArr,
                                    tInArray& vArr,
                                    stormReal_t& alpha,
                                    stormReal_t& beta,
                                    tInArray const& bArr,
                                    stormOperator<tInArray, tOutArray> const& linOp,
                                    stormPreconditioner<tInArray> const* preOp) {

  // ----------------------
  // Initialize the bidiagonalization procedure:
  // 𝛽 ← ‖𝒃‖, 𝒖 ← 𝒃/𝛽,
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝒔 ← 𝓐*𝒖, 𝒕 ← 𝓟*𝒔,
  // 𝗲𝗹𝘀𝗲: 𝒕 ← 𝓐*𝒖, 𝗲𝗻𝗱 𝗶𝗳
  // 𝛼 ← ‖𝒕‖, 𝒗 ← 𝒕/𝛼.
  // ----------------------
  beta = stormBlas::Norm2(bArr); stormBlas::Scale(uArr, bArr, 1.0/beta);
  if (preOp != nullptr) {
    linOp.ConjMatVec(sArr, uArr);
    preOp->ConjMatVec(tArr, sArr);
  } else {
    linOp.ConjMatVec(tArr, uArr);
  }
  alpha = stormBlas::Norm2(tArr); stormBlas::Scale(vArr, tArr, 1.0/alpha);

} // stormGolubKahanSolver<...>::InitBidiagonalization

template<class tInArray, class tOutArray>
void stormGolubKahanSolver<tInArray, tOutArray>::
          ContinueBidiagonalization(tInArray& sArr,
                                    tInArray& tArr,
                                    tInArray& uArr,
                                    tInArray& vArr,
                                    stormReal_t& alpha,
                                    stormReal_t& beta,
                                    stormOperator<tInArray, tOutArray> const& linOp,
                                    stormPreconditioner<tInArray> const* preOp) {

  // ----------------------
  // Continue the bidiagonalization procedure:
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
  stormBlas::Sub(tArr, tArr, uArr, alpha);
  beta = stormBlas::Norm2(tArr); stormBlas::Scale(uArr, tArr, 1.0/beta);
  if (preOp != nullptr) {
    linOp.ConjMatVec(sArr, uArr);
    preOp->ConjMatVec(tArr, sArr);
  } else {
    linOp.ConjMatVec(tArr, uArr);
  }
  stormBlas::Sub(tArr, tArr, vArr, beta);
  alpha = stormBlas::Norm2(tArr); stormBlas::Scale(vArr, tArr, 1.0/alpha);

} // stormGolubKahanSolver<...>::ContinueBidiagonalization

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Solve a right preconditioned linear least squares problem \
///   with the @c LSQR method.
///
/// @c LSQR is algebraically equivalent to applying @c CG
/// to the normal equations: (𝓐[𝓟])*𝓐[𝓟]𝒚 = (𝓐[𝓟])*𝒃, 𝒙 = [𝓟]𝒚,
/// (or, equivalently, [𝓟*]𝓐*𝓐[𝓟]𝒚 = [𝓟*]𝓐*𝒃, 𝒙 = [𝓟]𝒚),
/// but has better numerical properties.
///
/// @note The residual norm ‖𝓐[𝓟]𝒚 - 𝒃‖₂ decreases monotonically, \
///   while the normal equation's residual norm ‖(𝓐[𝓟])*(𝓐[𝓟]𝒚 - 𝒃)‖ \
///   is not guaranteed to decrease. Please make sure that the right \
///   stopping criterion is set.
///
/// @warning Using @c LSQR is not recommended in the \
///   self-adjoint case, please consider @c MINRES instead.
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
template<class tInArray, class tOutArray = tInArray>
class stormLsqrSolver final : public stormGolubKahanSolver<tInArray, tOutArray> {
private:
  stormReal_t alpha, beta, rho, rhoBar, theta, phi, phiBar, phiTilde, cs, sn;
  tInArray sArr, tArr, rArr, uArr, vArr, wArr, zArr;

  stormReal_t Init(tInArray& xArr,
                   tOutArray const& bArr,
                   stormOperator<tInArray, tOutArray> const& linOp,
                   stormPreconditioner<tInArray> const* preOp) override;

  stormReal_t Iterate(tInArray& xArr,
                      tOutArray const& bArr,
                      stormOperator<tInArray, tOutArray> const& linOp,
                      stormPreconditioner<tInArray> const* preOp) override;

  void Finalize(tInArray& xArr,
                tOutArray const& bArr,
                stormOperator<tInArray, tOutArray> const& linOp,
                stormPreconditioner<tInArray> const* preOp) override;

}; // class stormLsqrSolver<...>

template<class tInArray, class tOutArray>
stormReal_t stormLsqrSolver<tInArray, tOutArray>::
                                Init(tInArray& xArr,
                                     tOutArray const& bArr,
                                     stormOperator<tInArray, tOutArray> const& linOp,
                                     stormPreconditioner<tInArray> const* preOp) {

  // ----------------------
  // Allocate the intermediate arrays:
  // ----------------------
  /// @todo Some of these should be allocated like x, others like b.
  stormUtils::AllocLike(xArr, tArr, rArr, uArr, vArr, wArr, zArr);
  if (preOp != nullptr) {
    stormUtils::AllocLike(xArr, sArr);
  }

  // ----------------------
  // Utilize the initial guess
  // (solve 𝓐𝒛 = 𝒓 with zero initial guess, where 𝒓 = 𝒃 - 𝓐𝒙):
  // 𝒓 ← 𝓐𝒙,
  // 𝒓 ← 𝒃 - 𝒓,
  // 𝒛 ← {𝟢}ᵀ,
  // ----------------------
  linOp.MatVec(rArr, xArr);
  stormBlas::Sub(rArr, bArr, rArr);
  stormBlas::Fill(zArr, 0.0);

  // ----------------------
  // Initialize the bidiagonalization procedure:
  // 𝒖, 𝒗, 𝛼, 𝛽 ← 𝘉𝘪𝘋𝘪𝘢𝘨(𝒖, 𝒗, 𝛼, 𝛽, 𝒓, 𝓐[, 𝓟]).
  // ----------------------
  stormGolubKahanSolver<tInArray, tOutArray>::
    InitBidiagonalization(sArr, tArr, uArr, vArr, alpha, beta, rArr, linOp, preOp);

  // ----------------------
  // 𝜑̅ ← 𝛽, 𝜌̅ ← 𝛼.
  // 𝒘 ← 𝒗,
  // ----------------------
  phiBar = beta; rhoBar = alpha;
  stormBlas::Set(wArr, vArr);

  return phiBar;

} // stormLsqrSolver<...>::Init

template<class tInArray, class tOutArray>
stormReal_t stormLsqrSolver<tInArray, tOutArray>::
                            Iterate(tInArray& xArr,
                                    tOutArray const& bArr,
                                    stormOperator<tInArray, tOutArray> const& linOp,
                                    stormPreconditioner<tInArray> const* preOp) {

  // ----------------------
  // Continue the bidiagonalization procedure:
  // 𝒖, 𝒗, 𝛼, 𝛽 ← 𝘉𝘪𝘋𝘪𝘢𝘨(𝒖, 𝒗, 𝛼, 𝛽, 𝓐[, 𝓟]).
  // ----------------------
  stormGolubKahanSolver<tInArray, tOutArray>::
    ContinueBidiagonalization(sArr, tArr, uArr, vArr, alpha, beta, linOp, preOp);

  // ----------------------
  // Construct and apply rotation:
  // 𝜌 ← (𝜌̅² + 𝛽²)¹ᐟ²,
  // 𝑐𝑠 ← 𝜌̅/𝜌, 𝑠𝑛 ← 𝛽/𝜌,
  // 𝜃 ← 𝑠𝑛⋅𝛼, 𝜌̅ ← -𝑐𝑠⋅𝛼,
  // 𝜑 ← 𝑐𝑠⋅𝜑, 𝜑̅ ← 𝑠𝑛⋅𝜑̅.
  // ----------------------
  rho = std::hypot(rhoBar, beta);
  cs = rhoBar/rho, sn = beta/rho;
  theta = sn*alpha, rhoBar = -cs*alpha;
  phi = cs*phiBar, phiBar = sn*phiBar;

  // ----------------------
  // Update 𝒛-solution:
  // 𝒛 ← 𝒛 + (𝜑/𝜌)𝒘,
  // 𝒘 ← 𝒗 - (𝜃/𝜌)𝒘.
  // ----------------------
  stormBlas::Add(zArr, zArr, wArr, phi/rho);
  stormBlas::Sub(wArr, vArr, wArr, theta/rho);

  return phiBar;

} // stormLsqrSolver<...>::Iterate

template<class tInArray, class tOutArray>
void stormLsqrSolver<tInArray, tOutArray>::
                    Finalize(tInArray& xArr,
                             tOutArray const& bArr,
                             stormOperator<tInArray, tOutArray> const& linOp,
                             stormPreconditioner<tInArray> const* preOp) {

  // ----------------------
  // Compute 𝒙-solution:
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝒕 ← 𝓟𝒛, 𝒙 ← 𝒙 + 𝒕.
  // 𝗲𝗹𝘀𝗲: 𝒙 ← 𝒙 + 𝒛. 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  if (preOp != nullptr) {
    preOp->MatVec(tArr, zArr);
    stormBlas::Add(xArr, xArr, tArr);
  } else {
    stormBlas::Add(xArr, xArr, zArr);
  }

} // stormLsqrSolver<...>::Finalize

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Solve a right preconditioned linear least squares problem \
///   using the @c LSMR method.
///
/// @c LSMR is algebraically equivalent to applying @c MINRES
/// to the normal equations: (𝓐[𝓟])*𝓐[𝓟]𝒚 = (𝓐[𝓟])*𝒃, 𝒙 = [𝓟]𝒚,
/// (or, equivalently, [𝓟*]𝓐*𝓐[𝓟]𝒚 = [𝓟*]𝓐*𝒃, 𝒙 = [𝓟]𝒚),
/// but has better numerical properties.
///
/// @note The normal equation's residual norm ‖(𝓐[𝓟])*(𝓐[𝓟]𝒚 - 𝒃)‖ \
///   decreases monotonically, while the residual norm ‖𝓐[𝓟]𝒚 - 𝒃‖
///   is not guaranteed to decrease (but decreases on practice). \
///   Please make sure that the right stopping criterion is set.
///
/// @warning Using @c LSMR is not recommended in the \
///   self-adjoint case, please consider @c MINRES instead.
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
template<class tInArray, class tOutArray = tInArray>
class stormLsmrSolver final : public stormGolubKahanSolver<tInArray, tOutArray> {
private:
  stormReal_t alpha, alphaBar, beta, rho, rhoBar, cs, sn,
    theta, thetaBar, psi, psiBar, psiTilde, zeta, csBar, snBar;
  tInArray rArr, sArr, tArr, wArr, hArr, uArr, vArr, zArr;

  stormReal_t Init(tInArray& xArr,
                   tOutArray const& bArr,
                   stormOperator<tInArray, tOutArray> const& linOp,
                   stormPreconditioner<tInArray> const* preOp) override;

  stormReal_t Iterate(tInArray& xArr,
                      tOutArray const& bArr,
                      stormOperator<tInArray, tOutArray> const& linOp,
                      stormPreconditioner<tInArray> const* preOp) override;

  void Finalize(tInArray& xArr,
                tOutArray const& bArr,
                stormOperator<tInArray, tOutArray> const& linOp,
                stormPreconditioner<tInArray> const* preOp) override;

}; // class stormLsmrSolver<...>

template<class tInArray, class tOutArray>
stormReal_t stormLsmrSolver<tInArray, tOutArray>::
                                Init(tInArray& xArr,
                                     tOutArray const& bArr,
                                     stormOperator<tInArray, tOutArray> const& linOp,
                                     stormPreconditioner<tInArray> const* preOp) {

  // ----------------------
  // Allocate the intermediate arrays:
  // ----------------------
  /// @todo Some of these should be allocated like x, others like b.
  stormUtils::AllocLike(xArr, tArr, rArr, uArr, vArr, wArr, hArr, zArr);
  if (preOp != nullptr) {
    stormUtils::AllocLike(xArr, sArr);
  }

  // ----------------------
  // Utilize the initial guess
  // (solve 𝓐𝒛 = 𝒓 with zero initial guess, where 𝒓 = 𝒃 - 𝓐𝒙):
  // 𝒓 ← 𝓐𝒙,
  // 𝒓 ← 𝒃 - 𝒓,
  // 𝒛 ← {𝟢}ᵀ,
  // ----------------------
  linOp.MatVec(rArr, xArr);
  stormBlas::Sub(rArr, bArr, rArr);
  stormBlas::Fill(zArr, 0.0);

  // ----------------------
  // Initialize the bidiagonalization procedure:
  // 𝒖, 𝒗, 𝛼, 𝛽 ← 𝘉𝘪𝘋𝘪𝘢𝘨(𝒖, 𝒗, 𝛼, 𝛽, 𝒓, 𝓐[, 𝓟]).
  // ----------------------
  stormGolubKahanSolver<tInArray, tOutArray>::
    InitBidiagonalization(sArr, tArr, uArr, vArr, alpha, beta, rArr, linOp, preOp);

  // ----------------------
  // 𝛼̅ ← 𝛼, 𝜓̅ ← 𝛼𝛽,
  // 𝜁 ← 𝟣, 𝑐̅𝑠̅ ← 𝟣, 𝑠̅𝑛̅ ← 𝟢,
  // 𝒘 ← 𝒗, 𝒉 ← {𝟢}ᵀ.
  // ----------------------
  alphaBar = alpha; psiBar = alpha*beta;
  zeta = 1.0; csBar = 1.0; snBar = 0.0;
  stormBlas::Set(wArr, vArr); stormBlas::Fill(hArr, 0.0);

  return std::abs(psiBar);

} // stormLsmrSolver<...>::Init

template<class tInArray, class tOutArray>
stormReal_t stormLsmrSolver<tInArray, tOutArray>::
                            Iterate(tInArray& xArr,
                                    tOutArray const& bArr,
                                    stormOperator<tInArray, tOutArray> const& linOp,
                                    stormPreconditioner<tInArray> const* preOp) {

  // ----------------------
  // Continue the bidiagonalization procedure:
  // 𝒖, 𝒗, 𝛼, 𝛽 ← 𝘉𝘪𝘋𝘪𝘢𝘨(𝒖, 𝒗, 𝛼, 𝛽, 𝓐[, 𝓟]).
  // ----------------------
  stormGolubKahanSolver<tInArray, tOutArray>::
    ContinueBidiagonalization(sArr, tArr, uArr, vArr, alpha, beta, linOp, preOp);

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
  cs = alphaBar/rho, sn = beta/rho;
  theta = sn*alpha, alphaBar = cs*alpha;
  thetaBar = snBar*rho, rhoBar = std::hypot(csBar*rho, theta);
  csBar = csBar*rho/rhoBar, snBar = theta/rhoBar;
  psi = csBar*psiBar, psiBar = -snBar*psiBar;

  // ----------------------
  // Update 𝒛-solution:
  // 𝒉 ← 𝒘 - (𝜃𝜌/𝜁)𝒉, 𝜁 ← 𝜌𝜌̅,
  // 𝒛 ← 𝒛 + (𝜓/𝜁)𝒉,
  // 𝒘 ← 𝒗 - (𝜃/𝜌)𝒘.
  // ----------------------
  stormBlas::Sub(hArr, wArr, hArr, thetaBar*rho/zeta); zeta = rho*rhoBar;
  stormBlas::Add(zArr, zArr, hArr, psi/zeta);
  stormBlas::Sub(wArr, vArr, wArr, theta/rho);

  return std::abs(psiBar);

} // stormLsmrSolver<...>::Iterate

template<class tInArray, class tOutArray>
void stormLsmrSolver<tInArray, tOutArray>::
                    Finalize(tInArray& xArr,
                             tOutArray const& bArr,
                             stormOperator<tInArray, tOutArray> const& linOp,
                             stormPreconditioner<tInArray> const* preOp) {

  // ----------------------
  // Compute 𝒙-solution:
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝒕 ← 𝓟𝒛, 𝒙 ← 𝒙 + 𝒕.
  // 𝗲𝗹𝘀𝗲: 𝒙 ← 𝒙 + 𝒛. 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  if (preOp != nullptr) {
    preOp->MatVec(tArr, zArr);
    stormBlas::Add(xArr, xArr, tArr);
  } else {
    stormBlas::Add(xArr, xArr, zArr);
  }

} // stormLsmrSolver<...>::Finalize

#endif // ifndef _STORM_SOLVER_LSQR_HXX_
