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
#ifndef _STORM_SOLVER_GMRES_HXX_
#define _STORM_SOLVER_GMRES_HXX_

#include <algorithm>

#include <stormSolvers/stormSolver.hxx>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Solve a linear operator equation: [𝓟]𝓐𝒙 = [𝓟]𝒃, using
///   the monstrous @c GMRES (Generalized Minimal Residual) method.
///
/// Preconditioned residual norm, ‖[𝓟]𝒓‖, where 𝒓 = 𝒃 - 𝓐𝒙, is reported.
///
/// The classical GMRES(𝑚) implementation with restarts
/// after 𝑚 iterations is used.
/// 
/// @c GMRES may be applied to the singular problems, and the square
/// least squares problems: ‖(𝓐[𝓟]𝒚 - 𝒃)‖₂ → 𝘮𝘪𝘯, 𝒙 = [𝓟]𝒚,
/// although convergeance to minimum norm solution is not guaranteed
/// (is this true?).
///
/// References:
/// @verbatim
/// [1] Saad and M.H. Schultz,
///     "GMRES: A generalized minimal residual algorithm for solving
///      nonsymmetric linear systems",
///     SIAM J. Sci. Stat. Comput., 7:856–869, 1986.
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class tArray>
class stormGmresSolver final : public stormRestartableSolver<tArray> {
private:
  tArray rArr;
  std::vector<tArray> QArr;
  std::vector<stormReal_t> beta, cs, sn;
  std::vector<std::vector<stormReal_t>> H;

  void PreInit(tArray& xArr,
               const tArray& bArr, 
               bool hasPreOp) override;

  stormReal_t ReInit(tArray& xArr,
                     const tArray& bArr,
                     const stormOperator<tArray>& linOp,
                     const stormPreconditioner<tArray>* preOp) override;

  stormReal_t ReIterate(stormSize_t k,
                        tArray& xArr,
                        const tArray& bArr,
                        const stormOperator<tArray>& linOp,
                        const stormPreconditioner<tArray>* preOp) override;

  void ReFinalize(stormSize_t k,
                  tArray& xArr,
                  const tArray& bArr,
                  const stormOperator<tArray>& linOp,
                  const stormPreconditioner<tArray>* preOp) override;

}; // class stormGmresSolver<...>

template<class tArray>
void stormGmresSolver<tArray>::PreInit(tArray& xArr,
                                       const tArray& bArr, 
                                       bool hasPreOp) {

  // ----------------------
  // Allocate the intermediate arrays:
  // ----------------------
  const stormSize_t m = this->NumIterationsBeforeRestart;
  beta.resize(m + 1), cs.resize(m), sn.resize(m);
  H.assign(m + 1, std::vector<stormReal_t>(m, 0.0));
  stormUtils::AllocLike(xArr, rArr);
  QArr.resize(m + 1);
  for (stormSize_t i = 0; i <= m; ++i) {
    stormUtils::AllocLike(xArr, QArr[i]);
  }
  if (hasPreOp) {
    //stormUtils::AllocLike(xArr, wArr, yArr, zArr);
  }

} // stormGmresSolver<...>::Init

template<class tArray>
stormReal_t stormGmresSolver<tArray>::ReInit(tArray& xArr,
                                             const tArray& bArr,
                                             const stormOperator<tArray>& linOp,
                                             const stormPreconditioner<tArray>* preOp) {

  // ----------------------
  // Initialize:
  // 𝒓 ← 𝓐𝒙,
  // 𝒓 ← 𝒃 - 𝒓,
  // 𝜑 ← ‖𝒓‖,
  // ----------------------
  linOp.MatVec(rArr, xArr);
  stormUtils::Sub(rArr, bArr, rArr);
  const stormReal_t phi = stormUtils::Norm2(rArr);

  // ----------------------
  // 𝒄𝒔 ← {𝟢}ᵀ, 𝒔𝒏 ← {𝟢}ᵀ,
  // 𝜷 ← {𝜑,𝟢,…,𝟢}ᵀ,
  // 𝓠₁ ← 𝒓/𝜑. 
  // ----------------------
  std::fill(cs.begin(), cs.end(), 0.0);
  std::fill(sn.begin(), sn.end(), 0.0);
  beta[0] = phi; std::fill(beta.begin() + 1, beta.end(), 0.0);
  stormUtils::Scale(QArr[0], rArr, 1.0/phi);

  return phi;

} // stormGmresSolver<...>::ReInit

template<class tArray>
stormReal_t stormGmresSolver<tArray>::ReIterate(stormSize_t k,
                                                tArray& xArr,
                                                const tArray& bArr,
                                                const stormOperator<tArray>& linOp,
                                                const stormPreconditioner<tArray>* preOp) {
  std::cout << "ReI " << k << std::endl; 

  // ----------------------
  // Arnoldi iteration:
  // 𝓠ₖ₊₁ ← 𝓐𝓠ₖ,
  // 𝗳𝗼𝗿 𝑖 = 𝟣, 𝑘 𝗱𝗼:
  //   𝓗ᵢₖ ← <𝓠ₖ₊₁⋅𝓠ᵢ>,
  //   𝓠ₖ₊₁ ← 𝓠ₖ₊₁ - 𝓗ᵢₖ𝓠ᵢ,
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // 𝓗ₖ₊₁,ₖ ← ‖𝓠ₖ₊₁‖, 
  // 𝓠ₖ₊₁ ← 𝓠ₖ₊₁/𝓗ₖ₊₁,ₖ.  
  // ----------------------
  linOp.MatVec(QArr[k + 1], QArr[k]);
  for (stormSize_t i = 0; i <= k; ++i) {
    H[i][k] = stormUtils::Dot(QArr[k + 1], QArr[i]);
    stormUtils::Sub(QArr[k + 1], QArr[k + 1], QArr[i], H[i][k]);
  }
  H[k + 1][k] = stormUtils::Norm2(QArr[k + 1]); 
  stormUtils::Scale(QArr[k + 1], QArr[k + 1], 1.0/H[k + 1][k]);

  // ----------------------
  // Eliminate the last element in 𝓗
  // and and update the rotation matrix:
  // 𝗳𝗼𝗿 𝑖 = 𝟣, 𝑘 - 𝟣 𝗱𝗼:
  //   𝜒 ← 𝒄𝒔ᵢ⋅𝓗ᵢₖ + 𝒔𝒏ᵢ⋅𝓗ᵢ₊₁,ₖ,
  //   𝓗ᵢ₊₁,ₖ ← -𝒔𝒏ᵢ⋅𝓗ᵢₖ + 𝒄𝒔ᵢ⋅𝓗ᵢ₊₁,ₖ 
  //   𝓗ᵢₖ ← 𝜒,
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // 𝒄𝒔ₖ, 𝒔𝒏ₖ ← 𝘚𝘺𝘮𝘖𝘳𝘵𝘩𝘰(𝓗ₖₖ, 𝓗ₖ₊₁,ₖ),
  // 𝓗ₖₖ ← 𝒄𝒔ₖ⋅𝓗ₖₖ + 𝒔𝒏ₖ⋅𝓗ₖ₊₁,ₖ,
  // 𝓗ₖ₊₁,ₖ ← 𝟢.
  // ----------------------
  for (stormSize_t i = 0; i < k; ++i) {
    const stormReal_t chi = cs[i]*H[i][k] + sn[i]*H[i + 1][k];
    H[i + 1][k] = -sn[i]*H[i][k] + cs[i]*H[i+1][k];
    H[i][k] = chi;
  }
  std::tie(cs[k], sn[k], std::ignore) =
    stormUtils::SymOrtho(H[k][k], H[k + 1][k]);
  H[k][k] = cs[k]*H[k][k] + sn[k]*H[k + 1][k];
  H[k + 1][k] = 0.0;

  // ----------------------
  // Update the residual norm:
  // 𝜷ₖ₊₁ ← -𝒔𝒏ₖ⋅𝜷ₖ, 𝜷ₖ ← 𝒄𝒔ₖ⋅𝜷ₖ,
  // 𝜑 ← |𝜷ₖ₊₁|,
  // ----------------------
  beta[k + 1] = -sn[k]*beta[k], beta[k] *= cs[k];
  const stormReal_t phi = std::abs(beta[k + 1]);

  return phi;

} // stormGmresSolver<...>::ReIterate

template<class tArray>
void stormGmresSolver<tArray>::ReFinalize(stormSize_t k,
                                          tArray& xArr,
                                          const tArray& bArr,
                                          const stormOperator<tArray>& linOp,
                                          const stormPreconditioner<tArray>* preOp) {
  std::cout << "ReF " << k << std::endl; 

  // ----------------------
  // Compute 𝒙-solution:
  // 𝜷₁:ₖ ← (𝓗₁:ₖ,₁:ₖ)⁻¹𝜷₁:ₖ, 
  // 𝗳𝗼𝗿 𝑖 = 1, 𝑘 𝗱𝗼:
  //   𝒙 ← 𝒙 + 𝜷ᵢ𝓠ᵢ.
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // // Since 𝓗₁:ₖ is upper triangular, 
  // // operations can be combined:
  // 𝗳𝗼𝗿 𝑖 = 𝑘, 𝟣, -𝟣 𝗱𝗼:
  //   𝜷ᵢ ← (𝜷ᵢ - <𝓗ᵢ,ᵢ₊₁:ₖ⋅𝜷ᵢ₊₁:ₖ>)/𝓗ᵢᵢ,
  //   𝒙 ← 𝒙 + 𝜷ᵢ𝓠ᵢ.
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // ----------------------
  for (stormPtrDiff_t i = k; i >= 0; --i) {
    for (stormSize_t j = i + 1; j <= k; ++j) {
      beta[i] -= H[i][j]*beta[j];
    }
    beta[i] /= H[i][i];
    stormUtils::Add(xArr, xArr, QArr[i], beta[i]);
  }

} // stormGmresSolver<...>::ReFinalize

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Solve a linear operator equation: [𝓟]𝓐𝒙 = [𝓟]𝒃, using \
///   the yet more monstrous @c FGMRES (Flexible Generalized \
///   Minimal Residual) method.
///
/// Preconditioned residual norm, ‖[𝓟]𝒓‖, where 𝒓 = 𝒃 - 𝓐𝒙, is reported.
///
/// @c FGMRES is intended to be used with flexible preconditioners.
///
/// The classical FGMRES(𝑚) implementation with restarts
/// after 𝑚 iterations is used.
///
/// @c FGMRES may be applied to the singular problems, and the square
/// least squares problems: ‖(𝓐[𝓟]𝒚 - 𝒃)‖₂ → 𝘮𝘪𝘯, 𝒙 = [𝓟]𝒚,
/// although convergeance to minimum norm solution is not guaranteed
/// (is this true?).
///
/// References:
/// @verbatim
/// [1] Saad and M.H. Schultz,
///     "GMRES: A generalized minimal residual algorithm for solving
///      nonsymmetric linear systems",
///     SIAM J. Sci. Stat. Comput., 7:856–869, 1986.
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class tArray>
class stormFgmresSolver final : public stormRestartableSolver<tArray> {
private:

  void PreInit(tArray& xArr,
               const tArray& bArr, 
               bool hasPreOp) override;

  stormReal_t ReInit(tArray& xArr,
                     const tArray& bArr,
                     const stormOperator<tArray>& linOp,
                     const stormPreconditioner<tArray>* preOp) override;

  stormReal_t ReIterate(stormSize_t k,
                        tArray& xArr,
                        const tArray& bArr,
                        const stormOperator<tArray>& linOp,
                        const stormPreconditioner<tArray>* preOp) override;

  void ReFinalize(stormSize_t k,
                  tArray& xArr,
                  const tArray& bArr,
                  const stormOperator<tArray>& linOp,
                  const stormPreconditioner<tArray>* preOp) override;

}; // class stormFgmresSolver<...>

template<class tArray>
void stormFgmresSolver<tArray>::PreInit(tArray& xArr,
                                       const tArray& bArr, 
                                       bool hasPreOp) {
  // ----------------------
  // Allocate the intermediate arrays:
  // ----------------------
  //stormUtils::AllocLike(xArr, pArr, rArr, rTildeArr, sArr, tArr, vArr);
  if (hasPreOp) {
    //stormUtils::AllocLike(xArr, wArr, yArr, zArr);
  }

  _STORM_NOT_IMPLEMENTED_();

} // stormFgmresSolver<...>::Init

template<class tArray>
stormReal_t stormFgmresSolver<tArray>::ReInit(tArray& xArr,
                                              const tArray& bArr,
                                              const stormOperator<tArray>& linOp,
                                              const stormPreconditioner<tArray>* preOp) {

  _STORM_NOT_IMPLEMENTED_();

} // stormFgmresSolver<...>::ReInit

template<class tArray>
stormReal_t stormFgmresSolver<tArray>::ReIterate(stormSize_t k,
                                                 tArray& xArr,
                                                 const tArray& bArr,
                                                 const stormOperator<tArray>& linOp,
                                                 const stormPreconditioner<tArray>* preOp) {

  _STORM_NOT_IMPLEMENTED_();

} // stormFgmresSolver<...>::ReIterate

template<class tArray>
void stormFgmresSolver<tArray>::ReFinalize(stormSize_t k,
                                           tArray& xArr,
                                           const tArray& bArr,
                                           const stormOperator<tArray>& linOp,
                                           const stormPreconditioner<tArray>* preOp) {

  _STORM_NOT_IMPLEMENTED_();

} // stormFgmresSolver<...>::ReFinalize

#endif // ifndef _STORM_SOLVER_GMRES_HXX_
