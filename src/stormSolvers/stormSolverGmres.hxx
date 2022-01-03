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

#include <numeric>
#include <algorithm>

#include <stormSolvers/stormSolver.hxx>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Base class for @c GMRES and @c FGMRES.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<bool Flexible, class tArray>
class stormBaseGmresSolver : public stormRestartableSolver<tArray> {
private:
  std::vector<stormReal_t> beta, cs, sn;
  std::vector<std::vector<stormReal_t>> H;
  std::vector<tArray> QArr, ZArr;

  void PreInit(tArray& xArr,
               tArray const& bArr, 
               bool hasPreOp) override;

  stormReal_t ReInit(tArray& xArr,
                     tArray const& bArr,
                     stormOperator<tArray> const& linOp,
                     stormPreconditioner<tArray> const* preOp) override;

  stormReal_t ReIterate(stormSize_t k,
                        tArray& xArr,
                        tArray const& bArr,
                        stormOperator<tArray> const& linOp,
                        stormPreconditioner<tArray> const* preOp) override;

  void ReFinalize(stormSize_t k,
                  tArray& xArr,
                  tArray const& bArr,
                  stormOperator<tArray> const& linOp,
                  stormPreconditioner<tArray> const* preOp) override;

protected:

  stormBaseGmresSolver() = default;

}; // class stormBaseGmresSolver<...>

template<bool Flexible, class tArray>
void stormBaseGmresSolver<Flexible, tArray>::
                          PreInit(tArray& xArr,
                                  tArray const& bArr, 
                                  bool hasPreOp) {

  // ----------------------
  // Allocate the intermediate arrays:
  // ----------------------
  stormSize_t const m = this->NumIterationsBeforeRestart;
  beta.resize(m + 1), cs.resize(m), sn.resize(m);
  H.assign(m + 1, std::vector<stormReal_t>(m, 0.0));
  QArr.resize(m + 1);
  for (tArray& qArr : QArr) {
    stormUtils::AllocLike(xArr, qArr);
  }
  if (hasPreOp) {
    ZArr.resize(Flexible ? m : 1);
    for (tArray& zArr : ZArr) {
      stormUtils::AllocLike(xArr, zArr);
    }
  }

} // stormBaseGmresSolver<...>::Init

template<bool Flexible, class tArray>
stormReal_t stormBaseGmresSolver<Flexible, tArray>::
                                  ReInit(tArray& xArr,
                                         tArray const& bArr,
                                         stormOperator<tArray> const& linOp,
                                         stormPreconditioner<tArray> const* preOp) {

  // ----------------------
  // Initialize:
  // 𝓠₀ ← 𝓐𝒙,
  // 𝓠₀ ← 𝒃 - 𝓠₀,
  // 𝜑 ← ‖𝓠₀‖,
  // ----------------------
  linOp.MatVec(QArr[0], xArr);
  stormBlas::Sub(QArr[0], bArr, QArr[0]);
  stormReal_t const phi = stormBlas::Norm2(QArr[0]);

  // ----------------------
  // 𝒄𝒔 ← {𝟢}ᵀ, 𝒔𝒏 ← {𝟢}ᵀ,
  // 𝜷 ← {𝜑,𝟢,…,𝟢}ᵀ,
  // 𝓠₀ ← 𝓠₀/𝜑. 
  // ----------------------
  std::fill(cs.begin(), cs.end(), 0.0);
  std::fill(sn.begin(), sn.end(), 0.0);
  beta[0] = phi, std::fill(beta.begin() + 1, beta.end(), 0.0);
  stormBlas::Scale(QArr[0], QArr[0], 1.0/phi);

  return phi;

} // stormBaseGmresSolver<...>::ReInit

template<bool Flexible, class tArray>
stormReal_t stormBaseGmresSolver<Flexible, tArray>::
                              ReIterate(stormSize_t k,
                                        tArray& xArr,
                                        tArray const& bArr,
                                        stormOperator<tArray> const& linOp,
                                        stormPreconditioner<tArray> const* preOp) {

  // ----------------------
  // Continue the Arnoldi procedure:
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝑗 ← 𝘍𝘭𝘦𝘹𝘪𝘣𝘭𝘦 ? 𝑘 : 𝟢,
  //   𝓠ₖ₊₁, 𝓩ⱼ ← 𝓐𝓟𝓠ₖ, 𝓟𝓠ₖ,
  // 𝗲𝗹𝘀𝗲:
  //   𝓠ₖ₊₁ ← 𝓐𝓠ₖ,
  // 𝗲𝗻𝗱 𝗶𝗳
  // 𝗳𝗼𝗿 𝑖 = 𝟢, 𝑘 𝗱𝗼:
  //   𝓗ᵢₖ ← <𝓠ₖ₊₁⋅𝓠ᵢ>,
  //   𝓠ₖ₊₁ ← 𝓠ₖ₊₁ - 𝓗ᵢₖ𝓠ᵢ,
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // 𝓗ₖ₊₁,ₖ ← ‖𝓠ₖ₊₁‖, 
  // 𝓠ₖ₊₁ ← 𝓠ₖ₊₁/𝓗ₖ₊₁,ₖ.  
  // ----------------------
  if (preOp != nullptr) {
    stormSize_t const j = Flexible ? k : 0;
    stormUtils::MatVecRightPre(QArr[k + 1], ZArr[j], QArr[k], linOp, preOp);
  } else {
    linOp.MatVec(QArr[k + 1], QArr[k]);
  }
  for (stormSize_t i = 0; i <= k; ++i) {
    H[i][k] = stormBlas::Dot(QArr[k + 1], QArr[i]);
    stormBlas::Sub(QArr[k + 1], QArr[k + 1], QArr[i], H[i][k]);
  }
  H[k + 1][k] = stormBlas::Norm2(QArr[k + 1]); 
  stormBlas::Scale(QArr[k + 1], QArr[k + 1], 1.0/H[k + 1][k]);

  // ----------------------
  // Eliminate the last element in 𝓗
  // and and update the rotation matrix:
  // 𝗳𝗼𝗿 𝑖 = 𝟢, 𝑘 - 𝟣 𝗱𝗼:
  //   𝜒 ← 𝒄𝒔ᵢ⋅𝓗ᵢₖ + 𝒔𝒏ᵢ⋅𝓗ᵢ₊₁,ₖ,
  //   𝓗ᵢ₊₁,ₖ ← -𝒔𝒏ᵢ⋅𝓗ᵢₖ + 𝒄𝒔ᵢ⋅𝓗ᵢ₊₁,ₖ 
  //   𝓗ᵢₖ ← 𝜒,
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // 𝒄𝒔ₖ, 𝒔𝒏ₖ ← 𝘚𝘺𝘮𝘖𝘳𝘵𝘩𝘰(𝓗ₖₖ, 𝓗ₖ₊₁,ₖ),
  // 𝓗ₖₖ ← 𝒄𝒔ₖ⋅𝓗ₖₖ + 𝒔𝒏ₖ⋅𝓗ₖ₊₁,ₖ,
  // 𝓗ₖ₊₁,ₖ ← 𝟢.
  // ----------------------
  for (stormSize_t i = 0; i < k; ++i) {
    stormReal_t const chi = cs[i]*H[i][k] + sn[i]*H[i + 1][k];
    H[i + 1][k] = -sn[i]*H[i][k] + cs[i]*H[i+1][k];
    H[i][k] = chi;
  }
  std::tie(cs[k], sn[k], std::ignore) =
    stormBlas::SymOrtho(H[k][k], H[k + 1][k]);
  H[k][k] = cs[k]*H[k][k] + sn[k]*H[k + 1][k];
  H[k + 1][k] = 0.0;

  // ----------------------
  // Update the 𝜷-solution and residual norm:
  // 𝜷ₖ₊₁ ← -𝒔𝒏ₖ⋅𝜷ₖ, 𝜷ₖ ← 𝒄𝒔ₖ⋅𝜷ₖ,
  // 𝜑 ← |𝜷ₖ₊₁|,
  // ----------------------
  beta[k + 1] = -sn[k]*beta[k], beta[k] *= cs[k];
  stormReal_t const phi = std::abs(beta[k + 1]);

  return phi;

} // stormBaseGmresSolver<...>::ReIterate

template<bool Flexible, class tArray>
void stormBaseGmresSolver<Flexible, tArray>::
                      ReFinalize(stormSize_t k,
                                 tArray& xArr,
                                 tArray const& bArr,
                                 stormOperator<tArray> const& linOp,
                                 stormPreconditioner<tArray> const* preOp) {

  // ----------------------
  // Finalize the 𝜷-solution:
  // 𝜷ₖ ← 𝜷ₖ/𝓗ₖₖ,
  // 𝗳𝗼𝗿 𝑖 = 𝑘 - 𝟣, 𝟢, -𝟣 𝗱𝗼:
  //   𝜷ᵢ ← (𝜷ᵢ - <𝓗ᵢ,ᵢ₊₁:ₖ⋅𝜷ᵢ₊₁:ₖ>)/𝓗ᵢᵢ,
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // ----------------------
  beta[k] /= H[k][k];
  for (stormPtrDiff_t i = k - 1; i >= 0; --i) {
    beta[i] -= std::inner_product(
      beta.begin() + i + 1, beta.begin() + k + 1, H[i].begin() + i + 1, 0.0);
    beta[i] /= H[i][i];
  }

  // ----------------------
  // Compute 𝒙-solution:
  // 𝗶𝗳 𝓟 = 𝗻𝗼𝗻𝗲:
  //   𝗳𝗼𝗿 𝑖 = 𝟢, 𝑘 𝗱𝗼:
  //     𝒙 ← 𝒙 + 𝜷ᵢ⋅𝓠ᵢ,
  //   𝗲𝗻𝗱 𝗳𝗼𝗿
  // 𝗲𝗹𝘀𝗲 𝗶𝗳 𝘍𝘭𝘦𝘹𝘪𝘣𝘭𝘦:
  //   𝗳𝗼𝗿 𝑖 = 𝟢, 𝑘 𝗱𝗼:
  //     𝒙 ← 𝒙 + 𝜷ᵢ⋅𝓩ᵢ,
  //   𝗲𝗻𝗱 𝗳𝗼𝗿
  // 𝗲𝗹𝘀𝗲:
  //   𝓠₀ ← 𝜷₀𝓠₀,
  //   𝗳𝗼𝗿 𝑖 = 𝟣, 𝑘 𝗱𝗼:
  //     𝓠₀ ← 𝓠₀ + 𝜷ᵢ⋅𝓠ᵢ,
  //   𝗲𝗻𝗱 𝗳𝗼𝗿
  //   𝓩₀ ← 𝓟𝓠₀,
  //   𝒙 ← 𝒙 + 𝓩₀.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  if (preOp == nullptr) {
    for (stormSize_t i = 0; i <= k; ++i) {
      stormBlas::Add(xArr, xArr, QArr[i], beta[i]);
    }
  } else if constexpr (Flexible) {
    for (stormSize_t i = 0; i <= k; ++i) {
      stormBlas::Add(xArr, xArr, ZArr[i], beta[i]);
    }
  } else {
    stormBlas::Scale(QArr[0], QArr[0], beta[0]);
    for (stormSize_t i = 1; i <= k; ++i) {
      stormBlas::Add(QArr[0], QArr[0], QArr[i], beta[i]);
    }
    preOp->MatVec(ZArr[0], QArr[0]);
    stormBlas::Add(xArr, xArr, ZArr[0]);
  }

} // stormBaseGmresSolver<...>::ReFinalize

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Solve a linear operator equation with the
///   monstrous @c GMRES (Generalized Minimal Residual) method.
///
/// @c GMRES may be applied to the singular problems, and the square
/// least squares problems: ‖(𝓐[𝓟]𝒚 - 𝒃)‖₂ → 𝘮𝘪𝘯, 𝒙 = [𝓟]𝒚,
/// although convergeance to minimum norm solution is not guaranteed
/// (is this true?).
///
/// References:
/// @verbatim
/// [1] Saad, Yousef and Martin H. Schultz. 
///     “GMRES: A generalized minimal residual algorithm for solving 
///      nonsymmetric linear systems.” 
///     SIAM J. Sci. Stat. Comput., 7:856–869, 1986.
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class tArray>
class stormGmresSolver final : public stormBaseGmresSolver<false, tArray> {

}; // class stormGmresSolver<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Solve a linear operator equation, using \
///   the yet more monstrous @c FGMRES (Flexible Generalized \
///   Minimal Residual) method.
///
/// @c FGMRES allows usage of the variable (or flexible)
/// preconditioners with the price of doubleing of the memory 
/// requirements.
///
/// @c FGMRES may be applied to the singular problems, and the square
/// least squares problems: ‖(𝓐[𝓟]𝒚 - 𝒃)‖₂ → 𝘮𝘪𝘯, 𝒙 = [𝓟]𝒚,
/// although convergeance to minimum norm solution is not guaranteed
/// (is this true?).
///
/// References:
/// @verbatim
/// [1] Saad, Yousef. 
///     “A Flexible Inner-Outer Preconditioned GMRES Algorithm.” 
///     SIAM J. Sci. Comput. 14 (1993): 461-469.
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class tArray>
class stormFgmresSolver final : public stormBaseGmresSolver<true, tArray> {

}; // class stormFgmresSolver<...>

#endif // ifndef _STORM_SOLVER_GMRES_HXX_
