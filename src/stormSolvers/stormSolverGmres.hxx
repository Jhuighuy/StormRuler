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

#include <array>
#include <vector>
#include <type_traits>
#include <numeric>
#include <algorithm>

#include <stormSolvers/stormSolver.hxx>

/// ----------------------------------------------------------------- ///
/// @brief Base class for @c GMRES and @c FGMRES.
/// ----------------------------------------------------------------- ///
template<bool Flexible, class tArray>
class stormBaseGmresSolver : public stormInnerOuterIterativeSolver<tArray> {
private:
  std::vector<stormReal_t> beta, cs, sn;
  std::vector<std::vector<stormReal_t>> h;
  std::vector<tArray> qArr;
  std::conditional_t<Flexible, std::vector<tArray>, std::array<tArray, 1>> zArr;

  void OuterInit(tArray& xArr,
                 tArray const& bArr,
                 stormOperator<tArray> const& linOp,
                 stormPreconditioner<tArray> const* preOp) override;

  stormReal_t InnerInit(tArray& xArr,
                        tArray const& bArr,
                        stormOperator<tArray> const& linOp,
                        stormPreconditioner<tArray> const* preOp) override;

  stormReal_t InnerIterate(tArray& xArr,
                           tArray const& bArr,
                           stormOperator<tArray> const& linOp,
                           stormPreconditioner<tArray> const* preOp) override;

  void InnerFinalize(tArray& xArr,
                     tArray const& bArr,
                     stormOperator<tArray> const& linOp,
                     stormPreconditioner<tArray> const* preOp) override;

protected:

  stormBaseGmresSolver() = default;

}; // class stormBaseGmresSolver<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Solve a linear operator equation with the
///   monstrous @c GMRES (Generalized Minimal Residual) method.
///
/// @c GMRES is typically more robust than the @c BiCG type solvers, \
///   but it may be slower than the @c BiCG solvers for the \
///   well-conditioned moderate sized problems.
///
/// In the self-adjoint operator unpreconditioned case, \
///   @c GMRES, is algebraically equivalent to @c MINRES method, \
///   however, the need for restarts may lead to the much slower \
///   @c GMRES convergence rate. 
///
/// @c GMRES may be applied to the singular problems, and the square \
///   least squares problems, although, similarly to @c MINRES, \
///   convergeance to minimum norm solution is not guaranteed.
///
/// References:
/// @verbatim
/// [1] Saad, Yousef and Martin h. Schultz. 
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
/// @c FGMRES is typically more robust than the @c BiCG type solvers, \
///   but it may be slower than the @c BiCG solvers for the \
///   well-conditioned moderate sized problems.
///
/// @c FGMRES does the same amount of operations per iteration \
///   as @c GMRES, but also allows usage of the variable (or flexible) \
///   preconditioners with the price of doubleing of the memory \ 
///   requirements. For the static preconditioners, @c FGMRES requires \
///   one preconditioner-vector product less than @c GMRES. \
///   @c FGMRES supports only the right preconditioning.
///
/// @c FGMRES may be applied to the singular problems, and the square \
///   least squares problems, although, similarly to @c MINRES, \
///   convergeance to minimum norm solution is not guaranteed.
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

template<bool Flexible, class tArray>
void stormBaseGmresSolver<Flexible, tArray>::
                        OuterInit(tArray& xArr,
                                  tArray const& bArr,
                                  stormOperator<tArray> const& linOp,
                                  stormPreconditioner<tArray> const* preOp) {

  stormSize_t const m = this->NumInnerIterations;
  beta.resize(m + 1), cs.resize(m), sn.resize(m);
  h.assign(m + 1, std::vector<stormReal_t>(m, 0.0));
  qArr.resize(m + 1);
  for (tArray& qArr : qArr) {
    stormUtils::AllocLike(xArr, qArr);
  }
  if (preOp != nullptr) {
    if constexpr (Flexible) {
      zArr.resize(m);
    }
    for (tArray& zArr : zArr) {
      stormUtils::AllocLike(xArr, zArr);
    }
  }

} // stormBaseGmresSolver<...>::OuterInit

template<bool Flexible, class tArray>
stormReal_t stormBaseGmresSolver<Flexible, tArray>::
                               InnerInit(tArray& xArr,
                                         tArray const& bArr,
                                         stormOperator<tArray> const& linOp,
                                         stormPreconditioner<tArray> const* preOp) {

  bool const leftPre = (preOp != nullptr) && 
    (!Flexible) && (this->PreSide == stormPreconditionerSide::Left);

  // ----------------------
  // Initialize:
  // 𝒒₀ ← 𝓐𝒙,
  // 𝒒₀ ← 𝒃 - 𝒒₀,
  // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
  //   𝒛₀ ← 𝒒₀,
  //   𝒒₀ ← 𝓟𝒛₀.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  linOp.MatVec(qArr[0], xArr);
  stormBlas::Sub(qArr[0], bArr, qArr[0]);
  if (leftPre) {
    std::swap(zArr[0], qArr[0]);
    preOp->MatVec(qArr[0], zArr[0]);
  }

  // ----------------------
  // 𝑐𝑠 ← {𝟢}ᵀ, 𝑠𝑛 ← {𝟢}ᵀ,
  // 𝜑 ← ‖𝒒₀‖,
  // 𝛽 ← {𝜑,𝟢,…,𝟢}ᵀ,
  // 𝒒₀ ← 𝒒₀/𝜑. 
  // ----------------------
  std::fill(cs.begin(), cs.end(), 0.0);
  std::fill(sn.begin(), sn.end(), 0.0);
  stormReal_t const phi = stormBlas::Norm2(qArr[0]);
  beta[0] = phi, std::fill(beta.begin() + 1, beta.end(), 0.0);
  stormBlas::Scale(qArr[0], qArr[0], 1.0/phi);

  return phi;

} // stormBaseGmresSolver<...>::InnerInit

template<bool Flexible, class tArray>
stormReal_t stormBaseGmresSolver<Flexible, tArray>::
                            InnerIterate(tArray& xArr,
                                         tArray const& bArr,
                                         stormOperator<tArray> const& linOp,
                                         stormPreconditioner<tArray> const* preOp) {

  stormSize_t const k = this->InnerIteration;

  bool const leftPre = (preOp != nullptr) && 
    (!Flexible) && (this->PreSide == stormPreconditionerSide::Left);
  bool const rightPre = (preOp != nullptr) && 
    (Flexible || (this->PreSide == stormPreconditionerSide::Right));

  // ----------------------
  // Continue the Arnoldi procedure:
  // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
  //   𝒒ₖ₊₁ ← 𝓟(𝒛₀ ← 𝓐𝑞ₖ),
  // 𝗲𝗹𝘀𝗲 𝗶𝗳 𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦:
  //   𝑗 ← 𝘍𝘭𝘦𝘹𝘪𝘣𝘭𝘦 ? 𝑘 : 𝟢,
  //   𝒒ₖ₊₁ ← 𝓐(𝒛ⱼ ← 𝓟𝑞ₖ),
  // 𝗲𝗹𝘀𝗲:
  //   𝒒ₖ₊₁ ← 𝓐𝒒ₖ,
  // 𝗲𝗻𝗱 𝗶𝗳
  // 𝗳𝗼𝗿 𝑖 = 𝟢, 𝑘 𝗱𝗼:
  //   𝒉ᵢₖ ← <𝒒ₖ₊₁⋅𝒒ᵢ>,
  //   𝒒ₖ₊₁ ← 𝒒ₖ₊₁ - 𝒉ᵢₖ⋅𝒒ᵢ,
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // 𝒉ₖ₊₁,ₖ ← ‖𝒒ₖ₊₁‖, 
  // 𝒒ₖ₊₁ ← 𝒒ₖ₊₁/𝒉ₖ₊₁,ₖ.  
  // ----------------------
  if (leftPre) {
    stormBlas::MatVec(qArr[k + 1], *preOp, zArr[0], linOp, qArr[k]);
  } else if (rightPre) {
    stormSize_t const j = Flexible ? k : 0;
    stormBlas::MatVec(qArr[k + 1], linOp, zArr[j], *preOp, qArr[k]);
  } else {
    linOp.MatVec(qArr[k + 1], qArr[k]);
  }
  for (stormSize_t i = 0; i <= k; ++i) {
    h[i][k] = stormBlas::Dot(qArr[k + 1], qArr[i]);
    stormBlas::Sub(qArr[k + 1], qArr[k + 1], qArr[i], h[i][k]);
  }
  h[k + 1][k] = stormBlas::Norm2(qArr[k + 1]); 
  stormBlas::Scale(qArr[k + 1], qArr[k + 1], 1.0/h[k + 1][k]);

  // ----------------------
  // Eliminate the last element in {𝒉ᵢⱼ}
  // and and update the rotation matrix:
  // 𝗳𝗼𝗿 𝑖 = 𝟢, 𝑘 - 𝟣 𝗱𝗼:
  //   𝜒 ← 𝑐𝑠ᵢ⋅𝒉ᵢₖ + 𝑠𝑛ᵢ⋅𝒉ᵢ₊₁,ₖ,
  //   𝒉ᵢ₊₁,ₖ ← -𝑠𝑛ᵢ⋅𝒉ᵢₖ + 𝑐𝑠ᵢ⋅𝒉ᵢ₊₁,ₖ,
  //   𝒉ᵢₖ ← 𝜒,
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // 𝑐𝑠ₖ, 𝑠𝑛ₖ ← 𝘚𝘺𝘮𝘖𝘳𝘵𝘩𝘰(𝒉ₖₖ, 𝒉ₖ₊₁,ₖ),
  // 𝒉ₖₖ ← 𝑐𝑠ₖ⋅𝒉ₖₖ + 𝑠𝑛ₖ⋅𝒉ₖ₊₁,ₖ,
  // 𝒉ₖ₊₁,ₖ ← 𝟢.
  // ----------------------
  for (stormSize_t i = 0; i < k; ++i) {
    stormReal_t const chi = cs[i]*h[i][k] + sn[i]*h[i + 1][k];
    h[i + 1][k] = -sn[i]*h[i][k] + cs[i]*h[i+1][k];
    h[i][k] = chi;
  }
  std::tie(cs[k], sn[k], std::ignore) =
    stormBlas::SymOrtho(h[k][k], h[k + 1][k]);
  h[k][k] = cs[k]*h[k][k] + sn[k]*h[k + 1][k];
  h[k + 1][k] = 0.0;

  // ----------------------
  // Update the 𝛽-solution and residual norm:
  // 𝛽ₖ₊₁ ← -𝑠𝑛ₖ⋅𝛽ₖ, 𝛽ₖ ← 𝑐𝑠ₖ⋅𝛽ₖ,
  // 𝜑 ← |𝛽ₖ₊₁|.
  // ----------------------
  beta[k + 1] = -sn[k]*beta[k], beta[k] *= cs[k];
  stormReal_t const phi = std::abs(beta[k + 1]);

  return phi;

} // stormBaseGmresSolver<...>::InnerIterate

template<bool Flexible, class tArray>
void stormBaseGmresSolver<Flexible, tArray>::
                    InnerFinalize(tArray& xArr,
                                  tArray const& bArr,
                                  stormOperator<tArray> const& linOp,
                                  stormPreconditioner<tArray> const* preOp) {

  stormSize_t const k = this->InnerIteration;

  bool const rightPre = (preOp != nullptr) && 
    (Flexible || (this->PreSide == stormPreconditionerSide::Right));

  // ----------------------
  // Finalize the 𝛽-solution:
  // 𝛽ₖ ← 𝛽ₖ/𝒉ₖₖ,
  // 𝗳𝗼𝗿 𝑖 = 𝑘 - 𝟣, 𝟢, -𝟣 𝗱𝗼:
  //   𝛽ᵢ ← (𝛽ᵢ - <𝒉ᵢ,ᵢ₊₁:ₖ⋅𝛽ᵢ₊₁:ₖ>)/𝒉ᵢᵢ.
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // ----------------------
  beta[k] /= h[k][k];
  for (stormPtrDiff_t i = k - 1; i >= 0; --i) {
    beta[i] -= std::inner_product(
      beta.begin() + i + 1, beta.begin() + k + 1, h[i].begin() + i + 1, 0.0);
    beta[i] /= h[i][i];
  }

  // ----------------------
  // Compute 𝒙-solution:
  // 𝗶𝗳 𝗻𝗼𝘁 𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦:
  //   𝗳𝗼𝗿 𝑖 = 𝟢, 𝑘 𝗱𝗼:
  //     𝒙 ← 𝒙 + 𝛽ᵢ⋅𝒒ᵢ.
  //   𝗲𝗻𝗱 𝗳𝗼𝗿
  // 𝗲𝗹𝘀𝗲 𝗶𝗳 𝘍𝘭𝘦𝘹𝘪𝘣𝘭𝘦:
  //   𝗳𝗼𝗿 𝑖 = 𝟢, 𝑘 𝗱𝗼:
  //     𝒙 ← 𝒙 + 𝛽ᵢ⋅𝒛ᵢ.
  //   𝗲𝗻𝗱 𝗳𝗼𝗿
  // 𝗲𝗹𝘀𝗲:
  //   𝒒₀ ← 𝛽₀⋅𝒒₀,
  //   𝗳𝗼𝗿 𝑖 = 𝟣, 𝑘 𝗱𝗼:
  //     𝒒₀ ← 𝒒₀ + 𝛽ᵢ⋅𝒒ᵢ,
  //   𝗲𝗻𝗱 𝗳𝗼𝗿
  //   𝒛₀ ← 𝓟𝒒₀,
  //   𝒙 ← 𝒙 + 𝒛₀.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  if (!rightPre) {
    for (stormSize_t i = 0; i <= k; ++i) {
      stormBlas::Add(xArr, xArr, qArr[i], beta[i]);
    }
  } else if constexpr (Flexible) {
    for (stormSize_t i = 0; i <= k; ++i) {
      stormBlas::Add(xArr, xArr, zArr[i], beta[i]);
    }
  } else {
    stormBlas::Scale(qArr[0], qArr[0], beta[0]);
    for (stormSize_t i = 1; i <= k; ++i) {
      stormBlas::Add(qArr[0], qArr[0], qArr[i], beta[i]);
    }
    preOp->MatVec(zArr[0], qArr[0]);
    stormBlas::Add(xArr, xArr, zArr[0]);
  }

} // stormBaseGmresSolver<...>::InnerFinalize

#endif // ifndef _STORM_SOLVER_GMRES_HXX_
