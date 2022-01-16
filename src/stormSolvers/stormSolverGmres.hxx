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
#include <stormBlas/stormArray.hxx>
#include <stormBlas/stormTensor.hxx>

/// ----------------------------------------------------------------- ///
/// @brief Base class for @c GMRES and @c FGMRES.
/// ----------------------------------------------------------------- ///
template<bool Flexible, class Vector>
class stormBaseGmresSolver : public stormInnerOuterIterativeSolver<Vector> {
private:
  std::vector<stormReal_t> betaData, csData, snData;
  std::vector<Vector> qVec_;
  std::conditional_t<Flexible, std::vector<Vector>, std::array<Vector, 1>> zVec_;
  stormVectorView<stormReal_t> beta_, cs_, sn_;
  stormMatrix<stormReal_t> h_;

  void OuterInit(Vector& xVec,
                 Vector const& bVec,
                 stormOperator<Vector> const& linOp,
                 stormPreconditioner<Vector> const* preOp) override;

  stormReal_t InnerInit(Vector& xVec,
                        Vector const& bVec,
                        stormOperator<Vector> const& linOp,
                        stormPreconditioner<Vector> const* preOp) override;

  stormReal_t InnerIterate(Vector& xVec,
                           Vector const& bVec,
                           stormOperator<Vector> const& linOp,
                           stormPreconditioner<Vector> const* preOp) override;

  void InnerFinalize(Vector& xVec,
                     Vector const& bVec,
                     stormOperator<Vector> const& linOp,
                     stormPreconditioner<Vector> const* preOp) override;

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
/// @c GMRES is algebraically equivalent to @c MINRES method \
///   in the self-adjoint operator unpreconditioned case, \
///   however, the need for restarts may lead to the much slower \
///   @c GMRES convergence rate. 
///
/// @c GMRES may be applied to the singular problems, and the square \
///   least squares problems, although, similarly to @c MINRES, \
///   convergeance to minimum norm solution is not guaranteed.
///
/// References:
/// @verbatim
/// [1] Saad, Yousef and Martin h_. Schultz. 
///     “GMRES: A generalized minimal residual algorithm for solving 
///      nonsymmetric linear systems.” 
///     SIAM J. Sci. Stat. Comput., 7:856–869, 1986.
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Vector>
class stormGmresSolver final : public stormBaseGmresSolver<false, Vector> {

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
///   usage. For the static preconditioners, @c FGMRES requires \
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
template<class Vector>
class stormFgmresSolver final : public stormBaseGmresSolver<true, Vector> {

}; // class stormFgmresSolver<...>

template<bool Flexible, class Vector>
void stormBaseGmresSolver<Flexible, Vector>::
                        OuterInit(Vector& xVec,
                                  Vector const& bVec,
                                  stormOperator<Vector> const& linOp,
                                  stormPreconditioner<Vector> const* preOp) {

  stormSize_t const m = this->NumInnerIterations;
  betaData.resize(m + 1), csData.resize(m), snData.resize(m);
  h_.Assign(m + 1, m);
  qVec_.resize(m + 1);
  for (Vector& qVec_ : qVec_) {
    stormUtils::AllocLike(xVec, qVec_);
  }
  if (preOp != nullptr) {
    if constexpr (Flexible) {
      zVec_.resize(m);
    }
    for (Vector& zVec_ : zVec_) {
      stormUtils::AllocLike(xVec, zVec_);
    }
  }

  beta_.Assign(betaData.data(), m + 1); 
  cs_.Assign(csData.data(), m); 
  sn_.Assign(snData.data(), m);

} // stormBaseGmresSolver<...>::OuterInit

template<bool Flexible, class Vector>
stormReal_t stormBaseGmresSolver<Flexible, Vector>::
                               InnerInit(Vector& xVec,
                                         Vector const& bVec,
                                         stormOperator<Vector> const& linOp,
                                         stormPreconditioner<Vector> const* preOp) {

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
  linOp.MatVec(qVec_[0], xVec);
  stormBlas::Sub(qVec_[0], bVec, qVec_[0]);
  if (leftPre) {
    std::swap(zVec_[0], qVec_[0]);
    preOp->MatVec(qVec_[0], zVec_[0]);
  }

  // ----------------------
  // 𝑐𝑠 ← {𝟢}ᵀ, 𝑠𝑛 ← {𝟢}ᵀ,
  // 𝜑 ← ‖𝒒₀‖,
  // 𝛽 ← {𝜑,𝟢,…,𝟢}ᵀ,
  // 𝒒₀ ← 𝒒₀/𝜑. 
  // ----------------------
  std::fill(csData.begin(), csData.end(), 0.0);
  std::fill(snData.begin(), snData.end(), 0.0);
  stormReal_t const phi = stormBlas::Norm2(qVec_[0]);
  beta_(0) = phi, std::fill(betaData.begin() + 1, betaData.end(), 0.0);
  stormBlas::Scale(qVec_[0], qVec_[0], 1.0/phi);

  return phi;

} // stormBaseGmresSolver<...>::InnerInit

template<bool Flexible, class Vector>
stormReal_t stormBaseGmresSolver<Flexible, Vector>::
                            InnerIterate(Vector& xVec,
                                         Vector const& bVec,
                                         stormOperator<Vector> const& linOp,
                                         stormPreconditioner<Vector> const* preOp) {

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
    stormBlas::MatVec(qVec_[k + 1], *preOp, zVec_[0], linOp, qVec_[k]);
  } else if (rightPre) {
    stormSize_t const j = Flexible ? k : 0;
    stormBlas::MatVec(qVec_[k + 1], linOp, zVec_[j], *preOp, qVec_[k]);
  } else {
    linOp.MatVec(qVec_[k + 1], qVec_[k]);
  }
  for (stormSize_t i = 0; i <= k; ++i) {
    h_(i, k) = stormBlas::Dot(qVec_[k + 1], qVec_[i]);
    stormBlas::Sub(qVec_[k + 1], qVec_[k + 1], qVec_[i], h_(i, k));
  }
  h_(k + 1, k) = stormBlas::Norm2(qVec_[k + 1]); 
  stormBlas::Scale(qVec_[k + 1], qVec_[k + 1], 1.0/h_(k + 1, k));

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
    stormReal_t const chi = cs_(i)*h_(i, k) + sn_(i)*h_(i + 1, k);
    h_(i + 1, k) = -sn_(i)*h_(i, k) + cs_(i)*h_(i + 1, k);
    h_(i, k) = chi;
  }
  std::tie(cs_(k), sn_(k), std::ignore) =
    stormBlas::SymOrtho(h_(k, k), h_(k + 1, k));
  h_(k, k) = cs_(k)*h_(k, k) + sn_(k)*h_(k + 1, k);
  h_(k + 1, k) = 0.0;

  // ----------------------
  // Update the 𝛽-solution and residual norm:
  // 𝛽ₖ₊₁ ← -𝑠𝑛ₖ⋅𝛽ₖ, 𝛽ₖ ← 𝑐𝑠ₖ⋅𝛽ₖ,
  // 𝜑 ← |𝛽ₖ₊₁|.
  // ----------------------
  beta_(k + 1) = -sn_(k)*beta_(k), beta_(k) *= cs_(k);
  stormReal_t const phi = std::abs(beta_(k + 1));

  return phi;

} // stormBaseGmresSolver<...>::InnerIterate

template<bool Flexible, class Vector>
void stormBaseGmresSolver<Flexible, Vector>::
                    InnerFinalize(Vector& xVec,
                                  Vector const& bVec,
                                  stormOperator<Vector> const& linOp,
                                  stormPreconditioner<Vector> const* preOp) {

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
  beta_(k) /= h_(k, k);
  for (stormPtrDiff_t i = k - 1; i >= 0; --i) {
    //beta_(i) -= std::inner_product(
    //  beta_.begin() + i + 1, beta_.begin() + k + 1, h_[i].begin() + i + 1, 0.0);
    for (stormSize_t j = i + 1; j <= k + 1; ++j) beta_(i) -= h_(i, j)*beta_(j);
    beta_(i) /= h_(i, i);
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
      stormBlas::Add(xVec, xVec, qVec_[i], beta_(i));
    }
  } else if constexpr (Flexible) {
    for (stormSize_t i = 0; i <= k; ++i) {
      stormBlas::Add(xVec, xVec, zVec_[i], beta_(i));
    }
  } else {
    stormBlas::Scale(qVec_[0], qVec_[0], beta_(0));
    for (stormSize_t i = 1; i <= k; ++i) {
      stormBlas::Add(qVec_[0], qVec_[0], qVec_[i], beta_(i));
    }
    preOp->MatVec(zVec_[0], qVec_[0]);
    stormBlas::Add(xVec, xVec, zVec_[0]);
  }

} // stormBaseGmresSolver<...>::InnerFinalize

#endif // ifndef _STORM_SOLVER_GMRES_HXX_
