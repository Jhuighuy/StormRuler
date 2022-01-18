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
#include <algorithm>

#include <stormSolvers/stormSolver.hxx>
#include <stormBlas/stormTensor.hxx>

/// ----------------------------------------------------------------- ///
/// @brief Base class for @c GMRES and @c FGMRES.
/// ----------------------------------------------------------------- ///
template<bool Flexible, class Vector>
class stormBaseGmresSolver : public stormInnerOuterIterativeSolver<Vector> {
private:
  stormVector<stormReal_t> beta_, cs_, sn_;
  stormMatrix<stormReal_t> H_;
  std::vector<Vector> qVecs_;
  std::conditional_t<Flexible, std::vector<Vector>, std::array<Vector, 1>> zVecs_;

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
/// [1] Saad, Yousef and Martin H_. Schultz. 
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

  beta_.Assign(m + 1), cs_.Assign(m), sn_.Assign(m), H_.Assign(m + 1, m);
  qVecs_.resize(m + 1);
  for (Vector& qVec : qVecs_) {
    qVec.Assign(xVec, false);
  }
  if (preOp != nullptr) {
    if constexpr (Flexible) {
      zVecs_.resize(m);
    }
    for (Vector& zVec : zVecs_) {
      zVec.Assign(xVec, false);
    }
  }

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
  linOp.MatVec(qVecs_[0], xVec);
  stormBlas::Sub(qVecs_[0], bVec, qVecs_[0]);
  if (leftPre) {
    std::swap(zVecs_[0], qVecs_[0]);
    preOp->MatVec(qVecs_[0], zVecs_[0]);
  }

  // ----------------------
  // 𝛽₀ ← ‖𝒒₀‖,
  // 𝒒₀ ← 𝒒₀/𝛽₀. 
  // ----------------------
  beta_(0) = stormBlas::Norm2(qVecs_[0]);
  stormBlas::Scale(qVecs_[0], qVecs_[0], 1.0/beta_(0));

  return beta_(0);

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
  //   𝐻ᵢₖ ← <𝒒ₖ₊₁⋅𝒒ᵢ>,
  //   𝒒ₖ₊₁ ← 𝒒ₖ₊₁ - 𝐻ᵢₖ⋅𝒒ᵢ,
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // 𝐻ₖ₊₁,ₖ ← ‖𝒒ₖ₊₁‖, 
  // 𝒒ₖ₊₁ ← 𝒒ₖ₊₁/𝐻ₖ₊₁,ₖ.  
  // ----------------------
  if (leftPre) {
    stormBlas::MatVec(qVecs_[k + 1], *preOp, zVecs_[0], linOp, qVecs_[k]);
  } else if (rightPre) {
    stormSize_t const j = Flexible ? k : 0;
    stormBlas::MatVec(qVecs_[k + 1], linOp, zVecs_[j], *preOp, qVecs_[k]);
  } else {
    linOp.MatVec(qVecs_[k + 1], qVecs_[k]);
  }
  for (stormSize_t i = 0; i <= k; ++i) {
    H_(i, k) = stormBlas::Dot(qVecs_[k + 1], qVecs_[i]);
    stormBlas::Sub(qVecs_[k + 1], qVecs_[k + 1], qVecs_[i], H_(i, k));
  }
  H_(k + 1, k) = stormBlas::Norm2(qVecs_[k + 1]); 
  stormBlas::Scale(qVecs_[k + 1], qVecs_[k + 1], 1.0/H_(k + 1, k));

  // ----------------------
  // Eliminate the last element in {𝐻ᵢⱼ}
  // and and update the rotation matrix:
  // 𝗳𝗼𝗿 𝑖 = 𝟢, 𝑘 - 𝟣 𝗱𝗼:
  //   𝜒 ← 𝑐𝑠ᵢ⋅𝐻ᵢₖ + 𝑠𝑛ᵢ⋅𝐻ᵢ₊₁,ₖ,
  //   𝐻ᵢ₊₁,ₖ ← -𝑠𝑛ᵢ⋅𝐻ᵢₖ + 𝑐𝑠ᵢ⋅𝐻ᵢ₊₁,ₖ,
  //   𝐻ᵢₖ ← 𝜒,
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // 𝑐𝑠ₖ, 𝑠𝑛ₖ ← 𝘚𝘺𝘮𝘖𝘳𝘵𝘩𝘰(𝐻ₖₖ, 𝐻ₖ₊₁,ₖ),
  // 𝐻ₖₖ ← 𝑐𝑠ₖ⋅𝐻ₖₖ + 𝑠𝑛ₖ⋅𝐻ₖ₊₁,ₖ,
  // 𝐻ₖ₊₁,ₖ ← 𝟢.
  // ----------------------
  for (stormSize_t i = 0; i < k; ++i) {
    stormReal_t const chi = cs_(i)*H_(i, k) + sn_(i)*H_(i + 1, k);
    H_(i + 1, k) = -sn_(i)*H_(i, k) + cs_(i)*H_(i + 1, k);
    H_(i, k) = chi;
  }
  std::tie(cs_(k), sn_(k), std::ignore) =
    stormBlas::SymOrtho(H_(k, k), H_(k + 1, k));
  H_(k, k) = cs_(k)*H_(k, k) + sn_(k)*H_(k + 1, k);
  H_(k + 1, k) = 0.0;

  // ----------------------
  // Update the 𝛽-solution and residual norm:
  // 𝛽ₖ₊₁ ← -𝑠𝑛ₖ⋅𝛽ₖ, 𝛽ₖ ← 𝑐𝑠ₖ⋅𝛽ₖ.
  // ----------------------
  beta_(k + 1) = -sn_(k)*beta_(k), beta_(k) *= cs_(k);

  return std::abs(beta_(k + 1));

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
  // 𝛽ₖ ← 𝛽ₖ/𝐻ₖₖ,
  // 𝗳𝗼𝗿 𝑖 = 𝑘 - 𝟣, 𝟢, -𝟣 𝗱𝗼:
  //   𝛽ᵢ ← (𝛽ᵢ - <𝐻ᵢ,ᵢ₊₁:ₖ⋅𝛽ᵢ₊₁:ₖ>)/𝐻ᵢᵢ.
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // ----------------------
  /// @todo This should be replaced with a BLAS call that \
  ///   solves an equation with the upper-triangular part of 𝐻:
  ///   𝛽₀:ₖ ← (𝐻₀:ₖ,₀:ₖ)⁻¹𝛽₀:ₖ.
  beta_(k) /= H_(k, k);
  for (stormPtrDiff_t i = k - 1; i >= 0; --i) {
    for (stormSize_t j = i + 1; j <= k + 1; ++j) {
      beta_(i) -= H_(i, j)*beta_(j);
    }
    beta_(i) /= H_(i, i);
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
      stormBlas::Add(xVec, xVec, qVecs_[i], beta_(i));
    }
  } else if constexpr (Flexible) {
    for (stormSize_t i = 0; i <= k; ++i) {
      stormBlas::Add(xVec, xVec, zVecs_[i], beta_(i));
    }
  } else {
    stormBlas::Scale(qVecs_[0], qVecs_[0], beta_(0));
    for (stormSize_t i = 1; i <= k; ++i) {
      stormBlas::Add(qVecs_[0], qVecs_[0], qVecs_[i], beta_(i));
    }
    preOp->MatVec(zVecs_[0], qVecs_[0]);
    stormBlas::Add(xVec, xVec, zVecs_[0]);
  }

} // stormBaseGmresSolver<...>::InnerFinalize

#endif // ifndef _STORM_SOLVER_GMRES_HXX_
