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
#ifndef _STORM_SOLVER_IDRs_HXX_
#define _STORM_SOLVER_IDRs_HXX_

#include <vector>

#include <stormSolvers/stormSolver.hxx>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Solve a non-singular linear operator equation \
///   equation with the @c IDR(s) (Induced Dimension Reduction) method.
///
/// References:
/// @verbatim
/// [1] Peter Sonneveld, Martin B. van Gijzen. 
///     “IDR(s): A Family of Simple and Fast Algorithms for Solving 
///      Large Nonsymmetric Systems of Linear Equations.” 
///     SIAM J. Sci. Comput. 31 (2008): 1035-1062.
/// [2] Martin B. van Gijzen, Peter Sonneveld. 
///     “Algorithm 913: An elegant IDR(s) variant that efficiently 
///      exploits biorthogonality properties.” 
///     ACM Trans. Math. Softw. 38 (2011): 5:1-5:19.
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Vector>
class stormIdrsSolver final : public stormInnerOuterIterativeSolver<Vector> {
private:
  stormReal_t omega_;
  stormVector<stormReal_t> phi_, gamma_;
  stormMatrix<stormReal_t> mu_;
  Vector rVec_, vVec_, tVec_;
  std::vector<Vector> pVecs_, uVecs_, gVecs_;

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

public:

  stormIdrsSolver() {
    this->NumInnerIterations = 4;
  }

}; // class stormIdrsSolver<...>

template<class Vector>
void stormIdrsSolver<Vector>::OuterInit(Vector& xVec,
                                               Vector const& bVec,
                                               stormOperator<Vector> const& linOp,
                                               stormPreconditioner<Vector> const* preOp) {

  stormSize_t const s = this->NumInnerIterations;

  bool const leftPre = (preOp != nullptr) && 
    (this->PreSide == stormPreconditionerSide::Left);

  phi_.Assign(s);
  gamma_.Assign(s);
  mu_.Assign(s, s);

  rVec_.Assign(xVec, false);
  vVec_.Assign(xVec, false);
  tVec_.Assign(xVec, false);
  pVecs_.resize(s);
  for (Vector& pVec : pVecs_) {
    pVec.Assign(xVec, false);
  }
  uVecs_.resize(s);
  for (Vector& uVec : uVecs_) {
    uVec.Assign(xVec, false);
  }
  gVecs_.resize(s);
  for (Vector& gVec : gVecs_) {
    gVec.Assign(xVec, false);
  }

  // ----------------------
  // 𝒓 ← 𝓐𝒙,
  // 𝒓 ← 𝒃 - 𝒓,
  // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
  //   𝒗 ← 𝒓,
  //   𝒓 ← 𝓟𝒗.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  linOp.MatVec(rVec_, xVec);
  stormBlas::Sub(rVec_, bVec, rVec_);
  if (leftPre) {
    std::swap(vVec_, rVec_);
    preOp->MatVec(rVec_, vVec_);
  }

} // stormIdrsSolver<...>::OuterInit

template<class Vector>
stormReal_t stormIdrsSolver<Vector>::InnerInit(Vector& xVec,
                                               Vector const& bVec,
                                               stormOperator<Vector> const& linOp,
                                               stormPreconditioner<Vector> const* preOp) {

  stormSize_t const s = this->NumInnerIterations;

  // ----------------------
  // Build shadow space and initialize 𝜑:
  // 𝗶𝗳 𝘍𝘪𝘳𝘴𝘵𝘐𝘵𝘦𝘳𝘢𝘵𝘪𝘰𝘯:
  //   𝜔 ← 𝜇₀₀ ← 𝟣,
  //   𝜑₀ ← ‖𝒓‖,
  //   𝒑₀ ← 𝒓/𝜑₀,
  //   𝗳𝗼𝗿 𝑖 = 𝟣, 𝑠 - 𝟣 𝗱𝗼:
  //     𝜇ᵢᵢ ← 𝟣, 𝜑ᵢ ← 𝟢,
  //     𝒑ᵢ ← random, 
  //     𝗳𝗼𝗿 𝑗 = 𝟢, 𝑖 - 𝟣 𝗱𝗼:
  //       𝜇ᵢⱼ ← 𝟢,
  //       𝒑ᵢ ← 𝒑ᵢ - <𝒑ᵢ⋅𝒑ⱼ>⋅𝒑ⱼ,
  //     𝗲𝗻𝗱 𝗳𝗼𝗿
  //     𝒑ᵢ ← 𝒑ᵢ/‖𝒑ᵢ‖.
  //   𝗲𝗻𝗱 𝗳𝗼𝗿
  // 𝗲𝗹𝘀𝗲:
  //   𝗳𝗼𝗿 𝑖 = 𝟢, 𝑠 - 𝟣 𝗱𝗼:
  //     𝜑ᵢ ← <𝒑ᵢ⋅𝒓>.
  //   𝗲𝗻𝗱 𝗳𝗼𝗿
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  bool const firstIteration = this->Iteration == 0;
  if (firstIteration) {
    omega_ = mu_(0, 0) = 1.0;
    phi_(0) = stormBlas::Norm2(rVec_);
    stormBlas::Scale(pVecs_[0], rVec_, 1.0/phi_(0));
    for (stormSize_t i = 1; i < s; ++i) {
      mu_(i, i) = 1.0, phi_(i) = 0.0;
      stormBlas::RandFill(pVecs_[i]);
      for (stormSize_t j = 0; j < i; ++j) {
        mu_(i, j) = 0.0;
        stormBlas::Sub(pVecs_[i], pVecs_[i], 
          pVecs_[j], stormBlas::Dot(pVecs_[i], pVecs_[j]));
      }
      stormBlas::Scale(pVecs_[i], 
        pVecs_[i], 1.0/stormBlas::Norm2(pVecs_[i]));
    }
  } else {
    for (stormSize_t i = 0; i < s; ++i) {
      phi_(i) = stormBlas::Dot(pVecs_[i], rVec_);
    }
  }

  return phi_(0);

} // stormIdrsSolver<...>::InnerInit

template<class Vector>
stormReal_t stormIdrsSolver<Vector>::InnerIterate(Vector& xVec,
                                                  Vector const& bVec,
                                                  stormOperator<Vector> const& linOp,
                                                  stormPreconditioner<Vector> const* preOp) {

  stormSize_t const s = this->NumInnerIterations;
  stormSize_t const k = this->InnerIteration;

  bool const leftPre = (preOp != nullptr) && 
    (this->PreSide == stormPreconditionerSide::Left);
  bool const rightPre = (preOp != nullptr) && 
    (this->PreSide == stormPreconditionerSide::Right);

  // ----------------------
  // Compute 𝛾: 
  // {𝛾}ₖ:ₛ₋₁ ← ({𝜇}ₖ:ₛ₋₁,ₖ:ₛ₋₁)⁻¹{𝜑}ₖ:ₛ₋₁.
  // ----------------------
  for (stormSize_t i = k; i < s; ++i) {
    gamma_(i) = phi_(i);
  }
  for (stormSize_t i = k; i < s; ++i) {
    for (stormSize_t j = k; j < i; ++j) {
      gamma_(i) -= gamma_(j)*mu_(i, j);
    }
    gamma_(i) /= mu_(i, i);
  }

  // ----------------------
  // Compute new 𝒈ₖ and 𝒖ₖ:
  // 𝒗 ← 𝒓 - 𝛾ₖ⋅𝒈ₖ,
  // 𝗳𝗼𝗿 𝑖 = 𝑘 + 𝟣, 𝑠 - 𝟣 𝗱𝗼:
  //   𝒗 ← 𝒗 - 𝛾ᵢ⋅𝒈ᵢ.
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // 𝗶𝗳 𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦:
  //   𝒕 ← 𝒗,
  //   𝒗 ← 𝓟𝒕,
  // 𝗲𝗻𝗱 𝗶𝗳
  // 𝒖ₖ ← 𝜔⋅𝒗 + 𝛾ₖ⋅𝒖ₖ,
  // 𝗳𝗼𝗿 𝑖 = 𝑘 + 𝟣, 𝑠 - 𝟣 𝗱𝗼:
  //   𝒖ₖ ← 𝒖ₖ + 𝛾ᵢ⋅𝒖ᵢ.
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
  //   𝒈ₖ ← 𝓟(𝒕 ← 𝓐𝒖ₖ),
  // 𝗲𝗹𝘀𝗲:
  //   𝒈ₖ ← 𝓐𝒖ₖ,
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  stormBlas::Sub(vVec_, rVec_, gVecs_[k], gamma_(k));
  for (stormSize_t i = k + 1; i < s; ++i) {
    stormBlas::Sub(vVec_, vVec_, gVecs_[i], gamma_(i));
  }
  if (rightPre) {
    std::swap(tVec_, vVec_);
    preOp->MatVec(vVec_, tVec_);
  }
  stormBlas::Add(uVecs_[k], uVecs_[k], gamma_(k), vVec_, omega_);
  for (stormSize_t i = k + 1; i < s; ++i) {
    stormBlas::Add(uVecs_[k], uVecs_[k], uVecs_[i], gamma_(i));
  }
  if (leftPre) {
    stormBlas::MatVec(gVecs_[k], *preOp, tVec_, linOp, uVecs_[k]);
  } else {
    linOp.MatVec(gVecs_[k], uVecs_[k]);
  }

  // ----------------------
  // Bi-orthogonalize the new vectors:
  // 𝗳𝗼𝗿 𝑖 = 𝟢, 𝑘 - 𝟣 𝗱𝗼:
  //   𝛼 ← <𝒑ᵢ⋅𝒈ₖ>/𝜇ᵢᵢ,
  //   𝒖ₖ ← 𝒖ₖ - 𝛼⋅𝒖ᵢ,
  //   𝒈ₖ ← 𝒈ₖ - 𝛼⋅𝒈ᵢ.
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // ----------------------
  for (stormSize_t i = 0; i < k; ++i) {
    stormReal_t const alpha = 
      stormBlas::Dot(pVecs_[i], gVecs_[k])/mu_(i, i);
    stormBlas::Sub(uVecs_[k], uVecs_[k], uVecs_[i], alpha);
    stormBlas::Sub(gVecs_[k], gVecs_[k], gVecs_[i], alpha);
  }

  // ----------------------
  // Compute the new column of 𝜇: 
  // 𝗳𝗼𝗿 𝑖 = 𝑘, 𝑠 - 𝟣 𝗱𝗼:
  //   𝜇ᵢₖ ← <𝒑ᵢ⋅𝒈ₖ>.
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // ----------------------
  for (stormSize_t i = k; i < s; ++i) {
    mu_(i, k) = stormBlas::Dot(pVecs_[i], gVecs_[k]);
  }

  // ----------------------
  // Update the solution and the residual:
  // 𝛽 ← 𝜑ₖ/𝜇ₖₖ,
  // 𝒙 ← 𝒙 + 𝛽⋅𝒖ₖ,
  // 𝒓 ← 𝒓 - 𝛽⋅𝒈ₖ.
  // ----------------------
  stormReal_t const beta = phi_(k)/mu_(k, k);
  stormBlas::Add(xVec, xVec, uVecs_[k], beta);
  stormBlas::Sub(rVec_, rVec_, gVecs_[k], beta);

  // ----------------------
  // Update 𝜑:
  // 𝜑ₖ ← 𝟢,
  // 𝜑ₖ₊₁:ₛ₋₁ ← 𝜑ₖ₊₁:ₛ₋₁ - 𝛽⋅𝜇ₖ₊₁:ₛ₋₁,ₖ.
  // ----------------------
  phi_(k) = 0.0;
  for (stormSize_t i = k + 1; i < s; ++i) {
    phi_(i) -= beta*mu_(i, k);
  }

  // ----------------------
  // Enter the next 𝓖 subspace:
  // 𝗶𝗳 𝑘 = 𝑠 - 𝟣:
  //   𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
  //     𝒕 ← 𝓟(𝒗 ← 𝓐𝒓),
  //   𝗲𝗹𝘀𝗲 𝗶𝗳 𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦:
  //     𝒕 ← 𝓐(𝒗 ← 𝓟𝒓),
  //   𝗲𝗹𝘀𝗲:
  //     𝒕 ← 𝓐𝒓,
  //   𝗲𝗻𝗱 𝗶𝗳
  //   𝜔 ← <𝒕⋅𝒓>/<𝒕⋅𝒕>,
  //   𝒙 ← 𝒙 + 𝜔⋅(𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦 ? 𝒗 : 𝒓),
  //   𝒓 ← 𝒓 - 𝜔⋅𝒕.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  if (k == s - 1) {
    if (leftPre) {
      stormBlas::MatVec(tVec_, *preOp, vVec_, linOp, rVec_);
    } else if (rightPre) {
      stormBlas::MatVec(tVec_, linOp, vVec_, *preOp, rVec_);
    } else {
      linOp.MatVec(tVec_, rVec_);
    }
    omega_ = stormUtils::SafeDivide(
      stormBlas::Dot(tVec_, rVec_), stormBlas::Dot(tVec_, tVec_));
    stormBlas::Add(xVec, xVec, rightPre ? vVec_ : rVec_, omega_);
    stormBlas::Sub(rVec_, rVec_, tVec_, omega_);
  }

  return stormBlas::Norm2(rVec_);

} // stormIdrsSolver<...>::InnerIterate

#endif // ifndef _STORM_SOLVER_IDRs_HXX_