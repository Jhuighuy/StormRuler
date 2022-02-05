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

#include <stormBlas/stormTensor.hxx>
#include <stormBlas/stormSubspace.hxx>
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
///     “Algorithm 913: An Elegant IDR(s) Variant that Efficiently
///      Exploits Biorthogonality Properties.”
///     ACM Trans. Math. Softw. 38 (2011): 5:1-5:19.
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Vector>
class stormIdrsSolver final : public stormInnerOuterIterativeSolver<Vector> {
private:
  stormReal_t omega_;
  stormVector<stormReal_t> phi_, gamma_;
  stormMatrix<stormReal_t> mu_;
  Vector rVec_, vVec_, zVec_;
  stormSubspace<Vector> pVecs_, uVecs_, gVecs_;

  stormReal_t OuterInit(Vector const& xVec,
                        Vector const& bVec,
                        stormOperator<Vector> const& linOp,
                        stormPreconditioner<Vector> const* preOp) override;

  void InnerInit(Vector const& xVec,
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
stormReal_t stormIdrsSolver<Vector>::OuterInit(Vector const& xVec,
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
  if (preOp != nullptr) {
    zVec_.Assign(xVec, false);
  }

  pVecs_.Assign(s, xVec, false);
  uVecs_.Assign(s, xVec, false);
  gVecs_.Assign(s, xVec, false);

  // ----------------------
  // 𝒓 ← 𝓐𝒙,
  // 𝒓 ← 𝒃 - 𝒓,
  // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
  //   𝒛 ← 𝒓,
  //   𝒓 ← 𝓟𝒛.
  // 𝗲𝗻𝗱 𝗶𝗳
  // 𝜑₀ ← ‖𝒓‖.
  // ----------------------
  linOp.MatVec(rVec_, xVec);
  rVec_.Sub(bVec, rVec_);
  if (leftPre) {
    std::swap(zVec_, rVec_);
    preOp->MatVec(rVec_, zVec_);
  }
  phi_(0) = rVec_.Norm2();

  return phi_(0);

} // stormIdrsSolver<...>::OuterInit

template<class Vector>
void stormIdrsSolver<Vector>::InnerInit(Vector const& xVec,
                                        Vector const& bVec,
                                        stormOperator<Vector> const& linOp,
                                        stormPreconditioner<Vector> const* preOp) {

  stormSize_t const s = this->NumInnerIterations;

  // ----------------------
  // Build shadow space and initialize 𝜑:
  // 𝗶𝗳 𝘍𝘪𝘳𝘴𝘵𝘐𝘵𝘦𝘳𝘢𝘵𝘪𝘰𝘯:
  //   𝜔 ← 𝜇₀₀ ← 𝟣,
  //   𝒑₀ ← 𝒓/𝜑₀,
  //   𝗳𝗼𝗿 𝑖 = 𝟣, 𝑠 - 𝟣 𝗱𝗼:
  //     𝜇ᵢᵢ ← 𝟣, 𝜑ᵢ ← 𝟢,
  //     𝒑ᵢ ← 𝘙𝘢𝘯𝘥𝘰𝘮,
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
    pVecs_(0).Scale(rVec_, 1.0/phi_(0));
    for (stormSize_t i = 1; i < s; ++i) {
      mu_(i, i) = 1.0, phi_(i) = 0.0;
      pVecs_(i).RandFill();
      for (stormSize_t j = 0; j < i; ++j) {
        mu_(i, j) = 0.0;
        pVecs_(i).Sub(pVecs_(i), pVecs_(j), stormBlas::Dot(pVecs_(i), pVecs_(j)));
      }
      pVecs_(i).Scale(1.0/pVecs_(i).Norm2());
    }
  } else {
    for (stormSize_t i = 0; i < s; ++i) {
      phi_(i) = stormBlas::Dot(pVecs_(i), rVec_);
    }
  }

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
  // 𝛾ₖ:ₛ₋₁ ← (𝜇ₖ:ₛ₋₁,ₖ:ₛ₋₁)⁻¹⋅𝜑ₖ:ₛ₋₁.
  // ----------------------
  for (stormSize_t i = k; i < s; ++i) {
    gamma_(i) = phi_(i);
    for (stormSize_t j = k; j < i; ++j) {
      gamma_(i) -= mu_(i, j)*gamma_(j);
    }
    gamma_(i) /= mu_(i, i);
  }

  // ----------------------
  // Compute the new 𝒈ₖ and 𝒖ₖ vectors:
  // 𝒗 ← 𝒓 - 𝛾ₖ⋅𝒈ₖ,
  // 𝗳𝗼𝗿 𝑖 = 𝑘 + 𝟣, 𝑠 - 𝟣 𝗱𝗼:
  //   𝒗 ← 𝒗 - 𝛾ᵢ⋅𝒈ᵢ,
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // 𝗶𝗳 𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦:
  //   𝒛 ← 𝒗,
  //   𝒗 ← 𝓟𝒛,
  // 𝗲𝗻𝗱 𝗶𝗳
  // 𝒖ₖ ← 𝜔⋅𝒗 + 𝛾ₖ⋅𝒖ₖ,
  // 𝗳𝗼𝗿 𝑖 = 𝑘 + 𝟣, 𝑠 - 𝟣 𝗱𝗼:
  //   𝒖ₖ ← 𝒖ₖ + 𝛾ᵢ⋅𝒖ᵢ,
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
  //   𝒈ₖ ← 𝓟(𝒛 ← 𝓐𝒖ₖ).
  // 𝗲𝗹𝘀𝗲:
  //   𝒈ₖ ← 𝓐𝒖ₖ.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  vVec_.Sub(rVec_, gVecs_(k), gamma_(k));
  for (stormSize_t i = k + 1; i < s; ++i) {
    vVec_.Sub(vVec_, gVecs_(i), gamma_(i));
  }
  if (rightPre) {
    std::swap(zVec_, vVec_);
    preOp->MatVec(vVec_, zVec_);
  }
  stormBlas::Add(uVecs_(k), uVecs_(k), gamma_(k), vVec_, omega_);
  for (stormSize_t i = k + 1; i < s; ++i) {
    stormBlas::Add(uVecs_(k), uVecs_(k), uVecs_(i), gamma_(i));
  }
  if (leftPre) {
    stormBlas::MatVec(gVecs_(k), *preOp, zVec_, linOp, uVecs_(k));
  } else {
    linOp.MatVec(gVecs_(k), uVecs_(k));
  }

  // ----------------------
  // Biorthogonalize the new 𝒈ₖ and 𝒖ₖ vectors:
  // 𝗳𝗼𝗿 𝑖 = 𝟢, 𝑘 - 𝟣 𝗱𝗼:
  //   𝛼 ← <𝒑ᵢ⋅𝒈ₖ>/𝜇ᵢᵢ,
  //   𝒖ₖ ← 𝒖ₖ - 𝛼⋅𝒖ᵢ,
  //   𝒈ₖ ← 𝒈ₖ - 𝛼⋅𝒈ᵢ.
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // ----------------------
  for (stormSize_t i = 0; i < k; ++i) {
    stormReal_t const alpha =
      stormUtils::SafeDivide(stormBlas::Dot(pVecs_(i), gVecs_(k)), mu_(i, i));
    uVecs_(k).Sub(uVecs_(k), uVecs_(i), alpha);
    gVecs_(k).Sub(gVecs_(k), gVecs_(i), alpha);
  }

  // ----------------------
  // Compute the new column of 𝜇:
  // 𝗳𝗼𝗿 𝑖 = 𝑘, 𝑠 - 𝟣 𝗱𝗼:
  //   𝜇ᵢₖ ← <𝒑ᵢ⋅𝒈ₖ>.
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // ----------------------
  for (stormSize_t i = k; i < s; ++i) {
    mu_(i, k) = stormBlas::Dot(pVecs_(i), gVecs_(k));
  }

  // ----------------------
  // Update the solution and the residual:
  // 𝛽 ← 𝜑ₖ/𝜇ₖₖ,
  // 𝒙 ← 𝒙 + 𝛽⋅𝒖ₖ,
  // 𝒓 ← 𝒓 - 𝛽⋅𝒈ₖ.
  // ----------------------
  stormReal_t const beta = 
    stormUtils::SafeDivide(phi_(k), mu_(k, k));
  stormBlas::Add(xVec, xVec, uVecs_(k), beta);
  rVec_.Sub(rVec_, gVecs_(k), beta);

  // ----------------------
  // Update 𝜑:
  // 𝜑ₖ₊₁:ₛ₋₁ ← 𝜑ₖ₊₁:ₛ₋₁ - 𝛽⋅𝜇ₖ₊₁:ₛ₋₁,ₖ.
  // ----------------------
  for (stormSize_t i = k + 1; i < s; ++i) {
    phi_(i) -= beta*mu_(i, k);
  }

  if (k == s - 1) {
    // ----------------------
    // Enter the next 𝓖 subspace:
    // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
    //   𝒗 ← 𝓟(𝒛 ← 𝓐𝒓),
    // 𝗲𝗹𝘀𝗲 𝗶𝗳 𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦:
    //   𝒗 ← 𝓐(𝒛 ← 𝓟𝒓),
    // 𝗲𝗹𝘀𝗲:
    //   𝒗 ← 𝓐𝒓,
    // 𝗲𝗻𝗱 𝗶𝗳
    // 𝜔 ← <𝒗⋅𝒓>/<𝒗⋅𝒗>,
    // 𝒙 ← 𝒙 + 𝜔⋅(𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦 ? 𝒛 : 𝒓),
    // 𝒓 ← 𝒓 - 𝜔⋅𝒗.
    // ----------------------
    if (leftPre) {
      stormBlas::MatVec(vVec_, *preOp, zVec_, linOp, rVec_);
    } else if (rightPre) {
      stormBlas::MatVec(vVec_, linOp, zVec_, *preOp, rVec_);
    } else {
      linOp.MatVec(vVec_, rVec_);
    }
    omega_ = stormUtils::SafeDivide(
      stormBlas::Dot(vVec_, rVec_), stormBlas::Dot(vVec_, vVec_));
    stormBlas::Add(xVec, xVec, rightPre ? zVec_ : rVec_, omega_);
    rVec_.Sub(rVec_, vVec_, omega_);
  }

  return rVec_.Norm2();

} // stormIdrsSolver<...>::InnerIterate

#endif // ifndef _STORM_SOLVER_IDRs_HXX_
