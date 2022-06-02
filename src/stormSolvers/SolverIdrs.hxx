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

#pragma once

#include <stormBase.hxx>
#include <stormSolvers/LegacyTensor.hxx>
#include <stormSolvers/Solver.hxx>
#include <stormSolvers/Subspace.hxx>
#include <stormSolvers/Vector.hxx>

namespace Storm {

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief The @c IDR(s) (Induced Dimension Reduction)
///   linear operator equation solver.
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
template<vector_like Vector>
class IdrsSolver final : public InnerOuterIterativeSolver<Vector> {
private:

  real_t omega_;
  stormVector<real_t> phi_, gamma_;
  stormMatrix<real_t> mu_;
  Vector rVec_, vVec_, zVec_;
  Subspace<Vector> pVecs_, uVecs_, gVecs_;

  real_t OuterInit(Vector const& xVec, Vector const& bVec,
                   Operator<Vector> const& linOp,
                   Preconditioner<Vector> const* preOp) override;

  void InnerInit(Vector const& xVec, Vector const& bVec,
                 Operator<Vector> const& linOp,
                 Preconditioner<Vector> const* preOp) override;

  real_t InnerIterate(Vector& xVec, Vector const& bVec,
                      Operator<Vector> const& linOp,
                      Preconditioner<Vector> const* preOp) override;

public:

  IdrsSolver() {
    this->NumInnerIterations = 4;
  }

}; // class IdrsSolver

template<vector_like Vector>
real_t IdrsSolver<Vector>::OuterInit(Vector const& xVec, Vector const& bVec,
                                     Operator<Vector> const& linOp,
                                     Preconditioner<Vector> const* preOp) {
  size_t const s{this->NumInnerIterations};

  bool const leftPre{(preOp != nullptr) &&
                     (this->PreSide == PreconditionerSide::Left)};

  phi_.Assign(s);
  gamma_.Assign(s);
  mu_.Assign(s, s);

  rVec_.Assign(xVec, false);
  vVec_.Assign(xVec, false);
  if (preOp != nullptr) { zVec_.Assign(xVec, false); }

  pVecs_.Assign(s, xVec, false);
  uVecs_.Assign(s, xVec, false);
  gVecs_.Assign(s, xVec, false);

  // Initialize:
  // ----------------------
  // 𝒓 ← 𝒃 - 𝓐𝒙,
  // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
  //   𝒛 ← 𝒓,
  //   𝒓 ← 𝓟𝒛.
  // 𝗲𝗻𝗱 𝗶𝗳
  // 𝜑₀ ← ‖𝒓‖.
  // ----------------------
  linOp.Residual(rVec_, bVec, xVec);
  if (leftPre) {
    zVec_.Swap(rVec_);
    preOp->MatVec(rVec_, zVec_);
  }
  phi_(0) = rVec_.Norm2();

  return phi_(0);

} // IdrsSolver::OuterInit

template<vector_like Vector>
void IdrsSolver<Vector>::InnerInit(Vector const& xVec, Vector const& bVec,
                                   Operator<Vector> const& linOp,
                                   Preconditioner<Vector> const* preOp) {
  size_t const s{this->NumInnerIterations};

  // Build shadow space and initialize 𝜑:
  // ----------------------
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
  bool const firstIteration{this->Iteration == 0};
  if (firstIteration) {
    omega_ = mu_(0, 0) = 1.0;
    pVecs_(0).Scale(rVec_, 1.0 / phi_(0));
    for (size_t i{1}; i < s; ++i) {
      mu_(i, i) = 1.0, phi_(i) = 0.0;
      pVecs_(i).RandFill();
      for (size_t j{0}; j < i; ++j) {
        mu_(i, j) = 0.0;
        pVecs_(i).SubAssign(pVecs_(j), pVecs_(i).Dot(pVecs_(j)));
      }
      pVecs_(i).ScaleAssign(1.0 / pVecs_(i).Norm2());
    }
  } else {
    for (size_t i{0}; i < s; ++i) {
      phi_(i) = pVecs_(i).Dot(rVec_);
    }
  }

} // IdrsSolver::InnerInit

template<vector_like Vector>
real_t IdrsSolver<Vector>::InnerIterate(Vector& xVec, Vector const& bVec,
                                        Operator<Vector> const& linOp,
                                        Preconditioner<Vector> const* preOp) {
  size_t const s{this->NumInnerIterations};
  size_t const k{this->InnerIteration};

  bool const leftPre{(preOp != nullptr) &&
                     (this->PreSide == PreconditionerSide::Left)};
  bool const rightPre{(preOp != nullptr) &&
                      (this->PreSide == PreconditionerSide::Right)};

  // Compute 𝛾:
  // ----------------------
  // 𝛾ₖ:ₛ₋₁ ← (𝜇ₖ:ₛ₋₁,ₖ:ₛ₋₁)⁻¹⋅𝜑ₖ:ₛ₋₁.
  // ----------------------
  for (size_t i = k; i < s; ++i) {
    gamma_(i) = phi_(i);
    for (size_t j = k; j < i; ++j) {
      gamma_(i) -= mu_(i, j) * gamma_(j);
    }
    gamma_(i) /= mu_(i, i);
  }

  // Compute the new 𝒈ₖ and 𝒖ₖ vectors:
  // ----------------------
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
  for (size_t i{k + 1}; i < s; ++i) {
    vVec_.SubAssign(gVecs_(i), gamma_(i));
  }
  if (rightPre) {
    std::swap(zVec_, vVec_);
    preOp->MatVec(vVec_, zVec_);
  }
  uVecs_(k).Add(uVecs_(k), gamma_(k), vVec_, omega_);
  for (size_t i{k + 1}; i < s; ++i) {
    uVecs_(k).AddAssign(uVecs_(i), gamma_(i));
  }
  if (leftPre) {
    preOp->MatVec(gVecs_(k), zVec_, linOp, uVecs_(k));
  } else {
    linOp.MatVec(gVecs_(k), uVecs_(k));
  }

  // Biorthogonalize the new vectors 𝒈ₖ and 𝒖ₖ:
  // ----------------------
  // 𝗳𝗼𝗿 𝑖 = 𝟢, 𝑘 - 𝟣 𝗱𝗼:
  //   𝛼 ← <𝒑ᵢ⋅𝒈ₖ>/𝜇ᵢᵢ,
  //   𝒖ₖ ← 𝒖ₖ - 𝛼⋅𝒖ᵢ,
  //   𝒈ₖ ← 𝒈ₖ - 𝛼⋅𝒈ᵢ.
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // ----------------------
  for (size_t i{0}; i < k; ++i) {
    real_t const alpha{Utils::SafeDivide(pVecs_(i).Dot(gVecs_(k)), mu_(i, i))};
    uVecs_(k).SubAssign(uVecs_(i), alpha);
    gVecs_(k).SubAssign(gVecs_(i), alpha);
  }

  // Compute the new column of 𝜇:
  // ----------------------
  // 𝗳𝗼𝗿 𝑖 = 𝑘, 𝑠 - 𝟣 𝗱𝗼:
  //   𝜇ᵢₖ ← <𝒑ᵢ⋅𝒈ₖ>.
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // ----------------------
  for (size_t i{k}; i < s; ++i) {
    mu_(i, k) = pVecs_(i).Dot(gVecs_(k));
  }

  // Update the solution and the residual:
  // ----------------------
  // 𝛽 ← 𝜑ₖ/𝜇ₖₖ,
  // 𝒙 ← 𝒙 + 𝛽⋅𝒖ₖ,
  // 𝒓 ← 𝒓 - 𝛽⋅𝒈ₖ.
  // ----------------------
  real_t const beta{Utils::SafeDivide(phi_(k), mu_(k, k))};
  xVec.AddAssign(uVecs_(k), beta);
  rVec_.SubAssign(gVecs_(k), beta);

  // Update 𝜑:
  // ----------------------
  // 𝜑ₖ₊₁:ₛ₋₁ ← 𝜑ₖ₊₁:ₛ₋₁ - 𝛽⋅𝜇ₖ₊₁:ₛ₋₁,ₖ.
  // ----------------------
  for (size_t i{k + 1}; i < s; ++i) {
    phi_(i) -= beta * mu_(i, k);
  }

  if (k == s - 1) {
    // Enter the next 𝓖 subspace:
    // ----------------------
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
      preOp->MatVec(vVec_, zVec_, linOp, rVec_);
    } else if (rightPre) {
      linOp.MatVec(vVec_, zVec_, *preOp, rVec_);
    } else {
      linOp.MatVec(vVec_, rVec_);
    }
    omega_ = Utils::SafeDivide(vVec_.Dot(rVec_), vVec_.Dot(vVec_));
    xVec.AddAssign(rightPre ? zVec_ : rVec_, omega_);
    rVec_.SubAssign(vVec_, omega_);
  }

  return rVec_.Norm2();

} // IdrsSolver::InnerIterate

} // namespace Storm
