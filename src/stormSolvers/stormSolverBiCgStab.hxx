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
#ifndef _STORM_SOLVER_BICGSTAB_
#define _STORM_SOLVER_BICGSTAB_

#include <cmath>

#include <stormBlas/stormTensor.hxx>
#include <stormBlas/stormSubspace.hxx>
#include <stormSolvers/stormSolver.hxx>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Solve a linear operator equation with the good old \
///   @c BiCGStab (Biconjugate Gradients Stabilized) method.
///
/// @c BiCGStab, like the other @c BiCG type solvers, requires \
///   two operator multiplications per iteration.
///
/// @c BiCGStab typically converges much smoother, than \
///   @c CGS. @todo Breakdowns?
///
/// References:
/// @verbatim
/// [1] Henk A. van der Vorst.
///     “Bi-CGSTAB: A Fast and Smoothly Converging Variant of Bi-CG
///      for the Solution of Nonsymmetric Linear Systems.”
///     SIAM J. Sci. Comput. 13 (1992): 631-644.
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Vector>
class stormBiCgStabSolver final : public stormIterativeSolver<Vector> {
private:
  stormReal_t alpha_, rho_, omega_;
  Vector pVec_, rVec_, rTildeVec_, tVec_, vVec_, zVec_;

  stormReal_t Init(Vector const& xVec,
                   Vector const& bVec,
                   stormOperator<Vector> const& linOp,
                   stormPreconditioner<Vector> const* preOp) override;

  stormReal_t Iterate(Vector& xVec,
                      Vector const& bVec,
                      stormOperator<Vector> const& linOp,
                      stormPreconditioner<Vector> const* preOp) override;

}; // class stormBiCgStabSolver<...>

template<class Vector>
stormReal_t stormBiCgStabSolver<Vector>::Init(Vector const& xVec,
                                              Vector const& bVec,
                                              stormOperator<Vector> const& linOp,
                                              stormPreconditioner<Vector> const* preOp) {

  bool const leftPre = (preOp != nullptr) &&
    (this->PreSide == stormPreconditionerSide::Left);

  pVec_.Assign(xVec, false);
  rVec_.Assign(xVec, false);
  rTildeVec_.Assign(xVec, false);
  tVec_.Assign(xVec, false);
  vVec_.Assign(xVec, false);
  if (preOp != nullptr) {
    zVec_.Assign(xVec, false);
  }

  // ----------------------
  // 𝒓 ← 𝓐𝒙,
  // 𝒓 ← 𝒃 - 𝒓,
  // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
  //   𝒛 ← 𝒓,
  //   𝒓 ← 𝓟𝒛,
  // 𝗲𝗻𝗱 𝗶𝗳
  // 𝒓̃ ← 𝒓,
  // 𝜌 ← <𝒓̃⋅𝒓>.
  // ----------------------
  linOp.MatVec(rVec_, xVec);
  stormBlas::Sub(rVec_, bVec, rVec_);
  if (leftPre) {
    std::swap(zVec_, rVec_);
    preOp->MatVec(rVec_, zVec_);
  }
  stormBlas::Set(rTildeVec_, rVec_);
  rho_ = stormBlas::Dot(rTildeVec_, rVec_);

  return std::sqrt(rho_);

} // stormBiCgStabSolver<...>::Init

template<class Vector>
stormReal_t stormBiCgStabSolver<Vector>::Iterate(Vector& xVec,
                                                 Vector const& bVec,
                                                 stormOperator<Vector> const& linOp,
                                                 stormPreconditioner<Vector> const* preOp) {

  bool const leftPre = (preOp != nullptr) &&
    (this->PreSide == stormPreconditionerSide::Left);
  bool const rightPre = (preOp != nullptr) &&
    (this->PreSide == stormPreconditionerSide::Right);

  // ----------------------
  // Continue the iterations:
  // 𝗶𝗳 𝘍𝘪𝘳𝘴𝘵𝘐𝘵𝘦𝘳𝘢𝘵𝘪𝘰𝘯:
  //   𝒑 ← 𝒓.
  // 𝗲𝗹𝘀𝗲:
  //   𝜌̅ ← 𝜌,
  //   𝜌 ← <𝒓̃⋅𝒓>,
  //   𝛽 ← (𝜌/𝜌̅)⋅(𝛼/𝜔),
  //   𝒑 ← 𝒑 - 𝜔⋅𝒗,
  //   𝒑 ← 𝒓 + 𝛽⋅𝒑.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  bool const firstIteration = this->Iteration == 0;
  if (firstIteration) {
    stormBlas::Set(pVec_, rVec_);
  } else {
    stormReal_t const rhoBar = rho_;
    rho_ = stormBlas::Dot(rTildeVec_, rVec_);
    stormReal_t const beta =
      stormUtils::SafeDivide(rho_, rhoBar)*stormUtils::SafeDivide(alpha_, omega_);
    stormBlas::Sub(pVec_, pVec_, vVec_, omega_);
    stormBlas::Add(pVec_, rVec_, pVec_, beta);
  }

  // ----------------------
  // Update the solution and the residual:
  // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
  //   𝒗 ← 𝓟(𝒛 ← 𝓐𝒑),
  // 𝗲𝗹𝘀𝗲 𝗶𝗳 𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦:
  //   𝒗 ← 𝓐(𝒛 ← 𝓟𝒑),
  // 𝗲𝗹𝘀𝗲:
  //   𝒗 ← 𝓐𝒑,
  // 𝗲𝗻𝗱 𝗶𝗳
  // 𝛼 ← 𝜌/<𝒓̃⋅𝒗>,
  // 𝒙 ← 𝒙 + 𝛼⋅(𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦 ? 𝒛 : 𝒑),
  // 𝒓 ← 𝒓 - 𝛼⋅𝒗.
  // ----------------------
  if (leftPre) {
    stormBlas::MatVec(vVec_, *preOp, zVec_, linOp, pVec_);
  } else if (rightPre) {
    stormBlas::MatVec(vVec_, linOp, zVec_, *preOp, pVec_);
  } else {
    linOp.MatVec(vVec_, pVec_);
  }
  alpha_ = stormUtils::SafeDivide(rho_, stormBlas::Dot(rTildeVec_, vVec_));
  stormBlas::Add(xVec, xVec, rightPre ? zVec_ : pVec_, alpha_);
  stormBlas::Sub(rVec_, rVec_, vVec_, alpha_);

  // ----------------------
  // Update the solution and the residual again:
  // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
  //   𝒕 ← 𝓟(𝒛 ← 𝓐𝒓),
  // 𝗲𝗹𝘀𝗲 𝗶𝗳 𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦:
  //   𝒕 ← 𝓐(𝒛 ← 𝓟𝒓),
  // 𝗲𝗹𝘀𝗲:
  //   𝒕 ← 𝓐𝒓,
  // 𝗲𝗻𝗱 𝗶𝗳
  // 𝜔 ← <𝒕⋅𝒓>/<𝒕⋅𝒕>,
  // 𝒙 ← 𝒙 + 𝜔⋅(𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦 ? 𝒛 : 𝒓),
  // 𝒓 ← 𝒓 - 𝜔⋅𝒕.
  // ----------------------
  if (leftPre) {
    stormBlas::MatVec(tVec_, *preOp, zVec_, linOp, rVec_);
  } else if (rightPre) {
    stormBlas::MatVec(tVec_, linOp, zVec_, *preOp, rVec_);
  } else {
    linOp.MatVec(tVec_, rVec_);
  }
  omega_ = stormUtils::SafeDivide(
    stormBlas::Dot(tVec_, rVec_), stormBlas::Dot(tVec_, tVec_));
  stormBlas::Add(xVec, xVec, rightPre ? zVec_ : rVec_, omega_);
  stormBlas::Sub(rVec_, rVec_, tVec_, omega_);

  return stormBlas::Norm2(rVec_);

} // stormBiCgStabSolver<...>::Iterate

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Solve a linear operator equation with the \
///   @c BiCGStab(l) (Biconjugate Gradients Stabilized) method.
///
/// @c BiCGStab(l), like the other @c BiCG type solvers, requires \
///   two operator multiplications per iteration.
///
/// References:
/// @verbatim
/// [1] Gerard L. G. Sleijpen and Diederik R. Fokkema. 
///     “BiCGStab(l) for Linear Equations involving 
///      Unsymmetric Matrices with Complex Spectrum.” 
///     Electronic Transactions on Numerical Analysis 1 (1993): 11-32.
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Vector>
class stormBiCGStabLSolver final : public stormInnerOuterIterativeSolver<Vector> {
private:
  stormReal_t alpha_, rho_, omega_;
  stormVector<stormReal_t> gamma_, gammaBar_, gammaBarBar_, sigma_;
  stormMatrix<stormReal_t> tau_;
  Vector rTildeVec_, zVec_;
  stormSubspace<Vector> rVecs_, uVecs_;

  stormReal_t OuterInit(Vector const& xVec,
                        Vector const& bVec,
                        stormOperator<Vector> const& linOp,
                        stormPreconditioner<Vector> const* preOp) override;

  stormReal_t InnerIterate(Vector& xVec,
                           Vector const& bVec,
                           stormOperator<Vector> const& linOp,
                           stormPreconditioner<Vector> const* preOp) override;

public:

  stormBiCGStabLSolver() {
    this->NumInnerIterations = 2;
  }

}; // class stormBiCGStabLSolver<...>

template<class Vector>
stormReal_t stormBiCGStabLSolver<Vector>::
                OuterInit(Vector const& xVec,
                          Vector const& bVec,
                          stormOperator<Vector> const& linOp,
                          stormPreconditioner<Vector> const* preOp) {

  stormSize_t const l = this->NumInnerIterations;

  gamma_.Assign(l + 1);
  gammaBar_.Assign(l + 1);
  gammaBarBar_.Assign(l + 1);
  sigma_.Assign(l + 1);
  tau_.Assign(l + 1, l + 1);

  rTildeVec_.Assign(xVec, false);
  if (preOp != nullptr) {
    zVec_.Assign(xVec, false);
  }

  rVecs_.Assign(l + 1, xVec, false);
  uVecs_.Assign(l + 1, xVec, false);

  // ----------------------
  // 𝒖₀ ← {𝟢}ᵀ,
  // 𝒓₀ ← 𝓐𝒙,
  // 𝒓₀ ← 𝒃 - 𝒓₀,
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝒛 ← 𝒓₀,
  //   𝒓₀ ← 𝓟𝒛,
  // 𝗲𝗻𝗱 𝗶𝗳
  // 𝒓̃ ← 𝒓₀,
  // 𝜌 ← <𝒓̃⋅𝒓₀>.
  // ----------------------
  stormBlas::Fill(uVecs_(0), 0.0);
  linOp.MatVec(rVecs_(0), xVec);
  stormBlas::Sub(rVecs_(0), bVec, rVecs_(0));
  if (preOp != nullptr) {
    std::swap(zVec_, rVecs_(0));
    preOp->MatVec(rVecs_(0), zVec_);
  }
  stormBlas::Set(rTildeVec_, rVecs_(0));
  rho_ = stormBlas::Dot(rTildeVec_, rVecs_(0));

  return std::sqrt(rho_);

} // stormBiCGStabLSolver<...>::OuterInit

template<class Vector>
stormReal_t stormBiCGStabLSolver<Vector>::
                    InnerIterate(Vector& xVec,
                                 Vector const& bVec,
                                 stormOperator<Vector> const& linOp,
                                 stormPreconditioner<Vector> const* preOp) {

  stormSize_t const l = this->NumInnerIterations;
  stormSize_t const j = this->InnerIteration;

  // ----------------------
  // BiCG part:
  // 𝗶𝗳 𝘍𝘪𝘳𝘴𝘵𝘐𝘵𝘦𝘳𝘢𝘵𝘪𝘰𝘯:
  //   𝒖₀ ← 𝒓₀,
  // 𝗲𝗹𝘀𝗲:
  //   𝜌̅ ← 𝜌,
  //   𝜌 ← <𝒓̃⋅𝒓ⱼ>,
  //   𝛽 ← 𝛼⋅𝜌/𝜌̅,
  //   𝗳𝗼𝗿 𝑖 = 𝟢, 𝑗 𝗱𝗼:
  //     𝒖ᵢ ← 𝒓ᵢ - 𝛽⋅𝒖ᵢ,
  //   𝗲𝗻𝗱 𝗳𝗼𝗿
  // 𝗲𝗻𝗱 𝗶𝗳
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝒖ⱼ₊₁ ← 𝓟(𝒛 ← 𝓐𝒖ⱼ),
  // 𝗲𝗹𝘀𝗲:
  //   𝒖ⱼ₊₁ ← 𝓐𝒖ⱼ,
  // 𝗲𝗻𝗱 𝗶𝗳
  // 𝛼 ← 𝜌/<𝒓̃⋅𝒖ⱼ₊₁>,
  // 𝗳𝗼𝗿 𝑖 = 𝟢, 𝑗 𝗱𝗼:
  //   𝒓ᵢ ← 𝒓ᵢ - 𝛼⋅𝒖ᵢ₊₁.
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // ----------------------
  bool const firstIteration = this->Iteration == 0;
  if (firstIteration) {
    stormBlas::Set(uVecs_(0), rVecs_(0));
  } else {
    stormReal_t const rhoBar = rho_;
    rho_ = stormBlas::Dot(rTildeVec_, rVecs_(j));
    stormReal_t const beta =
      alpha_*stormUtils::SafeDivide(rho_, rhoBar);
    for (stormSize_t i = 0; i <= j; ++i) {
      stormBlas::Sub(uVecs_(i), rVecs_(i), uVecs_(i), beta);
    }
  }
  if (preOp != nullptr) {
    stormBlas::MatVec(uVecs_(j + 1), *preOp, zVec_, linOp, uVecs_(j));
  } else {
    linOp.MatVec(uVecs_(j + 1), uVecs_(j));
  }
  alpha_ = stormUtils::SafeDivide(rho_, stormBlas::Dot(rTildeVec_, uVecs_(j + 1)));
  for (stormSize_t i = 0; i <= j; ++i) {
    stormBlas::Sub(rVecs_(i), rVecs_(i), uVecs_(i + 1), alpha_);
  }

  // ----------------------
  // Update the solution and the residual:
  // 𝒙 ← 𝒙 + 𝛼⋅𝒖₀,
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝒓ⱼ₊₁ ← 𝓟(𝒛 ← 𝓐𝒓ⱼ).
  // 𝗲𝗹𝘀𝗲:
  //   𝒓ⱼ₊₁ ← 𝓐𝒓ⱼ.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  stormBlas::Add(xVec, xVec, uVecs_(0), alpha_);
  if (preOp != nullptr) {
    stormBlas::MatVec(rVecs_(j + 1), *preOp, zVec_, linOp, rVecs_(j));
  } else {
    linOp.MatVec(rVecs_(j + 1), rVecs_(j));
  }

  if (j == l - 1) {
    // ----------------------
    // Minimal residual part:
    // 𝗳𝗼𝗿 𝑗 = 𝟣, 𝑙 𝗱𝗼:
    //   𝗳𝗼𝗿 𝑖 = 𝟣, 𝑗 - 𝟣 𝗱𝗼:
    //     𝜏ᵢⱼ ← <𝒓ᵢ⋅𝒓ⱼ>/𝜎ᵢ,
    //     𝒓ⱼ ← 𝒓ⱼ - 𝜏ᵢⱼ⋅𝒓ᵢ,
    //   𝗲𝗻𝗱 𝗳𝗼𝗿
    //   𝜎ⱼ ← <𝒓ⱼ⋅𝒓ⱼ>,
    //   𝛾̅ⱼ ← <𝒓₀⋅𝒓ⱼ>/𝜎ⱼ,
    // 𝗲𝗻𝗱 𝗳𝗼𝗿
    // ----------------------
    for (stormSize_t j = 1; j <= l; ++j) {
      for (stormSize_t i = 1; i < j; ++i) {
        tau_(i, j) = 
          stormUtils::SafeDivide(stormBlas::Dot(rVecs_(i), rVecs_(j)), sigma_(i));
        stormBlas::Sub(rVecs_(j), rVecs_(j), rVecs_(i), tau_(i, j));
      }
      sigma_(j) = stormBlas::Dot(rVecs_(j), rVecs_(j));
      gammaBar_(j) = 
        stormUtils::SafeDivide(stormBlas::Dot(rVecs_(0), rVecs_(j)), sigma_(j));
    }

    // ----------------------
    // 𝜔 ← 𝛾ₗ ← 𝛾̅ₗ, 𝜌 ← -𝜔⋅𝜌, 
    // 𝗳𝗼𝗿 𝑗 = 𝑙 - 𝟣, 𝟣, -𝟣 𝗱𝗼:
    //   𝛾ⱼ ← 𝛾̅ⱼ,
    //   𝗳𝗼𝗿 𝑖 = 𝑗 + 𝟣, 𝑙 𝗱𝗼:
    //     𝛾ⱼ ← 𝛾ⱼ - 𝜏ⱼᵢ⋅𝛾ᵢ,
    //   𝗲𝗻𝗱 𝗳𝗼𝗿
    // 𝗲𝗻𝗱 𝗳𝗼𝗿
    // 𝗳𝗼𝗿 𝑗 = 𝟣, 𝑙 - 𝟣 𝗱𝗼:
    //   𝛾̿ⱼ ← 𝛾ⱼ₊₁,
    //   𝗳𝗼𝗿 𝑖 = 𝑗 + 𝟣, 𝑙 - 𝟣 𝗱𝗼:
    //     𝛾̿ⱼ ← 𝛾̿ⱼ + 𝜏ⱼᵢ⋅𝛾ᵢ₊₁.
    //   𝗲𝗻𝗱 𝗳𝗼𝗿
    // 𝗲𝗻𝗱 𝗳𝗼𝗿
    // ----------------------
    omega_ = gamma_(l) = gammaBar_(l), rho_ *= -omega_;
    for (stormSize_t j = l - 1; j != 0; --j) {
      gamma_(j) = gammaBar_(j);
      for (stormSize_t i = j + 1; i <= l; ++i) {
        gamma_(j) -= tau_(j, i)*gamma_(i);
      }
    }
    for (stormSize_t j = 1; j < l; ++j) {
      gammaBarBar_(j) = gamma_(j + 1);
      for (stormSize_t i = j + 1; i < l; ++i) {
        gammaBarBar_(j) += tau_(j, i)*gamma_(i + 1);
      }
    }

    // ----------------------
    // Update the solution and the residual again:
    // 𝒙 ← 𝒙 + 𝛾₁⋅𝒓₀,
    // 𝒓₀ ← 𝒓₀ - 𝛾̅ₗ⋅𝒓ₗ,
    // 𝒖₀ ← 𝒖₀ - 𝛾ₗ⋅𝒖ₗ,
    // 𝗳𝗼𝗿 𝑗 = 𝟣, 𝑙 - 𝟣 𝗱𝗼:
    //   𝒙 ← 𝒙 + 𝛾̿ⱼ⋅𝒓ⱼ,
    //   𝒓₀ ← 𝒓₀ - 𝛾̅ⱼ⋅𝒓ⱼ,
    //   𝒖₀ ← 𝒖₀ - 𝛾ⱼ⋅𝒖ⱼ.
    // 𝗲𝗻𝗱 𝗳𝗼𝗿
    // ----------------------
    stormBlas::Add(xVec, xVec, rVecs_(0), gamma_(1));
    stormBlas::Sub(rVecs_(0), rVecs_(0), rVecs_(l), gammaBar_(l));
    stormBlas::Sub(uVecs_(0), uVecs_(0), uVecs_(l), gamma_(l));
    for (stormSize_t j = 1; j < l; ++j) {
      stormBlas::Add(xVec, xVec, rVecs_(j), gammaBarBar_(j));
      stormBlas::Sub(rVecs_(0), rVecs_(0), rVecs_(j), gammaBar_(j));
      stormBlas::Sub(uVecs_(0), uVecs_(0), uVecs_(j), gamma_(j));
    }
  }

  return stormBlas::Norm2(rVecs_(0));

} // stormBiCGStabLSolver<...>::InnerIterate

#endif // ifndef _STORM_SOLVER_BICGSTAB_
