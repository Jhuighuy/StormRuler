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

#include <cmath>

#include <stormBase.hxx>
#include <stormSolvers/Solver.hxx>
#include <stormSolvers/Vector.hxx>

namespace Storm {

/// ----------------------------------------------------------------- ///
/// @brief Base class for @c TFQMR and @c TFQMR1.
/// ----------------------------------------------------------------- ///
template<vector_like Vector, bool L1>
class BaseTfqmrSolver : public IterativeSolver<Vector> {
private:

  real_t rho_, tau_;
  Vector dVec_, rTildeVec_, uVec_, vVec_, yVec_, sVec_, zVec_;

  real_t Init(Vector const& xVec, Vector const& bVec,
              Operator<Vector> const& linOp,
              Preconditioner<Vector> const* preOp) override;

  real_t Iterate(Vector& xVec, Vector const& bVec,
                 Operator<Vector> const& linOp,
                 Preconditioner<Vector> const* preOp) override;

protected:

  BaseTfqmrSolver() = default;

}; // class BaseTfqmrSolver

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief The @c TFQMR (Transpose-Free Quasi-Minimal Residual) \
///   linear operator equation solver.
///
/// @c TFQMR, like the other @c BiCG type methods, normally \
///   requires two operator-vector products per iteration. \
///   But, unlike the other @c BiCG type methods, @c TFQMR does not \
///   implicitly contain the residual norm estimate, only the rough \
///   upper bound is avariable, so at the latter iterations an extra \
///   operator-vector product per iteration may be required for the \
///   explicit residual estimation.
///
/// @c TFQMR typically converges much smoother, than \
///   @c CGS and @c BiCGStab. @todo Breakdowns?
///
/// References:
/// @verbatim
/// [1] Freund, Roland W.
///     “A Transpose-Free Quasi-Minimal Residual Algorithm
///      for Non-Hermitian Linear Systems.”
///     SIAM J. Sci. Comput. 14 (1993): 470-482.
/// [2] Freund, Roland W.
///     “Transpose-Free Quasi-Minimal Residual Methods
///      for Non-Hermitian Linear Systems.” (1994).
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<vector_like Vector>
class TfqmrSolver final : public BaseTfqmrSolver<Vector, false> {};

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief The @c TFQMR1 (Transpose-Free 1-norm \
///   Quasi-Minimal Residual) linear operator equation solver.
///
/// @c TFQMR1, like the other @c BiCG type solvers, requires \
///   two operator-vector products per iteration. Unlike @c TFQMR, \
///   @c TFQMR1 implicitly contains the residual norm estimate, so no \
///   extra operator-vector products are required.
///
/// @c TFQMR1 typically converges much smoother, than \
///   @c CGS and @c BiCGStab and is slightly faster than \
///   @c TFQMR. @todo Breakdowns?
///
/// References:
/// @verbatim
/// [1] H.M Bücker,
///     “A Transpose-Free 1-norm Quasi-Minimal Residual Algorithm
///      for Non-Hermitian Linear Systems.“, FZJ-ZAM-IB-9706.
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<vector_like Vector>
class Tfqmr1Solver final : public BaseTfqmrSolver<Vector, true> {};

template<vector_like Vector, bool L1>
real_t BaseTfqmrSolver<Vector, L1>::Init(Vector const& xVec, Vector const& bVec,
                                         Operator<Vector> const& linOp,
                                         Preconditioner<Vector> const* preOp) {
  bool const leftPre{(preOp != nullptr) &&
                     (this->PreSide == PreconditionerSide::Left)};

  dVec_.Assign(xVec, false);
  rTildeVec_.Assign(xVec, false);
  uVec_.Assign(xVec, false);
  vVec_.Assign(xVec, false);
  yVec_.Assign(xVec, false);
  sVec_.Assign(xVec, false);
  if (preOp != nullptr) { zVec_.Assign(xVec, false); }

  // Initialize:
  // ----------------------
  // 𝗶𝗳 𝘓₁:
  //   𝒅 ← 𝒙,
  // 𝗲𝗹𝘀𝗲:
  //   𝒅 ← {𝟢}ᵀ,
  // 𝗲𝗻𝗱 𝗶𝗳
  // 𝒚 ← 𝒃 - 𝓐𝒙,
  // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
  //   𝒛 ← 𝒚,
  //   𝒚 ← 𝓟𝒛,
  // 𝗲𝗻𝗱 𝗶𝗳
  // 𝒖 ← 𝒚,
  // 𝒓̃ ← 𝒖,
  // 𝜌 ← <𝒓̃⋅𝒓>, 𝜏 ← 𝜌¹ᐟ².
  // ----------------------
  if constexpr (L1) {
    dVec_.Set(xVec);
  } else {
    dVec_.Fill(0.0);
  }
  linOp.Residual(yVec_, bVec, xVec);
  if (leftPre) {
    zVec_.Swap(yVec_);
    preOp->MatVec(yVec_, zVec_);
  }
  uVec_.Set(yVec_);
  rTildeVec_.Set(uVec_);
  rho_ = rTildeVec_.Dot(uVec_), tau_ = std::sqrt(rho_);

  return tau_;

} // BaseTfqmrSolver::Init

template<vector_like Vector, bool L1>
real_t
BaseTfqmrSolver<Vector, L1>::Iterate(Vector& xVec, Vector const& bVec,
                                     Operator<Vector> const& linOp,
                                     Preconditioner<Vector> const* preOp) {
  bool const leftPre{(preOp != nullptr) &&
                     (this->PreSide == PreconditionerSide::Left)};
  bool const rightPre{(preOp != nullptr) &&
                      (this->PreSide == PreconditionerSide::Right)};

  // Continue the iterations:
  // ----------------------
  // 𝗶𝗳 𝘍𝘪𝘳𝘴𝘵𝘐𝘵𝘦𝘳𝘢𝘵𝘪𝘰𝘯:
  //   𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
  //     𝒔 ← 𝓟(𝒛 ← 𝓐𝒚),
  //   𝗲𝗹𝘀𝗲 𝗶𝗳 𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦:
  //     𝒔 ← 𝓐(𝒛 ← 𝓟𝒚),
  //   𝗲𝗹𝘀𝗲:
  //     𝒔 ← 𝓐𝒚.
  //   𝗲𝗻𝗱 𝗶𝗳
  //   𝒗 ← 𝒔,
  // 𝗲𝗹𝘀𝗲:
  //   𝜌̅ ← 𝜌,
  //   𝜌 ← <𝒓̃⋅𝒖>,
  //   𝛽 ← 𝜌/𝜌̅,
  //   𝒗 ← 𝒔 + 𝛽⋅𝒗,
  //   𝒚 ← 𝒖 + 𝛽⋅𝒚,
  //   𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
  //     𝒔 ← 𝓟(𝒛 ← 𝓐𝒚),
  //   𝗲𝗹𝘀𝗲 𝗶𝗳 𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦:
  //     𝒔 ← 𝓐(𝒛 ← 𝓟𝒚),
  //   𝗲𝗹𝘀𝗲:
  //     𝒔 ← 𝓐𝒚,
  //   𝗲𝗻𝗱 𝗶𝗳
  //   𝒗 ← 𝒔 + 𝛽⋅𝒗.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  bool const firstIteration{this->Iteration == 0};
  if (firstIteration) {
    if (leftPre) {
      preOp->MatVec(sVec_, zVec_, linOp, yVec_);
    } else if (rightPre) {
      linOp.MatVec(sVec_, zVec_, *preOp, yVec_);
    } else {
      linOp.MatVec(sVec_, yVec_);
    }
    vVec_.Set(sVec_);
  } else {
    real_t const rhoBar{rho_};
    rho_ = rTildeVec_.Dot(uVec_);
    real_t const beta{Utils::SafeDivide(rho_, rhoBar)};
    vVec_.Add(sVec_, vVec_, beta);
    yVec_.Add(uVec_, yVec_, beta);
    if (leftPre) {
      preOp->MatVec(sVec_, zVec_, linOp, yVec_);
    } else if (rightPre) {
      linOp.MatVec(sVec_, zVec_, *preOp, yVec_);
    } else {
      linOp.MatVec(sVec_, yVec_);
    }
    vVec_.Add(sVec_, vVec_, beta);
  }

  // Update the solution:
  // ----------------------
  // 𝛼 ← 𝜌/<𝒓̃⋅𝒗>,
  // 𝗳𝗼𝗿 𝑚 = 𝟢, 𝟣 𝗱𝗼:
  //   𝒖 ← 𝒖 - 𝛼⋅𝒔,
  //   𝒅 ← 𝒅 + 𝛼⋅(𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦 ? 𝒛 : 𝒚),
  //   𝜔 ← ‖𝒖‖,
  //   𝗶𝗳 𝘓₁:
  //     𝗶𝗳 𝜔 < 𝜏:
  //       𝜏 ← 𝜔, 𝒙 ← 𝒅,
  //     𝗲𝗻𝗱 𝗶𝗳
  //   𝗲𝗹𝘀𝗲:
  //     𝑐𝑠, 𝑠𝑛 ← 𝘚𝘺𝘮𝘖𝘳𝘵𝘩𝘰(𝜏, 𝜔),
  //     𝜏 ← 𝑐𝑠⋅𝜔,
  //     𝒙 ← 𝒙 + 𝑐𝑠²⋅𝒅,
  //     𝒅 ← 𝑠𝑛²⋅𝒅,
  //   𝗲𝗻𝗱 𝗶𝗳
  //   𝗶𝗳 𝑚 = 𝟢:
  //     𝒚 ← 𝒚 - 𝛼⋅𝒗,
  //     𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
  //       𝒔 ← 𝓟(𝒛 ← 𝓐𝒚).
  //     𝗲𝗹𝘀𝗲 𝗶𝗳 𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦:
  //       𝒔 ← 𝓐(𝒛 ← 𝓟𝒚).
  //     𝗲𝗹𝘀𝗲:
  //       𝒔 ← 𝓐𝒚.
  //     𝗲𝗻𝗱 𝗶𝗳
  //   𝗲𝗻𝗱 𝗶𝗳
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // ----------------------
  real_t const alpha{Utils::SafeDivide(rho_, rTildeVec_.Dot(vVec_))};
  for (size_t m{0}; m <= 1; ++m) {
    uVec_.SubAssign(sVec_, alpha);
    dVec_.AddAssign(rightPre ? zVec_ : yVec_, alpha);
    real_t const omega{uVec_.Norm2()};
    if constexpr (L1) {
      if (omega < tau_) { tau_ = omega, xVec.Set(dVec_); }
    } else {
      auto const [cs, sn, rr] = Utils::SymOrtho(tau_, omega);
      tau_ = omega * cs;
      xVec.AddAssign(dVec_, std::pow(cs, 2));
      dVec_.ScaleAssign(std::pow(sn, 2));
    }
    if (m == 0) {
      yVec_.SubAssign(vVec_, alpha);
      if (leftPre) {
        preOp->MatVec(sVec_, zVec_, linOp, yVec_);
      } else if (rightPre) {
        linOp.MatVec(sVec_, zVec_, *preOp, yVec_);
      } else {
        linOp.MatVec(sVec_, yVec_);
      }
    }
  }

  // Compute the residual norm
  // (or it's upper bound estimate in the ℒ₂ case):
  // ----------------------
  // 𝜏̃ ← 𝜏,
  // 𝗶𝗳 𝗻𝗼𝘁 𝘓₁:
  //   𝜏̃ ← 𝜏⋅(𝟤𝑘 + 𝟥)¹ᐟ².
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  real_t tauTilde{tau_};
  if constexpr (!L1) {
    size_t const k = this->Iteration;
    tauTilde *= std::sqrt(2.0 * k + 3.0);
  }

  return tauTilde;

} // BaseTfqmrSolver::Iterate

} // namespace Storm
