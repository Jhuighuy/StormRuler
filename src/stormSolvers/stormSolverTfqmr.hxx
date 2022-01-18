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
#ifndef _STORM_SOLVER_TFQMR_
#define _STORM_SOLVER_TFQMR_

#include <stormSolvers/stormSolver.hxx>

/// ----------------------------------------------------------------- ///
/// @brief Base class for @c TFQMR and @c TFQMR1.
/// ----------------------------------------------------------------- ///
template<bool L1, class Vector>
class stormBaseTfqmrSolver : public stormIterativeSolver<Vector> {
private:
  stormReal_t rho_, tau_;
  Vector dVec_, rTildeVec_, uVec_, vVec_, yVec_, sVec_, zVec_;

  stormReal_t Init(Vector& xVec,
                   Vector const& bVec,
                   stormOperator<Vector> const& linOp,
                   stormPreconditioner<Vector> const* preOp) override;

  stormReal_t Iterate(Vector& xVec,
                      Vector const& bVec,
                      stormOperator<Vector> const& linOp,
                      stormPreconditioner<Vector> const* preOp) override;

protected:

  stormBaseTfqmrSolver() = default;

}; // class stormBaseTfqmrSolver<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Solve a linear operator equation with the
///   @c TFQMR (Transpose-Free Quasi-Minimal Residual) method.
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
template<class Vector>
class stormTfqmrSolver final : public stormBaseTfqmrSolver<false, Vector> {

}; // class stormTfqmrSolver<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Solve a linear operator equation with the
///   @c TFQMR1 (Transpose-Free 1-norm Quasi-Minimal Residual) method.
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
template<class Vector>
class stormTfqmr1Solver final : public stormBaseTfqmrSolver<true, Vector> {

}; // class stormTfqmr1Solver<...>

template<bool L1, class Vector>
stormReal_t stormBaseTfqmrSolver<L1, Vector>::
                              Init(Vector& xVec,
                                   Vector const& bVec,
                                   stormOperator<Vector> const& linOp,
                                   stormPreconditioner<Vector> const* preOp) {

  bool const leftPre = (preOp != nullptr) && 
    (this->PreSide == stormPreconditionerSide::Left);

  dVec_.Assign(xVec, false);
  rTildeVec_.Assign(xVec, false);
  uVec_.Assign(xVec, false);
  vVec_.Assign(xVec, false);
  yVec_.Assign(xVec, false);
  sVec_.Assign(xVec, false);
  if (preOp != nullptr) {
    zVec_.Assign(xVec, false);
  }

  // ----------------------
  // Initialize:
  // 𝗶𝗳 𝘓₁:
  //   𝒅 ← 𝒙,
  // 𝗲𝗹𝘀𝗲:
  //   𝒅 ← {𝟢}ᵀ,
  // 𝗲𝗻𝗱 𝗶𝗳
  // 𝒚 ← 𝓐𝒙,
  // 𝒚 ← 𝒃 - 𝒚,
  // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
  //   𝒛 ← 𝒚,
  //   𝒚 ← 𝓟𝒛,
  // 𝗲𝗻𝗱 𝗶𝗳
  // 𝒖 ← 𝒚,
  // 𝒓̃ ← 𝒖,
  // 𝜌 ← <𝒓̃⋅𝒓>, 𝜏 ← 𝜌¹ᐟ².
  // ----------------------
  if constexpr (L1) {
    stormBlas::Set(dVec_, xVec);
  } else {
    stormBlas::Fill(dVec_, 0.0);
  }
  linOp.MatVec(yVec_, xVec);
  stormBlas::Sub(yVec_, bVec, yVec_);
  if (leftPre) {
    std::swap(zVec_, yVec_);
    preOp->MatVec(yVec_, zVec_);
  }
  stormBlas::Set(uVec_, yVec_);
  stormBlas::Set(rTildeVec_, uVec_);
  rho_ = stormBlas::Dot(rTildeVec_, uVec_), tau_ = std::sqrt(rho_);

  return tau_;

} // stormBaseTfqmrSolver<...>::Init

template<bool L1, class Vector>
stormReal_t stormBaseTfqmrSolver<L1, Vector>::
                           Iterate(Vector& xVec,
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
  bool const firstIteration = this->Iteration == 0;
  if (firstIteration) {
    if (leftPre) {
      stormBlas::MatVec(sVec_, *preOp, zVec_, linOp, yVec_);
    } else if (rightPre) {
      stormBlas::MatVec(sVec_, linOp, zVec_, *preOp, yVec_);
    } else {
      linOp.MatVec(sVec_, yVec_);
    }
    stormBlas::Set(vVec_, sVec_);
  } else {
    stormReal_t const rhoBar = rho_;
    rho_ = stormBlas::Dot(rTildeVec_, uVec_);
    stormReal_t const beta = rho_/rhoBar;
    stormBlas::Add(vVec_, sVec_, vVec_, beta);
    stormBlas::Add(yVec_, uVec_, yVec_, beta);
    if (leftPre) {
      stormBlas::MatVec(sVec_, *preOp, zVec_, linOp, yVec_);
    } else if (rightPre) {
      stormBlas::MatVec(sVec_, linOp, zVec_, *preOp, yVec_);
    } else {
      linOp.MatVec(sVec_, yVec_);
    }
    stormBlas::Add(vVec_, sVec_, vVec_, beta);
  }

  // ----------------------
  // Update the solution:
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
  stormReal_t const alpha =
    stormUtils::SafeDivide(rho_, stormBlas::Dot(rTildeVec_, vVec_));
  for (stormSize_t m = 0; m <= 1; ++m) {
    stormBlas::Sub(uVec_, uVec_, sVec_, alpha);
    stormBlas::Add(dVec_, dVec_, rightPre ? zVec_ : yVec_, alpha);
    stormReal_t const omega = stormBlas::Norm2(uVec_);
    if constexpr (L1) {
      if (omega < tau_) {
        tau_ = omega, stormBlas::Set(xVec, dVec_);
      }
    } else {
      auto const [cs, sn, rr] =
        stormBlas::SymOrtho(tau_, omega);
      tau_ = omega*cs;
      stormBlas::Add(xVec, xVec, dVec_, std::pow(cs, 2));
      stormBlas::Scale(dVec_, dVec_, std::pow(sn, 2));
    }
    if (m == 0) {
      stormBlas::Sub(yVec_, yVec_, vVec_, alpha);
      if (leftPre) {
        stormBlas::MatVec(sVec_, *preOp, zVec_, linOp, yVec_);
      } else if (rightPre) {
        stormBlas::MatVec(sVec_, linOp, zVec_, *preOp, yVec_);
      } else {
        linOp.MatVec(sVec_, yVec_);
      }
    }
  }

  // ----------------------
  // Compute the residual norm 
  // (or it's upper bound estimate in the ℒ₂ case):
  // 𝜏̃ ← 𝜏,
  // 𝗶𝗳 𝗻𝗼𝘁 𝘓₁:
  //   𝜏̃ ← 𝜏⋅(𝟤𝑘 + 𝟥)¹ᐟ².
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  stormReal_t tauTilde = tau_;
  if constexpr (!L1) {
    stormSize_t const k = this->Iteration;
    tauTilde *= std::sqrt(2.0*k + 3.0);
  }

  return tauTilde;

} // stormBaseTfqmrSolver<...>::Iterate

#endif // ifndef _STORM_SOLVER_TFQMR_
