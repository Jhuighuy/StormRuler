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

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief The @c CGS (Conjugate Gradients Squared)
///   linear operator equation solver.
///
/// @c CGS, like the other @c BiCG type solvers, requires
///   two operator multiplications per iteration.
///
/// @warning @c CGS convergence behavior may be very erratic.
///
/// References:
/// @verbatim
/// [1] Sonneveld, Peter.
///     “CGS, A Fast Lanczos-Type Solver for Nonsymmetric Linear systems.”
///     SIAM J. Sci. Stat. Comput., 10:36-52, 1989.
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<vector_like Vector>
class CgsSolver final : public IterativeSolver<Vector> {
private:

  real_t rho_;
  Vector pVec_, qVec_, rVec_, rTildeVec_, uVec_, vVec_;

  real_t Init(Vector const& xVec, Vector const& bVec,
              Operator<Vector> const& linOp,
              Preconditioner<Vector> const* preOp) override;

  real_t Iterate(Vector& xVec, Vector const& bVec,
                 Operator<Vector> const& linOp,
                 Preconditioner<Vector> const* preOp) override;

}; // class CgsSolver

template<vector_like Vector>
real_t CgsSolver<Vector>::Init(Vector const& xVec, Vector const& bVec,
                               Operator<Vector> const& linOp,
                               Preconditioner<Vector> const* preOp) {
  bool const leftPre{(preOp != nullptr) &&
                     (this->PreSide == PreconditionerSide::Left)};

  pVec_.Assign(xVec, false);
  qVec_.Assign(xVec, false);
  rVec_.Assign(xVec, false);
  rTildeVec_.Assign(xVec, false);
  uVec_.Assign(xVec, false);
  vVec_.Assign(xVec, false);

  // Initialize:
  // ----------------------
  // 𝒓 ← 𝒃 - 𝓐𝒙,
  // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
  //   𝒖 ← 𝒓,
  //   𝒓 ← 𝓟𝒖,
  // 𝗲𝗻𝗱 𝗶𝗳
  // 𝒓̃ ← 𝒓,
  // 𝜌 ← <𝒓̃⋅𝒓>.
  // ----------------------
  linOp.Residual(rVec_, bVec, xVec);
  if (leftPre) {
    uVec_.Swap(rVec_);
    preOp->MatVec(rVec_, uVec_);
  }
  rTildeVec_.Set(rVec_);
  rho_ = rTildeVec_.Dot(rVec_);

  return std::sqrt(rho_);

} // CgsSolver::Init

template<vector_like Vector>
real_t CgsSolver<Vector>::Iterate(Vector& xVec, Vector const& bVec,
                                  Operator<Vector> const& linOp,
                                  Preconditioner<Vector> const* preOp) {
  bool const leftPre{(preOp != nullptr) &&
                     (this->PreSide == PreconditionerSide::Left)};
  bool const rightPre{(preOp != nullptr) &&
                      (this->PreSide == PreconditionerSide::Right)};

  // Continue the iterations:
  // ----------------------
  // 𝗶𝗳 𝘍𝘪𝘳𝘴𝘵𝘐𝘵𝘦𝘳𝘢𝘵𝘪𝘰𝘯:
  //   𝒖 ← 𝒓,
  //   𝒑 ← 𝒖.
  // 𝗲𝗹𝘀𝗲:
  //   𝜌̅ ← 𝜌,
  //   𝜌 ← <𝒓̃⋅𝒓>,
  //   𝛽 ← 𝜌/𝜌̅,
  //   𝒖 ← 𝒓 + 𝛽⋅𝒒,
  //   𝒑 ← 𝒒 + 𝛽⋅𝒑,
  //   𝒑 ← 𝒖 + 𝛽⋅𝒑.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  bool const firstIteration{this->Iteration == 0};
  if (firstIteration) {
    uVec_.Set(rVec_);
    pVec_.Set(uVec_);
  } else {
    real_t const rhoBar{rho_};
    rho_ = rTildeVec_.Dot(rVec_);
    real_t const beta{Utils::SafeDivide(rho_, rhoBar)};
    uVec_.Add(rVec_, qVec_, beta);
    pVec_.Add(qVec_, pVec_, beta);
    pVec_.Add(uVec_, pVec_, beta);
  }

  // ----------------------
  // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
  //   𝒗 ← 𝓟(𝒒 ← 𝓐𝒑),
  // 𝗲𝗹𝘀𝗲 𝗶𝗳 𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦:
  //   𝒗 ← 𝓐(𝒒 ← 𝓟𝒑),
  // 𝗲𝗹𝘀𝗲:
  //   𝒗 ← 𝓐𝒑,
  // 𝗲𝗻𝗱 𝗶𝗳
  // 𝛼 ← 𝜌/<𝒓̃⋅𝒗>,
  // 𝒒 ← 𝒖 - 𝛼⋅𝒗,
  // 𝒗 ← 𝒖 + 𝒒.
  // ----------------------
  if (leftPre) {
    preOp->MatVec(vVec_, qVec_, linOp, pVec_);
  } else if (rightPre) {
    linOp.MatVec(vVec_, qVec_, *preOp, pVec_);
  } else {
    linOp.MatVec(vVec_, pVec_);
  }
  real_t const alpha{Utils::SafeDivide(rho_, rTildeVec_.Dot(vVec_))};
  qVec_.Sub(uVec_, vVec_, alpha);
  vVec_.Add(uVec_, qVec_);

  // Update the solution and the residual:
  // ----------------------
  // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
  //   𝒙 ← 𝒙 + 𝛼⋅𝒗,
  //   𝒗 ← 𝓟(𝒖 ← 𝓐𝒗),
  //   𝒓 ← 𝒓 - 𝛼⋅𝒗.
  // 𝗲𝗹𝘀𝗲 𝗶𝗳 𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦:
  //   𝒗 ← 𝓐(𝒖 ← 𝓟𝒗),
  //   𝒙 ← 𝒙 + 𝛼⋅𝒖,
  //   𝒓 ← 𝒓 - 𝛼⋅𝒗.
  // 𝗲𝗹𝘀𝗲:
  //   𝒖 ← 𝓐𝒗,
  //   𝒙 ← 𝒙 + 𝛼⋅𝒗,
  //   𝒓 ← 𝒓 - 𝛼⋅𝒖.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  if (leftPre) {
    xVec.AddAssign(vVec_, alpha);
    preOp->MatVec(vVec_, uVec_, linOp, vVec_);
    rVec_.SubAssign(vVec_, alpha);
  } else if (rightPre) {
    linOp.MatVec(vVec_, uVec_, *preOp, vVec_);
    xVec.AddAssign(uVec_, alpha);
    rVec_.SubAssign(vVec_, alpha);
  } else {
    linOp.MatVec(uVec_, vVec_);
    xVec.AddAssign(vVec_, alpha);
    rVec_.SubAssign(uVec_, alpha);
  }

  return rVec_.Norm2();

} // CgsSolver::Iterate

} // namespace Storm
