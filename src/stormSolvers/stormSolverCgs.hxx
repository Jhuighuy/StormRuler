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
#ifndef _STORM_SOLVER_CGS_HXX_
#define _STORM_SOLVER_CGS_HXX_

#include <cmath>

#include <stormBase.hxx>
#include <stormSolvers/stormSolver.hxx>

_STORM_NAMESPACE_BEGIN_

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Solve a non-singular operator equation \
///   equation with the @c CGS (Conjugate Gradients Squared) method.
///
/// @c CGS, like the other @c BiCG type solvers, requires \
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
template<class Vector>
class CgsSolver final : public IterativeSolver<Vector> {
private:
  Real_t rho_;
  Vector pVec_, qVec_, rVec_, rTildeVec_, uVec_, vVec_;

  Real_t Init(Vector const& xVec,
              Vector const& bVec,
              Operator<Vector> const& linOp,
              Preconditioner<Vector> const* preOp) override;

  Real_t Iterate(Vector& xVec,
                 Vector const& bVec,
                 Operator<Vector> const& linOp,
                 Preconditioner<Vector> const* preOp) override;

}; // class CgsSolver<...>

template<class Vector>
Real_t CgsSolver<Vector>::Init(Vector const& xVec,
                               Vector const& bVec,
                               Operator<Vector> const& linOp,
                               Preconditioner<Vector> const* preOp) {

  bool const leftPre = (preOp != nullptr) && 
    (this->PreSide == PreconditionerSide::Left);

  pVec_.Assign(xVec, false);
  qVec_.Assign(xVec, false); 
  rVec_.Assign(xVec, false); 
  rTildeVec_.Assign(xVec, false);
  uVec_.Assign(xVec, false); 
  vVec_.Assign(xVec, false);

  // ----------------------
  // 𝒓 ← 𝓐𝒙,
  // 𝒓 ← 𝒃 - 𝒓,
  // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
  //   𝒖 ← 𝒓,
  //   𝒓 ← 𝓟𝒖,
  // 𝗲𝗻𝗱 𝗶𝗳
  // 𝒓̃ ← 𝒓,
  // 𝜌 ← <𝒓̃⋅𝒓>.
  // ----------------------
  linOp.MatVec(rVec_, xVec);
  Blas::Sub(rVec_, bVec, rVec_);
  if (leftPre) {
    std::swap(uVec_, rVec_);
    preOp->MatVec(rVec_, uVec_);
  }
  Blas::Set(rTildeVec_, rVec_);
  rho_ = Blas::Dot(rTildeVec_, rVec_);

  return std::sqrt(rho_);

} // CgsSolver<...>::Init

template<class Vector>
Real_t CgsSolver<Vector>::Iterate(Vector& xVec,
                                  Vector const& bVec,
                                  Operator<Vector> const& linOp,
                                  Preconditioner<Vector> const* preOp) {

  bool const leftPre = (preOp != nullptr) && 
    (this->PreSide == PreconditionerSide::Left);
  bool const rightPre = (preOp != nullptr) && 
    (this->PreSide == PreconditionerSide::Right);

  // ----------------------
  // Continue the iterations:
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
  bool const firstIteration = this->Iteration == 0;
  if (firstIteration) {
    Blas::Set(uVec_, rVec_);
    Blas::Set(pVec_, uVec_);
  } else {
    Real_t const rhoBar = rho_; 
    rho_ = Blas::Dot(rTildeVec_, rVec_);
    Real_t const beta = Utils::SafeDivide(rho_, rhoBar);
    Blas::Add(uVec_, rVec_, qVec_, beta);
    Blas::Add(pVec_, qVec_, pVec_, beta);
    Blas::Add(pVec_, uVec_, pVec_, beta);
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
  Real_t const alpha = 
    Utils::SafeDivide(rho_, Blas::Dot(rTildeVec_, vVec_));
  Blas::Sub(qVec_, uVec_, vVec_, alpha);
  Blas::Add(vVec_, uVec_, qVec_);

  // ----------------------
  // Update the solution and the residual:
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
    Blas::Add(xVec, xVec, vVec_, alpha);
    preOp->MatVec(vVec_, uVec_, linOp, vVec_);
    Blas::Sub(rVec_, rVec_, vVec_, alpha);
  } else if (rightPre) {
    linOp.MatVec(vVec_, uVec_, *preOp, vVec_);
    Blas::Add(xVec, xVec, uVec_, alpha);
    Blas::Sub(rVec_, rVec_, vVec_, alpha);
  } else {
    linOp.MatVec(uVec_, vVec_);
    Blas::Add(xVec, xVec, vVec_, alpha);
    Blas::Sub(rVec_, rVec_, uVec_, alpha);
  }

  return Blas::Norm2(rVec_);

} // CgsSolver<...>::Iterate

_STORM_NAMESPACE_END_

#endif // ifndef _STORM_SOLVER_CGS_HXX_
