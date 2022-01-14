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
/// [1] van der Vorst, Henk A.
///     “Bi-CGSTAB: A Fast and Smoothly Converging Variant of Bi-CG
///      for the Solution of Nonsymmetric Linear Systems.”
///     SIAM J. Sci. Comput. 13 (1992): 631-644.
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Vector>
class stormBiCgStabSolver final : public stormIterativeSolver<Vector> {
private:
  stormReal_t alpha, rho, omega;
  Vector pVec, rVec, rTildeVec, tVec, vVec, zVec;

  stormReal_t Init(Vector& xVec,
                   Vector const& bVec,
                   stormOperator<Vector> const& linOp,
                   stormPreconditioner<Vector> const* preOp) override;

  stormReal_t Iterate(Vector& xVec,
                      Vector const& bVec,
                      stormOperator<Vector> const& linOp,
                      stormPreconditioner<Vector> const* preOp) override;

}; // class stormBiCgStabSolver<...>

template<class Vector>
stormReal_t stormBiCgStabSolver<Vector>::Init(Vector& xVec,
                                              Vector const& bVec,
                                              stormOperator<Vector> const& linOp,
                                              stormPreconditioner<Vector> const* preOp) {

  bool const leftPre = (preOp != nullptr) && 
    (this->PreSide == stormPreconditionerSide::Left);

  stormUtils::AllocLike(xVec, pVec, rVec, rTildeVec, tVec, vVec);
  if (preOp != nullptr) {
    stormUtils::AllocLike(xVec, zVec);
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
  linOp.MatVec(rVec, xVec);
  stormBlas::Sub(rVec, bVec, rVec);
  if (leftPre) {
    std::swap(zVec, rVec);
    preOp->MatVec(rVec, zVec);
  }
  stormBlas::Set(rTildeVec, rVec);
  rho = stormBlas::Dot(rTildeVec, rVec);

  return std::sqrt(rho);

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
    stormBlas::Set(pVec, rVec);
  } else {
    stormReal_t const rhoBar = rho; 
    rho = stormBlas::Dot(rTildeVec, rVec);
    stormReal_t const beta = 
      stormUtils::SafeDivide(rho, rhoBar)*stormUtils::SafeDivide(alpha, omega);
    stormBlas::Sub(pVec, pVec, vVec, omega);
    stormBlas::Add(pVec, rVec, pVec, beta);
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
    stormBlas::MatVec(vVec, *preOp, zVec, linOp, pVec);
  } else if (rightPre) {
    stormBlas::MatVec(vVec, linOp, zVec, *preOp, pVec);
  } else {
    linOp.MatVec(vVec, pVec);
  }
  alpha = stormUtils::SafeDivide(rho, stormBlas::Dot(rTildeVec, vVec));
  stormBlas::Add(xVec, xVec, rightPre ? zVec : pVec, alpha);
  stormBlas::Sub(rVec, rVec, vVec, alpha);

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
    stormBlas::MatVec(tVec, *preOp, zVec, linOp, rVec);
  } else if (rightPre) {
    stormBlas::MatVec(tVec, linOp, zVec, *preOp, rVec);
  } else {
    linOp.MatVec(tVec, rVec);
  }
  omega = stormUtils::SafeDivide(
    stormBlas::Dot(tVec, rVec), stormBlas::Dot(tVec, tVec));
  stormBlas::Add(xVec, xVec, rightPre ? zVec : rVec, omega);
  stormBlas::Sub(rVec, rVec, tVec, omega);

  return stormBlas::Norm2(rVec);

} // stormBiCgStabSolver<...>::Iterate

#endif // ifndef _STORM_SOLVER_BICGSTAB_
