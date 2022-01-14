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

#include <stormSolvers/stormSolver.hxx>

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
class stormCgsSolver final : public stormIterativeSolver<Vector> {
private:
  stormReal_t rho;
  Vector pVec, qVec, rVec, rTildeVec, uVec, vVec;

  stormReal_t Init(Vector& xVec,
                   Vector const& bVec,
                   stormOperator<Vector> const& linOp,
                   stormPreconditioner<Vector> const* preOp) override;

  stormReal_t Iterate(Vector& xVec,
                      Vector const& bVec,
                      stormOperator<Vector> const& linOp,
                      stormPreconditioner<Vector> const* preOp) override;

}; // class stormCgsSolver<...>

template<class Vector>
stormReal_t stormCgsSolver<Vector>::Init(Vector& xVec,
                                         Vector const& bVec,
                                         stormOperator<Vector> const& linOp,
                                         stormPreconditioner<Vector> const* preOp) {

  bool const leftPre = (preOp != nullptr) && 
    (this->PreSide == stormPreconditionerSide::Left);

  stormUtils::AllocLike(xVec, pVec, qVec, rVec, rTildeVec, uVec, vVec);

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
  linOp.MatVec(rVec, xVec);
  stormBlas::Sub(rVec, bVec, rVec);
  if (leftPre) {
    std::swap(uVec, rVec);
    preOp->MatVec(rVec, uVec);
  }
  stormBlas::Set(rTildeVec, rVec);
  rho = stormBlas::Dot(rTildeVec, rVec);

  return std::sqrt(rho);

} // stormCgsSolver<...>::Init

template<class Vector>
stormReal_t stormCgsSolver<Vector>::Iterate(Vector& xVec,
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
    stormBlas::Set(uVec, rVec);
    stormBlas::Set(pVec, uVec);
  } else {
    stormReal_t const rhoBar = rho; 
    rho = stormBlas::Dot(rTildeVec, rVec);
    stormReal_t const beta = 
      stormUtils::SafeDivide(rho, rhoBar);
    stormBlas::Add(uVec, rVec, qVec, beta);
    stormBlas::Add(pVec, qVec, pVec, beta);
    stormBlas::Add(pVec, uVec, pVec, beta);
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
    stormBlas::MatVec(vVec, *preOp, qVec, linOp, pVec);
  } else if (rightPre) {
    stormBlas::MatVec(vVec, linOp, qVec, *preOp, pVec);
  } else {
    linOp.MatVec(vVec, pVec);
  }
  stormReal_t const alpha = 
    stormUtils::SafeDivide(rho, stormBlas::Dot(rTildeVec, vVec));
  stormBlas::Sub(qVec, uVec, vVec, alpha);
  stormBlas::Add(vVec, uVec, qVec);

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
    stormBlas::Add(xVec, xVec, vVec, alpha);
    stormBlas::MatVec(vVec, *preOp, uVec, linOp, vVec);
    stormBlas::Sub(rVec, rVec, vVec, alpha);
  } else if (rightPre) {
    stormBlas::MatVec(vVec, linOp, uVec, *preOp, vVec);
    stormBlas::Add(xVec, xVec, uVec, alpha);
    stormBlas::Sub(rVec, rVec, vVec, alpha);
  } else {
    linOp.MatVec(uVec, vVec);
    stormBlas::Add(xVec, xVec, vVec, alpha);
    stormBlas::Sub(rVec, rVec, uVec, alpha);
  }

  return stormBlas::Norm2(rVec);

} // stormCgsSolver<...>::Iterate

#endif // ifndef _STORM_SOLVER_CGS_HXX_
