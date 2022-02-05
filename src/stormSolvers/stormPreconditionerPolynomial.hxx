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
#ifndef _STORM_PRECONDITIONER_POLYNOMIAL_
#define _STORM_PRECONDITIONER_POLYNOMIAL_

#include <cmath>

#include <StormRuler_API.h>
#include <stormSolvers/stormPreconditioner.hxx>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract polynomial preconditioner.
/// 
/// @todo Document me!
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Vector>
class stormPolynomialPreconditioner : public stormPreconditioner<Vector> {
}; // class stormPolynomialPreconditioner<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Chebyshev polynomial preconditioner.
/// 
/// @todo Document me!
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Vector>
class stormChebyshevPreconditioner final : 
    public stormPolynomialPreconditioner<Vector> {
private:
  stormSize_t NumIterations = 10;
  /// @todo: Estimate the true eigenvalue bounds!
  stormReal_t lambdaMin = 0.3*8000.0, lambdaMax = 1.2*8000.0;
  mutable Vector pVec, rVec;
  stormOperator<Vector> const* linOp;

  void Build(Vector const& xVec,
             Vector const& bVec,
             stormOperator<Vector> const& linOp) override;

  void MatVec(Vector& yVec,
              Vector const& xVec) const override;

  void ConjMatVec(Vector& xVec,
                  Vector const& yVec) const override {
    MatVec(xVec, yVec);
  }

}; // class stormIdentityPreconditioner<...>

template<class Vector>
void stormChebyshevPreconditioner<Vector>::Build(Vector const& xVec,
                                                 Vector const& bVec,
                                                 stormOperator<Vector> const& linOp) {

  pVec.Assign(xVec, false);
  rVec.Assign(xVec, false);
  this->linOp = &linOp;

} // stormChebyshevPreconditioner<...>::Build

template<class Vector>
void stormChebyshevPreconditioner<Vector>::MatVec(Vector& yVec,
                                                  Vector const& xVec) const {

  assert(linOp != nullptr && "Preconditioner was not built!");

  // ----------------------
  // Initialize the Chebyshev iterations:
  // 𝒓 ← 𝒙,
  // 𝒚 ← {𝟢}ᵀ,
  // 𝑐 ← ½(𝜆ₘₐₓ - 𝜆ₘᵢₙ),
  // 𝑑 ← ½(𝜆ₘₐₓ + 𝜆ₘᵢₙ).
  // ----------------------
  rVec.Assign(xVec);
  yVec.Fill(0.0);
  const stormReal_t c = 0.5*(lambdaMax - lambdaMin);
  const stormReal_t d = 0.5*(lambdaMax + lambdaMin);

  stormReal_t alpha;
  for (stormSize_t iteration = 0; iteration < NumIterations; ++iteration) {
    
    // ----------------------
    // Continue the Chebyshev iterations:
    // 𝗶𝗳 𝑘 = 𝟢:
    //   𝛼 ← 1/𝑑,
    //   𝒑 ← 𝒓,
    // 𝗲𝗹𝘀𝗲:
    //   𝗶𝗳 𝑘 = 2: 𝛽 ← ½(𝑐⋅𝛼)²,
    //   𝗲𝗹𝘀𝗲: 𝛽 ← (½⋅𝑐⋅𝛼)², 𝗲𝗻𝗱 𝗶𝗳
    //   𝛼 ← 𝛼/(𝑑⋅𝛼 - 𝛽),
    //   𝒑 ← 𝒓 + 𝛽⋅𝒑.
    // 𝗲𝗻𝗱 𝗶𝗳
    // ----------------------
    if (iteration == 0) {
      alpha = 1.0/d;
      pVec.Assign(rVec);
    } else {
      stormReal_t beta;
      if (iteration == 2) {
        beta = 0.5*std::pow(c*alpha, 2);
      } else {
        beta = std::pow(0.5*c*alpha, 2);
      }
      alpha /= (d*alpha - beta);
      stormBlas::Add(pVec, rVec, pVec, beta);
    }

    // ----------------------
    // Update the solution:
    // 𝒚 ← 𝒚 + 𝛼𝒑,
    // 𝒓 ← 𝓐𝒚,
    // 𝒓 ← 𝒙 - 𝒓.
    // ----------------------
    stormBlas::Add(yVec, yVec, pVec, alpha);
    linOp->MatVec(rVec, yVec);
    stormBlas::Sub(rVec, xVec, rVec);
  }

} // stormChebyshevPreconditioner<...>::MatVec

#endif // ifndef _STORM_PRECONDITIONER_POLYNOMIAL_
