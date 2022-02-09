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
#include <stormSolvers/Preconditioner.hxx>

namespace Storm {

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract polynomial preconditioner.
/// 
/// @todo Document me!
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Vector>
class PolynomialPreconditioner : public Preconditioner<Vector> {

}; // class PolynomialPreconditioner<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Chebyshev polynomial preconditioner.
/// 
/// @todo Document me!
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Vector>
class ChebyshevPreconditioner final : 
    public PolynomialPreconditioner<Vector> {
private:
  size_t NumIterations = 10;
  /// @todo: Estimate the true eigenvalue bounds!
  real_t lambdaMin = 0.3*8000.0, lambdaMax = 1.2*8000.0;
  mutable Vector pVec, rVec;
  Operator<Vector> const* linOp;

  void Build(Vector const& xVec,
             Vector const& bVec,
             Operator<Vector> const& linOp) override;

  void MatVec(Vector& yVec,
              Vector const& xVec) const override;

  void ConjMatVec(Vector& xVec,
                  Vector const& yVec) const override {
    MatVec(xVec, yVec);
  }

}; // class stormIdentityPreconditioner<...>

template<class Vector>
void ChebyshevPreconditioner<Vector>::Build(Vector const& xVec,
                                            Vector const& bVec,
                                            Operator<Vector> const& linOp) {

  pVec.Assign(xVec, false);
  rVec.Assign(xVec, false);
  this->linOp = &linOp;

} // ChebyshevPreconditioner<...>::Build

template<class Vector>
void ChebyshevPreconditioner<Vector>::MatVec(Vector& yVec,
                                             Vector const& xVec) const {

  assert(linOp != nullptr && "Preconditioner was not built!");

  // ----------------------
  // Initialize the Chebyshev iterations:
  // 𝒓 ← 𝒙,
  // 𝒚 ← {𝟢}ᵀ,
  // 𝑐 ← ½(𝜆ₘₐₓ - 𝜆ₘᵢₙ),
  // 𝑑 ← ½(𝜆ₘₐₓ + 𝜆ₘᵢₙ).
  // ----------------------
  Blas::Set(rVec, xVec);
  Blas::Fill(yVec, 0.0);
  real_t const c = 0.5*(lambdaMax - lambdaMin);
  real_t const d = 0.5*(lambdaMax + lambdaMin);

  real_t alpha;
  for (size_t iteration = 0; iteration < NumIterations; ++iteration) {
    
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
      Blas::Set(pVec, rVec);
    } else {
      real_t beta;
      if (iteration == 2) {
        beta = 0.5*std::pow(c*alpha, 2);
      } else {
        beta = std::pow(0.5*c*alpha, 2);
      }
      alpha /= (d*alpha - beta);
      Blas::Add(pVec, rVec, pVec, beta);
    }

    // ----------------------
    // Update the solution:
    // 𝒚 ← 𝒚 + 𝛼𝒑,
    // 𝒓 ← 𝓐𝒚,
    // 𝒓 ← 𝒙 - 𝒓.
    // ----------------------
    Blas::Add(yVec, yVec, pVec, alpha);
    linOp->MatVec(rVec, yVec);
    Blas::Sub(rVec, xVec, rVec);
  }

} // ChebyshevPreconditioner<...>::MatVec

} // namespace Storm
