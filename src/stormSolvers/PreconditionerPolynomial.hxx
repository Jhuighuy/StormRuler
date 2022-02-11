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
public:
  size_t Degree = 10;

private:
  real_t theta_, delta_;
  mutable Vector rVec_, pVec_;
  Operator<Vector> const* LinOp_;

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

  rVec_.Assign(xVec, false);
  pVec_.Assign(xVec, false);
  this->LinOp_ = &linOp;

  //PowerIterations<Vector> powerIterations;
  //lambdaMax =
  //  powerIterations.EstimateLargestEigenvalue(pVec, linOp);
  //lambdaMin = 0.01*lambdaMax;

  /// @todo: Estimate the true eigenvalue bounds!
  real_t const alpha = 0.95*1.046599390654509e+00;
  real_t const beta = 1.05*8.003575342439456e+02;

  theta_ = 0.5*(beta + alpha);
  delta_ = 0.5*(beta - alpha);

} // ChebyshevPreconditioner<...>::Build

template<class Vector>
void ChebyshevPreconditioner<Vector>::MatVec(Vector& yVec,
                                             Vector const& xVec) const {

  // ----------------------
  // Clear the solution:
  // 𝒚 ← {𝟢}ᵀ.
  // ----------------------
  Blas::Fill(yVec, 0.0);

  real_t alpha;
  for (size_t k = 0; k < Degree; ++k) {

    // ----------------------
    // Compute the residual:
    // 𝒓 ← 𝒙 - 𝓐𝒚.
    // ----------------------
    LinOp_->Residual(rVec_, xVec, yVec);

    // ----------------------
    // Update the solution:
    // 𝗶𝗳 𝑘 = 𝟢:
    //   𝒑 ← 𝒓/𝜃,
    // 𝗲𝗹𝘀𝗲:
    //   𝛼 ← 𝑘 = 𝟣 ? 𝟤⋅𝜃/(𝟤⋅𝜃² - 𝛿²) : 𝟣/(𝜃 - ¼⋅𝛼⋅𝛿²),
    //   𝛽 ← 𝛼⋅𝜃 - 𝟣,
    //   𝒑 ← 𝛼⋅𝒓 + 𝛽⋅𝒑,
    // 𝗲𝗻𝗱 𝗶𝗳
    // 𝒚 ← 𝒚 + 𝒑.
    // ----------------------
    if (k == 0) {
      Blas::Scale(pVec_, rVec_, 1.0/theta_);
    } else {
      alpha = k == 1 ?
        2.0*theta_/(2.0*std::pow(theta_, 2) - std::pow(delta_, 2)) :
        1.0/(theta_ - 0.25*alpha*std::pow(delta_, 2));
      real_t const beta = alpha*theta_ - 1.0;
      Blas::Add(pVec_, rVec_, alpha, pVec_, beta);
    }
    Blas::Add(yVec, yVec, pVec_);

  }

} // ChebyshevPreconditioner<...>::MatVec

} // namespace Storm
