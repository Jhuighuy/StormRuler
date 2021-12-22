/// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ///
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
/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ///
#ifndef _STORM_PRECONDITIONER_POLYNOMIAL_
#define _STORM_PRECONDITIONER_POLYNOMIAL_

#include <StormRuler_API.h>
#include <stormOperator.hxx>

#include <cmath>

/// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ///
/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ///

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract polynomial preconditioner.
/// 
/// @todo Document me!
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class tArray>
class stormPolynomialPreconditioner : public stormPreconditioner<tArray> {
}; // class stormPolynomialPreconditioner<...>

/// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ///
/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ///

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Chebyshev polynomial preconditioner.
/// 
/// @todo Document me!
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class tArray>
class stormChebyshevPreconditioner final : public stormPolynomialPreconditioner<tArray> {
private:
  stormSize_t MaxIterations = 10;
  /// @todo: Estimate the true eigenvalue bounds!
  stormReal_t lambdaMin = 1.0, lambdaMax = 800.0;
  const stormOperator<tArray>* linOp;

  void Build(const stormOperator<tArray>& linOp) override final {
    this->linOp = &linOp;
  }

  void MatVec(tArray& yArr,
              const tArray& xArr) const override final;

  void ConjMatVec(tArray& xArr,
                  const tArray& yArr) const override final {
    MatVec(xArr, yArr);
  }

}; // class stormIdentityPreconditioner<...>

template<class tArray>
void stormChebyshevPreconditioner<tArray>::MatVec(tArray& yArr,
                                                  const tArray& xArr) const {
  assert(linOp != nullptr && "Preconditioner was not built!");
                                                    
  tArray pArr, rArr;
  stormUtils::AllocLike(xArr, pArr, rArr);

  // ----------------------
  // Initialize the Chebyshev iterations:
  // 𝒓 ← 𝒙,
  // 𝒚 ← {𝟢}ᵀ,
  // 𝑐 ← ½(𝜆ₘₐₓ - 𝜆ₘᵢₙ),
  // 𝑑 ← ½(𝜆ₘₐₓ + 𝜆ₘᵢₙ).
  // ----------------------
  stormUtils::Set(rArr, xArr);
  stormUtils::Fill(yArr, 0.0);
  const stormReal_t c = 0.5*(lambdaMax - lambdaMin);
  const stormReal_t d = 0.5*(lambdaMax + lambdaMin);

  for (stormSize_t iteration = 1; iteration <= MaxIterations; ++iteration) {
    
    // ----------------------
    // Continue the Chebyshev iterations:
    // 𝗶𝗳 𝑘 == 1:
    //   𝛼 ← 1/𝑑,
    //   𝒑 ← 𝒓,
    // 𝗲𝗹𝘀𝗲:
    //   𝗶𝗳 𝑘 == 2: 𝛽 ← ½(𝑐⋅𝛼)²,
    //   𝗲𝗹𝘀𝗲: 𝛽 ← (½⋅𝑐⋅𝛼)², 𝗲𝗻𝗱 𝗶𝗳
    //   𝛼 ← 1/(𝑑 - 𝛽/𝛼),
    //   𝒑 ← 𝒓 + 𝛽𝒑.
    // 𝗲𝗻𝗱 𝗶𝗳
    // ----------------------
    stormReal_t alpha;
    if (iteration == 1) {
      alpha = 1.0/d;
      stormUtils::Set(pArr, rArr);
    } else {
      stormReal_t beta;
      if (iteration == 2) {
        beta = 0.5*std::pow(c*alpha, 2);
      } else {
        beta = std::pow(0.5*c*alpha, 2);
      }
      alpha = 1.0/(d - beta/alpha);
      stormUtils::Add(pArr, rArr, pArr, beta);
    }

    // ----------------------
    // 𝒚 ← 𝒚 + 𝛼𝒑,
    // 𝒓 ← 𝓐𝒚,
    // 𝒓 ← 𝒙 - 𝒓.
    // ----------------------
    stormUtils::Add(yArr, yArr, pArr, alpha);
    linOp->MatVec(rArr, yArr);
    stormUtils::Sub(rArr, xArr, rArr);
  }

} // stormChebyshevPreconditioner<...>::MatVec

/// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ///
/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ///

#endif // ifndef _STORM_PRECONDITIONER_POLYNOMIAL_
