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

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Solve a linear operator equation with the
///   @c TFQMR (Transpose-Free Quasi-Minimal Residual) method.
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
template<class tArray>
class stormTfqmrSolver final : public stormIterativeSolver<tArray> {
private:
  stormReal_t tau, rho, theta, eta;
  tArray dArr, xTildeArr, rTildeArr, uArr, vArr, yArr, sArr, zArr;

  stormReal_t Init(tArray& xArr,
                   tArray const& bArr,
                   stormOperator<tArray> const& linOp,
                   stormPreconditioner<tArray> const* preOp) override;

  stormReal_t Iterate(tArray& xArr,
                      tArray const& bArr,
                      stormOperator<tArray> const& linOp,
                      stormPreconditioner<tArray> const* preOp) override;

  void Finalize(tArray& xArr,
                tArray const& bArr,
                stormOperator<tArray> const& linOp,
                stormPreconditioner<tArray> const* preOp) override;

}; // class stormTfqmrSolver<...>

template<class tArray>
stormReal_t stormTfqmrSolver<tArray>::Init(tArray& xArr,
                                           tArray const& bArr,
                                           stormOperator<tArray> const& linOp,
                                           stormPreconditioner<tArray> const* preOp) {

  bool const leftPre = (preOp != nullptr) && 
    (this->PreSide == stormPreconditionerSide::Left);
  bool const rightPre = (preOp != nullptr) && 
    (this->PreSide == stormPreconditionerSide::Right);

  stormUtils::AllocLike(xArr, dArr, rTildeArr, uArr, vArr, yArr, sArr);
  if (preOp != nullptr) {
    stormUtils::AllocLike(xArr, zArr);
    if (rightPre) stormUtils::AllocLike(xArr, xTildeArr);
  }

  // ----------------------
  // Initialize:
  // 𝒚 ← 𝓐𝒙,
  // 𝒚 ← 𝒃 - 𝒚,
  // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
  //   𝒛 ← 𝒚,
  //   𝒚 ← 𝓟𝒛,
  //   𝒔 ← 𝓟(𝒛 ← 𝓐𝒚),
  // 𝗲𝗹𝘀𝗲 𝗶𝗳 𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦:
  //   𝒙̃, 𝒙 ← 𝒙̃, 𝒙,
  //   𝒙̃ ← {𝟢}ᵀ,
  //   𝒔 ← 𝓐(𝒛 ← 𝓟𝒚),
  // 𝗲𝗹𝘀𝗲:
  //   𝒔 ← 𝓐𝒚.
  // 𝗲𝗻𝗱 𝗶𝗳
  // 𝒖 ← 𝒚,
  // 𝒗 ← 𝒔,
  // 𝒅 ← {𝟢}ᵀ,
  // 𝒓̃ ← 𝒖,
  // 𝜌 ← <𝒓̃⋅𝒖>,
  // 𝜏 ← 𝜌¹ᐟ², 𝜗 ← 𝟢, 𝜂 ← 𝟢.
  // ----------------------
  linOp.MatVec(yArr, xArr);
  stormBlas::Sub(yArr, bArr, yArr);
  if (leftPre) {
    std::swap(zArr, yArr);
    preOp->MatVec(yArr, zArr);
    stormBlas::MatVec(sArr, *preOp, zArr, linOp, yArr);
  } else if (rightPre) {
    std::swap(xArr, xTildeArr);
    stormBlas::Fill(xArr, 0.0);
    stormBlas::MatVec(sArr, linOp, zArr, *preOp, yArr);
  } else {
    linOp.MatVec(sArr, yArr);
  }
  stormBlas::Set(uArr, yArr);
  stormBlas::Set(vArr, sArr);
  stormBlas::Fill(dArr, 0.0);
  stormBlas::Set(rTildeArr, uArr);
  rho = stormBlas::Dot(rTildeArr, uArr);
  tau = std::sqrt(rho), theta = 0.0, eta = 0.0;

  return tau;

} // stormTfqmrSolver<...>::Init

template<class tArray>
stormReal_t stormTfqmrSolver<tArray>::Iterate(tArray& xArr,
                                              tArray const& bArr,
                                              stormOperator<tArray> const& linOp,
                                              stormPreconditioner<tArray> const* preOp) {

  bool const leftPre = (preOp != nullptr) && 
    (this->PreSide == stormPreconditionerSide::Left);
  bool const rightPre = (preOp != nullptr) && 
    (this->PreSide == stormPreconditionerSide::Right);

  // ----------------------
  // Continue the iterations:
  // 𝜎 ← <𝒓̃⋅𝒗>, 𝛼 ← 𝜌/𝜎,
  // 𝗳𝗼𝗿 𝑚 = 𝟢, 𝟣 𝗱𝗼:
  //   𝒖 ← 𝒖 - 𝛼⋅𝒔,
  //   𝒅 ← 𝒚 + (𝜗²⋅𝜂/𝛼)⋅𝒅,
  //   𝜗 ← ‖𝒖‖/𝜏,
  //   𝑐𝑠 ← 𝟣/(𝟣 + 𝜗²)¹ᐟ²,
  //   𝜏 ← 𝜏⋅𝜗⋅𝑐𝑠, 𝜂 ← 𝛼⋅(𝑐𝑠)²,
  //   𝒙 ← 𝒙 + 𝜂⋅𝒅,
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
  stormReal_t const sigma =
    stormBlas::Dot(rTildeArr, vArr), alpha = rho/sigma;
  for (stormSize_t m = 0; m <= 1; ++m) {
    stormBlas::Sub(uArr, uArr, sArr, alpha);
    stormBlas::Add(dArr, yArr, dArr, std::pow(theta, 2)*eta/alpha);
    theta = stormBlas::Norm2(uArr)/tau;
    stormReal_t const cs = 1.0/std::hypot(1.0, theta);
    tau *= theta*cs, eta = alpha*std::pow(cs, 2);
    stormBlas::Add(xArr, xArr, dArr, eta);
    if (m == 0) {
      stormBlas::Sub(yArr, yArr, vArr, alpha);
      if (leftPre) {
        stormBlas::MatVec(sArr, *preOp, zArr, linOp, yArr);
      } else if (rightPre) {
        stormBlas::MatVec(sArr, linOp, zArr, *preOp, yArr);
      } else {
        linOp.MatVec(sArr, yArr);
      }
    }
  }

  // ----------------------
  // 𝜌̅ ← 𝜌,
  // 𝜌 ← <𝒓̃⋅𝒖>, 𝛽 ← 𝜌/𝜌̅,
  // 𝒗 ← 𝒔 + 𝛽⋅𝒗,
  // 𝒚 ← 𝒖 + 𝛽⋅𝒚,
  // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
  //   𝒔 ← 𝓟(𝒛 ← 𝓐𝒚),
  // 𝗲𝗹𝘀𝗲 𝗶𝗳 𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦:
  //   𝒔 ← 𝓐(𝒛 ← 𝓟𝒚),
  // 𝗲𝗹𝘀𝗲:
  //   𝒔 ← 𝓐𝒚,
  // 𝗲𝗻𝗱 𝗶𝗳
  // 𝒗 ← 𝒔 + 𝛽⋅𝒗.
  // ----------------------
  stormReal_t const rhoBar = rho;
  rho = stormBlas::Dot(rTildeArr, uArr);
  stormReal_t const beta = rho/rhoBar;
  stormBlas::Add(vArr, sArr, vArr, beta);
  stormBlas::Add(yArr, uArr, yArr, beta);
  if (leftPre) {
    stormBlas::MatVec(sArr, *preOp, zArr, linOp, yArr);
  } else if (rightPre) {
    stormBlas::MatVec(sArr, linOp, zArr, *preOp, yArr);
  } else {
    linOp.MatVec(sArr, yArr);
  }
  stormBlas::Add(vArr, sArr, vArr, beta);

  // ----------------------
  // Compute the residual upper bound:
  // 𝜑̃ ← 𝜏⋅(𝟤𝑘 + 𝟥)¹ᐟ².
  // ----------------------
  stormSize_t const k = this->Iteration;
  stormReal_t const phiTilde = tau*std::sqrt(2.0*k + 3.0);

  return phiTilde;

} // stormTfqmrSolver<...>::Iterate

template<class tArray>
void stormTfqmrSolver<tArray>::Finalize(tArray& xArr,
                                        tArray const& bArr,
                                        stormOperator<tArray> const& linOp,
                                        stormPreconditioner<tArray> const* preOp) {

  bool const rightPre = (preOp != nullptr) && 
    (this->PreSide == stormPreconditionerSide::Right);

  // ----------------------
  // Finalize the solution:
  // 𝗶𝗳 𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦:
  //   𝒛 ← 𝓟𝒙,
  //   𝒙 ← 𝒙̃ + 𝒛.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  if (rightPre) {
    preOp->MatVec(zArr, xArr);
    stormBlas::Add(xArr, xTildeArr, zArr);
  }

} // stormTfqmrSolver<...>::Finalize

#endif // ifndef _STORM_SOLVER_TFQMR_
