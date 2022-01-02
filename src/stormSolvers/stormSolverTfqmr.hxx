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
  tArray dArr, rTildeArr, uHatArr, vArr, yArr, yBarArr, zArr, zBarArr;

  stormReal_t Init(tArray& xArr,
                   tArray const& bArr,
                   stormOperator<tArray> const& linOp,
                   stormPreconditioner<tArray> const* preOp) override;

  stormReal_t Iterate(tArray& xArr,
                      tArray const& bArr,
                      stormOperator<tArray> const& linOp,
                      stormPreconditioner<tArray> const* preOp) override;

}; // class stormTfqmrSolver<...>

template<class tArray>
stormReal_t stormTfqmrSolver<tArray>::Init(tArray& xArr,
                                           tArray const& bArr,
                                           stormOperator<tArray> const& linOp,
                                           stormPreconditioner<tArray> const* preOp) {
  // ----------------------
  // Allocate the intermediate arrays:
  // ----------------------
  stormUtils::AllocLike(xArr, dArr, rTildeArr, uHatArr, vArr, 
    yArr, yBarArr, zArr, zBarArr);
  if (preOp != nullptr) {
    //stormUtils::AllocLike(xArr, wArr, yArr, zArr);
  }

  // ----------------------
  // Initialize:
  // 𝒚 ← 𝓐𝒙,
  // 𝒚 ← 𝒃 - 𝒚,
  // 𝒖̂ ← 𝒚,
  // 𝒛 ← 𝓐𝒚, 
  // 𝒗 ← 𝒛,
  // 𝒅 ← {𝟢}ᵀ,
  // 𝒓̃ ← 𝒖̂,
  // 𝜌 ← <𝒓̃⋅𝒖̂>,
  // 𝜏 ← 𝜌¹ᐟ², 𝜗 ← 𝟢, 𝜂 ← 𝟢.
  // ----------------------
  linOp.MatVec(yArr, xArr);
  stormBlas::Sub(yArr, bArr, yArr);
  stormBlas::Set(uHatArr, yArr);
  linOp.MatVec(zArr, yArr);
  stormBlas::Set(vArr, zArr);
  stormBlas::Fill(dArr, 0.0);
  stormBlas::Set(rTildeArr, uHatArr);
  rho = stormBlas::Dot(rTildeArr, uHatArr);
  tau = std::sqrt(rho), theta = 0.0, eta = 0.0;

  return tau;

} // stormTfqmrSolver<...>::Init

template<class tArray>
stormReal_t stormTfqmrSolver<tArray>::Iterate(tArray& xArr,
                                              tArray const& bArr,
                                              stormOperator<tArray> const& linOp,
                                              stormPreconditioner<tArray> const* preOp) {

  // ----------------------
  // Continue the iterations:
  // 𝜎 ← <𝒓̃⋅𝒗>, 𝛼 ← 𝜌/𝜎,
  // 𝒚̅, 𝒛̅ ← 𝒚, 𝒛,
  // 𝒚 ← 𝒚̅ - 𝛼⋅𝒗,
  // 𝒛 ← 𝓐𝒚.
  // ----------------------
  stormReal_t const sigma = 
    stormBlas::Dot(rTildeArr, vArr), alpha = rho/sigma;
  std::swap(yBarArr, yArr), std::swap(zBarArr, zArr);
  stormBlas::Sub(yArr, yBarArr, vArr, alpha);
  linOp.MatVec(zArr, yArr);

  // ----------------------
  // 𝗳𝗼𝗿 𝑚 = 𝟢, 𝟣 𝗱𝗼:
  //   𝒖̂ ← 𝒖̂ - 𝛼⋅𝒛̅,
  //   𝒅 ← 𝒚̅ + (𝜗²⋅𝜂/𝛼)⋅𝒅,
  //   𝜗 ← ‖𝒖̂‖/𝜏, 
  //   𝑐𝑠 ← 𝟣/(𝟣 + 𝜗²)¹ᐟ²,
  //   𝜏 ← 𝜏⋅𝜗⋅𝑐𝑠, 𝜂 ← 𝛼⋅(𝑐𝑠)²,
  //   𝒙 ← 𝒙 + 𝜂⋅𝒅,
  //   𝗶𝗳 𝑚 = 𝟢: 
  //     𝒚̅, 𝒛̅ ← 𝒚, 𝒛.
  //   𝗲𝗻𝗱 𝗶𝗳
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // ----------------------
  for (stormSize_t m = 0; m <= 1; ++m) {
    stormBlas::Sub(uHatArr, uHatArr, zBarArr, alpha);
    stormBlas::Add(dArr, yBarArr, dArr, std::pow(theta, 2)*eta/alpha);
    theta = stormBlas::Norm2(uHatArr)/tau;
    stormReal_t const cs = 1.0/std::hypot(1.0, theta);
    tau *= theta*cs, eta = alpha*std::pow(cs, 2);
    stormBlas::Add(xArr, xArr, dArr, eta);
    if (m == 0) {
      std::swap(yBarArr, yArr), std::swap(zBarArr, zArr);
    }
  }

  // ----------------------
  // 𝜌̅ ← 𝜌, 
  // 𝜌 ← <𝒓̃⋅𝒖̂>, 𝛽 ← 𝜌/𝜌̅,
  // 𝒚 ← 𝒖̂ + 𝛽⋅𝒚̅,
  // 𝒛 ← 𝓐𝒚,
  // 𝒗 ← 𝒛̅ + 𝛽⋅𝒗,
  // 𝒗 ← 𝒛 + 𝛽⋅𝒗.
  // ----------------------
  stormReal_t const rhoBar = rho;
  rho = stormBlas::Dot(rTildeArr, uHatArr);
  stormReal_t const beta = rho/rhoBar;
  stormBlas::Add(yArr, uHatArr, yBarArr, beta);
  linOp.MatVec(zArr, yArr);
  stormBlas::Add(vArr, zBarArr, vArr, beta);
  stormBlas::Add(vArr, zArr, vArr, beta);

  // ----------------------
  // Compute the residual upper bound:
  // 𝜑 ← 𝜏⋅(𝟤𝑘 + 𝟥)¹ᐟ².
  // ----------------------
  stormSize_t const k = this->Iteration;
  stormReal_t const phi = tau*std::sqrt(2.0*k + 3.0);

  return phi;

} // stormTfqmrSolver<...>::Iterate

#endif // ifndef _STORM_SOLVER_TFQMR_
