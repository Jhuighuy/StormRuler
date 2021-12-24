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
/// @brief Solve a linear operator equation: [𝓟]𝓐𝒙 = [𝓟]𝒃, using
///   the @c TFQMR (Transpose-Free Quasi-Minimal Residual) method.
///
/// Preconditioned residual norm, ‖[𝓟]𝒓‖, where 𝒓 = 𝒃 - 𝓐𝒙, is reported.
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

  stormReal_t Init(tArray& xArr,
                   const tArray& bArr,
                   const stormOperator<tArray>& linOp,
                   const stormPreconditioner<tArray>* preOp) override;

  stormReal_t Iterate(tArray& xArr,
                      const tArray& bArr,
                      const stormOperator<tArray>& linOp,
                      const stormPreconditioner<tArray>* preOp) override;

}; // class stormTfqmrSolver<...>

template<class tArray>
stormReal_t stormTfqmrSolver<tArray>::Init(tArray& xArr,
                                           const tArray& bArr,
                                           const stormOperator<tArray>& linOp,
                                           const stormPreconditioner<tArray>* preOp) {
  // ----------------------
  // Allocate the intermediate arrays:
  // ----------------------
  //stormUtils::AllocLike(xArr, pArr, rArr, rTildeArr, sArr, tArr, vArr);
  if (preOp != nullptr) {
    //stormUtils::AllocLike(xArr, wArr, yArr, zArr);
  }

  _STORM_NOT_IMPLEMENTED_();

} // stormTfqmrSolver<...>::Init

template<class tArray>
stormReal_t stormTfqmrSolver<tArray>::Iterate(tArray& xArr,
                                              const tArray& bArr,
                                              const stormOperator<tArray>& linOp,
                                              const stormPreconditioner<tArray>* preOp) {

  _STORM_NOT_IMPLEMENTED_();

} // stormTfqmrSolver<...>::Iterate

#endif // ifndef _STORM_SOLVER_TFQMR_
