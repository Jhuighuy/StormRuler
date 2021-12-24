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
#ifndef _STORM_SOLVER_GMRES_HXX_
#define _STORM_SOLVER_GMRES_HXX_

#include <stormSolvers/stormSolver.hxx>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Solve a linear operator equation: [𝓟]𝓐𝒙 = [𝓟]𝒃, using
///   the monstrous @c GMRES (Generalized Minimal Residual) method.
///
/// The classical @c GMRES(𝑚) implementation with restarts
/// after 𝑚 iterations is used.
///
/// Preconditioned residual norm, ‖[𝓟]𝒓‖, where 𝒓 = 𝒃 - 𝓐𝒙, is reported.
/// 
/// @c GMRES may be applied to the singular problems, and the square
/// least squares problems: ‖(𝓐[𝓟]𝒚 - 𝒃)‖₂ → 𝘮𝘪𝘯, 𝒙 = [𝓟]𝒚,
/// although convergeance to minimum norm solution is not guaranteed
/// (is this true?).
///
/// References:
/// @verbatim
/// [1] Saad and M.H. Schultz,
///     "GMRES: A generalized minimal residual algorithm for solving
///      nonsymmetric linear systems",
///     SIAM J. Sci. Stat. Comput., 7:856–869, 1986.
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class tArray>
class stormGmresSolver final : public stormRestartableSolver<tArray> {
private:

  void PreInit(tArray& xArr,
               const tArray& bArr, 
               bool hasPreOp) override;

  stormReal_t ReInit(tArray& xArr,
                     const tArray& bArr,
                     const stormOperator<tArray>& linOp,
                     const stormPreconditioner<tArray>* preOp) override;

  stormReal_t ReIterate(stormSize_t k,
                        tArray& xArr,
                        const tArray& bArr,
                        const stormOperator<tArray>& linOp,
                        const stormPreconditioner<tArray>* preOp) override;

  void ReFinalize(stormSize_t k,
                  tArray& xArr,
                  const tArray& bArr,
                  const stormOperator<tArray>& linOp,
                  const stormPreconditioner<tArray>* preOp) override;

}; // class stormGmresSolver<...>

template<class tArray>
void stormGmresSolver<tArray>::PreInit(tArray& xArr,
                                       const tArray& bArr, 
                                       bool hasPreOp) {
  // ----------------------
  // Allocate the intermediate arrays:
  // ----------------------
  //stormUtils::AllocLike(xArr, pArr, rArr, rTildeArr, sArr, tArr, vArr);
  if (hasPreOp) {
    //stormUtils::AllocLike(xArr, wArr, yArr, zArr);
  }

  _STORM_NOT_IMPLEMENTED_();

} // stormGmresSolver<...>::Init

template<class tArray>
stormReal_t stormGmresSolver<tArray>::ReInit(tArray& xArr,
                                             const tArray& bArr,
                                             const stormOperator<tArray>& linOp,
                                             const stormPreconditioner<tArray>* preOp) {

  _STORM_NOT_IMPLEMENTED_();

} // stormGmresSolver<...>::ReInit

template<class tArray>
stormReal_t stormGmresSolver<tArray>::ReIterate(stormSize_t k,
                                                tArray& xArr,
                                                const tArray& bArr,
                                                const stormOperator<tArray>& linOp,
                                                const stormPreconditioner<tArray>* preOp) {

  _STORM_NOT_IMPLEMENTED_();

} // stormGmresSolver<...>::Iterate

template<class tArray>
void stormGmresSolver<tArray>::ReFinalize(stormSize_t k,
                                          tArray& xArr,
                                          const tArray& bArr,
                                          const stormOperator<tArray>& linOp,
                                          const stormPreconditioner<tArray>* preOp) {

  _STORM_NOT_IMPLEMENTED_();

} // stormGmresSolver<...>::Finalize

#endif // ifndef _STORM_SOLVER_GMRES_HXX_
