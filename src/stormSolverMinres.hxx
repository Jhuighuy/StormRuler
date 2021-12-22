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
#ifndef _STORM_SOLVER_MINRES_HXX_
#define _STORM_SOLVER_MINRES_HXX_

#include <stormSolver.hxx>

#include <cmath>

/// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ///
/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ///

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Solve a linear self-adjoint indefinite operator equation 
///   [𝓜]𝓐[𝓜ᵀ]𝒚 = [𝓜]𝒃, [𝓜ᵀ]𝒚 = 𝒙, [𝓜𝓜ᵀ = 𝓟], using the @c MINRES method.
///
/// @c MINRES can be applied to the singular problems, and the self-adjoint
/// least squares problems: ‖[𝓜](𝓐[𝓜ᵀ]𝒚 - 𝒃)‖₂ → 𝘮𝘪𝘯, 𝒙 = [𝓜ᵀ]𝒚, 
/// although convergeance to minimum norm solution is not guaranteed.
///
/// @note Despite 𝓐 may be indefinite, a positive-definite \
///   preconditioner 𝓟 is explicitly required.
///
/// References:
/// @verbatim
/// [1] Paige, C. and M. Saunders. 
///     “Solution of Sparse Indefinite Systems of Linear Equations.” 
///     SIAM Journal on Numerical Analysis 12 (1975): 617-629.
/// [2] Choi, S.-C. T.
///     “Iterative Methods for Singular Linear Equations and 
///     Least-Squares Problems” PhD thesis, ICME, Stanford University.
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class tArray, class tOperator>
class stormMinresSolver final : public stormIterativeSolver<tArray, tOperator> {
private:

protected:

  /// @brief Initialize the @c MINRES solver.
  ///
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param linOp Self-adjoint linear operator, 𝓐(𝒙).
  /// @param preOp Self-adjoint positive definite linear preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Preconditioned residual norm, \
  ///   square root of <𝒓⋅𝒛>, where 𝒓 = 𝒃 - 𝓐𝒙  and 𝒛 = [𝓟]𝒓.
  stormReal_t Init(tArray& xArr,
                   const tArray& bArr,
                   const tOperator& linOp,
                   const tOperator* preOp) override final;

  /// @brief Iterate the @c MINRES solver.
  ///
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param linOp Self-adjoint linear operator, 𝓐(𝒙).
  /// @param preOp Self-adjoint positive definite linear preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Preconditioned residual norm, \
  ///   square root of <𝒓⋅𝒛>, where 𝒓 = 𝒃 - 𝓐𝒙  and 𝒛 = [𝓟]𝒓.
  stormReal_t Iterate(tArray& xArr,
                      const tArray& bArr,
                      const tOperator& linOp,
                      const tOperator* preOp) override final;

}; // class stormMinresSolver<...>

template<class tArray, class tOperator>
stormReal_t stormMinresSolver<tArray, tOperator>::Init(tArray& xArr,
                                                       const tArray& bArr,
                                                       const tOperator& linOp,
                                                       const tOperator* preOp) {
  // ----------------------
  // Allocate the intermediate arrays:
  // ----------------------
  //stormUtils::AllocLike(xArr, pArr, rArr, tArr);
  if (preOp != nullptr) {
    //stormUtils::AllocLike(xArr, zArr);
  }

  _STORM_NOT_IMPLEMENTED_();

} // stormMinresSolver<...>::Init

template<class tArray, class tOperator>
stormReal_t stormMinresSolver<tArray, tOperator>::Iterate(tArray& xArr,
                                                          const tArray& bArr,
                                                          const tOperator& linOp,
                                                          const tOperator* preOp) {

  _STORM_NOT_IMPLEMENTED_();

} // stormMinresSolver<...>::Iterate

/// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ///
/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ///

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Solve a linear operator equation: [𝓟]𝓐𝒙 = [𝓟]𝒃, using 
///   the monstrous Generalized Minimal Residual method (@c GMRES).
///
/// The classical @c GMRES(𝑚) implementation with restarts
/// after 𝑚 iterations is used.
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
template<class tArray, class tOperator>
class stormGmresSolver final : public stormIterativeSolver<tArray, tOperator> {
private:

protected:

  /// @brief Initialize the @c GMRES solver.
  ///
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param linOp Linear operator, 𝓐(𝒙).
  /// @param preOp Linear preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Preconditioned residual norm, ‖[𝓟]𝒓‖, where 𝒓 = 𝒃 - 𝓐𝒙.
  stormReal_t Init(tArray& xArr,
                   const tArray& bArr,
                   const tOperator& linOp,
                   const tOperator* preOp) override final;

  /// @brief Iterate the @c GMRES solver.
  ///
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param linOp Linear operator, 𝓐(𝒙).
  /// @param preOp Linear preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Preconditioned residual norm, ‖[𝓟]𝒓‖, where 𝒓 = 𝒃 - 𝓐𝒙.
  stormReal_t Iterate(tArray& xArr,
                      const tArray& bArr,
                      const tOperator& linOp,
                      const tOperator* preOp) override final;

  /// @brief Finalize the @c GMRES iterations.
  ///
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param linOp Linear operator, 𝓐(𝒙).
  /// @param preOp Linear preconditioner operator, 𝓟(𝒙).
  void Finalize(tArray& xArr,
                const tArray& bArr,
                const tOperator& linOp,
                const tOperator* preOp = nullptr) override final;

}; // class stormGmresSolver<...>

template<class tArray, class tOperator>
stormReal_t stormGmresSolver<tArray, tOperator>::Init(tArray& xArr,
                                                      const tArray& bArr,
                                                      const tOperator& linOp,
                                                      const tOperator* preOp) {
  // ----------------------
  // Allocate the intermediate arrays:
  // ----------------------
  //stormUtils::AllocLike(xArr, pArr, rArr, rTildeArr, sArr, tArr, vArr);
  if (preOp != nullptr) {
    //stormUtils::AllocLike(xArr, wArr, yArr, zArr);
  }

  _STORM_NOT_IMPLEMENTED_();

} // stormGmresSolver<...>::Init

template<class tArray, class tOperator>
stormReal_t stormGmresSolver<tArray, tOperator>::Iterate(tArray& xArr,
                                                         const tArray& bArr,
                                                         const tOperator& linOp,
                                                         const tOperator* preOp) {

  _STORM_NOT_IMPLEMENTED_();

} // stormGmresSolver<...>::Iterate

template<class tArray, class tOperator>
void stormGmresSolver<tArray, tOperator>::Finalize(tArray& xArr,
                                                  const tArray& bArr,
                                                  const tOperator& linOp,
                                                  const tOperator* preOp) {

  _STORM_NOT_IMPLEMENTED_();

} // stormGmresSolver<...>::Finalize

/// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ///
/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ///

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Solve a linear operator equation: [𝓟]𝓐𝒙 = [𝓟]𝒃, using 
///   the Transpose-Free Quasi-Minimal Residual method (@c TFQMR).
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
template<class tArray, class tOperator>
class stormTfqmrSolver final : public stormIterativeSolver<tArray, tOperator> {
private:

protected:

  /// @brief Initialize the @c TFQMR solver.
  ///
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param linOp Linear operator, 𝓐(𝒙).
  /// @param preOp Linear preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Preconditioned residual norm, ‖[𝓟]𝒓‖, where 𝒓 = 𝒃 - 𝓐𝒙.
  stormReal_t Init(tArray& xArr,
                   const tArray& bArr,
                   const tOperator& linOp,
                   const tOperator* preOp) override final;

  /// @brief Iterate the @c TFQMR solver.
  ///
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param linOp Linear operator, 𝓐(𝒙).
  /// @param preOp Linear preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Preconditioned residual norm, ‖[𝓟]𝒓‖, where 𝒓 = 𝒃 - 𝓐𝒙.
  stormReal_t Iterate(tArray& xArr,
                      const tArray& bArr,
                      const tOperator& linOp,
                      const tOperator* preOp) override final;

}; // class stormTfqmrSolver<...>

template<class tArray, class tOperator>
stormReal_t stormTfqmrSolver<tArray, tOperator>::Init(tArray& xArr,
                                                      const tArray& bArr,
                                                      const tOperator& linOp,
                                                      const tOperator* preOp) {
  // ----------------------
  // Allocate the intermediate arrays:
  // ----------------------
  //stormUtils::AllocLike(xArr, pArr, rArr, rTildeArr, sArr, tArr, vArr);
  if (preOp != nullptr) {
    //stormUtils::AllocLike(xArr, wArr, yArr, zArr);
  }

  _STORM_NOT_IMPLEMENTED_();

} // stormTfqmrSolver<...>::Init

template<class tArray, class tOperator>
stormReal_t stormTfqmrSolver<tArray, tOperator>::Iterate(tArray& xArr,
                                                         const tArray& bArr,
                                                         const tOperator& linOp,
                                                         const tOperator* preOp) {

  _STORM_NOT_IMPLEMENTED_();

} // stormTfqmrSolver<...>::Iterate

/// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ///
/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ///

#endif // ifndef _STORM_SOLVER_MINRES_HXX_
