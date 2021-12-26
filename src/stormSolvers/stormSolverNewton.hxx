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
#ifndef _STORM_SOLVER_NEWTON_HXX_
#define _STORM_SOLVER_NEWTON_HXX_

#include <limits>

#include <stormSolvers/stormSolver.hxx>

#if 0
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Solve a nonlinear operator equation with the Newton's method.
///
/// The classical Newton iterations are based on the linearization 
/// of 𝓐(𝒙) near 𝒙: 
///
/// 𝓐(𝒙̂) ≈ 𝓐(𝒙) + [∂𝓐(𝒙)/∂𝒙](𝒙̂ - 𝒙) = 𝒃, 
///
/// or, alternatively:
///
/// [∂𝓐(𝒙)/∂𝒙]𝒕 = 𝒓, 𝒕 = 𝒙̂ - 𝒙, 𝒓 = 𝒃 - 𝓐(𝒙)
///
/// where 𝒙 and 𝒙̂ are the current and updated solution vectors.
/// Therefore, a linear equation has to be solved on each iteration,
/// linear operator 𝓙(𝒙) ≈ ∂𝓐(𝒙)/∂𝒙 for computing Jacobian-vector 
/// products is required.
///
/// References:
/// @verbatim
/// [1] ???
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class tArray, class tOperator>
class stormNewtonSolver : public stormIterativeSolver<tArray, tOperator> {
private:
  stormReal_t Init(tArray& xArr,
                   const tArray& bArr,
                   const tOperator& anyOp,
                   const tOperator* preOp) override final;

  stormReal_t Iterate(tArray& xArr,
                      const tArray& bArr,
                      const tOperator& anyOp,
                      const tOperator* preOp) override final;

}; // class stormNewtonSolver<...>
#endif

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Solve a nonlinear operator equation with the \
///   first order @c JFNK (Jacobian free-Newton-Krylov) method.
///
/// For the Newton iterations, computing of the Jacobian-vector
/// products 𝒛 = 𝓙(𝒙)𝒚, where 𝓙(𝒙) ≈ ∂𝓐(𝒙)/∂𝒙 is required.
/// Consider the expansion:
///
/// 𝓐(𝒙 + 𝛿⋅𝒚) = 𝓐(𝒙) + 𝛿⋅[∂𝓐(𝒙)/∂𝒙]𝒚 + 𝓞(𝛿²),
///
/// where 𝛿 is some small number. Therefore,
///
/// 𝓙(𝒙)𝒚 = [𝓐(𝒙 + 𝛿⋅𝒚) - 𝓐(𝒙)]/𝛿 = [∂𝓐(𝒙)/∂𝒙]𝒚 + 𝓞(𝛿).
///
/// Expression above may be used as the formula for computing
/// the (approximate) Jacobian-vector products. Parameter 𝛿 is commonly 
/// defined as [1]:
///
/// 𝛿 = 𝜇⋅‖𝒚‖⁺, 𝜇 = (𝜀ₘ)¹ᐟ²⋅(1 + ‖𝒙‖)¹ᐟ²,
///
/// where 𝜀ₘ is the machine roundoff, ‖𝒚‖⁺ is the pseudo-inverse to ‖𝒚‖.
///
/// References:
/// @verbatim
/// [1] Liu, Wei, Lilun Zhang, Ying Zhong, Yongxian Wang, 
///     Yonggang Che, Chuanfu Xu and Xinghua Cheng. 
///     “CFD High-order Accurate Scheme JFNK Solver.” 
///     Procedia Engineering 61 (2013): 9-15.
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class tArray>
class stormJfnkSolver final : public stormIterativeSolver<tArray> {
private:
  tArray sArr, tArr, rArr, wArr;

  stormReal_t Init(tArray& xArr,
                   const tArray& bArr,
                   const stormOperator<tArray>& linOp,
                   const stormPreconditioner<tArray>* preOp) override;

  stormReal_t Iterate(tArray& xArr,
                      const tArray& bArr,
                      const stormOperator<tArray>& linOp,
                      const stormPreconditioner<tArray>* preOp) override;

}; // class stormJfnkSolver<...>

template<class tArray>
stormReal_t stormJfnkSolver<tArray>::Init(tArray& xArr,
                                          const tArray& bArr,
                                          const stormOperator<tArray>& linOp,
                                          const stormPreconditioner<tArray>* preOp) {

  // ----------------------
  // Allocate the intermediate arrays:
  // ----------------------
  stormUtils::AllocLike(xArr, sArr, tArr, rArr, wArr);

  // ----------------------
  // Compute residual:
  // 𝒘 ← 𝓐(𝒙),
  // 𝒓 ← 𝒃 - 𝒘.
  // ----------------------
  linOp.MatVec(wArr, xArr);
  stormUtils::Sub(rArr, bArr, wArr);

  return stormUtils::Norm2(rArr);  

} // stormJfnkSolver<...>::Init

template<class tArray>
stormReal_t stormJfnkSolver<tArray>::Iterate(tArray& xArr,
                                             const tArray& bArr,
                                             const stormOperator<tArray>& linOp,
                                             const stormPreconditioner<tArray>* preOp) {

  // ----------------------
  // Solve the Jacobian equation:
  // 𝜇 ← (𝜀ₘ)¹ᐟ²⋅(1 + ‖𝒙‖)]¹ᐟ²,
  // 𝒕 ← 𝒓,
  // 𝒕 ← 𝓙(𝒙)⁻¹𝒓,
  // 𝒙 ← 𝒙 + 𝒕.
  // ----------------------
  static const stormReal_t sqrtOfEpsilon = 
    std::sqrt(std::numeric_limits<stormReal_t>::epsilon());
  const stormReal_t mu = 
    sqrtOfEpsilon*std::sqrt(1.0 + stormUtils::Norm2(xArr));
  stormUtils::Set(tArr, rArr);
  {
    /// @todo Refactor me!
    /// @todo equation parameters!
    //call jacConvParams%Init(1e-8_dp, 1e-8_dp, 2000, 'Newton')
    //call LinSolve(mesh, 'GMRES', '', tArr, rArr, ApproxJacobianMatVecWithX, jacConvParams)
    auto solver = std::make_unique<stormBiCgStabSolver<tArray>>();
    solver->AbsoluteTolerance = 1.0e-8;
    solver->RelativeTolerance = 1.0e-8;
    auto op = stormMakeOperator<tArray>(
      [&](stormArray& zArr, const stormArray& yArr) {

        // ----------------------
        // Compute the Jacobian-vector product:
        // 𝛿 ← 𝜇⋅‖𝒚‖⁺,
        // 𝒔 ← 𝒙 + 𝛿⋅𝒚,
        // 𝒛 ← 𝓐(𝒔),
        // 𝒛 ← 𝛿⁺⋅𝒛 - 𝛿⁺⋅𝒘.
        // ----------------------
        const stormReal_t delta = 
          stormUtils::SafeDivide(mu, stormUtils::Norm2(yArr));
        stormUtils::Add(sArr, xArr, yArr, delta);
        linOp.MatVec(zArr, sArr);
        const stormReal_t deltaInverse = stormUtils::SafeDivide(1.0, delta);
        stormUtils::Sub(zArr, zArr, wArr, deltaInverse, deltaInverse);

      });
    solver->Solve(tArr, rArr, *op);
  }
  stormUtils::Add(xArr, xArr, tArr);

  // ----------------------
  // Compute residual:
  // 𝒘 ← 𝓐(𝒙),
  // 𝒓 ← 𝒃 - 𝒘.
  // ----------------------
  linOp.MatVec(wArr, xArr);
  stormUtils::Sub(rArr, bArr, wArr);

  return stormUtils::Norm2(rArr);  

} // stormJfnkSolver<...>::Iterate

#endif // ifndef _STORM_SOLVER_NEWTON_HXX_
