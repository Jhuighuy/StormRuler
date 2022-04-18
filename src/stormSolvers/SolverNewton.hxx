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

#include <limits>

#include <stormBase.hxx>
#include <stormSolvers/Solver.hxx>

namespace Storm {

#if 0
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief The @c Newton method nonlinear operator equation solver.
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
template<class Vector>
class NewtonSolver : public IterativeSolver<Vector, tOperator> {
private:
  real_t Init(Vector const& xVec,
              Vector const& bVec,
              Operator<Vector> const& anyOp,
              Preconditioner<Vector> const* preOp) override final;

  real_t Iterate(Vector& xVec,
                 Vector const& bVec,
                 Operator<Vector> const& anyOp,
                 Preconditioner<Vector> const* preOp) override final;

}; // class NewtonSolver<...>
#endif

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief The first-order @c JFNK (Jacobian free-Newton-Krylov) \
///   nonlinear operator equation solver.
///
/// For the @c Newton iterations, computing of the Jacobian-vector
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
template<VectorLike Vector>
class JfnkSolver final : public IterativeSolver<Vector> {
private:
  Vector sVec_, tVec_, rVec_, wVec_;

  real_t Init(Vector const& xVec,
              Vector const& bVec,
              Operator<Vector> const& linOp,
              Preconditioner<Vector> const* preOp) override;

  real_t Iterate(Vector& xVec,
                 Vector const& bVec,
                 Operator<Vector> const& linOp,
                 Preconditioner<Vector> const* preOp) override;

}; // class JfnkSolver<...>

template<VectorLike Vector>
real_t JfnkSolver<Vector>::Init(Vector const& xVec,
                                Vector const& bVec,
                                Operator<Vector> const& linOp,
                                Preconditioner<Vector> const* preOp) {

  sVec_.Assign(xVec, false);
  tVec_.Assign(xVec, false);
  rVec_.Assign(xVec, false);
  wVec_.Assign(xVec, false);

  // ----------------------
  // 𝒘 ← 𝓐(𝒙),
  // 𝒓 ← 𝒃 - 𝒘.
  // ----------------------
  linOp.MatVec(wVec_, xVec);
  rVec_.Sub(bVec, wVec_);

  return rVec_.Norm2();  

} // JfnkSolver<...>::Init

template<VectorLike Vector>
real_t JfnkSolver<Vector>::Iterate(Vector& xVec,
                                   Vector const& bVec,
                                   Operator<Vector> const& linOp,
                                   Preconditioner<Vector> const* preOp) {

  // ----------------------
  // Solve the Jacobian equation:
  // 𝜇 ← (𝜀ₘ)¹ᐟ²⋅(1 + ‖𝒙‖)]¹ᐟ²,
  // 𝒕 ← 𝒓,
  // 𝒕 ← 𝓙(𝒙)⁻¹𝒓.
  // ----------------------
  static real_t const sqrtOfEpsilon = 
    std::sqrt(std::numeric_limits<real_t>::epsilon());
  real_t const mu = 
    sqrtOfEpsilon*std::sqrt(1.0 + xVec.Norm2());
  tVec_.Set(rVec_);
  {
    /// @todo Refactor me!
    /// @todo equation parameters!
    //call jacConvParams%Init(1e-8_dp, 1e-8_dp, 2000, 'Newton')
    //call LinSolve(mesh, 'GMRES', '', tVec_, rVec_, ApproxJacobianMatVecWithX, jacConvParams)
    auto solver = std::make_unique<BiCgStabSolver<Vector>>();
    solver->AbsoluteTolerance = 1.0e-8;
    solver->RelativeTolerance = 1.0e-8;
    auto op = MakeOperator<Vector>(
      [&](Vector& zVec, Vector const& yVec) {

        // ----------------------
        // Compute the Jacobian-vector product:
        // 𝛿 ← 𝜇⋅‖𝒚‖⁺,
        // 𝒔 ← 𝒙 + 𝛿⋅𝒚,
        // 𝒛 ← 𝓐(𝒔),
        // 𝒛 ← 𝛿⁺⋅𝒛 - 𝛿⁺⋅𝒘.
        // ----------------------
        real_t const delta = 
          Utils::SafeDivide(mu, yVec.Norm2());
        sVec_.Add(xVec, yVec, delta);
        linOp.MatVec(zVec, sVec_);
        real_t const deltaInverse = Utils::SafeDivide(1.0, delta);
        zVec.Sub(zVec, deltaInverse, wVec_, deltaInverse);

      });
    solver->Solve(tVec_, rVec_, *op);
  }

  // ----------------------
  // Update the solution and the residual:
  // 𝒙 ← 𝒙 + 𝒕,
  // 𝒘 ← 𝓐(𝒙),
  // 𝒓 ← 𝒃 - 𝒘.
  // ----------------------
  xVec.AddAssign(tVec_);
  linOp.MatVec(wVec_, xVec);
  rVec_.Sub(bVec, wVec_);

  return rVec_.Norm2();  

} // JfnkSolver<...>::Iterate

} // namespace Storm
