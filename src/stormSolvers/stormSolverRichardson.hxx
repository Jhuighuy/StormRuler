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
#ifndef _STORM_SOLVERS_RICHARDSON_HXX_
#define _STORM_SOLVERS_RICHARDSON_HXX_

#include <stormSolvers/stormSolver.hxx>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Solve a non-singular operator equation \
///   equation with the @c Richardson method.
///
/// References:
/// @verbatim
/// [1] ???
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Vector>
class stormRichardsonSolver final : public stormIterativeSolver<Vector> {
public:

  stormReal_t RelaxationFactor = 1.0e-4;

private:
  Vector rVec_, zVec_;

  stormReal_t Init(Vector const& xVec,
                   Vector const& bVec,
                   stormOperator<Vector> const& linOp,
                   stormPreconditioner<Vector> const* preOp) override;

  stormReal_t Iterate(Vector& xVec,
                      Vector const& bVec,
                      stormOperator<Vector> const& linOp,
                      stormPreconditioner<Vector> const* preOp) override;

}; // class stormRichardsonSolver<...>

template<class Vector>
stormReal_t stormRichardsonSolver<Vector>::
                      Init(Vector const& xVec,
                           Vector const& bVec,
                           stormOperator<Vector> const& linOp,
                           stormPreconditioner<Vector> const* preOp) {

  rVec_.Assign(xVec, false);
  if (preOp != nullptr) {
    zVec_.Assign(xVec, false);
  }

  // ----------------------
  // Initialize:
  // 𝒓 ← 𝓐𝒙,
  // 𝒓 ← 𝒃 - 𝒓,
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝒛 ← 𝒓,
  //   𝒓 ← 𝓟𝒛.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  linOp.MatVec(rVec_, xVec);
  rVec_.Sub(bVec, rVec_);
  if (preOp != nullptr) {
    std::swap(zVec_, rVec_);
    preOp->MatVec(rVec_, zVec_);
  }

  return rVec_.Norm2();

} // stormRichardsonSolver<...>::Init

template<class Vector>
stormReal_t stormRichardsonSolver<Vector>::
                        Iterate(Vector& xVec,
                                Vector const& bVec,
                                stormOperator<Vector> const& linOp,
                                stormPreconditioner<Vector> const* preOp) {

  stormReal_t const& omega = RelaxationFactor;

  // ----------------------
  // Update the solution and the residual:
  // 𝒙 ← 𝒙 + 𝜔⋅𝒓,
  // 𝒓 ← 𝓐𝒙,
  // 𝒓 ← 𝒃 - 𝒓,
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝒛 ← 𝒓,
  //   𝒓 ← 𝓟𝒛.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  stormBlas::Add(xVec, xVec, rVec_, omega);
  linOp.MatVec(rVec_, xVec);
  rVec_.Sub(bVec, rVec_);
  if (preOp != nullptr) {
    std::swap(zVec_, rVec_);
    preOp->MatVec(rVec_, zVec_);
  }

  return rVec_.Norm2();

} // stormRichardsonSolver<...>::Iterate

#endif // ifndef _STORM_SOLVERS_RICHARDSON_HXX_
