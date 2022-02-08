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

#include <stormBase.hxx>
#include <stormSolvers/stormSolver.hxx>

namespace Storm {

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
class RichardsonSolver final : public IterativeSolver<Vector> {
public:
  real_t RelaxationFactor = 1.0e-4;

private:
  Vector rVec_, zVec_;

  real_t Init(Vector const& xVec,
              Vector const& bVec,
              Operator<Vector> const& linOp,
              Preconditioner<Vector> const* preOp) override;

  real_t Iterate(Vector& xVec,
                 Vector const& bVec,
                 Operator<Vector> const& linOp,
                 Preconditioner<Vector> const* preOp) override;

}; // class RichardsonSolver<...>

template<class Vector>
real_t RichardsonSolver<Vector>::Init(Vector const& xVec,
                                      Vector const& bVec,
                                      Operator<Vector> const& linOp,
                                      Preconditioner<Vector> const* preOp) {

  rVec_.Assign(xVec, false);
  if (preOp != nullptr) {
    zVec_.Assign(xVec, false);
  }

  // ----------------------
  // Initialize:
  // 𝒓 ← 𝒃 - 𝓐𝒙,
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝒛 ← 𝒓,
  //   𝒓 ← 𝓟𝒛.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  linOp.Residual(rVec_, bVec, xVec);
  if (preOp != nullptr) {
    std::swap(zVec_, rVec_);
    preOp->MatVec(rVec_, zVec_);
  }

  return Blas::Norm2(rVec_);

} // RichardsonSolver<...>::Init

template<class Vector>
real_t RichardsonSolver<Vector>::Iterate(Vector& xVec,
                                         Vector const& bVec,
                                         Operator<Vector> const& linOp,
                                         Preconditioner<Vector> const* preOp) {

  real_t const& omega = RelaxationFactor;

  // ----------------------
  // Update the solution and the residual:
  // 𝒙 ← 𝒙 + 𝜔⋅𝒓,
  // 𝒓 ← 𝒃 - 𝓐𝒙,
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝒛 ← 𝒓,
  //   𝒓 ← 𝓟𝒛.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  Blas::Add(xVec, xVec, rVec_, omega);
  linOp.Residual(rVec_, bVec, xVec);
  if (preOp != nullptr) {
    std::swap(zVec_, rVec_);
    preOp->MatVec(rVec_, zVec_);
  }

  return Blas::Norm2(rVec_);

} // RichardsonSolver<...>::Iterate

} // namespace Storm
