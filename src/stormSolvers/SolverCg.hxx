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

#include <cmath>

#include <stormBase.hxx>
#include <stormSolvers/Solver.hxx>
#include <stormSolvers/Vector.hxx>

namespace Storm {

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief The @c CG (Conjugate Gradients) linear self-adjoint
///   definite operator equation solver.
///
/// @c CG may be applied to the consistent singular problems,
/// it converges towards..
///
/// References:
/// @verbatim
/// [1] Hestenes, Magnus R. and Eduard Stiefel.
///     “Methods of conjugate gradients for solving linear systems.”
///     Journal of research of the National
///     Bureau of Standards 49 (1952): 409-435.
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<VectorLike Vector>
class CgSolver final : public IterativeSolver<Vector> {
private:

  real_t gamma_;
  Vector pVec_, rVec_, zVec_;

  real_t Init(Vector const& xVec, Vector const& bVec,
              Operator<Vector> const& linOp,
              Preconditioner<Vector> const* preOp) override;

  real_t Iterate(Vector& xVec, Vector const& bVec,
                 Operator<Vector> const& linOp,
                 Preconditioner<Vector> const* preOp) override;

}; // class CgSolver

template<VectorLike Vector>
real_t CgSolver<Vector>::Init(Vector const& xVec, Vector const& bVec,
                              Operator<Vector> const& linOp,
                              Preconditioner<Vector> const* preOp) {
  pVec_.Assign(xVec, false);
  rVec_.Assign(xVec, false);
  zVec_.Assign(xVec, false);

  // Initialize:
  // ----------------------
  // 𝒓 ← 𝒃 - 𝓐𝒙.
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝒛 ← 𝓟𝒓,
  //   𝒑 ← 𝒛,
  //   𝛾 ← <𝒓⋅𝒛>,
  // 𝗲𝗹𝘀𝗲:
  //   𝒑 ← 𝒓,
  //   𝛾 ← <𝒓⋅𝒓>.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  linOp.Residual(rVec_, bVec, xVec);
  if (preOp != nullptr) {
    preOp->MatVec(zVec_, rVec_);
    pVec_ <<= zVec_;
    gamma_ = Blas::Dot(rVec_, zVec_);
  } else {
    pVec_ <<= rVec_;
    gamma_ = Blas::Dot(rVec_, rVec_);
  }

  return (preOp != nullptr) ? Blas::Norm2(rVec_) : std::sqrt(gamma_);

} // CgSolver::Init

template<VectorLike Vector>
real_t CgSolver<Vector>::Iterate(Vector& xVec, Vector const& bVec,
                                 Operator<Vector> const& linOp,
                                 Preconditioner<Vector> const* preOp) {
  // Iterate:
  // ----------------------
  // 𝒛 ← 𝓐𝒑,
  // 𝛼 ← 𝛾/<𝒑⋅𝒛>,
  // 𝒙 ← 𝒙 + 𝛼⋅𝒑,
  // 𝒓 ← 𝒓 - 𝛼⋅𝒛.
  // ----------------------
  linOp.MatVec(zVec_, pVec_);
  real_t const alpha{Utils::SafeDivide(gamma_, Blas::Dot(pVec_, zVec_))};
  xVec += alpha * pVec_;
  rVec_ -= alpha * zVec_;

  // ----------------------
  // 𝛾̅ ← 𝛾,
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝒛 ← 𝓟𝒓,
  //   𝛾 ← <𝒓⋅𝒛>,
  // 𝗲𝗹𝘀𝗲:
  //   𝛾 ← <𝒓⋅𝒓>.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  real_t const gammaBar{gamma_};
  if (preOp != nullptr) {
    preOp->MatVec(zVec_, rVec_);
    gamma_ = Blas::Dot(rVec_, zVec_);
  } else {
    gamma_ = Blas::Dot(rVec_, rVec_);
  }

  // ----------------------
  // 𝛽 ← 𝛾/𝛾̅,
  // 𝒑 ← (𝓟 ≠ 𝗻𝗼𝗻𝗲 ? 𝒛 : 𝒓) + 𝛽⋅𝒑.
  // ----------------------
  real_t const beta = Utils::SafeDivide(gamma_, gammaBar);
  pVec_ <<= (preOp != nullptr ? zVec_ : rVec_) + beta * pVec_;

  return (preOp != nullptr) ? Blas::Norm2(rVec_) : std::sqrt(gamma_);

} // CgSolver::Iterate

} // namespace Storm
