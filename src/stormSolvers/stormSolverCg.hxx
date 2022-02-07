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
#ifndef _STORM_SOLVER_CG_HXX_
#define _STORM_SOLVER_CG_HXX_

#include <cmath>

#include <stormBase.hxx>
#include <stormSolvers/stormSolver.hxx>

_STORM_NAMESPACE_BEGIN_

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Solve a linear self-adjoint definite operator equation \
///   equation with the @c CG (Conjugate Gradients) method.
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
template<class Vector>
class CgSolver final : public IterativeSolver<Vector> {
private:
  Real_t alpha_;
  Vector pVec_, rVec_, zVec_;

  Real_t Init(Vector const& xVec,
              Vector const& bVec,
              Operator<Vector> const& linOp,
              Preconditioner<Vector> const* preOp) override;

  Real_t Iterate(Vector& xVec,
                 Vector const& bVec,
                 Operator<Vector> const& linOp,
                 Preconditioner<Vector> const* preOp) override;

}; // class CgSolver<...>

template<class Vector>
Real_t CgSolver<Vector>::Init(Vector const& xVec,
                              Vector const& bVec,
                              Operator<Vector> const& linOp,
                              Preconditioner<Vector> const* preOp) {

  pVec_.Assign(xVec, false);
  rVec_.Assign(xVec, false);
  zVec_.Assign(xVec, false);

  // ----------------------
  // Initialize:
  // 𝒓 ← 𝓐𝒙,
  // 𝒓 ← 𝒃 - 𝒓.
  // ----------------------
  linOp.MatVec(rVec_, xVec);
  stormBlas::Sub(rVec_, bVec, rVec_);

  // ----------------------
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝒛 ← 𝓟𝒓,
  //   𝒑 ← 𝒛,
  //   𝛼 ← <𝒓⋅𝒛>,
  // 𝗲𝗹𝘀𝗲:
  //   𝒑 ← 𝒓,
  //   𝛼 ← <𝒓⋅𝒓>.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  if (preOp != nullptr) {
    preOp->MatVec(zVec_, rVec_);
    stormBlas::Set(pVec_, zVec_);
    alpha_ = stormBlas::Dot(rVec_, zVec_);
  } else {
    stormBlas::Set(pVec_, rVec_);
    alpha_ = stormBlas::Dot(rVec_, rVec_);
  }

  return (preOp != nullptr) ? stormBlas::Norm2(rVec_) : std::sqrt(alpha_);

} // CgSolver<...>::Init

template<class Vector>
Real_t CgSolver<Vector>::Iterate(Vector& xVec,
                                 Vector const& bVec,
                                 Operator<Vector> const& linOp,
                                 Preconditioner<Vector> const* preOp) {

  // ----------------------
  // Iterate:
  // 𝒛 ← 𝓐𝒑,
  // 𝛼̅ ← 𝛼,
  // 𝛼 ← 𝛼/<𝒑⋅𝒛>,
  // 𝒙 ← 𝒙 + 𝛼⋅𝒑,
  // 𝒓 ← 𝒓 - 𝛼⋅𝒛,
  // ----------------------
  linOp.MatVec(zVec_, pVec_);
  Real_t const alphaBar = alpha_;
  stormUtils::SafeDivideEquals(alpha_, stormBlas::Dot(pVec_, zVec_));
  stormBlas::Add(xVec, xVec, pVec_, alpha_);
  stormBlas::Sub(rVec_, rVec_, zVec_, alpha_);

  // ----------------------
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝒛 ← 𝓟𝒓,
  //   𝛼 ← <𝒓⋅𝒛>,
  // 𝗲𝗹𝘀𝗲:
  //   𝛼 ← <𝒓⋅𝒓>.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  if (preOp != nullptr) {
    preOp->MatVec(zVec_, rVec_);
    alpha_ = stormBlas::Dot(rVec_, zVec_);
  } else {
    alpha_ = stormBlas::Dot(rVec_, rVec_);
  }

  // ----------------------
  // 𝛽 ← 𝛼/𝛼̅,
  // 𝒑 ← (𝓟 ≠ 𝗻𝗼𝗻𝗲 ? 𝒛 : 𝒓) + 𝛽⋅𝒑.
  // ----------------------
  Real_t const beta = stormUtils::SafeDivide(alpha_, alphaBar);
  stormBlas::Add(pVec_, (preOp != nullptr ? zVec_ : rVec_), pVec_, beta);

  return (preOp != nullptr) ? stormBlas::Norm2(rVec_) : std::sqrt(alpha_);

} // CgSolver<...>::Iterate

_STORM_NAMESPACE_END_

#endif // ifndef _STORM_SOLVER_CG_HXX_
