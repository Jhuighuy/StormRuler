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
  Vector p_vec_, r_vec_, z_vec_;

  real_t init(Vector const& x_vec, Vector const& b_vec,
              Operator<Vector> const& lin_op,
              Preconditioner<Vector> const* pre_op) override;

  real_t iterate(Vector& x_vec, Vector const& b_vec,
                 Operator<Vector> const& lin_op,
                 Preconditioner<Vector> const* pre_op) override;

}; // class CgSolver

template<VectorLike Vector>
real_t CgSolver<Vector>::init(Vector const& x_vec, Vector const& b_vec,
                              Operator<Vector> const& lin_op,
                              Preconditioner<Vector> const* pre_op) {
  p_vec_.assign(x_vec, false);
  r_vec_.assign(x_vec, false);
  z_vec_.assign(x_vec, false);

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
  lin_op.Residual(r_vec_, b_vec, x_vec);
  if (pre_op != nullptr) {
    pre_op->mul(z_vec_, r_vec_);
    p_vec_ <<= z_vec_;
    gamma_ = dot_product(r_vec_, z_vec_);
  } else {
    p_vec_ <<= r_vec_;
    gamma_ = dot_product(r_vec_, r_vec_);
  }

  return (pre_op != nullptr) ? norm_2(r_vec_) : std::sqrt(gamma_);

} // CgSolver::init

template<VectorLike Vector>
real_t CgSolver<Vector>::iterate(Vector& x_vec, Vector const& b_vec,
                                 Operator<Vector> const& lin_op,
                                 Preconditioner<Vector> const* pre_op) {
  // Iterate:
  // ----------------------
  // 𝒛 ← 𝓐𝒑,
  // 𝛼 ← 𝛾/<𝒑⋅𝒛>,
  // 𝒙 ← 𝒙 + 𝛼⋅𝒑,
  // 𝒓 ← 𝒓 - 𝛼⋅𝒛.
  // ----------------------
  lin_op.mul(z_vec_, p_vec_);
  real_t const alpha{safe_divide(gamma_, dot_product(p_vec_, z_vec_))};
  x_vec += alpha * p_vec_;
  r_vec_ -= alpha * z_vec_;

  // ----------------------
  // 𝛾̅ ← 𝛾,
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝒛 ← 𝓟𝒓,
  //   𝛾 ← <𝒓⋅𝒛>,
  // 𝗲𝗹𝘀𝗲:
  //   𝛾 ← <𝒓⋅𝒓>.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  real_t const gamma_bar{gamma_};
  if (pre_op != nullptr) {
    pre_op->mul(z_vec_, r_vec_);
    gamma_ = dot_product(r_vec_, z_vec_);
  } else {
    gamma_ = dot_product(r_vec_, r_vec_);
  }

  // ----------------------
  // 𝛽 ← 𝛾/𝛾̅,
  // 𝒑 ← (𝓟 ≠ 𝗻𝗼𝗻𝗲 ? 𝒛 : 𝒓) + 𝛽⋅𝒑.
  // ----------------------
  real_t const beta = safe_divide(gamma_, gamma_bar);
  p_vec_ <<= (pre_op != nullptr ? z_vec_ : r_vec_) + beta * p_vec_;

  return (pre_op != nullptr) ? norm_2(r_vec_) : std::sqrt(gamma_);

} // CgSolver::iterate

} // namespace Storm
