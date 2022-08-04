/**
 * Copyright (C) 2022 Oleg Butakov
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#pragma once

#include <utility>

#include <Storm/Base.hpp>

#include <Storm/Utils/Math.hpp>

#include <Storm/Blass/Vector.hpp>

#include <Storm/Solvers/Solver.hpp>

namespace Storm {

/**
 * @brief The CGS (Conjugate Gradients Squared) linear operator equation solver.
 *
 * CGS, like the other BiCG type solvers, requires two operator multiplications
 * per iteration.
 *
 * @warning CGS convergence behavior may be very erratic.
 *
 * References:
 * @verbatim
 * [1] Sonneveld, Peter.
 *     “CGS, A Fast Lanczos-Type Solver for Nonsymmetric Linear systems.”
 *     SIAM J. Sci. Stat. Comput., 10:36-52, 1989.
 * @endverbatim
 */
template<VectorLike Vector>
class CgsSolver final : public IterativeSolver<Vector> {
private:

  real_t rho_;
  Vector p_vec_, q_vec_, r_vec_, r_tilde_vec_, u_vec_, v_vec_;

  real_t init(const Vector& x_vec, const Vector& b_vec,
              const Operator<Vector>& lin_op,
              const Preconditioner<Vector>* pre_op) override;

  real_t iterate(Vector& x_vec, const Vector& b_vec,
                 const Operator<Vector>& lin_op,
                 const Preconditioner<Vector>* pre_op) override;

}; // class CgsSolver

template<VectorLike Vector>
real_t CgsSolver<Vector>::init(const Vector& x_vec, const Vector& b_vec,
                               const Operator<Vector>& lin_op,
                               const Preconditioner<Vector>* pre_op) {
  const bool left_pre{(pre_op != nullptr) &&
                      (this->pre_side == PreconditionerSide::Left)};

  p_vec_.assign(x_vec, false);
  q_vec_.assign(x_vec, false);
  r_vec_.assign(x_vec, false);
  r_tilde_vec_.assign(x_vec, false);
  u_vec_.assign(x_vec, false);
  v_vec_.assign(x_vec, false);

  // Initialize:
  // ----------------------
  // 𝒓 ← 𝒃 - 𝓐𝒙,
  // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
  //   𝒖 ← 𝒓,
  //   𝒓 ← 𝓟𝒖,
  // 𝗲𝗻𝗱 𝗶𝗳
  // 𝒓̃ ← 𝒓,
  // 𝜌 ← <𝒓̃⋅𝒓>.
  // ----------------------
  lin_op.Residual(r_vec_, b_vec, x_vec);
  if (left_pre) {
    std::swap(u_vec_, r_vec_);
    pre_op->mul(r_vec_, u_vec_);
  }
  r_tilde_vec_ <<= r_vec_;
  rho_ = dot_product(r_tilde_vec_, r_vec_);

  return math::sqrt(rho_);

} // CgsSolver::init

template<VectorLike Vector>
real_t CgsSolver<Vector>::iterate(Vector& x_vec, const Vector& b_vec,
                                  const Operator<Vector>& lin_op,
                                  const Preconditioner<Vector>* pre_op) {
  const bool left_pre{(pre_op != nullptr) &&
                      (this->pre_side == PreconditionerSide::Left)};
  const bool right_pre{(pre_op != nullptr) &&
                       (this->pre_side == PreconditionerSide::Right)};

  // Continue the iterations:
  // ----------------------
  // 𝗶𝗳 𝘍𝘪𝘳𝘴𝘵𝘐𝘵𝘦𝘳𝘢𝘵𝘪𝘰𝘯:
  //   𝒖 ← 𝒓,
  //   𝒑 ← 𝒖.
  // 𝗲𝗹𝘀𝗲:
  //   𝜌̅ ← 𝜌,
  //   𝜌 ← <𝒓̃⋅𝒓>,
  //   𝛽 ← 𝜌/𝜌̅,
  //   𝒖 ← 𝒓 + 𝛽⋅𝒒,
  //   𝒑 ← 𝒖 + 𝛽⋅(𝒒 + 𝛽⋅𝒑).
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  const bool first_iteration{this->iteration == 0};
  if (first_iteration) {
    u_vec_ <<= r_vec_;
    p_vec_ <<= u_vec_;
  } else {
    const real_t rho_bar{
        std::exchange(rho_, dot_product(r_tilde_vec_, r_vec_))};
    const real_t beta{math::safe_divide(rho_, rho_bar)};
    u_vec_ <<= r_vec_ + beta * q_vec_;
    p_vec_ <<= u_vec_ + beta * (q_vec_ + beta * p_vec_);
  }

  // ----------------------
  // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
  //   𝒗 ← 𝓟(𝒒 ← 𝓐𝒑),
  // 𝗲𝗹𝘀𝗲 𝗶𝗳 𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦:
  //   𝒗 ← 𝓐(𝒒 ← 𝓟𝒑),
  // 𝗲𝗹𝘀𝗲:
  //   𝒗 ← 𝓐𝒑,
  // 𝗲𝗻𝗱 𝗶𝗳
  // 𝛼 ← 𝜌/<𝒓̃⋅𝒗>,
  // 𝒒 ← 𝒖 - 𝛼⋅𝒗,
  // 𝒗 ← 𝒖 + 𝒒.
  // ----------------------
  if (left_pre) {
    pre_op->mul(v_vec_, q_vec_, lin_op, p_vec_);
  } else if (right_pre) {
    lin_op.mul(v_vec_, q_vec_, *pre_op, p_vec_);
  } else {
    lin_op.mul(v_vec_, p_vec_);
  }
  const real_t alpha{
      math::safe_divide(rho_, dot_product(r_tilde_vec_, v_vec_))};
  q_vec_ <<= u_vec_ - alpha * v_vec_;
  v_vec_ <<= u_vec_ + q_vec_;

  // Update the solution and the residual:
  // ----------------------
  // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
  //   𝒙 ← 𝒙 + 𝛼⋅𝒗,
  //   𝒗 ← 𝓟(𝒖 ← 𝓐𝒗),
  //   𝒓 ← 𝒓 - 𝛼⋅𝒗.
  // 𝗲𝗹𝘀𝗲 𝗶𝗳 𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦:
  //   𝒗 ← 𝓐(𝒖 ← 𝓟𝒗),
  //   𝒙 ← 𝒙 + 𝛼⋅𝒖,
  //   𝒓 ← 𝒓 - 𝛼⋅𝒗.
  // 𝗲𝗹𝘀𝗲:
  //   𝒖 ← 𝓐𝒗,
  //   𝒙 ← 𝒙 + 𝛼⋅𝒗,
  //   𝒓 ← 𝒓 - 𝛼⋅𝒖.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  if (left_pre) {
    x_vec += alpha * v_vec_;
    pre_op->mul(v_vec_, u_vec_, lin_op, v_vec_);
    r_vec_ -= alpha * v_vec_;
  } else if (right_pre) {
    lin_op.mul(v_vec_, u_vec_, *pre_op, v_vec_);
    x_vec += alpha * u_vec_;
    r_vec_ -= alpha * v_vec_;
  } else {
    lin_op.mul(u_vec_, v_vec_);
    x_vec += alpha * v_vec_;
    r_vec_ -= alpha * u_vec_;
  }

  return norm_2(r_vec_);

} // CgsSolver::iterate

} // namespace Storm
