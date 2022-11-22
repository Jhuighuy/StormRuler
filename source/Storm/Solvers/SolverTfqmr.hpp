// Copyright (C) 2020-2023 Oleg Butakov
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR Allocator PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT
// SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
// OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

#pragma once

#include <Storm/Base.hpp>

#include <Storm/Bittern/Math.hpp>
#include <Storm/Bittern/Matrix.hpp>

#include <Storm/Solvers/Solver.hpp>

#include <utility>

namespace Storm {

/// @brief Parametrized TFQMR solver.
template<legacy_vector_like Vector, bool L1>
class BaseTfqmrSolver : public IterativeSolver<Vector> {
private:

  real_t rho_, tau_;
  Vector d_vec_, r_tilde_vec_, u_vec_, v_vec_, y_vec_, s_vec_, z_vec_;

  real_t init(const Vector& x_vec, const Vector& b_vec,
              const Operator<Vector>& lin_op,
              const Preconditioner<Vector>* pre_op) override {
    const bool left_pre =
        (pre_op != nullptr) && (this->pre_side == PreconditionerSide::Left);

    d_vec_.assign(x_vec, false);
    r_tilde_vec_.assign(x_vec, false);
    u_vec_.assign(x_vec, false);
    v_vec_.assign(x_vec, false);
    y_vec_.assign(x_vec, false);
    s_vec_.assign(x_vec, false);
    if (pre_op != nullptr) z_vec_.assign(x_vec, false);

    // Initialize:
    // ----------------------
    // 𝗶𝗳 𝘓₁:
    //   𝒅 ← 𝒙,
    // 𝗲𝗹𝘀𝗲:
    //   𝒅 ← {𝟢}ᵀ,
    // 𝗲𝗻𝗱 𝗶𝗳
    // 𝒚 ← 𝒃 - 𝓐𝒙,
    // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
    //   𝒛 ← 𝒚,
    //   𝒚 ← 𝓟𝒛,
    // 𝗲𝗻𝗱 𝗶𝗳
    // 𝒖 ← 𝒚,
    // 𝒓̃ ← 𝒖,
    // 𝜌 ← <𝒓̃⋅𝒓>, 𝜏 ← 𝜌¹ᐟ².
    // ----------------------
    if constexpr (L1) {
      d_vec_ <<= x_vec;
    } else {
      fill_with(d_vec_, 0.0);
    }
    lin_op.Residual(y_vec_, b_vec, x_vec);
    if (left_pre) {
      std::swap(z_vec_, y_vec_);
      pre_op->mul(y_vec_, z_vec_);
    }
    u_vec_ <<= y_vec_;
    r_tilde_vec_ <<= u_vec_;
    rho_ = dot_product(r_tilde_vec_, u_vec_), tau_ = sqrt(rho_);

    return tau_;
  }

  real_t iterate(Vector& x_vec, const Vector& b_vec,
                 const Operator<Vector>& lin_op,
                 const Preconditioner<Vector>* pre_op) override {
    const bool left_pre =
        (pre_op != nullptr) && (this->pre_side == PreconditionerSide::Left);
    const bool right_pre =
        (pre_op != nullptr) && (this->pre_side == PreconditionerSide::Right);

    // Continue the iterations:
    // ----------------------
    // 𝗶𝗳 𝘍𝘪𝘳𝘴𝘵𝘐𝘵𝘦𝘳𝘢𝘵𝘪𝘰𝘯:
    //   𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
    //     𝒔 ← 𝓟(𝒛 ← 𝓐𝒚),
    //   𝗲𝗹𝘀𝗲 𝗶𝗳 𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦:
    //     𝒔 ← 𝓐(𝒛 ← 𝓟𝒚),
    //   𝗲𝗹𝘀𝗲:
    //     𝒔 ← 𝓐𝒚.
    //   𝗲𝗻𝗱 𝗶𝗳
    //   𝒗 ← 𝒔,
    // 𝗲𝗹𝘀𝗲:
    //   𝜌̅ ← 𝜌,
    //   𝜌 ← <𝒓̃⋅𝒖>,
    //   𝛽 ← 𝜌/𝜌̅,
    //   𝒗 ← 𝒔 + 𝛽⋅𝒗,
    //   𝒚 ← 𝒖 + 𝛽⋅𝒚,
    //   𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
    //     𝒔 ← 𝓟(𝒛 ← 𝓐𝒚),
    //   𝗲𝗹𝘀𝗲 𝗶𝗳 𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦:
    //     𝒔 ← 𝓐(𝒛 ← 𝓟𝒚),
    //   𝗲𝗹𝘀𝗲:
    //     𝒔 ← 𝓐𝒚,
    //   𝗲𝗻𝗱 𝗶𝗳
    //   𝒗 ← 𝒔 + 𝛽⋅𝒗.
    // 𝗲𝗻𝗱 𝗶𝗳
    // ----------------------
    const bool first_iteration = this->iteration == 0;
    if (first_iteration) {
      if (left_pre) pre_op->mul(s_vec_, z_vec_, lin_op, y_vec_);
      else if (right_pre) lin_op.mul(s_vec_, z_vec_, *pre_op, y_vec_);
      else lin_op.mul(s_vec_, y_vec_);
      v_vec_ <<= s_vec_;
    } else {
      const real_t rho_bar =
          std::exchange(rho_, dot_product(r_tilde_vec_, u_vec_));
      const real_t beta = safe_divide(rho_, rho_bar);
      v_vec_ <<= s_vec_ + beta * v_vec_;
      y_vec_ <<= u_vec_ + beta * y_vec_;
      if (left_pre) pre_op->mul(s_vec_, z_vec_, lin_op, y_vec_);
      else if (right_pre) lin_op.mul(s_vec_, z_vec_, *pre_op, y_vec_);
      else lin_op.mul(s_vec_, y_vec_);
      v_vec_ <<= s_vec_ + beta * v_vec_;
    }

    // Update the solution:
    // ----------------------
    // 𝛼 ← 𝜌/<𝒓̃⋅𝒗>,
    // 𝗳𝗼𝗿 𝑚 = 𝟢, 𝟣 𝗱𝗼:
    //   𝒖 ← 𝒖 - 𝛼⋅𝒔,
    //   𝒅 ← 𝒅 + 𝛼⋅(𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦 ? 𝒛 : 𝒚),
    //   𝜔 ← ‖𝒖‖,
    //   𝗶𝗳 𝘓₁:
    //     𝗶𝗳 𝜔 < 𝜏:
    //       𝜏 ← 𝜔, 𝒙 ← 𝒅,
    //     𝗲𝗻𝗱 𝗶𝗳
    //   𝗲𝗹𝘀𝗲:
    //     𝑐𝑠, 𝑠𝑛 ← 𝘚𝘺𝘮𝘖𝘳𝘵𝘩𝘰(𝜏, 𝜔),
    //     𝜏 ← 𝑐𝑠⋅𝜔,
    //     𝒙 ← 𝒙 + 𝑐𝑠²⋅𝒅,
    //     𝒅 ← 𝑠𝑛²⋅𝒅,
    //   𝗲𝗻𝗱 𝗶𝗳
    //   𝗶𝗳 𝑚 = 𝟢:
    //     𝒚 ← 𝒚 - 𝛼⋅𝒗,
    //     𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
    //       𝒔 ← 𝓟(𝒛 ← 𝓐𝒚).
    //     𝗲𝗹𝘀𝗲 𝗶𝗳 𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦:
    //       𝒔 ← 𝓐(𝒛 ← 𝓟𝒚).
    //     𝗲𝗹𝘀𝗲:
    //       𝒔 ← 𝓐𝒚.
    //     𝗲𝗻𝗱 𝗶𝗳
    //   𝗲𝗻𝗱 𝗶𝗳
    // 𝗲𝗻𝗱 𝗳𝗼𝗿
    // ----------------------
    const real_t alpha = safe_divide(rho_, dot_product(r_tilde_vec_, v_vec_));
    for (size_t m = 0; m <= 1; ++m) {
      u_vec_ -= alpha * s_vec_;
      d_vec_ += alpha * (right_pre ? z_vec_ : y_vec_);
      const real_t omega = norm_2(u_vec_);
      if constexpr (L1) {
        if (omega < tau_) tau_ = omega, x_vec <<= d_vec_;
      } else {
        const auto [cs, sn, rr] = sym_ortho(tau_, omega);
        tau_ = omega * cs;
        x_vec += std::pow(cs, 2) * d_vec_;
        d_vec_ *= std::pow(sn, 2);
      }
      if (m == 0) {
        y_vec_ -= alpha * v_vec_;
        if (left_pre) pre_op->mul(s_vec_, z_vec_, lin_op, y_vec_);
        else if (right_pre) lin_op.mul(s_vec_, z_vec_, *pre_op, y_vec_);
        else lin_op.mul(s_vec_, y_vec_);
      }
    }

    // Compute the residual norm
    // (or it's upper bound estimate in the ℒ₂ case):
    // ----------------------
    // 𝜏̃ ← 𝜏,
    // 𝗶𝗳 𝗻𝗼𝘁 𝘓₁:
    //   𝜏̃ ← 𝜏⋅(𝟤𝑘 + 𝟥)¹ᐟ².
    // 𝗲𝗻𝗱 𝗶𝗳
    // ----------------------
    real_t tau_tilde = tau_;
    if constexpr (!L1) {
      const size_t k{this->iteration};
      tau_tilde *= std::sqrt(2.0 * k + 3.0);
    }

    return tau_tilde;
  }

protected:

  BaseTfqmrSolver() = default;

}; // class BaseTfqmrSolver

// -----------------------------------------------------------------------------

/// @brief The TFQMR (Transpose-Free Quasi-Minimal Residual) linear operator
/// equation solver.
///
/// TFQMR, like the other BiCG type methods, normally requires two
/// operator-vector products per iteration. But, unlike the other BiCG type
/// methods, TFQMR does not implicitly contain the residual norm estimate,
/// only the rough upper bound is avariable, so at the latter iterations an
/// extra operator-vector product per iteration may be required for the explicit
/// residual estimation.
///
/// TFQMR typically converges much smoother, than CGS and BiCGStab.
/// @todo Breakdowns?
///
/// References:
/// @verbatim
/// [1] Freund, Roland W.
///     “A Transpose-Free Quasi-Minimal Residual Algorithm for Non-Hermitian
///      Linear Systems.”
///     SIAM J. Sci. Comput. 14 (1993): 470-482.
/// [2] Freund, Roland W.
///     “Transpose-Free Quasi-Minimal Residual Methods for Non-Hermitian Linear
///      Systems.” (1994).
/// @endverbatim
template<legacy_vector_like Vector>
class TfqmrSolver final : public BaseTfqmrSolver<Vector, false> {};

/// @brief The TFQMR1 (Transpose-Free 1-norm Quasi-Minimal Residual) linear
/// operator equation solver.
///
/// TFQMR1, like the other BiCG type solvers, requires two operator-vector
/// products per iteration. Unlike TFQMR, TFQMR1 implicitly contains the
/// residual norm estimate, so no extra operator-vector products are required.
///
/// TFQMR1 typically converges much smoother, than CGS and BiCGStab and is
/// slightly faster than TFQMR.
/// @todo Breakdowns?
///
/// References:
/// @verbatim
/// [1] H.M Bücker,
///     “A Transpose-Free 1-norm Quasi-Minimal Residual Algorithm
///      for Non-Hermitian Linear Systems.“, FZJ-ZAM-IB-9706.
/// @endverbatim
template<legacy_vector_like Vector>
class Tfqmr1Solver final : public BaseTfqmrSolver<Vector, true> {};

} // namespace Storm
