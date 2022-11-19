// Copyright (C) 2020 - 2023 Oleg Butakov
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

#include <Storm/Utils/Math.hpp>

#include <Storm/Bittern/Matrix.hpp>
#include <Storm/Bittern/MatrixDense.hpp>

#include <Storm/Solvers/Solver.hpp>

#include <utility>

namespace Storm {

/// @brief The BiCGStab (Biconjugate Gradients Stabilized) linear operator
/// equation solver.
///
/// BiCGStab, like the other BiCG type solvers, requires two operator
/// multiplications per iteration. BiCGStab typically converges much smoother,
/// than CGS.
///
/// References:
/// @verbatim
/// [1] Henk A. van der Vorst.
///     “Bi-CGSTAB: A Fast and Smoothly Converging Variant of Bi-CG for the
///      Solution of Nonsymmetric Linear Systems.”
///     SIAM J. Sci. Comput. 13 (1992): 631-644.
/// @endverbatim
template<legacy_vector_like Vector>
class BiCgStabSolver final : public IterativeSolver<Vector> {
private:

  real_t alpha_, rho_, omega_;
  Vector p_vec_, r_vec_, r_tilde_vec_, t_vec_, v_vec_, z_vec_;

  real_t init(const Vector& x_vec, const Vector& b_vec,
              const Operator<Vector>& lin_op,
              const Preconditioner<Vector>* pre_op) override {
    const bool left_pre =
        (pre_op != nullptr) && (this->pre_side == PreconditionerSide::Left);

    p_vec_.assign(x_vec, false);
    r_vec_.assign(x_vec, false);
    r_tilde_vec_.assign(x_vec, false);
    t_vec_.assign(x_vec, false);
    v_vec_.assign(x_vec, false);
    if (pre_op != nullptr) { z_vec_.assign(x_vec, false); }

    // Initialize:
    // ----------------------
    // 𝒓 ← 𝒃 - 𝓐𝒙,
    // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
    //   𝒛 ← 𝒓,
    //   𝒓 ← 𝓟𝒛,
    // 𝗲𝗻𝗱 𝗶𝗳
    // 𝒓̃ ← 𝒓,
    // 𝜌 ← <𝒓̃⋅𝒓>.
    // ----------------------
    lin_op.Residual(r_vec_, b_vec, x_vec);
    if (left_pre) {
      std::swap(z_vec_, r_vec_);
      pre_op->mul(r_vec_, z_vec_);
    }
    r_tilde_vec_ <<= r_vec_;
    rho_ = dot_product(r_tilde_vec_, r_vec_);

    return sqrt(rho_);
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
    //   𝒑 ← 𝒓.
    // 𝗲𝗹𝘀𝗲:
    //   𝜌̅ ← 𝜌,
    //   𝜌 ← <𝒓̃⋅𝒓>,
    //   𝛽 ← (𝜌/𝜌̅)⋅(𝛼/𝜔),
    //   𝒑 ← 𝒓 + 𝛽⋅(𝒑 - 𝜔⋅𝒗).
    // 𝗲𝗻𝗱 𝗶𝗳
    // ----------------------
    const bool first_iteration = this->iteration == 0;
    if (first_iteration) {
      p_vec_ <<= r_vec_;
    } else {
      const real_t rho_bar =
          std::exchange(rho_, dot_product(r_tilde_vec_, r_vec_));
      const real_t beta{safe_divide(alpha_ * rho_, omega_ * rho_bar)};
      p_vec_ <<= r_vec_ + beta * (p_vec_ - omega_ * v_vec_);
    }

    // Update the solution and the residual:
    // ----------------------
    // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
    //   𝒗 ← 𝓟(𝒛 ← 𝓐𝒑),
    // 𝗲𝗹𝘀𝗲 𝗶𝗳 𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦:
    //   𝒗 ← 𝓐(𝒛 ← 𝓟𝒑),
    // 𝗲𝗹𝘀𝗲:
    //   𝒗 ← 𝓐𝒑,
    // 𝗲𝗻𝗱 𝗶𝗳
    // 𝛼 ← 𝜌/<𝒓̃⋅𝒗>,
    // 𝒙 ← 𝒙 + 𝛼⋅(𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦 ? 𝒛 : 𝒑),
    // 𝒓 ← 𝒓 - 𝛼⋅𝒗.
    // ----------------------
    if (left_pre) {
      pre_op->mul(v_vec_, z_vec_, lin_op, p_vec_);
    } else if (right_pre) {
      lin_op.mul(v_vec_, z_vec_, *pre_op, p_vec_);
    } else {
      lin_op.mul(v_vec_, p_vec_);
    }
    alpha_ = safe_divide(rho_, dot_product(r_tilde_vec_, v_vec_));
    x_vec += alpha_ * (right_pre ? z_vec_ : p_vec_);
    r_vec_ -= alpha_ * v_vec_;

    // Update the solution and the residual again:
    // ----------------------
    // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
    //   𝒕 ← 𝓟(𝒛 ← 𝓐𝒓),
    // 𝗲𝗹𝘀𝗲 𝗶𝗳 𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦:
    //   𝒕 ← 𝓐(𝒛 ← 𝓟𝒓),
    // 𝗲𝗹𝘀𝗲:
    //   𝒕 ← 𝓐𝒓,
    // 𝗲𝗻𝗱 𝗶𝗳
    // 𝜔 ← <𝒕⋅𝒓>/<𝒕⋅𝒕>,
    // 𝒙 ← 𝒙 + 𝜔⋅(𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦 ? 𝒛 : 𝒓),
    // 𝒓 ← 𝒓 - 𝜔⋅𝒕.
    // ----------------------
    if (left_pre) {
      pre_op->mul(t_vec_, z_vec_, lin_op, r_vec_);
    } else if (right_pre) {
      lin_op.mul(t_vec_, z_vec_, *pre_op, r_vec_);
    } else {
      lin_op.mul(t_vec_, r_vec_);
    }
    omega_ =
        safe_divide(dot_product(t_vec_, r_vec_), dot_product(t_vec_, t_vec_));
    x_vec += omega_ * (right_pre ? z_vec_ : r_vec_);
    r_vec_ -= omega_ * t_vec_;

    return norm_2(r_vec_);
  }

}; // class BiCgStabSolver

// -----------------------------------------------------------------------------

/// @brief The BiCGStab(l) (Biconjugate Gradients Stabilized)
///   linear operator equation solver.
///
/// BiCGStab(l), like the other BiCG type solvers, requires two operator
/// multiplications per iteration.
///
/// References:
/// @verbatim
/// [1] Gerard L. G. Sleijpen and Diederik R. Fokkema.
///     “BiCGStab(l) for Linear Equations involving Unsymmetric Matrices with
///     Complex Spectrum.”
///     Electronic Transactions on Numerical Analysis 1 (1993): 11-32.
/// @endverbatim
template<legacy_vector_like Vector>
class BiCgStabLSolver final : public InnerOuterIterativeSolver<Vector> {
private:

  real_t alpha_, rho_, omega_;
  DenseVector<real_t> gamma_, gamma_bar_, gamma_bbar_, sigma_;
  DenseMatrix<real_t> tau_;
  Vector r_tilde_vec_, z_vec_;
  std::vector<Vector> r_vecs_, u_vecs_;

  real_t outer_init(const Vector& x_vec, const Vector& b_vec,
                    const Operator<Vector>& lin_op,
                    const Preconditioner<Vector>* pre_op) override {
    const size_t l = this->num_inner_iterations;

    gamma_.assign(l + 1);
    gamma_bar_.assign(l + 1);
    gamma_bbar_.assign(l + 1);
    sigma_.assign(l + 1);
    tau_.assign(l + 1, l + 1);

    r_tilde_vec_.assign(x_vec, false);
    if (pre_op != nullptr) { z_vec_.assign(x_vec, false); }

    r_vecs_.resize(l + 1);
    u_vecs_.resize(l + 1);
    for (Vector& r_vec : r_vecs_) {
      r_vec.assign(x_vec, false);
    }
    for (Vector& u_vec : u_vecs_) {
      u_vec.assign(x_vec, false);
    }

    // Initialize:
    // ----------------------
    // 𝒖₀ ← {𝟢}ᵀ,
    // 𝒓₀ ← 𝒃 - 𝓐𝒙,
    // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
    //   𝒛 ← 𝒓₀,
    //   𝒓₀ ← 𝓟𝒛,
    // 𝗲𝗻𝗱 𝗶𝗳
    // 𝒓̃ ← 𝒓₀,
    // 𝜌 ← <𝒓̃⋅𝒓₀>.
    // ----------------------
    fill_with(u_vecs_[0], 0.0);
    lin_op.Residual(r_vecs_[0], b_vec, x_vec);
    if (pre_op != nullptr) {
      std::swap(z_vec_, r_vecs_[0]);
      pre_op->mul(r_vecs_[0], z_vec_);
    }
    r_tilde_vec_ <<= r_vecs_[0];
    rho_ = dot_product(r_tilde_vec_, r_vecs_[0]);

    return sqrt(rho_);
  }

  real_t inner_iterate(Vector& x_vec, const Vector& b_vec,
                       const Operator<Vector>& lin_op,
                       const Preconditioner<Vector>* pre_op) override {
    const size_t l = this->num_inner_iterations;
    const size_t j = this->inner_iteration;

    // BiCG part:
    // ----------------------
    // 𝗶𝗳 𝘍𝘪𝘳𝘴𝘵𝘐𝘵𝘦𝘳𝘢𝘵𝘪𝘰𝘯:
    //   𝒖₀ ← 𝒓₀,
    // 𝗲𝗹𝘀𝗲:
    //   𝜌̅ ← 𝜌,
    //   𝜌 ← <𝒓̃⋅𝒓ⱼ>,
    //   𝛽 ← 𝛼⋅𝜌/𝜌̅,
    //   𝗳𝗼𝗿 𝑖 = 𝟢, 𝑗 𝗱𝗼:
    //     𝒖ᵢ ← 𝒓ᵢ - 𝛽⋅𝒖ᵢ,
    //   𝗲𝗻𝗱 𝗳𝗼𝗿
    // 𝗲𝗻𝗱 𝗶𝗳
    // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
    //   𝒖ⱼ₊₁ ← 𝓟(𝒛 ← 𝓐𝒖ⱼ),
    // 𝗲𝗹𝘀𝗲:
    //   𝒖ⱼ₊₁ ← 𝓐𝒖ⱼ,
    // 𝗲𝗻𝗱 𝗶𝗳
    // 𝛼 ← 𝜌/<𝒓̃⋅𝒖ⱼ₊₁>,
    // 𝗳𝗼𝗿 𝑖 = 𝟢, 𝑗 𝗱𝗼:
    //   𝒓ᵢ ← 𝒓ᵢ - 𝛼⋅𝒖ᵢ₊₁.
    // 𝗲𝗻𝗱 𝗳𝗼𝗿
    // ----------------------
    const bool first_iteration = this->iteration == 0;
    if (first_iteration) {
      u_vecs_[0] <<= r_vecs_[0];
    } else {
      const real_t rho_bar =
          std::exchange(rho_, dot_product(r_tilde_vec_, r_vecs_[j]));
      const real_t beta = safe_divide(alpha_ * rho_, rho_bar);
      for (size_t i = 0; i <= j; ++i) {
        u_vecs_[i] <<= r_vecs_[i] - beta * u_vecs_[i];
      }
    }
    if (pre_op != nullptr) {
      pre_op->mul(u_vecs_[j + 1], z_vec_, lin_op, u_vecs_[j]);
    } else {
      lin_op.mul(u_vecs_[j + 1], u_vecs_[j]);
    }
    alpha_ = safe_divide(rho_, dot_product(r_tilde_vec_, u_vecs_[j + 1]));
    for (size_t i = 0; i <= j; ++i) {
      r_vecs_[i] -= alpha_ * u_vecs_[i + 1];
    }

    // Update the solution and the residual:
    // ----------------------
    // 𝒙 ← 𝒙 + 𝛼⋅𝒖₀,
    // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
    //   𝒓ⱼ₊₁ ← 𝓟(𝒛 ← 𝓐𝒓ⱼ).
    // 𝗲𝗹𝘀𝗲:
    //   𝒓ⱼ₊₁ ← 𝓐𝒓ⱼ.
    // 𝗲𝗻𝗱 𝗶𝗳
    // ----------------------
    x_vec += alpha_ * u_vecs_[0];
    if (pre_op != nullptr) {
      pre_op->mul(r_vecs_[j + 1], z_vec_, lin_op, r_vecs_[j]);
    } else {
      lin_op.mul(r_vecs_[j + 1], r_vecs_[j]);
    }

    if (j == l - 1) {
      // Minimal residual part:
      // ----------------------
      // 𝗳𝗼𝗿 𝑗 = 𝟣, 𝑙 𝗱𝗼:
      //   𝗳𝗼𝗿 𝑖 = 𝟣, 𝑗 - 𝟣 𝗱𝗼:
      //     𝜏ᵢⱼ ← <𝒓ᵢ⋅𝒓ⱼ>/𝜎ᵢ,
      //     𝒓ⱼ ← 𝒓ⱼ - 𝜏ᵢⱼ⋅𝒓ᵢ,
      //   𝗲𝗻𝗱 𝗳𝗼𝗿
      //   𝜎ⱼ ← <𝒓ⱼ⋅𝒓ⱼ>,
      //   𝛾̅ⱼ ← <𝒓₀⋅𝒓ⱼ>/𝜎ⱼ,
      // 𝗲𝗻𝗱 𝗳𝗼𝗿
      // ----------------------
      for (size_t j = 1; j <= l; ++j) {
        for (size_t i = 1; i < j; ++i) {
          tau_(i, j) =
              safe_divide(dot_product(r_vecs_[i], r_vecs_[j]), sigma_(i));
          r_vecs_[j] -= tau_(i, j) * r_vecs_[i];
        }
        sigma_(j) = dot_product(r_vecs_[j], r_vecs_[j]);
        gamma_bar_(j) =
            safe_divide(dot_product(r_vecs_[0], r_vecs_[j]), sigma_(j));
      }

      // ----------------------
      // 𝜔 ← 𝛾ₗ ← 𝛾̅ₗ, 𝜌 ← -𝜔⋅𝜌,
      // 𝗳𝗼𝗿 𝑗 = 𝑙 - 𝟣, 𝟣, -𝟣 𝗱𝗼:
      //   𝛾ⱼ ← 𝛾̅ⱼ,
      //   𝗳𝗼𝗿 𝑖 = 𝑗 + 𝟣, 𝑙 𝗱𝗼:
      //     𝛾ⱼ ← 𝛾ⱼ - 𝜏ⱼᵢ⋅𝛾ᵢ,
      //   𝗲𝗻𝗱 𝗳𝗼𝗿
      // 𝗲𝗻𝗱 𝗳𝗼𝗿
      // 𝗳𝗼𝗿 𝑗 = 𝟣, 𝑙 - 𝟣 𝗱𝗼:
      //   𝛾̿ⱼ ← 𝛾ⱼ₊₁,
      //   𝗳𝗼𝗿 𝑖 = 𝑗 + 𝟣, 𝑙 - 𝟣 𝗱𝗼:
      //     𝛾̿ⱼ ← 𝛾̿ⱼ + 𝜏ⱼᵢ⋅𝛾ᵢ₊₁.
      //   𝗲𝗻𝗱 𝗳𝗼𝗿
      // 𝗲𝗻𝗱 𝗳𝗼𝗿
      // ----------------------
      omega_ = gamma_(l) = gamma_bar_(l), rho_ *= -omega_;
      for (size_t j = l - 1; j != 0; --j) {
        gamma_(j) = gamma_bar_(j);
        for (size_t i = j + 1; i <= l; ++i) {
          gamma_(j) -= tau_(j, i) * gamma_(i);
        }
      }
      for (size_t j = 1; j < l; ++j) {
        gamma_bbar_(j) = gamma_(j + 1);
        for (size_t i = j + 1; i < l; ++i) {
          gamma_bbar_(j) += tau_(j, i) * gamma_(i + 1);
        }
      }

      // Update the solution and the residual again:
      // ----------------------
      // 𝒙 ← 𝒙 + 𝛾₁⋅𝒓₀,
      // 𝒓₀ ← 𝒓₀ - 𝛾̅ₗ⋅𝒓ₗ,
      // 𝒖₀ ← 𝒖₀ - 𝛾ₗ⋅𝒖ₗ,
      // 𝗳𝗼𝗿 𝑗 = 𝟣, 𝑙 - 𝟣 𝗱𝗼:
      //   𝒙 ← 𝒙 + 𝛾̿ⱼ⋅𝒓ⱼ,
      //   𝒓₀ ← 𝒓₀ - 𝛾̅ⱼ⋅𝒓ⱼ,
      //   𝒖₀ ← 𝒖₀ - 𝛾ⱼ⋅𝒖ⱼ.
      // 𝗲𝗻𝗱 𝗳𝗼𝗿
      // ----------------------
      x_vec += gamma_(1) * r_vecs_[0];
      r_vecs_[0] -= gamma_bar_(l) * r_vecs_[l];
      u_vecs_[0] -= gamma_(l) * u_vecs_[l];
      for (size_t j = 1; j < l; ++j) {
        x_vec += gamma_bbar_(j) * r_vecs_[j];
        r_vecs_[0] -= gamma_bar_(j) * r_vecs_[j];
        u_vecs_[0] -= gamma_(j) * u_vecs_[j];
      }
    }

    return norm_2(r_vecs_[0]);
  }

public:

  BiCgStabLSolver() {
    this->num_inner_iterations = 2;
  }

}; // class BiCgStabLSolver

} // namespace Storm
