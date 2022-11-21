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

namespace Storm {

/// @brief The CG (Conjugate Gradients) linear self-adjoint definite operator
/// equation solver.
///
/// CG may be applied to the consistent singular problems, it converges
/// towards..
///
/// References:
/// @verbatim
/// [1] Hestenes, Magnus R. and Eduard Stiefel.
///     “Methods of conjugate gradients for solving linear systems.”
///     Journal of research of the National Bureau of Standards 49 (1952):
/// 409-435.
/// @endverbatim
template<legacy_vector_like Vector>
class CgSolver final : public IterativeSolver<Vector> {
private:

  real_t gamma_;
  Vector p_vec_, r_vec_, z_vec_;

  real_t init(const Vector& x_vec, const Vector& b_vec,
              const Operator<Vector>& lin_op,
              const Preconditioner<Vector>* pre_op) override {
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
  }

  real_t iterate(Vector& x_vec, const Vector& b_vec,
                 const Operator<Vector>& lin_op,
                 const Preconditioner<Vector>* pre_op) override {
    // Iterate:
    // ----------------------
    // 𝒛 ← 𝓐𝒑,
    // 𝛼 ← 𝛾/<𝒑⋅𝒛>,
    // 𝒙 ← 𝒙 + 𝛼⋅𝒑,
    // 𝒓 ← 𝒓 - 𝛼⋅𝒛.
    // ----------------------
    lin_op.mul(z_vec_, p_vec_);
    const real_t alpha = safe_divide(gamma_, dot_product(p_vec_, z_vec_));
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
    const real_t gamma_bar = gamma_;
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
    const real_t beta = safe_divide(gamma_, gamma_bar);
    p_vec_ <<= (pre_op != nullptr ? z_vec_ : r_vec_) + beta * p_vec_;

    return (pre_op != nullptr) ? norm_2(r_vec_) : sqrt(gamma_);
  }

}; // class CgSolver

} // namespace Storm
