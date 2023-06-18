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
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <Storm/Base.hpp>

#include <Storm/Bittern/Math.hpp>
#include <Storm/Bittern/Matrix.hpp>

#include <Storm/Solvers/Solver.hpp>

namespace Storm {

// -----------------------------------------------------------------------------

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

  real_t _gamma;
  Vector _p_vec, _r_vec, _z_vec;

  real_t init(const Vector& x_vec, const Vector& b_vec,
              const Operator<Vector>& lin_op,
              const Preconditioner<Vector>* pre_op) override {
    _p_vec.assign(x_vec, false);
    _r_vec.assign(x_vec, false);
    _z_vec.assign(x_vec, false);

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
    lin_op.Residual(_r_vec, b_vec, x_vec);
    if (pre_op != nullptr) {
      pre_op->mul(_z_vec, _r_vec);
      _p_vec <<= _z_vec;
      _gamma = dot_product(_r_vec, _z_vec);
    } else {
      _p_vec <<= _r_vec;
      _gamma = dot_product(_r_vec, _r_vec);
    }

    return (pre_op != nullptr) ? norm_2(_r_vec) : std::sqrt(_gamma);
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
    lin_op.mul(_z_vec, _p_vec);
    const real_t alpha = safe_divide(_gamma, dot_product(_p_vec, _z_vec));
    x_vec += alpha * _p_vec;
    _r_vec -= alpha * _z_vec;

    // ----------------------
    // 𝛾̅ ← 𝛾,
    // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
    //   𝒛 ← 𝓟𝒓,
    //   𝛾 ← <𝒓⋅𝒛>,
    // 𝗲𝗹𝘀𝗲:
    //   𝛾 ← <𝒓⋅𝒓>.
    // 𝗲𝗻𝗱 𝗶𝗳
    // ----------------------
    const real_t gamma_bar = _gamma;
    if (pre_op != nullptr) {
      pre_op->mul(_z_vec, _r_vec);
      _gamma = dot_product(_r_vec, _z_vec);
    } else {
      _gamma = dot_product(_r_vec, _r_vec);
    }

    // ----------------------
    // 𝛽 ← 𝛾/𝛾̅,
    // 𝒑 ← (𝓟 ≠ 𝗻𝗼𝗻𝗲 ? 𝒛 : 𝒓) + 𝛽⋅𝒑.
    // ----------------------
    const real_t beta = safe_divide(_gamma, gamma_bar);
    _p_vec <<= (pre_op != nullptr ? _z_vec : _r_vec) + beta * _p_vec;

    return (pre_op != nullptr) ? norm_2(_r_vec) : sqrt(_gamma);
  }

}; // class CgSolver

// -----------------------------------------------------------------------------

} // namespace Storm
