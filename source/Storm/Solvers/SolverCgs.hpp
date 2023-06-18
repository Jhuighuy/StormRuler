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

#include <utility>

namespace Storm {

// -----------------------------------------------------------------------------

/// @brief The CGS (Conjugate Gradients Squared) linear operator equation
/// solver.
///
/// CGS, like the other BiCG type solvers, requires two operator multiplications
/// per iteration.
///
/// @warning CGS convergence behavior may be very erratic.
///
/// References:
/// @verbatim
/// [1] Sonneveld, Peter.
///     “CGS, A Fast Lanczos-Type Solver for Nonsymmetric Linear systems.”
///     SIAM J. Sci. Stat. Comput., 10:36-52, 1989.
/// @endverbatim
template<legacy_vector_like Vector>
class CgsSolver final : public IterativeSolver<Vector> {
private:

  real_t _rho;
  Vector _p_vec, _q_vec, _r_vec, _r_tilde_vec, _u_vec, _v_vec;

  real_t init(const Vector& x_vec, const Vector& b_vec,
              const Operator<Vector>& lin_op,
              const Preconditioner<Vector>* pre_op) override {
    const bool left_pre =
        (pre_op != nullptr) && (this->pre_side == PreconditionerSide::Left);

    _p_vec.assign(x_vec, false);
    _q_vec.assign(x_vec, false);
    _r_vec.assign(x_vec, false);
    _r_tilde_vec.assign(x_vec, false);
    _u_vec.assign(x_vec, false);
    _v_vec.assign(x_vec, false);

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
    lin_op.Residual(_r_vec, b_vec, x_vec);
    if (left_pre) {
      std::swap(_u_vec, _r_vec);
      pre_op->mul(_r_vec, _u_vec);
    }
    _r_tilde_vec <<= _r_vec;
    _rho = dot_product(_r_tilde_vec, _r_vec);

    return sqrt(_rho);
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
    const bool first_iteration = this->iteration == 0;
    if (first_iteration) {
      _u_vec <<= _r_vec;
      _p_vec <<= _u_vec;
    } else {
      const real_t rho_bar =
          std::exchange(_rho, dot_product(_r_tilde_vec, _r_vec));
      const real_t beta = safe_divide(_rho, rho_bar);
      _u_vec <<= _r_vec + beta * _q_vec;
      _p_vec <<= _u_vec + beta * (_q_vec + beta * _p_vec);
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
    if (left_pre) pre_op->mul(_v_vec, _q_vec, lin_op, _p_vec);
    else if (right_pre) lin_op.mul(_v_vec, _q_vec, *pre_op, _p_vec);
    else lin_op.mul(_v_vec, _p_vec);
    const real_t alpha = safe_divide(_rho, dot_product(_r_tilde_vec, _v_vec));
    _q_vec <<= _u_vec - alpha * _v_vec;
    _v_vec <<= _u_vec + _q_vec;

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
      x_vec += alpha * _v_vec;
      pre_op->mul(_v_vec, _u_vec, lin_op, _v_vec);
      _r_vec -= alpha * _v_vec;
    } else if (right_pre) {
      lin_op.mul(_v_vec, _u_vec, *pre_op, _v_vec);
      x_vec += alpha * _u_vec;
      _r_vec -= alpha * _v_vec;
    } else {
      lin_op.mul(_u_vec, _v_vec);
      x_vec += alpha * _v_vec;
      _r_vec -= alpha * _u_vec;
    }

    return norm_2(_r_vec);
  }

}; // class CgsSolver

// -----------------------------------------------------------------------------

} // namespace Storm
