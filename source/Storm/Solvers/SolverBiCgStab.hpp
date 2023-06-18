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

#include <Storm/Bittern/Matrix.hpp>
#include <Storm/Crow/MathUtils.hpp>

#include <Storm/Solvers/MatrixDense.hpp>
#include <Storm/Solvers/Solver.hpp>

#include <utility>
#include <vector>

namespace Storm {

// -----------------------------------------------------------------------------

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

  real_t _alpha, _rho, _omega;
  Vector _p_vec, _r_vec, _r_tilde_vec, _t_vec, _v_vec, _z_vec;

  real_t init(const Vector& x_vec, const Vector& b_vec,
              const Operator<Vector>& lin_op,
              const Preconditioner<Vector>* pre_op) override {
    const bool left_pre =
        (pre_op != nullptr) && (this->pre_side == PreconditionerSide::Left);

    _p_vec.assign(x_vec, false);
    _r_vec.assign(x_vec, false);
    _r_tilde_vec.assign(x_vec, false);
    _t_vec.assign(x_vec, false);
    _v_vec.assign(x_vec, false);
    if (pre_op != nullptr) _z_vec.assign(x_vec, false);

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
    lin_op.Residual(_r_vec, b_vec, x_vec);
    if (left_pre) {
      std::swap(_z_vec, _r_vec);
      pre_op->mul(_r_vec, _z_vec);
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
      _p_vec <<= _r_vec;
    } else {
      const real_t rho_bar =
          std::exchange(_rho, dot_product(_r_tilde_vec, _r_vec));
      const real_t beta{safe_divide(_alpha * _rho, _omega * rho_bar)};
      _p_vec <<= _r_vec + beta * (_p_vec - _omega * _v_vec);
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
    if (left_pre) pre_op->mul(_v_vec, _z_vec, lin_op, _p_vec);
    else if (right_pre) lin_op.mul(_v_vec, _z_vec, *pre_op, _p_vec);
    else lin_op.mul(_v_vec, _p_vec);

    _alpha = safe_divide(_rho, dot_product(_r_tilde_vec, _v_vec));
    x_vec += _alpha * (right_pre ? _z_vec : _p_vec);
    _r_vec -= _alpha * _v_vec;

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
    if (left_pre) pre_op->mul(_t_vec, _z_vec, lin_op, _r_vec);
    else if (right_pre) lin_op.mul(_t_vec, _z_vec, *pre_op, _r_vec);
    else lin_op.mul(_t_vec, _r_vec);
    _omega =
        safe_divide(dot_product(_t_vec, _r_vec), dot_product(_t_vec, _t_vec));
    x_vec += _omega * (right_pre ? _z_vec : _r_vec);
    _r_vec -= _omega * _t_vec;

    return norm_2(_r_vec);
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

  real_t _alpha, _rho, _omega;
  DenseVector<real_t> _gamma, _gamma_bar, _gamma_bbar, _sigma;
  DenseMatrix<real_t> _tau;
  Vector _r_tilde_vec, _z_vec;
  std::vector<Vector> _r_vecs, _u_vecs;

  real_t outer_init(const Vector& x_vec, const Vector& b_vec,
                    const Operator<Vector>& lin_op,
                    const Preconditioner<Vector>* pre_op) override {
    const size_t l = this->num_inner_iterations;

    _gamma.assign(l + 1);
    _gamma_bar.assign(l + 1);
    _gamma_bbar.assign(l + 1);
    _sigma.assign(l + 1);
    _tau.assign(l + 1, l + 1);

    _r_tilde_vec.assign(x_vec, false);
    if (pre_op != nullptr) _z_vec.assign(x_vec, false);

    _r_vecs.resize(l + 1);
    _u_vecs.resize(l + 1);
    for (Vector& r_vec : _r_vecs) r_vec.assign(x_vec, false);
    for (Vector& u_vec : _u_vecs) u_vec.assign(x_vec, false);

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
    fill_with(_u_vecs[0], 0.0);
    lin_op.Residual(_r_vecs[0], b_vec, x_vec);
    if (pre_op != nullptr) {
      std::swap(_z_vec, _r_vecs[0]);
      pre_op->mul(_r_vecs[0], _z_vec);
    }
    _r_tilde_vec <<= _r_vecs[0];
    _rho = dot_product(_r_tilde_vec, _r_vecs[0]);

    return sqrt(_rho);
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
      _u_vecs[0] <<= _r_vecs[0];
    } else {
      const real_t rho_bar =
          std::exchange(_rho, dot_product(_r_tilde_vec, _r_vecs[j]));
      const real_t beta = safe_divide(_alpha * _rho, rho_bar);
      for (size_t i = 0; i <= j; ++i) {
        _u_vecs[i] <<= _r_vecs[i] - beta * _u_vecs[i];
      }
    }
    if (pre_op != nullptr) {
      pre_op->mul(_u_vecs[j + 1], _z_vec, lin_op, _u_vecs[j]);
    } else {
      lin_op.mul(_u_vecs[j + 1], _u_vecs[j]);
    }
    _alpha = safe_divide(_rho, dot_product(_r_tilde_vec, _u_vecs[j + 1]));
    for (size_t i = 0; i <= j; ++i) {
      _r_vecs[i] -= _alpha * _u_vecs[i + 1];
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
    x_vec += _alpha * _u_vecs[0];
    if (pre_op != nullptr) {
      pre_op->mul(_r_vecs[j + 1], _z_vec, lin_op, _r_vecs[j]);
    } else {
      lin_op.mul(_r_vecs[j + 1], _r_vecs[j]);
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
          _tau(i, j) =
              safe_divide(dot_product(_r_vecs[i], _r_vecs[j]), _sigma(i));
          _r_vecs[j] -= _tau(i, j) * _r_vecs[i];
        }
        _sigma(j) = dot_product(_r_vecs[j], _r_vecs[j]);
        _gamma_bar(j) =
            safe_divide(dot_product(_r_vecs[0], _r_vecs[j]), _sigma(j));
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
      _omega = _gamma(l) = _gamma_bar(l), _rho *= -_omega;
      for (size_t j = l - 1; j != 0; --j) {
        _gamma(j) = _gamma_bar(j);
        for (size_t i = j + 1; i <= l; ++i) {
          _gamma(j) -= _tau(j, i) * _gamma(i);
        }
      }
      for (size_t j = 1; j < l; ++j) {
        _gamma_bbar(j) = _gamma(j + 1);
        for (size_t i = j + 1; i < l; ++i) {
          _gamma_bbar(j) += _tau(j, i) * _gamma(i + 1);
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
      x_vec += _gamma(1) * _r_vecs[0];
      _r_vecs[0] -= _gamma_bar(l) * _r_vecs[l];
      _u_vecs[0] -= _gamma(l) * _u_vecs[l];
      for (size_t j = 1; j < l; ++j) {
        x_vec += _gamma_bbar(j) * _r_vecs[j];
        _r_vecs[0] -= _gamma_bar(j) * _r_vecs[j];
        _u_vecs[0] -= _gamma(j) * _u_vecs[j];
      }
    }

    return norm_2(_r_vecs[0]);
  }

public:

  BiCgStabLSolver() {
    this->num_inner_iterations = 2;
  }

}; // class BiCgStabLSolver

// -----------------------------------------------------------------------------

} // namespace Storm
