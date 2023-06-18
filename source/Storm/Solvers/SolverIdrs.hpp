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

/// @brief The IDR(s) (Induced Dimension Reduction) linear operator equation
/// solver.
///
/// References:
/// @verbatim
/// [1] Peter Sonneveld, Martin B. van Gijzen.
///     “IDR(s): A Family of Simple and Fast Algorithms for Solving Large
///      Nonsymmetric Systems of Linear Equations.”
///     SIAM J. Sci. Comput. 31 (2008): 1035-1062.
/// [2] Martin B. van Gijzen, Peter Sonneveld.
///     “Algorithm 913: An Elegant IDR(s) Variant that Efficiently Exploits
///      Biorthogonality Properties.”
///      ACM Trans. Math. Softw. 38 (2011): 5:1-5:19.
/// @endverbatim
template<legacy_vector_like Vector>
class IdrsSolver final : public InnerOuterIterativeSolver<Vector> {
private:

  real_t _omega;
  DenseVector<real_t> _phi, _gamma;
  DenseMatrix<real_t> _mu;
  Vector _r_vec, _v_vec, _z_vec;
  std::vector<Vector> _p_vecs, _u_vecs, _g_vecs;

  real_t outer_init(const Vector& x_vec, const Vector& b_vec,
                    const Operator<Vector>& lin_op,
                    const Preconditioner<Vector>* pre_op) override {
    const size_t s = this->num_inner_iterations;

    const bool left_pre =
        (pre_op != nullptr) && (this->pre_side == PreconditionerSide::Left);

    _phi.assign(s);
    _gamma.assign(s);
    _mu.assign(s, s);

    _r_vec.assign(x_vec, false);
    _v_vec.assign(x_vec, false);
    if (pre_op != nullptr) _z_vec.assign(x_vec, false);

    _p_vecs.resize(s);
    _u_vecs.resize(s);
    _g_vecs.resize(s);
    for (Vector& p_vec : _p_vecs) p_vec.assign(x_vec, false);
    for (Vector& u_vec : _u_vecs) u_vec.assign(x_vec, false);
    for (Vector& g_vec : _g_vecs) g_vec.assign(x_vec, false);

    // Initialize:
    // ----------------------
    // 𝒓 ← 𝒃 - 𝓐𝒙,
    // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
    //   𝒛 ← 𝒓,
    //   𝒓 ← 𝓟𝒛.
    // 𝗲𝗻𝗱 𝗶𝗳
    // 𝜑₀ ← ‖𝒓‖.
    // ----------------------
    lin_op.Residual(_r_vec, b_vec, x_vec);
    if (left_pre) {
      std::swap(_z_vec, _r_vec);
      pre_op->mul(_r_vec, _z_vec);
    }
    _phi(0) = norm_2(_r_vec);

    return _phi(0);
  }

  void inner_init(const Vector& x_vec, const Vector& b_vec,
                  const Operator<Vector>& lin_op,
                  const Preconditioner<Vector>* pre_op) override {
    const size_t s = this->num_inner_iterations;

    // Build shadow space and initialize 𝜑:
    // ----------------------
    // 𝗶𝗳 𝘍𝘪𝘳𝘴𝘵𝘐𝘵𝘦𝘳𝘢𝘵𝘪𝘰𝘯:
    //   𝜔 ← 𝜇₀₀ ← 𝟣,
    //   𝒑₀ ← 𝒓/𝜑₀,
    //   𝗳𝗼𝗿 𝑖 = 𝟣, 𝑠 - 𝟣 𝗱𝗼:
    //     𝜇ᵢᵢ ← 𝟣, 𝜑ᵢ ← 𝟢,
    //     𝒑ᵢ ← 𝘙𝘢𝘯𝘥𝘰𝘮,
    //     𝗳𝗼𝗿 𝑗 = 𝟢, 𝑖 - 𝟣 𝗱𝗼:
    //       𝜇ᵢⱼ ← 𝟢,
    //       𝒑ᵢ ← 𝒑ᵢ - <𝒑ᵢ⋅𝒑ⱼ>⋅𝒑ⱼ,
    //     𝗲𝗻𝗱 𝗳𝗼𝗿
    //     𝒑ᵢ ← 𝒑ᵢ/‖𝒑ᵢ‖.
    //   𝗲𝗻𝗱 𝗳𝗼𝗿
    // 𝗲𝗹𝘀𝗲:
    //   𝗳𝗼𝗿 𝑖 = 𝟢, 𝑠 - 𝟣 𝗱𝗼:
    //     𝜑ᵢ ← <𝒑ᵢ⋅𝒓>.
    //   𝗲𝗻𝗱 𝗳𝗼𝗿
    // 𝗲𝗻𝗱 𝗶𝗳
    // ----------------------
    const bool first_iteration = this->iteration == 0;
    if (first_iteration) {
      _omega = _mu(0, 0) = 1.0;
      _p_vecs[0] <<= _r_vec / _phi(0);
      for (size_t i = 1; i < s; ++i) {
        _mu(i, i) = 1.0, _phi(i) = 0.0;
        fill_randomly(_p_vecs[i]);
        for (size_t j = 0; j < i; ++j) {
          _mu(i, j) = 0.0;
          _p_vecs[i] -= dot_product(_p_vecs[i], _p_vecs[j]) * _p_vecs[j];
        }
        _p_vecs[i] /= norm_2(_p_vecs[i]);
      }
    } else {
      for (size_t i = 0; i < s; ++i) {
        _phi(i) = dot_product(_p_vecs[i], _r_vec);
      }
    }
  }

  real_t inner_iterate(Vector& x_vec, const Vector& b_vec,
                       const Operator<Vector>& lin_op,
                       const Preconditioner<Vector>* pre_op) override {
    const size_t s = this->num_inner_iterations;
    const size_t k = this->inner_iteration;

    const bool left_pre =
        (pre_op != nullptr) && (this->pre_side == PreconditionerSide::Left);
    const bool right_pre =
        (pre_op != nullptr) && (this->pre_side == PreconditionerSide::Right);

    // Compute 𝛾:
    // ----------------------
    // 𝛾ₖ:ₛ₋₁ ← (𝜇ₖ:ₛ₋₁,ₖ:ₛ₋₁)⁻¹⋅𝜑ₖ:ₛ₋₁.
    // ----------------------
    /// @todo:
    /// slice(_gamma, {k, s}) =
    ///    solve(slice(_mu, {k, s}, {k, s}), slice(_phi, {k, s}));
    for (size_t i = k; i < s; ++i) {
      _gamma(i) = _phi(i);
      for (size_t j = k; j < i; ++j) {
        _gamma(i) -= _mu(i, j) * _gamma(j);
      }
      _gamma(i) /= _mu(i, i);
    }

    // Compute the new 𝒈ₖ and 𝒖ₖ vectors:
    // ----------------------
    // 𝒗 ← 𝒓 - 𝛾ₖ⋅𝒈ₖ,
    // 𝗳𝗼𝗿 𝑖 = 𝑘 + 𝟣, 𝑠 - 𝟣 𝗱𝗼:
    //   𝒗 ← 𝒗 - 𝛾ᵢ⋅𝒈ᵢ,
    // 𝗲𝗻𝗱 𝗳𝗼𝗿
    // 𝗶𝗳 𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦:
    //   𝒛 ← 𝒗,
    //   𝒗 ← 𝓟𝒛,
    // 𝗲𝗻𝗱 𝗶𝗳
    // 𝒖ₖ ← 𝜔⋅𝒗 + 𝛾ₖ⋅𝒖ₖ,
    // 𝗳𝗼𝗿 𝑖 = 𝑘 + 𝟣, 𝑠 - 𝟣 𝗱𝗼:
    //   𝒖ₖ ← 𝒖ₖ + 𝛾ᵢ⋅𝒖ᵢ,
    // 𝗲𝗻𝗱 𝗳𝗼𝗿
    // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
    //   𝒈ₖ ← 𝓟(𝒛 ← 𝓐𝒖ₖ).
    // 𝗲𝗹𝘀𝗲:
    //   𝒈ₖ ← 𝓐𝒖ₖ.
    // 𝗲𝗻𝗱 𝗶𝗳
    // ----------------------
    _v_vec <<= _r_vec - _gamma(k) * _g_vecs[k];
    for (size_t i = k + 1; i < s; ++i) {
      _v_vec -= _gamma(i) * _g_vecs[i];
    }
    if (right_pre) {
      std::swap(_z_vec, _v_vec);
      pre_op->mul(_v_vec, _z_vec);
    }
    _u_vecs[k] <<= _omega * _v_vec + _gamma(k) * _u_vecs[k];
    for (size_t i = k + 1; i < s; ++i) {
      _u_vecs[k] += _gamma(i) * _u_vecs[i];
    }
    if (left_pre) {
      pre_op->mul(_g_vecs[k], _z_vec, lin_op, _u_vecs[k]);
    } else {
      lin_op.mul(_g_vecs[k], _u_vecs[k]);
    }

    // Biorthogonalize the new vectors 𝒈ₖ and 𝒖ₖ:
    // ----------------------
    // 𝗳𝗼𝗿 𝑖 = 𝟢, 𝑘 - 𝟣 𝗱𝗼:
    //   𝛼 ← <𝒑ᵢ⋅𝒈ₖ>/𝜇ᵢᵢ,
    //   𝒖ₖ ← 𝒖ₖ - 𝛼⋅𝒖ᵢ,
    //   𝒈ₖ ← 𝒈ₖ - 𝛼⋅𝒈ᵢ.
    // 𝗲𝗻𝗱 𝗳𝗼𝗿
    // ----------------------
    for (size_t i = 0; i < k; ++i) {
      const real_t alpha =
          safe_divide(dot_product(_p_vecs[i], _g_vecs[k]), _mu(i, i));
      _u_vecs[k] -= alpha * _u_vecs[i];
      _g_vecs[k] -= alpha * _g_vecs[i];
    }

    // Compute the new column of 𝜇:
    // ----------------------
    // 𝗳𝗼𝗿 𝑖 = 𝑘, 𝑠 - 𝟣 𝗱𝗼:
    //   𝜇ᵢₖ ← <𝒑ᵢ⋅𝒈ₖ>.
    // 𝗲𝗻𝗱 𝗳𝗼𝗿
    // ----------------------
    for (size_t i = k; i < s; ++i) {
      _mu(i, k) = dot_product(_p_vecs[i], _g_vecs[k]);
    }

    // Update the solution and the residual:
    // ----------------------
    // 𝛽 ← 𝜑ₖ/𝜇ₖₖ,
    // 𝒙 ← 𝒙 + 𝛽⋅𝒖ₖ,
    // 𝒓 ← 𝒓 - 𝛽⋅𝒈ₖ.
    // ----------------------
    const real_t beta = safe_divide(_phi(k), _mu(k, k));
    x_vec += beta * _u_vecs[k];
    _r_vec -= beta * _g_vecs[k];

    // Update 𝜑:
    // ----------------------
    // 𝜑ₖ₊₁:ₛ₋₁ ← 𝜑ₖ₊₁:ₛ₋₁ - 𝛽⋅𝜇ₖ₊₁:ₛ₋₁,ₖ.
    // ----------------------
    /// @todo:
    /// slice(_phi, {k + 1, s}) -= beta * slice(_mu, {k + 1, s}, k);
    for (size_t i = k + 1; i < s; ++i) {
      _phi(i) -= beta * _mu(i, k);
    }

    if (k == s - 1) {
      // Enter the next 𝓖 subspace:
      // ----------------------
      // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
      //   𝒗 ← 𝓟(𝒛 ← 𝓐𝒓),
      // 𝗲𝗹𝘀𝗲 𝗶𝗳 𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦:
      //   𝒗 ← 𝓐(𝒛 ← 𝓟𝒓),
      // 𝗲𝗹𝘀𝗲:
      //   𝒗 ← 𝓐𝒓,
      // 𝗲𝗻𝗱 𝗶𝗳
      // 𝜔 ← <𝒗⋅𝒓>/<𝒗⋅𝒗>,
      // 𝒙 ← 𝒙 + 𝜔⋅(𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦 ? 𝒛 : 𝒓),
      // 𝒓 ← 𝒓 - 𝜔⋅𝒗.
      // ----------------------
      if (left_pre) pre_op->mul(_v_vec, _z_vec, lin_op, _r_vec);
      else if (right_pre) lin_op.mul(_v_vec, _z_vec, *pre_op, _r_vec);
      else lin_op.mul(_v_vec, _r_vec);

      _omega =
          safe_divide(dot_product(_v_vec, _r_vec), dot_product(_v_vec, _v_vec));
      x_vec += _omega * (right_pre ? _z_vec : _r_vec);
      _r_vec -= _omega * _v_vec;
    }

    return norm_2(_r_vec);
  }

public:

  IdrsSolver() {
    this->num_inner_iterations = 4;
  }

}; // class IdrsSolver

// -----------------------------------------------------------------------------

} // namespace Storm
