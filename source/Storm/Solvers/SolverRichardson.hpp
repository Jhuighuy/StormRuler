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

#include <Storm/Solvers/Solver.hpp>

#include <utility>

namespace Storm {

// -----------------------------------------------------------------------------

/// @brief The Richardson iteration linear operator equation solver.
///
/// References:
/// @verbatim
/// [1] ???
/// @endverbatim
template<legacy_vector_like Vector>
class RichardsonSolver final : public IterativeSolver<Vector> {
public:

  real_t relaxation_factor = 1.0e-4;

private:

  Vector _r_vec, _z_vec;

  real_t init(const Vector& x_vec, const Vector& b_vec,
              const Operator<Vector>& lin_op,
              const Preconditioner<Vector>* pre_op) override {
    _r_vec.assign(x_vec, false);
    if (pre_op != nullptr) _z_vec.assign(x_vec, false);

    // Initialize:
    // ----------------------
    // 𝒓 ← 𝒃 - 𝓐𝒙,
    // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
    //   𝒛 ← 𝒓,
    //   𝒓 ← 𝓟𝒛.
    // 𝗲𝗻𝗱 𝗶𝗳
    // ----------------------
    lin_op.Residual(_r_vec, b_vec, x_vec);
    if (pre_op != nullptr) {
      std::swap(_z_vec, _r_vec);
      pre_op->mul(_r_vec, _z_vec);
    }

    return norm_2(_r_vec);
  }

  real_t iterate(Vector& x_vec, const Vector& b_vec,
                 const Operator<Vector>& lin_op,
                 const Preconditioner<Vector>* pre_op) override {
    const real_t& omega = relaxation_factor;

    // Update the solution and the residual:
    // ----------------------
    // 𝒙 ← 𝒙 + 𝜔⋅𝒓,
    // 𝒓 ← 𝒃 - 𝓐𝒙,
    // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
    //   𝒛 ← 𝒓,
    //   𝒓 ← 𝓟𝒛.
    // 𝗲𝗻𝗱 𝗶𝗳
    // ----------------------
    x_vec += omega * _r_vec;
    lin_op.Residual(_r_vec, b_vec, x_vec);
    if (pre_op != nullptr) {
      std::swap(_z_vec, _r_vec);
      pre_op->mul(_r_vec, _z_vec);
    }

    return norm_2(_r_vec);
  }

}; // class RichardsonSolver

// -----------------------------------------------------------------------------

} // namespace Storm
