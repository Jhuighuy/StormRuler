/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// Copyright (C) 2020-2023 Oleg Butakov
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

#include <Storm/Base.hpp>

#include <Storm/Solvers/Operator.hpp>

namespace Storm {

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Largest eigenvalue estimator based on the Power Iterations.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<legacy_vector_like Vector>
class PowerIterations final {
public:

  /// @brief Estimate the largest eigenvalue of
  ///   the linear operator 𝓐 using the Power Iterations method.
  ///
  /// @param x_vec On input: a non-zero vector that is used as
  ///   the initial guess for the Power iterations; on output:
  ///   estimate of the eigenvector, corresponding to the largest eigenvalue.
  /// @param lin_op Linear operator, 𝓐𝒙.
  /// @param max_iters Maximum number of the iterations.
  /// @param relative_tolerance Relative error tolerance
  ///   to terminate the iterations before the maximum number is reached.
  ///
  /// @returns Estimate the largest eigenvalue of 𝓐.
  static real_t estimate_largest_eigenvalue(Vector& x_vec,
                                            const Operator<Vector>& lin_op,
                                            size_t max_iters = 20,
                                            real_t relative_tolerance = 1.0e-8);

}; // class PowerIterations

template<legacy_vector_like Vector>
real_t PowerIterations<Vector>::estimate_largest_eigenvalue(
    Vector& x_vec, const Operator<Vector>& lin_op, size_t max_iters,
    real_t relative_tolerance) {
  Vector y_vec;
  y_vec.assign(x_vec, false);

  // Initialize the power iterations:
  // ----------------------
  // 𝜆 ← 𝟣,
  // 𝒙 ← 𝘙𝘢𝘯𝘥𝘰𝘮(),
  // 𝒙 ← 𝒙/‖𝒙‖.
  // ----------------------
  real_t lambda{1.0};
  x_vec.RandFill();
  x_vec.ScaleAssign(1.0 / x_vec.Norm2());

  for (size_t k{0}; k < max_iters; ++k) {
    // Continue the power iterations:
    // ----------------------
    // 𝒚 ← 𝓐𝒙,
    // 𝜆̅ ← 𝜆, 𝜆 ← <𝒙⋅𝒚>,
    // 𝒙 ← 𝒚/‖𝒚‖.
    // ----------------------
    lin_op.mul(y_vec, x_vec);
    const real_t lambda_bar{lambda};
    lambda = x_vec.Dot(y_vec);
    x_vec.Scale(y_vec, 1.0 / y_vec.Norm2());

    // Check for the convergence on 𝜆 and 𝜆̅:
    if (std::abs((lambda - lambda_bar) / lambda_bar) < relative_tolerance) {
      break;
    }
  }

  return lambda;

} // PowerIterations::estimate_largest_eigenvalue

} // namespace Storm
