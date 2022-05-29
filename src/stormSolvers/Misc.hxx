/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// Copyright (C) 2022 Oleg Butakov
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

#include <stormBase.hxx>
#include <stormSolvers/Operator.hxx>

namespace Storm {

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Largest eigenvalue estimator based on the Power Iterations.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<VectorLike Vector>
class PowerIterations final : public Object {
public:

  /// @brief Estimate the largest eigenvalue of \
  ///   the linear operator 𝓐 using the Power Iterations method.
  ///
  /// @param xVec On input: a non-zero vector that is used as \
  ///   the initial guess for the Power iterations; on output: \
  ///   estimate of the eigenvector, corresponding to the largest eigenvalue.
  /// @param linOp Linear operator, 𝓐𝒙.
  /// @param maxIterations Maximum number of the iterations.
  /// @param relativeTolerance Relative error tolerance \
  ///   to terminate the iterations before the maximum number is reached.
  ///
  /// @returns Estimate the largest eigenvalue of 𝓐.
  static real_t
    EstimateLargestEigenvalue(Vector& xVec,
                              Operator<Vector> const& linOp,
                              size_t maxIterations = 20,
                              real_t relativeTolerance = 1.0e-8);

}; // class PowerIterations

template<VectorLike Vector>
real_t PowerIterations<Vector>::
    EstimateLargestEigenvalue(Vector& xVec,
                              Operator<Vector> const& linOp,
                              size_t maxIterations,
                              real_t relativeTolerance) {

  Vector yVec;
  yVec.Assign(xVec, false);

  // ----------------------
  // Initialize the Power Iterations:
  // 𝜆 ← 𝟣,
  // 𝒙 ← 𝘙𝘢𝘯𝘥𝘰𝘮(),
  // 𝒙 ← 𝒙/‖𝒙‖.
  // ----------------------
  real_t lambda = 1.0;
  xVec.RandFill();
  xVec.ScaleAssign(1.0/xVec.Norm2());

  for (size_t iteration = 0; iteration < maxIterations; ++iteration) {

    // ----------------------
    // Continue the Power Iterations:
    // 𝒚 ← 𝓐𝒙,
    // 𝜆̅ ← 𝜆, 𝜆 ← <𝒙⋅𝒚>,
    // 𝒙 ← 𝒚/‖𝒚‖.
    // ----------------------
    linOp.MatVec(yVec, xVec);
    //real_t const lambdaBar = lambda;
    lambda = xVec.Dot(yVec);
    xVec.Scale(yVec, 1.0/yVec.Norm2());

    // ----------------------
    // Check for the convergence on 𝜆 and 𝜆̅:
    // ----------------------
    //if (std::abs((lambda - lambdaBar)/lambdaBar) < relativeTolerance) {
    //  break;
    //}
  }

  return lambda;

} // PowerIterations::EstimateLargestEigenvalue

} // namespace Storm
