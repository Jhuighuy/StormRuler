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

#include <Storm/Solvers/Preconditioner.hpp>

#include <cmath>

namespace Storm {

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief @c Chebyshev polynomial preconditioner.
///
/// @c Chebyshev preconditioner can operate in the matrix-free mode.
///
/// @verbatim
/// [1] Saad, Yousef.
///     “Iterative methods for sparse linear systems.” (2003).
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<legacy_vector_like Vector>
class ChebyshevPreconditioner final : public Preconditioner<Vector> {
public:

  size_t degree{5};

private:

  real_t theta_, delta_;
  mutable Vector r_vec_, p_vec_;
  const Operator<Vector>* LinOp_;

  void build(const Vector& x_vec, const Vector& b_vec,
             const Operator<Vector>& lin_op) override;

  void mul(Vector& y_vec, const Vector& x_vec) const override;

  void conj_mul(Vector& x_vec, const Vector& y_vec) const override {
    mul(x_vec, y_vec);
  }

}; // class ChebyshevPreconditioner

template<legacy_vector_like Vector>
void ChebyshevPreconditioner<Vector>::build(const Vector& x_vec,
                                            const Vector& b_vec,
                                            const Operator<Vector>& lin_op) {
  r_vec_.assign(x_vec, false);
  p_vec_.assign(x_vec, false);
  this->LinOp_ = &lin_op;

  // PowerIterations<Vector> powerIterations;
  // const real_t beta =
  //   powerIterations.EstimateLargestEigenvalue(p_vec_, lin_op);
  // const real_t alpha = 0.01*beta;

  /// @todo: Estimate the true eigenvalue bounds!
  const real_t alpha = 0.95 * 1.046599390654509e+00;
  const real_t beta = 1.05 * 8.003575342439456e+02;

  // Initialize the center and the semi-major
  // axis of ellipse containing the eigenvalues:
  // ----------------------
  // 𝜃 ← ½⋅(𝛽 + 𝛼),
  // 𝛿 ← ½⋅(𝛽 - 𝛼).
  // ----------------------
  theta_ = 0.5 * (beta + alpha);
  delta_ = 0.5 * (beta - alpha);

} // ChebyshevPreconditioner::Build

template<legacy_vector_like Vector>
void ChebyshevPreconditioner<Vector>::mul(Vector& y_vec,
                                          const Vector& x_vec) const {
  // Initialize the solution:
  // ----------------------
  // 𝛼 ← 𝟤/𝜃,
  // 𝒑 ← 𝒙/𝜃,
  // 𝒚 ← 𝒑.
  // ----------------------
  real_t alpha{2.0 / theta_};
  p_vec_ <<= x_vec / theta_;
  y_vec <<= p_vec_;

  for (size_t k = 0; k < degree; ++k) {
    // Compute the residual and update the solution:
    // ----------------------
    // 𝛼 ← 𝟣/(𝜃 - ¼⋅𝛼⋅𝛿²),
    // 𝛽 ← 𝛼⋅𝜃 - 𝟣,
    // 𝒓 ← 𝒙 - 𝓐𝒚,
    // 𝒑 ← 𝛼⋅𝒓 + 𝛽⋅𝒑,
    // 𝒚 ← 𝒚 + 𝒑.
    // ----------------------
    alpha = 1.0 / (theta_ - 0.25 * alpha * std::pow(delta_, 2));
    const real_t beta = alpha * theta_ - 1.0;
    LinOp_->Residual(r_vec_, x_vec, y_vec);
    p_vec_ <<= alpha * r_vec_ + beta * p_vec_;
    y_vec += p_vec_;
  }

} // ChebyshevPreconditioner::mul

} // namespace Storm
