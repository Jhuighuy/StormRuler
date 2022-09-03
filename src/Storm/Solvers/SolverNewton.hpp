/// Copyright (C) 2022 Oleg Butakov
///
/// Permission is hereby granted, free of charge, to any person obtaining a copy
/// of this software and associated documentation files (the "Software"), to
/// deal in the Software without restriction, including without limitation the
/// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
/// sell copies of the Software, and to permit persons to whom the Software is
/// furnished to do so, subject to the following conditions:
///
/// The above copyright notice and this permission notice shall be included in
/// all copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
/// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
/// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
/// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
/// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
/// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
/// IN THE SOFTWARE.

#pragma once

#include <Storm/Base.hpp>

#include <Storm/Utils/Math.hpp>

#include <Storm/Blass/Vector.hpp>

#include <Storm/Solvers/Solver.hpp>

#include <limits>

namespace Storm {

/// @brief The Newton method nonlinear operator equation solver.
///
/// The classical Newton iterations are based on the linearization of 𝓐(𝒙)
/// near 𝒙:
///
/// 𝓐(𝒙̂) ≈ 𝓐(𝒙) + [∂𝓐(𝒙)/∂𝒙](𝒙̂ - 𝒙) = 𝒃,
///
/// or, alternatively:
///
/// [∂𝓐(𝒙)/∂𝒙]𝒕 = 𝒓, 𝒕 = 𝒙̂ - 𝒙, 𝒓 = 𝒃 - 𝓐(𝒙)
///
/// where 𝒙 and 𝒙̂ are the current and updated solution vectors. Therefore, a
/// linear equation has to be solved on each iteration, linear operator
/// 𝓙(𝒙) ≈ ∂𝓐(𝒙)/∂𝒙 for computing Jacobian-vector products is required.
///
/// References:
/// @verbatim
/// [1] ???
/// @endverbatim
template<VectorLike Vector>
class NewtonSolver : public IterativeSolver<Vector> {
private:

  real_t init(const Vector& x_vec, const Vector& b_vec,
              const Operator<Vector>& any_op,
              const Preconditioner<Vector>* pre_op) override final;

  real_t iterate(Vector& x_vec, const Vector& b_vec,
                 const Operator<Vector>& any_op,
                 const Preconditioner<Vector>* pre_op) override final;

}; // class NewtonSolver

template<VectorLike Vector>
real_t NewtonSolver<Vector>::init(const Vector& x_vec, const Vector& b_vec,
                                  const Operator<Vector>& any_op,
                                  const Preconditioner<Vector>* pre_op) {
  STORM_TERMINATE_("Newton solver is not implemented yet!");
} // NewtonSolver::init

template<VectorLike Vector>
real_t NewtonSolver<Vector>::iterate(Vector& x_vec, const Vector& b_vec,
                                     const Operator<Vector>& any_op,
                                     const Preconditioner<Vector>* pre_op) {
  STORM_TERMINATE_("Newton solver is not implemented yet!");
} // NewtonSolver::iterate

/// @brief The first-order JFNK (Jacobian free-Newton-Krylov) nonlinear operator
/// equation solver.
///
/// For the @c Newton iterations, computing of the Jacobian-vector products
/// 𝒛 = 𝓙(𝒙)𝒚, where 𝓙(𝒙) ≈ ∂𝓐(𝒙)/∂𝒙 is required. Consider the expansion:
///
/// 𝓐(𝒙 + 𝛿⋅𝒚) = 𝓐(𝒙) + 𝛿⋅[∂𝓐(𝒙)/∂𝒙]𝒚 + 𝓞(𝛿²),
///
/// where 𝛿 is some small number. Therefore,
///
/// 𝓙(𝒙)𝒚 = [𝓐(𝒙 + 𝛿⋅𝒚) - 𝓐(𝒙)]/𝛿 = [∂𝓐(𝒙)/∂𝒙]𝒚 + 𝓞(𝛿).
///
/// Expression above may be used as the formula for computing the (approximate)
/// Jacobian-vector products. Parameter 𝛿 is commonly defined as [1]:
///
/// 𝛿 = 𝜇⋅‖𝒚‖⁺, 𝜇 = (𝜀ₘ)¹ᐟ²⋅(1 + ‖𝒙‖)¹ᐟ²,
///
/// where 𝜀ₘ is the machine roundoff, ‖𝒚‖⁺ is the pseudo-inverse to ‖𝒚‖.
///
/// References:
/// @verbatim
/// [1] Liu, Wei, Lilun Zhang, Ying Zhong, Yongxian Wang,
///     Yonggang Che, Chuanfu Xu and Xinghua Cheng.
///     “CFD High-order Accurate Scheme JFNK Solver.”
///     Procedia Engineering 61 (2013): 9-15.
/// @endverbatim
template<VectorLike Vector>
class JfnkSolver final : public IterativeSolver<Vector> {
private:

  Vector s_vec_, t_vec_, r_vec_, w_vec_;

  real_t init(const Vector& x_vec, const Vector& b_vec,
              const Operator<Vector>& any_op,
              const Preconditioner<Vector>* pre_op) override;

  real_t iterate(Vector& x_vec, const Vector& b_vec,
                 const Operator<Vector>& any_op,
                 const Preconditioner<Vector>* pre_op) override;

}; // class JfnkSolver

template<VectorLike Vector>
real_t JfnkSolver<Vector>::init(const Vector& x_vec, const Vector& b_vec,
                                const Operator<Vector>& any_op,
                                const Preconditioner<Vector>* pre_op) {
  s_vec_.assign(x_vec, false);
  t_vec_.assign(x_vec, false);
  r_vec_.assign(x_vec, false);
  w_vec_.assign(x_vec, false);

  // Initialize:
  // ----------------------
  // 𝒘 ← 𝓐(𝒙),
  // 𝒓 ← 𝒃 - 𝒘.
  // ----------------------
  any_op.mul(w_vec_, x_vec);
  r_vec_ <<= b_vec - w_vec_;

  return norm_2(r_vec_);

} // JfnkSolver::init

template<VectorLike Vector>
real_t JfnkSolver<Vector>::iterate(Vector& x_vec, const Vector& b_vec,
                                   const Operator<Vector>& any_op,
                                   const Preconditioner<Vector>* pre_op) {
  // Solve the Jacobian equation:
  // ----------------------
  // 𝜇 ← (𝜀ₘ)¹ᐟ²⋅(1 + ‖𝒙‖)]¹ᐟ²,
  // 𝒕 ← 𝒓,
  // 𝒕 ← 𝓙(𝒙)⁻¹𝒓.
  // ----------------------
  static const real_t sqrt_of_epsilon =
      std::sqrt(std::numeric_limits<real_t>::epsilon());
  const real_t mu = sqrt_of_epsilon * math::sqrt(1.0 + norm_2(x_vec));
  t_vec_ <<= r_vec_;
  {
    auto solver = std::make_unique<BiCgStabSolver<Vector>>();
    solver->absolute_error_tolerance = 1.0e-8;
    solver->relative_error_tolerance = 1.0e-8;
    auto op = make_operator<Vector>([&](Vector& z_vec, const Vector& y_vec) {
      // Compute the Jacobian-vector product:
      // ----------------------
      // 𝛿 ← 𝜇⋅‖𝒚‖⁺,
      // 𝒔 ← 𝒙 + 𝛿⋅𝒚,
      // 𝒛 ← 𝓐(𝒔),
      // 𝒛 ← 𝛿⁺⋅𝒛 - 𝛿⁺⋅𝒘.
      // ----------------------
      const real_t delta = math::safe_divide(mu, norm_2(y_vec));
      s_vec_ <<= x_vec + delta * y_vec;
      any_op.mul(z_vec, s_vec_);
      const real_t delta_inverse = math::safe_divide(1.0, delta);
      z_vec <<= delta_inverse * (z_vec - w_vec_);
    });
    solver->solve(t_vec_, r_vec_, *op);
  }

  // Update the solution and the residual:
  // ----------------------
  // 𝒙 ← 𝒙 + 𝒕,
  // 𝒘 ← 𝓐(𝒙),
  // 𝒓 ← 𝒃 - 𝒘.
  // ----------------------
  x_vec += t_vec_;
  any_op.mul(w_vec_, x_vec);
  r_vec_ <<= b_vec - w_vec_;

  return norm_2(r_vec_);

} // JfnkSolver::iterate

} // namespace Storm
