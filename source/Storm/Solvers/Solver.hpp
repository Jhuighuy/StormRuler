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
#include <Storm/Solvers/Preconditioner.hpp>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace Storm {

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract operator equation solver.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<legacy_vector_like InVector, legacy_vector_like OutVector = InVector>
class Solver : public Object {
public:

  /// @brief Solve the operator equation 𝓐(𝒙) = 𝒃.
  ///
  /// @param x_vec Solution vector, 𝒙.
  /// @param b_vec Right-hand-side vector, 𝒃.
  /// @param any_op Equation operator, 𝓐(𝒙).
  ///
  /// @returns Status of operation.
  virtual bool solve(InVector& x_vec, const OutVector& b_vec,
                     const Operator<InVector, OutVector>& any_op) = 0;

}; // class Solver

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract operator equation iterative solver.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<legacy_vector_like InVector, legacy_vector_like OutVector = InVector>
class IterativeSolver : public Solver<InVector, OutVector> {
public:

  size_t iteration{0};
  size_t num_iterations{2000};
  real_t absolute_error{0.0};
  real_t relative_error{0.0};

  real_t absolute_error_tolerance{1.0e-6};
  real_t relative_error_tolerance{1.0e-6};

  PreconditionerSide pre_side{PreconditionerSide::Right};
  std::unique_ptr<Preconditioner<InVector>> pre_op{nullptr};
  std::string name;

protected:

  /// @brief Initialize the iterative solver.
  ///
  /// @param x_vec Initial guess for the solution vector, 𝒙.
  /// @param b_vec Right-hand-side vector, 𝒃.
  /// @param any_op Equation operator, 𝓐(𝒙).
  /// @param pre_op Preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Residual norm of the initial guess, ‖𝒃 - 𝓐(𝒙)‖.
  virtual real_t init(const InVector& x_vec, const OutVector& b_vec,
                      const Operator<InVector, OutVector>& any_op,
                      const Preconditioner<InVector>* pre_op) = 0;

  /// @brief Iterate the solver.
  ///
  /// @param x_vec Solution vector, 𝒙.
  /// @param b_vec Right-hand-side vector, 𝒃.
  /// @param any_op Equation operator, 𝓐(𝒙).
  /// @param pre_op Preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Residual norm, ‖𝒃 - 𝓐(𝒙)‖.
  virtual real_t iterate(InVector& x_vec, const OutVector& b_vec,
                         const Operator<InVector, OutVector>& any_op,
                         const Preconditioner<InVector>* pre_op) = 0;

  /// @brief Finalize the iterations.
  ///
  /// @param x_vec Solution vector, 𝒙.
  /// @param b_vec Right-hand-side vector, 𝒃.
  /// @param any_op Equation operator, 𝓐(𝒙).
  /// @param pre_op Preconditioner operator, 𝓟(𝒙).
  virtual void finalize(InVector& x_vec, const OutVector& b_vec,
                        const Operator<InVector, OutVector>& any_op,
                        const Preconditioner<InVector>* pre_op) {}

public:

  bool solve(InVector& x_vec, const OutVector& b_vec,
             const Operator<InVector, OutVector>& any_op) final {
    // Initialize the solver.
    if (pre_op != nullptr) {
      pre_op->build(x_vec, b_vec, any_op);
    }
    const real_t initial_error = init(x_vec, b_vec, any_op, pre_op.get());
    absolute_error = initial_error;
    if (absolute_error_tolerance > 0.0 &&
        absolute_error < absolute_error_tolerance) {
      finalize(x_vec, b_vec, any_op, pre_op.get());
      return true;
    }

    // Iterate the solver:
    bool converged = false;
    for (iteration = 0; !converged && (iteration < num_iterations);
         ++iteration) {
      absolute_error = iterate(x_vec, b_vec, any_op, pre_op.get());
      relative_error = absolute_error / initial_error;
      converged |= (absolute_error_tolerance > 0.0) &&
                   (absolute_error < absolute_error_tolerance);
      converged |= (relative_error_tolerance > 0.0) &&
                   (relative_error < relative_error_tolerance);
    }

    // Exit the solver.
    finalize(x_vec, b_vec, any_op, pre_op.get());
    STORM_INFO("n_iter: {:>4d}, abs_err: {:>-12e}, rel_err: {:>-12e}", //
               iteration, absolute_error, relative_error);
    return converged;
  }

}; // class IterativeSolver

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract inner-outer iterative solver.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<legacy_vector_like InVector, legacy_vector_like OutVector = InVector>
class InnerOuterIterativeSolver : public IterativeSolver<InVector, OutVector> {
public:

  size_t inner_iteration{0};
  size_t num_inner_iterations{50};

protected:

  /// @brief Initialize the outer iterations.
  ///
  /// This function is used invoked only once,
  ///   in the initialization phase.
  ///
  /// @param x_vec Initial guess for the solution vector, 𝒙.
  /// @param b_vec Right-hand-side vector, 𝒃.
  /// @param any_op Equation operator, 𝓐(𝒙).
  /// @param pre_op Preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Residual norm of the initial guess, ‖𝒃 - 𝓐(𝒙)‖.
  virtual real_t outer_init(const InVector& x_vec, const OutVector& b_vec,
                            const Operator<InVector, OutVector>& any_op,
                            const Preconditioner<InVector>* pre_op) = 0;

  /// @brief Initialize the inner iterations.
  ///
  /// This function is invoked before the each inner iteration loop.
  ///
  /// @param x_vec Solution vector, 𝒙.
  /// @param b_vec Right-hand-side vector, 𝒃.
  /// @param any_op Equation operator, 𝓐(𝒙).
  /// @param pre_op Preconditioner operator, 𝓟(𝒙).
  virtual void inner_init(const InVector& x_vec, const OutVector& b_vec,
                          const Operator<InVector, OutVector>& any_op,
                          const Preconditioner<InVector>* pre_op) {}

  /// @brief Perform the inner iteration.
  ///
  /// @param x_vec Solution vector, 𝒙.
  /// @param b_vec Right-hand-side vector, 𝒃.
  /// @param any_op Equation operator, 𝓐(𝒙).
  /// @param pre_op Preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Residual norm, ‖𝒃 - 𝓐(𝒙)‖.
  virtual real_t inner_iterate(InVector& x_vec, const OutVector& b_vec,
                               const Operator<InVector, OutVector>& any_op,
                               const Preconditioner<InVector>* pre_op) = 0;

  /// @brief Finalize the inner iterations.
  ///
  /// This function is called in order to finalize
  ///   the inner iterations or if some stopping criterion is met.
  ///
  /// @param x_vec Solution vector, 𝒙.
  /// @param b_vec Right-hand-side vector, 𝒃.
  /// @param any_op Equation operator, 𝓐(𝒙).
  /// @param pre_op Preconditioner operator, 𝓟(𝒙).
  virtual void inner_finalize(InVector& x_vec, const OutVector& b_vec,
                              const Operator<InVector, OutVector>& any_op,
                              const Preconditioner<InVector>* pre_op) {}

  /// @brief Finalize the outer iterations.
  ///
  /// This function is used invoked only once,
  ///   when some stopping criterion is met.
  ///
  /// @param x_vec Solution vector, 𝒙.
  /// @param b_vec Right-hand-side vector, 𝒃.
  /// @param any_op Equation operator, 𝓐(𝒙).
  /// @param pre_op Preconditioner operator, 𝓟(𝒙).
  virtual void outer_finalize(InVector& x_vec, const OutVector& b_vec,
                              const Operator<InVector, OutVector>& any_op,
                              const Preconditioner<InVector>* pre_op) {}

private:

  real_t init(const InVector& x_vec, const OutVector& b_vec,
              const Operator<InVector, OutVector>& any_op,
              const Preconditioner<InVector>* pre_op) final {
    return outer_init(x_vec, b_vec, any_op, pre_op);
  }

  real_t iterate(InVector& x_vec, const OutVector& b_vec,
                 const Operator<InVector, OutVector>& any_op,
                 const Preconditioner<InVector>* pre_op) final {
    inner_iteration = this->iteration % num_inner_iterations;
    if (inner_iteration == 0) {
      inner_init(x_vec, b_vec, any_op, pre_op);
    }
    const real_t residual_norm{inner_iterate(x_vec, b_vec, any_op, pre_op)};
    if (inner_iteration == num_inner_iterations - 1) {
      inner_finalize(x_vec, b_vec, any_op, pre_op);
    }
    return residual_norm;
  }

  void finalize(InVector& x_vec, const OutVector& b_vec,
                const Operator<InVector, OutVector>& any_op,
                const Preconditioner<InVector>* pre_op) final {
    if (inner_iteration != num_inner_iterations - 1) {
      inner_finalize(x_vec, b_vec, any_op, pre_op);
    }
    outer_finalize(x_vec, b_vec, any_op, pre_op);
  }

}; // class InnerOuterIterativeSolver

template<template<class> class Solver, class Vector>
bool solve(Vector& x_vec, const Vector& b_vec, const Operator<Vector>& any_op) {
  Solver<Vector> solver{};
  return solver.solve(x_vec, b_vec, any_op);
}

/// ----------------------------------------------------------------- ///
/// @brief Solve an operator equation 𝓐(𝒙) = 𝒃,
///   when 𝓐(𝒙) is a non-uniform operator (𝓐(𝟢) ≠ 𝟢).
/// ----------------------------------------------------------------- ///
template<legacy_vector_like Vector>
bool solve_non_uniform(Solver<Vector>& solver, Vector& x_vec,
                       const Vector& b_vec, const Operator<Vector>& any_op) {
  Vector z_vec, f_vec;

  z_vec.assign(x_vec, false);
  f_vec.assign(b_vec, false);

  // Solve an equation with the "uniformed" operator:
  // 𝓐(𝒙) - 𝓐(𝟢) = 𝒃 - 𝓐(𝟢).
  fill_with(f_vec, 0.0);
  any_op.mul(z_vec, f_vec);
  f_vec <<= b_vec - z_vec;

  const auto uni_op =
      make_operator<Vector>([&](Vector& y_vec, const Vector& x_vec) {
        any_op.mul(y_vec, x_vec);
        y_vec -= z_vec;
      });

  return solver.solve(x_vec, f_vec, *uni_op);
}

} // namespace Storm
