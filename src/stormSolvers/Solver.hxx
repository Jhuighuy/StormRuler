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

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include <stormBase.hxx>
#include <stormSolvers/Operator.hxx>
#include <stormSolvers/Preconditioner.hxx>
#include <stormUtils/Object.hxx>

namespace Storm {

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract operator equation solver.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<VectorLike InVector, VectorLike OutVector = InVector>
class Solver : public Object {
public:

  /// @brief Solve the operator equation 𝓐(𝒙) = 𝒃.
  ///
  /// @param x_vec Solution vector, 𝒙.
  /// @param b_vec Right-hand-side vector, 𝒃.
  /// @param any_op Equation operator, 𝓐(𝒙).
  ///
  /// @returns Status of operation.
  virtual bool Solve(InVector& x_vec, OutVector const& b_vec,
                     Operator<InVector, OutVector> const& any_op) = 0;

}; // class Solver

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract operator equation iterative solver.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<VectorLike InVector, VectorLike OutVector = InVector>
class IterativeSolver : public Solver<InVector, OutVector> {
public:

  size_t Iteration{0};
  size_t NumIterations{2000};
  real_t AbsoluteError{0.0};
  real_t RelativeError{0.0};

  real_t AbsoluteTolerance{1.0e-6};
  real_t RelativeTolerance{1.0e-6};
  bool VerifySolution{false};

  PreconditionerSide PreSide{PreconditionerSide::Right};
  std::unique_ptr<Preconditioner<InVector>> PreOp{nullptr};

protected:

  /// @brief Initialize the iterative solver.
  ///
  /// @param x_vec Initial guess for the solution vector, 𝒙.
  /// @param b_vec Right-hand-side vector, 𝒃.
  /// @param any_op Equation operator, 𝓐(𝒙).
  /// @param pre_op Preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Residual norm of the initial guess, ‖𝒃 - 𝓐(𝒙)‖.
  virtual real_t Init(InVector const& x_vec, OutVector const& b_vec,
                      Operator<InVector, OutVector> const& any_op,
                      Preconditioner<InVector> const* pre_op) = 0;

  /// @brief Iterate the solver.
  ///
  /// @param x_vec Solution vector, 𝒙.
  /// @param b_vec Right-hand-side vector, 𝒃.
  /// @param any_op Equation operator, 𝓐(𝒙).
  /// @param pre_op Preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Residual norm, ‖𝒃 - 𝓐(𝒙)‖.
  virtual real_t Iterate(InVector& x_vec, OutVector const& b_vec,
                         Operator<InVector, OutVector> const& any_op,
                         Preconditioner<InVector> const* pre_op) = 0;

  /// @brief Finalize the iterations.
  ///
  /// @param x_vec Solution vector, 𝒙.
  /// @param b_vec Right-hand-side vector, 𝒃.
  /// @param any_op Equation operator, 𝓐(𝒙).
  /// @param pre_op Preconditioner operator, 𝓟(𝒙).
  virtual void Finalize(InVector& x_vec, OutVector const& b_vec,
                        Operator<InVector, OutVector> const& any_op,
                        Preconditioner<InVector> const* pre_op) {}

public:

  bool Solve(InVector& x_vec, OutVector const& b_vec,
             Operator<InVector, OutVector> const& any_op) override final;

}; // class IterativeSolver

template<VectorLike InVector, VectorLike OutVector>
bool IterativeSolver<InVector, OutVector>::Solve(
    InVector& x_vec, OutVector const& b_vec,
    Operator<InVector, OutVector> const& any_op) {
  // Initialize the solver.
  if (PreOp != nullptr) { PreOp->Build(x_vec, b_vec, any_op); }
  real_t const initialError{
      (AbsoluteError = Init(x_vec, b_vec, any_op, PreOp.get()))};
  if (AbsoluteTolerance > 0.0 && AbsoluteError < AbsoluteTolerance) {
    Finalize(x_vec, b_vec, any_op, PreOp.get());
    return true;
  }

  // Iterate the solver:
  bool converged = false;
  for (Iteration = 0; !converged && (Iteration < NumIterations); ++Iteration) {
    AbsoluteError = Iterate(x_vec, b_vec, any_op, PreOp.get());
    RelativeError = AbsoluteError / initialError;

    converged |=
        (AbsoluteTolerance > 0.0) && (AbsoluteError < AbsoluteTolerance);
    converged |=
        (RelativeTolerance > 0.0) && (RelativeError < RelativeTolerance);
  }

  // Exit the solver.
  Finalize(x_vec, b_vec, any_op, PreOp.get());
  return converged;

} // IterativeSolver::Solve

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract inner-outer iterative solver.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<VectorLike InVector, VectorLike OutVector = InVector>
class InnerOuterIterativeSolver : public IterativeSolver<InVector, OutVector> {
public:

  size_t InnerIteration{0};
  size_t NumInnerIterations{50};

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
  virtual real_t OuterInit(InVector const& x_vec, OutVector const& b_vec,
                           Operator<InVector, OutVector> const& any_op,
                           Preconditioner<InVector> const* pre_op) = 0;

  /// @brief Initialize the inner iterations.
  ///
  /// This function is invoked before the each inner iteration loop.
  ///
  /// @param x_vec Solution vector, 𝒙.
  /// @param b_vec Right-hand-side vector, 𝒃.
  /// @param any_op Equation operator, 𝓐(𝒙).
  /// @param pre_op Preconditioner operator, 𝓟(𝒙).
  virtual void InnerInit(InVector const& x_vec, OutVector const& b_vec,
                         Operator<InVector, OutVector> const& any_op,
                         Preconditioner<InVector> const* pre_op) {}

  /// @brief Perform the inner iteration.
  ///
  /// @param x_vec Solution vector, 𝒙.
  /// @param b_vec Right-hand-side vector, 𝒃.
  /// @param any_op Equation operator, 𝓐(𝒙).
  /// @param pre_op Preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Residual norm, ‖𝒃 - 𝓐(𝒙)‖.
  virtual real_t InnerIterate(InVector& x_vec, OutVector const& b_vec,
                              Operator<InVector, OutVector> const& any_op,
                              Preconditioner<InVector> const* pre_op) = 0;

  /// @brief Finalize the inner iterations.
  ///
  /// This function is called in order to finalize
  ///   the inner iterations or if some stopping criterion is met.
  ///
  /// @param x_vec Solution vector, 𝒙.
  /// @param b_vec Right-hand-side vector, 𝒃.
  /// @param any_op Equation operator, 𝓐(𝒙).
  /// @param pre_op Preconditioner operator, 𝓟(𝒙).
  virtual void InnerFinalize(InVector& x_vec, OutVector const& b_vec,
                             Operator<InVector, OutVector> const& any_op,
                             Preconditioner<InVector> const* pre_op) {}

  /// @brief Finalize the outer iterations.
  ///
  /// This function is used invoked only once,
  ///   when some stopping criterion is met.
  ///
  /// @param x_vec Solution vector, 𝒙.
  /// @param b_vec Right-hand-side vector, 𝒃.
  /// @param any_op Equation operator, 𝓐(𝒙).
  /// @param pre_op Preconditioner operator, 𝓟(𝒙).
  virtual void OuterFinalize(InVector& x_vec, OutVector const& b_vec,
                             Operator<InVector, OutVector> const& any_op,
                             Preconditioner<InVector> const* pre_op) {}

private:

  real_t Init(InVector const& x_vec, OutVector const& b_vec,
              Operator<InVector, OutVector> const& any_op,
              Preconditioner<InVector> const* pre_op) override final {
    return OuterInit(x_vec, b_vec, any_op, pre_op);
  }

  real_t Iterate(InVector& x_vec, OutVector const& b_vec,
                 Operator<InVector, OutVector> const& any_op,
                 Preconditioner<InVector> const* pre_op) override final {
    InnerIteration = this->Iteration % NumInnerIterations;
    if (InnerIteration == 0) { InnerInit(x_vec, b_vec, any_op, pre_op); }
    real_t const residualNorm{InnerIterate(x_vec, b_vec, any_op, pre_op)};
    if (InnerIteration == NumInnerIterations - 1) {
      InnerFinalize(x_vec, b_vec, any_op, pre_op);
    }
    return residualNorm;
  }

  void Finalize(InVector& x_vec, OutVector const& b_vec,
                Operator<InVector, OutVector> const& any_op,
                Preconditioner<InVector> const* pre_op) override final {
    if (InnerIteration != NumInnerIterations - 1) {
      InnerFinalize(x_vec, b_vec, any_op, pre_op);
    }
    OuterFinalize(x_vec, b_vec, any_op, pre_op);
  }

}; // class InnerOuterIterativeSolver

/// ----------------------------------------------------------------- ///
/// @brief Solve an operator equation 𝓐(𝒙) = 𝒃,
///   when 𝓐(𝒙) is a non-uniform operator (𝓐(𝟢) ≠ 𝟢).
/// ----------------------------------------------------------------- ///
template<VectorLike Vector>
bool SolveNonUniform(Solver<Vector>& solver, Vector& x_vec, Vector const& b_vec,
                     Operator<Vector> const& any_op) {
  Vector z_vec, f_vec;

  z_vec.assign(x_vec, false);
  f_vec.assign(b_vec, false);

  // Solve an equation with the "uniformed" operator:
  // 𝓐(𝒙) - 𝓐(𝟢) = 𝒃 - 𝓐(𝟢).
  Blas::Fill(f_vec, 0.0);
  any_op.MatVec(z_vec, f_vec);
  f_vec <<= b_vec - z_vec;

  auto const uniOp =
      MakeOperator<Vector>([&](Vector& y_vec, Vector const& x_vec) {
        any_op.MatVec(y_vec, x_vec);
        y_vec -= z_vec;
      });

  return solver.Solve(x_vec, f_vec, *uniOp);

} // SolveNonUniform

} // namespace Storm
