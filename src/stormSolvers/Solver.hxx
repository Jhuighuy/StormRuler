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

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdexcept>

#include <stormBase.hxx>
#include <stormUtils/Class.hxx>

#include <stormSolvers/Operator.hxx>
#include <stormSolvers/Preconditioner.hxx>

namespace Storm {

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract operator equation solver.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<VectorLike InVector, VectorLike OutVector = InVector>
class Solver : public Object {
public:
  STORM_CLASS_(Solver, Object)

  /// @brief Solve the operator equation 𝓐(𝒙) = 𝒃.
  ///
  /// @param xVec Solution vector, 𝒙.
  /// @param bVec Right-hand-side vector, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  ///
  /// @returns Status of operation.
  virtual bool Solve(InVector& xVec,
                     OutVector const& bVec,
                     Operator<InVector, OutVector> const& anyOp) = 0;

}; // class Solver<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract operator equation iterative solver.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<VectorLike InVector, VectorLike OutVector = InVector>
class IterativeSolver : public Solver<InVector, OutVector> {
public:
  STORM_CLASS_(IterativeSolver, Solver<InVector, OutVector>)

  STORM_FIELD_(size_t, Iteration, {0})
  STORM_FIELD_(size_t, NumIterations, {2000})
  STORM_FIELD_(real_t, AbsoluteError, {0.0}) 
  STORM_FIELD_(real_t, RelativeError, {0.0}) 

  STORM_FIELD_(real_t, AbsoluteTolerance, {1.0e-6}) 
  STORM_FIELD_(real_t, RelativeTolerance, {1.0e-6}) 
  STORM_FIELD_(bool, VerifySolution, {true})

  STORM_FIELD_(PreconditionerSide, PreSide, {PreconditionerSide::Right})
  STORM_FIELD_(std::unique_ptr<Preconditioner<InVector>>, PreOp, {nullptr})

protected:

  /// @brief Initialize the iterative solver.
  ///
  /// @param xVec Initial guess for the solution vector, 𝒙.
  /// @param bVec Right-hand-side vector, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  /// @param preOp Preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Residual norm of the initial guess, ‖𝒃 - 𝓐(𝒙)‖.
  virtual real_t Init(InVector const& xVec,
                      OutVector const& bVec,
                      Operator<InVector, OutVector> const& anyOp,
                      Preconditioner<InVector> const* preOp) = 0;

  /// @brief Iterate the solver.
  ///
  /// @param xVec Solution vector, 𝒙.
  /// @param bVec Right-hand-side vector, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  /// @param preOp Preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Residual norm, ‖𝒃 - 𝓐(𝒙)‖.
  virtual real_t Iterate(InVector& xVec,
                         OutVector const& bVec,
                         Operator<InVector, OutVector> const& anyOp,
                         Preconditioner<InVector> const* preOp) = 0;

  /// @brief Finalize the iterations.
  ///
  /// @param xVec Solution vector, 𝒙.
  /// @param bVec Right-hand-side vector, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  /// @param preOp Preconditioner operator, 𝓟(𝒙).
  virtual void Finalize(InVector& xVec,
                        OutVector const& bVec,
                        Operator<InVector, OutVector> const& anyOp,
                        Preconditioner<InVector> const* preOp) {}

public:

  bool Solve(InVector& xVec,
             OutVector const& bVec,
             Operator<InVector, OutVector> const& anyOp) override final;

}; // class IterativeSolver<...>

template<VectorLike InVector, VectorLike OutVector>
bool IterativeSolver<InVector, OutVector>::
                          Solve(InVector& xVec,
                                OutVector const& bVec,
                                Operator<InVector, OutVector> const& anyOp) {

  // ----------------------
  // Initialize the solver.
  // ----------------------
  if (PreOp != nullptr) {
    PreOp->Build(xVec, bVec, anyOp);
  }
  real_t const initialError =
    (AbsoluteError = Init(xVec, bVec, anyOp, PreOp.get()));
  std::cout << std::fixed << std::scientific << std::setprecision(15);
  std::cout << "\tI\t" << initialError << std::endl;
  if (AbsoluteTolerance > 0.0 && AbsoluteError < AbsoluteTolerance) {
    Finalize(xVec, bVec, anyOp, PreOp.get());
    return true;
  }

  // ----------------------
  // Iterate the solver:
  // ----------------------
  bool converged = false;
  for (Iteration = 0; !converged && (Iteration < NumIterations); ++Iteration) {
    AbsoluteError = Iterate(xVec, bVec, anyOp, PreOp.get());
    RelativeError = AbsoluteError/initialError;

    converged |=
      (AbsoluteTolerance > 0.0) && (AbsoluteError < AbsoluteTolerance);
    converged |=
      (RelativeTolerance > 0.0) && (RelativeError < RelativeTolerance);

    if (Iteration % 20 == 0 || converged) {
      std::cout << "\t" << (Iteration + 1) << "\t"
        << AbsoluteError << "\t" << RelativeError << std::endl;
    }
  }

  Finalize(xVec, bVec, anyOp, PreOp.get());
  std::cout << "\t----------------------" << std::endl;

  if (VerifySolution) {
    OutVector rVec;
    rVec.Assign(bVec, false);
    anyOp.Residual(rVec, bVec, xVec);
    real_t const
      trueAbsoluteError = rVec.Norm2(),
      trueRelativeError = trueAbsoluteError/initialError;
    std::cout << "\tT\t"
      << trueAbsoluteError << "\t" << trueRelativeError << std::endl;
    std::cout << "\t----------------------" << std::endl;
  }

  return converged;

} // IterativeSolver<...>::Solve

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract inner-outer iterative solver.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<VectorLike InVector, VectorLike OutVector = InVector>
class InnerOuterIterativeSolver : public IterativeSolver<InVector, OutVector> {
public:
  STORM_CLASS_(InnerOuterIterativeSolver, IterativeSolver<InVector, OutVector>)

  STORM_FIELD_(size_t, InnerIteration, {0})
  STORM_FIELD_(size_t, NumInnerIterations, {50})

protected:

  /// @brief Initialize the outer iterations.
  ///
  /// This function is used invoked only once, \
  ///   in the initialization phase.
  ///
  /// @param xVec Initial guess for the solution vector, 𝒙.
  /// @param bVec Right-hand-side vector, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  /// @param preOp Preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Residual norm of the initial guess, ‖𝒃 - 𝓐(𝒙)‖.
  virtual real_t OuterInit(InVector const& xVec,
                           OutVector const& bVec,
                           Operator<InVector, OutVector> const& anyOp,
                           Preconditioner<InVector> const* preOp) = 0;

  /// @brief Initialize the inner iterations.
  ///
  /// This function is invoked before the each inner iteration loop.
  ///
  /// @param xVec Solution vector, 𝒙.
  /// @param bVec Right-hand-side vector, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  /// @param preOp Preconditioner operator, 𝓟(𝒙).
  virtual void InnerInit(InVector const& xVec,
                         OutVector const& bVec,
                         Operator<InVector, OutVector> const& anyOp,
                         Preconditioner<InVector> const* preOp) {}

  /// @brief Perform the inner iteration.
  ///
  /// @param xVec Solution vector, 𝒙.
  /// @param bVec Right-hand-side vector, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  /// @param preOp Preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Residual norm, ‖𝒃 - 𝓐(𝒙)‖.
  virtual real_t InnerIterate(InVector& xVec,
                              OutVector const& bVec,
                              Operator<InVector, OutVector> const& anyOp,
                              Preconditioner<InVector> const* preOp) = 0;

  /// @brief Finalize the inner iterations.
  ///
  /// This function is called in order to finalize \
  ///   the inner iterations or if some stopping criterion is met.
  ///
  /// @param xVec Solution vector, 𝒙.
  /// @param bVec Right-hand-side vector, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  /// @param preOp Preconditioner operator, 𝓟(𝒙).
  virtual void InnerFinalize(InVector& xVec,
                             OutVector const& bVec,
                             Operator<InVector, OutVector> const& anyOp,
                             Preconditioner<InVector> const* preOp) {}

  /// @brief Finalize the outer iterations.
  ///
  /// This function is used invoked only once, \
  ///   when some stopping criterion is met.
  ///
  /// @param xVec Solution vector, 𝒙.
  /// @param bVec Right-hand-side vector, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  /// @param preOp Preconditioner operator, 𝓟(𝒙).
  virtual void OuterFinalize(InVector& xVec,
                             OutVector const& bVec,
                             Operator<InVector, OutVector> const& anyOp,
                             Preconditioner<InVector> const* preOp) {}

private:

  real_t Init(InVector const& xVec,
              OutVector const& bVec,
              Operator<InVector, OutVector> const& anyOp,
              Preconditioner<InVector> const* preOp) override final {
    return OuterInit(xVec, bVec, anyOp, preOp);
  }

  real_t Iterate(InVector& xVec,
                 OutVector const& bVec,
                 Operator<InVector, OutVector> const& anyOp,
                 Preconditioner<InVector> const* preOp) override final {
    InnerIteration = this->Iteration % NumInnerIterations;
    if (InnerIteration == 0) {
      InnerInit(xVec, bVec, anyOp, preOp);
    }
    real_t const residualNorm = InnerIterate(xVec, bVec, anyOp, preOp);
    if (InnerIteration == NumInnerIterations - 1) {
      InnerFinalize(xVec, bVec, anyOp, preOp);
    }
    return residualNorm;
  }

  void Finalize(InVector& xVec,
                OutVector const& bVec,
                Operator<InVector, OutVector> const& anyOp,
                Preconditioner<InVector> const* preOp) override final {
    if (InnerIteration != NumInnerIterations - 1) {
      InnerFinalize(xVec, bVec, anyOp, preOp);
    }
    OuterFinalize(xVec, bVec, anyOp, preOp);
  }

}; // class InnerOuterIterativeSolver<...>

/// ----------------------------------------------------------------- ///
/// @brief Solve the operator equation 𝓐(𝒙) = 𝒃, \
///   when 𝓐(𝒙) is the non-uniform operator (𝓐(𝟢) ≠ 𝟢).
/// ----------------------------------------------------------------- ///
template<VectorLike Vector>
bool SolveNonUniform(Solver<Vector>& solver,
                     Vector& xVec,
                     Vector const& bVec,
                     Operator<Vector> const& anyOp) {

  Vector zVec, fVec;

  zVec.Assign(xVec, false);
  fVec.Assign(bVec, false);

  // ----------------------
  // Solve an equation with the "uniformed" operator:
  // 𝓐(𝒙) - 𝓐(𝟢) = 𝒃 - 𝓐(𝟢).
  // ----------------------
  fVec.Fill(0.0);
  anyOp.MatVec(zVec, fVec);
  fVec.Sub(bVec, zVec);
  
  auto const uniOp = MakeOperator<Vector>(
    [&](Vector& yVec, Vector const& xVec) {
      anyOp.MatVec(yVec, xVec);
      yVec.SubAssign(zVec);
    });

  return solver.Solve(xVec, fVec, *uniOp);

} // SolveNonUniform<...>

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

}; // class PowerIterations<...>

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

} // PowerIterations<...>::EstimateLargestEigenvalue

} // namespace Storm
