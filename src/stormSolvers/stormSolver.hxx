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
#ifndef _STORM_SOLVER_HXX_
#define _STORM_SOLVER_HXX_

#include <iomanip>
#include <iostream>
#include <stdexcept>

#include <stormSolvers/stormOperator.hxx>
#include <stormSolvers/stormPreconditioner.hxx>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract operator equation solver.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class InVector, class OutVector = InVector>
class stormSolver : public stormBaseObject {
public:

  /// @brief Solve the operator equation 𝓐(𝒙) = 𝒃.
  ///
  /// @param xVec Solution vector, 𝒙.
  /// @param bVec Right-hand-side vector, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  ///
  /// @returns Status of operation.
  virtual bool Solve(InVector& xVec,
                     OutVector const& bVec,
                     stormOperator<InVector, OutVector> const& anyOp) = 0;

}; // class stormSolver<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract operator equation iterative solver.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class InVector, class OutVector = InVector>
class stormIterativeSolver : public stormSolver<InVector, OutVector> {
public:
  stormSize_t Iteration = 0;
  stormSize_t NumIterations = 2000;
  stormReal_t AbsoluteError = 0.0, RelativeError = 0.0;
  stormReal_t AbsoluteTolerance = 1.0e-6, RelativeTolerance = 1.0e-6;
  bool VerifySolution = false;

public:
  stormPreconditionerSide PreSide = stormPreconditionerSide::Right;
  std::unique_ptr<stormPreconditioner<InVector>> PreOp = nullptr;

protected:

  /// @brief Initialize the iterative solver.
  ///
  /// @param xVec Solution vector, 𝒙.
  /// @param bVec Right-hand-side vector, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  /// @param preOp Preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Residual norm, ‖𝒃 - 𝓐(𝒙)‖.
  virtual stormReal_t Init(InVector& xVec,
                           OutVector const& bVec,
                           stormOperator<InVector, OutVector> const& anyOp,
                           stormPreconditioner<InVector> const* preOp) = 0;

  /// @brief Iterate the solver.
  ///
  /// @param xVec Solution vector, 𝒙.
  /// @param bVec Right-hand-side vector, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  /// @param preOp Preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Residual norm, ‖𝒃 - 𝓐(𝒙)‖.
  virtual stormReal_t Iterate(InVector& xVec,
                              OutVector const& bVec,
                              stormOperator<InVector, OutVector> const& anyOp,
                              stormPreconditioner<InVector> const* preOp) = 0;

  /// @brief Finalize the iterations.
  ///
  /// @param xVec Solution vector, 𝒙.
  /// @param bVec Right-hand-side vector, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  /// @param preOp Preconditioner operator, 𝓟(𝒙).
  virtual void Finalize(InVector& xVec,
                        OutVector const& bVec,
                        stormOperator<InVector, OutVector> const& anyOp,
                        stormPreconditioner<InVector> const* preOp) {}

public:

  bool Solve(InVector& xVec,
             OutVector const& bVec,
             stormOperator<InVector, OutVector> const& anyOp) override final;

}; // class stormIterativeSolver<...>

template<class InVector, class OutVector>
bool stormIterativeSolver<InVector, OutVector>::
                            Solve(InVector& xVec,
                                  OutVector const& bVec,
                                  stormOperator<InVector, OutVector> const& anyOp) {
  // ----------------------
  // Initialize the solver.
  // ----------------------
  if (PreOp != nullptr) {
    PreOp->Build(xVec, bVec, anyOp);
  }
  stormReal_t const initialError =
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
  for (Iteration = 0; Iteration < NumIterations; ++Iteration) {
    AbsoluteError = Iterate(xVec, bVec, anyOp, PreOp.get());
    RelativeError = AbsoluteError/initialError;
    std::cout << "\t" << (Iteration + 1) << "\t"
      << AbsoluteError << "\t" << RelativeError << std::endl;

    if (AbsoluteTolerance > 0.0 && AbsoluteError < AbsoluteTolerance) {
      converged = true;
      break;
    }
    if (RelativeTolerance > 0.0 && RelativeError < RelativeTolerance) {
      converged = true;
      break;
    }
  }

  Finalize(xVec, bVec, anyOp, PreOp.get());
  std::cout << "\t----------------------" << std::endl;

  if (VerifySolution) {
    OutVector rVec;
    rVec.Assign(bVec, false);
    anyOp.MatVec(rVec, xVec);
    stormBlas::Sub(rVec, bVec, rVec);
    stormReal_t const 
      trueAbsoluteError = stormBlas::Norm2(rVec),
      trueRelativeError = trueAbsoluteError/initialError; 
    std::cout << "\tT\t"
      << trueAbsoluteError << "\t" << trueRelativeError << std::endl;
    std::cout << "\t----------------------" << std::endl;
  }

  return converged;

} // stormIterativeSolver<...>::Solve

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract inner-outer iterative solver.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class InVector, class OutVector = InVector>
class stormInnerOuterIterativeSolver : 
    public stormIterativeSolver<InVector, OutVector> {
public:
  stormSize_t InnerIteration = 0;
  stormSize_t NumInnerIterations = 50;

protected:

  /// @brief Initialize the outer iterations.
  /// 
  /// This function is used invoked only once, \
  ///   in the initialization phase.
  /// 
  /// @param xVec Solution vector, 𝒙.
  /// @param bVec Right-hand-side vector, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  /// @param preOp Preconditioner operator, 𝓟(𝒙).
  virtual void OuterInit(InVector& xVec,
                         OutVector const& bVec,
                         stormOperator<InVector, OutVector> const& anyOp,
                         stormPreconditioner<InVector> const* preOp) = 0;

  /// @brief Initialize the inner iterations.
  ///
  /// This function is invoked before the each inner iteration loop.
  ///
  /// @param xVec Solution vector, 𝒙.
  /// @param bVec Right-hand-side vector, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  /// @param preOp Preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Residual norm, ‖𝒃 - 𝓐𝒙‖.
  virtual stormReal_t InnerInit(InVector& xVec,
                                OutVector const& bVec,
                                stormOperator<InVector, OutVector> const& anyOp,
                                stormPreconditioner<InVector> const* preOp) = 0;

  /// @brief Perform the inner iteration.
  ///
  /// @param xVec Solution vector, 𝒙.
  /// @param bVec Right-hand-side vector, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  /// @param preOp Preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Residual norm, ‖𝒃 - 𝓐𝒙‖.
  virtual stormReal_t InnerIterate(InVector& xVec,
                                   OutVector const& bVec,
                                   stormOperator<InVector, OutVector> const& anyOp,
                                   stormPreconditioner<InVector> const* preOp) = 0;

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
                             stormOperator<InVector, OutVector> const& anyOp,
                             stormPreconditioner<InVector> const* preOp) {}

  
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
                             stormOperator<InVector, OutVector> const& anyOp,
                             stormPreconditioner<InVector> const* preOp) {}

private:

  stormReal_t Init(InVector& xVec,
                   OutVector const& bVec,
                   stormOperator<InVector, OutVector> const& anyOp,
                   stormPreconditioner<InVector> const* preOp) override final {
    OuterInit(xVec, bVec, anyOp, preOp);
    return InnerInit(xVec, bVec, anyOp, preOp);
  }

  stormReal_t Iterate(InVector& xVec,
                      OutVector const& bVec,
                      stormOperator<InVector, OutVector> const& anyOp,
                      stormPreconditioner<InVector> const* preOp) override final {
    InnerIteration = this->Iteration % NumInnerIterations;
    if (InnerIteration == 0 && this->Iteration != 0) {
      InnerIteration = NumInnerIterations - 1;
      InnerFinalize(xVec, bVec, anyOp, preOp);
      InnerIteration = 0;
      InnerInit(xVec, bVec, anyOp, preOp);
    }
    return InnerIterate(xVec, bVec, anyOp, preOp);
  }

  void Finalize(InVector& xVec,
                OutVector const& bVec,
                stormOperator<InVector, OutVector> const& anyOp,
                stormPreconditioner<InVector> const* preOp) override final {
    InnerFinalize(xVec, bVec, anyOp, preOp);
    OuterFinalize(xVec, bVec, anyOp, preOp);
  }

}; // class stormInnerOuterIterativeSolver<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Largest eigenvalue estimator based on the Power Iterations.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class tArray>
class stormPowerIterations final : public stormBaseObject {
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
  static stormReal_t
    EstimateLargestEigenvalue(tArray& xVec,
                              stormOperator<tArray> const& linOp,
                              stormSize_t maxIterations = 20,
                              stormReal_t relativeTolerance = 1.0e-8);

}; // class stormPowerIterations<...>

template<class tArray>
stormReal_t stormPowerIterations<tArray>::
    EstimateLargestEigenvalue(tArray& xVec,
                              stormOperator<tArray> const& linOp,
                              stormSize_t maxIterations,
                              stormReal_t relativeTolerance) {

  stormArray yVec;
  yVec.Assign(xVec, false);

  // ----------------------
  // Initialize the Power Iterations:
  // 𝜆 ← 𝟣,
  // 𝒙 ← 𝘙𝘢𝘯𝘥𝘰𝘮(),
  // 𝒙 ← 𝒙/‖𝒙‖.
  // ----------------------
  stormReal_t lambda = 1.0;
  stormBlas::RandFill(xVec);
  stormBlas::Scale(xVec, xVec, 1.0/stormBlas::Norm2(xVec));

  for (stormSize_t iteration = 0; iteration < maxIterations; ++iteration) {

    // ----------------------
    // Continue the Power Iterations:
    // 𝒚 ← 𝓐𝒙,
    // 𝜆̅ ← 𝜆, 𝜆 ← <𝒙⋅𝒚>,
    // 𝒙 ← 𝒚/‖𝒚‖.
    // ----------------------
    linOp.MatVec(yVec, xVec);
    stormReal_t const lambdaBar = lambda;
    lambda = stormBlas::Dot(xVec, yVec);
    stormBlas::Scale(xVec, yVec, 1.0/stormBlas::Norm2(yVec));

    // ----------------------
    // Check for the convergence on 𝜆 and 𝜆̅:
    // ----------------------
    if (std::abs((lambda - lambdaBar)/lambdaBar) < relativeTolerance) {
      break;
    }
  }

  return lambda;

} // stormPowerIterations<...>::EstimateLargestEigenvalue

#endif // ifndef _STORM_SOLVER_HXX_
