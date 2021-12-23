/// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ///
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
/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ///
#ifndef _STORM_SOLVER_HXX_
#define _STORM_SOLVER_HXX_

#include <StormRuler_API.h>
#include <stormSolvers/stormOperator.hxx>
#include <stormSolvers/stormPreconditioner.hxx>

#include <iostream>
#include <stdexcept>

/// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ///
/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ///

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract operator equation 𝓐(𝒙) = 𝒃 solver.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class tInArray, class tOutArray = tInArray>
class stormSolver : public stormBaseObject {
public:

  /// @brief Solve the operator equation.
  ///
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  ///
  /// @returns Status of operation.
  virtual bool Solve(tInArray& xArr,
                     const tOutArray& bArr,
                     const stormOperator<tInArray, tOutArray>& anyOp) = 0;

}; // class stormSolver<...>

/// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ///
/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ///

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract operator equation 𝓟(𝓐(𝒙)) = 𝓟(𝒃) iterative solver.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class tInArray, class tOutArray = tInArray>
class stormIterativeSolver : public stormSolver<tInArray, tOutArray> {
private:
  stormSize_t NumIterations = 0;
  stormSize_t MaxNumIterations = 2000;
  stormReal_t AbsoluteError = 0.0, RelativeError = 0.0;
  stormReal_t AbsoluteTolerance = 1.0e-6, RelativeTolerance = 1.0e-6;

public:
  stormPreconditioner<tInArray>* PreOp = nullptr;

public:

  /// @brief Solve the operator equation.
  ///
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  ///
  /// @returns Status of operation.
  bool Solve(tInArray& xArr,
             const tOutArray& bArr,
             const stormOperator<tInArray, tOutArray>& anyOp) override final;

protected:

  /// @brief Initialize the iterative solver.
  ///
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  /// @param preOp Preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Residual norm-like value.
  virtual stormReal_t Init(tInArray& xArr,
                           const tOutArray& bArr,
                           const stormOperator<tInArray, tOutArray>& anyOp,
                           const stormPreconditioner<tInArray>* preOp) = 0;

  /// @brief Iterate the solver.
  ///
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  /// @param preOp Preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Residual norm-like value.
  virtual stormReal_t Iterate(tInArray& xArr,
                              const tOutArray& bArr,
                              const stormOperator<tInArray, tOutArray>& anyOp,
                              const stormPreconditioner<tInArray>* preOp) = 0;

  /// @brief Finalize the iterations.
  ///
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  /// @param preOp Preconditioner operator, 𝓟(𝒙).
  virtual void Finalize(tInArray& xArr,
                        const tOutArray& bArr,
                        const stormOperator<tInArray, tOutArray>& anyOp,
                        const stormPreconditioner<tInArray>* preOp) {}

}; // class stormIterativeSolver<...>

template<class tInArray, class tOutArray>
bool stormIterativeSolver<tInArray, tOutArray>::
                            Solve(tInArray& xArr,
                                  const tOutArray& bArr,
                                  const stormOperator<tInArray, tOutArray>& anyOp) {
  // ----------------------
  // Initialize the solver.
  // ----------------------
  if (PreOp != nullptr) {
    PreOp->Build(anyOp);
  }
  const stormReal_t initialError =
    (AbsoluteError = Init(xArr, bArr, anyOp, PreOp));
  std::cout << "\t1 " << initialError << std::endl;
  if (AbsoluteTolerance > 0.0 && AbsoluteError < AbsoluteTolerance) {
    Finalize(xArr, bArr, anyOp, PreOp);
    return true;
  }

  // ----------------------
  // Iterate the solver:
  // ----------------------
  for (NumIterations = 1; NumIterations <= MaxNumIterations; ++NumIterations) {
    AbsoluteError = Iterate(xArr, bArr, anyOp, PreOp);
    RelativeError = AbsoluteError/initialError;
    std::cout << "\t" << NumIterations << " "
      << AbsoluteError << " " << RelativeError << std::endl;

    if (AbsoluteTolerance > 0.0 && AbsoluteError < AbsoluteTolerance) {
      Finalize(xArr, bArr, anyOp, PreOp);
      return true;
    }
    if (RelativeTolerance > 0.0 && RelativeError < RelativeTolerance) {
      Finalize(xArr, bArr, anyOp, PreOp);
      return true;
    }
  }

  Finalize(xArr, bArr, anyOp, PreOp);
  return false;

} // stormIterativeSolver<...>::Solve

/// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ///
/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ///

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Largest eigenvalue estimator based on the Power Iterations.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class tArray>
class stormPowerIterations final : public stormBaseObject {
public:

  /// @brief Estimate the largest eigenvalue of \
  ///   the linear operator 𝓐 using the Power Iterations method.
  ///
  /// @param xArr On input: a non-zero vector that is used as \
  ///   the initial guess for the Power iterations; on output: \
  ///   estimate of the eigenvector, corresponding to the largest eigenvalue.
  /// @param linOp Linear operator, 𝓐(𝒙).
  /// @param maxIterations Maximum number of the iterations.
  /// @param relativeTolerance Relative error tolerance \
  ///   to terminate the iterations before the maximum number is reached.
  ///
  /// @returns Estimate the largest eigenvalue of 𝓐.
  static stormReal_t
    EstimateLargestEigenvalue(tArray& xArr,
                              const stormOperator<tArray>& linOp,
                              stormSize_t maxIterations = 20,
                              stormReal_t relativeTolerance = 1.0e-8);

}; // class stormPowerIterations<...>

template<class tArray>
stormReal_t stormPowerIterations<tArray>::
    EstimateLargestEigenvalue(tArray& xArr,
                              const stormOperator<tArray>& linOp,
                              stormSize_t maxIterations,
                              stormReal_t relativeTolerance) {

  stormArray yArr;
  stormUtils::AllocLike(xArr, yArr);

  // ----------------------
  // Initialize the Power Iterations:
  // 𝜆 ← 𝟣,
  // 𝒙 ← 𝘙𝘢𝘯𝘥𝘰𝘮(),
  // 𝒙 ← 𝒙/‖𝒙‖.
  // ----------------------
  stormReal_t lambda = 1.0;
  stormUtils::RandFill(xArr);
  stormUtils::Scale(xArr, xArr, 1.0/stormUtils::Norm2(xArr));

  for (stormSize_t iteration = 1; iteration <= maxIterations; ++iteration) {

    // ----------------------
    // Continue the Power Iterations:
    // 𝒚 ← 𝓐𝒙,
    // 𝜆̅ ← 𝜆, 𝜆 ← <𝒙⋅𝒚>,
    // 𝒙 ← 𝒚/‖𝒚‖.
    // ----------------------
    linOp.MatVec(yArr, xArr);
    const stormReal_t lambdaBar = lambda;
    lambda = stormUtils::Dot(xArr, yArr);
    stormUtils::Scale(xArr, yArr, 1.0/stormUtils::Norm2(yArr));

    // ----------------------
    // Check for the convergence on 𝜆 and 𝜆̅:
    // ----------------------
    if (std::abs((lambda - lambdaBar)/lambdaBar) < relativeTolerance) {
      break;
    }
  }

  return lambda;

} // stormPowerIterations<...>::EstimateLargestEigenvalue

/// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ///
/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ///

#endif // ifndef _STORM_SOLVER_HXX_
