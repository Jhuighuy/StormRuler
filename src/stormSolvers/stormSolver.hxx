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
template<class tInArray, class tOutArray = tInArray>
class stormSolver : public stormBaseObject {
public:

  /// @brief Solve the operator equation 𝓐(𝒙) = 𝒃.
  ///
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  ///
  /// @returns Status of operation.
  virtual bool Solve(tInArray& xArr,
                     tOutArray const& bArr,
                     stormOperator<tInArray, tOutArray> const& anyOp) = 0;

}; // class stormSolver<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract operator equation iterative solver.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class tInArray, class tOutArray = tInArray>
class stormIterativeSolver : public stormSolver<tInArray, tOutArray> {
public:
  stormSize_t Iteration = 0;
  stormSize_t NumIterations = 2000;
  stormReal_t AbsoluteError = 0.0, RelativeError = 0.0;
  stormReal_t AbsoluteTolerance = 1.0e-6, RelativeTolerance = 1.0e-6;

public:
  std::unique_ptr<stormPreconditioner<tInArray>> PreOp = nullptr;

protected:

  /// @brief Initialize the iterative solver.
  ///
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  /// @param preOp Preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Residual norm, ‖𝒃 - 𝓐𝒙‖.
  virtual stormReal_t Init(tInArray& xArr,
                           tOutArray const& bArr,
                           stormOperator<tInArray, tOutArray> const& anyOp,
                           stormPreconditioner<tInArray> const* preOp) = 0;

  /// @brief Iterate the solver.
  ///
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  /// @param preOp Preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Residual norm, ‖𝒃 - 𝓐𝒙‖.
  virtual stormReal_t Iterate(tInArray& xArr,
                              tOutArray const& bArr,
                              stormOperator<tInArray, tOutArray> const& anyOp,
                              stormPreconditioner<tInArray> const* preOp) = 0;

  /// @brief Finalize the iterations.
  ///
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  /// @param preOp Preconditioner operator, 𝓟(𝒙).
  virtual void Finalize(tInArray& xArr,
                        tOutArray const& bArr,
                        stormOperator<tInArray, tOutArray> const& anyOp,
                        stormPreconditioner<tInArray> const* preOp) {}

public:

  bool Solve(tInArray& xArr,
             tOutArray const& bArr,
             stormOperator<tInArray, tOutArray> const& anyOp) override final;

}; // class stormIterativeSolver<...>

template<class tInArray, class tOutArray>
bool stormIterativeSolver<tInArray, tOutArray>::
                            Solve(tInArray& xArr,
                                  tOutArray const& bArr,
                                  stormOperator<tInArray, tOutArray> const& anyOp) {
  // ----------------------
  // Initialize the solver.
  // ----------------------
  if (PreOp != nullptr) {
    PreOp->Build(anyOp);
  }
  stormReal_t const initialError =
    (AbsoluteError = Init(xArr, bArr, anyOp, PreOp.get()));
  std::cout << std::fixed << std::scientific << std::setprecision(15);
  std::cout << "\tI\t" << initialError << std::endl;
  if (AbsoluteTolerance > 0.0 && AbsoluteError < AbsoluteTolerance) {
    Finalize(xArr, bArr, anyOp, PreOp.get());
    return true;
  }

  // ----------------------
  // Iterate the solver:
  // ----------------------
  bool converged = false;
  for (Iteration = 0; Iteration < NumIterations; ++Iteration) {
    AbsoluteError = Iterate(xArr, bArr, anyOp, PreOp.get());
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

  Finalize(xArr, bArr, anyOp, PreOp.get());
  std::cout << "\t----------------------" << std::endl;
  return converged;

} // stormIterativeSolver<...>::Solve

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract restartable iterative solver.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class tInArray, class tOutArray = tInArray>
class stormRestartableSolver : public stormIterativeSolver<tInArray, tOutArray> {
public:
  stormSize_t NumIterationsBeforeRestart = 50;

protected:

  /// @brief Preinitialize the iterative solver \
  ///   (primarily allocate the intermediate data structures).
  /// 
  /// This function is used invoked only once, \
  ///   in the initialization phase.
  /// 
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param hasPreOp True if the preconditioned solver is used. 
  virtual void PreInit(tInArray& xArr,
                       tOutArray const& bArr, 
                       bool hasPreOp) = 0;

  /// @brief Reinitialize the iterative solver after the restart.
  ///
  /// This function is invoked in the initialization phase \
  ///   and after the each restart.
  ///
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  /// @param preOp Preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Residual norm-like value.
  virtual stormReal_t ReInit(tInArray& xArr,
                             tOutArray const& bArr,
                             stormOperator<tInArray, tOutArray> const& anyOp,
                             stormPreconditioner<tInArray> const* preOp) = 0;

  /// @brief Iterate the solver between the restarts.
  ///
  /// @param k Current iteration number between the restarts.
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  /// @param preOp Preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Residual norm-like value.
  virtual stormReal_t ReIterate(stormSize_t k,
                                tInArray& xArr,
                                tOutArray const& bArr,
                                stormOperator<tInArray, tOutArray> const& anyOp,
                                stormPreconditioner<tInArray> const* preOp) = 0;

  /// @brief Finalize the iterations before the restart.
  ///
  /// This function is called in order to finalize \
  ///   the intermediate computation results before the restart happens \
  ///   or some stopping criterion is met.
  ///
  /// @param k Last iteration number.
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  /// @param preOp Preconditioner operator, 𝓟(𝒙).
  virtual void ReFinalize(stormSize_t k,
                          tInArray& xArr,
                          tOutArray const& bArr,
                          stormOperator<tInArray, tOutArray> const& anyOp,
                          stormPreconditioner<tInArray> const* preOp) {}

private:

  stormReal_t Init(tInArray& xArr,
                   tOutArray const& bArr,
                   stormOperator<tInArray, tOutArray> const& anyOp,
                   stormPreconditioner<tInArray> const* preOp) override final {
    PreInit(xArr, bArr, preOp != nullptr);
    return ReInit(xArr, bArr, anyOp, preOp);
  }

  stormReal_t Iterate(tInArray& xArr,
                      tOutArray const& bArr,
                      stormOperator<tInArray, tOutArray> const& anyOp,
                      stormPreconditioner<tInArray> const* preOp) override final {
    const stormSize_t k = this->Iteration % NumIterationsBeforeRestart;
    if (k == 0 && this->Iteration != 0) {
      ReFinalize(NumIterationsBeforeRestart - 1, xArr, bArr, anyOp, preOp);
      ReInit(xArr, bArr, anyOp, preOp);
    }
    return ReIterate(k, xArr, bArr, anyOp, preOp);
  }

  void Finalize(tInArray& xArr,
                tOutArray const& bArr,
                stormOperator<tInArray, tOutArray> const& anyOp,
                stormPreconditioner<tInArray> const* preOp) override final {
    const stormSize_t k = this->Iteration % NumIterationsBeforeRestart;
    ReFinalize(k, xArr, bArr, anyOp, preOp);
  }

}; // class stormRestartableSolver<...>

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
                              stormOperator<tArray> const& linOp,
                              stormSize_t maxIterations = 20,
                              stormReal_t relativeTolerance = 1.0e-8);

}; // class stormPowerIterations<...>

template<class tArray>
stormReal_t stormPowerIterations<tArray>::
    EstimateLargestEigenvalue(tArray& xArr,
                              stormOperator<tArray> const& linOp,
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

  for (stormSize_t iteration = 0; iteration < maxIterations; ++iteration) {

    // ----------------------
    // Continue the Power Iterations:
    // 𝒚 ← 𝓐𝒙,
    // 𝜆̅ ← 𝜆, 𝜆 ← <𝒙⋅𝒚>,
    // 𝒙 ← 𝒚/‖𝒚‖.
    // ----------------------
    linOp.MatVec(yArr, xArr);
    stormReal_t const lambdaBar = lambda;
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

#endif // ifndef _STORM_SOLVER_HXX_
