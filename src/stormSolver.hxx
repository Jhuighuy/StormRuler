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
#include <stormOperator.hxx>

#include <iostream>

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
                           const stormOperator<tInArray, tInArray>* preOp = nullptr) = 0;

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
                              const stormOperator<tInArray, tInArray>* preOp = nullptr) = 0;
  
  /// @brief Finalize the iterations.
  ///
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  /// @param preOp Preconditioner operator, 𝓟(𝒙).
  virtual void Finalize(tInArray& xArr,
                        const tOutArray& bArr,
                        const stormOperator<tInArray, tOutArray>& anyOp,
                        const stormOperator<tInArray, tInArray>* preOp = nullptr) {}

}; // class stormIterativeSolver<...>

template<class tInArray, class tOutArray>
bool stormIterativeSolver<tInArray, tOutArray>::
                            Solve(tInArray& xArr,
                                  const tOutArray& bArr,
                                  const stormOperator<tInArray, tOutArray>& anyOp) {
  // ----------------------
  // Initialize the solver.
  // ----------------------
  const stormReal_t initialError = 
    (AbsoluteError = Init(xArr, bArr, anyOp));
  std::cout << "\t1 " << initialError << std::endl;
  if (AbsoluteTolerance > 0.0 && AbsoluteError < AbsoluteTolerance) {
    Finalize(xArr, bArr, anyOp);
    return true;
  }

  // ----------------------
  // Iterate the solver:
  // ----------------------
  for (NumIterations = 1; NumIterations < MaxNumIterations; ++NumIterations) {
    AbsoluteError = Iterate(xArr, bArr, anyOp);
    RelativeError = AbsoluteError/initialError;
    std::cout << "\t" << NumIterations << " " 
      << AbsoluteError << " " << RelativeError << std::endl;

    if (AbsoluteTolerance > 0.0 && AbsoluteError < AbsoluteTolerance) {
      Finalize(xArr, bArr, anyOp);
      return true;
    }
    if (RelativeTolerance > 0.0 && RelativeError < RelativeTolerance) {
      Finalize(xArr, bArr, anyOp);
      return true;
    }
  }

  Finalize(xArr, bArr, anyOp);
  return false;

} // stormIterativeSolver<...>::Solve

/// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ///
/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ///

#endif // ifndef _STORM_SOLVER_HXX_
