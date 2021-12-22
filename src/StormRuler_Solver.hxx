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
#ifndef _STORM_RULER_SOLVER_HXX_
#define _STORM_RULER_SOLVER_HXX_

#include <StormRuler_API.h>

#include <iostream>
#include <functional>

/// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ///
/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ///

#define _STORM_NOT_IMPLEMENTED_() do { \
  std::cerr << "not implemented" << std::endl; exit(1); } while(false)

template<class tArray>
class stormLinearOperator {
public:
  std::function<void(tArray&, const tArray&)> MatVec;
  std::function<void(tArray&, const tArray&)> ConjMatVec;
};

class stormBaseObject {
public:
  virtual ~stormBaseObject() = default;
};

class stormArray {
public:
  stormMesh_t Mesh;
  stormArray_t Array;
  bool Copy = false;

  ~stormArray() {
    if (Copy) stormFree(Array);
  }
};

namespace stormUtils {
  stormReal_t SafeDivide(stormReal_t x, stormReal_t y) {
    return (y == 0.0) ? 0.0 : (x/y);
  }

  stormReal_t Norm2(const stormArray& z) {
    return stormNorm2(z.Mesh, z.Array);
  }
  stormReal_t Dot(const stormArray& z, const stormArray& y) {
    return stormDot(z.Mesh, z.Array, y.Array);
  }

  void Set(stormArray& z, const stormArray& y) {
    stormSet(z.Mesh, z.Array, y.Array);
  }

  void Fill(stormArray& z, stormReal_t a) {
    stormFill(z.Mesh, z.Array, a);
  }

  void Scale(stormArray& z, const stormArray& y, stormReal_t a) {
    stormScale(z.Mesh, z.Array, y.Array, a);
  }

  void Add(stormArray& z, const stormArray& y, const stormArray& x, 
          stormReal_t a = 1.0, stormReal_t b = 1.0) {
    stormAdd(z.Mesh, z.Array, y.Array, x.Array, a, b);
  }
  void Sub(stormArray& z, const stormArray& y, const stormArray& x, 
          stormReal_t a = 1.0, stormReal_t b = 1.0) {
    stormSub(z.Mesh, z.Array, y.Array, x.Array, a, b);
  }

  void AllocLike(const stormArray& like, stormArray& z) {
    z.Mesh = like.Mesh;
    z.Array = stormAllocLike(like.Array);
    z.Copy = true;
  }
  template<class... tArray>
  void AllocLike(const stormArray& like, stormArray& z, tArray&... zz) {
    AllocLike(like, z);
    AllocLike(like, zz...);
  }
}

/// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ///
/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ///

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract operator equation 𝓐(𝒙) = 𝒃 solver.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class tArray, class tOperator>
class stormSolver : public stormBaseObject {
public:

  /// @brief Solve the operator equation. 
  ///
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  ///
  /// @returns Status of operation.
  virtual bool Solve(tArray& xArr,
                     const tArray& bArr,
                     const tOperator& anyOp) = 0;

}; // class stormSolver<...>

/// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ///
/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ///

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract operator equation 𝓟(𝓐(𝒙)) = 𝓟(𝒃) iterative solver.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class tArray, class tOperator>
class stormIterativeSolver : public stormSolver<tArray, tOperator> {
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
  bool Solve(tArray& xArr,
             const tArray& bArr,
             const tOperator& anyOp) override final;

protected:

  /// @brief Initialize the iterative solver.
  ///
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  /// @param preOp Preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Residual norm-like value.
  virtual stormReal_t Init(tArray& xArr,
                           const tArray& bArr,
                           const tOperator& anyOp,
                           const tOperator* preOp = nullptr) = 0;

  /// @brief Iterate the solver.
  ///
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  /// @param preOp Preconditioner operator, 𝓟(𝒙).
  ///
  /// @returns Residual norm-like value.
  virtual stormReal_t Iterate(tArray& xArr,
                              const tArray& bArr,
                              const tOperator& anyOp,
                              const tOperator* preOp = nullptr) = 0;
  
  /// @brief Finalize the iterations.
  ///
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param anyOp Equation operator, 𝓐(𝒙).
  /// @param preOp Preconditioner operator, 𝓟(𝒙).
  virtual void Finalize(tArray& xArr,
                        const tArray& bArr,
                        const tOperator& anyOp,
                        const tOperator* preOp = nullptr) {}

}; // class stormIterativeSolver<...>

template<class tArray, class tOperator>
bool stormIterativeSolver<tArray, tOperator>::Solve(tArray& xArr,
                                                    const tArray& bArr,
                                                    const tOperator& anyOp) {
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

#endif // ifndef _STORM_RULER_SOLVER_HXX_
