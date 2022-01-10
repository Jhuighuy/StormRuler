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
#ifndef _STORM_PRECONDITIONER_HXX_
#define _STORM_PRECONDITIONER_HXX_

#include <iostream>

#include <stormSolvers/stormOperator.hxx>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Preconditioner side.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
enum class stormPreconditionerSide {

  /// Left preconditioned equation is solved, 𝓟(𝓐(𝒙)) = 𝓟(𝒃).
  Left,

  /// Right preconditioned equation is solved, 𝓐(𝓟(𝒙̃)) = 𝒃, 𝓟(𝒙̃) = 𝒙.
  Right,

  /// Symmetric preconditioned equation is solved, \
  ///   𝓜(𝓐(𝓝(𝒙̃))) = 𝓜(𝒃), 𝓝(𝒙̃) = 𝒙, 𝓟 = 𝓝.
  Symmetric,

}; // enum class stormPreconditionerSide

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract preconditioner operator.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class tArray>
class stormPreconditioner : public stormOperator<tArray> {
public:

  /// @brief Build the preconditioner.
  ///
  /// @param xArr Solution (block-)array, 𝒙.
  /// @param bArr Right-hand-side (block-)array, 𝒃.
  /// @param anyOp Operator to build the preconditioner upon.
  virtual void Build(tArray const& xArr,
                     tArray const& bArr,
                     stormOperator<tArray> const& anyOp) {}

}; // class stormPreconditioner<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Identity preconditioner, \
///   intended to be used for debugging only.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class tArray>
class stormIdentityPreconditioner final : public stormPreconditioner<tArray> {
private:

  void MatVec(tArray& yArr,
              tArray const& xArr) const override {
    std::cout << "`stormIdentityPreconditioner<...>::MatVec`!" << std::endl;
    stormBlas::Set(yArr, xArr);
  }

  void ConjMatVec(tArray& xArr,
                  tArray const& yArr) const override {
    std::cout << "`stormIdentityPreconditioner<...>::ConjMatVec`!" << std::endl;
    stormBlas::Set(xArr, yArr);
  }

}; // class stormIdentityPreconditioner<...>

#endif // ifndef _STORM_PRECONDITIONER_HXX_
