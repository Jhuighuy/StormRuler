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

#include <stormBase.hxx>
#include <stormBlas/stormOperator.hxx>

_STORM_NAMESPACE_BEGIN_

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Preconditioner side.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
enum class PreconditionerSide {

  /// @brief Left preconditioned equation is solved, 𝓟𝓐𝒙 = 𝓟𝒃.
  ///
  /// When the left preconditioning is used, iterative solver tracks \
  ///   convergence by the left preconditioned residual norm, ‖𝓟(𝒃 - 𝓐𝒙)‖.
  Left,

  /// Right preconditioned equation is solved, 𝓐𝓟𝒙̃ = 𝒃, 𝓟𝒙̃ = 𝒙.
  ///
  /// When the right preconditioning is used, iterative solver tracks \
  ///   convergence by the unpreconditioned residual norm, ‖𝒃 - 𝓐𝒙‖.
  Right,

  /// Symmetric preconditioned equation is solved, \
  ///   𝓜𝓐𝓝𝒙̃ = 𝓜𝒃, 𝓝𝒙̃ = 𝒙, 𝓟 = 𝓜𝓝.
  ///
  /// When the symmetric preconditioning is used, iterative solver tracks \
  ///   convergence by the partially preconditioned residual norm, ‖𝓜(𝒃 - 𝓐𝒙)‖.
  Symmetric,

}; // enum class PreconditionerSide

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract preconditioner operator.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Vector>
class Preconditioner : public Operator<Vector> {
public:

  /// @brief Build the preconditioner.
  ///
  /// @param xVec Solution vector, 𝒙.
  /// @param bVec Right-hand-side vector, 𝒃.
  /// @param anyOp Operator to build the preconditioner upon.
  virtual void Build(Vector const& xVec,
                     Vector const& bVec,
                     Operator<Vector> const& anyOp) {}

}; // class Preconditioner<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Identity preconditioner, \
///   intended to be used for debugging only.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Vector>
class IdentityPreconditioner final : public Preconditioner<Vector> {
private:

  void MatVec(Vector& yVec,
              Vector const& xVec) const override {
    std::cout << "`IdentityPreconditioner<...>::MatVec`!" << std::endl;
    Blas::Set(yVec, xVec);
  }

  void ConjMatVec(Vector& xVec,
                  Vector const& yVec) const override {
    std::cout << "`IdentityPreconditioner<...>::ConjMatVec`!" << std::endl;
    Blas::Set(xVec, yVec);
  }

}; // class IdentityPreconditioner<...>

_STORM_NAMESPACE_END_

#endif // ifndef _STORM_PRECONDITIONER_HXX_
