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

#include <iostream>

#include <Storm/Base.hpp>

#include <Storm/Utils/Enum.hpp>

#include <Storm/Solvers/Operator.hpp>

namespace Storm {

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Preconditioner side.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
class PreconditionerSide final : public Enum<PreconditionerSide> {
  // clang-format off

  StormEnum_(PreconditionerSide)

  /// @brief Left preconditioned equation is solved, 𝓟𝓐𝒙 = 𝓟𝒃.
  ///
  /// When the left preconditioning is used, iterative solver tracks 
  ///   convergence by the left preconditioned residual norm, ‖𝓟(𝒃 - 𝓐𝒙)‖.
  StormEnumValue_(Left)

  /// Right preconditioned equation is solved, 𝓐𝓟𝒙̃ = 𝒃, 𝓟𝒙̃ = 𝒙.
  ///
  /// When the right preconditioning is used, iterative solver tracks 
  ///   convergence by the unpreconditioned residual norm, ‖𝒃 - 𝓐𝒙‖.
  StormEnumValue_(Right)

  /// Symmetric preconditioned equation is solved, 
  ///   𝓜𝓐𝓝𝒙̃ = 𝓜𝒃, 𝓝𝒙̃ = 𝒙, 𝓟 = 𝓜𝓝.
  ///
  /// When the symmetric preconditioning is used, iterative solver tracks 
  ///   convergence by the partially preconditioned residual norm, ‖𝓜(𝒃 - 𝓐𝒙)‖.
  StormEnumValue_(Symmetric)

  // clang-format on

}; // enum class PreconditionerSide

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract preconditioner operator.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<VectorLike Vector>
class Preconditioner : public Operator<Vector> {
public:

  /// @brief Build the preconditioner.
  ///
  /// @param x_vec Solution vector, 𝒙.
  /// @param b_vec Right-hand-side vector, 𝒃.
  /// @param any_op Operator to build the preconditioner upon.
  virtual void Build(const Vector& x_vec, const Vector& b_vec,
                     const Operator<Vector>& any_op) {}

}; // class Preconditioner

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Identity preconditioner,
///   intended to be used for debugging only.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<VectorLike Vector>
class IdentityPreconditioner final : public Preconditioner<Vector> {
private:

  void mul(Vector& y_vec, const Vector& x_vec) const override {
    std::clog << "IdentityPreconditioner::mul called" << std::endl;
    y_vec <<= x_vec;
  }

  void conj_mul(Vector& x_vec, const Vector& y_vec) const override {
    std::clog << "IdentityPreconditioner::conj_mul called" << std::endl;
    x_vec <<= y_vec;
  }

}; // class IdentityPreconditioner

} // namespace Storm
