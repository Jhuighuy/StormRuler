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

#include <stormBase.hxx>

namespace Storm {

// clang-format off

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Vector-like concept.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Vector>
concept VectorLike = 
  requires(Vector& zVec, Vector const& yVec, bool copy) {
    // Require assignment:
    // 𝒛 ← 𝒚.
    { zVec.Assign(yVec) };
    { zVec.Assign(yVec, copy) };
  } &&
  requires(Vector const& zVec, Vector const& yVec) {
    // Require norm and operator computation:
    // <𝒛⋅𝒚> → ℙ, ‖𝒛‖ → ℝ.
    { zVec.Dot(yVec) } -> std::same_as<real_t>;
    { zVec.Norm2() } -> std::same_as<real_t>;
} &&
requires(Vector& zVec, Vector const& yVec, Vector const& xVec) {
  { zVec.Fill(0.0) };
  { zVec.RandFill() };
  { zVec.Scale(yVec, 1.0) };
  { zVec.ScaleAssign(1.0) };
};

// clang-format on

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Operator-like concept.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Operator, class InVector, class OutVector = InVector>
concept OperatorLike = requires(Operator& anyOp, OutVector& yVec,
                                InVector const& xVec) {
  anyOp(yVec, xVec);
};

} // namespace Storm
