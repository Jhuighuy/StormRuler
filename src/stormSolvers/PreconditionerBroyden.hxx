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
#include <stormBlas/stormSubspace.hxx>
#include <stormSolvers/Preconditioner.hxx>

namespace Storm {

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief @c Broyden's method preconditioner.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Vector>
class BroydenPreconditioner final : public Preconditioner<Vector> {
public:
  size_t Rank = 10;

private:
  real_t omega_ = 1.0*7800.0;
  Subspace<Vector> uVecs_, vVecs_;

  void Build(Vector const& xVec,
             Vector const& bVec,
             Operator<Vector> const& linOp) override;

  void MatVec(Vector& yVec,
              Vector const& xVec) const override;

  void MatVec(size_t k,
              Vector& yVec,
              Vector const& xVec) const;

}; // class BroydenPreconditioner<...>

template<class Vector>
void BroydenPreconditioner<Vector>::Build(Vector const& xVec,
                                          Vector const& bVec,
                                          Operator<Vector> const& linOp) {

  size_t const m = Rank;

  Vector rVec, yVec;

  rVec.Assign(xVec, false);
  yVec.Assign(xVec, false);
  uVecs_.Assign(m, xVec, false);
  vVecs_.Assign(m, xVec, false);

  // ----------------------
  // Compute the initial residual:
  // 𝒓 ← 𝒃 - 𝓐𝒙.
  // ----------------------
  linOp.Residual(rVec, bVec, xVec);

  for (size_t k = 0; k < m; ++k) {

    // ----------------------
    // Build the new vector 𝒗ₖ and update the residual:
    // 𝗶𝗳 𝑘 = 𝟢:
    //   𝒗ₖ ← 𝒓/𝜔,
    // 𝗲𝗹𝘀𝗲:
    //   𝒗ₖ ← (𝓗ₖ₋₁)𝒓,
    // 𝗲𝗻𝗱 𝗶𝗳
    // 𝒖ₖ ← 𝓐𝒗ₖ,
    // 𝒓 ← 𝒓 - 𝒖ₖ.
    // ----------------------
    if (k == 0) {
      Blas::Scale(vVecs_(k), rVec, 1.0/omega_);
    } else {
      MatVec(k - 1, vVecs_(k), rVec);
    }
    linOp.MatVec(uVecs_(k), vVecs_(k));
    Blas::Sub(rVec, rVec, uVecs_(k));

    // ----------------------
    // Build the new vector 𝒖ₖ and matrix 𝓗ₖ:
    // 𝗶𝗳 𝑘 = 𝟢:
    //   𝒖ₖ ← 𝒖ₖ/𝜔,
    // 𝗲𝗹𝘀𝗲:
    //   𝒖ₖ ← (𝓗ₖ₋₁)𝒖ₖ,
    // 𝗲𝗻𝗱 𝗶𝗳
    // 𝒖ₖ ← (𝒗ₖ - 𝒖ₖ)/<𝒗ₖ⋅𝒖ₖ>.
    // ----------------------
    if (k != 0) {
      MatVec(k - 1, uVecs_(k), uVecs_(k));
    } else {
      Blas::Scale(uVecs_(k), uVecs_(k), 1.0/omega_);
    }
    real_t const alpha = 1.0/Blas::Dot(uVecs_(k), vVecs_(k));
    Blas::Sub(uVecs_(k), vVecs_(k), alpha, uVecs_(k), alpha);
  }

} // BroydenPreconditioner<...>::Build

template<class Vector>
void BroydenPreconditioner<Vector>::MatVec(Vector& yVec,
                                           Vector const& xVec) const {

  size_t const m = Rank;

  // ----------------------
  // Compute the inverse approximation product:
  // 𝒚 ← (𝓗ₘ₋₁)𝒙.
  // ----------------------
  if (m == 0) {
    Blas::Scale(yVec, xVec, 1.0/omega_);
  } else {
    MatVec(m - 1, yVec, xVec);
  }

} // BroydenPreconditioner<...>::MatVec

template<class Vector>
void BroydenPreconditioner<Vector>::MatVec(size_t k,
                                           Vector& yVec,
                                           Vector const& xVec) const {

  // ----------------------
  // Compute the (𝓗ₖ₋₁)𝒙: 
  // 𝗶𝗳 𝑘 = 𝟢:
  //   𝒚 ← 𝒙/𝜔.
  // 𝗲𝗹𝘀𝗲:
  //   𝒚 ← (𝓗ₖ₋₁)𝒙.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  if (k == 0) {
    Blas::Scale(yVec, xVec, 1.0/omega_);
  } else {
    MatVec(k - 1, yVec, xVec);
  }

  // ----------------------
  // Compute the rank-1 correction for (𝓗ₖ)𝒙:
  // 𝒚 ← 𝒚 + <𝒚⋅𝒗ₖ>⋅𝒖ₖ.
  // ----------------------
  Blas::Add(yVec, yVec, uVecs_(k), Blas::Dot(vVecs_(k), yVec));

} // BroydenPreconditioner<...>::MatVec

} // namespace Storm
