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
#include <stormBlas/stormTensor.hxx>
#include <stormBlas/stormSubspace.hxx>
#include <stormSolvers/Preconditioner.hxx>

namespace Storm {

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief @c Broyden's method preconditioner.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<VectorLike Vector>
class BroydenPreconditioner final : public Preconditioner<Vector> {
public:
  size_t Rank = 40;

private:
  real_t omega_ = 1.0+0.75*8000.0, sigma_;
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

template<VectorLike Vector>
void BroydenPreconditioner<Vector>::Build(Vector const& xVec,
                                          Vector const& bVec,
                                          Operator<Vector> const& linOp) {

#if 0
  size_t const m = Rank;

  Vector rVec, qVec, yVec;

  rVec.Assign(xVec, false);
  qVec.Assign(xVec, false);
  yVec.Assign(xVec, false);
  uVecs_.Assign(m, xVec, false);
  vVecs_.Assign(m, xVec, false);

  // ----------------------
  // Compute the initial residual:
  // 𝒚 ← 𝒙,
  // 𝒓 ← 𝒃 - 𝓐𝒙.
  // ----------------------
  Blas::Set(yVec, xVec);
  //Blas::RandFill(yVec);
  linOp.Residual(rVec, bVec, yVec);

  std::cout << "Building LGB:" << std::endl;
  std::cout << Blas::Norm2(rVec) << std::endl;

  //real_t hui = 2.0/8000.0;
  for (size_t k = 0; k < m; ++k) {

    // ----------------------
    // 𝒗ₖ ← (𝓗ₖ)𝒓,
    // 𝒒 ← 𝓐𝒗ₖ,
    // 𝒖ₖ ← (𝓗ₖ)𝒒,
    // ----------------------
    MatVec(k, vVecs_(k), rVec);
    linOp.MatVec(qVec, vVecs_(k));
    MatVec(k, uVecs_(k), qVec);

    // ----------------------
    // 𝜏 ← <𝒗ₖ⋅𝒗ₖ>/<𝒗ₖ⋅𝒖ₖ>,
    // 𝒚 ← 𝒚 + 𝜏⋅𝒗ₖ,
    // 𝒓 ← 𝒓 - 𝜏⋅𝒒.
    // ----------------------
    //real_t const tau = Blas::Dot(vVecs_(k), vVecs_(k))/Blas::Dot(vVecs_(k), uVecs_(k));
    real_t const tau = Blas::Dot(rVec, qVec)/Blas::Dot(qVec, qVec);
    //real_t const tau = Blas::Dot(vVecs_(k), uVecs_(k))/Blas::Dot(uVecs_(k), uVecs_(k));
    Blas::Add(yVec, yVec, vVecs_(k), tau);
    Blas::Sub(rVec, rVec, qVec, tau);

    std::cout << k << ' ' << Blas::Norm2(rVec) << std::endl;

    // ----------------------
    // 𝒖ₖ ← (𝒗ₖ - 𝒖ₖ)/<𝒗ₖ⋅𝒖ₖ>.
    // ----------------------
    real_t const alpha = 1.0/Blas::Dot(uVecs_(k), vVecs_(k));
    Blas::Sub(uVecs_(k), vVecs_(k), alpha, uVecs_(k), alpha);

  }

  //auto leftPreOp = MakeOperator<Vector>(
  //  [&](Vector& sVec, Vector const& tVec) {
  //    this->Operator<Vector>::MatVec(sVec, rVec, linOp, tVec);
  //    //linOp.MatVec(sVec, rVec, *this, tVec);
  //  });
  //PowerIterations<Vector> powerIterations;
  //sigma_ = 1.0;
  //sigma_ =
  //  powerIterations.EstimateLargestEigenvalue(qVec, *leftPreOp/**this*/, 200);
  //std::cout << "My lambda max = " << sigma_ << " " << 1.0/sigma_ << std::endl;
  //abort();
#endif
  abort();

} // BroydenPreconditioner<...>::Build

template<VectorLike Vector>
void BroydenPreconditioner<Vector>::MatVec(Vector& yVec,
                                           Vector const& xVec) const {

  size_t const m = Rank;

  // ----------------------
  // Compute the inverse approximation product:
  // 𝒚 ← (𝓗ₘ)𝒙.
  // ----------------------
  MatVec(m, yVec, xVec);
  //Blas::Scale(yVec, yVec, omega_);

} // BroydenPreconditioner<...>::MatVec

template<VectorLike Vector>
void BroydenPreconditioner<Vector>::MatVec(size_t k,
                                           Vector& yVec,
                                           Vector const& xVec) const {

#if 0
  Blas::Scale(yVec, xVec, 1.0/omega_);
  for (size_t i = 0; i < k; ++i) {
    Blas::Add(yVec, yVec, uVecs_(i), Blas::Dot(vVecs_(i), yVec));
  }

#if 0
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
#endif
#endif

} // BroydenPreconditioner<...>::MatVec

template<VectorLike Vector>
class BfgsPreconditioner final : public Preconditioner<Vector> {
public:
  size_t Rank = 20;

private:
  real_t omega_ = 1.0, sigma_;
  stormVector<real_t> rho_;
  Subspace<Vector> sVecs_, yVecs_;

  void Build(Vector const& xVec,
             Vector const& bVec,
             Operator<Vector> const& linOp) override;

  void MatVec(Vector& yVec,
              Vector const& xVec) const override;

  void MatVec(size_t k,
              Vector& yVec,
              Vector const& xVec) const;

}; // class BfgsPreconditioner<...>

template<VectorLike Vector>
void BfgsPreconditioner<Vector>::Build(Vector const& xVec,
                                       Vector const& bVec,
                                       Operator<Vector> const& linOp) {

#if 0
  size_t const m = Rank;

  real_t alpha;
  Vector uVec, rVec, pVec, qVec;

  rho_.Assign(m);
  sVecs_.Assign(m, xVec, false);
  yVecs_.Assign(m, xVec, false);

  uVec.Assign(xVec, false);
  rVec.Assign(xVec, false);
  pVec.Assign(xVec, false);
  qVec.Assign(xVec, false);

  // ----------------------
  // Initialize the CG solver:
  // 𝒖 ← 𝒙,
  // 𝒓 ← 𝒃 - 𝓐𝒖,
  // 𝒑 ← 𝒓,
  // 𝛾 ← <𝒓⋅𝒓>.
  // ----------------------
  Blas::Fill(uVec, 0.0);
  Blas::RandFill(qVec);
  Blas::Scale(qVec, qVec, 1.0/Blas::Norm2(qVec));
  Blas::Add(qVec, bVec, qVec, 0.001);
  linOp.Residual(rVec, qVec, uVec);
  //Blas::Set(uVec, xVec);
  //linOp.Residual(rVec, bVec, uVec);
  Blas::Set(pVec, rVec);
  alpha = Blas::Dot(rVec, rVec);

  size_t const N = 400;
  for (size_t k = 0; k < N; ++k) {

    // ----------------------
    // Iterate:
    // 𝒒 ← 𝓐𝒑,
    // 𝛼̅ ← 𝛼,
    // 𝛼 ← 𝛼̅/<𝒑⋅𝒒>,
    // 𝒖 ← 𝒖 + 𝛼⋅𝒑,
    // 𝒓 ← 𝒓 - 𝛼⋅𝒒,
    // ----------------------
    linOp.MatVec(qVec, pVec);
    real_t const alphaBar = alpha;
    alpha /= Blas::Dot(pVec, qVec);
    Blas::Add(uVec, uVec, pVec, alpha);
    Blas::Sub(rVec, rVec, qVec, alpha);

    // ----------------------
    // 𝒔ₖ ← 𝛼⋅𝒑,
    // 𝒚ₖ ← 𝛼⋅𝒒,
    // 𝜌ₖ ← 𝟣/<𝒔ₖ⋅𝒚ₖ>.
    // ----------------------
    if (k % 20 == 0) {
      size_t const l = k/20;
      Blas::Scale(sVecs_(l), pVec, alpha);
      Blas::Scale(yVecs_(l), qVec, alpha);
      rho_(l) = 1.0/Blas::Dot(sVecs_(l), yVecs_(l));
    }

    // ----------------------
    // 𝛼 ← <𝒓⋅𝒓>,
    // 𝛽 ← 𝛼/𝛼̅,
    // 𝒑 ← 𝒓 + 𝛽⋅𝒑.
    // ----------------------
    alpha = Blas::Dot(rVec, rVec);
    real_t const beta = alpha/alphaBar;
    Blas::Add(pVec, rVec, pVec, beta);

    std::cout << std::sqrt(alpha) << std::endl;

  }

  omega_ = 1.0/(rho_(m - 1)*Blas::Dot(yVecs_(m - 1), yVecs_(m - 1)));
#endif

} // BfgsPreconditioner<...>::Build

template<VectorLike Vector>
void BfgsPreconditioner<Vector>::MatVec(Vector& uVec,
                                        Vector const& xVec) const {

  size_t const m = Rank;

  // ----------------------
  // Compute the approximate inverse Hessian product:
  // 𝒖 ← (𝓗ₘ)𝒙.
  // ----------------------
  MatVec(m, uVec, xVec);

} // BfgsPreconditioner<...>::MatVec

template<VectorLike Vector>
void BfgsPreconditioner<Vector>::MatVec(size_t k,
                                        Vector& uVec,
                                        Vector const& xVec) const {

#if 0
  // ----------------------
  // Compute 𝒖 ← (𝓗ₖ)𝒙, where:
  // 𝓗₀ = 𝜔⋅𝓘,
  // 𝓗ₖ₊₁ = (𝓥ₖ)*(𝓗ₖ)(𝓥ₖ) + 𝜌ₖ⋅𝒔ₖ⋅(𝒔ₖ)*,
  // 𝓥ₖ = 𝓘 - 𝜌ₖ⋅𝒚ₖ⋅(𝒔ₖ)*.
  // ----------------------

  if (k == 0) {

    // ----------------------
    // Compute 𝒖 ← (𝓗₀)𝒙:
    // 𝒖 ← 𝜔⋅𝒙.
    // ----------------------
    Blas::Scale(uVec, xVec, omega_);

  } else {

    // ----------------------
    // Compute 𝒖 ← (𝓗ₖ)𝒙 = 
    //   (𝓥ₖ₋₁)*(𝓗ₖ₋₁)(𝓥ₖ₋₁) + 𝜌ₖ₋₁⋅𝒔ₖ₋₁⋅(𝒔ₖ₋₁)*
    //   with 𝓥ₖ₋₁ = 𝓘 - 𝜌ₖ₋₁⋅𝒚ₖ₋₁⋅(𝒔ₖ₋₁)*:
    // 𝛼 ← <𝒔ₖ₋₁⋅𝒙>,
    // 𝒖 ← 𝒙 - 𝛼⋅𝜌ₖ₋₁⋅𝒚ₖ₋₁,
    // 𝒖 ← (𝓗ₖ₋₁)𝒖,
    // 𝛽 ← 𝛼 - <𝒚ₖ₋₁⋅𝒖>,
    // 𝒖 ← 𝒖 + 𝛽⋅𝜌ₖ₋₁⋅𝒔ₖ₋₁.
    // ----------------------
    real_t const alpha = Blas::Dot(sVecs_(k - 1), xVec);
    Blas::Sub(uVec, xVec, yVecs_(k - 1), alpha*rho_(k - 1));
    MatVec(k - 1, uVec, uVec);
    real_t const beta = alpha - Blas::Dot(yVecs_(k - 1), uVec);
    Blas::Add(uVec, uVec, sVecs_(k - 1), beta*rho_(k - 1));
  }
#endif

} // BfgsPreconditioner<...>::MatVec

} // namespace Storm
