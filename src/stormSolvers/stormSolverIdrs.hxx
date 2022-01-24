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
#ifndef _STORM_SOLVER_IDRs_HXX_
#define _STORM_SOLVER_IDRs_HXX_

#include <vector>

#include <stormSolvers/stormSolver.hxx>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Solve a non-singular operator equation \
///   equation with the @c IDR(s) (Induced Dimension Reduction) method.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Vector>
class stormIdrsSolver final : public stormInnerOuterIterativeSolver<Vector> {
private:
  stormReal_t psi_, omega_;
  stormVector<stormReal_t> phi_, gamma_;
  stormMatrix<stormReal_t> mu_;
  Vector rVec_, vVec_, tVec_;
  std::vector<Vector> pVecs_, gVecs_, uVecs_;

  void OuterInit(Vector& xVec,
                 Vector const& bVec,
                 stormOperator<Vector> const& linOp,
                 stormPreconditioner<Vector> const* preOp) override;

  stormReal_t InnerInit(Vector& xVec,
                        Vector const& bVec,
                        stormOperator<Vector> const& linOp,
                        stormPreconditioner<Vector> const* preOp) override;

  stormReal_t InnerIterate(Vector& xVec,
                           Vector const& bVec,
                           stormOperator<Vector> const& linOp,
                           stormPreconditioner<Vector> const* preOp) override;

}; // class stormIdrsSolver<...>

template<class Vector>
void stormIdrsSolver<Vector>::OuterInit(Vector& xVec,
                                               Vector const& bVec,
                                               stormOperator<Vector> const& linOp,
                                               stormPreconditioner<Vector> const* preOp) {

  stormSize_t const s = (this->NumInnerIterations = 4);

  phi_.Assign(s);
  gamma_.Assign(s);
  mu_.Assign(s, s);

  rVec_.Assign(xVec, false);
  vVec_.Assign(xVec, false);
  tVec_.Assign(xVec, false);
  pVecs_.resize(s);
  for (Vector& pVec : pVecs_) {
    pVec.Assign(xVec, false);
  }
  gVecs_.resize(s);
  for (Vector& gVec : gVecs_) {
    gVec.Assign(xVec, false);
  }
  uVecs_.resize(s);
  for (Vector& uVec : uVecs_) {
    uVec.Assign(xVec, false);
  }

  // ----------------------
  // 𝒓 ← 𝓐𝒙,
  // 𝒓 ← 𝒃 - 𝒓,
  // 𝜓 ← ‖𝒓‖,
  // ----------------------
  linOp.MatVec(rVec_, xVec);
  stormBlas::Sub(rVec_, bVec, rVec_);
  psi_ = stormBlas::Norm2(rVec_);

} // stormIdrsSolver<...>::OuterInit

template<class Vector>
stormReal_t stormIdrsSolver<Vector>::InnerInit(Vector& xVec,
                                               Vector const& bVec,
                                               stormOperator<Vector> const& linOp,
                                               stormPreconditioner<Vector> const* preOp) {

  stormSize_t const s = this->NumInnerIterations;

  // ----------------------
  // Build shadow space:
  // 𝗶𝗳 𝘍𝘪𝘳𝘴𝘵𝘐𝘵𝘦𝘳𝘢𝘵𝘪𝘰𝘯:
  //   𝜔 ← 𝟣,
  //   𝒑₀ ← 𝒓,
  //   𝗳𝗼𝗿 𝑖 = 𝟣, 𝑠 - 𝟣 𝗱𝗼:
  //     𝒑ᵢ ← random, 
  //   𝗲𝗻𝗱 𝗳𝗼𝗿
  //   𝗳𝗼𝗿 𝑖 = 𝟢, 𝑠 - 𝟣 𝗱𝗼:
  //     𝗳𝗼𝗿 𝑘 = 𝟢, 𝑖 - 𝟣 𝗱𝗼:
  //       𝜇ᵢₖ ← 𝟢,
  //       𝛼 ← <𝒑ᵢ⋅𝒑ₖ>,
  //       𝒑ₖ ← 𝒑ₖ - 𝛼⋅𝒑ᵢ,
  //     𝗲𝗻𝗱 𝗳𝗼𝗿
  //     𝜇ᵢᵢ ← 𝟣,
  //     𝛼 ← ‖𝒑ᵢ‖,
  //     𝒑ᵢ ← 𝒑ᵢ/𝛼,
  //   𝗲𝗻𝗱 𝗳𝗼𝗿
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  bool const firstIteration = this->Iteration == 0;
  if (firstIteration) {
    omega_ = 1.0;
    stormBlas::Set(pVecs_[0], rVec_);
    for (stormSize_t i = 1; i < s; ++i) {
      stormBlas::RandFill(pVecs_[i]);
    }
    for (stormSize_t i = 0; i < s; ++i) {
      stormReal_t alpha;
      for (stormSize_t k = 0; k < i; ++k) {
        mu_(i, k) = 0.0;
        alpha = stormBlas::Dot(pVecs_[i], pVecs_[k]);
        stormBlas::Sub(pVecs_[i], pVecs_[i], pVecs_[k], alpha);
      }
      mu_(i, i) = 1.0;
      alpha = stormBlas::Norm2(pVecs_[i]);
      stormBlas::Scale(pVecs_[i], pVecs_[i], 1.0/alpha);
    }
  }

  // ----------------------
  // 𝗳𝗼𝗿 𝑖 = 𝟢, 𝑠 - 𝟣 𝗱𝗼:
  //   𝜑ᵢ ← <𝒑ᵢ⋅𝒓>.
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // ----------------------
  /// @todo Merge with the upper step. 
  for (stormSize_t k = 0; k < s; ++k) {
    phi_(k) = stormBlas::Dot(pVecs_[k], rVec_);
  }

  return psi_;

} // stormIdrsSolver<...>::InnerInit

template<class Vector>
stormReal_t stormIdrsSolver<Vector>::InnerIterate(Vector& xVec,
                                                  Vector const& bVec,
                                                  stormOperator<Vector> const& linOp,
                                                  stormPreconditioner<Vector> const* preOp) {

  stormSize_t const s = this->NumInnerIterations;
  stormSize_t const k = this->InnerIteration;

  // ----------------------
  // 𝛄 ← 𝑀⁻¹𝞿.
  // ----------------------
  for (stormSize_t i = 0; i < s; ++i) {
    gamma_(i) = phi_(i);
  }
  for (stormSize_t i = 0; i < s; ++i) {
    for (stormSize_t j = 0; j < i; ++j) {
      gamma_(i) -= gamma_(j)*mu_(i, j);
    }
    gamma_(i) /= mu_(i, i);
  }

  // ----------------------
  // 𝒗 ← 𝒓,
  // 𝗳𝗼𝗿 𝑖 = 𝑘, 𝑠 - 𝟣 𝗱𝗼:
  //   𝒗 ← 𝒗 - 𝛾ᵢ⋅𝒈ᵢ,
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // ----------------------
  stormBlas::Set(vVec_, rVec_);
  for (stormSize_t i = k; i < s; ++i) {
    stormBlas::Sub(vVec_, vVec_, gVecs_[i], gamma_(i));
  }

  // ----------------------
  // 𝒗 ← 𝓟𝒗,
  // 𝒖ₖ ← 𝛾ₖ⋅𝒖ₖ,
  // 𝗳𝗼𝗿 𝑖 = 𝑘 + 𝟣, 𝑠 - 𝟣 𝗱𝗼:
  //   𝒖ₖ ← 𝒖ₖ + 𝛾ᵢ⋅𝒖ᵢ,
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // 𝒖ₖ ← 𝒖ₖ + 𝜔⋅𝒗,
  // ----------------------
  /// @todo Apply preconditioning!
  stormBlas::Scale(uVecs_[k], uVecs_[k], gamma_(k));
  for (stormSize_t i = k + 1; i < s; ++i) {
    stormBlas::Add(uVecs_[k], uVecs_[k], uVecs_[i], gamma_(i));
  }
  stormBlas::Add(uVecs_[k], uVecs_[k], vVec_, omega_);

  // ----------------------
  // 𝒈ₖ ← 𝓐𝒖ₖ,
  // 𝗳𝗼𝗿 𝑖 = 𝟢, 𝑘 - 𝟣 𝗱𝗼:
  //   𝛼 ← <𝒑ᵢ⋅𝒈ₖ>/𝜇ᵢᵢ,
  //   𝒈ₖ ← 𝒈ₖ - 𝛼⋅𝒈ᵢ,
  //   𝒖ₖ ← 𝒖ₖ - 𝛼⋅𝒖ᵢ,
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // 𝗳𝗼𝗿 𝑖 = 𝑘, 𝑠 - 𝟣 𝗱𝗼:
  //   𝜇ᵢₖ ← <𝒑ᵢ⋅𝒈ₖ>,
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // ----------------------
  linOp.MatVec(gVecs_[k], uVecs_[k]);
  for (stormSize_t i = 0; i < k; ++i) {
    stormReal_t const alpha = 
      stormBlas::Dot(pVecs_[i], gVecs_[k])/mu_(i, i);
    stormBlas::Sub(gVecs_[k], gVecs_[k], gVecs_[i], alpha);
    stormBlas::Sub(uVecs_[k], uVecs_[k], uVecs_[i], alpha);
  }
  for (stormSize_t i = k; i < s; ++i) {
    mu_(i, k) = stormBlas::Dot(pVecs_[i], gVecs_[k]);
  }

  // ----------------------
  // Update the solution and the residual:
  // 𝛽 ← 𝜑ₖ/𝜇ₖₖ,
  // 𝒓 ← 𝒓 - 𝛽⋅𝒈ₖ,
  // 𝒙 ← 𝒙 + 𝛽⋅𝒖ₖ,
  // ----------------------
  stormReal_t const beta = phi_(k)/mu_(k, k);
  stormBlas::Sub(rVec_, rVec_, gVecs_[k], beta);
  stormBlas::Add(xVec, xVec, uVecs_[k], beta);

  // ----------------------
  // Update 𝞿:
  // 𝗶𝗳 𝑘 < 𝑠 - 𝟣:
  //   𝗳𝗼𝗿 𝑖 = 𝟢, 𝑘 𝗱𝗼:
  //     𝜑ᵢ ← 𝟢,
  //   𝗲𝗻𝗱 𝗳𝗼𝗿
  //   𝗳𝗼𝗿 𝑖 = 𝑘 + 𝟣, 𝑠 - 𝟣 𝗱𝗼:
  //     𝜑ᵢ ← 𝜑ᵢ - 𝛽⋅𝜇ᵢₖ,
  //   𝗲𝗻𝗱 𝗳𝗼𝗿
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  if (k < s - 1) {
    for (stormSize_t i = 0; i <= k; ++i) {
      phi_(i) = 0.0;
    }
    for (stormSize_t i = k + 1; i < s; ++i) {
      phi_(i) -= beta*mu_(i, k);
    }
  }

  // ----------------------
  // Enter the next 𝓖 subspace:
  // 𝗶𝗳 𝑘 = 𝑠 - 𝟣:
  //   𝒕 ← 𝓐(𝒗 ← 𝓟𝒓),
  //   𝜔 ← <𝒕⋅𝒓>/<𝒕⋅𝒕>,
  //   𝒓 ← 𝒓 - 𝜔⋅𝒕,
  //   𝒙 ← 𝒙 + 𝜔⋅𝒗.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  if (k == s - 1) {
    stormBlas::Set(vVec_, rVec_);
    linOp.MatVec(tVec_, vVec_);
    omega_ = stormUtils::SafeDivide(
      stormBlas::Dot(tVec_, rVec_), stormBlas::Dot(tVec_, tVec_));
    stormBlas::Sub(rVec_, rVec_, tVec_, omega_);
    stormBlas::Add(xVec, xVec, vVec_, omega_);
  }

  return stormBlas::Norm2(rVec_);

} // stormIdrsSolver<...>::InnerIterate

#endif // ifndef _STORM_SOLVER_IDRs_HXX_
