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
  //   𝗳𝗼𝗿 𝑘 = 𝟣, 𝑠 - 𝟣 𝗱𝗼:
  //     𝒑ₖ ← random, 
  //   𝗲𝗻𝗱 𝗳𝗼𝗿
  //   𝗳𝗼𝗿 𝑘 = 𝟢, 𝑠 - 𝟣 𝗱𝗼:
  //     𝗳𝗼𝗿 𝑗 = 𝟢, 𝑘 - 𝟣 𝗱𝗼:
  //       𝜇ₖⱼ ← 𝟢,
  //       𝛼 ← <𝒑ₖ⋅𝒑ⱼ>,
  //       𝒑ₖ ← 𝒑ₖ - 𝛼⋅𝒑ⱼ,
  //     𝗲𝗻𝗱 𝗳𝗼𝗿
  //     𝜇ₖₖ ← 𝟣,
  //     𝛼 ← ‖𝒑ₖ‖,
  //     𝒑ₖ ← 𝒑ₖ/𝛼,
  //   𝗲𝗻𝗱 𝗳𝗼𝗿
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  bool const firstIteration = this->Iteration == 0;
  if (firstIteration) {
    omega_ = 1.0;
    stormBlas::Set(pVecs_[0], rVec_);
    for (stormSize_t k = 1; k < s; ++k) {
      stormBlas::RandFill(pVecs_[k]);
    }
    for (stormSize_t k = 0; k < s; ++k) {
      stormReal_t alpha;
      for (stormSize_t j = 0; j < k; ++j) {
        mu_(k, j) = 0.0;
        alpha = stormBlas::Dot(pVecs_[k], pVecs_[j]);
        stormBlas::Sub(pVecs_[k], pVecs_[k], pVecs_[j], alpha);
      }
      mu_(k, k) = 1.0;
      alpha = stormBlas::Norm2(pVecs_[k]);
      stormBlas::Scale(pVecs_[k], pVecs_[k], 1.0/alpha);
    }
  }

  // ----------------------
  // 𝗳𝗼𝗿 𝑘 = 𝟢, 𝑠 - 𝟣 𝗱𝗼:
  //   𝜑ₖ ← <𝒑ₖ⋅𝒓>.
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
  for (stormSize_t j = 0; j < s; ++j) {
    gamma_(j) = phi_(j);
  }
  for (stormSize_t j = 0; j < s; ++j) {
    for (stormSize_t i = 0; i < j; ++i) {
      gamma_(j) -= gamma_(i)*mu_(j, i);
    }
    gamma_(j) /= mu_(j, j);
  }

  // ----------------------
  // 𝒗 ← 𝒓,
  // 𝗳𝗼𝗿 𝑗 = 𝑘, 𝑠 - 𝟣 𝗱𝗼:
  //   𝒗 ← 𝒗 - 𝛾ⱼ⋅𝒈ⱼ,
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // ----------------------
  stormBlas::Set(vVec_, rVec_);
  for (stormSize_t j = k; j < s; ++j) {
    stormBlas::Sub(vVec_, vVec_, gVecs_[j], gamma_(j));
  }

  // ----------------------
  // 𝒗 ← 𝓟𝒗,
  // 𝒖ₖ ← 𝛾ₖ⋅𝒖ₖ,
  // 𝗳𝗼𝗿 𝑗 = 𝑘 + 𝟣, 𝑠 - 𝟣 𝗱𝗼:
  //   𝒖ₖ ← 𝒖ₖ + 𝛾ⱼ⋅𝒖ⱼ,
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // 𝒖ₖ ← 𝒖ₖ + 𝜔⋅𝒗,
  // ----------------------
  /// @todo Apply preconditioning!
  stormBlas::Scale(uVecs_[k], uVecs_[k], gamma_(k));
  for (stormSize_t j = k + 1; j < s; ++j) {
    stormBlas::Add(uVecs_[k], uVecs_[k], uVecs_[j], gamma_(j));
  }
  stormBlas::Add(uVecs_[k], uVecs_[k], vVec_, omega_);

  // ----------------------
  // 𝒈ₖ ← 𝓐𝒖ₖ,
  // 𝗳𝗼𝗿 𝑗 = 𝟢, 𝑘 - 𝟣 𝗱𝗼:
  //   𝛼 ← <𝒈ₖ⋅𝒑ⱼ>/𝜇ⱼⱼ,
  //   𝒈ₖ ← 𝒈ₖ - 𝛼⋅𝒈ⱼ,
  //   𝒖ₖ ← 𝒖ₖ - 𝛼⋅𝒖ⱼ,
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // 𝗳𝗼𝗿 𝑗 = 𝑘, 𝑠 - 𝟣 𝗱𝗼:
  //   𝜇ⱼₖ ← <𝒈ₖ⋅𝒑ⱼ>,
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // ----------------------
  linOp.MatVec(gVecs_[k], uVecs_[k]);
  for (stormSize_t j = 0; j < k; ++j) {
    stormReal_t const alpha = 
      stormBlas::Dot(gVecs_[k], pVecs_[j])/mu_(j, j);
    stormBlas::Sub(gVecs_[k], gVecs_[k], gVecs_[j], alpha);
    stormBlas::Sub(uVecs_[k], uVecs_[k], uVecs_[j], alpha);
  }
  for (stormSize_t j = k; j < s; ++j) {
    mu_(j, k) = stormBlas::Dot(gVecs_[j], pVecs_[k]);
  }

  // ----------------------
  // 𝛽 ← 𝜑ₖ/𝜇ₖₖ,
  // 𝒓 ← 𝒓 - 𝛽⋅𝒈ₖ,
  // 𝒙 ← 𝒙 + 𝛽⋅𝒖ₖ,
  // ----------------------
  stormReal_t const beta = phi_(k)/mu_(k, k);
  stormBlas::Sub(rVec_, rVec_, gVecs_[k], beta);
  stormBlas::Add(xVec, xVec, uVecs_[k], beta);

  // ----------------------
  // 𝗶𝗳 𝑘 < 𝑠 - 𝟣:
  //   𝗳𝗼𝗿 𝑗 = 𝟢, 𝑘 𝗱𝗼:
  //     𝜑ⱼ ← 𝟢,
  //   𝗲𝗻𝗱 𝗳𝗼𝗿
  //   𝗳𝗼𝗿 𝑗 = 𝑘 + 𝟣, 𝑠 - 𝟣 𝗱𝗼:
  //     𝜑ⱼ ← 𝜑ⱼ - 𝛽⋅𝜇ⱼₖ,
  //   𝗲𝗻𝗱 𝗳𝗼𝗿
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  if (k < s - 1) {
    for (stormSize_t j = 0; j <= k; ++j) {
      phi_(j) = 0.0;
    }
    for (stormSize_t j = k + 1; j < s; ++j) {
      phi_(j) -= beta*mu_(j, k);
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
