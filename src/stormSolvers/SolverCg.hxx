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

#include <cmath>
#if 0
#include <fstream>
#include <mkl.h>
#endif

#include <stormBase.hxx>
#include <stormSolvers/Solver.hxx>

namespace Storm {

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief The @c CG (Conjugate Gradients) linear self-adjoint \
///   definite operator equation solver.
///
/// @c CG may be applied to the consistent singular problems,
/// it converges towards..
///
/// References:
/// @verbatim
/// [1] Hestenes, Magnus R. and Eduard Stiefel.
///     “Methods of conjugate gradients for solving linear systems.”
///     Journal of research of the National 
///     Bureau of Standards 49 (1952): 409-435.
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Vector>
class CgSolver final : public IterativeSolver<Vector> {
private:
  real_t gamma_;
  Vector pVec_, rVec_, zVec_;
  std::vector<real_t> Diagonal_, SubDiagonal_;

  real_t Init(Vector const& xVec,
              Vector const& bVec,
              Operator<Vector> const& linOp,
              Preconditioner<Vector> const* preOp) override;

  real_t Iterate(Vector& xVec,
                 Vector const& bVec,
                 Operator<Vector> const& linOp,
                 Preconditioner<Vector> const* preOp) override;

}; // class CgSolver<...>

template<class Vector>
real_t CgSolver<Vector>::Init(Vector const& xVec,
                              Vector const& bVec,
                              Operator<Vector> const& linOp,
                              Preconditioner<Vector> const* preOp) {

  pVec_.Assign(xVec, false);
  rVec_.Assign(xVec, false);
  zVec_.Assign(xVec, false);

  // ----------------------
  // Initialize:
  // 𝒓 ← 𝒃 - 𝓐𝒙.
  // ----------------------
  linOp.Residual(rVec_, bVec, xVec);

  // ----------------------
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝒛 ← 𝓟𝒓,
  //   𝒑 ← 𝒛,
  //   𝛾 ← <𝒓⋅𝒛>,
  // 𝗲𝗹𝘀𝗲:
  //   𝒑 ← 𝒓,
  //   𝛾 ← <𝒓⋅𝒓>.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  if (preOp != nullptr) {
    preOp->MatVec(zVec_, rVec_);
    Blas::Set(pVec_, zVec_);
    gamma_ = Blas::Dot(rVec_, zVec_);
  } else {
    Blas::Set(pVec_, rVec_);
    gamma_ = Blas::Dot(rVec_, rVec_);
  }

  return (preOp != nullptr) ? Blas::Norm2(rVec_) : std::sqrt(gamma_);

} // CgSolver<...>::Init

template<class Vector>
real_t CgSolver<Vector>::Iterate(Vector& xVec,
                                 Vector const& bVec,
                                 Operator<Vector> const& linOp,
                                 Preconditioner<Vector> const* preOp) {

  // ----------------------
  // Iterate:
  // 𝒛 ← 𝓐𝒑,
  // 𝛾̅ ← 𝛾,
  // 𝛼 ← 𝛾/<𝒑⋅𝒛>,
  // 𝒙 ← 𝒙 + 𝛼⋅𝒑,
  // 𝒓 ← 𝒓 - 𝛼⋅𝒛,
  // ----------------------
  linOp.MatVec(zVec_, pVec_);
  real_t const gammaBar = gamma_;
  real_t const alpha = 
    Utils::SafeDivide(gamma_, Blas::Dot(pVec_, zVec_));
  Blas::Add(xVec, xVec, pVec_, alpha);
  Blas::Sub(rVec_, rVec_, zVec_, alpha);

  // ----------------------
  // 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  //   𝒛 ← 𝓟𝒓,
  //   𝛾 ← <𝒓⋅𝒛>,
  // 𝗲𝗹𝘀𝗲:
  //   𝛾 ← <𝒓⋅𝒓>.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  if (preOp != nullptr) {
    preOp->MatVec(zVec_, rVec_);
    gamma_ = Blas::Dot(rVec_, zVec_);
  } else {
    gamma_ = Blas::Dot(rVec_, rVec_);
  }

  // ----------------------
  // 𝛽 ← 𝛾/𝛾̅,
  // 𝒑 ← (𝓟 ≠ 𝗻𝗼𝗻𝗲 ? 𝒛 : 𝒓) + 𝛽⋅𝒑.
  // ----------------------
  real_t const beta = Utils::SafeDivide(gamma_, gammaBar);
  Blas::Add(pVec_, (preOp != nullptr ? zVec_ : rVec_), pVec_, beta);

#if 0
  bool const computeEigenvalues = true;
  if (computeEigenvalues) {

    // ----------------------
    // Update the tridiagonal matrix:
    // 𝗶𝗳 𝘍𝘪𝘳𝘴𝘵𝘐𝘵𝘦𝘳𝘢𝘵𝘪𝘰𝘯:
    //  𝑇ₖₖ ← 𝟣/𝛼,
    // 𝗲𝗹𝘀𝗲:
    //  𝑇ₖₖ ← 𝑇ₖₖ + 𝟣/𝛼,
    // 𝗲𝗻𝗱 𝗶𝗳
    // 𝑇ₖ₊₁,ₖ₊₁ ← 𝛽/𝛼,
    // 𝑇ₖ₊₁,ₖ ← 𝛽¹ᐟ²/𝛼.
    // ----------------------
    real_t const alphaInverse = 1.0/alpha;
    bool const firstIteration = this->Iteration == 0;
    if (firstIteration) {
      Diagonal_.push_back(alphaInverse);
    } else {
      Diagonal_.back() += alphaInverse;
    }
    Diagonal_.push_back(beta*alphaInverse);
    SubDiagonal_.push_back(std::sqrt(beta)*alphaInverse);

    lapack_int const n = this->Iteration + 1;
    if (n == 200) {
      std::vector<real_t> eigenvalues(n);
      std::vector<lapack_int> iblock(n), isplit(n);

      lapack_int m, nsplit, info;
      info = LAPACKE_dstebz('A', 'E', n, 
        0.0, 0.0, 
        0, 0,
        0.0,
        Diagonal_.data(),
        SubDiagonal_.data(),
        &m, &nsplit, eigenvalues.data(), iblock.data(), isplit.data());

      std::cout << "emax = " << eigenvalues.front() << std::endl;
      std::cout << "emin = " << eigenvalues.back() << std::endl;
      std::cout << info << std::endl;

      std::ofstream file("eigenvalues.txt");
      for (real_t const& ev : eigenvalues) {
        file << ev << " 0" << std::endl;
      }
      file.close();

      abort();
    }

  }
#endif

  return (preOp != nullptr) ? Blas::Norm2(rVec_) : std::sqrt(gamma_);

} // CgSolver<...>::Iterate

} // namespace Storm
