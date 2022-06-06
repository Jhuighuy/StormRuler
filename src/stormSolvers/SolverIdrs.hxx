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

#include <utility>

#include <stormBase.hxx>
#include <stormSolvers/LegacyTensor.hxx>
#include <stormSolvers/Solver.hxx>
#include <stormSolvers/Subspace.hxx>
#include <stormSolvers/Vector.hxx>

namespace Storm {

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief The @c IDR(s) (Induced Dimension Reduction)
///   linear operator equation solver.
///
/// References:
/// @verbatim
/// [1] Peter Sonneveld, Martin B. van Gijzen.
///     “IDR(s): A Family of Simple and Fast Algorithms for Solving
///      Large Nonsymmetric Systems of Linear Equations.”
///     SIAM J. Sci. Comput. 31 (2008): 1035-1062.
/// [2] Martin B. van Gijzen, Peter Sonneveld.
///     “Algorithm 913: An Elegant IDR(s) Variant that Efficiently
///      Exploits Biorthogonality Properties.”
///     ACM Trans. Math. Softw. 38 (2011): 5:1-5:19.
/// @endverbatim
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<VectorLike Vector>
class IdrsSolver final : public InnerOuterIterativeSolver<Vector> {
private:

  real_t omega_;
  stormVector<real_t> phi_, gamma_;
  stormMatrix<real_t> mu_;
  Vector r_vec_, v_vec_, z_vec_;
  Subspace<Vector> p_vecs_, u_vecs_, g_vecs_;

  real_t outer_init(Vector const& x_vec, Vector const& b_vec,
                    Operator<Vector> const& lin_op,
                    Preconditioner<Vector> const* pre_op) override;

  void inner_init(Vector const& x_vec, Vector const& b_vec,
                  Operator<Vector> const& lin_op,
                  Preconditioner<Vector> const* pre_op) override;

  real_t inner_iterate(Vector& x_vec, Vector const& b_vec,
                       Operator<Vector> const& lin_op,
                       Preconditioner<Vector> const* pre_op) override;

public:

  IdrsSolver() {
    this->num_inner_iterations = 4;
  }

}; // class IdrsSolver

template<VectorLike Vector>
real_t IdrsSolver<Vector>::outer_init(Vector const& x_vec, Vector const& b_vec,
                                      Operator<Vector> const& lin_op,
                                      Preconditioner<Vector> const* pre_op) {
  size_t const s{this->num_inner_iterations};

  bool const left_pre{(pre_op != nullptr) &&
                      (this->pre_side == PreconditionerSide::Left)};

  phi_.assign(s);
  gamma_.assign(s);
  mu_.assign(s, s);

  r_vec_.assign(x_vec, false);
  v_vec_.assign(x_vec, false);
  if (pre_op != nullptr) { z_vec_.assign(x_vec, false); }

  p_vecs_.assign(s, x_vec, false);
  u_vecs_.assign(s, x_vec, false);
  g_vecs_.assign(s, x_vec, false);

  // Initialize:
  // ----------------------
  // 𝒓 ← 𝒃 - 𝓐𝒙,
  // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
  //   𝒛 ← 𝒓,
  //   𝒓 ← 𝓟𝒛.
  // 𝗲𝗻𝗱 𝗶𝗳
  // 𝜑₀ ← ‖𝒓‖.
  // ----------------------
  lin_op.Residual(r_vec_, b_vec, x_vec);
  if (left_pre) {
    std::swap(z_vec_, r_vec_);
    pre_op->mul(r_vec_, z_vec_);
  }
  phi_(0) = Blas::Norm2(r_vec_);

  return phi_(0);

} // IdrsSolver::outer_init

template<VectorLike Vector>
void IdrsSolver<Vector>::inner_init(Vector const& x_vec, Vector const& b_vec,
                                    Operator<Vector> const& lin_op,
                                    Preconditioner<Vector> const* pre_op) {
  size_t const s{this->num_inner_iterations};

  // Build shadow space and initialize 𝜑:
  // ----------------------
  // 𝗶𝗳 𝘍𝘪𝘳𝘴𝘵𝘐𝘵𝘦𝘳𝘢𝘵𝘪𝘰𝘯:
  //   𝜔 ← 𝜇₀₀ ← 𝟣,
  //   𝒑₀ ← 𝒓/𝜑₀,
  //   𝗳𝗼𝗿 𝑖 = 𝟣, 𝑠 - 𝟣 𝗱𝗼:
  //     𝜇ᵢᵢ ← 𝟣, 𝜑ᵢ ← 𝟢,
  //     𝒑ᵢ ← 𝘙𝘢𝘯𝘥𝘰𝘮,
  //     𝗳𝗼𝗿 𝑗 = 𝟢, 𝑖 - 𝟣 𝗱𝗼:
  //       𝜇ᵢⱼ ← 𝟢,
  //       𝒑ᵢ ← 𝒑ᵢ - <𝒑ᵢ⋅𝒑ⱼ>⋅𝒑ⱼ,
  //     𝗲𝗻𝗱 𝗳𝗼𝗿
  //     𝒑ᵢ ← 𝒑ᵢ/‖𝒑ᵢ‖.
  //   𝗲𝗻𝗱 𝗳𝗼𝗿
  // 𝗲𝗹𝘀𝗲:
  //   𝗳𝗼𝗿 𝑖 = 𝟢, 𝑠 - 𝟣 𝗱𝗼:
  //     𝜑ᵢ ← <𝒑ᵢ⋅𝒓>.
  //   𝗲𝗻𝗱 𝗳𝗼𝗿
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  bool const first_iteration{this->iteration == 0};
  if (first_iteration) {
    omega_ = mu_(0, 0) = 1.0;
    p_vecs_(0) <<= r_vec_ / phi_(0);
    for (size_t i{1}; i < s; ++i) {
      mu_(i, i) = 1.0, phi_(i) = 0.0;
      Blas::RandFill(p_vecs_(i));
      for (size_t j{0}; j < i; ++j) {
        mu_(i, j) = 0.0;
        p_vecs_(i) -= Blas::Dot(p_vecs_(i), p_vecs_(j)) * p_vecs_(j);
      }
      p_vecs_(i) /= Blas::Norm2(p_vecs_(i));
    }
  } else {
    for (size_t i{0}; i < s; ++i) {
      phi_(i) = Blas::Dot(p_vecs_(i), r_vec_);
    }
  }

} // IdrsSolver::inner_init

template<VectorLike Vector>
real_t IdrsSolver<Vector>::inner_iterate(Vector& x_vec, Vector const& b_vec,
                                         Operator<Vector> const& lin_op,
                                         Preconditioner<Vector> const* pre_op) {
  size_t const s{this->num_inner_iterations};
  size_t const k{this->inner_iteration};

  bool const left_pre{(pre_op != nullptr) &&
                      (this->pre_side == PreconditionerSide::Left)};
  bool const right_pre{(pre_op != nullptr) &&
                       (this->pre_side == PreconditionerSide::Right)};

  // Compute 𝛾:
  // ----------------------
  // 𝛾ₖ:ₛ₋₁ ← (𝜇ₖ:ₛ₋₁,ₖ:ₛ₋₁)⁻¹⋅𝜑ₖ:ₛ₋₁.
  // ----------------------
  for (size_t i = k; i < s; ++i) {
    gamma_(i) = phi_(i);
    for (size_t j = k; j < i; ++j) {
      gamma_(i) -= mu_(i, j) * gamma_(j);
    }
    gamma_(i) /= mu_(i, i);
  }

  // Compute the new 𝒈ₖ and 𝒖ₖ vectors:
  // ----------------------
  // 𝒗 ← 𝒓 - 𝛾ₖ⋅𝒈ₖ,
  // 𝗳𝗼𝗿 𝑖 = 𝑘 + 𝟣, 𝑠 - 𝟣 𝗱𝗼:
  //   𝒗 ← 𝒗 - 𝛾ᵢ⋅𝒈ᵢ,
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // 𝗶𝗳 𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦:
  //   𝒛 ← 𝒗,
  //   𝒗 ← 𝓟𝒛,
  // 𝗲𝗻𝗱 𝗶𝗳
  // 𝒖ₖ ← 𝜔⋅𝒗 + 𝛾ₖ⋅𝒖ₖ,
  // 𝗳𝗼𝗿 𝑖 = 𝑘 + 𝟣, 𝑠 - 𝟣 𝗱𝗼:
  //   𝒖ₖ ← 𝒖ₖ + 𝛾ᵢ⋅𝒖ᵢ,
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
  //   𝒈ₖ ← 𝓟(𝒛 ← 𝓐𝒖ₖ).
  // 𝗲𝗹𝘀𝗲:
  //   𝒈ₖ ← 𝓐𝒖ₖ.
  // 𝗲𝗻𝗱 𝗶𝗳
  // ----------------------
  v_vec_ <<= r_vec_ - gamma_(k) * g_vecs_(k);
  for (size_t i{k + 1}; i < s; ++i) {
    v_vec_ -= gamma_(i) * g_vecs_(i);
  }
  if (right_pre) {
    std::swap(z_vec_, v_vec_);
    pre_op->mul(v_vec_, z_vec_);
  }
  u_vecs_(k) <<= omega_ * v_vec_ + gamma_(k) * u_vecs_(k);
  for (size_t i{k + 1}; i < s; ++i) {
    u_vecs_(k) += gamma_(i) * u_vecs_(i);
  }
  if (left_pre) {
    pre_op->mul(g_vecs_(k), z_vec_, lin_op, u_vecs_(k));
  } else {
    lin_op.mul(g_vecs_(k), u_vecs_(k));
  }

  // Biorthogonalize the new vectors 𝒈ₖ and 𝒖ₖ:
  // ----------------------
  // 𝗳𝗼𝗿 𝑖 = 𝟢, 𝑘 - 𝟣 𝗱𝗼:
  //   𝛼 ← <𝒑ᵢ⋅𝒈ₖ>/𝜇ᵢᵢ,
  //   𝒖ₖ ← 𝒖ₖ - 𝛼⋅𝒖ᵢ,
  //   𝒈ₖ ← 𝒈ₖ - 𝛼⋅𝒈ᵢ.
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // ----------------------
  for (size_t i{0}; i < k; ++i) {
    real_t const alpha{
        utils::safe_div(Blas::Dot(p_vecs_(i), g_vecs_(k)), mu_(i, i))};
    u_vecs_(k) -= alpha * u_vecs_(i);
    g_vecs_(k) -= alpha * g_vecs_(i);
  }

  // Compute the new column of 𝜇:
  // ----------------------
  // 𝗳𝗼𝗿 𝑖 = 𝑘, 𝑠 - 𝟣 𝗱𝗼:
  //   𝜇ᵢₖ ← <𝒑ᵢ⋅𝒈ₖ>.
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // ----------------------
  for (size_t i{k}; i < s; ++i) {
    mu_(i, k) = Blas::Dot(p_vecs_(i), g_vecs_(k));
  }

  // Update the solution and the residual:
  // ----------------------
  // 𝛽 ← 𝜑ₖ/𝜇ₖₖ,
  // 𝒙 ← 𝒙 + 𝛽⋅𝒖ₖ,
  // 𝒓 ← 𝒓 - 𝛽⋅𝒈ₖ.
  // ----------------------
  real_t const beta{utils::safe_div(phi_(k), mu_(k, k))};
  x_vec += beta * u_vecs_(k);
  r_vec_ -= beta * g_vecs_(k);

  // Update 𝜑:
  // ----------------------
  // 𝜑ₖ₊₁:ₛ₋₁ ← 𝜑ₖ₊₁:ₛ₋₁ - 𝛽⋅𝜇ₖ₊₁:ₛ₋₁,ₖ.
  // ----------------------
  for (size_t i{k + 1}; i < s; ++i) {
    phi_(i) -= beta * mu_(i, k);
  }

  if (k == s - 1) {
    // Enter the next 𝓖 subspace:
    // ----------------------
    // 𝗶𝗳 𝘓𝘦𝘧𝘵𝘗𝘳𝘦:
    //   𝒗 ← 𝓟(𝒛 ← 𝓐𝒓),
    // 𝗲𝗹𝘀𝗲 𝗶𝗳 𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦:
    //   𝒗 ← 𝓐(𝒛 ← 𝓟𝒓),
    // 𝗲𝗹𝘀𝗲:
    //   𝒗 ← 𝓐𝒓,
    // 𝗲𝗻𝗱 𝗶𝗳
    // 𝜔 ← <𝒗⋅𝒓>/<𝒗⋅𝒗>,
    // 𝒙 ← 𝒙 + 𝜔⋅(𝘙𝘪𝘨𝘩𝘵𝘗𝘳𝘦 ? 𝒛 : 𝒓),
    // 𝒓 ← 𝒓 - 𝜔⋅𝒗.
    // ----------------------
    if (left_pre) {
      pre_op->mul(v_vec_, z_vec_, lin_op, r_vec_);
    } else if (right_pre) {
      lin_op.mul(v_vec_, z_vec_, *pre_op, r_vec_);
    } else {
      lin_op.mul(v_vec_, r_vec_);
    }
    omega_ =
        utils::safe_div(Blas::Dot(v_vec_, r_vec_), Blas::Dot(v_vec_, v_vec_));
    x_vec += omega_ * (right_pre ? z_vec_ : r_vec_);
    r_vec_ -= omega_ * v_vec_;
  }

  return Blas::Norm2(r_vec_);

} // IdrsSolver::inner_iterate

} // namespace Storm
