/*
 *  ______  ______   ______   ______  __  __   ______   ______   ______
 * /\  ___\/\  ___\ /\  __ \ /\__  _\/\ \_\ \ /\  ___\ /\  __ \ /\  ___\
 * \ \  __\\ \  _\  \ \  __ \\/_/\ \/\ \  __ \\ \  __\ \ \  __/ \ \___  \
 *  \ \_\   \ \_____\\ \_\ \_\  \ \_\ \ \_\ \_\\ \_____\\ \_\ \_\\/\_____\
 *   \/_/    \/_____/ \/_/\/_/   \/_/  \/_/\/_/ \/_____/ \/_/ /_/ \/_____/
 *
 * Copyright (c) 2021 Oleg Butakov
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#pragma once
#ifndef GRADIENT_SCHEME_HH_
#define GRADIENT_SCHEME_HH_

#include "SkunkBase.hh"
#include <stormMesh/Mesh.hxx>

namespace Storm {

/** Abstract cell-centered gradient scheme. */
class iGradientScheme : public tObject<iGradientScheme> {
public:

  /** Compute cell-centered gradients. */
  virtual void get_gradients(size_t num_vars, tVectorField& grad_u,
                             const tScalarField& u) const = 0;
}; // class iGradientScheme

/**
 * Weighted Least-Squares gradient estimation scheme, cell-based:
 * computes cell-centered gradients based on the cell-centered values.
 *
 * This gradient scheme is a second-order scheme for any meshes.
 * Also, this gradient scheme is by far the fastest one.
 */
class cLeastSquaresGradientScheme final : public iGradientScheme {
private:

  std::shared_ptr<const Mesh> m_mesh;
  tMatrixField m_inverse_matrices;

public:

  /** Initialize the gradient scheme. */
  explicit cLeastSquaresGradientScheme(std::shared_ptr<const Mesh> mesh)
      : m_mesh(std::move(mesh)), m_inverse_matrices(1, m_mesh->cells().size()) {
    init_gradients_();
  }

private:

  void init_gradients_();

public:

  /** Compute cell-centered gradients. */
  void get_gradients(size_t num_vars, tVectorField& grad_u,
                     const tScalarField& u) const final;
}; // class cLeastSquaresGradientScheme

} // namespace Storm

#endif // GRADIENT_SCHEME_HH_
