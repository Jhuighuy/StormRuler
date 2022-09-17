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
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
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
#ifndef CONVECTION_SCHEME_HH_
#define CONVECTION_SCHEME_HH_

#include "SkunkBase.hh"
#include "FluxScheme.hh"
#include "GradientLimiterScheme.hh"
#include "GradientScheme.hh"

namespace Storm {

/**
 * Abstract convection scheme.
 */
class iConvectionScheme : public tObject<iConvectionScheme> {
public:
    /** Compute the nonlinear convection. */
    virtual void get_cell_convection(size_t num_vars,
                                     tScalarField& conv_u,
                                     const tScalarField& u) const = 0;
}; // class iConvectionScheme

/**
 * Piecewise-constant upwind convection scheme.
 * This is a first-order scheme.
 */
class cUpwindConvectionScheme final : public iConvectionScheme {
public:
    std::shared_ptr<const Mesh> m_mesh;
    std::shared_ptr<iFluxScheme> m_flux;

public:
    explicit cUpwindConvectionScheme(std::shared_ptr<const Mesh> mesh):
        m_mesh(std::move(mesh)),
        m_flux(new tLaxFriedrichsFluxScheme<tGasPhysics>()) {
    }

    /** Compute the first-order upwind nonlinear convection. */
    void get_cell_convection(size_t num_vars,
                             tScalarField& div_f,
                             const tScalarField& u) const final;
}; // class cUpwindConvectionScheme

/**
 * Piecewise-linear upwind convection scheme.
 * This is a second-order scheme.
 */
class cUpwind2ConvectionScheme final : public iConvectionScheme {
public:
    std::shared_ptr<const Mesh> m_mesh;
    std::shared_ptr<iFluxScheme> m_flux;
    std::shared_ptr<iGradientScheme> m_gradient_scheme;
    std::shared_ptr<iGradientLimiterScheme> m_gradient_limiter_scheme;

public:
    explicit cUpwind2ConvectionScheme(std::shared_ptr<const Mesh> mesh):
        m_mesh(std::move(mesh)),
        m_flux(new tHllFluxScheme<tGasPhysics>()),
        m_gradient_scheme(new cLeastSquaresGradientScheme(m_mesh)),
        m_gradient_limiter_scheme(new cCubic2GradientLimiterScheme(m_mesh)) {
    }

    /** Compute the second-order upwind nonlinear convection. */
    void get_cell_convection(size_t num_vars,
                             tScalarField& div_f,
                             const tScalarField& u) const final;
}; // class cUpwindConvectionScheme

} // namespace feathers

#endif // CONVECTION_SCHEME_HH_
