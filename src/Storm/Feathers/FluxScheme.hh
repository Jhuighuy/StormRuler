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
#ifndef FLUX_SCHEME_HH_
#define FLUX_SCHEME_HH_

#include "SkunkBase.hh"
#include "SkunkHydro.hh"
#include <stormMesh/Field.hh>
//#include "SkunkFluidPhysics.hh"

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

namespace Storm {

/**
 * Abstract numerical flux.
 */
class iFluxScheme : public tObject<iFluxScheme> {
public:
    /** Compute the numerical flux. */
    virtual void get_numerical_flux(size_t num_vars,
                                    const vec3_t& n,
                                    tScalarConstSubField cons_r,
                                    tScalarConstSubField cons_l,
                                    tScalarSubField flux) const = 0;
}; // class iFluxScheme

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Local Lax-Friedrichs (Rusanov) numerical flux.
 *
 * Use this numerical flux if all other fails. 
 * It should always work.
 */
/** @{ */
template<typename tPhysics>
class tLaxFriedrichsFluxScheme;
template<>
class tLaxFriedrichsFluxScheme<tGasPhysics> final : public iFluxScheme {
public:
    /** Compute the numerical flux. */
    void get_numerical_flux(size_t num_vars,
                            const vec3_t& n,
                            tScalarConstSubField cons_r,
                            tScalarConstSubField cons_l,
                            tScalarSubField flux) const final;
}; // class tLaxFriedrichsFluxScheme<tGasPhysics>
/** @} */

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Harten-Lax-van Leer-Einfeldt numerical flux.
 *
 * Use this numerical flux if HLLC fails. 
 * It should (almost) always work.
 */
/** @{ */
template<typename tPhysics>
class tHllFluxScheme;
template<>
class tHllFluxScheme<tGasPhysics> : public iFluxScheme {
public:
    /** Compute the numerical flux. */
    void get_numerical_flux(size_t num_vars,
                            const vec3_t& n,
                            tScalarConstSubField cons_r,
                            tScalarConstSubField cons_l,
                            tScalarSubField flux) const final;
}; // class tHllFluxScheme<tGasPhysics>
/** @} */

/**
 * Harten-Lax-van Leer-Contact numerical flux.
 *
 * Optimal choice for both gas and plasma physics.
 * In plasma physics case may be a bit more dissipative, but more consistent than HLLD/Roe.
 */
/** @{ */
template<typename tPhysics>
class tHllcFluxScheme;
template<>
class tHllcFluxScheme<tGasPhysics> : public iFluxScheme {
public:
    /** Compute the numerical flux. */
    void get_numerical_flux(size_t num_vars,
                            const vec3_t& n,
                            tScalarConstSubField cons_r,
                            tScalarConstSubField cons_l,
                            tScalarSubField flux) const override;
}; // class tHllcFluxScheme<tGasPhysics>
/** @} */

} // namespace feathers

#endif  // ifndef FLUX_SCHEME_HH_
