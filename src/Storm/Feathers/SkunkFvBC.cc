/**
 *    ______             __     __  _____ _____
 *   / __/ /____ _____  / /__  /  |/  / // / _ \
 *  _\ \/  '_/ // / _ \/  '_/ / /|_/ / _  / // /
 * /___/_/\_\\_,_/_//_/_/\_\ /_/  /_/_//_/____/
 *
 * Copyright (c) 2019 Oleg Butakov
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

#include "SkunkFvBC.hh"

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

/**
 * @brief Compute the ghost state.
 * Scalars are copied, velocity components are set to zero.
 */
template<typename MhdPhysicsT>
void MhdFvBcFarFieldT<MhdPhysicsT>::get_ghost_state_(const Storm::vec3_t& n, const Storm::vec3_t& r, const Storm::vec3_t& r_ghost,
                                                     const MhdFluidStateT& u,
                                                     MhdFluidStateT& u_ghost) const {
    u_ghost = u;
}   // MhdFvBcFarFieldT::get_ghost_state_

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

/**
 * @brief Compute the ghost state. 
 * Scalars are copied, velocity components are set to zero.
 */
template<typename MhdPhysicsT>
void MhdFvBcNoSlipT<MhdPhysicsT>::get_ghost_state_(const Storm::vec3_t& n, const Storm::vec3_t& r, const Storm::vec3_t& r_ghost,
                                                   const MhdFluidStateT& u,
                                                   MhdFluidStateT& u_ghost) const {
    u_ghost = u;
    if (vfunc != nullptr) {
        u_ghost.vel = vfunc(r_ghost);
    } else {
        u_ghost.vel = {};
    }
}   // MhdFvBcNoSlipT::get_ghost_state_

/**
 * @brief Compute the ghost state. 
 * Scalars are copied, velocity components are set to zero.
 */
template<typename MhdPhysicsT>
void MhdFvBcSlipT<MhdPhysicsT>::get_ghost_state_(const Storm::vec3_t& n, const Storm::vec3_t& r, const Storm::vec3_t& r_ghost,
                                                 const MhdFluidStateT& u,
                                                 MhdFluidStateT& u_ghost) const {
    u_ghost = u;
    u_ghost.vel -= u.vel_n * n;
}   // MhdFvBcSlipT::get_ghost_state_

// ************************************************************************************ //
// ************************************************************************************ //
// ************************************************************************************ //

template class MhdFvBcPT<tGasPhysics>;
template class MhdFvBcFarFieldT<tGasPhysics>;
template class MhdFvBcNoSlipT<tGasPhysics>;
template class MhdFvBcSlipT<tGasPhysics>;
