//  ______  ______   ______   ______  __  __   ______   ______   ______
// /\  ___\/\  ___\ /\  __ \ /\__  _\/\ \_\ \ /\  ___\ /\  __ \ /\  ___\ 
// \ \  __\\ \  _\  \ \  __ \\/_/\ \/\ \  __ \\ \  __\ \ \  __/ \ \___  \ 
//  \ \_\   \ \_____\\ \_\ \_\  \ \_\ \ \_\ \_\\ \_____\\ \_\ \_\\/\_____\ 
//   \/_/    \/_____/ \/_/\/_/   \/_/  \/_/\/_/ \/_____/ \/_/ /_/ \/_____/
//
// Copyright (C) 2020 - 2023 Oleg Butakov
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR Allocator PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT
// SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
// OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

#pragma once

#include "Field.hpp"
#include "SkunkHydro.hpp"
#include <functional>

namespace Storm::Feathers {

/**
 * @brief Abstract boundary condition.
 */
template<typename MhdPhysicsT>
class MhdFvBcPT : public std::enable_shared_from_this<MhdFvBcPT<MhdPhysicsT>> {
public:

  using MhdFluidStateT = typename MhdPhysicsT::MhdFluidStateT;
  static constexpr size_t num_vars = MhdPhysicsT::num_vars;

public:

  /** @brief Compute the ghost states. */
  void get_ghost_state(const vec2_t& n, const vec2_t& r, //
                       const real_t* u, real_t* u_ghost) const {
    get_ghost_state(vec3_t(n, 0.0), vec3_t(r, 0.0), u, u_ghost);
  }
  void get_ghost_state(const vec3_t& n, const vec3_t& r, //
                       const real_t* u, real_t* u_ghost) const {
    MhdFluidStateT u_state(n, u), u_ghost_state;
    get_ghost_state_(n, r, u_state, u_ghost_state);
    u_ghost_state.make_cons(5, u_ghost);
  }

private:

  /** @brief Compute the ghost state. */
  virtual void get_ghost_state_(const vec3_t& n, const vec3_t& r,
                                const MhdFluidStateT& u,
                                MhdFluidStateT& u_ghost) const = 0;
}; // MhdFvBcT

/**
 * @brief Far field boundary condition.
 * Sets ghost state values to infinity state.
 */
template<typename MhdPhysicsT>
class MhdFvBcFarFieldT : public MhdFvBcPT<MhdPhysicsT> {
public:

  using typename MhdFvBcPT<MhdPhysicsT>::MhdFluidStateT;
  using MhdFvBcPT<MhdPhysicsT>::num_vars;

private:

  /** @brief Compute the ghost state. */
  void get_ghost_state_(const vec3_t& n, const vec3_t& r,
                        const MhdFluidStateT& u,
                        MhdFluidStateT& u_ghost) const override {
    u_ghost = u;
  }

}; // class MhdFvBcFarFieldT

/**
 * @brief No-slip wall boundary condition.
 */
template<typename MhdPhysicsT>
class MhdFvBcNoSlipT : public MhdFvBcPT<MhdPhysicsT> {
public:

  using typename MhdFvBcPT<MhdPhysicsT>::MhdFluidStateT;
  using MhdFvBcPT<MhdPhysicsT>::num_vars;
  std::function<vec3_t(vec3_t)> vfunc;

  MhdFvBcNoSlipT(std::function<vec3_t(vec3_t)> vfunc = nullptr)
      : vfunc(std::move(vfunc)) {}

private:

  /** @brief Compute the ghost state. */
  void get_ghost_state_(const vec3_t& n, const vec3_t& r,
                        const MhdFluidStateT& u,
                        MhdFluidStateT& u_ghost) const override {
    u_ghost = u;
    if (vfunc != nullptr) {
      u_ghost.vel = vfunc(r);
    } else {
      u_ghost.vel = {};
    }
  }

}; // class MhdFvBcNoSlipT

/**
 * @brief Slip wall boundary condition.
 */
template<typename MhdPhysicsT>
class MhdFvBcSlipT : public MhdFvBcPT<MhdPhysicsT> {
public:

  using typename MhdFvBcPT<MhdPhysicsT>::MhdFluidStateT;
  using MhdFvBcPT<MhdPhysicsT>::num_vars;

private:

  /** @brief Compute the ghost state. */
  void get_ghost_state_(const vec3_t& n, const vec3_t& r,
                        const MhdFluidStateT& u,
                        MhdFluidStateT& u_ghost) const override {
    u_ghost = u;
    u_ghost.vel -= u.vel_n * n;
  }

}; // class MhdFluidBcSlipTMhdFluidBcSlipT

} // namespace Storm::Feathers
