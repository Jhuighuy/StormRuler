// Orchid/Skunk -- 2D / 3D Euler / MagnetoHydroDynamics solver.
// Copyright(C) Butakov Oleg and Co. 2019.

#pragma once

#include "SkunkBase.hh"
#include <cstring>

// ************************************************************************************
// //
// ************************************************************************************
// //
// ************************************************************************************
// //

static const real_t Gamma = 1.4;
static const real_t Gamma1 = Gamma - 1.0;

using namespace Storm;

class MhdHydroVars {
public:

  real_t rho = 0.0; /**< Fluid density, 𝜌. */
  real_t p = 0.0;   /**< Fluid pressure, 𝑝. */
  vec3_t vel = {};  /**< Fluid velocity, 𝒗. */
  real_t vel_n = 0.0; /**< Fluid velocity normal component, 𝒗ₙ = 𝒗⋅𝒏. */
  real_t eps = 0.0; /**< Fluid internal energy, 𝜀. */
  real_t nrg = 0.0; /**< Fluid specific total energy, 𝐸 = ½𝒗² + 𝜀. */
  real_t ent = 0.0; /**< Fluid specific enthalpy, 𝐻 = 𝐸 + 𝑝/𝜌. */
  real_t c_snd = 0.0; /**< Fluid sound speed, 𝑐 = (∂𝑝/∂𝜌)¹ᐟ². */
  const real_t* rest_prim = nullptr; /**< Additional advected scalars, 𝑞ᵢ. */
  const real_t* rest_cons =
      nullptr; /**< Additional advected scalars, 𝑢ᵢ = 𝜌𝑞ᵢ. */

public:

  explicit MhdHydroVars() = default;
  explicit MhdHydroVars(const vec3_t& n, const real_t* q_cons,
                        const real_t* q_prim = nullptr);

public:

  /** make primitive variables, 𝑸 = (𝜌,𝑝,𝒗,𝑞ᵢ,…)ᵀ. */
  real_t* make_prim(uint_t num_vars, real_t* prim) const {
    *reinterpret_cast<std::array<real_t, 5>*>(prim) = {rho, p, vel.x, vel.y,
                                                       vel.z};
    return prim;
  }

  /** make conserved variables, 𝑼 = (𝜌,𝜌𝐸,𝜌𝒗,𝜌𝑞ᵢ,…)ᵀ. */
  real_t* make_cons(uint_t num_vars, real_t* cons) const {
    *reinterpret_cast<std::array<real_t, 5>*>(cons) = {
        rho, rho * nrg, rho * vel.x, rho * vel.y, rho * vel.z};
    for (uint_t i = 5; i < num_vars; ++i) {
      if (rest_cons != nullptr) {
        cons[i] = rest_cons[i - 5];
      } else if (rest_prim != nullptr) {
        cons[i] = rho * rest_prim[i - 5];
      }
    }
    return cons;
  }

  /** make flux variables, 𝑭ₙ = (𝜌𝒗ₙ,𝜌𝐻𝒗ₙ,𝜌𝒗𝒗ₙ,𝜌𝑞ᵢ𝒗ₙ,…)ᵀ. */
  real_t* make_flux(uint_t num_vars, const vec3_t& n, real_t* flux) const {
    *reinterpret_cast<std::array<real_t, 5>*>(flux) = {
        rho * vel_n, rho * vel_n * ent, rho * vel_n * vel.x + p * n.x,
        rho * vel_n * vel.y + p * n.y, rho * vel_n * vel.z + p * n.z};
    for (uint_t i = 5; i < num_vars; ++i) {
      if (rest_cons != nullptr) {
        flux[i] = vel_n * rest_cons[i - 5];
      } else if (rest_prim != nullptr) {
        flux[i] = rho * vel_n * rest_prim[i - 5];
      }
    }
    return flux;
  }
};

inline MhdHydroVars::MhdHydroVars(const vec3_t& n, const real_t* q_cons,
                                  const real_t* q_prim)
    : MhdHydroVars() {
  if (q_cons != nullptr) {
    rho = q_cons[0];
    nrg = q_cons[1] / rho;
    vel.x = q_cons[2] / rho;
    vel.y = q_cons[3] / rho;
    vel.z = q_cons[4] / rho;
    vel_n = glm::dot(vel, n);
    eps = nrg - 0.5 * glm::dot(vel, vel);
    p = Gamma1 * rho * eps;
    ent = nrg + p / rho;
  } else if (q_prim != nullptr) {
    rho = q_prim[0];
    p = q_prim[1];
    vel.x = q_prim[2];
    vel.y = q_prim[3];
    vel.z = q_prim[4];
    vel_n = glm::dot(vel, n);
    eps = p / rho / Gamma1;
    nrg = eps + 0.5 * glm::dot(vel, vel);
    ent = nrg + p / rho;
  }
  c_snd = std::sqrt(Gamma * p / rho);
}

typedef class MhdHydroVars MhdFluidVarsIdealGas;

class tGasPhysics {
public:

  static constexpr int_t num_vars = 5;
  typedef MhdFluidVarsIdealGas MhdFluidStateT;
  typedef MhdFluidVarsIdealGas tFluidState;
}; // class tGasPhysics
