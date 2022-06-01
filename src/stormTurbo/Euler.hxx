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

#include <stormBase.hxx>
#include <stormSolvers/Mat.hxx>

namespace Storm::Turbo {

static constexpr real_t Gamma = 1.4;

template<class Value>
class GasState {
public:

  /// @brief Density.
  Value rho; 

  /// @brief Pressure.
  Value p;

  /// @brief Velocity vector.
  Vec3D<Value> vel;

  /// @brief Internal energy.
  Value eps;

  /// @brief Specific total energy.
  Value nrg;

  /// @brief Specific entalpy.
  Value ent;

  /// @brief Sound speed.
  Value snd;

  /// @brief Load the gas state from the \
  ///   conservative variables vector @p cons or \
  ///   the primitive variables vector @p prim. 
  GasState(Vec3D<real_t> const& n,
           Vec<Value, 5> const* cons,
           Vec<Value, 5> const* prim = nullptr);

  /// @brief Make the conserved variables vector.
  Vec<Value, 5> ToCons() const noexcept;

  /// @brief Make the primitve variables vector.
  Vec<Value, 5> ToPrim() const noexcept;

  /// @brief Make the fux variables vector.
  Vec<Value, 5> ToFlux(Vec3D<real_t> const& n) const noexcept;

}; // class GasState

template<class Value>
GasState<Value>::GasState(Vec3D<real_t> const& n,
                          Vec<Value, 5> const* consPtr,
                          Vec<Value, 5> const* primPtr) {

  if (consPtr != nullptr) {
    auto const& cons = *consPtr;
    rho = cons(0);
    nrg = cons(1)/rho;
    vel = {cons(2), cons(3), cons(4)};
    vel = vel/rho;
    eps = nrg - 0.5*dot(vel, vel);
    p   = (Gamma - 1.0)*rho*eps;
  } else if (primPtr != nullptr) {
    auto const& prim = *primPtr;
    rho = prim(0);
    p   = prim(1);
    vel = {prim(2), prim(2), prim(3)};
    eps = p/rho/(Gamma - 1.0);
    nrg = eps + 0.5*dot(vel, vel);
  }
  ent = nrg + p/rho;
  snd = std::sqrt(Gamma*p/rho);

} // GasState::GasState

template<class Value>
Vec<Value, 5> GasState<Value>::ToCons() const noexcept {

  return { rho, rho*nrg, rho*vel(0), rho*vel(1), rho*vel(2) };

} // GasState::ToCons

template<class Value>
Vec<Value, 5> GasState<Value>::ToPrim() const noexcept {

  return { rho, p, vel(0), vel(1), vel(2) };

} // GasState::ToPrim

template<class Value>
Vec<Value, 5> GasState<Value>::ToFlux(Vec3D<real_t> const& n) const noexcept {

  real_t const vn = rho*dot(vel, n);
  return { vn, vn*ent, vn*vel(0) + p*n(0), vn*vel(1) + p*n(1), vn*vel(2) + p*n(2) };

} // GasState::ToPrim

template<class Value>
Vec<Value, 5> LaxFriedrichsFlux(Vec3D<real_t> const& n,
                                Vec<Value, 5> const& rCons,
                                Vec<Value, 5> const& lCons) {

  GasState<Value> ur(n, &rCons), ul(n, &lCons);

  Value const ss = std::max(
    std::abs(dot(n, ur.vel)) + ur.snd,
    std::abs(dot(n, ul.vel)) + ul.snd);

  return 0.5*(ur.ToFlux(n) + ul.ToFlux(n)) - ss*(rCons - lCons);

} // LaxFriedrichsFlux<..>

} // namespace Storm::Turbo
