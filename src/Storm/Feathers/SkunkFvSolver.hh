
#pragma once

#include "ConvectionScheme.hh"
#include "Field.hh"
#include "FluxScheme.hh"
#include "GradientLimiterScheme.hh"
#include "SkunkFvBC.hh"

#include <map>

namespace Storm::Feathers {

/**
 * @brief A Finite volume solver.
 */
template<typename MhdPhysicsT>
class MhdFvSolverT :
    public std::enable_shared_from_this<MhdFvSolverT<MhdPhysicsT>> {
public:

  using MhdFluidStateT = typename MhdPhysicsT::MhdFluidStateT;
  static constexpr int_t num_vars = MhdPhysicsT::num_vars;

private:

  std::shared_ptr<const Mesh> m_mesh;
  std::shared_ptr<iConvectionScheme> m_conv;
  std::map<int_t, std::shared_ptr<MhdFvBcPT<MhdPhysicsT>>> m_bcs;

public:

  explicit MhdFvSolverT(std::shared_ptr<const Mesh> mesh);

public:

  void calc_func(tScalarField& u, tScalarField& u_out) const;
  void calc_step(real_t& dt, tScalarField& u, tScalarField& u_hat) const;
}; // class MhdFvSolverT

} // namespace Storm::Feathers