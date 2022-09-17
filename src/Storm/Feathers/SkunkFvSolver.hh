// ************************************************************************************
// // Orchid/Skunk -- 2D / 3D Euler / MagnetoHydroDynamics solver. Copyright(C)
// Butakov Oleg and Co. 2019.
// ************************************************************************************
// //

#pragma once

#include "ConvectionScheme.hh"
#include "FluxScheme.hh"
#include "GradientLimiterScheme.hh"
#include "SkunkBase.hh"
#include "SkunkFvBC.hh"

#include <map>

// ************************************************************************************
// //
// ************************************************************************************
// //
// ************************************************************************************
// //

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

  std::shared_ptr<const cMesh> m_mesh;
  std::shared_ptr<Storm::iConvectionScheme> m_conv;
  std::map<int_t, std::shared_ptr<MhdFvBcPT<MhdPhysicsT>>> m_bcs;

public:

  explicit MhdFvSolverT(std::shared_ptr<const cMesh> mesh);

public:

  void calc_func(Storm::tScalarField& u, Storm::tScalarField& u_out) const;
  void calc_step(real_t& dt, Storm::tScalarField& u,
                 Storm::tScalarField& u_hat) const;
}; // class MhdFvSolverT

// ************************************************************************************
// //
// ************************************************************************************
// //
// ************************************************************************************
// //
