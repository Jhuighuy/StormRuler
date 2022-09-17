// ************************************************************************************
// // Orchid/Skunk -- 2D / 3D Euler / MagnetoHydroDynamics solver. Copyright(C)
// Butakov Oleg and Co. 2019.
// ************************************************************************************
// //

#include "SkunkFvSolver.hh"

// ************************************************************************************
// //
// ************************************************************************************
// //
// ************************************************************************************
// //

template<typename MhdPhysicsT>
MhdFvSolverT<MhdPhysicsT>::MhdFvSolverT(std::shared_ptr<const cMesh> mesh)
    : m_mesh(mesh), m_conv(new Storm::cUpwind2ConvectionScheme(mesh)) {
  m_bcs[1] = std::make_shared<MhdFvBcFarFieldT<MhdPhysicsT>>();
  m_bcs[2] = std::make_shared<MhdFvBcSlipT<MhdPhysicsT>>();
}

// ------------------------------------------------------------------------------------
// //
// ------------------------------------------------------------------------------------
// //

/**
 * @brief Compute spacial discretization.
 */
template<typename MhdPhysicsT>
void MhdFvSolverT<MhdPhysicsT>::calc_func(Storm::tScalarField& u,
                                          Storm::tScalarField& u_out) const {
  using namespace Storm;

  /*
   * Clear fields and apply boundary conditions.
   */

  ForEach(cell_views(*m_mesh), [&](CellView cell) { u_out[cell].fill(0.0); });
  for (size_t mark = 1; mark < m_mesh->num_face_marks(); ++mark) {
    const auto& bc = m_bcs.at(mark);
    ForEach(face_views(*m_mesh, FaceMark(mark)), [&](FaceView face) {
      bc->get_ghost_state(
          face.normal(), face.inner_cell().center(), face.outer_cell().center(),
          u[face.inner_cell()].data(), u[face.outer_cell()].data());
    });
  }

  m_conv->get_cell_convection(5, u_out, u);

} // MhdFvSolverT::calc_func

template<typename MhdPhysicsT>
void MhdFvSolverT<MhdPhysicsT>::calc_step(real_t& dt, Storm::tScalarField& u,
                                          Storm::tScalarField& u_hat) const {
  /*
   * Compute.
   */
  calc_func(u, u_hat);
  ForEach(int_cell_views(*m_mesh), [&](CellView cell) {
    for (uint_t i = 0; i < num_vars; ++i) {
      u_hat[cell][i] = u[cell][i] - dt * u_hat[cell][i];
    }
  });
} // MhdFvSolverT::calc_step

// ************************************************************************************
// //
// ************************************************************************************
// //
// ************************************************************************************
// //

template class MhdFvSolverT<tGasPhysics>;
