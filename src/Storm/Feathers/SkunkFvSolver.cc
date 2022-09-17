#include "SkunkFvSolver.hh"

namespace Storm::Feathers {

template<typename MhdPhysicsT>
MhdFvSolverT<MhdPhysicsT>::MhdFvSolverT(std::shared_ptr<const Mesh> mesh)
    : m_mesh(mesh), m_conv(new cUpwind2ConvectionScheme(mesh)) {
  m_bcs[1] = std::make_shared<MhdFvBcFarFieldT<MhdPhysicsT>>();
  m_bcs[2] = std::make_shared<MhdFvBcSlipT<MhdPhysicsT>>();
}

/**
 * @brief Compute spacial discretization.
 */
template<typename MhdPhysicsT>
void MhdFvSolverT<MhdPhysicsT>::calc_func(tScalarField& u,
                                          tScalarField& u_out) const {
  using namespace Storm;

  /*
   * Clear fields and apply boundary conditions.
   */

  ForEach(m_mesh->cells(), [&](CellView<Mesh> cell) { u_out[cell].fill(0.0); });
  for (size_t mark = 1; mark < m_mesh->num_face_labels(); ++mark) {
    const auto& bc = m_bcs.at(mark);
    ForEach(m_mesh->faces(Label(mark)), [&](FaceView<Mesh> face) {
      bc->get_ghost_state(
          face.normal(), face.inner_cell().center(), face.outer_cell().center(),
          u[face.inner_cell()].data(), u[face.outer_cell()].data());
    });
  }

  m_conv->get_cell_convection(5, u_out, u);

} // MhdFvSolverT::calc_func

template<typename MhdPhysicsT>
void MhdFvSolverT<MhdPhysicsT>::calc_step(real_t& dt, tScalarField& u,
                                          tScalarField& u_hat) const {
  /*
   * Compute.
   */
  calc_func(u, u_hat);
  ForEach(m_mesh->unlabeled_cells(), [&](CellView<Mesh> cell) {
    for (uint_t i = 0; i < num_vars; ++i) {
      u_hat[cell][i] = u[cell][i] - dt * u_hat[cell][i];
    }
  });
} // MhdFvSolverT::calc_step

template class MhdFvSolverT<tGasPhysics>;

} // namespace Storm::Feathers
