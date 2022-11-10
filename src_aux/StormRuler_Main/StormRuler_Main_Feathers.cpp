// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// Copyright (C) 2021 Oleg Butakov
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights  to use,
// copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following
// conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

#define _USE_MATH_DEFINES 1

#include <Storm/Utils/Banner.hpp>

#include <Storm/Solvers/Operator.hpp>
#include <Storm/Solvers/SolverCg.hpp>

#include <Storm/Mallard/IoTetgen.hpp>
#include <Storm/Mallard/IoVtk.hpp>
#include <Storm/Mallard/MeshUnstructured.hpp>
#include <Storm/Mallard/Shape.hpp>
// #include <Storm/Bittern/Matrix.hpp>
#include <Storm/Feathers/SkunkFvSolver.hpp>

#include <algorithm>
#include <cstring>
#include <fstream>
#include <math.h>
#include <ranges>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

using namespace Storm;
using namespace Storm::Feathers;

inline std::string my_to_string(size_t i) {
  std::string s = std::to_string(i);
  std::string z(5 - s.size(), '0');
  return z + s;
}

template<std::ranges::input_range Range>
auto ForEachSum(Range&& r, auto x, auto f) {
  std::ranges::for_each(r, [&](auto y) { x += f(y); });
  return x;
}

template<size_t N = 5>
void save_vtk(auto& mesh, const char* path,
              const std::vector<sFieldDesc<N>>& fields) {
  std::ofstream file(path);
  file << std::setprecision(std::numeric_limits<real_t>::digits10 + 1);
  file << "# vtk DataFile Version 2.0" << std::endl;
  file << "# Generated by Feathers/StormRuler/Mesh2VTK" << std::endl;
  file << "ASCII" << std::endl;
  file << "DATASET UNSTRUCTURED_GRID" << std::endl;

  file << "POINTS " << mesh.num_nodes() << " double" << std::endl;
  std::ranges::for_each(mesh.nodes(), [&](auto node) {
    const auto& pos = node.position();
    file << pos.x << " " << pos.y << " " << 0.0 << std::endl;
  });
  file << std::endl;

  size_t const sumNumCellAdjNodes =
      ForEachSum(mesh.interior_cells(), size_t(0),
                 [](auto cell) { return cell.nodes().size() + 1; });
  file << "CELLS " << mesh.num_cells({}) << " " << sumNumCellAdjNodes
       << std::endl;
  std::ranges::for_each(mesh.interior_cells(), [&](auto cell) {
    file << cell.nodes().size() << " ";
    cell.for_each_node(
        [&](NodeIndex node_index) { file << node_index << " "; });
    file << std::endl;
  });
  file << std::endl;

  file << "CELL_TYPES " << mesh.num_cells({}) << std::endl;
  std::ranges::for_each(mesh.interior_cells(),
                        [&](auto cell) { file << "5" << std::endl; });
  file << std::endl;

  file << "CELL_DATA " << mesh.num_cells({}) << std::endl;
  for (const sFieldDesc<N>& field : fields) {
    file << "SCALARS " << field.name << " double 1" << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;
    std::ranges::for_each(mesh.interior_cells(), [&](auto cell) {
      file << (*field.scalar)[cell] /*[field.var_index]*/ << std::endl;
    });
  }
  file << std::endl;
} // Mesh::save_vtk

#if 1
namespace {
static double tau = 1.0e-3, Gamma = 1.0e-4, sigma = 2.0;

template<mesh Mesh>
static void stormDivGrad(const Mesh& mesh, //
                         CellField<Mesh, real_t>& u, real_t dt,
                         const CellField<Mesh, real_t>& c) {
  std::ranges::for_each(mesh.interior_faces(), [&](FaceView<Mesh> face) {
    const CellView<Mesh> cell_inner = face.inner_cell();
    const CellView<Mesh> cell_outer = face.outer_cell();

    // Reconstruct the face values.

    // Compute the flux.
    const auto flux = dt * (c[cell_outer] - c[cell_inner]) /
                      glm::length(cell_outer.center() - cell_inner.center());
    u[cell_inner] += (face.area() / cell_inner.volume()) * flux;
    u[cell_outer] -= (face.area() / cell_outer.volume()) * flux;
  });
}

template<mesh Mesh>
static void cahn_hilliard_step(const Mesh& mesh, //
                               const CellField<Mesh, real_t>& c,
                               // const StormArray<Vec2D<real_t>>& v, //
                               CellField<Mesh, real_t>& c_hat,
                               CellField<Mesh, real_t>& w_hat) {
  // SetBCs_c(mesh, c, c);
  // SetBCs_v(mesh, v);

  constexpr auto dF_dc = [](real_t c) noexcept {
    return 2.0 * c * (c - 1.0) * (2.0 * c - 1.0);
  };

  CellField<Mesh, real_t> f{mesh};

  f <<= map(dF_dc, c);

  c_hat <<= c;
  solve<CgSolver>(c_hat, c,
                  *make_operator<CellField<Mesh, real_t>>(
                      [&](CellField<Mesh, real_t>& c_hat,
                          const CellField<Mesh, real_t>& c_in) {
                        // w_hat <<= f + sigma * (c_in - c) - Gamma *
                        // DIVGRAD(c_in);
                        w_hat <<= f + sigma * (c_in - c);
                        // SetBCs_c(mesh, c_in, c);
                        stormDivGrad(mesh, w_hat, -Gamma, c_in);

                        // c_hat <<= c_in + tau * CONV(v, c_in) - tau *
                        // DIVGRAD(w_hat)
                        c_hat <<= c_in;
                        // SetBCs_w(mesh, w_hat);
                        // stormConvection(mesh, c_hat, -tau, c_in, v);
                        stormDivGrad(mesh, c_hat, -tau, w_hat);
                      }));

  // w_hat = dF_dc(c_hat) - Gamma * DIVGRAD(c_hat)
  // w_hat <<= map(dF_dc, c_hat);
  // SetBCs_c(mesh, c_hat, c);
  // stormDivGrad(mesh, w_hat, -Gamma, c_hat);

} // cahn_hilliard_step

template<mesh Mesh>
static void cahn_hilliard_solve(const Mesh& mesh) {
  CellField<Mesh, real_t> c{mesh};
  CellField<Mesh, real_t> c_hat{mesh};
  CellField<Mesh, real_t> w_hat{mesh};

  std::ranges::for_each(mesh.interior_cells(), [&](CellView<Mesh> cell) {
    c[cell] = (1.0 * rand()) / RAND_MAX;
  });

  double total_time = 0.0;
  for (int time = 0; time <= 200000; ++time) {
    for (int frac = 0; time != 0 && frac < 1; ++frac) {
      struct timespec start, finish;
      clock_gettime(CLOCK_MONOTONIC, &start);

      cahn_hilliard_step(mesh, c, c_hat, w_hat);

      clock_gettime(CLOCK_MONOTONIC, &finish);
      double elapsed = (finish.tv_sec - start.tv_sec);
      elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
      total_time += elapsed;

      std::swap(c, c_hat);
    }

    printf("time = %f\n", total_time);
    save_vtk<1>(mesh, ("out/fields-" + my_to_string(time) + ".vtk").c_str(),
                {{"c", 0, &c}});
  }
}
} // namespace
#endif

#if 0
template<mesh Mesh>
static void euler_solve(const std::shared_ptr<Mesh>& mesh) {
  CellField<Mesh, real_t, 5> uc(mesh->num_cells());
  CellField<Mesh, real_t, 5> up(mesh->num_cells());
  for (size_t cell_ind = 0; cell_ind < mesh->num_cells(); ++cell_ind) {
    std::array<real_t, 5> q{2.0, 1.0, 1.0, 0.0, 0.0};
    // std::array<real_t, 5> q{1.4, 1.0, 3.0, 0.0, 0.0};
    MhdHydroVars v({}, nullptr, q.data());
    v.make_cons(5, uc[CellIndex<Mesh>{cell_ind}].data());
  }
  real_t dt = 1e-4;

  MhdFvSolverT<tGasPhysics> solver(mesh);
  const size_t freq = 200;
  save_vtk(*mesh, ("out/fields-" + my_to_string(0) + ".vtk").c_str(),
           {{"rho", 0, &uc}});
  real_t tt = 0.0;
  {
    for (size_t l = 1; l <= 2000000; ++l) {
      solver.calc_step(dt, uc, up);
      if (l % freq == 0) {
        std::cout << l / freq << "\t" << tt << "\t" << std::endl;
        save_vtk(*mesh,
                 ("out/fields-" + my_to_string(l / freq) + ".vtk").c_str(),
                 {{"rho", 0, &uc}});
      }

      uc.swap(up);
    }
  }
} // euler_solve
#endif

int main(int argc, char** argv) {
  print_banner();

  UnstructuredMesh<2, 2, VovTable> mesh1{};
  read_mesh_from_tetgen(mesh1, "tests/_data/mesh/step.1.");

  auto mesh = std::make_shared<UnstructuredMesh<2, 2, CsrTable>>();
  mesh->assign(std::move(mesh1));

  STORM_INFO_("mesh has {} edges", mesh->num_edges());
  STORM_INFO_("mesh has {} faces", mesh->num_faces());
  STORM_INFO_("mesh has {} cells", mesh->num_cells());
  STORM_INFO_("mesh loaded");

#if 1
  cahn_hilliard_solve(*mesh);
#endif

#if 0
  euler_solve(mesh);
#endif

  return 0;
}
