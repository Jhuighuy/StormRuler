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
#define _GNU_SOURCE 1
#define STORM_HEADER_ONLY 1

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define YURI 0

#include <StormRuler_API.h>

#include <Storm/Utils/Banner.hpp>

#include <Storm/Bittern/FastVector.hpp>
#include <Storm/Bittern/Mat.hpp>

#if 1
#include <Storm/Mallard/IoTetgen.hpp>
#include <Storm/Mallard/IoVtk.hpp>
#include <Storm/Mallard/MeshUnstructured.hpp>
#include <Storm/Mallard/Shape.hpp>
#include <Storm/Mallard/Visualizer.hpp>
// #include <Storm/Bittern/Matrix.hpp>
#endif

#include "NVT/Nvt.hpp"
#include <Storm/Solvers/PreconditionerFactory.hpp>
#include <Storm/Solvers/SolverFactory.hpp>
#include <Storm/Utils/Table.hpp>

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

template<class T>
void stormLinSolve2(const SolverType& method,
                    const PreconditionerType& preMethod, //
                    StormArray<T>& x, const StormArray<T>& b, auto matVec,
                    bool uniform = true) {
  size_t numMatVecs = 0;
  auto symOp = make_symmetric_operator<StormArray<T>>(
      [&](StormArray<T>& y, const StormArray<T>& x) {
        numMatVecs += 1;
        matVec(y, x);
      });

  auto solver = make_iterative_solver<StormArray<T>>(method);
  solver->pre_op = make_preconditioner<StormArray<T>>(preMethod);

  if (uniform) {
    solver->solve(x, b, *symOp);
  } else {
    solve_non_uniform(*solver, x, b, *symOp);
  }

} // stormLinSolve2

template<class T>
void stormNonlinSolve2(const SolverType& method, //
                       StormArray<T>& x, const StormArray<T>& b, auto matVec) {
  auto op = make_symmetric_operator<StormArray<T>>(matVec);
  auto solver = std::make_unique<JfnkSolver<StormArray<T>>>();
  solver->absolute_error_tolerance = 1.0e-4;
  solver->relative_error_tolerance = 1.0e-4;
  solver->solve(x, b, *op);

} // stormNonLinSolve2

static double tau = 1.0e-2, Gamma = 1.0e-4, sigma = 1.0, Sigma = 10.0;

static void SetBCs_c(stormMesh_t mesh, stormArray_t c_hat, stormArray_t c) {
  stormApplyBCs(mesh, c_hat, SR_ALL, SR_PURE_NEUMANN);
  // stormApplyBCs_CosWall(
  //     mesh, c_hat, c,
  //     0.01 / std::sqrt(Gamma) * 6.0 * std::cos(M_PI / 2 - M_PI / 18), 1);
  //  stormApplyBCs(mesh, c_hat, 4, SR_DIRICHLET(1.0));
} // SetBCs_c

static void SetBCs_w(stormMesh_t mesh, stormArray_t w) {
  stormApplyBCs(mesh, w, SR_ALL, SR_PURE_NEUMANN);
} // SetBCs_w

static void SetBCs_p(stormMesh_t mesh, stormArray_t p) {
  stormApplyBCs(mesh, p, SR_ALL, SR_PURE_NEUMANN);
  // stormApplyBCs(mesh, p, 2, SR_DIRICHLET(0.0));
  // stormApplyBCs(mesh, p, 4,
  //               SR_DIRICHLET(0.75 * Sigma * std::cos(M_PI / 2 - M_PI / 18) /
  //                            (2.0 * 0.01 * 26)));
} // SetBCs_p

static void SetBCs_v(stormMesh_t mesh, stormArray_t v) {
  stormApplyBCs(mesh, v, SR_ALL, SR_PURE_DIRICHLET);
  stormApplyBCs_SlipWall(mesh, v, 3);
  // stormApplyBCs(mesh, v, 2, SR_PURE_NEUMANN);
  // stormApplyBCs_InOutLet(mesh, v, 2);
  // stormApplyBCs(mesh, v, 4, SR_PURE_NEUMANN);
  stormApplyBCs_InOutLet(mesh, v, 4);
} // SetBCs_v

static void CahnHilliard_Step(stormMesh_t mesh, //
                              const StormArray<real_t>& c,
                              const StormArray<Vec2D<real_t>>& v, //
                              StormArray<real_t>& c_hat,
                              StormArray<real_t>& w_hat) {
  SetBCs_c(mesh, c, c);
  SetBCs_v(mesh, v);

  constexpr auto dF_dc = [](real_t c) noexcept {
    return 2.0 * c * (c - 1.0) * (2.0 * c - 1.0);
  };

  StormArray<real_t> f;
  f.assign(c, false);

  f <<= map(dF_dc, c);

  c_hat <<= c;
  stormLinSolve2(
      Storm::SolverType::Gmres, Storm::PreconditionerType::None /*"extr"*/,
      c_hat, c,
      [&](StormArray<real_t>& c_hat, const StormArray<real_t>& c_in) {
        // w_hat <<= f + 2.0 * sigma * (c_in - c) - Gamma * DIVGRAD(c_in);
        w_hat <<= f + 2.0 * sigma * (c_in - c);
        SetBCs_c(mesh, c_in, c);
        stormDivGrad(mesh, w_hat, -Gamma, c_in);

        // c_hat <<= c_in + tau * CONV(v, c_in) - tau * DIVGRAD(w_hat)
        c_hat <<= c_in;
        SetBCs_w(mesh, w_hat);
        stormConvection(mesh, c_hat, -tau, c_in, v);
        stormDivGrad(mesh, c_hat, -tau, w_hat);
      },
      false);

  // c_hat <<= map([](real_t c) { return std::clamp(c, 0.0, 1.0); }, c_hat);

  // w_hat = dF_dc(c_hat) - Gamma * DIVGRAD(c_hat)
  w_hat <<= map(dF_dc, c_hat);
  SetBCs_c(mesh, c_hat, c);
  stormDivGrad(mesh, w_hat, -Gamma, c_hat);

} // CahnHilliard_Step

static double mu_1 = 0.08, mu_2 = 0.08;
static double rho_1 = 1.0, rho_2 = 50.0;

static void NavierStokes_Step(stormMesh_t mesh, //
                              const StormArray<real_t>& p,
                              const StormArray<Vec2D<real_t>>& v, //
                              const StormArray<real_t>& c,
                              const StormArray<real_t>& w,
                              StormArray<real_t>& p_hat,
                              StormArray<Vec2D<real_t>>& v_hat,
                              StormArray<real_t>& rho) {
  StormArray<real_t> mu, rho_inv;
  mu.assign(rho, false);
  rho_inv.assign(rho, false);

  rho <<= map([](real_t c) noexcept { return rho_1 + (rho_2 - rho_1) * c; }, c);
  mu <<= map([](real_t c) noexcept { return mu_1 + (mu_2 - mu_1) * c; }, c);
  rho_inv <<= map([](real_t rho) noexcept { return 1.0 / rho; }, rho);

  {
    StormArray<Vec2D<real_t>> rhs;
    rhs.assign(v, false);

    // rhs <<= v + (tau * Sigma / math::sqrt(Gamma)) * rho_inv * c * GRAD(w);
    SetBCs_w(mesh, w);
    SetBCs_v(mesh, v);
    fill_with(rhs, Vec2D<real_t>{0.0, 0.0});
    stormGradient(mesh, rhs, tau, w);
    rhs *= (Sigma / std::sqrt(Gamma)) * c * rho_inv;
    rhs += v;

    v_hat <<= v;
    stormNonlinSolve2( //
        Storm::SolverType::Jfnk, v_hat, v,
        [&](StormArray<Vec2D<real_t>>& v_hat,
            const StormArray<Vec2D<real_t>>& v) {
          StormArray<Vec2D<real_t>> tmp;
          tmp.assign(v, false);

          // v_hat <<= v + tau * CONV(v, v) - tau * mu * rho_inv * DIVGRAD(v);
          v_hat <<= v;
          SetBCs_v(mesh, v);
          stormConvection(mesh, v_hat, -tau, v, v);
          fill_with(tmp, Vec2D<real_t>{0.0, 0.0});
          stormDivGrad(mesh, tmp, tau, v);
          tmp *= mu * rho_inv;
          v_hat -= tmp;
        });
  }

  v_hat <<= map(
      [](const Vec2D<real_t>& v) noexcept {
        return v + 0.3 * tau * Vec2D<real_t>{0.0, -1.0};
      },
      v_hat);

  {
    StormArray<real_t> rhs;
    rhs.assign(p, false);

    // rhs <<= p - DIV(v_hat) + ???
    rhs <<= p;
    stormDivergence(mesh, rhs, 1.0, v_hat);
    SetBCs_w(mesh, rho_inv);

    p_hat <<= p;
    stormLinSolve2(
        Storm::SolverType::Cg, Storm::PreconditionerType::None /*"extr"*/, //
        p_hat, rhs,
        [&](StormArray<real_t>& p_hat, const StormArray<real_t>& p) {
          // p_hat <<= p - tau * DIVGRAD(rho_inv, p)
          p_hat <<= p;
          SetBCs_p(mesh, p);
          stormDivWGrad(mesh, p_hat, -tau, rho_inv, p);
        },
        false);
  }

  {
    StormArray<Vec2D<real_t>> tmp;
    tmp.assign(v, false);

    // v_hat -= tau * rho_inv * GRAD(p_hat)
    fill_with(tmp, Vec2D<real_t>{0.0, 0.0});
    SetBCs_p(mesh, p_hat);
    stormGradient(mesh, tmp, tau, p_hat);
    v_hat += rho_inv * tmp;
  }

} // NavierStokes_Step

void Initial_Data(stormSize_t dim, const stormReal_t* r, stormSize_t size,
                  stormReal_t* c, const stormReal_t* _, void* env) {
  const static stormReal_t L = 1.0;
  bool in = false;
  if (fabs(r[0] - 0 * L) <= L * 0.101 && fabs(2 * L - r[1]) <= L * 0.665) {
    in = true;
  }
  // if (fabs(2 * L - r[1]) <= L * 0.4) { in = true; }

  *c = 0.0;
  if (in) { *c = 1.0; }

} // Initial_Data

void Init_For_NVT(Nvt& NVT_obj) {
  // Количество компонент
  //  int n_comp = 2;
  // Критическая температура
  std::vector<double> T_crit{304.2, 617.6};
  // Критическое давление
  std::vector<double> P_crit{73.7646e5, 21.0756e5};
  // Ацентрические факторы
  std::vector<double> ac_factors{0.225, 0.49};
  // Параметры парного взаимодействия
  std::vector<std::vector<double>> k_ij{{0, 0.1141}, {0.1141, 0}};
  // Количества компонент
  std::vector<double> N_parts{2549.336, 159.335375};
  // Объем
  double V = 1.0;
  // Температура
  double T = 295.15;
  // Давление
  double P_init = 4e6;
  // Количество узлов по phi
  int N_mesh = 101;

  NVT_obj.NVT_set_param(T_crit, P_crit, ac_factors, k_ij, N_parts, V, T, P_init,
                        N_mesh);
}

int main(int argc, char** argv) {
  print_banner();

#if 1
  UnstructuredMesh<2, 2, VovTable> mesh1{};
  read_mesh_from_tetgen(mesh1, "test/mesh/step.1.");
  // UnstructuredMesh<3, 3> mesh1{};
  // read_mesh_from_tetgen(mesh1, "test/mesh/mesh.1.");
  //     write_mesh_to_vtk(mesh1, "out/test.vtk");

  STORM_INFO_("mesh has {} edges", mesh1.num_edges());
  STORM_INFO_("mesh has {} faces", mesh1.num_faces());
  STORM_INFO_("mesh has {} cells", mesh1.num_cells());
  STORM_INFO_("mesh loaded");

  UnstructuredMesh<2, 2, CsrTable> mesh2{};
  // UnstructuredMesh<3, 3, CsrTable> mesh2{};
  mesh2.assign(std::move(mesh1));

  Vulture::visualize_mesh(mesh2);
  return 0;
#endif

#if 0
  {
    using namespace Storm;
    DenseMatrix<double> A{{10.0, 1.0, 2.0, 3.0},
                          {4.0, 11.0, 1.0, 5.0},
                          {1.0, 3.0, 13.0, 2.0},
                          {2.0, 4.0, 5.0, 25.0}};
    DenseMatrix<double> B(4, 4);
#if 0
    {
      B <<= matmul(A, A);
      A <<= matmul(A, A);
      std::cout << A << std::endl << std::endl;
      std::cout << B << std::endl << std::endl;
      std::cout << "=======" << std::endl << std::endl;
    }
#endif
    {
      inplace_inverse_lu(matmul(A, A), B);
      std::cout << A << std::endl << std::endl;
      std::cout << B << std::endl << std::endl;
      std::cout << matmul(B, matmul(A, A)) << std::endl << std::endl;
      std::cout << matmul(matmul(A, A), B) << std::endl << std::endl;
      std::cout << "=======" << std::endl << std::endl;
    }
    {
      auto AA = select_rows(A, 0, 1, 3);
      AA(0, 0) = 1.0;
      select_rows(AA, 0) -= select_rows(AA, 1);
      std::cout << AA << std::endl << std::endl;
      std::cout << "=======" << std::endl << std::endl;
    }
#if 1
    {
      auto AA = slice_rows(select_cols(A, 0, 1, 3), 0, 3, 1);
      auto BB = select_rows(select_cols(B, 0, 1, 3), 0, 1, 3);
      AA *= 1.4880;
      inplace_inverse_lu(AA, BB);
      std::cout << AA << std::endl << std::endl;
      std::cout << BB << std::endl << std::endl;
      std::cout << matmul(BB, AA) << std::endl << std::endl;
      std::cout << matmul(AA, BB) << std::endl << std::endl;
      std::cout << "=======" << std::endl << std::endl;
    }
#endif
  }
  return 0;
#endif

#if 0
  Nvt nvt_class = Nvt();

  Init_For_NVT(nvt_class);

  stormMesh_t mesh = SR_InitMesh();

  StormArray<real_t> c(mesh, SR_Alloc(mesh, 1, 0)),
      c_hat(mesh, SR_Alloc(mesh, 1, 0));
  StormArray<real_t> w_hat(mesh, SR_Alloc(mesh, 1, 0));
  StormArray<real_t> p(mesh, SR_Alloc(mesh, 1, 0)),
      p_hat(mesh, SR_Alloc(mesh, 1, 0));
  StormArray<Vec2D<real_t>> v(mesh, SR_Alloc(mesh, 1, 1)),
      v_hat(mesh, SR_Alloc(mesh, 1, 1));
  StormArray<real_t> rho(mesh, SR_Alloc(mesh, 1, 0));

  stormSpFuncProd(mesh, c, c, Initial_Data, STORM_NULL);
  fill_with(v, Vec2D<real_t>{0.0, 0.0});
  fill_with(p, 0.0);

  // Calculate Nvt beforehand


  double total_time = 0.0;

  for (int time = 0; time <= 200000; ++time) {
    for (int frac = 0; time != 0 && frac < 10; ++frac) {
      struct timespec start, finish;
      clock_gettime(CLOCK_MONOTONIC, &start);

      CahnHilliard_Step(mesh, c, v, c_hat, w_hat);
      NavierStokes_Step(mesh, p, v, c_hat, w_hat, p_hat, v_hat, rho);

      clock_gettime(CLOCK_MONOTONIC, &finish);
      double elapsed = (finish.tv_sec - start.tv_sec);
      elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
      total_time += elapsed;

      std::swap(c, c_hat);
      std::swap(p, p_hat);
      std::swap(v, v_hat);
    }

    char filename[256];
    printf("time = %f\n", total_time);
    sprintf(filename, "out/fld-%d.vti", time);
    stormIOList_t io = SR_IO_Begin();
    SR_IO_Add(io, v, "velocity");
    SR_IO_Add(io, c, "phase");
    SR_IO_Add(io, p, "pressure");
    SR_IO_Add(io, rho, "density");
    SR_IO_Flush(io, mesh, filename);
  }
#endif

  return 0;
}
