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

#include <algorithm>
#include <cstring>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define YURI 0

#include <StormRuler_API.h>

#include <Storm/Blass/Mat.hpp>
//#include <Storm/Blass/Matrix.hpp>

#include <Storm/NVT/Nvt.hpp>
#include <Storm/Solvers/PreconditionerFactory.hpp>
#include <Storm/Solvers/SolverFactory.hpp>

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

  auto solver = MakeIterativeSolver<StormArray<T>>(method);
  solver->pre_op = make_preconditioner<StormArray<T>>(preMethod);

  if (uniform) {
    solver->solve(x, b, *symOp);
  } else {
    solve_non_uniform(*solver, x, b, *symOp);
  }

  std::cout << "num matvecs = " << numMatVecs << " " << method.ToString()
            << std::endl;

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

// tau - time step
// Gamma - lambda^2
// sigma - const stab (num meth)
// Sigma - ??
// Try Sigma = 1??
static double tau = 1.0e-2, Gamma = 1.0e-4, sigma = 1.0, Sigma = 10.0; // 10.0;
static double Mobility;

static void SetBCs_c(stormMesh_t mesh, stormArray_t c_hat, stormArray_t c,
                     const std::vector<double>& BV) {
  stormApplyBCs(mesh, c_hat, SR_ALL, SR_PURE_NEUMANN);
  // stormApplyBCs_CosWall(
  //     mesh, c_hat, c,
  //     0.01 / std::sqrt(Gamma) * 6.0 * std::cos(M_PI / 2 - M_PI / 18), 1);
  //  stormApplyBCs(mesh, c_hat, 4, SR_DIRICHLET(1.0));
} // SetBCs_c

static void SetBCs_w(stormMesh_t mesh, stormArray_t w,
                     const std::vector<double>& BV) {
  stormApplyBCs(mesh, w, SR_ALL, SR_PURE_NEUMANN);
} // SetBCs_w

static void SetBCs_p(stormMesh_t mesh, stormArray_t p,
                     const std::vector<double>& BV) {
  stormApplyBCs(mesh, p, SR_ALL, SR_PURE_NEUMANN);
  // stormApplyBCs(mesh, p, 2, SR_DIRICHLET(0.0));
  // stormApplyBCs(mesh, p, 4,
  //               SR_DIRICHLET(0.75 * Sigma * std::cos(M_PI / 2 - M_PI / 18) /
  //                            (2.0 * 0.01 * 26)));
} // SetBCs_p

static void SetBCs_v(stormMesh_t mesh, stormArray_t v,
                     const std::vector<double>& BV) {
  stormApplyBCs(mesh, v, SR_ALL, SR_PURE_DIRICHLET);
  stormApplyBCs_SlipWall(mesh, v, 3);
  // stormApplyBCs(mesh, v, 2, SR_PURE_NEUMANN);
  // stormApplyBCs_InOutLet(mesh, v, 2);
  // stormApplyBCs(mesh, v, 4, SR_PURE_NEUMANN);

  [[maybe_unused]] stormReal_t R = 0.2 * 0.5;
  [[maybe_unused]] stormReal_t qFlux = 10.0 * 0.5 * M_PI * R * R * R * R;
  stormApplyBCs_InOutLet(mesh, v, 4, BV[4]);
} // SetBCs_v

static void CahnHilliard_Step(stormMesh_t mesh, //
                              const StormArray<real_t>& c,
                              const StormArray<Vec2D<real_t>>& v, //
                              StormArray<real_t>& c_hat,
                              StormArray<real_t>& w_hat,
                              const std::vector<double>& BV) {
  SetBCs_c(mesh, c, c, BV);
  SetBCs_v(mesh, v, BV);

  constexpr auto dF_dc = [](real_t c) {
    return 2.0 * c * (c - 1.0) * (2.0 * c - 1.0);
  };

  StormArray<real_t> f;
  f.assign(c, false);

  f <<= map(dF_dc, c);

  c_hat <<= c;
  stormLinSolve2(
      Storm::SolverType::BiCgStab, Storm::PreconditionerType::None /*"extr"*/,
      c_hat, c,
      [&](StormArray<real_t>& c_hat, const StormArray<real_t>& c_in) {
        // w_hat <<= f + 2.0 * sigma * (c_in - c) - Gamma * DIVGRAD(c_in);
        w_hat <<= f + 2.0 * sigma * (c_in - c);
        SetBCs_c(mesh, c_in, c, BV);
        stormDivGrad(mesh, w_hat, -Gamma, c_in);

        // c_hat <<= c_in + tau * CONV(v, c_in) - tau * DIVGRAD(w_hat)
        c_hat <<= c_in;
        SetBCs_w(mesh, w_hat, BV);
        stormConvection(mesh, c_hat, -tau, c_in, v);
        stormDivGrad(mesh, c_hat, -tau * Mobility, w_hat);
      },
      false);

  c_hat <<= map([](real_t c) { return std::clamp(c, 0.0, 1.0); }, c_hat);

  // w_hat = dF_dc(c_hat) - Gamma * DIVGRAD(c_hat)
  w_hat <<= map(dF_dc, c_hat);
  SetBCs_c(mesh, c_hat, c, BV);
  stormDivGrad(mesh, w_hat, -Gamma, c_hat);

} // CahnHilliard_Step

static double mu_1 = 0.008, mu_2 = 0.008;
static double rho_1 = 1.0, rho_2 = 50.0;

static void
NavierStokes_Step(stormMesh_t mesh, //
                  const StormArray<real_t>& p,
                  const StormArray<Vec2D<real_t>>& v, //
                  const StormArray<real_t>& c, const StormArray<real_t>& w,
                  StormArray<real_t>& p_hat, StormArray<Vec2D<real_t>>& v_hat,
                  StormArray<real_t>& rho, const std::vector<double>& BV) {
  StormArray<real_t> mu, rho_inv;
  mu.assign(rho, false);
  rho_inv.assign(rho, false);

  rho <<= map([](real_t c) { return rho_1 + (rho_2 - rho_1) * c; }, c);
  mu <<= map([](real_t c) { return mu_1 + (mu_2 - mu_1) * c; }, c);
  rho_inv <<= map([](real_t rho) { return 1.0 / rho; }, rho);

  {
    StormArray<Vec2D<real_t>> rhs;
    rhs.assign(v, false);

    // rhs <<= v + (tau * Sigma / math::sqrt(Gamma)) * rho_inv * c * GRAD(w);
    SetBCs_w(mesh, w, BV);
    SetBCs_v(mesh, v, BV);
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
          SetBCs_v(mesh, v, BV);
          stormConvection(mesh, v_hat, -tau, v, v);
          fill_with(tmp, Vec2D<real_t>{0.0, 0.0});
          stormDivGrad(mesh, tmp, tau, v);
          tmp *= mu * rho_inv;
          v_hat -= tmp;
        });

    // std::cout<<"NonlinSolve2 passed"<<std::endl;
  }

  v_hat <<= map(
      [](const Vec2D<real_t>& v) {
        return v + 0.0 * tau * Vec2D<real_t>{0.0, -1.0}; // no gravity!
      },
      v_hat);

  {
    StormArray<real_t> rhs;
    rhs.assign(p, false);

    // rhs <<= p - DIV(v_hat) + ???
    rhs <<= p;
    stormDivergence(mesh, rhs, 1.0, v_hat);
    SetBCs_w(mesh, rho_inv, BV);

    p_hat <<= p;
    stormLinSolve2(
        Storm::SolverType::Cg, Storm::PreconditionerType::None /*"extr"*/, //
        p_hat, rhs,
        [&](StormArray<real_t>& p_hat, const StormArray<real_t>& p) {
          // p_hat <<= p - tau * DIVGRAD(rho_inv, p)
          p_hat <<= p;
          SetBCs_p(mesh, p, BV);
          stormDivWGrad(mesh, p_hat, -tau, rho_inv, p);
        },
        false);
    // std::cout<<"LinSolve2 passed"<<std::endl;
  }

  {
    StormArray<Vec2D<real_t>> tmp;
    tmp.assign(v, false);

    // v_hat -= tau * rho_inv * GRAD(p_hat)
    fill_with(tmp, Vec2D<real_t>{0.0, 0.0});
    SetBCs_p(mesh, p_hat, BV);
    stormGradient(mesh, tmp, tau, p_hat);
    v_hat += rho_inv * tmp;
  }

} // NavierStokes_Step

void Initial_Data(stormSize_t dim, const stormReal_t* r, stormSize_t size,
                  stormReal_t* c, const stormReal_t* _, void* env) {
  // const static stormReal_t L = 1.0;
  const static stormReal_t L = 1e-2;
  const static stormReal_t l = 1e-3;
  const static stormReal_t h_small = 5e-3 / 120.0;
  bool in = false;
  if ((fabs(r[0]) <= 0.101 * 5e-3) && (r[1] >= L + h_small) &&
      (r[1] <= L + l + 2 * h_small)) {
    // if (fabs(r[0] - 0 * L) <= L * 0.101 && fabs(2 * L - r[1]) <= L * 0.665) {
    in = true;
  }
  // if (fabs(2 * L - r[1]) <= L * 0.4) { in = true; }

  *c = 0.0;
  if (in) { *c = 1.0; }

} // Initial_Data

void Init_For_NVT(Nvt& NVT_obj) {
  //Количество компонент
  // 0 - CO2, 1 - C10
  // int n_comp = 2;
  //Критическая температура
  std::vector<double> T_crit{304.2, 617.6};
  //Критическое давление
  std::vector<double> P_crit{73.7646e5, 21.0756e5};
  //Ацентрические факторы
  std::vector<double> ac_factors{0.225, 0.49};
  //Параметры парного взаимодействия
  std::vector<std::vector<double>> k_ij{{0, 0.1141}, {0.1141, 0}};
  //Количества компонент
  std::vector<double> N_parts{2549.336, 159.335375};
  //Объем
  double V = 1.0;
  //Температура
  double T = 295.15;
  //Давление
  double P_init = 5e6;
  //Количество узлов по phi
  int N_mesh = 101;

  NVT_obj.NVT_set_param(T_crit, P_crit, ac_factors, k_ij, N_parts, V, T, P_init,
                        N_mesh);
}

int main(int argc, char** argv) {
  // //mesh file name
  const char* mesh_filename = "./test/Domain-100-Tube_brandnew.ppm";


  // global parameters of the task
  [[maybe_unused]] const stormReal_t R_area = 5e-3; // main reservoir radius
  [[maybe_unused]] const stormReal_t L_to_R =
      2.0; // Length to radius ration (main reservoir)
  [[maybe_unused]] const stormReal_t R_to_r =
      10.0; // Main reservoir radius to small reservoir radius
  [[maybe_unused]] const stormReal_t TaskTime = 10000.0; // Total soultion time
  [[maybe_unused]] const stormReal_t P_form =
      5e6; // Formation pressure (initial?)
  [[maybe_unused]] const stormReal_t Ving_to_V =
      0.1; // V_injection to total volume
  [[maybe_unused]] const stormReal_t TimesScale =
      100.0; // Time factor ratio (= TaskTime / TimeRelax = TimeRelax / dt)
  [[maybe_unused]] const stormReal_t Lambda_to_h =
      5.0; // Interface width in cells
  [[maybe_unused]] const stormInt_t N_mesh_R =
      120; // Number of cells for the main radius
  [[maybe_unused]] const stormInt_t Wall_N_r =
      2; // Wall mesh size in r dimension
  [[maybe_unused]] const stormInt_t Wall_N_l =
      8; // Wall mesh size in l dimension

  // Calculate dependable parameters
  [[maybe_unused]] stormReal_t R_small = R_area / R_to_r;
  [[maybe_unused]] stormReal_t L_area = R_area * L_to_R;
  [[maybe_unused]] stormReal_t L_small = L_area / R_to_r;
  [[maybe_unused]] stormReal_t V_main = M_PI * R_area * R_area * L_area;
  [[maybe_unused]] stormReal_t V_small = M_PI * R_small * R_small * L_small;
  [[maybe_unused]] stormReal_t V_ing = V_main * Ving_to_V;
  [[maybe_unused]] stormReal_t tau_relax = TaskTime / TimesScale;
  tau = tau_relax / TimesScale;
  [[maybe_unused]] stormReal_t Qflux = V_ing / TaskTime;
  [[maybe_unused]] stormReal_t dR = R_area / N_mesh_R;
  // stormReal_t dR = 0.01;
  [[maybe_unused]] stormReal_t Lambda = Lambda_to_h * dR;
  /*static stormReal_t */ Mobility = Lambda * Lambda / tau_relax;

  // Qflux /= 10.0;

  // double u_mean = Qflux / M_PI / (R_small * R_small);
  // double Re = u_mean * R_small * rho_2 / mu_2;
  // std::cout << Re << std::endl;

  Gamma = Lambda * Lambda;

  // Values for boundaries; index -> value for boundary with <<index>> number
  std::vector<double> BoundaryKeyValues(6);
  BoundaryKeyValues[0] = 0.0;
  BoundaryKeyValues[1] = 0.0;
  BoundaryKeyValues[2] = 0.0;
  BoundaryKeyValues[3] = 0.0;
  BoundaryKeyValues[4] = Qflux; //поток на границе типа вход
  BoundaryKeyValues[5] = 0.0;

  // Initialize NVT
  Nvt nvt_class = Nvt();

  Init_For_NVT(nvt_class);

  // Calculate initial concetrations and quantities
  // CO2
  [[maybe_unused]] stormReal_t n_conc_CO2 = nvt_class.Concentration(0, P_form);
  [[maybe_unused]] stormReal_t N_CO2 = n_conc_CO2 * V_main;
  //[[maybe_unused]] //C10
  [[maybe_unused]] stormReal_t n_conc_C10 = nvt_class.Concentration(1, P_form);
  [[maybe_unused]] stormReal_t N_C10 = n_conc_C10 * V_small;

  // need to give these N to NVT!!!
  std::vector<double> N_load;
  N_load.push_back(N_CO2 / (V_main + V_small));
  N_load.push_back(N_C10 / (V_main + V_small));
  nvt_class.Load_N(N_load);

  // Formation pressure

  // make name a function parameter!
  stormMesh_t mesh = SR_InitMesh(mesh_filename, dR, dR);

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
  double time_of_task = 0.0;
  int num_write = 0;

  // for (int time = 0; time <= 200000; ++time) {
  while (time_of_task < TaskTime) {
    for (int frac = 0; num_write != 0 && frac < 1; ++frac) {
      struct timespec start, finish;
      clock_gettime(CLOCK_MONOTONIC, &start);

      CahnHilliard_Step(mesh, c, v, c_hat, w_hat, BoundaryKeyValues);
      NavierStokes_Step(mesh, p, v, c_hat, w_hat, p_hat, v_hat, rho,
                        BoundaryKeyValues);

      clock_gettime(CLOCK_MONOTONIC, &finish);
      double elapsed = (finish.tv_sec - start.tv_sec);
      elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
      total_time += elapsed;
      time_of_task += tau;

      std::swap(c, c_hat);
      std::swap(p, p_hat);
      std::swap(v, v_hat);
    }

    char filename[256];
    printf("time = %f\n", total_time);
    printf("time of task = %f\n", time_of_task);
    sprintf(filename, "out/fld-%d.vti", num_write);
    stormIOList_t io = SR_IO_Begin();
    SR_IO_Add(io, v, "velocity");
    SR_IO_Add(io, c, "phase");
    SR_IO_Add(io, p, "pressure");
    SR_IO_Add(io, rho, "density");
    SR_IO_Flush(io, mesh, filename);

    num_write++;
  }

  return 0;
}
