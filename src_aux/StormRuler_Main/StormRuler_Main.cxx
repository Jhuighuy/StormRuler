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

#include <StormRuler_API.h>
#include <stormBlas/stormTensor.hxx>
//#include <stormBlas/stormMatrixExtraction.hxx>
#include <stormSolvers/SolverFactory.hxx>
#include <stormSolvers/PreconditionerFactory.hxx>

//#include <stormMesh/Mesh.hxx>

#include <stormUtils/JsonValue.hxx>

#include <cstring>

template<typename stormMatVecFuncT_t>
STORM_INL void stormLinSolve2(stormMesh_t mesh,
                              Storm::SolverType const& method,
                              Storm::PreconditionerType const& preMethod,
                              stormArray_t x,
                              stormArray_t b,
                              stormMatVecFuncT_t matVec,
                              bool uniform = true) {
  using namespace Storm;

#if 0
  Array<real_t> xVec(14674, 15906), bVec(xVec);

  Storm::stormArray xx = {mesh, x}, bb = {mesh, b};
  ToArray(xVec, xx);
  ToArray(bVec, bb);

  stormSize_t numMatVecs = 0;
  Solve(
    method, preMethod,
    xVec, bVec,
    [&](Array<real_t>& yVec, Array<real_t> const& zVec) {
     numMatVecs += 1;
      stormArray_t y = stormAllocLike(x), z = stormAllocLike(b);
      Storm::stormArray yy = {mesh, y}, zz = {mesh, z};
      FromArray(yy, yVec), FromArray(zz, zVec);
      matVec(mesh, y, z);
      ToArray(yVec, yy);
      stormFree(y), stormFree(z);
    });
#else
  Storm::stormArray xx = {mesh, x}, bb = {mesh, b};
  stormSize_t numMatVecs = 0;
  auto symOp = Storm::MakeSymmetricOperator<Storm::stormArray>(
    [&](Storm::stormArray& yy, const Storm::stormArray& xx) {
      numMatVecs += 1;
      matVec(yy.Mesh, yy.Array, xx.Array);
    });

  auto solver = Storm::MakeIterativeSolver<Storm::stormArray>(method);
  solver->PreOp = Storm::MakePreconditioner<Storm::stormArray>(preMethod);

  if (uniform) {
    solver->Solve(xx, bb, *symOp);
  } else {
    Storm::SolveNonUniform(*solver, xx, bb, *symOp);
  }
#endif
  std::cout << "num matvecs = " << numMatVecs << ' ' <<  method.ToString() << std::endl;

} // stormLinSolve2<...>

template<typename stormMatVecFuncT_t>
STORM_INL void stormNonlinSolve2(stormMesh_t mesh,
                                 Storm::SolverType const& method,
                                 stormArray_t x,
                                 stormArray_t b, 
                                 stormMatVecFuncT_t matVec) {

  Storm::stormArray xx = {mesh, x}, bb = {mesh, b};
  auto op = Storm::MakeSymmetricOperator<Storm::stormArray>(
      [&](Storm::stormArray& yy, const Storm::stormArray& xx) {
        matVec(yy.Mesh, yy.Array, xx.Array);
      });

  std::unique_ptr<Storm::IterativeSolver<Storm::stormArray>> solver = 
    std::make_unique<Storm::JfnkSolver<Storm::stormArray>>();
  solver->AbsoluteTolerance = 1.0e-4;
  solver->RelativeTolerance = 1.0e-4;
  solver->Solve(xx, bb, *op);

} // stormNonLinSolve2<...>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <fstream>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define YURI 0

#define min(x, y) ( (x) < (y) ? (x) : (y) )
#define max(x, y) ( (x) > (y) ? (x) : (y) )

static double tau = 1.0e-2, Gamma = 1.0e-4, sigma = 1.0;

static void SetBCs_c(stormMesh_t mesh, stormArray_t c_hat, stormArray_t c) {
  stormApplyBCs(mesh, c_hat, SR_ALL, SR_PURE_NEUMANN);
  stormApplyBCs_CosWall(mesh, c_hat, c, 
  
    0.01/std::sqrt(Gamma)*6.0*std::cos(M_PI/2 - M_PI/16)
    
    , 1);
  stormApplyBCs(mesh, c_hat, 4, SR_DIRICHLET(1.0));
} // SetBCs_c

static void SetBCs_w(stormMesh_t mesh, stormArray_t w) {
  stormApplyBCs(mesh, w, SR_ALL, SR_PURE_NEUMANN);
} // SetBCs_w

static void SetBCs_p(stormMesh_t mesh, stormArray_t p) {
  stormApplyBCs(mesh, p, SR_ALL, SR_PURE_NEUMANN);
  //stormApplyBCs(mesh, p, 2, SR_DIRICHLET(1.0));
  //stormApplyBCs(mesh, p, 4, SR_DIRICHLET(2.0));
} // SetBCs_p

static void SetBCs_v(stormMesh_t mesh, stormArray_t v) {
  stormApplyBCs(mesh, v, SR_ALL, SR_PURE_DIRICHLET);
  stormApplyBCs_SlipWall(mesh, v, 3);
#if !YURI
  stormApplyBCs_InOutLet(mesh, v, 2);
#endif
  stormApplyBCs_InOutLet(mesh, v, 4);
} // SetBCs_v

#if YURI
stormSize_t slice = 983 + 1, num_slices = 1000;
stormVector<double> phi_set;
stormVector<double> alpha_sandwich;
stormMatrix<double> D1_W_vs_phi_sandwich;
stormTensor2R<double> N_part_sandwich;
stormTensor3R<double> n_part_vs_phi_sandwich;
double N2_min, N2_max, N2_cur;
static double const mol_mass[2] = {0.0440098, 0.1422853};

template<class Tensor>
void ReadTensor(Tensor& tensor, std::string path) {
  path = "/home/jhuighuy/Downloads/decane_droplet/" + path;
  std::ifstream ifstream(path);
  std::for_each_n(tensor.Data(), tensor.Size(), [&](auto& value) {
    ifstream >> value;
  });
}
#endif

void dWdC(stormSize_t size, stormReal_t* Wc, const stormReal_t* c, void* env) {
  const stormReal_t x = *c;

#if !YURI

  // [-1,+1] CH.
  //*Wc = (x < -1.0) ? 2.0*(1.0+x) : ( (x > 1.0) ? (2.0*(x-1.0)) : x*(x*x - 1.0) );
  // [0,1] CH.
  *Wc = (x < -1.0) ? 2.0*x : ( (x > 1.0) ? (2.0*(x-1.0)) : 2.0*x*(x - 1.0)*(2.0*x - 1.0) );

#else

  auto const D1_W_vs_phi = [&](auto i) { return D1_W_vs_phi_sandwich(slice, i) / 32.0; };

  // [0,1] MCH.
  const double h = 1.0/(phi_set.Size() - 1.0);
  if (x < 0.0) {
    *Wc = ( D1_W_vs_phi(0) + x*(D1_W_vs_phi(1) - D1_W_vs_phi(0))/h );
    return;
  }
  if (x > 1.0) {
    *Wc = ( D1_W_vs_phi(100) + (x - 1.0)*(D1_W_vs_phi(100) - D1_W_vs_phi(99))/h );
    return;
  }
  const int il = floor(x/h);
  const int ir = ceil(x/h);
  if (il == ir) {
    *Wc = ( D1_W_vs_phi(il) );  
  }

  *Wc = ( D1_W_vs_phi(il) + (x - h*il)*(D1_W_vs_phi(ir) - D1_W_vs_phi(il))/h );
  return;

#endif

} // dWdC

void Vol(stormSize_t size, stormReal_t* Ic, const stormReal_t* c, void* env) {
  stormReal_t x = *c;

  *Ic = x;

} // Vol

#if 1
static void CahnHilliard_Step(stormMesh_t mesh,
    stormArray_t c, stormArray_t v,
    stormArray_t c_hat, stormArray_t w_hat
#if YURI
    , double& V, double& V_hat
#endif
    ) {

  SetBCs_c(mesh, c, c);
  SetBCs_v(mesh, v);

  stormArray_t tmp = stormAllocLike(c);
  stormFuncProd(mesh, tmp, c, dWdC, STORM_NULL);

  stormSet(mesh, c_hat, c);
  stormLinSolve2(mesh,
#if YURI
      Storm::SolverType::BiCgStab,
      Storm::PreconditionerType::None/*"extr"*/, 
#else
      Storm::SolverType::Tfqmr1,
      Storm::PreconditionerType::None/*"extr"*/,
#endif
    c_hat, c,
    [&](stormMesh_t mesh, stormArray_t c_out, stormArray_t c_in) {

      stormSet(mesh, w_hat, tmp);
      stormAdd(mesh, w_hat, w_hat, c_in, 2.0*sigma);
      stormSub(mesh, w_hat, w_hat, c, 2.0*sigma);

      SetBCs_c(mesh, c_in, c);
      stormDivGrad(mesh, w_hat, -Gamma, c_in);

      stormSet(mesh, c_out, c_in);
      stormConvection(mesh, c_out, -tau, c_in, v);

      SetBCs_w(mesh, w_hat);
      stormDivGrad(mesh, c_out, -tau, w_hat);

    }, false);
    abort();

  stormFree(tmp);

  SetBCs_c(mesh, c_hat, c);
  stormFuncProd(mesh, w_hat, c_hat, dWdC, STORM_NULL);
  stormDivGrad(mesh, w_hat, -Gamma, c_hat);

#if YURI
  V = V_hat;
  V_hat = stormIntegrate(mesh, c_hat, Vol, STORM_NULL);
#endif

} // CahnHilliard_Step
#else
static void CahnHilliard_Step(stormMesh_t mesh,
    stormArray_t c, stormArray_t v,
    stormArray_t c_hat, stormArray_t w_hat
#if YURI
    , double& V, double& V_hat
#endif
    ) {

  // 
  // Compute a single time step of the
  // Cahn-Hilliard equation with convection term:
  //
  // ∂𝑐/∂𝑡 + ∇⋅𝑐𝒗 = Δ𝑤, 
  // 𝑤 = 𝑊'(𝑐) - 𝛾Δ𝑐, 𝑊(𝑐) = ¼(1 - 𝑐²)²,
  //
  // with the linear-implicit scheme:
  // 
  // 𝓠𝑐̂ = 𝑐̂ + 𝜏∇⋅̂𝑐𝒗 - 𝜏Δ(2𝜎𝑐̂ - 𝛾Δ𝑐̂) = 𝑐 + 𝜏Δ[𝑊'(𝑐) - 2𝜎𝑐],
  // 𝑤̂ = 𝑊'(𝑐̂) - 𝛾Δ𝑐̂.
  //

  SetBCs_c(mesh, c);
  SetBCs_v(mesh, v);

  stormArray_t rhs = stormAllocLike(c);
  stormSet(mesh, rhs, c);
  stormFuncProd(mesh, w_hat, c, dWdC, STORM_NULL);
  stormSub(mesh, w_hat, w_hat, c, 2.0*sigma);
  SetBCs_w(mesh, w_hat);
  stormDivGrad(mesh, rhs, tau, w_hat);

  stormSet(mesh, c_hat, c);
  stormLinSolve2(mesh,
#if YURI
      Storm::SolverType::BiCgStab,
      Storm::PreconditionerType::None/*"extr"*/, 
#else
      Storm::SolverType::BiCgStab,
      Storm::PreconditionerType::None/*"extr"*/,
#endif
    c_hat, rhs,
    [&](stormMesh_t mesh, stormArray_t Qc, stormArray_t c) {
      SetBCs_c(mesh, c);
      SetBCs_v(mesh, v);

      stormArray_t tmp = stormAllocLike(c);

      stormScale(mesh, tmp, c, 2.0*sigma);
      stormDivGrad(mesh, tmp, -Gamma, c);

      SetBCs_w(mesh, tmp);
      
      stormSet(mesh, Qc, c);
      stormConvection(mesh, Qc, -tau, c, v);
      stormDivGrad(mesh, Qc, -tau, tmp);

      stormFree(tmp);
    });
    //abort();
  stormFree(rhs);

  SetBCs_c(mesh, c_hat);
  stormFuncProd(mesh, w_hat, c_hat, dWdC, STORM_NULL);
  stormDivGrad(mesh, w_hat, -Gamma, c_hat);

#if YURI
  V = V_hat;
  V_hat = stormIntegrate(mesh, c_hat, Vol, STORM_NULL);
#endif

} // CahnHilliard_Step
#endif

static double mu_1 = 0.08, mu_2 = 0.08;
#if !YURI
static double rho_1 = 1.0, rho_2 = 50.0;
#endif

void InvRho(stormSize_t size, stormReal_t* inv_rho, const stormReal_t* rho, void* env) {
  *inv_rho = 1.0/(*rho);
} // InvRho

#if YURI
static int II;
void NVsC(stormSize_t size, stormReal_t* n, const stormReal_t* c, void* env) {
  stormReal_t x = *c;

  auto const nPart_vs_phi = [&](auto i, auto j) { return n_part_vs_phi_sandwich(slice, i, j); };

  x = max(0.0, min(1.0, x));

  const int i = II;
  double dd;
  const double h = 0.01;
  if (x < 0.0) {
    dd = ( nPart_vs_phi(0, i) + x*(nPart_vs_phi(1, i) - nPart_vs_phi(0, i))/h );
  }
  else if (x > 1.0) {
    dd = ( nPart_vs_phi(100, i) + (x - 1.0)*(nPart_vs_phi(100, i) - nPart_vs_phi(99, i))/h );
  } else {
    const int il = floor(x/0.01);
    const int ir = ceil(x/0.01);
    if (il == ir) {
      dd = ( nPart_vs_phi(il, i) );  
    } else {
      dd = ( nPart_vs_phi(il, i) + (x - il*h)*(nPart_vs_phi(ir, i) - nPart_vs_phi(il, i))/h );
    }
  }

  *n = dd; 
  return;

} // NVsC
#endif

void AddGravity(stormSize_t size, stormReal_t* rhs, const stormReal_t* rhs_, void* env) {
#if YURI
  //rhs[1] += tau*(-0.3);
#else
  //rhs[1] += tau*(-1.0);
#endif
}

static void NavierStokes_VaD_Step(stormMesh_t mesh,
  stormArray_t p, stormArray_t v,
  stormArray_t c, stormArray_t w,
  stormArray_t p_hat, stormArray_t v_hat, stormArray_t d, stormArray_t rho
#if YURI
  , stormArray_t n1, stormArray_t n2, double V, double V_hat
#endif
) {

  // 
  // Compute a single time step of the incompressible
  // Navier-Stokes equation:
  //
  // 𝜌(∂𝒗/∂𝑡 + 𝒗(∇⋅𝒗)) + ∇𝑝 = 𝜇Δ𝒗 + 𝙛,
  // 𝜌 = 𝜌₁(1 - 𝑐) + 𝜌₂*𝑐,
  // 𝜇 = 𝜇₁(1 - 𝑐) + 𝜇₂*𝑐,
  // ∇⋅𝒗 = 0, 𝙛 = -𝑐∇𝑤,
  //
  // with the semi-implicit scheme:
  // 
  // 𝜌 ← 𝜌₁ + (𝜌₂ - 𝜌₁)𝑐,
  // 𝜇 ← 𝜇₁ + (𝜇₂ - 𝜇₁)𝑐,
  // 𝒗̂ + 𝜏𝒗̂(∇⋅𝒗̂) - (𝜏𝜇/𝜌)Δ𝒗̂ = 𝒗 + (𝜏/𝜌)𝙛,
  // 𝑝̂ - 𝜏∇⋅(∇𝑝̂/𝜌) = 𝑝 - ∇⋅𝒗,
  // 𝒗̂ ← 𝒗̂ - (𝜏/𝜌)∇𝑝̂.
  //

#if YURI && 0
  static bool first = true;
  if (first) {
    first = false;
    std::ofstream slice_info;
    double alpha = alpha_sandwich(slice);
    slice_info.open("slices.txt");
    slice_info << "Starting from slice " << slice << std::endl;
    slice_info << "V_hat = " << V_hat << " alpha = " << alpha << std::endl; 
    slice_info << std::endl;
  }
  //double alpha_dot = (V_hat - V)/tau;
  double A = 0.0133/4.55;
  N2_cur = N2_cur + tau*n_part_vs_phi_sandwich(slice, 100, 1)*A;
  N2_cur = max(min(N2_cur, N2_max), N2_min);
  auto new_slice = (stormSize_t)std::round((N2_cur - N2_min)/(N2_max - N2_min)*(num_slices - 1));
  new_slice = num_slices - new_slice - 1;
  if (new_slice != slice) {
    std::ofstream slice_info;
    double alpha = alpha_sandwich(new_slice);
    slice_info.open("slices.txt", std::ios_base::app);
    slice_info << "Change from slice " << slice << " to " << new_slice << std::endl;
    slice_info << "V_hat = " << V_hat << " alpha = " << alpha << std::endl;
    slice_info << std::endl;
  }
  slice = new_slice;
#endif

  //
  // Compute 𝜌, 𝜇, 1/𝜌.
  //
#if YURI
  II = 0; stormFuncProd(mesh, n1, c, NVsC, STORM_NULL);
  II = 1; stormFuncProd(mesh, n2, c, NVsC, STORM_NULL);
  stormAdd(mesh, rho, n1, n2, mol_mass[1], mol_mass[0]);
#else
  stormFill(mesh, rho, rho_1);
  stormAdd(mesh, rho, rho, c, rho_2 - rho_1);
#endif

  stormArray_t rho_inv = stormAllocLike(rho);
  stormFuncProd(mesh, rho_inv, rho, InvRho, STORM_NULL);

  stormArray_t mu = stormAllocLike(c);
  stormFill(mesh, mu, mu_1);
  stormAdd(mesh, mu, mu, c, mu_2 - mu_1);

  //
  // Compute 𝒗̂ prediction.
  //
  SetBCs_w(mesh, w);
  SetBCs_v(mesh, v);

  stormArray_t rhs = stormAllocLike(v);

  stormFill(mesh, rhs, 0.0);
  stormGradient(mesh, rhs, tau, w);
  stormMul(mesh, rhs, c, rhs);
  stormMul(mesh, rhs, rho_inv, rhs);

  stormAdd(mesh, rhs, rhs, v);
  
  stormSet(mesh, v_hat, v);
  stormNonlinSolve2(mesh, Storm::SolverType::Jfnk, v_hat, v, 
    [&](stormMesh_t mesh, stormArray_t Av, stormArray_t v) {
      SetBCs_v(mesh, v);

      stormSet(mesh, Av, v);
      stormConvection(mesh, Av, -tau, v, v);

      stormArray_t tmp = stormAllocLike(v);

      stormFill(mesh, tmp, 0.0);
      stormDivGrad(mesh, tmp, tau, v);
      stormMul(mesh, tmp, mu, tmp);
      stormMul(mesh, tmp, rho_inv, tmp);
      stormSub(mesh, Av, Av, tmp, 1.0, 1.0);

      stormFree(tmp);
    });
  
  stormFree(rhs);

  stormFuncProd(mesh, v_hat, v_hat, AddGravity, STORM_NULL);

  //
  // Solve pressure equation and correct 𝒗̂.
  // 
  rhs = stormAllocLike(p);
  
  stormSet(mesh, rhs, p);
  stormDivergence(mesh, rhs, 1.0, v_hat);
  SetBCs_w(mesh, rho);
  SetBCs_w(mesh, rho_inv);
  stormRhieChowCorrection(mesh, rhs, 1.0, tau, p, rho);

  stormSet(mesh, p_hat, p);
  stormLinSolve2(mesh, 
    Storm::SolverType::Cg,
    Storm::PreconditionerType::None/*"extr"*/,
    p_hat, rhs,
    [&](stormMesh_t mesh, stormArray_t Lp, stormArray_t p) {
      SetBCs_p(mesh, p);

      stormSet(mesh, Lp, p);
      stormDivWGrad(mesh, Lp, -tau, rho_inv, p);
    }, false);
    //abort();

  stormFree(rhs);

  SetBCs_p(mesh, p_hat);

  stormArray_t tmp = stormAllocLike(v);
  stormFill(mesh, tmp, 0.0);
  stormGradient(mesh, tmp, tau, p_hat);
  stormMul(mesh, tmp, rho_inv, tmp);
  stormAdd(mesh, v_hat, v_hat, tmp);
  stormFree(tmp);

  if (d != STORM_NULL) {
    SetBCs_v(mesh, v_hat);
    stormFill(mesh, d, 0.0);
    stormDivergence(mesh, d, -1.0, v);
  }

  stormFree(rho_inv);
  stormFree(mu);

} // NavierStokes_VaD_Step

void Initial_Data(stormSize_t dim, const stormReal_t* r,
    stormSize_t size, stormReal_t* c, const stormReal_t* _, void* env) {

  static const stormReal_t L = 1.0;
  bool in = false;
  //if (fabs(r[0]-0*L) <= L*0.101 && fabs(2*L-r[1]) <= L*0.665) {
  //  in = true;
  //}
  if (fabs(2*L-r[1]) <= L*0.2) {
    in = true;
  }

  *c = 0.0;
  if (in) {
    *c = 1.0;
  }

} // Initial_Data

int main_turbo();

int main(int argc, char** argv) {

  if (argc > 1 && strcmp(argv[1], "turbo") == 0) {
    return main_turbo();
  }

#if YURI
  phi_set.Assign(101);
  ReadTensor(phi_set, "phi.csv");
  alpha_sandwich.Assign(num_slices);
  ReadTensor(alpha_sandwich, "alpha_sandwich.csv");
  D1_W_vs_phi_sandwich.Assign(num_slices, 101);
  ReadTensor(D1_W_vs_phi_sandwich, "D1_W_vs_phi_sandwich.csv");
  N_part_sandwich.Assign(num_slices, 2);
  ReadTensor(N_part_sandwich, "N_part_sandwich.csv");
  n_part_vs_phi_sandwich.Assign(num_slices, 101, 2);
  ReadTensor(n_part_vs_phi_sandwich, "n_part_vs_phi_sandwich.csv");
  N2_max = N_part_sandwich(0, 1);
  N2_min = N_part_sandwich(num_slices - 1, 1);
  N2_cur = N_part_sandwich(slice, 1);
#endif

  //FILE* dWdCF = fopen("dWdC.txt", "w");
  //for (double x = -0.1; x <= 1.1; x += 0.0001) {
  //  double y;
  //  dWdC(1, &y, &x, STORM_NULL);
  //  fprintf(dWdCF, "%f %f\n", x, y);
  //}
  //fclose(dWdCF);
  //return 1;

#if YURI
  FILE* volFile = fopen("vol.txt", "w");
#endif

  stormMesh_t mesh = SR_InitMesh();

  stormArray_t c, p, v, c_hat, w_hat, p_hat, v_hat, d, rho;
  c = SR_Alloc(mesh, 1, 0);
  p = SR_Alloc(mesh, 1, 0);
  v = SR_Alloc(mesh, 1, 1);
  d = SR_Alloc(mesh, 1, 0);
  rho = SR_Alloc(mesh, 1, 0);
  c_hat = stormAllocLike(c);
  w_hat = stormAllocLike(c);
  p_hat = stormAllocLike(p);
  v_hat = stormAllocLike(v);
#if YURI
  stormArray_t n1, n2;
  n1 = stormAllocLike(c);
  n2 = stormAllocLike(c);
  double V, V_hat;
#endif

  //stormFill(mesh, c, 1.0);
  stormSpFuncProd(mesh, c, c, Initial_Data, STORM_NULL);
  //SR_SFuncProd(mesh, v, v, Initial_Data, STORM_NULL);
  //stormFillRandom(mesh, c, -1.0, +1.0);
  stormFill(mesh, v, 0.0);
  stormFill(mesh, p, 0.0);

  double total_time = 0.0;

  for (int time = 0; time <= 0*20+1*200000; ++time) {

    for (int frac = 0; time != 0 && frac < 10; ++frac) {

      struct timespec start, finish;
      clock_gettime(CLOCK_MONOTONIC, &start);

      CahnHilliard_Step(mesh, c, v, c_hat, w_hat
#if YURI
        , V, V_hat
#endif
        );
      NavierStokes_VaD_Step(mesh, p, v, c_hat, w_hat, p_hat, v_hat, d, rho
#if YURI
        , n1, n2, V, V_hat
#endif
        );

      clock_gettime(CLOCK_MONOTONIC, &finish);
      double elapsed = (finish.tv_sec - start.tv_sec);
      elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
      total_time += elapsed;

#if YURI
      fprintf(volFile, "%f\n", V_hat), fflush(volFile);
#endif

      stormSwap(c, c_hat);
      stormSwap(p, p_hat);
      stormSwap(v, v_hat);
    }

    char filename[256];
    printf("time = %f\n", total_time);
    sprintf(filename, "out/fld-%d.vti", time);
    stormIOList_t io = SR_IO_Begin();
    SR_IO_Add(io, v, "velocity");
    SR_IO_Add(io, c, "phase");
    SR_IO_Add(io, p, "pressure");
    SR_IO_Add(io, rho, "density");
    SR_IO_Add(io, d, "divV");
#if YURI
    SR_IO_Add(io, n1, "n1");
    SR_IO_Add(io, n2, "n2");
#endif
    SR_IO_Flush(io, mesh, filename);
  }

#if YURI
  fclose(volFile);
#endif
  return 0;
}
