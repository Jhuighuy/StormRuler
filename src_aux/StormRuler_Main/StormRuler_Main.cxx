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
#include <stormSolvers/stormSolverFactory.hxx>
#include <stormSolvers/stormPreconditionerFactory.hxx>

#include <cstring>

template<typename stormMatVecFuncT_t>
STORM_INL void stormLinSolve2(stormMesh_t mesh,
                              stormString_t method,
                              stormString_t preMethod,
                              stormArray_t x,
                              stormArray_t b,
                              stormMatVecFuncT_t matVec) {
  stormArray xx = {mesh, x}, bb = {mesh, b};
  auto symOp = stormMakeSymmetricOperator<stormArray>(
    [&](stormArray& yy, const stormArray& xx) {
      matVec(yy.Mesh, yy.Array, xx.Array);
    });

  auto solver = stormMakeIterativeSolver<stormArray>(method);
  solver->PreOp = stormMakePreconditioner<stormArray>(preMethod);
  solver->Solve(xx, bb, *symOp);

} // stormLinSolve2<...>

template<typename stormMatVecFuncT_t>
STORM_INL void stormNonlinSolve2(stormMesh_t mesh,
                                 stormString_t method,
                                 stormArray_t x,
                                 stormArray_t b,
                                 stormMatVecFuncT_t matVec) {
  stormArray xx = {mesh, x}, bb = {mesh, b};
  std::unique_ptr<stormOperator<stormArray>> op =
    stormMakeSymmetricOperator<stormArray>(
      [&](stormArray& yy, const stormArray& xx) {
        matVec(yy.Mesh, yy.Array, xx.Array);
      });

  std::unique_ptr<stormIterativeSolver<stormArray>> solver = 
    std::make_unique<stormJfnkSolver<stormArray>>();
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

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define YURI 1

#define min(x, y) ( (x) < (y) ? (x) : (y) )
#define max(x, y) ( (x) > (y) ? (x) : (y) )

static double tau = 1.0e-2, Gamma = 1.0e-4, sigma = 1.0;

static void SetBCs_c(stormMesh_t mesh, stormArray_t c) {
  stormApplyBCs(mesh, c, SR_ALL, SR_PURE_NEUMANN);
  //stormApplyBCs(mesh, c, 4, SR_DIRICHLET(1.0));
} // SetBCs_c

static void SetBCs_w(stormMesh_t mesh, stormArray_t w) {
  stormApplyBCs(mesh, w, SR_ALL, SR_PURE_NEUMANN);
} // SetBCs_w

static void SetBCs_p(stormMesh_t mesh, stormArray_t p) {
  stormApplyBCs(mesh, p, SR_ALL, SR_PURE_NEUMANN);
  //stormApplyBCs(mesh, p, 2, SR_DIRICHLET(1.0));
  //stormApplyBCs(mesh, p, 4, SR_DIRICHLET(50.0));
} // SetBCs_p

static void SetBCs_v(stormMesh_t mesh, stormArray_t v) {
  stormApplyBCs(mesh, v, SR_ALL, SR_PURE_DIRICHLET);
  stormApplyBCs_SlipWall(mesh, v, 3);
  stormApplyBCs_InOutLet(mesh, v, 2);
  stormApplyBCs_InOutLet(mesh, v, 4);
} // SetBCs_v

stormSize_t slice = 0;
stormVector<double> phi_set;
stormMatrix<double> D1_W_vs_phi_sandwitch;
stormTensor3R<double> nPart_vs_phi_sandwitch;
static double const mol_mass[2] = {0.0440098, 0.1422853};

template<class Tensor>
void ReadTensor(Tensor& tensor, std::string path) {
  path = "/home/jhuighuy/Downloads/decane_droplet/" + path;
  std::ifstream ifstream(path);
  std::for_each_n(tensor.Data(), tensor.Size(), [&](auto& value) {
    ifstream >> value;
  });
}

void dWdC(stormSize_t size, stormReal_t* Wc, const stormReal_t* c, void* env) {
  const stormReal_t x = *c;

#if !YURI

  // [-1,+1] CH.
  //*Wc = (x < -1.0) ? 2.0*(1.0+x) : ( (x > 1.0) ? (2.0*(x-1.0)) : x*(x*x - 1.0) );
  // [0,1] CH.
  *Wc = (x < -1.0) ? 2.0*x : ( (x > 1.0) ? (2.0*(x-1.0)) : 2.0*x*(x - 1.0)*(2.0*x - 1.0) );

#else

  auto const D1_W_vs_phi = [&](auto i) { return D1_W_vs_phi_sandwitch(slice, i); };

  // [0,1] MCH.
  const double h = 1.0/(phi_set.Size() - 1.0);
  if (x < 0.0) {
    *Wc = ( D1_W_vs_phi(0) + x*(D1_W_vs_phi(1) - D1_W_vs_phi(0))/h )/32.0;
    return;
  }
  if (x > 1.0) {
    *Wc = ( D1_W_vs_phi(100) + (x - 1.0)*(D1_W_vs_phi(100) - D1_W_vs_phi(99))/h )/32.0;
    return;
  }
  const int il = floor(x/h);
  const int ir = ceil(x/h);
  if (il == ir) {
    *Wc = ( D1_W_vs_phi(il) )/32.0;  
  }
  *Wc = ( D1_W_vs_phi(il) + (x - h*il)*(D1_W_vs_phi(ir) - D1_W_vs_phi(il))/h )/32.0;
  return;

#endif
} // dWdC

void Vol(stormSize_t size, stormReal_t* Ic, const stormReal_t* c, void* env) {
  stormReal_t x = *c;

  //x = min(+1.0, max(-1.0, x));
  //x = 0.5*(x + 1.0);
  //*Ic = round(x);
  *Ic = 1.0 - x;
}

static stormReal_t CahnHilliard_Step(stormMesh_t mesh,
    stormArray_t c, stormArray_t v,
    stormArray_t c_hat, stormArray_t w_hat) {

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
  stormLinSolve2(mesh, STORM_KSP_GMRES, STORM_PRE_NONE/*"extr"*/, c_hat, rhs,
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

  return stormIntegrate(mesh, c_hat, Vol, STORM_NULL);
} // CahnHilliard_Step

static double mu_1 = 0.08, mu_2 = 0.08;
#if !YURI
static double rho_1 = 1.0, rho_2 = 1.0;
#endif

void InvRho(stormSize_t size, stormReal_t* inv_rho, const stormReal_t* rho, void* env) {
  *inv_rho = 1.0/(*rho);
} // InvRho

static int II;

void NVsC(stormSize_t size, stormReal_t* n, const stormReal_t* c, void* env) {
  stormReal_t x = *c;

  auto const nPart_vs_phi = [&](auto i, auto j) { return nPart_vs_phi_sandwitch(slice, i, j); };

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

static void NavierStokes_VaD_Step(stormMesh_t mesh,
  stormArray_t p, stormArray_t v,
  stormArray_t c, stormArray_t w,
  stormArray_t p_hat, stormArray_t v_hat, stormArray_t d, stormArray_t rho
#if YURI
  , stormArray_t n1, stormArray_t n2
#endif
) {

  // 
  // Compute a single time step of the incompressible
  // Navier-Stokes equation:
  //
  // 𝜌(∂𝒗/∂𝑡 + 𝒗(∇⋅𝒗)) + ∇𝑝 = 𝜇Δ𝒗 + 𝙛,
  // 𝜌 = ½𝜌₁(1 - 𝑐) + ½𝜌₂(1 + 𝑐),
  // 𝜇 = ½𝜇₁(1 - 𝑐) + ½𝜇₂(1 + 𝑐),
  // ∇⋅𝒗 = 0, 𝙛 = -𝑐∇𝑤,
  //
  // with the semi-implicit scheme:
  // 
  // 𝜌 ← ½(𝜌₁ + 𝜌₂) + (½𝜌₂ - ½𝜌₁)𝑐,
  // 𝜇 ← ½(𝜇₁ + 𝜇₂) + (½𝜇₂ - ½𝜇₁)𝑐,
  // 𝒗̂ + 𝜏𝒗̂(∇⋅𝒗̂) - (𝜏𝜇/𝜌)Δ𝒗̂ = 𝒗 + (𝜏/𝜌)𝙛,
  // 𝑝̂ - 𝜏∇⋅(∇𝑝̂/𝜌) = 𝑝 - ∇⋅𝒗,
  // 𝒗̂ ← 𝒗̂ - (𝜏/𝜌)∇𝑝̂.
  //

  //
  // Compute 𝜌, 𝜇, 1/𝜌.
  //
#if YURI
  II = 0; stormFuncProd(mesh, n1, c, NVsC, STORM_NULL);
  II = 1; stormFuncProd(mesh, n2, c, NVsC, STORM_NULL);
  stormAdd(mesh, rho, n1, n2, mol_mass[1], mol_mass[0]);
#else
  stormFill(mesh, rho, 0.5*(rho_1 + rho_2));
  stormAdd(mesh, rho, rho, c, 0.5*(rho_2 - rho_1));
#endif

  stormArray_t rho_inv = stormAllocLike(rho);
  stormFuncProd(mesh, rho_inv, rho, InvRho, STORM_NULL);

  stormArray_t mu = stormAllocLike(c);
  stormFill(mesh, mu, 0.5*(mu_1 + mu_2));
  stormAdd(mesh, mu, mu, c, 0.5*(mu_2 - mu_1));

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
  stormNonlinSolve2(mesh, STORM_JFNK, v_hat, v, 
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
  stormLinSolve2(mesh, STORM_KSP_CG, STORM_PRE_NONE/*"extr"*/, p_hat, rhs,
    [&](stormMesh_t mesh, stormArray_t Lp, stormArray_t p) {
      SetBCs_p(mesh, p);

      stormSet(mesh, Lp, p);
      stormDivWGrad(mesh, Lp, -tau, rho_inv, p);
    });

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
  int in = 0;
  if (fabs(r[0]-0*L) <= L*0.101 && 
      fabs(2*L-r[1]) <= L*0.665) {
    in = 1.0;
  }

  *c = 0.0;
  if (in) {
    *c = 1.0;
  }
} // Initial_Data

int main() {

  phi_set.Assign(101);
  ReadTensor(phi_set, "phi.csv");
  D1_W_vs_phi_sandwitch.Assign(1000, 101);
  ReadTensor(D1_W_vs_phi_sandwitch, "D1_W_vs_phi.csv");
  nPart_vs_phi_sandwitch.Assign(1000, 101, 2);
  ReadTensor(nPart_vs_phi_sandwitch, "n_part_vs_phi_sandwich.csv");

  //FILE* dWdCF = fopen("dWdC.txt", "w");
  //for (double x = -0.1; x <= 1.1; x += 0.0001) {
  //  double y;
  //  dWdC(1, &y, &x, STORM_NULL);
  //  fprintf(dWdCF, "%f %f\n", x, y);
  //}
  //fclose(dWdCF);
  //return 1;

  FILE* volFile = fopen("vol.txt", "w");

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

      stormReal_t vol = CahnHilliard_Step(mesh, c, v, c_hat, w_hat);
      NavierStokes_VaD_Step(mesh, p, v, c_hat, w_hat, p_hat, v_hat, d, rho
#if YURI
        , n1, n2
#endif
      );

      clock_gettime(CLOCK_MONOTONIC, &finish);
      double elapsed = (finish.tv_sec - start.tv_sec);
      elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
      total_time += elapsed;

      fprintf(volFile, "%f\n", vol), fflush(volFile);

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

  fclose(volFile);
  return 0;
}
