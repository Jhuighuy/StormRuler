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

#define SR_MATLAB 0
#include "StormRuler_API.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define min(x, y) ( (x) < (y) ? (x) : (y) )
#define max(x, y) ( (x) > (y) ? (x) : (y) )

static double tau = 1.0e-2, Gamma = 1.0e-4, sigma = 1.0;

static void SetBCs_c(SR_tMesh mesh, SR_tFieldR c) {
  SR_ApplyBCs(mesh, c, SR_ALL, SR_PURE_NEUMANN);
  //SR_ApplyBCs(mesh, c, 4, SR_DIRICHLET(1.0));
} // SetBCs_c

static void SetBCs_w(SR_tMesh mesh, SR_tFieldR w) {
  SR_ApplyBCs(mesh, w, SR_ALL, SR_PURE_NEUMANN);
} // SetBCs_w

static void SetBCs_p(SR_tMesh mesh, SR_tFieldR p) {
  SR_ApplyBCs(mesh, p, SR_ALL, SR_PURE_NEUMANN);
  ////SR_ApplyBCs(mesh, p, 2, SR_DIRICHLET(1.0));
  ////SR_ApplyBCs(mesh, p, 4, SR_DIRICHLET(50.0));
} // SetBCs_p

static void SetBCs_v(SR_tMesh mesh, SR_tFieldR v) {
  SR_ApplyBCs(mesh, v, SR_ALL, SR_PURE_DIRICHLET);
  SR_ApplyBCs_SlipWall(mesh, v, 3);
  SR_ApplyBCs_InOutLet(mesh, v, 2);
  SR_ApplyBCs_InOutLet(mesh, v, 4);
} // SetBCs_v

void dWdC(int size, SR_REAL* Wc, const SR_REAL* c, void* env) {
  const SR_REAL x = *c;

  *Wc = (x < -1.0) ? 2.0*(1.0+x) : ( (x > 1.0) ? (2.0*(x-1.0)) : x*(x*x - 1.0) );
} // dWdC

void Vol(int size, SR_REAL* Ic, const SR_REAL* c, void* env) {
  SR_REAL x = *c;

  x = min(+1.0, max(-1.0, x));
  x = 0.5*(x + 1.0);
  *Ic = round(x);
}

SR_tFieldR vvv;

static void CahnHilliard_MatVec(SR_tMesh mesh,
    SR_tFieldR Qc, SR_tFieldR c, void* env) {

  SR_tFieldR v = vvv;
  SetBCs_c(mesh, c);
  SetBCs_v(mesh, v);

  SR_tFieldR tmp = SR_Alloc_Mold(c);
  SR_Scale(mesh, tmp, c, 2.0*sigma);
  SR_DivGrad(mesh, tmp, -Gamma, c);

  SetBCs_w(mesh, tmp);
  
  SR_Set(mesh, Qc, c);
  SR_Conv(mesh, Qc, -tau, c, vvv);
  SR_DivGrad(mesh, Qc, -tau, tmp);

  SR_Free(tmp);
} // CahnHilliard_MatVec

static SR_REAL CahnHilliard_Step(SR_tMesh mesh,
    SR_tFieldR c, SR_tFieldR v,
    SR_tFieldR c_hat, SR_tFieldR w_hat) {

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

  SR_tFieldR rhs = SR_Alloc_Mold(c);
  SR_Set(mesh, rhs, c);
  SR_FuncProd(mesh, w_hat, c, dWdC, NULL);
  SR_Sub(mesh, w_hat, w_hat, c, 2.0*sigma, 1.0);
  SetBCs_w(mesh, w_hat);
  SR_DivGrad(mesh, rhs, tau, w_hat);

  SR_Set(mesh, c_hat, c);
  SR_LinSolve(mesh, "BiCGStab", "", c_hat, rhs, CahnHilliard_MatVec, vvv=v, NULL, NULL);
  SR_Free(rhs);

  SetBCs_c(mesh, c_hat);

  SR_FuncProd(mesh, w_hat, c_hat, dWdC, NULL);
  SR_DivGrad(mesh, w_hat, -Gamma, c_hat);

  return SR_Integrate(mesh, c_hat, Vol, NULL);
} // CahnHilliard_Step

static double beta = 0.1;
static double rho_1 = 1.0, rho_2 = 300.0, mu_1 = 0.1, mu_2 = 0.1;

void InvRho(int size, SR_REAL* inv_rho, const SR_REAL* rho, void* env) {
  *inv_rho = 1.0/(*rho);
} // InvRho

void Fix_Mu(int dim, const SR_REAL* r,
    int size, SR_REAL* fix_mu, const SR_REAL* mu, void* env) {
  *fix_mu = *mu;

  static const SR_REAL L = 2.0*M_PI;
  if (fabs(2*L-r[1]) <= L*0.300 || fabs(r[1] < L*0.5)) {
    *fix_mu = 0.1;
  }

} // Fix_Mu

static SR_tFieldR rho_inv_, mu_;

static void Poisson_VaD_MatVec(SR_tMesh mesh,
    SR_tFieldR Lp, SR_tFieldR p, void* env) {
  
  SR_tFieldR rho_inv = rho_inv_;

  SetBCs_p(mesh, p);

  SR_Set(mesh, Lp, p);
  SR_DivKGrad(mesh, Lp, -tau, rho_inv, p);
  
} // Poisson_MatVec

static void NavierStokes_VaD_MatVec(SR_tMesh mesh,
    SR_tFieldR Av, SR_tFieldR v, void* env) {

  SetBCs_v(mesh, v);

  SR_Set(mesh, Av, v);
  SR_Conv(mesh, Av, -tau, v, v);

  SR_tFieldR tmp = SR_Alloc_Mold(v);
  SR_Fill(mesh, tmp, 0.0, 0.0);
  SR_DivGrad(mesh, tmp, tau, v);
  SR_Mul(mesh, tmp, mu_, tmp);
  SR_Mul(mesh, tmp, rho_inv_, tmp);
  SR_Sub(mesh, Av, Av, tmp, 1.0, 1.0);
  SR_Free(tmp);

} // NavierStokes_VaD_MatVec

static void NavierStokes_VaD_Step(SR_tMesh mesh,
  SR_tFieldR p, SR_tFieldR v,
  SR_tFieldR c, SR_tFieldR w,
  SR_tFieldR p_hat, SR_tFieldR v_hat, SR_tFieldR d) {

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
  SR_tFieldR rho = SR_Alloc_Mold(c);
  SR_Fill(mesh, rho, 0.5*(rho_1 + rho_2), 0.0);
  SR_Add(mesh, rho, rho, c, 0.5*(rho_2 - rho_1), 1.0);

  SR_tFieldR rho_inv = SR_Alloc_Mold(rho);
  SR_FuncProd(mesh, rho_inv, rho, InvRho, NULL);

  SR_tFieldR mu = SR_Alloc_Mold(c);
  SR_Fill(mesh, mu, 0.5*(mu_1 + mu_2), 0.0);
  SR_Add(mesh, mu, mu, c, 0.5*(mu_2 - mu_1), 1.0);


  //
  // Compute 𝒗̂ prediction.
  //
  SetBCs_w(mesh, w);
  SetBCs_v(mesh, v);

  SR_tFieldR rhs = SR_Alloc_Mold(p);
  SR_Fill(mesh, rhs, 0.0, 0.0);
  SR_Grad(mesh, rhs, tau, w);
  SR_Mul(mesh, rhs, c, rhs);
  SR_Mul(mesh, rhs, rho_inv, rhs);

  SR_Add(mesh, rhs, rhs, v, 1.0, 1.0);
  
  SR_Set(mesh, v_hat, v);
  rho_inv_ = rho_inv, mu_ = mu;
  SR_Solve_JFNK(mesh, v_hat, v, NavierStokes_VaD_MatVec, NULL);

  //
  // Solve pressure equation and correct 𝒗̂.
  // 
  SR_Set(mesh, rhs, p);
  SR_Div(mesh, rhs, 1.0, v_hat);

  SR_Set(mesh, p_hat, p);
  SR_LinSolve(mesh, "CG", "", p_hat, rhs, 
    Poisson_VaD_MatVec, rho_inv_=rho_inv, NULL, NULL);
  SR_Free(rhs);

  SetBCs_p(mesh, p_hat);

  SR_tFieldR tmp = SR_Alloc_Mold(v);
  SR_Fill(mesh, tmp, 0.0, 0.0);
  SR_Grad(mesh, tmp, tau, p_hat);
  SR_Mul(mesh, tmp, rho_inv, tmp);
  SR_Add(mesh, v_hat, v_hat, tmp, 1.0, 1.0);
  SR_Free(tmp);

  if (d != NULL) {
    SetBCs_v(mesh, v_hat);
    SR_Fill(mesh, d, 0.0, 0.0);
    SR_Div(mesh, d, -1.0, v);
  }

  SR_Free(rho);
  SR_Free(rho_inv);
  SR_Free(mu);

} // NavierStokes_VaD_Step

void Initial_Data(int dim, const SR_REAL* r,
    int size, SR_REAL* c, const SR_REAL* _, void* env) {

  static const SR_REAL L = 1.0;
  int in = 0;
  if (fabs(r[0]-0*L) <= L*0.101 && 
      fabs(2*L-r[1]) <= L*0.665) {
    in = 1.0;
  }

  *c = -1.0;
  if (in) {
    *c = 1.0;
  }
} // Initial_Data

void main() {

  FILE* volFile = fopen("vol.txt", "w");

  SR_tMesh mesh = SR_InitMesh();

  SR_tFieldR c, p, v, c_hat, w_hat, p_hat, v_hat, d;
  c = SR_AllocR(mesh, 1, 0);
  p = SR_AllocR(mesh, 1, 0);
  v = SR_AllocR(mesh, 1, 1);
  d = SR_AllocR(mesh, 1, 0);
  c_hat = SR_Alloc_Mold(c);
  w_hat = SR_Alloc_Mold(c);
  p_hat = SR_Alloc_Mold(p);
  v_hat = SR_Alloc_Mold(v);

  //SR_Fill(mesh, c, 1.0, 0.0);
  SR_SFuncProd(mesh, c, c, Initial_Data, NULL);
  //SR_SFuncProd(mesh, v, v, Initial_Data, NULL);
  //SR_Fill_Random(mesh, c, -1.0, +1.0);
  SR_Fill(mesh, v, 0.0, 0.0);
  SR_Fill(mesh, p, 0.0, 0.0);

  for (int time = 0; time <= 200000; ++time) {

    for (int frac = 0; time != 0 && frac < 1; ++ frac) {

      SR_REAL vol = CahnHilliard_Step(mesh, c, v, c_hat, w_hat);
      NavierStokes_VaD_Step(mesh, p, v, c_hat, w_hat, p_hat, v_hat, d);
      fprintf(volFile, "%f\n", vol), fflush(volFile);

      SR_Swap(&c, &c_hat);
      SR_Swap(&p, &p_hat);
      SR_Swap(&v, &v_hat);
    }

    char filename[256];
    sprintf(filename, "out/fld-%d.vtk", time);
    SR_tIOList io = SR_IO_Begin();
    SR_IO_Add(io, v, "velocity");
    SR_IO_Add(io, c, "phase");
    SR_IO_Add(io, p, "pressure");
    SR_IO_Add(io, c, "density");
    SR_IO_Add(io, d, "divV");
    SR_IO_Flush(io, mesh, filename);
  }

  fclose(volFile);
  exit(0);
}
