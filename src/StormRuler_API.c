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

#include "StormRuler_API.h"

#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double tau = (M_PI/50)*(M_PI/50), Gamma = 0.01, sigma = 1.0;

static void Advection_MatVec(SR_tMesh mesh,
    SR_tFieldR Ac, SR_tFieldR c, void* env) {
  SR_tFieldR v = (SR_tFieldR)env;

  SR_ApplyBCs(mesh, v, SR_ALL, SR_PURE_DIRICHLET);
  SR_ApplyBCs(mesh, c, SR_ALL, SR_PURE_NEUMANN);
  SR_Set(mesh, Ac, c);
  SR_Conv(mesh, Ac, tau, c, v);
} // Advection_MatVec

static void Advection_Step(SR_tMesh mesh,
    SR_tFieldR c, SR_tFieldR v,
    SR_tFieldR c_hat) {

  // 
  // Compute a single time step of the
  // advection equation.
  //
  // ∂𝑐/∂𝑡 + ∇⋅𝑐𝒗 = 0.
  //
  // with the implicit scheme:
  // 
  // 𝑐̂ + 𝜏∇⋅𝑐̂𝒗 = 𝑐,
  //

  SR_LinSolve(mesh, c_hat, c, Advection_MatVec, v, 
    SR_BiCGStab, SR_Precond_None, NULL, NULL);
} // Advection_Step

void dWdC(int size, SR_REAL* Wc, const SR_REAL* c, void* env) {
  //*Wc = 0.5*(1.0 - (*c)*(*c));
  //*Wc = (*Wc)*(*Wc);
  *Wc = -(*c)*(1.0 - (*c)*(*c));
} // dWdC

static void CahnHilliard_MatVec(SR_tMesh mesh,
    SR_tFieldR Qc, SR_tFieldR c, void* env) {
  
  SR_tFieldR tmp = SR_Alloc_Mold(c);
  SR_ApplyBCs(mesh, c, SR_ALL, SR_PURE_NEUMANN);
  SR_Fill(mesh, tmp, 0.0, 0.0);
  SR_DivGrad(mesh, tmp, 1.0, c);
  
  SR_Scale(mesh, Qc, c, 1.0 - 2.0*tau*sigma);
  SR_ApplyBCs(mesh, tmp, SR_ALL, SR_PURE_NEUMANN);
  SR_DivGrad(mesh, Qc, tau*Gamma, tmp);

  SR_Free(tmp);
} // CahnHilliard_MatVec

static void CahnHilliard_Step(SR_tMesh mesh,
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
  // 𝑐̃ + 𝜏∇⋅𝑐̃𝒗 = 𝑐,
  // 𝓠𝑐̂ = (1 - 2𝜏𝜎)𝑐̂ + 𝜏𝛾Δ²𝑐̂ = (1 - 2𝜏𝜎)𝑐̃ + 𝜏∇⋅𝑐̃𝒗 + 𝜏Δ𝑊'(𝑐̃),
  // 𝑤̂ = 𝑑𝑊(𝑐̂)/𝑑𝑐 - 𝛾Δ𝑐̂.
  //

  SR_ApplyBCs(mesh, c, SR_ALL, SR_PURE_NEUMANN);
  SR_ApplyBCs(mesh, v, SR_ALL, SR_PURE_DIRICHLET);

  SR_FuncProdR(mesh, w_hat, c, dWdC, NULL);
  SR_ApplyBCs(mesh, w_hat, SR_ALL, SR_PURE_NEUMANN);

  SR_tFieldR RHS = SR_Alloc_Mold(c);
  SR_Scale(mesh, RHS, c, 1.0 - 2.0*tau*sigma);
  SR_Conv(mesh, RHS, tau, c, v);
  SR_DivGrad(mesh, RHS, tau, w_hat);

  SR_Set(mesh, c_hat, c);
  SR_LinSolve(mesh, c_hat, RHS, CahnHilliard_MatVec, NULL,
    SR_CG, SR_Precond_None, NULL, NULL);

  SR_Free(RHS);

  SR_ApplyBCs(mesh, c_hat, SR_ALL, SR_PURE_NEUMANN);

  SR_FuncProdR(mesh, w_hat, c_hat, dWdC, NULL);
  SR_DivGrad(mesh, w_hat, -Gamma, c_hat);
} // CahnHilliard_Step

static double rho = 1.0, mu = 0.1;

static void Poisson_MatVec(SR_tMesh mesh,
    SR_tFieldR Lp, SR_tFieldR p, void* env) {
  
  SR_ApplyBCs(mesh, p, SR_ALL, SR_PURE_NEUMANN);

  SR_Fill(mesh, Lp, 0.0, 0.0);
  SR_DivGrad(mesh, Lp, 1.0, p);
  
} // Poisson_MatVec

static void NavierStokes_Step(SR_tMesh mesh,
  SR_tFieldR p, SR_tFieldR v,
  SR_tFieldR c, SR_tFieldR w,
  SR_tFieldR p_hat, SR_tFieldR v_hat) {

  // 
  // Compute a single time step of the incompressible
  // Navier-Stokes equation with convection term:
  //
  // 𝜌(∂𝒗/∂𝑡 + 𝒗(∇⋅𝒗)) + ∇𝑝 = 𝜇Δ𝒗 + 𝙛,
  // ∇⋅𝒗 = 0, 𝙛 = 𝑐∇𝑤,
  //
  // with the semi-implicit scheme:
  // 
  // 𝒗̂ ← 𝒗 - 𝜏𝒗(∇⋅𝒗) + (𝜏𝜇/𝜌)Δ𝒗 + (𝜏/𝜌)𝙛,
  // Δ𝑝̂ = 𝜌/𝜏∇⋅𝒗,
  // 𝒗̂ ← 𝒗̂ - (𝜏/𝜌)∇𝑝̂.
  // 

  //
  // Compute 𝒗̂ prediction.
  //
  SR_ApplyBCs(mesh, w, SR_ALL, SR_PURE_NEUMANN);
  SR_ApplyBCs(mesh, p, SR_ALL, SR_PURE_NEUMANN);
  SR_ApplyBCs(mesh, v, SR_ALL, SR_PURE_DIRICHLET);

  SR_Set(mesh, v_hat, v);
  SR_Conv(mesh, v_hat, tau, v, v);
  SR_DivGrad(mesh, v_hat, tau*mu/rho, v);

  //
  // Compute 𝙛 = -𝑐∇𝑤, 𝒗̂ ← 𝒗̂ + 𝙛.
  //
  SR_tFieldR f = SR_Alloc_Mold(v);

  SR_Fill(mesh, f, 0.0, 0.0);
  SR_Grad(mesh, f, 1.0, w);
  SR_Mul(mesh, f, c, f);

  SR_Add(mesh, v_hat, v_hat, f, tau/rho, 1.0);

  SR_Free(f);

  SR_ApplyBCs(mesh, v_hat, SR_ALL, SR_PURE_DIRICHLET);

  //
  // Solve pressure equation and correct 𝒗̂.
  // 
  SR_tFieldR RHS = SR_Alloc_Mold(p);
  SR_Fill(mesh, RHS, 0.0, 0.0);
  SR_Div(mesh, RHS, -rho/tau, v_hat);

  SR_Set(mesh, p_hat, p);
  SR_LinSolve(mesh, p_hat, RHS, Poisson_MatVec, NULL,
    SR_MINRES, SR_Precond_None, NULL, NULL);

  SR_ApplyBCs(mesh, p_hat, SR_ALL, SR_PURE_NEUMANN);
  SR_Grad(mesh, v_hat, tau/rho, p_hat);

  SR_Free(RHS);

} // NavierStokes_Step

static double rho_1 = 1.0, rho_2 = 2.0, mu_1 = 0.1, mu_2 = 0.2;

static void NavierStokes_VaryingDensity_Step(SR_tMesh mesh,
  SR_tFieldR p, SR_tFieldR v,
  SR_tFieldR c, SR_tFieldR w,
  SR_tFieldR rho_hat, SR_tFieldR mu_hat, 
  SR_tFieldR p_hat, SR_tFieldR v_hat) {

  // 
  // Compute a single time step of the incompressible
  // Navier-Stokes equation with convection term:
  //
  // 𝜌(∂𝒗/∂𝑡 + 𝒗(∇⋅𝒗)) + ∇𝑝 = 𝜇Δ𝒗 + 𝙛,
  // ∇⋅𝒗 = 0, 𝙛 = 𝑐∇𝑤,
  // 𝜌 = ½𝜌₁(1 - 𝑐) + ½𝜌₂(1 + 𝑐),
  // 𝜇 = ½𝜇₁(1 - 𝑐) + ½𝜇₂(1 + 𝑐),
  //
  // with the semi-implicit scheme:
  // 
  // 𝜌̂ ← ½(𝜌₁ + 𝜌₂) - ½(𝜌₁ - 𝜌₂)𝑐,
  // 𝜇̂ ← ½(𝜇₁ + 𝜇₂) - ½(𝜇₁ - 𝜇₂)𝑐,
  // 𝒗̂ ← 𝒗 - 𝜏𝒗(∇⋅𝒗) + (𝜏𝜇̂Δ𝒗 + 𝜏𝙛)/𝜌̂,
  // ∇⋅(1/𝜌̂∇𝑝̂) = 1/𝜏∇⋅𝒗,
  // 𝒗̂ ← 𝒗̂ - (𝜏/𝜌)∇𝑝̂.
  // 

  //
  // Compute 𝜌̂, 𝜇̂.
  //
  SR_Fill(mesh, rho_hat, 0.5*(rho_1 + rho_2), 0.0);
  SR_Sub(mesh, rho_hat, rho_hat, c, 0.5*(rho_1 - rho_2), 1.0);
  SR_Fill(mesh, mu_hat, 0.5*(mu_1 + mu_2), 0.0);
  SR_Sub(mesh, mu_hat, mu_hat, c, 0.5*(rho_1 - rho_2), 1.0);

  //
  // Compute 𝒗̂ prediction.
  //
  SR_ApplyBCs(mesh, w, SR_ALL, SR_PURE_NEUMANN);
  SR_ApplyBCs(mesh, p, SR_ALL, SR_PURE_NEUMANN);
  SR_ApplyBCs(mesh, v, SR_ALL, SR_PURE_DIRICHLET);

  SR_Set(mesh, v_hat, v);
  SR_Conv(mesh, v_hat, tau, v, v);

  SR_tFieldR TMP = SR_Alloc_Mold(v);

  SR_Fill(mesh, TMP, 0.0, 0.0);
  SR_DivGrad(mesh, TMP, 1.0, v);
  SR_Mul(mesh, TMP, mu_hat, TMP);
  //SR_Mul(mesh, TMP, rho_hat, TMP, -1);
  SR_Add(mesh, v_hat, v_hat, TMP, tau, 1.0);

  SR_Free(TMP);

} // NavierStokes_VaryingDensity_Step

void pure_c_main() {

  lua_State *L = luaL_newstate();
  luaL_openlibs(L);
  lua_close(L);

  SR_tMesh mesh = SR_InitMesh();

  SR_tFieldR c, p, v, c_hat, w_hat, p_hat, v_hat;
  c = SR_AllocR(mesh, 1, 0);
  p = SR_AllocR(mesh, 1, 0);
  v = SR_AllocR(mesh, 1, 1);
  c_hat = SR_Alloc_Mold(c);
  w_hat = SR_Alloc_Mold(c);
  p_hat = SR_Alloc_Mold(p);
  v_hat = SR_Alloc_Mold(v);

  SR_Fill_Random(mesh, c, -1.0, +1.0);
  SR_Fill(mesh, v, 0.0, 0.0);

  SR_tIOList io = SR_IO_Begin();
  SR_IO_Add(io, v, "velocity");
  SR_IO_Add(io, c, "phase");
  SR_IO_Add(io, p, "pressure");
  SR_IO_Add(io, c, "density");
  SR_IO_Flush(io, mesh, "out/fld-0.vtk");

  for (int time = 0; time <= 200; ++time) {

    if (time == 0) continue;

    for (int frac = 0; frac < 10; ++ frac) {

      CahnHilliard_Step(mesh, c, v, c_hat, w_hat);
      NavierStokes_Step(mesh, p, v, c_hat, w_hat, p_hat, v_hat);

      SR_Swap(&c, &c_hat);
      SR_Swap(&p, &p_hat);
      SR_Swap(&v, &v_hat);
    }

    char filename[256];
    sprintf(filename, "out/fld-%d.vtk", time);
    io = SR_IO_Begin();
    SR_IO_Add(io, v, "velocity");
    SR_IO_Add(io, c, "phase");
    SR_IO_Add(io, p, "pressure");
    SR_IO_Add(io, c, "density");
    SR_IO_Flush(io, mesh, filename);
  }

  exit(0);
}
