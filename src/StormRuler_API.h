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

#pragma once

#include "StormRuler_Params.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

SR_OPAQUE_STRUCT(SR_SYMBOL);

SR_OPAQUE_STRUCT(SR_tMesh);

SR_OPAQUE_STRUCT(SR_tIOList);

SR_OPAQUE_STRUCT(SR_tFieldR);
SR_OPAQUE_STRUCT(SR_tFieldC);
SR_OPAQUE_STRUCT(SR_tFieldS);

typedef union {
  void* P;
  SR_tFieldR R;
  SR_tFieldC C;
  SR_tFieldS S;
} SR_tFieldA;

#define SR_NULL_A ((SR_tFieldA){NULL})

typedef enum {
  SR_Done,
  SR_Request_MatVec,
  SR_Request_MatVec_H,
} SR_eRequest;

SR_OPAQUE_STRUCT(SR_tRequestEnv);

#define SR_FIELD_GENERIC(x, func) \
  _Generic((x), SR_tFieldR: func##R, SR_tFieldC: func##C)
#define SR_FIELD_GENERIC_EXT(x, func) \
  _Generic((x), SR_tFieldR: func##R, SR_tFieldC: func##C, SR_tFieldS: func##S)

typedef struct SR_ANY {} SR_ANY;

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

SR_API SR_tMesh SR_InitMesh(void);

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

/// @{
SR_API SR_tFieldR SR_AllocR(
    SR_tMesh mesh, SR_INTEGER numVars, SR_INTEGER rank);
SR_API SR_tFieldC SR_AllocC(
    SR_tMesh mesh, SR_INTEGER numVars, SR_INTEGER rank);
SR_API SR_tFieldS SR_AllocS(
    SR_tMesh mesh, SR_INTEGER numVars, SR_INTEGER rank);
/// @}

/// @{
SR_API SR_tFieldR SR_Alloc_MoldR(SR_tFieldR mold);
SR_API SR_tFieldC SR_Alloc_MoldC(SR_tFieldC mold);
SR_API SR_tFieldS SR_Alloc_MoldS(SR_tFieldS mold);
#if SR_C11
#define SR_Alloc_Mold(mold) \
  SR_FIELD_GENERIC_EXT(mold, SR_Alloc_Mold)(mold)
#endif
/// @}

/// @{
SR_API void SR_FreeR(SR_tFieldR x);
SR_API void SR_FreeC(SR_tFieldC x);
SR_API void SR_FreeS(SR_tFieldS x);
#if SR_C11
#define SR_Free(x) \
  SR_FIELD_GENERIC_EXT(x, SR_Free)(x)
#endif
/// @}

/// @{
SR_INL void SR_SwapP(void** pX, void** pY) {
  void* z = *pX; *pX = *pY, *pY = z;
}
#if SR_C11
#define SR_Swap(pX, pY) SR_SwapP((void**)pX, (void**)pY)
#endif
/// @}

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

SR_API SR_tIOList SR_IO_Begin();

SR_API void SR_IO_Add(SR_tIOList IO, 
    SR_tFieldR x, const char* name);

SR_API void SR_IO_Flush(SR_tIOList IO, 
    SR_tMesh mesh, const char* filename);

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

/// @{
SR_API void SR_FillR(SR_tMesh mesh, 
    SR_tFieldR x, SR_REAL alpha, SR_REAL beta);
SR_API void SR_FillC(SR_tMesh mesh, 
    SR_tFieldC x, SR_REAL alpha, SR_COMPLEX beta);
SR_API void SR_FillS(SR_tMesh mesh, 
    SR_tFieldS x, SR_REAL alpha, SR_SYMBOL beta);
#if SR_C11
#define SR_Fill(mesh, x, alpha, beta) \
  SR_FIELD_GENERIC_EXT(x, SR_Fill)(mesh, x, alpha, beta)
#endif
/// @}

/// @{
SR_API void SR_Fill_RandomR(SR_tMesh mesh, 
    SR_tFieldR x, SR_REAL alpha, SR_REAL beta);
SR_API void SR_Fill_RandomC(SR_tMesh mesh, 
    SR_tFieldC x, SR_REAL alpha, SR_COMPLEX beta);
#if SR_C11
#define SR_Fill_Random(mesh, x, alpha, beta) \
  SR_FIELD_GENERIC(x, SR_Fill_Random)(mesh, x, alpha, beta)
#endif
/// @}

/// @{
SR_API void SR_SetR(SR_tMesh mesh, SR_tFieldR y, SR_tFieldR x);
SR_API void SR_SetC(SR_tMesh mesh, SR_tFieldC y, SR_tFieldC x);
SR_API void SR_SetS(SR_tMesh mesh, SR_tFieldS y, SR_tFieldS x);
#if SR_C11
#define SR_Set(mesh, y, x) \
  SR_FIELD_GENERIC_EXT(y, SR_Set)(mesh, y, x)
#endif
/// @}

/// @{
SR_API void SR_ScaleR(SR_tMesh mesh, 
    SR_tFieldR y, SR_tFieldR x, SR_REAL alpha);
SR_API void SR_ScaleC(SR_tMesh mesh, 
    SR_tFieldC y, SR_tFieldC x, SR_COMPLEX alpha);
SR_API void SR_ScaleS(SR_tMesh mesh, 
    SR_tFieldS y, SR_tFieldS x, SR_SYMBOL alpha);
#if SR_C11
#define SR_Scale(mesh, y, x, alpha) \
  SR_FIELD_GENERIC_EXT(y, SR_Scale)(mesh, y, x, alpha)
#endif
/// @}

/// @{
SR_API void SR_AddR(SR_tMesh mesh, SR_tFieldR z,
    SR_tFieldR y, SR_tFieldR x, SR_REAL alpha, SR_REAL beta);
SR_API void SR_AddC(SR_tMesh mesh, SR_tFieldC z,
    SR_tFieldC y, SR_tFieldC x, SR_COMPLEX alpha, SR_COMPLEX beta);
SR_API void SR_AddS(SR_tMesh mesh, SR_tFieldS z,
    SR_tFieldS y, SR_tFieldS x, SR_SYMBOL alpha, SR_SYMBOL beta);
#if SR_C11
#define SR_Add(mesh, z, y, x, alpha, beta) \
  SR_FIELD_GENERIC_EXT(y, SR_Add)(mesh, z, y, x, alpha, beta)
#endif
/// @}

/// @{
// No need to import from Fortran.
SR_INL void SR_SubR(SR_tMesh mesh, SR_tFieldR z,
    SR_tFieldR y, SR_tFieldR x, SR_REAL alpha, SR_REAL beta) {
  SR_AddR(mesh, z, y, x, -alpha, beta);
}
SR_INL void SR_SubC(SR_tMesh mesh, SR_tFieldC z,
    SR_tFieldC y, SR_tFieldC x, SR_COMPLEX alpha, SR_COMPLEX beta) {
  SR_AddC(mesh, z, y, x, -alpha, beta);
}
SR_INL void SR_SubS(SR_tMesh mesh, SR_tFieldS z,
    SR_tFieldS y, SR_tFieldS x, SR_SYMBOL alpha, SR_SYMBOL beta) {
  assert(0); //SR_AddS(mesh, z, y, x, -alpha, beta);
}
#if SR_C11
#define SR_Sub(mesh, z, y, x, alpha, beta) \
  SR_FIELD_GENERIC_EXT(y, SR_Sub)(mesh, z, y, x, alpha, beta)
#endif
/// @}

SR_API void SR_MulR(SR_tMesh mesh, SR_tFieldR z,
    SR_tFieldR y, SR_tFieldR x);
#define SR_Mul SR_MulR

/// @{
typedef void(*SR_tMapFuncR)(SR_INTEGER size, 
    SR_REAL* Fx, const SR_REAL* x, void* env);
typedef void(*SR_tMapFuncC)(SR_INTEGER size, 
    SR_COMPLEX* Fx, const SR_COMPLEX* x, void* env);
typedef void(*SR_tMapFuncS)(SR_INTEGER size, 
    SR_SYMBOL* Fx, const SR_SYMBOL* x, void* env);
/// @}

/// @{
SR_API void SR_FuncProdR(SR_tMesh mesh,
    SR_tFieldR y, SR_tFieldR x, SR_tMapFuncR f, void* env);
SR_API void SR_FuncProdC(SR_tMesh mesh,
    SR_tFieldC y, SR_tFieldC x, SR_tMapFuncC f, void* env);
SR_API void SR_FuncProdS(SR_tMesh mesh,
    SR_tFieldS y, SR_tFieldS x, SR_tMapFuncS f, void* env);
#if SR_C11
#define SR_FuncProd(mesh, y, x, f, env) \
  SR_FIELD_GENERIC_EXT(y, SR_FuncProd)(mesh, y, x, f, env)
#endif
/// @}

/// @{
typedef void(*SR_tSMapFuncR)(SR_INTEGER dim, const SR_REAL* r,
    SR_INTEGER size, SR_REAL* Fx, const SR_REAL* x, void* env);
typedef void(*SR_tSMapFuncC)(SR_INTEGER dim, const SR_REAL* r,
    SR_INTEGER size, SR_COMPLEX* Fx, const SR_COMPLEX* x, void* env);
typedef void(*SR_tSMapFuncS)(SR_INTEGER dim, const SR_REAL* r,
    SR_INTEGER size, SR_SYMBOL* Fx, const SR_SYMBOL* x, void* env);
/// @}

/// @{
SR_API void SR_SFuncProdR(SR_tMesh mesh,
    SR_tFieldR y, SR_tFieldR x, SR_tSMapFuncR f, void* env);
SR_API void SR_SFuncProdC(SR_tMesh mesh,
    SR_tFieldC y, SR_tFieldC x, SR_tSMapFuncC f, void* env);
SR_API void SR_SFuncProdS(SR_tMesh mesh,
    SR_tFieldS y, SR_tFieldS x, SR_tSMapFuncS f, void* env);
#if SR_C11
#define SR_SFuncProd(mesh, y, x, f, env) \
  SR_FIELD_GENERIC_EXT(y, SR_SFuncProd)(mesh, y, x, f, env)
#endif
/// @}

/// @{
typedef void(*SR_tMatVecFuncR)(SR_tMesh mesh,
    SR_tFieldR Ax, SR_tFieldR x, void* env);
typedef void(*SR_tMatVecFuncC)(SR_tMesh mesh,
    SR_tFieldC Ax, SR_tFieldC x, void* env);
typedef void(*SR_tMatVecFuncS)(SR_tMesh mesh,
    SR_tFieldS Ax, SR_tFieldS x, void* env);
/// @}

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

typedef enum {
  SR_Mat_SymmDefinite = 100,
  SR_Mat_SymmSemiDefinite,
  SR_Mat_Symm,
  SR_Mat_General,
  SR_Mat_General_Singular
} SR_eMatClass;

typedef enum {
  SR_Auto = 200,
  SR_CG, 
  SR_BiCGStab, 
  SR_Cheby, 
  SR_ChebyCG, 
  SR_MINRES, 
  SR_GMRES,
  SR_CG_MKL, 
  SR_FGMRES_MKL, 
  SR_LSQR, 
  SR_LSMR,
} SR_eSolver;

typedef enum {
  SR_Precond_None = 300,
  SR_Precond_Jacobi,
  SR_Precond_LU_SGS,
} SR_ePrecond;

/// @{
SR_API void SR_LinSolveR(SR_tMesh mesh,
    SR_tFieldR x, SR_tFieldR b, 
    SR_tMatVecFuncR MatVec, void* env,
    SR_eSolver solver, SR_ePrecond precond, 
    SR_tMatVecFuncR MatVec_H, void* env_H);
SR_API void SR_LinSolveC(SR_tMesh mesh,
    SR_tFieldC x, SR_tFieldC b, 
    SR_tMatVecFuncR MatVec, void* env,
    SR_eSolver solver, SR_ePrecond precond, 
    SR_tMatVecFuncR MatVec_H, void* env_H);
#if SR_C11
#define SR_LinSolve(mesh, x, b, MatVec, env, Solver, Precond, ...) \
  SR_FIELD_GENERIC(x, SR_LinSolve)( \
    mesh, x, b, MatVec, env, Solver, Precond, ##__VA_ARGS__)
#endif
/// @}

/// @{
SR_API SR_eRequest SR_RCI_LinSolveR(SR_tMesh mesh,
    SR_tFieldR x, SR_tFieldR b, 
    SR_eSolver Solver, SR_ePrecond Precond,
    SR_tFieldR* pAy, SR_tFieldR* pY);
#define SR_RCI_LinSolve SR_RCI_LinSolveR 
/// @}

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

/// @{
SR_API void SR_ApplyBCsR(SR_tMesh mesh, SR_tFieldR u,
    SR_INTEGER iBC, SR_REAL alpha, SR_REAL beta, SR_REAL gamma);
SR_API void SR_ApplyBCsC(SR_tMesh mesh, SR_tFieldR u,
    SR_INTEGER iBC, SR_COMPLEX alpha, SR_COMPLEX beta, SR_COMPLEX gamma);
SR_API void SR_ApplyBCsS(SR_tMesh mesh, SR_tFieldR u,
    SR_INTEGER iBC, SR_SYMBOL alpha, SR_SYMBOL beta, SR_SYMBOL gamma);
#if SR_C11
#define SR_ApplyBCs(mesh, u, iBC, ...) \
  SR_FIELD_GENERIC_EXT(u, SR_ApplyBCs)(mesh, u, iBC, ##__VA_ARGS__)
#endif
/// @}

#define SR_ALL 0
#define SR_PURE_DIRICHLET 1.0, 0.0, 0.0
#define SR_PURE_NEUMANN 0.0, 1.0, 0.0

/// @{
SR_API void SR_GradR(SR_tMesh mesh, 
    SR_tFieldR vVec, SR_REAL lambda, SR_tFieldR u);
SR_API void SR_GradC(SR_tMesh mesh,
    SR_tFieldC vVec, SR_COMPLEX lambda, SR_tFieldC u);
SR_API void SR_GradS(SR_tMesh mesh,
    SR_tFieldS vVec, SR_SYMBOL lambda, SR_tFieldS u);
#if SR_C11
#define SR_Grad(mesh, vVec, lambda, u) \
  SR_FIELD_GENERIC_EXT(vVec, SR_Grad)(mesh, vVec, lambda, u)
#endif
/// @}

/// @{
SR_API void SR_DivR(SR_tMesh mesh, 
    SR_tFieldR v, SR_REAL lambda, SR_tFieldR uVec);
SR_API void SR_DivC(SR_tMesh mesh,
    SR_tFieldC v, SR_COMPLEX lambda, SR_tFieldC uVec);
SR_API void SR_DivS(SR_tMesh mesh,
    SR_tFieldS v, SR_SYMBOL lambda, SR_tFieldS uVec);
#if SR_C11
#define SR_Div(mesh, v, lambda, uVec) \
  SR_FIELD_GENERIC_EXT(v, SR_Div)(mesh, v, lambda, uVec)
#endif
/// @}

/// @{
SR_API void SR_ConvR(SR_tMesh mesh, 
    SR_tFieldR v, SR_REAL lambda, SR_tFieldR u, SR_tFieldR a);
SR_API void SR_ConvC(SR_tMesh mesh,
    SR_tFieldC v, SR_REAL lambda, SR_tFieldC u, SR_tFieldR a);
SR_API void SR_ConvS(SR_tMesh mesh,
    SR_tFieldS v, SR_REAL lambda, SR_tFieldS u, SR_tFieldR a);
#if SR_C11
#define SR_Conv(mesh, v, lambda, u, a) \
  SR_FIELD_GENERIC_EXT(v, SR_Conv)(mesh, v, lambda, u, a)
#endif
/// @}

/// @{
SR_API void SR_DivGradR(SR_tMesh mesh,
    SR_tFieldR v, SR_REAL lambda, SR_tFieldR u);
SR_API void SR_DivGradC(SR_tMesh mesh,
    SR_tFieldC v, SR_REAL lambda, SR_tFieldC u);
SR_API void SR_DivGradS(SR_tMesh mesh,
    SR_tFieldS v, SR_REAL lambda, SR_tFieldS u);
#if SR_C11
#define SR_DivGrad(mesh, v, lambda, u) \
  SR_FIELD_GENERIC_EXT(v, SR_DivGrad)(mesh, v, lambda, u)
#endif
/// @}

/// @{
SR_API void SR_DivKGradR(SR_tMesh mesh,
    SR_tFieldR v, SR_REAL lambda, SR_tFieldR k, SR_tFieldR u);
SR_API void SR_DivKGradC(SR_tMesh mesh,
    SR_tFieldC v, SR_REAL lambda, SR_tFieldR k, SR_tFieldC u);
SR_API void SR_DivKGradS(SR_tMesh mesh,
    SR_tFieldS v, SR_REAL lambda, SR_tFieldR k, SR_tFieldS u);
#if SR_C11
#define SR_DivKGrad(mesh, v, lambda, k, u) \
  SR_FIELD_GENERIC_EXT(v, SR_DivGrad)(mesh, v, lambda, k, u)
#endif
/// @}
