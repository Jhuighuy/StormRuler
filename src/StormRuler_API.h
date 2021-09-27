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

#include <assert.h>
#include <stdlib.h>
#include <complex.h>

#define SR_C11 (__STDC_VERSION__ >= 201112L)

#ifdef __cplusplus
#define SR_API extern "C"
#define SR_INL inline

#else
#define SR_API extern
#define SR_INL static
#endif

#define SR_REAL double
#define SR_COMPLEX complex double

#define SR_OPAQUE_STRUCT struct { int _; }

typedef SR_OPAQUE_STRUCT* SR_tSymbol;

typedef SR_OPAQUE_STRUCT* SR_tMesh;

typedef SR_OPAQUE_STRUCT* SR_tIOList;

typedef SR_OPAQUE_STRUCT* SR_tFieldR;
typedef SR_OPAQUE_STRUCT* SR_tFieldC;
typedef SR_OPAQUE_STRUCT* SR_tFieldS;

#define SR_FIELD_GENERIC(x, R, C) \
  _Generic((x), SR_tFieldR: R, SR_tFieldC: C)
#define SR_FIELD_GENERIC_EX(x, R, C, S) \
  _Generic((x), SR_tFieldR: R, SR_tFieldC: C, SR_tFieldS: S)

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

SR_API SR_tMesh SR_InitMesh();

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

/// @{
SR_API SR_tFieldR SR_AllocR(SR_tMesh mesh, int numVars, int rank);
SR_API SR_tFieldC SR_AllocC(SR_tMesh mesh, int numVars, int rank);
SR_API SR_tFieldS SR_AllocS(SR_tMesh mesh, int numVars, int rank);
/// @}

/// @{
SR_API SR_tFieldR SR_Alloc_MoldR(SR_tFieldR mold);
SR_API SR_tFieldC SR_Alloc_MoldC(SR_tFieldC mold);
SR_API SR_tFieldS SR_Alloc_MoldS(SR_tFieldS mold);
#if SR_C11
#define SR_Alloc_Mold(mold) \
  SR_FIELD_GENERIC_EX(mold, SR_Alloc_MoldR, \
    SR_Alloc_MoldC, SR_Alloc_MoldS)(mold)
#endif
/// @}

/// @{
SR_API void SR_FreeR(SR_tFieldR x);
SR_API void SR_FreeC(SR_tFieldC x);
SR_API void SR_FreeS(SR_tFieldS x);
#if SR_C11
#define SR_Free(x) \
  SR_FIELD_GENERIC_EX(x, SR_FreeR, SR_FreeC, SR_FreeS)(x)
#endif
/// @}

/// @{
SR_INL void SR_SwapR(SR_tFieldR* pX, SR_tFieldR* pY) {
  SR_tFieldR z = *pX; *pX = *pY, *pY = z;
}
SR_INL void SR_SwapC(SR_tFieldC* pX, SR_tFieldC* pY) {
  SR_tFieldC z = *pX; *pX = *pY, *pY = z;
}
SR_INL void SR_SwapS(SR_tFieldS* pX, SR_tFieldS* pY) {
  SR_tFieldS z = *pX; *pX = *pY, *pY = z;
}
#if SR_C11
#define SR_Swap(pX, pY) \
  SR_FIELD_GENERIC_EX(*pX, SR_SwapR, SR_SwapC, SR_SwapS)(pX, pY)
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
    SR_tFieldS x, SR_REAL alpha, SR_tSymbol beta);
#if SR_C11
#define SR_Fill(mesh, x, alpha, beta) \
  SR_FIELD_GENERIC_EX(x, SR_FillR, SR_FillC, SR_FillS)( \
    mesh, x, alpha, beta)
#endif
/// @}

/// @{
SR_API void SR_Fill_RandomR(SR_tMesh mesh, 
    SR_tFieldR x, SR_REAL alpha, SR_REAL beta);
SR_API void SR_Fill_RandomC(SR_tMesh mesh, 
    SR_tFieldC x, SR_REAL alpha, SR_COMPLEX beta);
#if SR_C11
#define SR_Fill_Random(mesh, x, alpha, beta) \
  SR_FIELD_GENERIC(x, SR_Fill_RandomR, SR_Fill_RandomC)( \
    mesh, x, alpha, beta)
#endif
/// @}

/// @{
SR_API void SR_SetR(SR_tMesh mesh, SR_tFieldR y, SR_tFieldR x);
SR_API void SR_SetC(SR_tMesh mesh, SR_tFieldC y, SR_tFieldC x);
SR_API void SR_SetS(SR_tMesh mesh, SR_tFieldS y, SR_tFieldS x);
#if SR_C11
#define SR_Set(mesh, y, x) \
  SR_FIELD_GENERIC_EX(y, SR_SetR, SR_SetC, SR_SetS)( \
    mesh, y, x)
#endif
/// @}

/// @{
SR_API void SR_ScaleR(SR_tMesh mesh, 
    SR_tFieldR y, SR_tFieldR x, SR_REAL alpha);
SR_API void SR_ScaleC(SR_tMesh mesh, 
    SR_tFieldC y, SR_tFieldC x, SR_COMPLEX alpha);
SR_API void SR_ScaleS(SR_tMesh mesh, 
    SR_tFieldS y, SR_tFieldS x, SR_tSymbol alpha);
#if SR_C11
#define SR_Scale(mesh, y, x, alpha) \
  SR_FIELD_GENERIC_EX(y, SR_ScaleR, SR_ScaleC, SR_ScaleS)( \
    mesh, y, x, alpha)
#endif
/// @}

/// @{
SR_API void SR_AddR(SR_tMesh mesh, SR_tFieldR z,
    SR_tFieldR y, SR_tFieldR x, SR_REAL alpha, SR_REAL beta);
SR_API void SR_AddC(SR_tMesh mesh, SR_tFieldC z,
    SR_tFieldC y, SR_tFieldC x, SR_COMPLEX alpha, SR_COMPLEX beta);
SR_API void SR_AddS(SR_tMesh mesh, SR_tFieldS z,
    SR_tFieldS y, SR_tFieldS x, SR_tSymbol alpha, SR_tSymbol beta);
#if SR_C11
#define SR_Add(mesh, z, y, x, alpha, beta) \
  SR_FIELD_GENERIC_EX(y, SR_AddR, SR_AddC, SR_AddS)( \
    mesh, z, y, x, alpha, beta)
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
    SR_tFieldS y, SR_tFieldS x, SR_tSymbol alpha, SR_tSymbol beta) {
  assert(0); //SR_AddS(mesh, z, y, x, -alpha, beta);
}
#if SR_C11
#define SR_Sub(mesh, z, y, x, alpha, beta) \
  SR_FIELD_GENERIC_EX(y, SR_SubR, SR_SubC, SR_SubS)( \
    mesh, z, y, x, alpha, beta)
#endif
/// @}

SR_API void SR_MulR(SR_tMesh mesh, SR_tFieldR z,
    SR_tFieldR y, SR_tFieldR x);
#define SR_Mul SR_MulR

/// @{
typedef void(*SR_tMapFuncR)(int size, 
    SR_REAL* Fx, const SR_REAL* x, void* env);
typedef void(*SR_tMapFuncC)(int size, 
    SR_COMPLEX* Fx, const SR_COMPLEX* x, void* env);
typedef void(*SR_tMapFuncS)(int size, 
    SR_tSymbol* Fx, const SR_tSymbol* x, void* env);
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
  SR_FIELD_GENERIC_EX(y, SR_FuncProdR, SR_FuncProdC, SR_FuncProdS)( \
    mesh, y, x, f, env)
#endif
/// @}

/// @{
typedef void(*SR_tSMapFuncR)(int dim, const SR_REAL* r,
    int size, SR_REAL* Fx, const SR_REAL* x, void* env);
typedef void(*SR_tSMapFuncC)(int dim, const SR_REAL* r,
    int size, SR_COMPLEX* Fx, const SR_COMPLEX* x, void* env);
typedef void(*SR_tSMapFuncS)(int dim, const SR_REAL* r,
    int size, SR_tSymbol* Fx, const SR_tSymbol* x, void* env);
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
  SR_FIELD_GENERIC_EX(y, SR_SFuncProdR, SR_SFuncProdC, SR_SFuncProdS)( \
    mesh, y, x, f, env)
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
} SR_tMatClass;

typedef enum {
  SR_Auto = 200,
  SR_CG, SR_BiCGStab, 
  SR_Cheby, SR_ChebyCG, 
  SR_MINRES, SR_GMRES,
  SR_CG_MKL, SR_FGMRES_MKL, 
  SR_LSQR, SR_LSMR,
} SR_tSolver;

typedef enum {
  SR_Precond_None = 300,
  SR_Precond_Jacobi,
  SR_Precond_LU_SGS,
} SR_tPrecond;

/// @{
SR_API void SR_LinSolveR(SR_tMesh mesh,
    SR_tFieldR x, SR_tFieldR b, SR_tMatVecFuncR MatVec, void* env,
    SR_tSolver Solver, SR_tPrecond Precond, SR_tMatVecFuncR MatVec_H, void* env_H);
SR_API void SR_LinSolveC(SR_tMesh mesh,
    SR_tFieldC x, SR_tFieldC b, SR_tMatVecFuncR MatVec, void* env,
    SR_tSolver Solver, SR_tPrecond Precond, SR_tMatVecFuncR MatVec_H, void* env_H);
#if SR_C11
#define SR_LinSolve(mesh, x, b, MatVec, env, Solver, Precond, ...) \
  SR_FIELD_GENERIC(x, SR_LinSolveR, SR_LinSolveC)( \
    mesh, x, b, MatVec, env, Solver, Precond, ##__VA_ARGS__)
#endif
/// @}

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

/// @{
SR_API void SR_ApplyBCsR(SR_tMesh mesh, SR_tFieldR u,
    int iBC, SR_REAL alpha, SR_REAL beta, SR_REAL gamma);
SR_API void SR_ApplyBCsC(SR_tMesh mesh, SR_tFieldR u,
    int iBC, SR_COMPLEX alpha, SR_COMPLEX beta, SR_COMPLEX gamma);
SR_API void SR_ApplyBCsS(SR_tMesh mesh, SR_tFieldR u,
    int iBC, SR_tSymbol alpha, SR_tSymbol beta, SR_tSymbol gamma);
#if SR_C11
#define SR_ApplyBCs(mesh, u, iBC, ...) \
  SR_FIELD_GENERIC_EX(u, SR_ApplyBCsR, SR_ApplyBCsC, \
    SR_ApplyBCsS)(mesh, u, iBC, ##__VA_ARGS__)
#endif
/// @}

#define SR_ALL 0
#define SR_PURE_DIRICHLET 1.0, 0.0, 0.0
#define SR_PURE_NEUMANN 0.0, 1.0, 0.0

/// @{
SR_API void SR_GradR(SR_tMesh mesh, 
    SR_tFieldR vVec, SR_REAL lambda, SR_tFieldR u);
SR_API void SR_GradC(SR_tMesh mesh,
    SR_tFieldC vVec, SR_REAL lambda, SR_tFieldC u);
SR_API void SR_GradS(SR_tMesh mesh,
    SR_tFieldS vVec, SR_REAL lambda, SR_tFieldS u);
#if SR_C11
#define SR_Grad(mesh, vVec, lambda, u) \
  SR_FIELD_GENERIC_EX(vVec, SR_GradR, SR_GradC, SR_GradS)( \
    mesh, vVec, lambda, u)
#endif
/// @}

/// @{
SR_API void SR_DivR(SR_tMesh mesh, 
    SR_tFieldR v, SR_REAL lambda, SR_tFieldR uVec);
SR_API void SR_DivC(SR_tMesh mesh,
    SR_tFieldC v, SR_REAL lambda, SR_tFieldC uVec);
SR_API void SR_DivS(SR_tMesh mesh,
    SR_tFieldS v, SR_REAL lambda, SR_tFieldS uVec);
#if SR_C11
#define SR_Div(mesh, v, lambda, uVec) \
  SR_FIELD_GENERIC_EX(v, SR_DivR, SR_DivC, SR_DivS)( \
    mesh, v, lambda, uVec)
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
  SR_FIELD_GENERIC_EX(v, SR_ConvR, SR_ConvC, SR_ConvS)( \
    mesh, v, lambda, u, a)
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
  SR_FIELD_GENERIC_EX(v, SR_DivGradR, \
    SR_DivGradC, SR_DivGradS)(mesh, v, lambda, u)
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
  SR_FIELD_GENERIC_EX(v, SR_DivGradR, \
    SR_DivGradC, SR_DivGradS)(mesh, v, lambda, k, u)
#endif
/// @}
