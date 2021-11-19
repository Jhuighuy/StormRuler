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
#ifndef STORM_RULER_API_H
#define STORM_RULER_API_H

#include <stdint.h>
#include <assert.h>

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

// Detect C11.
#if (__STDC_VERSION__ >= 201112L)
#define STORM_C11_ 1
#else
#define STORM_C11_ 0
#endif

// Detect C++.
#if defined(__cplusplus)
#define STORM_CXX17_ 1
#else
#define STORM_CXX17_ 0
#endif

#if STORM_C11_ && STORM_CXX17_
#error StormRuler API: both C11 or C++17 targets found.
#endif

#if !STORM_C11_ && !STORM_CXX17_
#error StormRuler API: neither C11 nor C++17 targets found.
#endif

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

#define _STORM_CAT2_(x, y) x##y
#define _STORM_CAT_(x, y) _STORM_CAT2_(x, y)

#define _STORM_STR2_(x) #x
#define _STORM_STR_(x) _STORM_STR2_(x)

#define _STORM_OPAQUE_(type) \
  struct type; typedef struct type* type##_t

// Enable default arguments for C++.
#if STORM_C11_
#define _STORM_DEFAULT_(decl, ...) decl
#elif STORM_CXX17_
#define _STORM_DEFAULT_(decl, ...) decl = __VA_ARGS__
#endif

#if STORM_C11_
#define STORM_API extern
#define STORM_INL static
#elif STORM_CXX17_
#define STORM_API extern "C"
#define STORM_INL inline
#endif

#if STORM_C11_
#define STORM_NULL NULL
#elif STORM_CXX17_
#define STORM_NULL nullptr
#endif

// Basic types.
typedef int stormInt_t;
typedef double stormReal_t;
typedef void* stormOpaque_t;
typedef const char* stormString_t;

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

_STORM_OPAQUE_(stormSymbol);

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

_STORM_OPAQUE_(stormMesh);

STORM_API stormMesh_t SR_InitMesh(void);

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

_STORM_OPAQUE_(stormArray);

STORM_API stormArray_t SR_Alloc(
  stormMesh_t mesh, stormInt_t numVars, stormInt_t rank);

STORM_API stormArray_t stormAllocLike(stormArray_t array);

STORM_API void stormFree(stormArray_t x);

/// @{
#if STORM_C11_
STORM_INL void stormSwapP(stormOpaque_t* pX, stormOpaque_t* pY) {
  stormOpaque_t z = *pX; *pX = *pY, *pY = z;
}
#define stormSwap(x, y) stormSwapP((stormOpaque_t*)&(x), (stormOpaque_t*)&(y))
#elif STORM_CXX17_
#define stormSwap(x, y) std::swap(x, y)
#endif
/// @}

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

_STORM_OPAQUE_(stormIOList);

STORM_API stormIOList_t SR_IO_Begin();

STORM_API void SR_IO_Add(stormIOList_t IO,
  stormArray_t x, const char* name);

STORM_API void SR_IO_Flush(stormIOList_t IO,
  stormMesh_t mesh, const char* filename);

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

STORM_API void stormFill(stormMesh_t mesh,
                         stormArray_t x,
                         stormReal_t alpha);

STORM_API void stormRandFill(stormMesh_t mesh,
                             stormArray_t x,
                             _STORM_DEFAULT_(stormReal_t a, -1.0),
                             _STORM_DEFAULT_(stormReal_t b, +1.0));

STORM_API void stormSet(stormMesh_t mesh,
                        stormArray_t y,
                        stormArray_t x);

STORM_API void stormScale(stormMesh_t mesh,
                          stormArray_t y,
                          stormArray_t x,
                          stormReal_t alpha);

STORM_API void stormAdd(stormMesh_t mesh,
                        stormArray_t z,
                        stormArray_t y,
                        stormArray_t x,
                        _STORM_DEFAULT_(stormReal_t alpha, 1.0),
                        _STORM_DEFAULT_(stormReal_t  beta, 1.0));

// No need to import these from Fortran.
STORM_INL void stormSub(stormMesh_t mesh,
                        stormArray_t z,
                        stormArray_t y,
                        stormArray_t x,
                        _STORM_DEFAULT_(stormReal_t alpha, 1.0),
                        _STORM_DEFAULT_(stormReal_t  beta, 1.0)) {
  stormAdd(mesh, z, y, x, -alpha, beta);
}

STORM_API void stormMul(stormMesh_t mesh,
                        stormArray_t z,
                        stormArray_t y,
                        stormArray_t x);

typedef void(*stormMapFunc_t)(stormInt_t size,
                              stormReal_t* Fx,
                              const stormReal_t* x,
                              stormOpaque_t env);

STORM_API stormReal_t stormIntegrate(stormMesh_t mesh,
                                     stormArray_t x,
                                     stormMapFunc_t f,
                                     stormOpaque_t env);

STORM_API void stormFuncProd(stormMesh_t mesh,
                             stormArray_t y,
                             stormArray_t x,
                             stormMapFunc_t f,
                             _STORM_DEFAULT_(stormOpaque_t env, STORM_NULL));

typedef void(*stormSpMapFunc_t)(stormInt_t dim,
                                const stormReal_t* r,
                                stormInt_t size,
                                stormReal_t* Fx,
                                const stormReal_t* x,
                                stormOpaque_t env);

STORM_API void stormSpFuncProd(stormMesh_t mesh,
                               stormArray_t y,
                               stormArray_t x,
                               stormSpMapFunc_t f,
                               _STORM_DEFAULT_(stormOpaque_t env, STORM_NULL));

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

typedef void(*stormMatVecFunc_t)(stormMesh_t mesh,
                                 stormArray_t Ax,
                                 stormArray_t x,
                                 stormOpaque_t env);

/// @{
#define STORM_CG       "CG"
#define STORM_BiCGStab "BiCGStab"
#define STORM_MINRES   "MINRES"
#define STORM_GMRES    "GMRES"
#define STORM_QMR      "QMR"
#define STORM_TFQMR    "TFQMR"
#define STORM_LSQR     "LSQR"
#define STORM_LSMR     "LSMR"
/// @}

/// @{
#define STORM_NONE   ""
#define STORM_JACOBI "Jacobi"
#define STORM_LU_SGS "LU_SGS"
/// @}

STORM_API void stormLinSolve(stormMesh_t mesh,
                             stormString_t method,
                             stormString_t preMethod,
                             stormArray_t x,
                             stormArray_t b,
                             stormMatVecFunc_t matVec,
                             _STORM_DEFAULT_(stormOpaque_t env, STORM_NULL),
                             _STORM_DEFAULT_(stormMatVecFunc_t conjMatVec, STORM_NULL),
                             _STORM_DEFAULT_(stormOpaque_t conjEnv, STORM_NULL));

#if STORM_CXX17_
/// @{
template<typename stormMatVecFuncT_t>
STORM_INL void stormLinSolve(stormMesh_t mesh,
                             stormString_t method,
                             stormString_t preMethod,
                             stormArray_t x,
                             stormArray_t b,
                             stormMatVecFuncT_t matVec) {
  stormLinSolve(mesh, method, preMethod, x, b,
    [](stormMesh_t mesh, stormArray_t Ax, stormArray_t x, stormOpaque_t env) {
      (*static_cast<stormMatVecFuncT_t*>(env))(mesh, Ax, x);
    }, &matVec, STORM_NULL, STORM_NULL);
} // stormLinSolve
template<typename stormMatVecFuncT1_t, typename stormMatVecFuncT2_t>
STORM_INL void stormLinSolve(stormMesh_t mesh,
                             stormString_t method,
                             stormString_t preMethod,
                             stormArray_t x,
                             stormArray_t b,
                             stormMatVecFuncT1_t matVec,
                             stormMatVecFuncT2_t conjMatVec) {
  stormLinSolve(mesh, method, preMethod, x, b,
    [](stormMesh_t mesh, stormArray_t Ax, stormArray_t x, stormOpaque_t env) {
      (*static_cast<stormMatVecFuncT1_t*>(env))(mesh, Ax, x);
    }, &matVec,
    [](stormMesh_t mesh, stormArray_t Ax, stormArray_t x, stormOpaque_t env) {
      (*static_cast<stormMatVecFuncT2_t*>(env))(mesh, Ax, x);
    }, &conjMatVec);
} // stormLinSolve
/// @}
#endif

/// @{
#define STORM_NEWTON "Newton"
#define STORM_JFNK   "JFNK"
/// @}

STORM_API void stormNonlinSolve(stormMesh_t mesh,
                                stormString_t method,
                                stormArray_t x,
                                stormArray_t b,
                                stormMatVecFunc_t matVec,
                                _STORM_DEFAULT_(stormOpaque_t env, STORM_NULL));

#if STORM_CXX17_
template<typename stormMatVecFuncT_t>
STORM_INL void stormNonlinSolve(stormMesh_t mesh,
                                stormString_t method,
                                stormArray_t x,
                                stormArray_t b,
                                stormMatVecFuncT_t matVec) {
  stormNonlinSolve(mesh, method, x, b,
    [](stormMesh_t mesh, stormArray_t Ax, stormArray_t x, stormOpaque_t env) {
      (*static_cast<stormMatVecFuncT_t*>(env))(mesh, Ax, x);
    }, &matVec);
} // stormNonlinSolve
#endif

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

STORM_API void stormApplyBCs(stormMesh_t mesh,
                             stormArray_t u,
                             stormInt_t iBC,
                             stormReal_t alpha,
                             stormReal_t beta,
                             stormReal_t gamma);

#define SR_ALL 0
#define SR_PURE_DIRICHLET 1.0, 0.0, 0.0
#define SR_PURE_NEUMANN   0.0, 1.0, 0.0

#define SR_DIRICHLET(gamma) 1.0, 0.0, (gamma)
#define SR_NEUMANN(gamma)   0.0, 1.0, (gamma)

STORM_API void stormApplyBCs_SlipWall(stormMesh_t mesh,
                                      stormArray_t u,
                                      stormInt_t iBC);

STORM_API void stormApplyBCs_InOutLet(stormMesh_t mesh,
                                      stormArray_t u,
                                      stormInt_t iBC);

STORM_API void stormGradient(stormMesh_t mesh,
                             stormArray_t vVec,
                             stormReal_t lambda,
                             stormArray_t u);

STORM_API void stormDivergence(stormMesh_t mesh,
                               stormArray_t v,
                               stormReal_t lambda,
                               stormArray_t uVec);

STORM_API void stormDivGrad(stormMesh_t mesh,
                            stormArray_t v,
                            stormReal_t lambda,
                            stormArray_t u);

STORM_API void stormDivWGrad(stormMesh_t mesh,
                             stormArray_t v,
                             stormReal_t lambda,
                             stormArray_t k,
                             stormArray_t u);

STORM_API void stormConvection(stormMesh_t mesh,
                               stormArray_t v,
                               stormReal_t lambda,
                               stormArray_t u,
                               stormArray_t aVec);

STORM_API void stormRhieChowCorrection(stormMesh_t mesh,
                                       stormArray_t v,
                                       stormReal_t lambda,
                                       stormReal_t tau,
                                       stormArray_t p,
                                       stormArray_t rho);

#endif // ifndef STORM_RULER_API_H
