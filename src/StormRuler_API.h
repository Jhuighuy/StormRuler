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
#define STORM_CXX_ 1
#else
#define STORM_CXX_ 0
#endif

#if STORM_C11_ && STORM_CXX_
#error StormRuler API: both C11 or C++ targets found.
#endif

#if !STORM_C11_ && !STORM_CXX_
#error StormRuler API: neither C11 nor C++ targets found.
#endif

#define SR_DEFINE_OPAQUE_(type) \
  typedef struct type##Struct type##Struct; \
  typedef type##Struct* type

#define SR_DEFINE_FUNCPTR_(type, name, ...) typedef type(*name)(__VA_ARGS__)

#if STORM_C11_
#define SR_API extern
#define SR_INL static
#elif STORM_CXX_
#define SR_API extern "C"
#define SR_INL inline
#endif

typedef int stormInt_t;
typedef double stormReal_t;
typedef const char* stormString_t;
#if STORM_C11_
#include <complex.h>
typedef stormReal_t complex stormComplex_t;
#elif STORM_CXX_
#include <complex>
typedef std::complex<stormReal_t> stormComplex_t;
#endif

#if STORM_C11_
#define STORM_DEFAULT_(decl, ...) decl
#elif STORM_CXX_
#define STORM_DEFAULT_(decl, ...) decl = __VA_ARGS__
#endif

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

SR_DEFINE_OPAQUE_(stormSymbol_t);

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

SR_DEFINE_OPAQUE_(stormMesh_t);

SR_API stormMesh_t SR_InitMesh(void);

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

SR_DEFINE_OPAQUE_(stormArrayR_t);
SR_DEFINE_OPAQUE_(stormArrayC_t);
SR_DEFINE_OPAQUE_(stormArrayS_t);

#if STORM_C11_
#define STORM_GENERIC_(x, func) \
  _Generic((x), stormArrayR_t: func##R, stormArrayC_t: func##C)
#define STORM_GENERIC_EXT_(x, func) \
  _Generic((x), stormArrayR_t: func##R, stormArrayC_t: func##C, stormArrayS_t: func##S)
#elif STORM_CXX_
#if 1
#define SR__Generic(x, f, ...) f
#else
template<typename R, typename C, typename S=void*>
inline auto SR__Generic(stormArrayR_t, R r, C, S=nullptr) { return r; }
template<typename R, typename C, typename S=void*>
inline auto SR__Generic(stormArrayC_t, R, C c, S=nullptr) { return c; }
template<typename R, typename C, typename S>
inline auto SR__Generic(stormArrayS_t, R, C, S s) { return s; }
#endif
#define STORM_GENERIC_(x, func) \
  ( SR__Generic(x, func##R, func##C) )
#define STORM_GENERIC_EXT_(x, func, ...) \
  ( SR__Generic(x, func##R, func##C, func##S) )
#endif

/// @{
SR_API stormArrayR_t SR_AllocR(
  stormMesh_t mesh, stormInt_t numVars, stormInt_t rank);
SR_API stormArrayC_t SR_AllocC(
  stormMesh_t mesh, stormInt_t numVars, stormInt_t rank);
SR_API stormArrayS_t SR_AllocS(
  stormMesh_t mesh, stormInt_t numVars, stormInt_t rank);
/// @}

/// @{
SR_API stormArrayR_t SR_Alloc_MoldR(stormArrayR_t mold);
SR_API stormArrayC_t SR_Alloc_MoldC(stormArrayC_t mold);
SR_API stormArrayS_t SR_Alloc_MoldS(stormArrayS_t mold);
#define SR_Alloc_Mold(mold) \
  STORM_GENERIC_EXT_(mold, SR_Alloc_Mold)(mold)
/// @}

/// @{
SR_API void stormFreeR(stormArrayR_t x);
SR_API void stormFreeC(stormArrayC_t x);
SR_API void stormFreeS(stormArrayS_t x);
#define stormFree(x) STORM_GENERIC_EXT_(x, stormFree)(x)
/// @}

/// @{
#if STORM_C11_
SR_INL void stormSwapP(void** pX, void** pY) {
  void* z = *pX; *pX = *pY, *pY = z;
}
#define stormSwap(x, y) stormSwapP((void**)(&x), (void**)(&y))
#elif STORM_CXX_
#define stormSwap(x, y) std::swap(x, y)
#endif
/// @}

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

SR_DEFINE_OPAQUE_(SR_tIOList);

SR_API SR_tIOList SR_IO_Begin();

SR_API void SR_IO_Add(SR_tIOList IO, 
  stormArrayR_t x, const char* name);

SR_API void SR_IO_Flush(SR_tIOList IO, 
  stormMesh_t mesh, const char* filename);

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

/// @{
SR_API void stormFillR(stormMesh_t mesh, 
                       stormArrayR_t x, 
                       stormReal_t alpha, 
                       STORM_DEFAULT_(stormReal_t beta, 0.0));
SR_API void stormFillC(stormMesh_t mesh, 
                       stormArrayC_t x, 
                       stormReal_t alpha, 
                       STORM_DEFAULT_(stormComplex_t beta, 0.0));
SR_API void stormFillS(stormMesh_t mesh, 
                       stormArrayS_t x, 
                       stormReal_t alpha, 
                       STORM_DEFAULT_(stormSymbol_t beta, {}));
#define stormFill(mesh, x, alpha, beta) \
  STORM_GENERIC_EXT_(x, stormFill)(mesh, x, alpha, beta)
/// @}

/// @{
SR_API void stormRandFillR(stormMesh_t mesh, 
                           stormArrayR_t x, 
                           STORM_DEFAULT_(stormReal_t a, -1.0), 
                           STORM_DEFAULT_(stormReal_t b, +1.0));
SR_API void stormRandFillC(stormMesh_t mesh, 
                           stormArrayC_t x, 
                           STORM_DEFAULT_(stormReal_t a, -1.0), 
                           STORM_DEFAULT_(stormReal_t b, +1.0));
#define stormRandFill(mesh, x, alpha, beta) \
  STORM_GENERIC_(x, stormRandFill)(mesh, x, alpha, beta)
/// @}

/// @{
SR_API void stormSetR(stormMesh_t mesh, 
                      stormArrayR_t y, 
                      stormArrayR_t x);
SR_API void stormSetC(stormMesh_t mesh, 
                      stormArrayC_t y,
                      stormArrayC_t x);
SR_API void stormSetS(stormMesh_t mesh, 
                      stormArrayS_t y, 
                      stormArrayS_t x);
#define stormSet(mesh, y, x) \
  STORM_GENERIC_EXT_(y, stormSet)(mesh, y, x)
/// @}

/// @{
SR_API void stormScaleR(stormMesh_t mesh, 
                        stormArrayR_t y, 
                        stormArrayR_t x, 
                        stormReal_t alpha);
SR_API void stormScaleC(stormMesh_t mesh, 
                        stormArrayC_t y, 
                        stormArrayC_t x, 
                        stormComplex_t alpha);
SR_API void stormScaleS(stormMesh_t mesh, 
                        stormArrayS_t y, 
                        stormArrayS_t x, 
                        stormSymbol_t alpha);
#define stormScale(mesh, y, x, alpha) \
  STORM_GENERIC_EXT_(y, stormScale)(mesh, y, x, alpha)
/// @}

/// @{
SR_API void stormAddR(stormMesh_t mesh, 
                      stormArrayR_t z, 
                      stormArrayR_t y,
                      stormArrayR_t x, 
                      STORM_DEFAULT_(stormReal_t alpha, 1.0), 
                      STORM_DEFAULT_(stormReal_t  beta, 1.0));
SR_API void stormAddC(stormMesh_t mesh, 
                      stormArrayC_t z, 
                      stormArrayC_t y,
                      stormArrayC_t x, 
                      STORM_DEFAULT_(stormComplex_t alpha, 1.0), 
                      STORM_DEFAULT_(stormComplex_t  beta, 1.0));
SR_API void stormAddS(stormMesh_t mesh, 
                      stormArrayS_t z, 
                      stormArrayS_t y,
                      stormArrayS_t x, 
                      STORM_DEFAULT_(stormSymbol_t alpha, {}), 
                      STORM_DEFAULT_(stormSymbol_t  beta, {}));
#define stormAdd(mesh, z, y, x, alpha, beta) \
  STORM_GENERIC_EXT_(y, stormAdd)(mesh, z, y, x, alpha, beta)
/// @}

/// @{
// No need to import these from Fortran.
SR_INL void stormSubR(stormMesh_t mesh, 
                      stormArrayR_t z,
                      stormArrayR_t y, 
                      stormArrayR_t x, 
                      STORM_DEFAULT_(stormReal_t alpha, 1.0),
                      STORM_DEFAULT_(stormReal_t  beta, 1.0)) {
  stormAddR(mesh, z, y, x, -alpha, beta);
}
SR_INL void stormSubC(stormMesh_t mesh, 
                      stormArrayC_t z,
                      stormArrayC_t y, 
                      stormArrayC_t x, 
                      STORM_DEFAULT_(stormComplex_t alpha, 1.0), 
                      STORM_DEFAULT_(stormComplex_t  beta, 1.0)) {
  stormAddC(mesh, z, y, x, -alpha, beta);
}
SR_INL void stormSubS(stormMesh_t mesh, 
                      stormArrayS_t z,
                      stormArrayS_t y, 
                      stormArrayS_t x, 
                      STORM_DEFAULT_(stormSymbol_t alpha, {}), 
                      STORM_DEFAULT_(stormSymbol_t  beta, {})) {
  assert(0); //stormAddS(mesh, z, y, x, -alpha, beta);
}
#define stormSub(mesh, z, y, x, alpha, beta) \
  STORM_GENERIC_EXT_(y, stormSub)(mesh, z, y, x, alpha, beta)
/// @}

/// @{
SR_API void stormMulR(stormMesh_t mesh,
                      stormArrayR_t z, 
                      stormArrayR_t y, 
                      stormArrayR_t x);
SR_API void stormMulC(stormMesh_t mesh, 
                      stormArrayC_t z, 
                      stormArrayC_t y, 
                      stormArrayC_t x);
SR_API void stormMulS(stormMesh_t mesh, 
                      stormArrayS_t z, 
                      stormArrayS_t y, 
                      stormArrayS_t x);
#define stormMul(mesh, z, y, x) \
  STORM_GENERIC_EXT_(y, stormMul)(mesh, z, y, x)
/// @}

/// @{
typedef void(*stormMapFuncR_t)(stormInt_t size, 
                               stormReal_t* Fx, 
                               const stormReal_t* x, 
                               void* env);
typedef void(*stormMapFuncC_t)(stormInt_t size, 
                               stormComplex_t* Fx, 
                               const stormComplex_t* x, 
                               void* env);
typedef void(*stormMapFuncS_t)(stormInt_t size, 
                               stormSymbol_t* Fx, 
                               const stormSymbol_t* x, 
                               void* env);
/// @}

/// @{
SR_API stormReal_t SR_IntegrateR(stormMesh_t mesh,
                                 stormArrayR_t x, 
                                 stormMapFuncR_t f, 
                                 void* env);
SR_API stormComplex_t SR_IntegrateC(stormMesh_t mesh,
                                    stormArrayC_t x, 
                                    stormMapFuncC_t f, 
                                    void* env);
#define SR_Integrate(mesh, x, f, env) \
  STORM_GENERIC_(x, SR_Integrate)(mesh, x, f, env)
/// @}

/// @{
SR_API void SR_FuncProdR(stormMesh_t mesh,
                         stormArrayR_t y, 
                         stormArrayR_t x, 
                         stormMapFuncR_t f, 
                         STORM_DEFAULT_(void* env, NULL));
SR_API void SR_FuncProdC(stormMesh_t mesh,
                         stormArrayC_t y, 
                         stormArrayC_t x, 
                         stormMapFuncC_t f, 
                         STORM_DEFAULT_(void* env, NULL));
SR_API void SR_FuncProdS(stormMesh_t mesh,
                         stormArrayS_t y, 
                         stormArrayS_t x, 
                         stormMapFuncS_t f, 
                         STORM_DEFAULT_(void* env, NULL));
#define SR_FuncProd(mesh, y, x, f, env) \
  STORM_GENERIC_EXT_(y, SR_FuncProd)(mesh, y, x, f, env)
/// @}

/// @{
typedef void(*stormSpMapFuncR_t)(stormInt_t dim, 
                                 const stormReal_t* r,
                                 stormInt_t size, 
                                 stormReal_t* Fx, 
                                 const stormReal_t* x, 
                                 void* env);
typedef void(*stormSpMapFuncC_t)(stormInt_t dim, 
                                 const stormReal_t* r,
                                 stormInt_t size, 
                                 stormComplex_t* Fx, 
                                 const stormComplex_t* x, 
                                 void* env);
typedef void(*stormSpMapFuncS_t)(stormInt_t dim, 
                                 const stormReal_t* r,
                                 stormInt_t size, 
                                 stormSymbol_t* Fx, 
                                 const stormSymbol_t* x, 
                                 void* env);
/// @}

/// @{
SR_API void SR_SFuncProdR(stormMesh_t mesh,
                          stormArrayR_t y, 
                          stormArrayR_t x, 
                          stormSpMapFuncR_t f, 
                          STORM_DEFAULT_(void* env, NULL));
SR_API void SR_SFuncProdC(stormMesh_t mesh,
                          stormArrayC_t y, 
                          stormArrayC_t x, 
                          stormSpMapFuncC_t f, 
                          STORM_DEFAULT_(void* env, NULL));
SR_API void SR_SFuncProdS(stormMesh_t mesh,
                          stormArrayS_t y, 
                          stormArrayS_t x, 
                          stormSpMapFuncS_t f, 
                          STORM_DEFAULT_(void* env, NULL));
#define SR_SFuncProd(mesh, y, x, f, env) \
  STORM_GENERIC_EXT_(y, SR_SFuncProd)(mesh, y, x, f, env)
/// @}

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

/// @{
typedef void(*stormMatVecFuncR_t)(stormMesh_t mesh, 
                                  stormArrayR_t Ax, 
                                  stormArrayR_t x, 
                                  void* env);
typedef void(*stormMatVecFuncC_t)(stormMesh_t mesh, 
                                  stormArrayC_t Ax, 
                                  stormArrayC_t x, 
                                  void* env);
typedef void(*stormMatVecFuncS_t)(stormMesh_t mesh, 
                                  stormArrayS_t Ax, 
                                  stormArrayS_t x, 
                                  void* env);
/// @}

/// @{
#define stormCG       "CG"
#define stormBiCGStab "BiCGStab"
#define stormMINRES   "MINRES"
#define stormGMRES    "GMRES"
#define stormTFQMR    "TFQMR"
#define stormLSQR     "LSQR"
#define stormLSMR     "LSMR"
/// @}

/// @{
#define stormNone   ""
#define stormJacobi "Jacobi"
#define stormLU_SGS "LU_SGS"
/// @}

/// @{
SR_API void stormLinSolveR(stormMesh_t mesh,
                           stormString_t method,
                           stormString_t preMethod,
                           stormArrayR_t x, 
                           stormArrayR_t b, 
                           stormMatVecFuncR_t MatVec, 
                           STORM_DEFAULT_(void* env, NULL),
                           STORM_DEFAULT_(stormMatVecFuncR_t MatVec_H, NULL), 
                           STORM_DEFAULT_(void* env_H, NULL));
SR_API void stormLinSolveC(stormMesh_t mesh,
                           stormString_t method, 
                           stormString_t preMethod,
                           stormArrayC_t x, 
                           stormArrayC_t b, 
                           stormMatVecFuncC_t MatVec, 
                           STORM_DEFAULT_(void* env, NULL),
                           STORM_DEFAULT_(stormMatVecFuncC_t MatVec_H, NULL), 
                           STORM_DEFAULT_(void* env_H, NULL));
#define stormLinSolve(mesh, method, preMethod, x, b, MatVec, env, ...) \
  STORM_GENERIC_(x, stormLinSolve)( \
    mesh, method, preMethod, x, b, MatVec, env, ##__VA_ARGS__)
/// @}

/// @{
#define stormNewton "Newton"
#define stormJFNK   "JFNK"
/// @}

/// @{
SR_API void stormNonlinSolveR(stormMesh_t mesh,
                              stormString_t method,
                              stormArrayR_t x,
                              stormArrayR_t b, 
                              stormMatVecFuncR_t MatVec, 
                              STORM_DEFAULT_(void* env, NULL));
SR_API void stormNonlinSolveC(stormMesh_t mesh,
                              stormString_t method,
                              stormArrayC_t x,
                              stormArrayC_t b, 
                              stormMatVecFuncC_t MatVec, 
                              STORM_DEFAULT_(void* env, NULL));
#define stormNonlinSolve(mesh, method, x, b, MatVec, env, ...) \
  STORM_GENERIC_(x, stormNonlinSolve)( \
    mesh, method, x, b, MatVec, env, ##__VA_ARGS__)
/// @}

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

/// @{
SR_API void stormApplyBCsR(stormMesh_t mesh, 
                           stormArrayR_t u,
                           stormInt_t iBC, 
                           stormReal_t alpha, 
                           stormReal_t beta, 
                           stormReal_t gamma);
SR_API void stormApplyBCsC(stormMesh_t mesh, 
                           stormArrayR_t u,
                           stormInt_t iBC, 
                           stormComplex_t alpha, 
                           stormComplex_t beta, 
                           stormComplex_t gamma);
SR_API void stormApplyBCsS(stormMesh_t mesh, 
                           stormArrayR_t u,
                           stormInt_t iBC, 
                           stormSymbol_t alpha, 
                           stormSymbol_t beta, 
                           stormSymbol_t gamma);
#define stormApplyBCs(mesh, u, iBC, ...) \
  STORM_GENERIC_EXT_(u, stormApplyBCs)(mesh, u, iBC, ##__VA_ARGS__)
/// @}

#define SR_ALL 0
#define SR_PURE_DIRICHLET 1.0, 0.0, 0.0
#define SR_PURE_NEUMANN   0.0, 1.0, 0.0

#define SR_DIRICHLET(gamma) 1.0, 0.0, (gamma)
#define SR_NEUMANN(gamma)   0.0, 1.0, (gamma)

/// @{
SR_API void stormApplyBCs_SlipWallR(stormMesh_t mesh, 
                                    stormArrayR_t u,
                                    stormInt_t iBC);
#define stormApplyBCs_SlipWall stormApplyBCs_SlipWallR
/// @}

/// @{
SR_API void stormApplyBCs_InOutLetR(stormMesh_t mesh, 
                                  stormArrayR_t u,
                                  stormInt_t iBC);
#define stormApplyBCs_InOutLet stormApplyBCs_InOutLetR
/// @}

/// @{
SR_API void stormGradientR(stormMesh_t mesh, 
                           stormArrayR_t vVec, 
                           stormReal_t lambda, 
                           stormArrayR_t u);
SR_API void stormGradientC(stormMesh_t mesh,
                           stormArrayC_t vVec, 
                           stormComplex_t lambda, 
                           stormArrayC_t u);
SR_API void stormGradientS(stormMesh_t mesh,
                           stormArrayS_t vVec, 
                           stormSymbol_t lambda, 
                           stormArrayS_t u);
#define stormGradient(mesh, vVec, lambda, u) \
  STORM_GENERIC_EXT_(vVec, stormGradient)(mesh, vVec, lambda, u)
/// @}

/// @{
SR_API void stormDivergenceR(stormMesh_t mesh, 
                             stormArrayR_t v, 
                             stormReal_t lambda, 
                             stormArrayR_t uVec);
SR_API void stormDivergenceC(stormMesh_t mesh,
                             stormArrayC_t v, 
                             stormComplex_t lambda, 
                             stormArrayC_t uVec);
SR_API void stormDivergenceS(stormMesh_t mesh,
                             stormArrayS_t v, 
                             stormSymbol_t lambda, 
                             stormArrayS_t uVec);
#define stormDivergence(mesh, v, lambda, uVec) \
  STORM_GENERIC_EXT_(v, stormDivergence)(mesh, v, lambda, uVec)
/// @}

/// @{
SR_API void stormRhieChowCorrection(stormMesh_t mesh, 
                                    stormArrayR_t v,
                                    stormReal_t lambda, 
                                    stormReal_t tau, 
                                    stormArrayR_t p, 
                                    stormArrayR_t rho);
/// @}

/// @{
SR_API void stormConvectionR(stormMesh_t mesh, 
                             stormArrayR_t v, 
                             stormReal_t lambda, 
                             stormArrayR_t u, 
                             stormArrayR_t a);
SR_API void stormConvectionC(stormMesh_t mesh,
                             stormArrayC_t v, 
                             stormComplex_t lambda, 
                             stormArrayC_t u, 
                             stormArrayR_t a);
SR_API void stormConvectionS(stormMesh_t mesh,
                             stormArrayS_t v, 
                             stormSymbol_t lambda, 
                             stormArrayS_t u, 
                             stormArrayR_t a);
#define stormConvection(mesh, v, lambda, u, a) \
  STORM_GENERIC_EXT_(v, stormConvection)(mesh, v, lambda, u, a)
/// @}

/// @{
SR_API void stormDivGradR(stormMesh_t mesh,
                          stormArrayR_t v, 
                          stormReal_t lambda, 
                          stormArrayR_t u);
SR_API void stormDivGradC(stormMesh_t mesh,
                          stormArrayC_t v, 
                          stormComplex_t lambda, 
                          stormArrayC_t u);
SR_API void stormDivGradS(stormMesh_t mesh,
                          stormArrayS_t v, 
                          stormSymbol_t lambda, 
                          stormArrayS_t u);
#define stormDivGrad(mesh, v, lambda, u) \
  STORM_GENERIC_EXT_(v, stormDivGrad)(mesh, v, lambda, u)
/// @}

/// @{
SR_API void stormDivWGradR(stormMesh_t mesh,
                           stormArrayR_t v, 
                           stormReal_t lambda, 
                           stormArrayR_t k, 
                           stormArrayR_t u);
SR_API void stormDivWGradC(stormMesh_t mesh,
                           stormArrayC_t v, 
                           stormComplex_t lambda, 
                           stormArrayC_t k, 
                           stormArrayC_t u);
SR_API void stormDivWGradS(stormMesh_t mesh,
                           stormArrayS_t v, 
                           stormSymbol_t lambda, 
                           stormArrayR_t k, 
                           stormArrayS_t u);
#define stormDivWGrad(mesh, v, lambda, k, u) \
  STORM_GENERIC_EXT_(v, stormDivWGrad)(mesh, v, lambda, k, u)
/// @}

#endif // ifndef STORM_RULER_API_H
