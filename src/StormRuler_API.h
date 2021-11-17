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

#define _STORM_CAT2_(x, y) x##y
#define _STORM_CAT_(x, y) _STORM_CAT2_(x, y)

#define _STORM_STR2_(x) #x
#define _STORM_STR_(x) _STORM_STR2_(x)

#define _STORM_OPAQUE_(type) \
  struct type; \
  typedef struct type* type##_t

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

typedef int stormInt_t;
typedef double stormReal_t;
typedef const char* stormString_t;
#if STORM_C11_
#include <complex.h>
typedef stormReal_t complex stormComplex_t;
#elif STORM_CXX17_
#include <complex>
typedef std::complex<stormReal_t> stormComplex_t;
#endif

#if STORM_C11_
#define _STORM_DEFAULT_(decl, ...) decl
#elif STORM_CXX17_
#define _STORM_DEFAULT_(decl, ...) decl = __VA_ARGS__
#endif

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

_STORM_OPAQUE_(stormSymbol);

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

_STORM_OPAQUE_(stormMesh);

STORM_API stormMesh_t SR_InitMesh(void);

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

_STORM_OPAQUE_(stormArrayR);
_STORM_OPAQUE_(stormArrayC);
_STORM_OPAQUE_(stormArrayS);

#if STORM_C11_

#define _STORM_GENERIC_(selectOn, func) \
  ( _Generic((selectOn), \
             stormArrayR_t: _STORM_CAT_(func, R), \
             stormArrayC_t: _STORM_CAT_(func, C)) )

#define _STORM_GENERIC_EXT_(selectOn, func)(...) \
  ( _Generic((selectOn), \
             stormArrayR_t: _STORM_CAT_(func, R), \
             stormArrayC_t: _STORM_CAT_(func, C), \
             stormArrayS_t: _STORM_CAT_(func, S)) )

#elif STORM_CXX17_

#include <type_traits>

#define _STORM_GENERIC_(selectOn, func) \
  [&](auto... args) -> auto { \
    /* Kind of _Generic implementation. */ \
    typedef std::remove_cv_t<decltype(selectOn)> T; \
    if constexpr (std::is_same_v<T, stormArrayR_t>) { \
      return _STORM_CAT_(func, R)(args...); \
    } else if constexpr (std::is_same_v<T, stormArrayC_t>) { \
      return _STORM_CAT_(func, C)(args...); \
    } else { \
      static_assert(sizeof...(args) != sizeof...(args), \
        "Unexpected type of argument `" _STORM_STR2_(x) "` passed, " \
        "expected `stormArrayR_t` or `stormArrayC_t`."); \
    } \
  }

#define _STORM_GENERIC_EXT_(selectOn, func) \
  [&](auto... args) -> auto { \
    /* Kind of _Generic implementation. */ \
    typedef std::remove_cv_t<decltype(selectOn)> T; \
    if constexpr (std::is_same_v<T, stormArrayR_t>) { \
      return _STORM_CAT_(func, R)(args...); \
    } else if constexpr (std::is_same_v<T, stormArrayC_t>) { \
      return _STORM_CAT_(func, C)(args...); \
    } else if constexpr (std::is_same_v<T, stormArrayS_t>) { \
      return _STORM_CAT_(func, S)(args...); \
    } else { \
      static_assert(sizeof...(args) != sizeof...(args), \
        "Unexpected type of argument `" _STORM_STR2_(x) "` passed, " \
        "expected `stormArrayR_t`, `stormArrayC_t` or `stormArrayS_t`."); \
    } \
  }

#endif // STORM_C11_ or STORM_CXX17_

/// @{
STORM_API stormArrayR_t SR_AllocR(
  stormMesh_t mesh, stormInt_t numVars, stormInt_t rank);
STORM_API stormArrayC_t SR_AllocC(
  stormMesh_t mesh, stormInt_t numVars, stormInt_t rank);
STORM_API stormArrayS_t SR_AllocS(
  stormMesh_t mesh, stormInt_t numVars, stormInt_t rank);
/// @}

/// @{
STORM_API stormArrayR_t stormAllocLikeR(stormArrayR_t array);
STORM_API stormArrayC_t stormAllocLikeC(stormArrayC_t array);
STORM_API stormArrayS_t stormAllocLikeS(stormArrayS_t array);
#define stormAllocLike(array) \
  _STORM_GENERIC_EXT_(array, stormAllocLike)(array)
/// @}

/// @{
STORM_API void stormFreeR(stormArrayR_t x);
STORM_API void stormFreeC(stormArrayC_t x);
STORM_API void stormFreeS(stormArrayS_t x);
#define stormFree(x) _STORM_GENERIC_EXT_(x, stormFree)(x)
/// @}

/// @{
#if STORM_C11_
STORM_INL void stormSwapP(void** pX, void** pY) {
  void* z = *pX; *pX = *pY, *pY = z;
}
#define stormSwap(x, y) stormSwapP((void**)&(x), (void**)&(y))
#elif STORM_CXX17_
#define stormSwap(x, y) std::swap(x, y)
#endif
/// @}

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

_STORM_OPAQUE_(stormIOList);

STORM_API stormIOList_t SR_IO_Begin();

STORM_API void SR_IO_Add(stormIOList_t IO,
  stormArrayR_t x, const char* name);

STORM_API void SR_IO_Flush(stormIOList_t IO,
  stormMesh_t mesh, const char* filename);

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

/// @{
STORM_API void stormFillR(stormMesh_t mesh,
                          stormArrayR_t x,
                          stormReal_t alpha,
                          _STORM_DEFAULT_(stormReal_t beta, 0.0));
STORM_API void stormFillC(stormMesh_t mesh,
                          stormArrayC_t x,
                          stormReal_t alpha,
                          _STORM_DEFAULT_(stormComplex_t beta, 0.0));
STORM_API void stormFillS(stormMesh_t mesh,
                          stormArrayS_t x,
                          stormReal_t alpha,
                          _STORM_DEFAULT_(stormSymbol_t beta, {}));
#define stormFill(mesh, x, alpha, ...) \
  _STORM_GENERIC_EXT_(x, stormFill)(mesh, x, alpha, ##__VA_ARGS__)
/// @}

/// @{
STORM_API void stormRandFillR(stormMesh_t mesh,
                              stormArrayR_t x,
                              _STORM_DEFAULT_(stormReal_t a, -1.0),
                              _STORM_DEFAULT_(stormReal_t b, +1.0));
STORM_API void stormRandFillC(stormMesh_t mesh,
                              stormArrayC_t x,
                              _STORM_DEFAULT_(stormReal_t a, -1.0),
                              _STORM_DEFAULT_(stormReal_t b, +1.0));
#define stormRandFill(mesh, x, ...) \
  _STORM_GENERIC_(x, stormRandFill)(mesh, x, ##__VA_ARGS__)
/// @}

/// @{
STORM_API void stormSetR(stormMesh_t mesh,
                         stormArrayR_t y,
                         stormArrayR_t x);
STORM_API void stormSetC(stormMesh_t mesh,
                         stormArrayC_t y,
                         stormArrayC_t x);
STORM_API void stormSetS(stormMesh_t mesh,
                         stormArrayS_t y,
                         stormArrayS_t x);
#define stormSet(mesh, y, x) \
  _STORM_GENERIC_EXT_(y, stormSet)(mesh, y, x)
/// @}

/// @{
STORM_API void stormScaleR(stormMesh_t mesh,
                           stormArrayR_t y,
                           stormArrayR_t x,
                           stormReal_t alpha);
STORM_API void stormScaleC(stormMesh_t mesh,
                           stormArrayC_t y,
                           stormArrayC_t x,
                           stormComplex_t alpha);
STORM_API void stormScaleS(stormMesh_t mesh,
                           stormArrayS_t y,
                           stormArrayS_t x,
                           stormSymbol_t alpha);
#define stormScale(mesh, y, x, alpha) \
  _STORM_GENERIC_EXT_(y, stormScale)(mesh, y, x, alpha)
/// @}

/// @{
STORM_API void stormAddR(stormMesh_t mesh,
                         stormArrayR_t z,
                         stormArrayR_t y,
                         stormArrayR_t x,
                         _STORM_DEFAULT_(stormReal_t alpha, 1.0),
                         _STORM_DEFAULT_(stormReal_t  beta, 1.0));
STORM_API void stormAddC(stormMesh_t mesh,
                         stormArrayC_t z,
                         stormArrayC_t y,
                         stormArrayC_t x,
                         _STORM_DEFAULT_(stormComplex_t alpha, 1.0),
                         _STORM_DEFAULT_(stormComplex_t  beta, 1.0));
STORM_API void stormAddS(stormMesh_t mesh,
                         stormArrayS_t z,
                         stormArrayS_t y,
                         stormArrayS_t x,
                         _STORM_DEFAULT_(stormSymbol_t alpha, {}),
                         _STORM_DEFAULT_(stormSymbol_t  beta, {}));
#define stormAdd(mesh, z, y, x, ...) \
  _STORM_GENERIC_EXT_(y, stormAdd)(mesh, z, y, x, ##__VA_ARGS__)
/// @}

/// @{
// No need to import these from Fortran.
STORM_INL void stormSubR(stormMesh_t mesh,
                         stormArrayR_t z,
                         stormArrayR_t y,
                         stormArrayR_t x,
                         _STORM_DEFAULT_(stormReal_t alpha, 1.0),
                         _STORM_DEFAULT_(stormReal_t  beta, 1.0)) {
  stormAddR(mesh, z, y, x, -alpha, beta);
}
STORM_INL void stormSubC(stormMesh_t mesh,
                         stormArrayC_t z,
                         stormArrayC_t y,
                         stormArrayC_t x,
                         _STORM_DEFAULT_(stormComplex_t alpha, 1.0),
                         _STORM_DEFAULT_(stormComplex_t  beta, 1.0)) {
  stormAddC(mesh, z, y, x, -alpha, beta);
}
STORM_INL void stormSubS(stormMesh_t mesh,
                         stormArrayS_t z,
                         stormArrayS_t y,
                         stormArrayS_t x,
                         _STORM_DEFAULT_(stormSymbol_t alpha, {}),
                         _STORM_DEFAULT_(stormSymbol_t  beta, {})) {
  assert(0); //stormAddS(mesh, z, y, x, -alpha, beta);
}
#define stormSub(mesh, z, y, x, ...) \
  _STORM_GENERIC_EXT_(y, stormSub)(mesh, z, y, x, ##__VA_ARGS__)
/// @}

/// @{
STORM_API void stormMulR(stormMesh_t mesh,
                         stormArrayR_t z,
                         stormArrayR_t y,
                         stormArrayR_t x);
STORM_API void stormMulC(stormMesh_t mesh,
                         stormArrayC_t z,
                         stormArrayC_t y,
                         stormArrayC_t x);
STORM_API void stormMulS(stormMesh_t mesh,
                         stormArrayS_t z,
                         stormArrayS_t y,
                         stormArrayS_t x);
#define stormMul(mesh, z, y, x) \
  _STORM_GENERIC_EXT_(y, stormMul)(mesh, z, y, x)
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
STORM_API stormReal_t SR_IntegrateR(stormMesh_t mesh,
                                    stormArrayR_t x,
                                    stormMapFuncR_t f,
                                    void* env);
STORM_API stormComplex_t SR_IntegrateC(stormMesh_t mesh,
                                       stormArrayC_t x,
                                       stormMapFuncC_t f,
                                       void* env);
#define SR_Integrate(mesh, x, f, env) \
  _STORM_GENERIC_(x, SR_Integrate)(mesh, x, f, env)
/// @}

/// @{
STORM_API void stormFuncProdR(stormMesh_t mesh,
                              stormArrayR_t y,
                              stormArrayR_t x,
                              stormMapFuncR_t f,
                              _STORM_DEFAULT_(void* env, STORM_NULL));
STORM_API void stormFuncProdC(stormMesh_t mesh,
                              stormArrayC_t y,
                              stormArrayC_t x,
                              stormMapFuncC_t f,
                              _STORM_DEFAULT_(void* env, STORM_NULL));
STORM_API void stormFuncProdS(stormMesh_t mesh,
                              stormArrayS_t y,
                              stormArrayS_t x,
                              stormMapFuncS_t f,
                              _STORM_DEFAULT_(void* env, STORM_NULL));
#define stormFuncProd(mesh, y, x, f, env) \
  _STORM_GENERIC_EXT_(y, stormFuncProd)(mesh, y, x, f, env)
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
STORM_API void stormSpFuncProdR(stormMesh_t mesh,
                                stormArrayR_t y,
                                stormArrayR_t x,
                                stormSpMapFuncR_t f,
                                _STORM_DEFAULT_(void* env, STORM_NULL));
STORM_API void stormSpFuncProdC(stormMesh_t mesh,
                                stormArrayC_t y,
                                stormArrayC_t x,
                                stormSpMapFuncC_t f,
                                _STORM_DEFAULT_(void* env, STORM_NULL));
STORM_API void stormSpFuncProdS(stormMesh_t mesh,
                                stormArrayS_t y,
                                stormArrayS_t x,
                                stormSpMapFuncS_t f,
                                _STORM_DEFAULT_(void* env, STORM_NULL));
#define stormSpFuncProd(mesh, y, x, f, env) \
  _STORM_GENERIC_EXT_(y, stormSpFuncProd)(mesh, y, x, f, env)
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
#define STORM_CG       "CG"
#define STORM_BiCGStab "BiCGStab"
#define STORM_MINRES   "MINRES"
#define STORM_GMRES    "GMRES"
#define STORM_TFQMR    "TFQMR"
#define STORM_LSQR     "LSQR"
#define STORM_LSMR     "LSMR"
/// @}

/// @{
#define STORM_NONE   ""
#define STORM_JACOBI "Jacobi"
#define STORM_LU_SGS "LU_SGS"
/// @}

/// @{
STORM_API void stormLinSolveR(stormMesh_t mesh,
                              stormString_t method,
                              stormString_t preMethod,
                              stormArrayR_t x,
                              stormArrayR_t b,
                              stormMatVecFuncR_t matVec,
                              _STORM_DEFAULT_(void* env, STORM_NULL),
                              _STORM_DEFAULT_(stormMatVecFuncR_t matVec_H, STORM_NULL),
                              _STORM_DEFAULT_(void* env_H, STORM_NULL));
STORM_API void stormLinSolveC(stormMesh_t mesh,
                              stormString_t method,
                              stormString_t preMethod,
                              stormArrayC_t x,
                              stormArrayC_t b,
                              stormMatVecFuncC_t matVec,
                              _STORM_DEFAULT_(void* env, STORM_NULL),
                              _STORM_DEFAULT_(stormMatVecFuncC_t matVec_H, STORM_NULL),
                              _STORM_DEFAULT_(void* env_H, STORM_NULL));
#define stormLinSolve(mesh, method, preMethod, x, b, matVec, ...) \
  _STORM_GENERIC_(x, stormLinSolve)(mesh, method, preMethod, x, b, matVec, ##__VA_ARGS__)
/// @}

#if STORM_CXX17_
/// @{
template<typename stormArrayT_t, typename stormMatVecFuncT_t>
STORM_INL void stormLinSolveT(stormMesh_t mesh,
                              stormString_t method,
                              stormString_t preMethod,
                              stormArrayT_t x,
                              stormArrayT_t b,
                              stormMatVecFuncT_t matVec) {
  stormLinSolve(mesh, method, preMethod, x, b,
    [](stormMesh_t mesh, stormArrayT_t Ax, stormArrayT_t x, void* env) {
      (*static_cast<stormMatVecFuncT_t*>(env))(mesh, Ax, x);
    }, &matVec);
} // stormLinSolveT
template<typename stormArrayT_t,
  typename stormMatVecFuncT1_t, typename stormMatVecFuncT2_t>
STORM_INL void stormLinSolveT(stormMesh_t mesh,
                              stormString_t method,
                              stormString_t preMethod,
                              stormArrayT_t x,
                              stormArrayT_t b,
                              stormMatVecFuncT1_t matVec,
                              stormMatVecFuncT2_t matVec_H) {
  stormLinSolve(mesh, method, preMethod, x, b,
    [](stormMesh_t mesh, stormArrayT_t Ax, stormArrayT_t x, void* env) {
      (*static_cast<stormMatVecFuncT1_t*>(env))(mesh, Ax, x);
    }, &matVec,
    [](stormMesh_t mesh, stormArrayT_t Ax, stormArrayT_t x, void* env) {
      (*static_cast<stormMatVecFuncT2_t*>(env))(mesh, Ax, x);
    }, &matVec_H);
} // stormLinSolveT
/// @}
#endif

/// @{
#define STORM_NEWTON "Newton"
#define STORM_JFNK   "JFNK"
/// @}

/// @{
STORM_API void stormNonlinSolveR(stormMesh_t mesh,
                                 stormString_t method,
                                 stormArrayR_t x,
                                 stormArrayR_t b,
                                 stormMatVecFuncR_t matVec,
                                 _STORM_DEFAULT_(void* env, STORM_NULL));
STORM_API void stormNonlinSolveC(stormMesh_t mesh,
                                 stormString_t method,
                                 stormArrayC_t x,
                                 stormArrayC_t b,
                                 stormMatVecFuncC_t matVec,
                                 _STORM_DEFAULT_(void* env, STORM_NULL));
#define stormNonlinSolve(mesh, method, x, b, matVec, ...) \
  _STORM_GENERIC_(x, stormNonlinSolve)(mesh, method, x, b, matVec, ##__VA_ARGS__)
/// @}

#if STORM_CXX17_
template<typename stormArrayT_t, typename stormMatVecFuncT_t>
STORM_INL void stormNonlinSolveT(stormMesh_t mesh,
                                 stormString_t method,
                                 stormArrayT_t x,
                                 stormArrayT_t b,
                                 stormMatVecFuncT_t matVec) {
  stormNonlinSolve(mesh, method, x, b,
    [](stormMesh_t mesh, stormArrayT_t Ax, stormArrayT_t x, void* env) {
      (*static_cast<stormMatVecFuncT_t*>(env))(mesh, Ax, x);
    }, &matVec);
} // stormNonlinSolveT
#endif

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

/// @{
STORM_API void stormApplyBCsR(stormMesh_t mesh,
                              stormArrayR_t u,
                              stormInt_t iBC,
                              stormReal_t alpha,
                              stormReal_t beta,
                              stormReal_t gamma);
STORM_API void stormApplyBCsC(stormMesh_t mesh,
                              stormArrayR_t u,
                              stormInt_t iBC,
                              stormComplex_t alpha,
                              stormComplex_t beta,
                              stormComplex_t gamma);
STORM_API void stormApplyBCsS(stormMesh_t mesh,
                              stormArrayR_t u,
                              stormInt_t iBC,
                              stormSymbol_t alpha,
                              stormSymbol_t beta,
                              stormSymbol_t gamma);
#define stormApplyBCs(mesh, u, iBC, ...) \
  _STORM_GENERIC_EXT_(u, stormApplyBCs)(mesh, u, iBC, ##__VA_ARGS__)
/// @}

#define SR_ALL 0
#define SR_PURE_DIRICHLET 1.0, 0.0, 0.0
#define SR_PURE_NEUMANN   0.0, 1.0, 0.0

#define SR_DIRICHLET(gamma) 1.0, 0.0, (gamma)
#define SR_NEUMANN(gamma)   0.0, 1.0, (gamma)

/// @{
STORM_API void stormApplyBCs_SlipWallR(stormMesh_t mesh,
                                       stormArrayR_t u,
                                       stormInt_t iBC);
#define stormApplyBCs_SlipWall stormApplyBCs_SlipWallR
/// @}

/// @{
STORM_API void stormApplyBCs_InOutLetR(stormMesh_t mesh,
                                       stormArrayR_t u,
                                       stormInt_t iBC);
#define stormApplyBCs_InOutLet stormApplyBCs_InOutLetR
/// @}

/// @{
STORM_API void stormGradientR(stormMesh_t mesh,
                              stormArrayR_t vVec,
                              stormReal_t lambda,
                              stormArrayR_t u);
STORM_API void stormGradientC(stormMesh_t mesh,
                              stormArrayC_t vVec,
                              stormComplex_t lambda,
                              stormArrayC_t u);
STORM_API void stormGradientS(stormMesh_t mesh,
                              stormArrayS_t vVec,
                              stormSymbol_t lambda,
                              stormArrayS_t u);
#define stormGradient(mesh, vVec, lambda, u) \
  _STORM_GENERIC_EXT_(vVec, stormGradient)(mesh, vVec, lambda, u)
/// @}

/// @{
STORM_API void stormDivergenceR(stormMesh_t mesh,
                                stormArrayR_t v,
                                stormReal_t lambda,
                                stormArrayR_t uVec);
STORM_API void stormDivergenceC(stormMesh_t mesh,
                                stormArrayC_t v,
                                stormComplex_t lambda,
                                stormArrayC_t uVec);
STORM_API void stormDivergenceS(stormMesh_t mesh,
                                stormArrayS_t v,
                                stormSymbol_t lambda,
                                stormArrayS_t uVec);
#define stormDivergence(mesh, v, lambda, uVec) \
  _STORM_GENERIC_EXT_(v, stormDivergence)(mesh, v, lambda, uVec)
/// @}

/// @{
STORM_API void stormRhieChowCorrection(stormMesh_t mesh,
                                       stormArrayR_t v,
                                       stormReal_t lambda,
                                       stormReal_t tau,
                                       stormArrayR_t p,
                                       stormArrayR_t rho);
/// @}

/// @{
STORM_API void stormConvectionR(stormMesh_t mesh,
                                stormArrayR_t v,
                                stormReal_t lambda,
                                stormArrayR_t u,
                                stormArrayR_t aVec);
STORM_API void stormConvectionC(stormMesh_t mesh,
                                stormArrayC_t v,
                                stormComplex_t lambda,
                                stormArrayC_t u,
                                stormArrayR_t aVec);
STORM_API void stormConvectionS(stormMesh_t mesh,
                                stormArrayS_t v,
                                stormSymbol_t lambda,
                                stormArrayS_t u,
                                stormArrayR_t aVec);
#define stormConvection(mesh, v, lambda, u, a) \
  _STORM_GENERIC_EXT_(v, stormConvection)(mesh, v, lambda, u, a)
/// @}

/// @{
STORM_API void stormDivGradR(stormMesh_t mesh,
                             stormArrayR_t v,
                             stormReal_t lambda,
                             stormArrayR_t u);
STORM_API void stormDivGradC(stormMesh_t mesh,
                             stormArrayC_t v,
                             stormComplex_t lambda,
                             stormArrayC_t u);
STORM_API void stormDivGradS(stormMesh_t mesh,
                             stormArrayS_t v,
                             stormSymbol_t lambda,
                             stormArrayS_t u);
#define stormDivGrad(mesh, v, lambda, u) \
  _STORM_GENERIC_EXT_(v, stormDivGrad)(mesh, v, lambda, u)
/// @}

/// @{
STORM_API void stormDivWGradR(stormMesh_t mesh,
                              stormArrayR_t v,
                              stormReal_t lambda,
                              stormArrayR_t k,
                              stormArrayR_t u);
STORM_API void stormDivWGradC(stormMesh_t mesh,
                              stormArrayC_t v,
                              stormComplex_t lambda,
                              stormArrayC_t k,
                              stormArrayC_t u);
STORM_API void stormDivWGradS(stormMesh_t mesh,
                              stormArrayS_t v,
                              stormSymbol_t lambda,
                              stormArrayR_t k,
                              stormArrayS_t u);
#define stormDivWGrad(mesh, v, lambda, k, u) \
  _STORM_GENERIC_EXT_(v, stormDivWGrad)(mesh, v, lambda, k, u)
/// @}

#endif // ifndef STORM_RULER_API_H
