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

#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

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
#define STORM_CXX14_ 1
#else
#define STORM_CXX14_ 0
#endif

#if STORM_C11_ && STORM_CXX14_
#error StormRuler API: both C11 or C++17 targets found.
#endif

#if !STORM_C11_ && !STORM_CXX14_
#error StormRuler API: neither C11 nor C++17 targets found.
#endif

#if STORM_CXX14_
#include <vector>
#endif

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

#define _STORM_CAT2_(x, y) x##y
#define _STORM_CAT_(x, y) _STORM_CAT2_(x, y)

#define _STORM_STR2_(x) #x
#define _STORM_STR_(x) _STORM_STR2_(x)

#define _STORM_OPAQUE_(type) \
  struct type;               \
  typedef struct type* type##_t

// Enable default arguments for C++.
#if STORM_C11_
#define _STORM_DEFAULT_(decl, ...) decl
#elif STORM_CXX14_
#define _STORM_DEFAULT_(decl, ...) decl = __VA_ARGS__
#endif

#if STORM_C11_
#define STORM_API extern
#define STORM_INL static
#elif STORM_CXX14_
#define STORM_API extern "C"
#define STORM_INL inline
#endif

#if STORM_C11_
#define STORM_NULL NULL
#elif STORM_CXX14_
#define STORM_NULL nullptr
#endif

// Basic types.
typedef int stormInt_t;
typedef size_t stormSize_t;
typedef ptrdiff_t stormPtrDiff_t;
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

STORM_API stormArray_t SR_Alloc(stormMesh_t mesh, stormInt_t numVars,
                                stormInt_t rank);

STORM_API stormArray_t stormAllocOnMesh(stormMesh_t mesh, stormSize_t size,
                                        const stormSize_t* shape);
#if STORM_CXX14_
STORM_INL stormArray_t stormAllocOnMesh(stormMesh_t mesh,
                                        const std::vector<stormSize_t>& shape) {
  return stormAllocOnMesh(mesh, shape.size(), shape.data());
} // stormAllocOnMesh
#endif

STORM_API stormArray_t stormAllocLike(stormArray_t array);

STORM_API void stormFree(stormArray_t x);

STORM_API void stormArrayUnwrap(stormMesh_t mesh, stormArray_t x,
                                stormReal_t** data, stormSize_t* size);

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

_STORM_OPAQUE_(stormIOList);

STORM_API stormIOList_t SR_IO_Begin();

STORM_API void SR_IO_Add(stormIOList_t IO, stormArray_t x, const char* name);

STORM_API void SR_IO_Flush(stormIOList_t IO, stormMesh_t mesh,
                           const char* filename);

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

STORM_API void stormMul(stormMesh_t mesh, stormArray_t z, stormArray_t y,
                        stormArray_t x);

typedef void (*stormSpMapFunc_t)(stormSize_t dim, const stormReal_t* r,
                                 stormSize_t size, stormReal_t* y,
                                 const stormReal_t* x, stormOpaque_t env);

STORM_API void stormSpFuncProd(stormMesh_t mesh, stormArray_t y, stormArray_t x,
                               stormSpMapFunc_t f,
                               _STORM_DEFAULT_(stormOpaque_t env, STORM_NULL));

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

STORM_API void stormApplyBCs(stormMesh_t mesh, stormArray_t u, stormInt_t iBC,
                             stormReal_t alpha, stormReal_t beta,
                             stormReal_t gamma);

#define SR_ALL 0
#define SR_PURE_DIRICHLET 1.0, 0.0, 0.0
#define SR_PURE_NEUMANN 0.0, 1.0, 0.0

#define SR_DIRICHLET(gamma) 1.0, 0.0, (gamma)
#define SR_NEUMANN(gamma) 0.0, 1.0, (gamma)

STORM_API void stormApplyBCs_SlipWall(stormMesh_t mesh, stormArray_t u,
                                      stormInt_t iBC);

STORM_API void stormApplyBCs_CosWall(stormMesh_t mesh, stormArray_t phiBar,
                                     stormArray_t phi, stormReal_t a,
                                     stormInt_t iBC);

STORM_API void stormApplyBCs_InOutLet(stormMesh_t mesh, stormArray_t u,
                                      stormInt_t iBC);

STORM_API void stormGradient(stormMesh_t mesh, stormArray_t vVec,
                             stormReal_t lambda, stormArray_t u);

STORM_API void stormDivergence(stormMesh_t mesh, stormArray_t v,
                               stormReal_t lambda, stormArray_t uVec);

STORM_API void stormDivGrad(stormMesh_t mesh, stormArray_t v,
                            stormReal_t lambda, stormArray_t u);

STORM_API void stormDivWGrad(stormMesh_t mesh, stormArray_t v,
                             stormReal_t lambda, stormArray_t k,
                             stormArray_t u);

STORM_API void stormConvection(stormMesh_t mesh, stormArray_t v,
                               stormReal_t lambda, stormArray_t u,
                               stormArray_t aVec);

STORM_API void stormRhieChowCorrection(stormMesh_t mesh, stormArray_t v,
                                       stormReal_t lambda, stormReal_t tau,
                                       stormArray_t p, stormArray_t rho);

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

#endif // ifndef STORM_RULER_API_H
