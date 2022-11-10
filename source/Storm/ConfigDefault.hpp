/// Copyright (C) 2022 Oleg Butakov
///
/// Permission is hereby granted, free of charge, to any person obtaining a copy
/// of this software and associated documentation files (the "Software"), to
/// deal in the Software without restriction, including without limitation the
/// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
/// sell copies of the Software, and to permit persons to whom the Software is
/// furnished to do so, subject to the following conditions:
///
/// The above copyright notice and this permission notice shall be included in
/// all copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
/// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
/// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
/// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
/// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
/// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
/// IN THE SOFTWARE.

#pragma once

// Threading libraries configutation.
#ifndef STORM_OpenMP_ENABLED
#define STORM_OpenMP_ENABLED 0
#endif
#ifndef STORM_TBB_ENABLED
#define STORM_TBB_ENABLED 0
#endif

// Algebra libraries configutation.
#ifndef STORM_MKL_ENABLED
#define STORM_MKL_ENABLED 0
#endif
#ifndef STORM_BLAS_ENABLED
#define STORM_BLAS_ENABLED 0
#endif
#ifndef STORM_LAPACK_ENABLED
#define STORM_LAPACK_ENABLED 0
#endif

// Metis configutation.
#ifndef STORM_METIS_ENABLED
#define STORM_METIS_ENABLED 0
#endif

// Zlib configutation.
#ifndef STORM_ZLIB_ENABLED
#define STORM_ZLIB_ENABLED 0
#endif

// OpenGL configutation.
#ifndef STORM_OpenGL_ENABLED
#define STORM_OpenGL_ENABLED 0
#endif
#ifndef STORM_GLEW_ENABLED
#define STORM_GLEW_ENABLED 0
#endif
#ifndef STORM_GLFW_ENABLED
#define STORM_GLFW_ENABLED 0
#endif
#ifndef STORM_GLM_ENABLED
#define STORM_GLM_ENABLED 0
#endif
