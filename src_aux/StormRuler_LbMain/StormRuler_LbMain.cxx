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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <random>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//static double tau = 1.0e-2, Gamma = 1.0e-4, sigma = 1.0;

void InitF(int dim, const stormReal_t* r,
    int size, stormReal_t* f, const stormReal_t* _, void* __) {

  f[0] += 1;

} // Initial_Data

int main() {

  stormMesh_t mesh = SR_InitMesh();

  stormArray_t f, fHat, p, v, rho;
  f = stormAllocOnMesh(mesh, { 9 });
  fHat = stormAllocLike(f);

  p = SR_Alloc(mesh, 1, 0);
  v = SR_Alloc(mesh, 1, 1);
  rho = SR_Alloc(mesh, 1, 0);

  stormRandFill(mesh, f, 1.0-1.0e-2, 1.0+1.0e-2);
  //stormFill(mesh, f, 1.0);
  stormSpFuncProd(mesh, f, f, InitF, STORM_NULL);

  double total_time = 0.0;

  for (int time = 0; time <= 0*20+1*200000; ++time) {
    for (int frac = 0; time != 0 && frac < 500; ++frac) {

      struct timespec start, finish;

      clock_gettime(CLOCK_MONOTONIC, &start);

      ////////

      stormLbmMacroscopics(mesh, rho, v, {1.0}, f);
      stormLbmCollisionBGK(mesh, f, 0.6, rho, v);
      stormLbmStream(mesh, fHat, f);

      stormSwap(f, fHat);
      ////////

      clock_gettime(CLOCK_MONOTONIC, &finish);
      double elapsed = (finish.tv_sec - start.tv_sec);
      elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
      total_time += elapsed;

    }

    stormLbmMacroscopics(mesh, rho, v, {1.0}, f);

    char filename[256];
    printf("time = %f\n", total_time);
    sprintf(filename, "out/fld-%d.vtk", time);
    stormIOList_t io = SR_IO_Begin();
    SR_IO_Add(io, f, "distr");
    SR_IO_Add(io, v, "velocity");
    SR_IO_Add(io, p, "pressure");
    SR_IO_Add(io, rho, "density");
    SR_IO_Flush(io, mesh, filename);
  }

  return 0;
}
