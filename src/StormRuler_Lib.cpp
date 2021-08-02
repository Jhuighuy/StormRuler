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

#include <iostream>

#include "StormRuler_Lib.hpp"
#include "StormRuler_Lib.inc"

int main() {
  using namespace StormRuler;
  std::cout << "Hello from C++" << std::endl;

  Lib_InitializeMesh();

  tField<1> u = AllocateField<1>();
  tField<1> v = AllocateField<1>();

  BLAS_Add(u, u, v, 2.0, 1.0);

  BLAS_SFuncProd(v, u, [&](double* x, double* in, double* out) {
    std::cout << x[0] << ' ' << x[1] << ' ' << in[0] << ' ' << in[1] << std::endl;
  });

  //FDM_Gradient(v, 1.0, u, 'c');

  std::cout << "Dosvidanya from C++" << std::endl;

  return 0;
}
