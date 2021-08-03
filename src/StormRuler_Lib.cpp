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

#include <cmath>
#include <iostream>
#include <typeinfo>

template <typename T>
constexpr auto type_name(T&&) noexcept {
  std::string name;
#ifdef __clang__
  name = __PRETTY_FUNCTION__;
#elif defined(__GNUC__)
  name = __PRETTY_FUNCTION__;
#elif defined(_MSC_VER)
  name = __FUNCSIG__;
#endif
  return name;
}

#include "StormRuler_Lib.hpp"
using namespace StormRuler;

double dWdC(double c) {
  return -c*(1.0 - c*c);
}

double dt = (M_PI/50)*(M_PI/50), Gamma = 0.01;

void CahnHilliard_Step(tField<0> c, tField<1> v, 
                       tField<0> c_hat, tField<0> w_hat) {
  tField<0> rhs = AllocateField<0>();
  rhs << MAP(&dWdC, c);
  rhs << c - dt*CONV(c, v) + dt*DIVGRAD(rhs);

  SOLVE_BiCGSTAB([&](tField<0> in, tField<0> out) {
    out << 0.0;
    out += DIVGRAD(in);
    out << in + Gamma*dt*DIVGRAD(out);
  }, c_hat, rhs);

  w_hat << MAP(&dWdC, c_hat);
  w_hat -= Gamma*DIVGRAD(c_hat);
}

double rho = 1.0, nu = 0.1, beta = 0.0;

void NavierStokes_Step(tField<0> p, tField<1> v, tField<0> c, tField<0> w,
                       tField<0> p_hat, tField<1> v_hat) {
  tField<1> f = AllocateField<1>();
  // We want to write this:
  //f << dt/rho*c*GRAD(w);
  f << 0.0;
  f -= 1.0/rho*GRAD(w);
  BLAS_Mul(f, c, f);

  v_hat << v + dt*f - dt*CONV(v, v) - dt*beta/rho*GRAD(p) + dt*nu*DIVGRAD(v);

  tField<0> rhs = AllocateField<0>();
  rhs << 0.0;
  rhs += rho/dt*DIV(v_hat);

  p_hat << 0.0;
  SOLVE_BiCGSTAB([&](tField<0> in, tField<0> out) {
    out << 0.0;
    out += DIVGRAD(in);
  }, p_hat, rhs);

  v_hat -= dt/rho*GRAD(p_hat);
  p_hat << p_hat + beta*p;
}

int main() {
  std::cout << "Hello from C++" << std::endl;

  Lib_InitializeMesh();

  {
    auto c = AllocateField<0>(), 
      c_hat = AllocateField<0>(), w_hat = AllocateField<0>();

    BLAS_SFuncProd(c, c, [&](double* coords, double* in, double* out) {
      double x = coords[0], y = coords[1];
      out[0] = (abs(x - M_PI) < 1.0 && abs(y - M_PI) < 1.0) ? -1.0 : 1.0;
    });
    auto p = AllocateField<0>(), p_hat = AllocateField<0>();
    p << 1.0;
    auto v = AllocateField<1>(), v_hat = AllocateField<1>();
    v << 0.0;

    _Lib_IO_Begin();
    _Lib_IO_Add(v);
    _Lib_IO_Add(p);
    _Lib_IO_Add(c);
    _Lib_IO_End();

    for (int L = 1; L <= 200; ++L) {
      for (int M = 0; M < 10; ++M) {
        CahnHilliard_Step(c, v, c_hat, w_hat);
        NavierStokes_Step(p, v, c_hat, w_hat, p_hat, v_hat);

        std::swap(c, c_hat);
        std::swap(p, p_hat);
        std::swap(v, v_hat);
      }

      _Lib_IO_Begin();
      _Lib_IO_Add(v);
      _Lib_IO_Add(p);
      _Lib_IO_Add(c);
      _Lib_IO_End();
    }
  }
/*
  {
    double dt, gamma;
    tField<0> rhs = AllocateField<0>(), c = AllocateField<0>(), u = AllocateField<0>(); 
    
    rhs << MAP(&dWdC, c);
    rhs += dt*DIVGRAD(rhs);
    
    u << 0.0;
    u += DIVGRAD(u);
    u << c + gamma*dt*DIVGRAD(u);

    BLAS_SFuncProd(rhs, rhs, [&](double* x, double* in, double* out) {
      std::cout << x[0] << ' ' << x[1] << ' ' << in[0] << std::endl;
    });
  }
*/

/*
  {
    tField<1> u = AllocateField<1>();
    tField<1> v = AllocateField<1>();
    tField<1> w = AllocateField<1>();

    w << 1.0*u + 2.0*v + 3.0*w;

    tField<0> q = AllocateField<0>();
    w << 1.0*u + 2.0*v - 3.0*15.0*GRAD(q) - 3.0*w;
    //w << v - GRAD(q);

    BLAS_SFuncProd(w, w, [&](double* x, double* in, double* out) {
      std::cout << x[0] << ' ' << x[1] << ' ' << in[0] << ' ' << in[1] << std::endl;
    });
  }
*/

/*
  {
    auto u = AllocateField<0>();
    auto v = AllocateField<0>();
    auto w = AllocateField<0>();
    auto q = AllocateField<1>();

    w << 1.0*u + 2.0*v - 3.0*15.0*DIV(q) - 4.0*w;
      BLAS_SFuncProd(w, w, [&](double* x, double* in, double* out) {
      std::cout << x[0] << ' ' << x[1] << ' ' << in[0] << std::endl;
    });
  }
*/

/*
  {
    auto u = AllocateField<0>();
    auto v = AllocateField<0>();
    auto w = AllocateField<0>();
    auto q = AllocateField<0>();

    w << 1.0*u + 2.0*v - 3.0*15.0*DIVGRAD(q) - 4.0*w;
      BLAS_SFuncProd(w, w, [&](double* x, double* in, double* out) {
      std::cout << x[0] << ' ' << x[1] << ' ' << in[0] << std::endl;
    });
  }
*/

/*
  BLAS_Add(u, u, v, 2.0, 1.0);

  BLAS_SFuncProd(v, u, [&](double* x, double* in, double* out) {
    std::cout << x[0] << ' ' << x[1] << ' ' << in[0] << ' ' << in[1] << std::endl;
    out[0] = 228.0;
    out[1] = 282.0;
  });

  BLAS_SFuncProd(u, v, [&](double* x, double* in, double* out) {
    std::cout << x[0] << ' ' << x[1] << ' ' << in[0] << ' ' << in[1] << std::endl;
    out[0] = 228.0;
  });
*/

  std::cout << "Dosvidanya from C++" << std::endl;

  return 0;
}
