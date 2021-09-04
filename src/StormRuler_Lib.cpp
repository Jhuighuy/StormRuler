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

#include <cmath>
#include <iostream>
#include <typeinfo>

#include <array>
#include <vector>
#include <fstream>
#include <algorithm>

template <typename T>
auto type_name(T&&) noexcept {
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

//
// dx = 0.5*10^-9.
// sigma = from ML
// lambda = from ML.
// 

#include "StormRuler_Lib.hpp"
using namespace StormRuler;

#if 1

double lerp(double x, const std::vector<double>& ivals) {
  x *= (ivals.size() - 1);
  double x0 = std::floor(x), x1 = std::ceil(x);
  size_t i0 = x0, i1 = x1;
  if (i0 == i1) {
    return ivals[i0];
  }
  return (ivals[i0]*(x1 - x) + ivals[i1]*(x - x0))/(x1 - x0);
}

#if 1
double dWdC(double c) {
  return -c*(1.0 - c*c);
}
#else
std::vector<double> ivals
  {0.0,1.5906018654950425, 2.6188272215722872, 3.1827287996745524, 
   3.3685758977822746, 3.2493980037514945, 2.8881855486953696, 
   2.3403368866745748, 1.6556122214974556, 0.87976958254585613, 
   0.056010611492066174, -0.77366104347571574, -1.5670760775132615, 
   -2.2804359862391315, -2.8666406117308516, -3.2732681306319122, 
   -3.4400177858916012, -3.2953421222288592, -2.7518601347379161, 
   -1.6999219735344779,0.0 };
double dWdC(double c) {
  c = (c + 1.0)/2.0;
  return lerp(c, ivals)/4.0;
}
#endif

std::array<std::vector<double>, 5> npart_vs_phi{
  std::vector<double>
  {6710.9163537992381, 7016.0968832367416, 7267.4708047934419, 
   7467.6794466816755, 7619.4891958037397, 7725.8293934008398, 
   7789.6834271923535, 7814.0084041099281, 7801.6775529130664, 
   7755.4403442952816, 7677.8962798720177, 7571.4791735516073, 
   7438.4494621150980, 7280.8926552982284, 7100.7225141285435, 
   6899.6879735582252, 6679.3832544587422, 6441.2610998312812, 
   6186.6497111619237, 5916.7749020434767, 5633.1430155648895,},

  std::vector<double>
  {596.38096304074429, 726.26627185612199, 864.14571296483416, 
   1007.9730504827404, 1155.8052126901985, 1305.807682866187, 
   1456.2872423453753, 1605.7091814146702, 1752.7029587286388, 
   1896.0596994865075, 2034.7243272716003, 2167.7845968026477, 
   2294.4588613351207, 2414.0840808604371, 2526.1053790409492, 
   2630.0684142624577, 2725.6159965604561, 2812.4908503791771, 
   2890.5473612722299, 2959.7768720218924, 3020.3738137534788},

  std::vector<double>
  {58.046923614810076, 80.674177544752823, 108.10339135108698, 
   140.43119759948448, 177.65657240620095, 219.68756157301900, 
   266.35923728186754, 317.45121147741793, 372.70363011265459, 
   431.83110357622098, 494.53442140680886, 560.51017725587315, 
   629.45862285040687, 701.09020729514123, 775.13137062862745, 
   851.33028587631338, 929.46342264034206, 1009.3440990136767, 
   1090.8346865090803, 1173.8649915140859, 1258.4101959299621,},

  std::vector<double>
  {6.3071248575747028, 9.9736394172950753, 14.99068515221488, 
   21.581746187020453, 29.952862353277645, 40.28614968366314, 
   52.737496196885765, 67.436462031346466, 84.487814962479831, 
   103.97414496285877, 125.9590758070204, 150.49069161199168, 
   177.60489636391833, 207.32851225863752, 239.68198605401204, 
   274.68160790654986, 312.34114516262161, 352.67274397761679, 
   395.68682286951486, 441.39041435567265, 489.75411630323185,},

  std::vector<double>
  {0.54311042015136468, 1.0384935915258402, 1.8464351149816978, 
   3.0861019688128151, 4.8929367053538808, 7.4137944581949444, 
   10.801811533224829, 15.211423237191939, 20.793937650212985, 
   27.693875544909236, 36.046103645549799, 45.973635705055223, 
   57.585855956924433, 70.976820268110856, 86.223192382459786, 
   103.38124588382844, 122.48216978367047, 143.52459683432221, 
   166.46274264594953, 191.18764610702181, 217.4972982277452,},
};
double nPart(size_t i, double c) {
  c = (c + 1.0)/2.0;
  return lerp(c, npart_vs_phi[i]);
}

double mol_mass[] {
  0.016043, 0.058124, 0.1002, 0.14229, 0.1984
};

double dt = (M_PI/50)*(M_PI/50), Gamma = 0.01;

void CahnHilliard_Step(tField<0> c, tField<1> v, 
                       tField<0> c_hat, tField<0> w_hat) {
  tField<0> rhs = AllocField<0>();
  rhs << MAP(&dWdC, c);
  rhs << c - dt*CONV(c, v) + dt*DIVGRAD(rhs);

  c_hat << c;
  SOLVE_BiCGSTAB([&](tField<0> in, tField<0> out) {
    out << 0.0;
    out += DIVGRAD(in);
    out << in + Gamma*dt*DIVGRAD(out);
  }, c_hat, rhs);

  w_hat << MAP(&dWdC, c_hat);
  w_hat -= Gamma*DIVGRAD(c_hat);
}

double rho = 1.0, nu = 0.1, beta = 0.0;

void NavierStokes_Step(tField<0> p, tField<1> v, 
                       tField<0> c, tField<0> w,
                       tField<0> p_hat, tField<1> v_hat) {
  v_hat << v - dt*CONV(v, v) - dt*beta/rho*GRAD(p) - dt/rho*(c*GRAD(w)) + dt/rho*nu*DIVGRAD(v);

  tField<0> rhs = AllocField<0>();
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

double rho0 = 1.0, rho1 = 2.0;
 
#if 1
void NavierStokes_VaDensity_Step(tField<0> p, tField<1> v, 
                                 tField<0> c, tField<0> w,
                                 tField<0> p_hat, tField<1> v_hat) {
  tField<0> rho = AllocField<0>();
  rho << MAP([](double c) {
    c = std::min(1.0, std::max(c, -1.0)); 
    return rho0*(1.0-c)/2 + rho1*(1.0+c)/2; 
  }, c);
  tField<0> rho_inv = AllocField<0>();
  tField<0> one = AllocField<0>();
  one << 1.0;
  BLAS_Mul(rho_inv, rho, one, -1);

  //v_hat << v - dt*CONV(v, v) - dt*beta*GRAD(p) - dt*(c*GRAD(w)) + dt*nu*DIVGRAD(v);
  v_hat << v - dt*CONV(v, v) - 
    dt*beta*(GRAD(p)/rho) - dt*((c*GRAD(w))/rho) + dt*nu*(DIVGRAD(v)/rho);

  tField<0> rhs = AllocField<0>();
  rhs << 0.0;
  rhs += (1.0/dt*DIV(v_hat));

  p_hat << p;
  SOLVE_BiCGSTAB([&](tField<0> in, tField<0> out) {
    out << 0.0;
    FDM_DivWGrad(out, 1.0, rho_inv, in);
    //out += DIVGRAD(in);
  }, p_hat, rhs, 'L');

  v_hat << v_hat - dt*(GRAD(p_hat)/rho);
  p_hat << p_hat + beta*p;
}
#else
#define NavierStokes_VaDensity_Step NavierStokes_Step
#endif

double mu = 0.1e-3;

template<size_t M>
void CustomNavierStokes_Step(tField<0> p, tField<1> v, 
                             tField<0> c, tField<0> w, 
                             std::array<tField<0>, M> nPart_hat, tField<0> rho_hat,
                             tField<0> p_hat, tField<1> v_hat) {
  /*rho_hat << MAP([](double c) {
    c = std::min(1.0, std::max(c, -1.0)); 
    return rho0*(1.0-c)/2 + rho1*(1.0+c)/2; 
  }, c);*/

  rho_hat << 0.0;
  for (size_t i = 0; i < M; ++i) {
    nPart_hat[i] << MAP([&](double c){ return nPart(i, c); }, c);
    rho_hat << rho_hat + mol_mass[i]*nPart_hat[i];
  }
  rho_hat << 0.0*rho_hat + 0.1*rho_hat;
  tField<0> rho_inv = AllocField<0>();
  tField<0> one = AllocField<0>();
  one << 1.0;
  BLAS_Mul(rho_inv, rho_hat, one, -1);

  v_hat << v - dt*CONV(v, v) - 
    dt*beta*(GRAD(p)/rho_hat) - dt*((c*GRAD(w))/rho_hat) + dt*nu*(DIVGRAD(v)/rho_hat);

  tField<0> rhs = AllocField<0>();
  rhs << 0.0;
  rhs += (1.0/dt*DIV(v_hat));
  //rho_hat << rhs;
  //return;

  p_hat << 0.0;
  SOLVE_BiCGSTAB([&](tField<0> in, tField<0> out) {
    out << 0.0;
    FDM_DivWGrad(out, 1.0, rho_inv, in);
    //out += DIVGRAD(in);
  }, p_hat, rhs);

  v_hat << v_hat - dt*(GRAD(p_hat)/rho_hat);
  p_hat << p_hat + beta*p;
}

///////////////////////////////////////////////

int main() {
  std::cout << "Hello from C++ " << sizeof(int) << std::endl;

  Lib_InitializeMesh();

  {
    auto c = AllocField<0>(), 
      c_hat = AllocField<0>(), w_hat = AllocField<0>();

    BLAS_SFuncProd(c, c, [&](double* coords, double* in, double* out) {
      double x = coords[0] - M_PI, y = coords[1] - M_PI;
      out[0] = 1.0;
      //if (hypot(x, y) < 1.0) out[0] = -1.0;
      // Square.
      //if ((fabs(x) < 1.0) && (fabs(y) < 1.0)) out[0] = -1.0;
      // Cross
      if ((fabs(x) < 0.2) && (fabs(y) < 2.0)) out[0] = -1.0;
      if ((fabs(y) < 0.2) && (fabs(x) < 2.0)) out[0] = -1.0;
      // Uncomment the following lines to get some very interesting shape:
      //if ((fabs(x - 1.8) < 0.2) && (fabs(y + 1.0) < 1.0)) out[0] = -1.0;
      //if ((fabs(x + 1.8) < 0.2) && (fabs(y - 1.0) < 1.0)) out[0] = -1.0;
      //if ((fabs(y + 1.8) < 0.2) && (fabs(x + 1.0) < 1.0)) out[0] = -1.0;
      //if ((fabs(y - 1.8) < 0.2) && (fabs(x - 1.0) < 1.0)) out[0] = -1.0;

      //double x = coords[0], y = coords[1];
      //out[0] = (abs(x - M_PI) < 1.0 && abs(y - M_PI) < 1.0) ? -1.0 : 1.0;
    });
    auto p = AllocField<0>(), p_hat = AllocField<0>();
    p << 1.0;
    auto v = AllocField<1>(), v_hat = AllocField<1>();
    v << 0.0;
    std::array<tField<0>, 5> nPart_hat = {
      AllocField<0>(), AllocField<0>(), AllocField<0>(), AllocField<0>(), AllocField<0>()
    };
    auto rho_hat = AllocField<0>();

    _Lib_IO_Begin();
    _Lib_IO_Add(v, "velocity");
    _Lib_IO_Add(p, "pressure");
    _Lib_IO_Add(c, "phase");
    _Lib_IO_Add(rho_hat, "rho");
    _Lib_IO_End();

    for (int L = 1; L <= 2000; ++L) {
      for (int M = 0; M < 1; ++M) {
        CahnHilliard_Step(c, v, c_hat, w_hat);
        NavierStokes_VaDensity_Step(p, v, c_hat, w_hat, p_hat, v_hat);
        //CustomNavierStokes_Step(p, v, c_hat, w_hat, nPart_hat, rho_hat, p_hat, v_hat);

        std::swap(c, c_hat);
        std::swap(p, p_hat);
        std::swap(v, v_hat);
      }

      _Lib_IO_Begin();
      _Lib_IO_Add(v, "velocity");
      _Lib_IO_Add(p, "pressure");
      _Lib_IO_Add(c, "phase");
      _Lib_IO_Add(rho_hat, "rho");
      _Lib_IO_End();
    }
  }

  std::cout << "Dosvidanya from C++" << std::endl;

  return 0;
}

#endif
