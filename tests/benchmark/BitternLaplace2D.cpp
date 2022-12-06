/// Copyright (C) 2020-2023 Oleg Butakov
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

#include "../unit/_UnitTests.hpp" // TODO: share CHECK_NEAR somehow.
#include "./_Benchmarks.hpp"

#include <Storm/Base.hpp> // TODO: some more relevent include file
#include <Storm/Utils/Meta.hpp>

#include <doctest/doctest.h>
#include <nanobench.h>

#include <cmath>
#include <numbers>
#include <type_traits>

// -----------------------------------------------------------------------------

// Based on:
// https://github.com/romeric/expression_templates_benchmark

namespace Storm::Benchmarks
{

// TODO: our implementation!

} // namespace Storm::Benchmarks

// -----------------------------------------------------------------------------

#if STORM_BENCH_ARMADILLO_ENABLED

#include <armadillo>

namespace Storm::Benchmarks
{

template<class T, size_t N>
using arma_vec = typename arma::Col<T>::template fixed<N>;
template<class T, size_t N>
using arma_mat = typename arma::Mat<T>::template fixed<N, N>;

template<class T, size_t N, size_t NumIterations>
class Laplace2D_Armadillo
{
public:

  T operator()() const
  {
    static constexpr T pi = std::numbers::pi_v<T>;

    const arma_vec<T, N> x = arma::linspace(T{0.0}, pi, N);
    arma_mat<T, N> u;
    u.fill(T{0.0});
    u.col(0) = arma::sin(x);
    u.col(N - 1) = arma::sin(x) * std::exp(-pi);

    T error;
    for (size_t iteration = 0; iteration < NumIterations; iteration++) {
      const arma_mat<T, N> u_old = u;
      u(arma::span(1, N - 2), arma::span(1, N - 2)) =
          (T{4.0} * (u_old(arma::span(0, N - 3), arma::span(1, N - 2)) +
                     u_old(arma::span(2, N - 1), arma::span(1, N - 2)) +
                     u_old(arma::span(1, N - 2), arma::span(0, N - 3)) +
                     u_old(arma::span(1, N - 2), arma::span(2, N - 1))) +
           T{1.0} * (u_old(arma::span(0, N - 3), arma::span(0, N - 3)) +
                     u_old(arma::span(0, N - 3), arma::span(2, N - 1)) +
                     u_old(arma::span(2, N - 1), arma::span(0, N - 3)) +
                     u_old(arma::span(2, N - 1), arma::span(2, N - 1)))) /
          T{20.0};

      // `arma::norm` is notoriously slow.
      const T square_norm = arma::accu(arma::pow(u - u_old, 2));
      error = std::sqrt(square_norm);
      ankerl::nanobench::doNotOptimizeAway(error);
    }

    return error;
  }

}; // class Laplace2D_Armadillo

} // namespace Storm::Benchmarks

#endif // if STORM_BENCH_ARMADILLO_ENABLED

// -----------------------------------------------------------------------------

#if STORM_BENCH_BLAZE_ENABLED

#include <blaze/Blaze.h>

namespace Storm::Benchmarks
{

template<class T, size_t N, size_t NumIterations>
class Laplace2D_Blaze
{
public:

  T operator()() const
  {
    static constexpr T pi = std::numbers::pi_v<T>;

    const blaze::StaticVector<T, N> x = blaze::linspace(N, T{0.0}, pi);
    blaze::StaticMatrix<T, N, N> u(T{0.0});
    blaze::column(u, 0) = sin(x);
    blaze::column(u, N - 1) = sin(x) * std::exp(-pi);

    T error;
    for (size_t iteration = 0; iteration < NumIterations; iteration++) {
      const blaze::StaticMatrix<T, N, N> u_old = u;
      blaze::submatrix(u, 1, 1, N - 2, N - 2) =
          (T{4.0} * (blaze::submatrix(u_old, 0, 1, N - 2, N - 2) +
                     blaze::submatrix(u_old, 2, 1, N - 2, N - 2) +
                     blaze::submatrix(u_old, 1, 0, N - 2, N - 2) +
                     blaze::submatrix(u_old, 1, 2, N - 2, N - 2)) +
           T{1.0} * (blaze::submatrix(u_old, 0, 0, N - 2, N - 2) +
                     blaze::submatrix(u_old, 0, 2, N - 2, N - 2) +
                     blaze::submatrix(u_old, 2, 0, N - 2, N - 2) +
                     blaze::submatrix(u_old, 2, 2, N - 2, N - 2))) /
          T{20.0};

      error = norm(u - u_old);
      ankerl::nanobench::doNotOptimizeAway(error);
    }

    return error;
  }

}; // class Laplace2D_Blaze

} // namespace Storm::Benchmarks

#endif // if STORM_BENCH_BLAZE_ENABLED

// -----------------------------------------------------------------------------

#if STORM_BENCH_EIGEN_ENABLED

#define EIGEN_STACK_ALLOCATION_LIMIT (1ull << 32ull)
#include <eigen3/Eigen/Core>

namespace Storm::Benchmarks
{

template<class T, size_t N, size_t NumIterations>
class Laplace2D_Eigen
{
public:

  T operator()() const
  {
    static constexpr T pi = std::numbers::pi_v<T>;

    const auto x = Eigen::Matrix<T, N, 1>::LinSpaced(N, T{0.0}, pi);
    Eigen::Matrix<T, N, N> u;
    u.setZero();
    u.col(0) = x.array().sin();
    u.col(N - 1) = x.array().sin() * std::exp(-pi);

    T error;
    for (size_t iteration = 0; iteration < NumIterations; iteration++) {
      const Eigen::Matrix<T, N, N> u_old = u;
      u.block(1, 1, N - 2, N - 2) = //
          (T{4.0} * (u_old.block(0, 1, N - 2, N - 2) +
                     u_old.block(2, 1, N - 1, N - 2) +
                     u_old.block(1, 0, N - 2, N - 3) +
                     u_old.block(1, 2, N - 2, N - 1)) +
           T{1.0} * (u_old.block(0, 0, N - 3, N - 3) +
                     u_old.block(0, 2, N - 3, N - 1) +
                     u_old.block(2, 0, N - 1, N - 3) +
                     u_old.block(2, 2, N - 1, N - 1))) /
          T{20.0};

      error = (u - u_old).norm();
      ankerl::nanobench::doNotOptimizeAway(error);
    }

    return error;
  }

}; // class Laplace2D_Eigen

} // namespace Storm::Benchmarks

#endif // if STORM_BENCH_EIGEN_ENABLED

// -----------------------------------------------------------------------------

#if STORM_BENCH_XTENSOR_ENABLED

#include <xtensor/xfixed.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xview.hpp>

namespace Storm::Benchmarks
{

template<class T, size_t N, size_t NumIterations>
class Laplace2D_XTensor
{
public:

  T operator()() const
  {
    static constexpr T pi = std::numbers::pi_v<T>;

    const xt::xtensor_fixed<T, xt::xshape<N>> x = xt::linspace(T{0.0}, pi, N);
    xt::xtensor_fixed<T, xt::xshape<N, N>> u = xt::zeros<T>({N, N});
    xt::col(u, 0) = xt::sin(x);
    xt::col(u, N - 1) = xt::sin(x) * std::exp(-pi);

    T error;
    for (size_t iteration = 0; iteration < NumIterations; iteration++) {
      const xt::xtensor_fixed<T, xt::xshape<N, N>> u_old = u;
      xt::view(u, xt::range(1, N - 1), xt::range(1, N - 1)) =
          (T{4.0} *
               (xt::view(u_old, xt::range(0, N - 2), xt::range(1, N - 1)) +
                xt::view(u_old, xt::range(2, N - 0), xt::range(1, N - 1)) +
                xt::view(u_old, xt::range(1, N - 1), xt::range(0, N - 2)) +
                xt::view(u_old, xt::range(1, N - 1), xt::range(2, N - 0))) +
           T{1.0} *
               (xt::view(u_old, xt::range(0, N - 2), xt::range(0, N - 2)) +
                xt::view(u_old, xt::range(0, N - 2), xt::range(2, N - 0)) +
                xt::view(u_old, xt::range(2, N - 0), xt::range(0, N - 2)) +
                xt::view(u_old, xt::range(2, N - 0), xt::range(2, N - 0)))) /
          T{20.0};

      const T square_norm = xt::sum(xt::pow(u - u_old, 2))();
      error = std::sqrt(square_norm);
      ankerl::nanobench::doNotOptimizeAway(error);
    }

    return error;
  }

}; // class Laplace2D_XTensor

} // namespace Storm::Benchmarks

#endif // if STORM_BENCH_XTENSOR_ENABLED

// -----------------------------------------------------------------------------

namespace Storm::Benchmarks
{

template<class T, size_t N, size_t NumIterations, //
         template<class, size_t, size_t> class Laplace2D>
T run_laplace_2D(const char* library_name)
{
  const auto benchmark_name =
      fmt::format("BitternBenchmarks/Laplace2D({}, T={}, N={}, I={})", //
                  library_name, meta::type_name_v<T>, N, NumIterations);

  const auto laplace_2D = Laplace2D<T, N, NumIterations>{};

  T error;
  ankerl::nanobench::Bench().run(benchmark_name, [&] {
    error = laplace_2D();
    ankerl::nanobench::doNotOptimizeAway(error);
  });

  return error;
}

TEST_CASE("BitternBenchmarks/Laplace2D")
{
  auto run_subcase = [&]<class T, size_t N, size_t NumIterations>(
                         meta::type<T>, size_t_constant<N>,
                         size_t_constant<NumIterations>, //
                         T expected_error, T eps) {
    SUBCASE("Bittern")
    {
      // TODO: our implementation!
    }

#if STORM_BENCH_BLAZE_ENABLED
    SUBCASE("Armadillo")
    {
      const double error =
          run_laplace_2D<T, N, NumIterations, Laplace2D_Armadillo>("Armadillo");
      CHECK_NEAR(error, expected_error, eps);
    }
#endif

#if STORM_BENCH_BLAZE_ENABLED
    SUBCASE("Blaze")
    {
      const double error =
          run_laplace_2D<T, N, NumIterations, Laplace2D_Blaze>("Blaze");
      CHECK_NEAR(error, expected_error, eps);
    }
#endif

#if STORM_BENCH_EIGEN_ENABLED
    SUBCASE("Eigen")
    {
      const double error =
          run_laplace_2D<T, N, NumIterations, Laplace2D_Eigen>("Eigen");
      CHECK_NEAR(error, expected_error, eps);
    }
#endif

#if STORM_BENCH_XTENSOR_ENABLED
    SUBCASE("XTensor")
    {
      const double error =
          run_laplace_2D<T, N, NumIterations, Laplace2D_XTensor>("XTensor");
      CHECK_NEAR(error, expected_error, eps);
    }
#endif
  };

  static constexpr size_t NumIterations = 1000;

  SUBCASE("T=double")
  {
    using T = double;
    static constexpr T eps = 1.0e-6;

    SUBCASE("N=100")
    {
      static constexpr size_t N = 100;
      static constexpr T expected_error = 0.0069143;
      run_subcase(meta::type_v<T>, size_t_constant<N>{},
                  size_t_constant<NumIterations>{}, expected_error, eps);
    }

    SUBCASE("N=150")
    {
      static constexpr size_t N = 150;
      static constexpr T expected_error = 0.00994008;
      run_subcase(meta::type_v<T>, size_t_constant<N>{},
                  size_t_constant<NumIterations>{}, expected_error, eps);
    }

    SUBCASE("N=200")
    {
      static constexpr size_t N = 200;
      static constexpr T expected_error = 0.0121789;
      run_subcase(meta::type_v<T>, size_t_constant<N>{},
                  size_t_constant<NumIterations>{}, expected_error, eps);
    }
  }
}

} // namespace Storm::Benchmarks

// -----------------------------------------------------------------------------
