// Copyright (C) 2020-2023 Oleg Butakov
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "../unit/_UnitTests.hpp"
#include "./_Benchmarks.hpp"

#include <Storm/Base.hpp> /// @todo some more relevent include file
#include <Storm/Utils/Meta.hpp>

#include <doctest/doctest.h>
#include <nanobench.h>

#include <cmath>
#include <numbers>
#include <type_traits>

// Based on:
// https://shorturl.at/efBKL
// https://github.com/romeric/expression_templates_benchmark

// -----------------------------------------------------------------------------

namespace Storm::Benchmarks {

/// @todo Our implementation!

} // namespace Storm::Benchmarks

// -----------------------------------------------------------------------------

#if STORM_BENCH_ARMADILLO_ENABLED
#  ifdef NDEBUG
#    define ARMA_NO_DEBUG 1
#  else
#    define ARMA_NO_DEBUG 0
#  endif
#  define ARMA_MAT_PREALLOC 100000ull
#  include <armadillo>

namespace Storm::Benchmarks {

// Laplace2D benchmark implementation (Armadillo).
template<class T, size_t N, size_t NumIterations>
class Laplace2D_Armadillo final {
private:

  using arma_vec = typename arma::Col<T>::template fixed<N>;
  using arma_mat = typename arma::Mat<T>::template fixed<N, N>;

public:

  [[nodiscard]] T operator()() const noexcept {
    static constexpr T pi = std::numbers::pi_v<T>;

    const arma_vec x = arma::linspace(T{0.0}, pi, N);
    arma_mat u;
    u.fill(T{0.0});
    u.col(0) = arma::sin(x);
    u.col(N - 1) = arma::sin(x) * std::exp(-pi);

    T error;
    for (size_t iteration = 0; iteration < NumIterations; iteration++) {
      const arma_mat u_old = u;
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

      // Note: `arma::norm` computes operator norm, not Frobenius.
      error = std::sqrt(arma::accu(arma::pow(u - u_old, 2)));
      nanobench::doNotOptimizeAway(error);
    }

    return error;
  }

}; // class Laplace2D_Armadillo

} // namespace Storm::Benchmarks

#endif // if STORM_BENCH_ARMADILLO_ENABLED

// -----------------------------------------------------------------------------

#if STORM_BENCH_BLAZE_ENABLED
#  include <blaze/Blaze.h>

namespace Storm::Benchmarks {

// Laplace2D benchmark implementation (Blaze).
template<class T, size_t N, size_t NumIterations>
class Laplace2D_Blaze final {
public:

  [[nodiscard]] T operator()() const noexcept {
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

      error = blaze::norm(u - u_old);
      nanobench::doNotOptimizeAway(error);
    }

    return error;
  }

}; // class Laplace2D_Blaze

} // namespace Storm::Benchmarks

#endif // if STORM_BENCH_BLAZE_ENABLED

// -----------------------------------------------------------------------------

#if STORM_BENCH_EIGEN_ENABLED
#  define EIGEN_STACK_ALLOCATION_LIMIT (1ull << 32ull)
#  include <eigen3/Eigen/Core>

namespace Storm::Benchmarks {

// Laplace2D benchmark implementation (Eigen).
template<class T, size_t N, size_t NumIterations>
class Laplace2D_Eigen final {
public:

  [[nodiscard]] T operator()() const noexcept {
    static constexpr T pi = std::numbers::pi_v<T>;

    const auto x = Eigen::Matrix<T, N, 1>::LinSpaced(N, T{0.0}, pi);
    Eigen::Matrix<T, N, N> u;
    u.setZero();
    u.col(0) = x.array().sin();
    u.col(N - 1) = x.array().sin() * std::exp(-pi);

    T error;
    for (size_t iteration = 0; iteration < NumIterations; iteration++) {
      const Eigen::Matrix<T, N, N> u_old = u;
      u.block(1, 1, N - 2, N - 2) =
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
      nanobench::doNotOptimizeAway(error);
    }

    return error;
  }

}; // class Laplace2D_Eigen

} // namespace Storm::Benchmarks

#endif // if STORM_BENCH_EIGEN_ENABLED

// -----------------------------------------------------------------------------

#if STORM_BENCH_XTENSOR_ENABLED
#  include <xtensor/xfixed.hpp>
#  include <xtensor/xnorm.hpp>
#  include <xtensor/xview.hpp>

namespace Storm::Benchmarks {

// Laplace2D benchmark implementation (XTensor).
template<class T, size_t N, size_t NumIterations>
class Laplace2D_XTensor final {
public:

  [[nodiscard]] T operator()() const noexcept {
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

      error = xt::norm_l2(u - u_old)();
      nanobench::doNotOptimizeAway(error);
    }

    return error;
  }

}; // class Laplace2D_XTensor

} // namespace Storm::Benchmarks

#endif // if STORM_BENCH_XTENSOR_ENABLED

// -----------------------------------------------------------------------------

#if STORM_BENCH_NUMPY_ENABLED
#  include <pybind11/embed.h>
#  if STORM_COMPILER_GCC_
#    pragma GCC visibility push(hidden)
#  endif

namespace Storm::Benchmarks {

/// @todo Move this line to a better place.
namespace py = pybind11;

// Laplace2D benchmark implementation (NumPy).
template<class T, size_t N, size_t NumIterations>
class Laplace2D_NumPy final {
private:

  py::object laplace2D_;

public:

  Laplace2D_NumPy() {
    py::dict globals = py::globals();
    globals["N"] = N, globals["NumIterations"] = NumIterations;
    py::exec(R"""(
      from math import pi, sqrt
      import numpy

      def laplace2D():
        x = numpy.linspace(0.0, pi, N)
        u = numpy.zeros((N, N))
        u[:, 0] = numpy.sin(x)
        u[:, N - 1] = numpy.sin(x) * numpy.exp(-pi)

        for _ in range(NumIterations):
          u_old = numpy.copy(u)
          u[1:-1, 1:-1] = (
              4.0 * (u_old[0:-2, 1:-1] + u_old[2:  , 1:-1]  +
                     u_old[1:-1, 0:-2] + u_old[1:-1, 2:  ]) +
              1.0 * (u_old[0:-2, 0:-2] + u_old[0:-2, 2:  ]  +
                     u_old[2:  , 0:-2] + u_old[2:  , 2:  ])
          ) / 20.0

          error = sqrt(numpy.sum((u - u_old) ** 2))

        return error)""",
             /*globals=*/globals, /*locals=*/globals);
    laplace2D_ = globals["laplace2D"];
  }

  [[nodiscard]] T operator()() const {
    return laplace2D_().cast<T>();
  }

}; // class Laplace2D_NumPy

} // namespace Storm::Benchmarks

#  if STORM_COMPILER_GCC_
#    pragma GCC visibility pop
#  endif
#endif // if STORM_BENCH_NUMPY_ENABLED

// -----------------------------------------------------------------------------

namespace Storm::Benchmarks {

TEST_CASE("Bittern/Laplace2D") {
  const auto run_benchmarks = []<class T, size_t N, size_t NumIterations>(
                                  T expected_error, T tolerance,
                                  std::index_sequence<N, NumIterations>) {
    const auto run_impl = [&](const auto& impl, const char* library_name) {
      SUBCASE(library_name) {
        const auto benchmark_name =
            STORM_FORMAT("Laplace2D({}, T={}, N={}, I={})", library_name,
                         meta::type_name_v<T>, N, NumIterations);
        T error;
        nanobench::Bench{}.run(benchmark_name, [&] { error = impl(); });
        CHECK_NEAR(error, expected_error, tolerance);
      }
    };

    /// @todo Our implementation!
    static_cast<void>(run_impl);

#if STORM_BENCH_ARMADILLO_ENABLED
    run_impl(Laplace2D_Armadillo<T, N, NumIterations>{}, "Armadillo");
#endif

#if STORM_BENCH_BLAZE_ENABLED
    run_impl(Laplace2D_Blaze<T, N, NumIterations>{}, "Blaze");
#endif

#if STORM_BENCH_EIGEN_ENABLED
    run_impl(Laplace2D_Eigen<T, N, NumIterations>{}, "Eigen");
#endif

#if STORM_BENCH_XTENSOR_ENABLED
    run_impl(Laplace2D_XTensor<T, N, NumIterations>{}, "XTensor");
#endif

#if STORM_BENCH_NUMPY_ENABLED
    run_impl(Laplace2D_NumPy<T, N, NumIterations>{}, "NumPy");
#endif
  };

  static constexpr size_t NumIterations = 1000;

  SUBCASE("T=double") {
    using T = double;
    static constexpr T tolerance = 1.0e-4;

    SUBCASE("N=100") {
      static constexpr size_t N = 100;
      static constexpr T expected_error = 0.0069143;
      run_benchmarks(expected_error, tolerance,
                     std::index_sequence<N, NumIterations>{});
    }

    SUBCASE("N=150") {
      static constexpr size_t N = 150;
      static constexpr T expected_error = 0.00994008;
      run_benchmarks(expected_error, tolerance,
                     std::index_sequence<N, NumIterations>{});
    }

    SUBCASE("N=200") {
      static constexpr size_t N = 200;
      static constexpr T expected_error = 0.0121789;
      run_benchmarks(expected_error, tolerance,
                     std::index_sequence<N, NumIterations>{});
    }
  }
}

} // namespace Storm::Benchmarks

// -----------------------------------------------------------------------------
