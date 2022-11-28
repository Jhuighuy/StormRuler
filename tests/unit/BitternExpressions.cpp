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

#include "_UnitTests.hpp"

#include <Storm/Bittern/Mat.hpp>

#include <gtest/gtest.h>

#include <cmath>

namespace Storm
{

using namespace std::complex_literals;

TEST(BitternTestSuite, CompareMatrixExpressions)
{
  Mat2x2 mat1{+1.0, -2.0, //
              +3.0, -4.0};
  Mat2x2 mat2{+1.0, -3.0, //
              +2.0, -5.0};

  // Exact comparisons.
  EXPECT_TRUE(all(mat1 == mat1));
  EXPECT_TRUE(any(mat1 != mat2));
  EXPECT_TRUE(all(mat1 >= mat2));
  EXPECT_TRUE(any(mat1 > mat2));
  EXPECT_FALSE(all(mat1 > mat2));
  EXPECT_TRUE(any(mat1 <= mat2));
  EXPECT_FALSE(any(mat1 < mat2));

  Mat2x2<real_t> approx_mat1{+1.001, -2.001, //
                             +3.001, -4.001};

  // Approximate comparisons.
  EXPECT_TRUE(all(approx_equal(mat1, approx_mat1, 0.002)));
  EXPECT_FALSE(any(approx_equal(mat1, approx_mat1, 0.0001)));
}

TEST(BitternTestSuite, LogicalMatrixExpressions)
{
  Mat2x2 mat1{true, true, //
              true, true};
  Mat2x2 mat2{false, false, //
              false, false};
  Mat2x2 mat3{true, false, //
              false, true};

  EXPECT_FALSE(any(mat1 && mat2));
  EXPECT_TRUE(all(mat1 || mat2));
  EXPECT_TRUE(any(mat1 && mat3));
}

#if 0
TEST(BitternTestSuite, ArithmeticMatrixExpressions) {
  Mat2x2<real_t> mat1{+1.0, -2.0, //
                      +3.0, -4.0};
  Mat2x2<real_t> mat2{+5.0, -6.0, //
                      +7.0, -8.0};
  Mat2x2<real_t> mat3{+3.0, +9.0, //
                      +2.0, -1.0};
  EXPECT_TRUE(all((mat1 + 2.0 * mat2) * mat3 == //
                  mat2 * 2.0 * mat3 + mat3 * mat1));
  EXPECT_TRUE(all((mat1 / mat2 + mat2 / mat3) * mat3 * mat2 == //
                  mat2 * mat2 + mat1 * mat2));
}
#endif

/// @todo SignMatrixExpressions

// Assuming the expressions are working,
// testing the power functions to output the correct result.
TEST(BitternTestSuite, PowerMatrixExpressions)
{
  Mat2x2 mat{2.0, 3.0, //
             3.0, 2.0};

  Mat2x2 sqrt_mat{sqrt2, sqrt3, //
                  sqrt3, sqrt2};

  EXPECT_TRUE(all(approx_equal(sqrt(mat), sqrt_mat, EPS)));
  EXPECT_TRUE(all(approx_equal(pow(mat, 0.5), sqrt_mat, EPS)));
  EXPECT_TRUE(all(approx_equal(pow(sqrt_mat, 2), mat, EPS)));

  Mat2x2 cbrt_mat{cbrt(2.0), cbrt(3.0), //
                  cbrt(3.0), cbrt(2.0)};

  EXPECT_TRUE(all(approx_equal(cbrt(mat), cbrt_mat, EPS)));
  EXPECT_TRUE(all(approx_equal(pow(mat, 1.0 / 3.0), cbrt_mat, EPS)));
  EXPECT_TRUE(all(approx_equal(pow(cbrt_mat, 3), mat, EPS)));

  Mat2x2 two_pow_mat{4.0, 8.0, //
                     8.0, 4.0};

  EXPECT_TRUE(all(approx_equal(pow(2, mat), two_pow_mat, EPS)));
}

/// @todo ExponentialMatrixExpressions

// Assuming the expressions are working,
// testing the trigonometric functions to output the correct result.
TEST(BitternTestSuite, TrigonometricMatrixExpressions)
{
  Mat2x2 mat{+pi / 6.0, -pi / 4.0, //
             -pi / 4.0, +pi / 6.0};

  Mat2x2 sin_mat{+0.5, -1.0 / sqrt2, //
                 -1.0 / sqrt2, +0.5};

  EXPECT_TRUE(all(approx_equal(sin(mat), sin_mat, EPS)));
  EXPECT_TRUE(all(approx_equal(asin(sin_mat), mat, EPS)));

  Mat2x2 cos_mat{sqrt3 / 2.0, 1.0 / sqrt2, //
                 1.0 / sqrt2, sqrt3 / 2.0};

  EXPECT_TRUE(all(approx_equal(cos(mat), cos_mat, EPS)));
  EXPECT_TRUE(all(approx_equal(acos(cos_mat), abs(mat), EPS)));

  Mat2x2 tan_mat{1.0 / sqrt3, -1.0, //
                 -1.0, 1.0 / sqrt3};

  EXPECT_TRUE(all(approx_equal(tan(mat), tan_mat, EPS)));
  EXPECT_TRUE(all(approx_equal(atan(tan_mat), mat, EPS)));

  /// @todo atan2 tests!
}

// Assuming the expressions are working,
// testing the hyperbolic functions to output the correct result.
TEST(BitternTestSuite, HyperbolicMatrixExpressions)
{
  Mat2x2 mat{+pi / 6.0i, -pi / 4.0i, //
             -pi / 4.0i, +pi / 6.0i};

  Mat2x2 sinh_mat{-0.5i, +1.0i / sqrt2, //
                  +1.0i / sqrt2, -0.5i};

  EXPECT_TRUE(all(approx_equal(sinh(mat), sinh_mat, EPS)));
  EXPECT_TRUE(all(approx_equal(asinh(sinh_mat), mat, EPS)));

  Mat2x2 cosh_mat{sqrt3 / 2.0, 1.0 / sqrt2, //
                  1.0 / sqrt2, sqrt3 / 2.0};

  EXPECT_TRUE(all(approx_equal(cosh(mat), cosh_mat, EPS)));
  EXPECT_TRUE(all(approx_equal(acosh(cosh_mat), 1.0i * abs(mat), EPS)));

  Mat2x2 tanh_mat{-1.0i / sqrt3, 1.0i, //
                  1.0i, 1.0i / -sqrt3};

  EXPECT_TRUE(all(approx_equal(tanh(mat), tanh_mat, EPS)));
  EXPECT_TRUE(all(approx_equal(atanh(tanh_mat), mat, EPS)));
}

TEST(BitternTestSuite, CrossProductExpressions)
{
  Vec3D vec1{1.0, 0.0, 0.0}, vec2{0.0, 1.0, 0.0}, vec3{0.0, 0.0, 1.0};

  EXPECT_TRUE(all(cross_product(vec1, vec2) == vec3));
}

} // namespace Storm
