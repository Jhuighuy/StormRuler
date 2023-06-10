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

#include "./_UnitTests.hpp"

#include <Storm/Bittern/Mat.hpp>

#include <doctest/doctest.h>

namespace Storm::UnitTests {

// -----------------------------------------------------------------------------

TEST_CASE("Bittern/BoolMatrixReductions") {
  SUBCASE("true-mat") {
    Mat2x2<bool> true_mat{true, true, //
                          true, true};

    CHECK(all(true_mat));
    CHECK(any(true_mat));
  }

  SUBCASE("false-mat") {
    Mat2x2<bool> false_mat{false, false, //
                           false, false};

    CHECK_FALSE(all(false_mat));
    CHECK_FALSE(any(false_mat));
  }

  SUBCASE("mixed-mat") {
    Mat2x2<bool> mixed_mat{true, false, //
                           false, true};

    CHECK_FALSE(all(mixed_mat));
    CHECK(any(mixed_mat));
  }
}

// -----------------------------------------------------------------------------

TEST_CASE("Bittern/RealMatrixReductions") {
  Mat2x2<real_t> mat{+1.0_dp, -2.0_dp, //
                     +3.0_dp, -4.0_dp};

  SUBCASE("simple-reductions") {
    CHECK_EQ(sum(mat), -2.0_dp);
    CHECK_EQ(min_element(mat), -4.0_dp);
    CHECK_EQ(max_element(mat), +3.0_dp);
  }

  SUBCASE("norms") {
    const real_t eps = 1.0e-5_dp;
    CHECK_EQ(norm_1(mat), 10.0_dp);
    CHECK_NEAR(norm_2(mat), 5.47723_dp, eps);
    CHECK_NEAR(norm_p(mat, 3), 4.64158_dp, eps);
    CHECK_EQ(norm_inf(mat), 4.0_dp);
  }
}

// -----------------------------------------------------------------------------

TEST_CASE("Bittern/ComplexMatrixReductions") {
  Mat2x2<complex_t> mat{+1.0_dp - i * 2.0_dp, -2.0_dp + i * 3.0_dp, //
                        +3.0_dp - i * 4.0_dp, -4.0_dp + i * 5.0_dp};

  SUBCASE("simple-reductions") {
    CHECK_EQ(sum(mat), -2.0_dp + i * 2.0_dp);
  }

  SUBCASE("norms") {
    const real_t eps = 1.0e-5_dp;
    CHECK_NEAR(norm_1(mat), 17.24474_dp, eps);
    CHECK_NEAR(norm_2(mat), 9.16515_dp, eps);
    CHECK_NEAR(norm_p(mat, 3), 7.63792_dp, eps);
    CHECK_NEAR(norm_inf(mat), 6.40312_dp, eps);
  }
}

// -----------------------------------------------------------------------------

TEST_CASE("Bittern/DotProductReductions") {
  Mat2x2<real_t> mat1{+1.0_dp, -2.0_dp, //
                      +3.0_dp, -4.0_dp};
  Mat2x2<real_t> mat2{+5.0_dp, -6.0_dp, //
                      +7.0_dp, -8.0_dp};
  Mat2x2<complex_t> mat3{+1.0_dp - i * 2.0_dp, -2.0_dp + i * 3.0_dp, //
                         +3.0_dp - i * 4.0_dp, -4.0_dp + i * 5.0_dp};

  CHECK_EQ(dot_product(mat1, mat2), 70.0_dp);

  // We should get dot(mat2, mat3) == conj(dot(mat3, mat2)).
  CHECK_EQ(dot_product(mat2, mat3), 70.0_dp + i * 96.0_dp);
  CHECK_EQ(dot_product(mat3, mat2), 70.0_dp - i * 96.0_dp);

  // We should get a real number.
  CHECK_EQ(dot_product(mat3, mat3), 84.0_dp);
}

// -----------------------------------------------------------------------------

// Assuming the expressions are working.
TEST_CASE("Bittern/MatrixExpressionReductions") {
  Mat2x2<real_t> mat1{+1.0_dp, -2.0_dp, //
                      +3.0_dp, -4.0_dp};
  Mat2x2<complex_t> mat2{+1.0_dp - i * 2.0_dp, -2.0_dp + i * 3.0_dp, //
                         +3.0_dp - i * 4.0_dp, -4.0_dp + i * 5.0_dp};
  const auto mat_expr = 2.0_dp * mat1 + i * 3.0_dp * mat2;

  const real_t eps = 1.0e-5_dp;
  CHECK_EQ(sum(mat_expr), -10.0_dp - i * 6.0_dp);
  CHECK_NEAR(norm_inf(mat_expr), 25.94224_dp, eps);
}

} // namespace Storm::UnitTests
