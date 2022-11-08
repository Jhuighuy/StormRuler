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

#include "_Main.hpp"

#include <Storm/Bittern/Mat.hpp>

#include <gtest/gtest.h>

namespace Storm {

using namespace std::complex_literals;

TEST(BitternTestSuite, BoolMatrixReductions) {
  Mat2x2<bool> mat_true{true, true, //
                        true, true};

  EXPECT_TRUE(all(mat_true));
  EXPECT_TRUE(any(mat_true));

  Mat2x2<bool> mat_false{false, false, //
                         false, false};

  EXPECT_FALSE(all(mat_false));
  EXPECT_FALSE(any(mat_false));

  Mat2x2<bool> mat_mixed{true, false, //
                         false, true};

  EXPECT_FALSE(all(mat_mixed));
  EXPECT_TRUE(any(mat_mixed));
}

TEST(BitternTestSuite, RealMatrixReductions) {
  Mat2x2<real_t> mat{+1.0, -2.0, //
                     +3.0, -4.0};

  EXPECT_EQ(sum(mat), -2.0);
  EXPECT_EQ(min_element(mat), -4.0);
  EXPECT_EQ(max_element(mat), +3.0);
  EXPECT_EQ(norm_1(mat), 10.0);
  EXPECT_NEAR(norm_2(mat), 5.47723, EPS);
  EXPECT_NEAR(norm_p(mat, 3), 4.64158, EPS);
  EXPECT_EQ(norm_inf(mat), 4.0);
}

TEST(BitternTestSuite, ComplexMatrixReductions) {
  Mat2x2<complex_t> mat{+1.0 - 2.0i, -2.0 + 3.0i, //
                        +3.0 - 4.0i, -4.0 + 5.0i};

  EXPECT_EQ(sum(mat), -2.0 + 2.0i);
  EXPECT_NEAR(norm_1(mat), 17.24474, EPS);
  EXPECT_NEAR(norm_2(mat), 9.16515, EPS);
  EXPECT_NEAR(norm_p(mat, 3), 7.63792, EPS);
  EXPECT_NEAR(norm_inf(mat), 6.40312, EPS);
}

TEST(BitternTestSuite, DotProductReductions) {
  Mat2x2<real_t> mat1{+1.0, -2.0, //
                      +3.0, -4.0};
  Mat2x2<real_t> mat2{+5.0, -6.0, //
                      +7.0, -8.0};
  Mat2x2<complex_t> mat3{+1.0 - 2.0i, -2.0 + 3.0i, //
                         +3.0 - 4.0i, -4.0 + 5.0i};

  EXPECT_EQ(dot_product(mat1, mat2), 70.0);

  // We should get dot(mat2, mat3) == conj(dot(mat3, mat2)).
  EXPECT_EQ(dot_product(mat2, mat3), 70.0 + 96.0i);
  EXPECT_EQ(dot_product(mat3, mat2), 70.0 - 96.0i);

  // We should get a real number.
  EXPECT_EQ(dot_product(mat3, mat3), 84.0);
}

TEST(BitternTestSuite, MatrixExpressionReductions) {
  // Assuming the expressions are working.
  Mat2x2<real_t> mat1{+1.0, -2.0, //
                      +3.0, -4.0};
  Mat2x2<complex_t> mat2{+1.0 - 2.0i, -2.0 + 3.0i, //
                         +3.0 - 4.0i, -4.0 + 5.0i};
  const auto mat = 2.0 * mat1 + 3.0i * mat2;

  EXPECT_EQ(sum(mat), -10.0 - 6.0i);
  EXPECT_NEAR(norm_inf(mat), 25.94224, EPS);
}

#if 0
TEST(BitternTestSuite, BlockMatrixReduction) {
  Mat2x2<real_t> A{1.0, 2.0, //
                   3.0, 4.0};
  Mat2x2<real_t> B{5.0, 6.0, //
                   7.0, 8.0};
  Mat2x2<Mat2x2<real_t>> mat{A, B, //
                             B, A};
  EXPECT_EQ(sum(mat), 1.0);
}
#endif

} // namespace Storm
