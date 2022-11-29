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

namespace Storm
{

using namespace std::complex_literals; // NOLINT

// -----------------------------------------------------------------------------

TEST_CASE("Bittern/BoolMatrixReductions")
{
  Mat2x2<bool> mat_true{true, true, //
                        true, true};

  CHECK(all(mat_true));
  CHECK(any(mat_true));

  Mat2x2<bool> mat_false{false, false, //
                         false, false};

  CHECK_FALSE(all(mat_false));
  CHECK_FALSE(any(mat_false));

  Mat2x2<bool> mat_mixed{true, false, //
                         false, true};

  CHECK_FALSE(all(mat_mixed));
  CHECK(any(mat_mixed));
}

// -----------------------------------------------------------------------------

TEST_CASE("Bittern/RealMatrixReductions")
{
  Mat2x2<real_t> mat{+1.0, -2.0, //
                     +3.0, -4.0};

  CHECK_EQ(sum(mat), -2.0);
  CHECK_EQ(min_element(mat), -4.0);
  CHECK_EQ(max_element(mat), +3.0);
  CHECK_EQ(norm_1(mat), 10.0);
  CHECK_NEAR(norm_2(mat), 5.47723, EPS);
  CHECK_NEAR(norm_p(mat, 3), 4.64158, EPS);
  CHECK_EQ(norm_inf(mat), 4.0);
}

// -----------------------------------------------------------------------------

TEST_CASE("Bittern/ComplexMatrixReductions")
{
  Mat2x2<complex_t> mat{+1.0 - 2.0i, -2.0 + 3.0i, //
                        +3.0 - 4.0i, -4.0 + 5.0i};

  CHECK_EQ(sum(mat), -2.0 + 2.0i);
  CHECK_NEAR(norm_1(mat), 17.24474, EPS);
  CHECK_NEAR(norm_2(mat), 9.16515, EPS);
  CHECK_NEAR(norm_p(mat, 3), 7.63792, EPS);
  CHECK_NEAR(norm_inf(mat), 6.40312, EPS);
}

// -----------------------------------------------------------------------------

TEST_CASE("Bittern/DotProductReductions")
{
  Mat2x2<real_t> mat1{+1.0, -2.0, //
                      +3.0, -4.0};
  Mat2x2<real_t> mat2{+5.0, -6.0, //
                      +7.0, -8.0};
  Mat2x2<complex_t> mat3{+1.0 - 2.0i, -2.0 + 3.0i, //
                         +3.0 - 4.0i, -4.0 + 5.0i};

  CHECK_EQ(dot_product(mat1, mat2), 70.0);

  // We should get dot(mat2, mat3) == conj(dot(mat3, mat2)).
  CHECK_EQ(dot_product(mat2, mat3), 70.0 + 96.0i);
  CHECK_EQ(dot_product(mat3, mat2), 70.0 - 96.0i);

  // We should get a real number.
  CHECK_EQ(dot_product(mat3, mat3), 84.0);
}

// -----------------------------------------------------------------------------

TEST_CASE("Bittern/MatrixExpressionReductions")
{
  // Assuming the expressions are working.
  Mat2x2<real_t> mat1{+1.0, -2.0, //
                      +3.0, -4.0};
  Mat2x2<complex_t> mat2{+1.0 - 2.0i, -2.0 + 3.0i, //
                         +3.0 - 4.0i, -4.0 + 5.0i};
  const auto mat = 2.0 * mat1 + 3.0i * mat2;

  CHECK_EQ(sum(mat), -10.0 - 6.0i);
  CHECK_NEAR(norm_inf(mat), 25.94224, EPS);
}

} // namespace Storm
