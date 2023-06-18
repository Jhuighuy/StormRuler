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

/// @todo cast expressions!

// -----------------------------------------------------------------------------

TEST_CASE("Bittern/CompareMatrixExpressions") {
  const Mat2x2<real_t> mat1{+1.0_dp, -2.0_dp, //
                            +3.0_dp, -4.0_dp};

  SUBCASE("exact-comparisons") {
    const Mat2x2<real_t> mat2{+1.0_dp, -3.0_dp, //
                              +2.0_dp, -5.0_dp};

    CHECK(all(mat1 == mat1));
    CHECK(any(mat1 != mat2));
    CHECK(all(mat1 >= mat2));
    CHECK(any(mat1 > mat2));
    CHECK_FALSE(all(mat1 > mat2));
    CHECK(any(mat1 <= mat2));
    CHECK_FALSE(any(mat1 < mat2));
  }

  SUBCASE("approximate-comparisons") {
    const Mat2x2<real_t> approx_mat1{+1.001_dp, -2.001_dp, //
                                     +3.001_dp, -4.001_dp};

    CHECK(all(approx_equal(mat1, approx_mat1, 0.002_dp)));
    CHECK_FALSE(any(approx_equal(mat1, approx_mat1, 0.0001_dp)));
  }
}

/// @todo Min/Max!

// -----------------------------------------------------------------------------

TEST_CASE("Bittern/LogicalMatrixExpressions") {
  Mat2x2<bool> mat1{true, true, //
                    true, true};
  Mat2x2<bool> mat2{false, false, //
                    false, false};
  Mat2x2<bool> mat3{true, false, //
                    false, true};

  SUBCASE("logical-not") {
    CHECK_FALSE(any(!mat1));
    CHECK(all(!mat2));
  }

  SUBCASE("logical-and") {
    CHECK_FALSE(any(matrix_and(mat1, mat2)));
    CHECK(any(matrix_and(mat1, mat3)));
    /// @todo Check the multiargument versions.
  }

  SUBCASE("logical-or") {
    CHECK(all(matrix_or(mat1, mat2)));
    /// @todo Check the multiargument versions.
  }
}

/// @todo Merge!

// -----------------------------------------------------------------------------

TEST_CASE("Bittern/ArithmeticMatrixExpressions") {
  const Mat2x2<real_t> mat1{+1.0_dp, -2.0_dp, //
                            +3.0_dp, -4.0_dp};

  SUBCASE("real-expressions") {
    const Mat2x2<real_t> mat2{+5.0_dp, -6.0_dp, //
                              +7.0_dp, -8.0_dp};
    const Mat2x2<real_t> mat3{+3.0_dp, +9.0_dp, //
                              +2.0_dp, -1.0_dp};

    SUBCASE("expr-1") {
      const Mat2x2<real_t> result{+21.0_dp, -152.0_dp, //
                                  +53.0_dp, -74.0_dp};

      CHECK(all(mat1 + 10.0 * (mat2 - mat3) == result));
    }

    SUBCASE("expr-2") {
      const Mat2x2<real_t> result{+297.5_dp, +894.0_dp, //
                                  +189.5_dp, -116.0_dp};

      CHECK(all(-mat1 * mat2 * 0.5_dp + mat3 / 0.01_dp == result));
    }

    SUBCASE("expr-3") {
      const Mat2x2<real_t> result{+58.0_dp, +174.0_dp, //
                                  +16.0_dp, +8.0_dp};

      CHECK(all(24.0_dp / +mat1 + (18.0_dp * mat3 - 4.0_dp * mat2) == result));
    }

    SUBCASE("expr-4") {
      const Mat2x2<real_t> result{-4.0_dp, +8.0_dp, //
                                  +13.0_dp, +88.0_dp};
      CHECK(all(2.0_dp * ((9.0_dp * mat1 / mat3) - mat2) == result));
    }
  }
}

TEST_CASE("Bittern/NormalizeMatrixExpressions") {
  SUBCASE("normalize") {
    const Mat2x2<real_t> mat{0.0_dp, 1.0_dp, //
                             2.0_dp, 2.0_dp};

    CHECK(all(approx_equal(normalize(mat), mat / 3.0_dp)));
  }

  SUBCASE("normalize-zero") {
    const Mat2x2<real_t> mat{0.0_dp, 0.0_dp, //
                             0.0_dp, 0.0_dp};

    CHECK(all(normalize(mat) == mat));
  }
}

// -----------------------------------------------------------------------------

TEST_CASE("Bittern/ComplexMatrixExpressions") {
  const Mat2x2<real_t> mat1{+1.0_dp, -2.0_dp, //
                            +3.0_dp, -4.0_dp};
  const Mat2x2<complex_t> mat2{+3.0_dp - i * 4.0_dp, -5.0_dp + i * 12.0_dp, //
                               -8.0_dp - i * 15.0_dp, +7.0_dp + i * 24.0_dp};

  SUBCASE("real") {
    const Mat2x2<real_t> re_mat2{+3.0_dp, -5.0_dp, //
                                 -8.0_dp, +7.0_dp};

    CHECK(all(real(mat1) == mat1));
    CHECK(all(real(mat2) == re_mat2));
  }

  SUBCASE("imag") {
    const Mat2x2<real_t> im_mat1{0.0_dp, 0.0_dp, //
                                 0.0_dp, 0.0_dp};
    const Mat2x2<real_t> im_mat2{-4.0_dp, 12.0_dp, //
                                 -15.0_dp, 24.0_dp};

    CHECK(all(imag(mat1) == im_mat1));
    CHECK(all(imag(mat2) == im_mat2));
  }

  SUBCASE("conj") {
    const Mat2x2<complex_t> conj_mat2{
        +3.0_dp + i * 4.0_dp, -5.0_dp - i * 12.0_dp, //
        -8.0_dp + i * 15.0_dp, +7.0_dp - i * 24.0_dp};

    CHECK(all(conj(mat1) == mat1));
    CHECK(all(conj(mat2) == conj_mat2));
  }
}

// -----------------------------------------------------------------------------

TEST_CASE("Bittern/AbsMatrixExpressions") {
  SUBCASE("real-abs-sign") {
    const Mat2x2<real_t> mat{+1.0_dp, -2.0_dp, //
                             +3.0_dp, -4.0_dp};

    SUBCASE("abs") {
      const Mat2x2<real_t> abs_mat{1.0_dp, 2.0_dp, //
                                   3.0_dp, 4.0_dp};

      CHECK(all(abs(mat) == abs_mat));
    }

    SUBCASE("sign") {
      const Mat2x2<int> sign_mat{+1, -1, //
                                 +1, -1};

      CHECK(all(sign(mat) == sign_mat));
    }
  }

  SUBCASE("complex-abs") {
    const Mat2x2<complex_t> mat{+3.0_dp - i * 4.0_dp, -5.0_dp + i * 12.0_dp, //
                                -8.0_dp - i * 15.0_dp, +7.0_dp + i * 24.0_dp};

    const Mat2x2<real_t> abs_mat{5.0_dp, 13.0_dp, //
                                 17.0_dp, 25.0_dp};

    CHECK(all(abs(mat) == abs_mat));
  }
}

// -----------------------------------------------------------------------------

TEST_CASE("Bittern/PowerMatrixExpressions") {
  const Mat2x2<real_t> mat{2.0_dp, 3.0_dp, //
                           3.0_dp, 2.0_dp};

  SUBCASE("sqrt/pow(_, 2)/pow(_, 1/2)") {
    const Mat2x2<real_t> sqrt_mat{sqrt2, sqrt3, //
                                  sqrt3, sqrt2};

    CHECK(all(approx_equal(sqrt(mat), sqrt_mat)));
    CHECK(all(approx_equal(pow(mat, 0.5_dp), sqrt_mat)));
    CHECK(all(approx_equal(pow(sqrt_mat, 2), mat)));
  }

  SUBCASE("cbrt/pow(_, 3)/pow(_, 1/3)") {
    const Mat2x2<real_t> cbrt_mat{cbrt(2.0_dp), cbrt(3.0_dp), //
                                  cbrt(3.0_dp), cbrt(2.0_dp)};

    CHECK(all(approx_equal(cbrt(mat), cbrt_mat)));
    CHECK(all(approx_equal(pow(mat, 1.0_dp / 3.0_dp), cbrt_mat)));
    CHECK(all(approx_equal(pow(cbrt_mat, 3), mat)));
  }

  SUBCASE("pow(_, mat)") {
    const Mat2x2<real_t> two_pow_mat{4.0_dp, 8.0_dp, //
                                     8.0_dp, 4.0_dp};

    CHECK(all(approx_equal(pow(2, mat), two_pow_mat)));
  }

  SUBCASE("pow(mat, mat)") {
    const Mat2x2<real_t> mat_pow_mat{4.0_dp, 27.0_dp, //
                                     27.0_dp, 4.0_dp};

    CHECK(all(approx_equal(pow(mat, mat), mat_pow_mat)));
  }
}

// -----------------------------------------------------------------------------

TEST_CASE("Bittern/ExponentialMatrixExpression") {
  const Mat2x2<real_t> mat{0.0_dp, 1.0_dp, //
                           2.0_dp, 3.0_dp};

  SUBCASE("base-e") {
    const Mat2x2<real_t> exp_mat{1.0_dp, e, //
                                 e * e, e * e * e};

    CHECK(all(approx_equal(exp(mat), exp_mat)));
    CHECK(all(approx_equal(log(exp_mat), mat)));
  }

  SUBCASE("base-2") {
    const Mat2x2<real_t> exp2_mat{1.0_dp, 2.0_dp, //
                                  4.0_dp, 8.0_dp};

    CHECK(all(approx_equal(exp2(mat), exp2_mat)));
    CHECK(all(approx_equal(log2(exp2_mat), mat)));
  }

  SUBCASE("base-10") {
    const Mat2x2<real_t> exp10_mat{1.0_dp, 10.0_dp, //
                                   100.0_dp, 1000.0_dp};

    CHECK(all(approx_equal(log10(exp10_mat), mat)));
  }
}

// -----------------------------------------------------------------------------

TEST_CASE("Bittern/TrigonometricMatrixExpressions") {
  const Mat2x2<real_t> mat{+pi / 6.0_dp, -pi / 4.0_dp, //
                           -pi / 4.0_dp, +pi / 6.0_dp};

  SUBCASE("sin/asin") {
    const Mat2x2<real_t> sin_mat{+0.5_dp, -1.0_dp / sqrt2, //
                                 -1.0_dp / sqrt2, +0.5_dp};

    CHECK(all(approx_equal(sin(mat), sin_mat)));
    CHECK(all(approx_equal(asin(sin_mat), mat)));
  }

  SUBCASE("cos/acos") {
    const Mat2x2<real_t> cos_mat{sqrt3 / 2.0_dp, 1.0_dp / sqrt2, //
                                 1.0_dp / sqrt2, sqrt3 / 2.0_dp};

    CHECK(all(approx_equal(cos(mat), cos_mat)));
    CHECK(all(approx_equal(acos(cos_mat), abs(mat))));
  }

  SUBCASE("tan/atan") {
    const Mat2x2<real_t> tan_mat{+1.0_dp / sqrt3, -1.0_dp, //
                                 -1.0_dp, +1.0_dp / sqrt3};

    CHECK(all(approx_equal(tan(mat), tan_mat)));
    CHECK(all(approx_equal(atan(tan_mat), mat)));
  }
}

// -----------------------------------------------------------------------------

/// @todo What is wrong with you, GCC?
#define GCC_COSH_BUG \
  STORM_ARCH_X64&& STORM_COMPILER_GCC && (__GNUC__ == 12) && STORM_NDEBUG
#if GCC_COSH_BUG
#  pragma GCC push_options
#  pragma GCC optimize("O2")
#endif

TEST_CASE("Bittern/HyperbolicMatrixExpressions") {
  const Mat2x2<complex_t> mat{+pi / (6.0_dp * i), -pi / (4.0_dp * i), //
                              -pi / (4.0_dp * i), +pi / (6.0_dp * i)};

  SUBCASE("sinh/asinh") {
    const Mat2x2<complex_t> sinh_mat{-0.5_dp * i, +1.0_dp * i / sqrt2, //
                                     +1.0_dp * i / sqrt2, -0.5_dp * i};

    CHECK(all(approx_equal(sinh(mat), sinh_mat)));
    CHECK(all(approx_equal(asinh(sinh_mat), mat)));
  }

  SUBCASE("cosh/acosh") {
    const Mat2x2<complex_t> cosh_mat{sqrt3 / 2.0_dp, 1.0_dp / sqrt2, //
                                     1.0_dp / sqrt2, sqrt3 / 2.0_dp};

    CHECK(all(approx_equal(cosh(mat), cosh_mat)));
    CHECK(all(approx_equal(acosh(cosh_mat), (1.0_dp * i) * abs(mat))));
  }

  SUBCASE("tanh/atanh") {
    const Mat2x2<complex_t> tanh_mat{-1.0_dp * i / sqrt3, +1.0_dp * i, //
                                     +1.0_dp * i, 1.0_dp * i / -sqrt3};

    CHECK(all(approx_equal(tanh(mat), tanh_mat)));
    CHECK(all(approx_equal(atanh(tanh_mat), mat)));
  }
}

#if GCC_COSH_BUG
#  pragma GCC pop_options
#endif
#undef GCC_COSH_BUG

// -----------------------------------------------------------------------------

} // namespace Storm::UnitTests
