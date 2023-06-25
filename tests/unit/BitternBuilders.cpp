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

#include <Storm/Bittern/MatrixDense.hpp>

#include <doctest/doctest.h>

namespace Storm::UnitTests {

// -----------------------------------------------------------------------------

TEST_CASE("Bittern/CustomBuilders") {
  SUBCASE("weird") {
    const DenseMatrix mat{{3_sz, 4_sz}, //
                          {5_sz, 6_sz}};

    auto gen_func = [](size_t i, size_t j) { return 2_sz * ++i + ++j; };
    CHECK(all(build(std::tuple{2_sz, 2_sz}, gen_func) == mat));
  }
}

// -----------------------------------------------------------------------------

TEST_CASE("Bittern/ValuesBuilders") {
  SUBCASE("values") {
    const DenseMatrix threes_mat{{3.0_dp, 3.0_dp}, //
                                 {3.0_dp, 3.0_dp}};

    CHECK(all(values(std::array{2_sz, 2_sz}, 3.0_dp) == threes_mat));
  }

  SUBCASE("zeroes-1D") {
    const DenseMatrix zeroes_vec{0.0_dp, 0.0_dp};

    CHECK(all(zeroes(2_sz) == zeroes_vec));
  }

  SUBCASE("zeroes-2D") {
    const DenseMatrix zeroes_mat{{0.0_dp, 0.0_dp}, //
                                 {0.0_dp, 0.0_dp}};

    CHECK(all(zeroes(2_sz, 2_sz) == zeroes_mat));
  }

  SUBCASE("ones") {
    const DenseMatrix ones_mat{{1.0_dp, 1.0_dp}, //
                               {1.0_dp, 1.0_dp}};

    CHECK(all(ones(2_sz, 2_sz) == ones_mat));
  }
}

// -----------------------------------------------------------------------------

TEST_CASE("Bittern/EyeBuilders") {
  SUBCASE("eye-1D") {
    const DenseMatrix eye_vec{1.0_dp, 0.0_dp, 0.0_dp};

    CHECK(all(eye(3_sz) == eye_vec));
  }

  SUBCASE("eye-2D") {
    const DenseMatrix eye_mat{{1.0_dp, 0.0_dp, 0.0_dp}, //
                              {0.0_dp, 1.0_dp, 0.0_dp}, //
                              {0.0_dp, 0.0_dp, 1.0_dp}};

    CHECK(all(eye(3_sz, 3_sz) == eye_mat));
  }

  SUBCASE("eye-2D-custom-values") {
    const DenseMatrix custom_eye_mat{{0.5_dp, 1.0_dp, 1.0_dp}, //
                                     {1.0_dp, 0.5_dp, 1.0_dp}, //
                                     {1.0_dp, 1.0_dp, 0.5_dp}};

    CHECK(all(eye(std::array{3_sz, 3_sz}, 0.5_dp, 1.0_dp) == custom_eye_mat));
  }

  SUBCASE("eye-custom-type") {
    using element_t =
        matrix_element_t<decltype(eye(std::array{1_sz}, 1_sz, 0_sz))>;
    static_assert(std::is_same_v<size_t, element_t>);
  }
}

// -----------------------------------------------------------------------------

} // namespace Storm::UnitTests
