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

TEST_CASE("Bittern/TransposeVector") {
  const DenseMatrix vec{+3.0_dp - i * 4.0_dp, //
                        -5.0_dp + i * 6.0_dp};

  SUBCASE("transpose") {
    const FixedMatrix<complex_t, 1, 2> transposed_vec{
        {+3.0_dp - i * 4.0_dp, -5.0_dp + i * 6.0_dp}};

    // Transpose:
    CHECK(all(transpose(vec) == transposed_vec));

    // Identity transpose:
    CHECK(all(transpose<0>(vec) == vec));
  }

  SUBCASE("conj-transpose") {
    const FixedMatrix<complex_t, 1, 2> conj_transposed_vec{
        {+3.0_dp + i * 4.0_dp, -5.0_dp - i * 6.0_dp}};

    // Conjugate transpose:
    CHECK(all(conj_transpose(vec) == conj_transposed_vec));

    // Conjugate identity transpose:
    CHECK(all(conj_transpose<0>(vec) == conj(vec)));
  }
}

// -----------------------------------------------------------------------------

TEST_CASE("Bittern/TransposeMatrix") {
  const DenseMatrix mat{{+1.0_dp - i * 2.0_dp, -3.0_dp + i * 4.0_dp},
                        {+3.0_dp - i * 4.0_dp, -5.0_dp + i * 6.0_dp},
                        {+5.0_dp - i * 6.0_dp, -7.0_dp + i * 8.0_dp}};

  SUBCASE("transpose") {
    const DenseMatrix transposed_mat{
        {+1.0_dp - i * 2.0_dp, +3.0_dp - i * 4.0_dp, +5.0_dp - i * 6.0_dp},
        {-3.0_dp + i * 4.0_dp, -5.0_dp + i * 6.0_dp, -7.0_dp + i * 8.0_dp}};

    // Transpose:
    CHECK(all(transpose(mat) == transposed_mat));
    CHECK(all(transpose<1, 0>(mat) == transposed_mat));

    // Identity transpose:
    CHECK(all(transpose<0, 1>(mat) == mat));
  }

  SUBCASE("conj-transpose") {
    const DenseMatrix conj_transposed_mat{
        {+1.0_dp + i * 2.0_dp, +3.0_dp + i * 4.0_dp, +5.0_dp + i * 6.0_dp},
        {-3.0_dp - i * 4.0_dp, -5.0_dp - i * 6.0_dp, -7.0_dp - i * 8.0_dp}};

    // Conjugate transpose:
    CHECK(all(conj_transpose(mat) == conj_transposed_mat));
    CHECK(all(conj_transpose<1, 0>(mat) == conj_transposed_mat));

    // Conjugate identity transpose:
    CHECK(all(conj_transpose<0, 1>(mat) == conj(mat)));
  }
}

// -----------------------------------------------------------------------------

TEST_CASE("Bittern/RollMatrix") {
  SUBCASE("roll-1D") {
    const DenseMatrix mat{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

    CHECK(all(roll(mat, +2) == DenseMatrix{8, 9, 0, 1, 2, 3, 4, 5, 6, 7}));
    CHECK(all(roll(mat, -2) == DenseMatrix{2, 3, 4, 5, 6, 7, 8, 9, 0, 1}));
  }

  SUBCASE("roll-2D") {
    const DenseMatrix mat{{0, 1, 2, 3, 4}, //
                          {5, 6, 7, 8, 9}};

    CHECK(all(roll(mat, +1) == DenseMatrix{{5, 6, 7, 8, 9}, //
                                           {0, 1, 2, 3, 4}}));
    CHECK(all(roll(mat, -1) == DenseMatrix{{5, 6, 7, 8, 9}, //
                                           {0, 1, 2, 3, 4}}));

    CHECK(all(roll<1>(mat, +1) == DenseMatrix{{4, 0, 1, 2, 3}, //
                                              {9, 5, 6, 7, 8}}));
    CHECK(all(roll<1>(mat, -1) == DenseMatrix{{1, 2, 3, 4, 0}, //
                                              {6, 7, 8, 9, 5}}));

    CHECK(all(roll<0, 1>(mat, 1, 1) == DenseMatrix{{9, 5, 6, 7, 8}, //
                                                   {4, 0, 1, 2, 3}}));
    CHECK(all(roll<1, 0>(mat, 2, 1) == DenseMatrix{{8, 9, 5, 6, 7}, //
                                                   {3, 4, 0, 1, 2}}));
  }
}

// -----------------------------------------------------------------------------

} // namespace Storm::UnitTests
