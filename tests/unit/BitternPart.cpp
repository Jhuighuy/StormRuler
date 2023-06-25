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
#include <Storm/Utils/Meta.hpp>

#include <doctest/doctest.h>
#include <ranges>

namespace Storm::UnitTests {

// -----------------------------------------------------------------------------

TEST_CASE("Bittern/MatrixSlice") {
  const DenseMatrix mat{{1, 2, 3}, //
                        {4, 5, 6}};

  SUBCASE("slice-basic") {
    CHECK(at(mat, 0, 0) == 1);
    CHECK(all(at(mat, _, 0) == DenseMatrix{1, 4}));
    CHECK(all(at(mat, 0, _) == DenseMatrix{1, 2, 3}));
    CHECK(all(at(mat, _, _) == mat));
  }

  SUBCASE("slice-padded") {
    CHECK(all(at(mat) == mat));
    CHECK(all(at(mat, 0) == DenseMatrix{1, 2, 3}));
  }

  SUBCASE("slice-iota") {
    const DenseMatrix slice{{1, 2}, //
                            {4, 5}};

    CHECK(all(at(mat, _, std::views::iota(0, 2)) == slice));
  }
}

// -----------------------------------------------------------------------------

TEST_CASE("Bittern/MatrixSelect") {
  const DenseMatrix mat{{1, 2, 3}, //
                        {4, 5, 6}};

  SUBCASE("select") {
    CHECK(all(at(mat, _, std::array{0, 1}) == DenseMatrix{{1, 2}, {4, 5}}));
    CHECK(all(at(mat, _, std::array{1, 0}) == DenseMatrix{{2, 1}, {5, 4}}));
  }
}

// -----------------------------------------------------------------------------

TEST_CASE("Bittern/MatrixFilter") {
  const DenseMatrix mat{{1, 2, 3}, //
                        {4, 5, 6}};

  SUBCASE("filter") {
    const DenseMatrix mask1{{false, true, false}, //
                            {true, false, true}};
    const DenseMatrix mask2{false, true, false};

    CHECK(all(at(mat, mask1) == DenseMatrix{2, 4, 6}));
    CHECK(all(at(mat, _, mask2) == DenseMatrix{{2}, {5}}));
  }
}

// -----------------------------------------------------------------------------

TEST_CASE("Bittern/MatrixNewAxis") {
  SUBCASE("new-axis") {
    const DenseMatrix mat{1, 2, 3};

    CHECK(all(at(mat, new_axis, _) == FixedMatrix<int, 1, 3>{{1, 2, 3}}));
    CHECK(all(at(mat, _, new_axis) == DenseMatrix{{1}, {2}, {3}}));
  }
}

// -----------------------------------------------------------------------------

} // namespace Storm::UnitTests
