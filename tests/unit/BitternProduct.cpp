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

TEST_CASE("Bittern/CrossProduct") {
  const DenseMatrix vec1{1.0_dp, 0.0_dp, 0.0_dp};

  SUBCASE("cross") {
    const DenseMatrix vec2{0.0_dp, 1.0_dp, 0.0_dp};
    const DenseMatrix vec3{0.0_dp, 0.0_dp, 1.0_dp};

    CHECK(all(cross_product(vec1, vec2) == vec3));
    CHECK(all(cross_product(vec2, vec1) == -vec3));
  }

  SUBCASE("cross-linear-dependent") {
    CHECK(all(cross_product(2.0_dp * vec1, 3.0_dp * vec1) == zeroes(3_sz)));
  }
}

// -----------------------------------------------------------------------------

} // namespace Storm::UnitTests
