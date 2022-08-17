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

#pragma once

#include <Storm/Base.hpp>

#include <array>
#include <glm/glm.hpp>

namespace Storm {

template<class, size_t>
using FastVector = glm::dvec3;

#if 0
template<class Element, size_t NumRows>
struct FastVectorData;

#if STORM_COMPILER_GCC_
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#elif STORM_COMPILER_CLANG_
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wgnu-anonymous-struct"
#pragma clang diagnostic ignored "-Wnested-anon-types"
#elif STORM_COMPILER_MSVC_
#pragma warning(push)
#pragma warning(disable : 4201)
#endif

/// @brief Data for a fast 1D vector.
template<class Element>
struct FastVectorData<Element, 1> {
  union {
    std::array<Element, 1> data;
    Element x;
    Element r;
  };
}; // struct FastVectorData<Element, 1>

/// @brief Data for a fast 2D vector.
template<class Element>
struct FastVectorData<Element, 2> {
  union {
    std::array<Element, 2> data;
    struct {
      Element x, y;
    };
    struct {
      Element r, g;
    };
  };
}; // struct FastVectorData<Element, 2>

/// @brief Data for a fast 3D vector.
template<class Element>
struct FastVectorData<Element, 3> {
  union {
    std::array<Element, 3> data;
    struct {
      Element x, y, z;
    };
    struct {
      Element r, g, b;
    };
  };
}; // struct FastVectorData<Element, 3>

/// @brief Data for a 4D fast vector.
template<class Element>
struct FastVectorData<Element, 4> {
  union {
    std::array<Element, 4> data;
    struct {
      Element x, y, z, w;
    };
    struct {
      Element r, g, b, a;
    };
  };
}; // struct FastVectorData<Element, 4>

/// @brief Data for a 5D+ fast vector.
template<class Element, size_t NumRows>
struct FastVectorData {
  union {
    std::array<Element, NumRows> data;
  };
}; // struct FastVectorData<Element, 5+>

#if STORM_COMPILER_GCC_
#pragma GCC diagnostic pop
#elif STORM_COMPILER_CLANG_
#pragma clang diagnostic pop
#elif STORM_COMPILER_MSVC_
#pragma warning(pop)
#endif

/// @brief Fixed-sized vector.
template<class Element, size_t NumRows>
class FastVector : public FastVectorData<Element, NumRows> {
public:

  constexpr FastVector() = default;

}; // class FastVector
#endif

} // namespace Storm
