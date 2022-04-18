/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// Copyright (C) 2022 Oleg Butakov
///
/// Permission is hereby granted, free of charge, to any person
/// obtaining a copy of this software and associated documentation
/// files (the "Software"), to deal in the Software without
/// restriction, including without limitation the rights  to use,
/// copy, modify, merge, publish, distribute, sublicense, and/or
/// sell copies of the Software, and to permit persons to whom the
/// Software is furnished to do so, subject to the following
/// conditions:
///
/// The above copyright notice and this permission notice shall be
/// included in all copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
/// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
/// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
/// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
/// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
/// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
/// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
/// OTHER DEALINGS IN THE SOFTWARE.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///

#pragma once

#include <stormBase.hxx>

namespace Storm::Turbo {

template<class Function>
void ParFor(size_t beginX, size_t endX,
            size_t beginY, size_t endY, Function const& func) {

  for (size_t ix = beginX; ix != endX; ++ix) {
    for (size_t iy = beginY; iy != endY; ++iy) {
      func(ix, iy);
    }
  }

} // ParFor

template<class Function>
void ParFor(size_t beginX, size_t endX,
            size_t beginY, size_t endY,
            size_t beginZ, size_t endZ, Function const& func) {

  for (size_t ix = beginX; ix != endX; ++ix) {
    for (size_t iy = beginY; iy != endY; ++iy) {
      for (size_t iz = beginZ; iz != endZ; ++iz) {
        func(ix, iy, iz);
      }
    }
  }

} // ParFor

} // namespace Storm::Turbo
