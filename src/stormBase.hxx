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
#ifndef _STORM_BASE_HXX_
#define _STORM_BASE_HXX_

#include <cstddef>
#include <array>
#include <vector>
#include <type_traits>
#include <iostream>

#include <StormRuler_API.h>

#define _STORM_NOT_IMPLEMENTED_() do { \
  std::cerr << __FUNCTION__ << " not implemented" << std::endl; std::exit(1); } while(false)

#define StormEnabledAssert(x) assert(x)
#define StormDisabledAssert(x) static_cast<void>(x)
#define StormAssert stormAssert
#define _STORM_MESH_DEBUG_

namespace Storm {

using size_t = std::size_t;
using ptrdiff_t = std::ptrdiff_t;
using real_t = double;

static size_t const DynamicExtent = SIZE_MAX;

template<size_t Extent>
using ExtentSize = 
  std::conditional_t<Extent == DynamicExtent, 
                     size_t, 
                     std::integral_constant<size_t, Extent>>;

template<class Value, size_t Extent>
using ExtentArray =
  std::conditional_t<Extent == DynamicExtent,
                     std::vector<Value>,
                     std::array<Value, Extent>>;

} // namespace Storm

#endif // ifndef _STORM_BASE_HXX_
