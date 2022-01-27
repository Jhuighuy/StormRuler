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
#ifndef _STORM_SUBSPACE_HXX_
#define _STORM_SUBSPACE_HXX_

#include <array>
#include <vector>
#include <type_traits>

#include <stormBlas/stormShape.hxx>

/// ----------------------------------------------------------------- ///
/// @brief Vector subspace.
/// ----------------------------------------------------------------- ///
template<class Vector, stormSize_t Extent = stormDynamicExtent>
class stormSubspace {
private:
  [[no_unique_address]] stormShape<Extent> Shape_ = {};
  [[no_unique_address]] std::conditional_t<
    Extent == stormDynamicExtent,
    std::vector<Vector>, 
    std::array<Vector, Extent>> Vectors_; 

public:

  /// @brief Assign a vector to the each of the subspace vectors.
  /// @{
  void Assign(Vector const& like, bool copy) 
      requires(Extent != stormDynamicExtent) {
    for (Vector& vector : Vectors_) {
      vector.Assign(like, copy);
    }
  }
  void Assign(stormSize_t size, Vector const& like, bool copy) 
      requires(Extent == stormDynamicExtent) {
    Vectors_.resize(size);
    for (Vector& vector : Vectors_) {
      vector.Assign(like, copy);
    }
  }
  /// @}

  /// @brief Access vector at index.
  /// @{
  Vector& operator()(stormSize_t index) noexcept {
    return Vectors_[Shape_(index)];
  }
  Vector const& operator()(stormSize_t index) const noexcept {
    return Vectors_[Shape_(index)];
  }
  /// @}

}; // class stormSubspace

#endif // ifndef _STORM_SUBSPACE_HXX_
