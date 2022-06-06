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

#include <span>

#include <stormBase.hxx>
#include <stormSolvers/Vector.hxx>

namespace Storm {

/// ----------------------------------------------------------------- ///
/// @brief Vector subspace.
/// ----------------------------------------------------------------- ///
template<VectorLike Vector, size_t Extent = std::dynamic_extent>
class Subspace {
private:

  [[no_unique_address]] std::conditional_t<Extent == std::dynamic_extent,
                                           std::vector<Vector>,
                                           std::array<Vector, Extent>>
      Vectors_;

public:

  /// @brief Assign a vector to the each of the subspace vectors.
  /// @{
  void assign(Vector const& like,
              bool copy) requires(Extent != std::dynamic_extent) {
    for (Vector& vector : Vectors_) {
      vector.assign(like, copy);
    }
  }
  void assign(size_t size, Vector const& like,
              bool copy) requires(Extent == std::dynamic_extent) {
    Vectors_.resize(size);
    for (Vector& vector : Vectors_) {
      vector.assign(like, copy);
    }
  }
  /// @}

  /// @brief Access vector at @p index.
  /// @{
  Vector& operator()(size_t index) noexcept {
    StormAssert(index < Vectors_.size());
    return Vectors_[index];
  }
  Vector const& operator()(size_t index) const noexcept {
    StormAssert(index < Vectors_.size());
    return Vectors_[index];
  }
  /// @}

}; // class Subspace

} // namespace Storm
