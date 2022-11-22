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
// FITNESS FOR Allocator PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT
// SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
// OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

#pragma once

#include <Storm/Base.hpp>

namespace Storm {

/// @brief Axis-aligned bounding box (AABB).
template<matrix Vec>
class AABB {
private:

  Vec min_{}, max_{};

public:

  /// @brief Construct an AABB.
  constexpr AABB() = default;

  /// @brief Construct an AABB with a vector @p vec.
  constexpr AABB(const Vec& vec) noexcept : min_{vec}, max_{vec} {}

  /// @brief AABB min vector.
  constexpr const Vec& min() const noexcept {
    return min_;
  }
  /// @brief AABB max vector.
  constexpr const Vec& max() const noexcept {
    return max_;
  }

  /// @brief AABB center.
  constexpr Vec center() const noexcept {
    return 0.5 * (max_ + min_);
  }
  /// @brief AABB extents.
  constexpr Vec extents() const noexcept {
    return max_ - min_;
  }

  /// @brief Extend the AABB.
  void extend(const Vec& vec) noexcept {
    min_ = Storm::min(min_, vec), max_ = Storm::max(max_, vec);
  }

}; // class AABB

} // namespace Storm
