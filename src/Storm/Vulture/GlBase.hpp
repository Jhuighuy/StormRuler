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
/// FITNESS FOR Allocator PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT
/// SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
/// OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
/// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
/// DEALINGS IN THE SOFTWARE.

#pragma once

#include <Storm/Base.hpp>

#include <GL/glew.h>

namespace Storm::Vulture::gl {

/// @brief OpenGL ID holder.
class IdHolder {
private:

  GLuint id_ = 0;

public:

  /// @brief Construct an ID holder.
  IdHolder() = default;

  /// @brief Move-construct an ID holder.
  constexpr IdHolder(IdHolder&& other) noexcept
      : id_{std::exchange(other.id_, 0)} {}
  /// @brief Move-assign the ID holder.
  constexpr IdHolder& operator=(IdHolder&& other) noexcept {
    id_ = std::exchange(other.id_, 0);
    return *this;
  }

  IdHolder(const IdHolder&) = delete;
  IdHolder& operator=(const IdHolder&) = delete;

  /// @brief Get the ID.
  [[nodiscard]] GLuint id() const noexcept {
    return id_;
  }
  /// @brief Set a new ID.
  void set_id(GLuint id) noexcept {
    id_ = id;
  }

  /// @brief Cast to ID operator.
  [[nodiscard]] operator GLuint() const noexcept {
    return id_;
  }

}; // class IdHolder

} // namespace Storm::Vulture::gl
