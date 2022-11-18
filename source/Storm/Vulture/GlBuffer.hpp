// Copyright (C) 2020 - 2023 Oleg Butakov
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

#include <algorithm>
#include <concepts>
#include <ranges>
#include <vector>

#include <GL/glew.h>
#include <glm/glm.hpp>

namespace Storm::Vulture::gl {

/// @brief OpenGL buffer usage.
enum class BufferUsage : GLenum {
  stream_draw = GL_STREAM_DRAW,
  stream_read = GL_STREAM_READ,
  stream_copy = GL_STREAM_COPY,
  static_draw = GL_STATIC_DRAW,
  static_read = GL_STATIC_READ,
  static_copy = GL_STATIC_COPY,
  dynamic_draw = GL_DYNAMIC_DRAW,
  dynamic_read = GL_DYNAMIC_READ,
  dynamic_copy = GL_DYNAMIC_COPY,
}; // enum class BufferUsage

/// @brief OpenGL buffer binding target.
enum class BufferTarget : GLenum {
  array_buffer = GL_ARRAY_BUFFER,
  element_array_buffer = GL_ELEMENT_ARRAY_BUFFER,
  pixel_pack_buffer = GL_PIXEL_PACK_BUFFER,
  pixel_unpack_buffer = GL_PIXEL_UNPACK_BUFFER,
  texture_buffer = GL_TEXTURE_BUFFER,
}; // enum class BufferTarget

/// @brief OpenGL buffer.
template<class Type>
class Buffer : detail_::noncopyable_ {
private:

  GLuint buffer_id_;
  GLsizei buffer_size_ = 0;

public:

  /// @brief Construct a buffer.
  Buffer() {
    glGenBuffers(1, &buffer_id_);
  }

  /// @brief Construct a buffer with a size.
  /// @param usage Intended buffer usage.
  explicit Buffer(GLsizei size, BufferUsage usage = BufferUsage::static_draw)
      : Buffer{} {
    assign(size, usage);
  }

  /// @brief Construct a buffer with a @p size copies of @p value.
  /// @param usage Intended buffer usage.
  explicit Buffer(GLsizei size, const Type& value,
                  BufferUsage usage = BufferUsage::static_draw)
      : Buffer{} {
    assign(size, value, usage);
  }

  /// @brief Construct a buffer with a range.
  /// @param usage Intended buffer usage.
  template<std::ranges::input_range Range>
    requires std::same_as<std::ranges::range_value_t<Range>, Type>
  Buffer(Range&& values, BufferUsage usage = BufferUsage::static_draw)
      : Buffer{} {
    assign(std::forward<Range>(values), usage);
  }

  /// @brief Move-construct a buffer.
  Buffer(Buffer&&) = default;
  /// @brief Move-assign the buffer.
  Buffer& operator=(Buffer&&) = default;

  /// @brief Destruct the buffer.
  ~Buffer() {
    glDeleteBuffers(1, &buffer_id_);
  }

  /// @brief Cast to buffer ID.
  [[nodiscard]] constexpr operator GLuint() const noexcept {
    return buffer_id_;
  }

  /// @brief Buffer size.
  [[nodiscard]] constexpr GLsizei size() const noexcept {
    return buffer_size_;
  }

  /// @brief Bind the buffer to @p target.
  void bind(BufferTarget target) const {
    glBindBuffer(static_cast<GLenum>(target), buffer_id_);
  }

  /// @brief Assign the buffer @p size.
  /// @param usage Intended buffer usage.
  void assign(GLsizei size, BufferUsage usage = BufferUsage::static_draw) {
    glBindBuffer(GL_ARRAY_BUFFER, buffer_id_);
    buffer_size_ = size;
    const auto num_bytes = static_cast<GLsizeiptr>(size * sizeof(Type));
    glBufferData(GL_ARRAY_BUFFER, num_bytes, /*data*/ nullptr,
                 static_cast<GLenum>(usage));
  }

  /// @brief Assign the buffer @p size copies of @p value.
  void assign(GLsizei size, const Type& value,
              BufferUsage usage = BufferUsage::static_draw) {
    /// @todo C++23:
    // assign(std::views::repeat(value, size));
    assign(std::views::iota(size) |
           std::views::transform([&](auto) { return value; }));
  }

  /// @brief Assign the buffer @p values.
  /// @param usage Intended buffer usage.
  template<std::ranges::input_range Range>
    requires std::same_as<std::ranges::range_value_t<Range>, Type>
  void assign(Range&& values, BufferUsage usage = BufferUsage::static_draw) {
    glBindBuffer(GL_ARRAY_BUFFER, buffer_id_);
    if constexpr (std::ranges::sized_range<Range>) {
      assign(static_cast<GLsizei>(values.size()));
      const auto pointer = glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
      if (pointer == nullptr) { STORM_THROW_GL_("Failed to map the buffer!"); }
      std::ranges::copy(values, static_cast<Type*>(pointer));
      if (glUnmapBuffer(GL_ARRAY_BUFFER) != GL_TRUE) {
        STORM_THROW_GL_("Failed to unmap the buffer!");
      }
    } else {
      std::vector<Type> temp{};
      std::ranges::copy(values, std::back_inserter(temp));
      buffer_size_ = static_cast<GLsizei>(temp.size());
      const auto num_bytes =
          static_cast<GLsizeiptr>(buffer_size_ * sizeof(Type));
      glBufferData(GL_ARRAY_BUFFER, num_bytes, temp.data(),
                   static_cast<GLenum>(usage));
    }
  }

  /// @todo Refactor as `operator[]`.
  /// @brief Get the buffer value at @p index.
  Type get(size_t index) const {
    STORM_ASSERT_(index < static_cast<size_t>(buffer_size_),
                  "Index is out of range!");
    glBindBuffer(GL_ARRAY_BUFFER, buffer_id_);
    const auto pointer = glMapBufferRange(
        GL_ARRAY_BUFFER, static_cast<GLintptr>(index * sizeof(Type)),
        static_cast<GLsizeiptr>(sizeof(Type)), GL_MAP_READ_BIT);
    if (pointer == nullptr) { STORM_THROW_GL_("Failed to map the buffer!"); }
    Type value = *static_cast<const Type*>(pointer);
    if (glUnmapBuffer(GL_ARRAY_BUFFER) != GL_TRUE) {
      STORM_THROW_GL_("Failed to unmap the buffer!");
    }
    return value;
  }

  /// @brief Set the buffer @p value at @p index.
  void set(size_t index, const Type& value) {
    STORM_ASSERT_(index < static_cast<size_t>(buffer_size_),
                  "Index is out of range!");
    glBindBuffer(GL_ARRAY_BUFFER, buffer_id_);
    const auto offset = static_cast<GLintptr>(index * sizeof(Type));
    glBufferSubData(GL_ARRAY_BUFFER, offset, sizeof(Type), &value);
  }

}; // class Buffer

template<std::ranges::input_range Range, class... T>
Buffer(Range&&, T...) -> Buffer<std::ranges::range_value_t<Range>>;

} // namespace Storm::Vulture::gl
