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
  texture_buffer = GL_TEXTURE_BUFFER,
}; // enum class BufferTarget

/// @brief OpenGL buffer.
template<class Type>
class Buffer : detail_::noncopyable_ {
private:

  GLuint buffer_id_;

public:

  /// @brief Construct a buffer.
  Buffer() {
    glGenBuffers(1, &buffer_id_);
  }

  /// @brief Construct a buffer with a range.
  /// @param usage Intended buffer usage.
  /// @param target Buffer binding target.
  // clang-format off
  template<std::ranges::input_range Range>
    requires std::same_as<std::ranges::range_value_t<Range>, Type>
  Buffer(Range&& values, BufferUsage usage = BufferUsage::static_draw) : Buffer{} {
    // clang-format on
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

  /// @brief Bind the buffer to @p target.
  void bind(BufferTarget target) const {
    glBindBuffer(static_cast<GLenum>(target), buffer_id_);
  }

  /// @brief Assign the buffer @p values.
  /// @param usage Intended buffer usage.
  // clang-format off
  template<std::ranges::input_range Range>
    requires std::same_as<std::ranges::range_value_t<Range>, Type>
  void assign(Range&& values, BufferUsage usage = BufferUsage::static_draw) {
    // clang-format on
    glBindBuffer(GL_ARRAY_BUFFER, buffer_id_);
    if constexpr (std::ranges::sized_range<Range>) {
      const auto num_bytes =
          static_cast<GLsizeiptr>(values.size() * sizeof(Type));
      glBufferData(GL_ARRAY_BUFFER, num_bytes, /*data*/ nullptr,
                   static_cast<GLenum>(usage));
      const auto pointer =
          static_cast<Type*>(glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY));
      if (pointer == nullptr) { STORM_THROW_GL_("Failed to map the buffer!"); }
      std::ranges::copy(values, pointer);
      if (glUnmapBuffer(GL_ARRAY_BUFFER) != GL_TRUE) {
        STORM_THROW_GL_("Failed to unmap the buffer!");
      }
    } else {
      std::vector<Type> temp{};
      std::ranges::copy(values, std::back_inserter(temp));
      const auto num_bytes =
          static_cast<GLsizeiptr>(temp.size() * sizeof(Type));
      glBufferData(GL_ARRAY_BUFFER, num_bytes, temp.data(),
                   static_cast<GLenum>(usage));
    }
  }

  /// @brief Set the buffer @p value at @p index.
  void set(size_t index, const Type& value) {
    glBindBuffer(GL_ARRAY_BUFFER, buffer_id_);
    const auto offset = static_cast<GLintptr>(index * sizeof(Type));
    glBufferSubData(GL_ARRAY_BUFFER, offset, sizeof(Type), &value);
  }

}; // class Buffer

template<std::ranges::input_range Range, class... T>
Buffer(Range&&, T...) -> Buffer<std::ranges::range_value_t<Range>>;

/// @brief Vertex attribute type traits.
template<class>
class VertexAttribTypeTraits {
public:

  /// @brief Number of the components in the type.
  [[nodiscard]] static consteval GLint length() noexcept = delete;

  /// @brief Type enumeration.
  [[nodiscard]] static consteval GLenum type_enum() noexcept = delete;

}; // class VertexAttribTypeTraits

/// @brief Vertex attribute type.
// clang-format off
template<class Type>
concept vertex_attrib_type =
    requires {
      { VertexAttribTypeTraits<Type>::length() } -> std::same_as<GLint>;
      { VertexAttribTypeTraits<Type>::type_enum() } -> std::same_as<GLenum>;
    };
// clang-format on

/// @brief Number of the components in type.
template<class Type>
inline constexpr size_t
    vertex_attrib_length_v = VertexAttribTypeTraits<Type>::length();

/// @brief OpenGL type enumeration.
template<class Type>
inline constexpr GLenum
    vertex_attrib_type_enum_v = VertexAttribTypeTraits<Type>::type_enum();

#define MAKE_SCALAR_BASE_TYPE_(GLtype, GL_TYPE)                  \
  template<>                                                     \
  struct VertexAttribTypeTraits<GLtype> {                        \
    [[nodiscard]] static consteval GLint length() noexcept {     \
      return 1;                                                  \
    }                                                            \
    [[nodiscard]] static consteval GLenum type_enum() noexcept { \
      return GL_TYPE;                                            \
    }                                                            \
  };

// Scalar integral types.
MAKE_SCALAR_BASE_TYPE_(GLbyte, GL_BYTE)
MAKE_SCALAR_BASE_TYPE_(GLubyte, GL_UNSIGNED_BYTE)
MAKE_SCALAR_BASE_TYPE_(GLshort, GL_SHORT)
MAKE_SCALAR_BASE_TYPE_(GLushort, GL_UNSIGNED_SHORT)
MAKE_SCALAR_BASE_TYPE_(GLint, GL_INT)
MAKE_SCALAR_BASE_TYPE_(GLuint, GL_UNSIGNED_INT)

// Scalar floating-point types.
MAKE_SCALAR_BASE_TYPE_(GLfloat, GL_FLOAT)
MAKE_SCALAR_BASE_TYPE_(GLdouble, GL_DOUBLE)

#undef MAKE_SCALAR_BASE_TYPE_

// GLM vector type.
template<glm::length_t Length, vertex_attrib_type Type>
struct VertexAttribTypeTraits<glm::vec<Length, Type, glm::packed>> {
  [[nodiscard]] static consteval GLint length() noexcept {
    return static_cast<GLint>(Length * vertex_attrib_length_v<Type>);
  }
  [[nodiscard]] static consteval GLenum type_enum() noexcept {
    return vertex_attrib_type_enum_v<Type>;
  }
}; // VertexAttribTypeTraits<glm::vec>

// GLM matrix type.
template<glm::length_t Rows, glm::length_t Cols, vertex_attrib_type Type>
struct VertexAttribTypeTraits<glm::mat<Rows, Cols, Type, glm::packed>> {
  [[nodiscard]] static consteval GLint length() noexcept {
    return static_cast<GLint>(Rows * Cols * vertex_attrib_length_v<Type>);
  }
  [[nodiscard]] static consteval GLenum type_enum() noexcept {
    return vertex_attrib_type_enum_v<Type>;
  }
}; // VertexAttribTypeTraits<glm::mat>

/// @brief OpenGL vertex array.
class VertexArray final : detail_::noncopyable_ {
public:

  GLuint vertex_array_id_;

public:

  /// @brief Construct a vertex array.
  VertexArray() {
    glGenVertexArrays(1, &vertex_array_id_);
  }

  /// @brief Move-construct a vertex array.
  VertexArray(VertexArray&&) = default;
  /// @brief Move-assign the vertex array.
  VertexArray& operator=(VertexArray&&) = default;

  /// @brief Destruct the vertex array.
  ~VertexArray() {
    glDeleteVertexArrays(1, &vertex_array_id_);
  }

  /// @brief Build the vertex array buffer.
  /// @{
  template<vertex_attrib_type... Types>
  void build(const Buffer<Types>&... vertex_buffers) {
    glBindVertexArray(vertex_array_id_);
    (attach_vertex_attribs_(vertex_buffers), ...);
  }
  template<vertex_attrib_type... Types>
  void build_indexed(const Buffer<GLuint>& index_buffer,
                     const Buffer<Types>&... vertex_buffers) {
    glBindVertexArray(vertex_array_id_);
    index_buffer.bind(BufferTarget::element_array_buffer);
    (attach_vertex_attribs_(vertex_buffers), ...);
  }
  /// @}

private:

  template<vertex_attrib_type... Types>
  static void attach_vertex_attribs_(const Buffer<Types>&... vertex_buffers) {
    GLuint index{};
    const auto attach_single_attrib =
        [&]<class Type>(const Buffer<Type>& vertex_buffer) {
          vertex_buffer.bind(BufferTarget::array_buffer);
          glVertexAttribPointer(index, //
                                vertex_attrib_length_v<Type>,
                                vertex_attrib_type_enum_v<Type>,
                                /*normalized*/ GL_FALSE,
                                /*stride*/ 0,
                                /*offset*/ nullptr);
          glEnableVertexAttribArray(index);
          index += 1;
        };
    (attach_single_attrib(vertex_buffers), ...);
  }

}; // class VertexArray

} // namespace Storm::Vulture::gl
