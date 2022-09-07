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

#include <Storm/Vulture/GlBase.hpp>

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef> // std::byte
#include <ranges>
#include <vector>

#include <GL/glew.h>
#include <glm/glm.hpp>

namespace Storm::Vulture::gl {

/// @brief OpenGL device buffer.
class DeviceBuffer : public IdHolder {
public:

  /// @brief Construct a buffer.
  DeviceBuffer() {
    GLuint buffer_id;
    glGenBuffers(1, &buffer_id);
    set_id(buffer_id);
  }

  /// @brief Move-construct a buffer.
  DeviceBuffer(DeviceBuffer&&) = default;

  /// @brief Move-assign the buffer.
  DeviceBuffer& operator=(DeviceBuffer&&) = default;

  /// @brief Destruct the buffer.
  ~DeviceBuffer() {
    const GLuint buffer_id = id();
    glDeleteBuffers(1, &buffer_id);
  }

}; // class DeviceBuffer

/// @brief OpenGL vertex array.
class VertexArray final : public IdHolder {
public:

  /// @brief Construct a vertex array.
  VertexArray() {
    GLuint vertex_array_id;
    glGenVertexArrays(1, &vertex_array_id);
    set_id(vertex_array_id);
  }

  /// @brief Move-construct a vertex array.
  VertexArray(VertexArray&&) = default;

  /// @brief Move-assign the vertex array.
  VertexArray& operator=(VertexArray&&) = default;

  /// @brief Destruct the vertex array.
  ~VertexArray() {
    glBindVertexArray(0);
    const GLuint vertex_array_id = id();
    glDeleteVertexArrays(1, &vertex_array_id);
  }

}; // class VertexArray

/// @brief RAII binder of an OpenGL vertex array.
class BindVertexArray final {
public:

  /// @brief Bind the @p vertex_array.
  explicit BindVertexArray(const VertexArray& vertex_array) {
    STORM_ASSERT_(vertex_array != 0, "Invalid vertex array!");
    glBindVertexArray(vertex_array);
  }

  /// @brief Unbind the vertex_array.
  ~BindVertexArray() {
    glBindVertexArray(0);
  }

}; // class BindProgram

/// @brief Vertex atttribute type traits.
template<class>
class VertexAttribTypeTraits {
public:

  /// @brief Type that would be stored in an OpenGL buffer.
  struct storage_type;

  /// @brief OpenGL type enumeration.
  [[nodiscard]] static constexpr GLenum type() noexcept = delete;

  /// @brief Number of the components in type.
  [[nodiscard]] static constexpr GLint length() noexcept = delete;

}; // class VertexAttribTypeTraits

/// @brief Vertex attribute type.
// clang-format off
template<class Type>
concept vertex_attrib_valid_type =
    requires {
      typename VertexAttribTypeTraits<Type>::storage_type;
      { VertexAttribTypeTraits<Type>::type() } -> std::same_as<GLenum>;
      { VertexAttribTypeTraits<Type>::length() } -> std::same_as<GLint>;
    } &&
    std::convertible_to<Type,
      typename VertexAttribTypeTraits<Type>::storage_type>;
// clang-format on

/// @brief Vertex attribute range.
template<class Range>
concept vertex_attrib_type_range = std::ranges::range<Range> &&
    vertex_attrib_valid_type<std::ranges::range_value_t<Range>>;

/// @brief Type that would be stored in an OpenGL buffer.
template<class Type>
using vertex_attrib_storage_type_t =
    typename VertexAttribTypeTraits<Type>::storage_type;

/// @brief OpenGL type enumeration.
template<class Type>
inline constexpr GLenum
    vertex_attrib_type = VertexAttribTypeTraits<Type>::type();

/// @brief Number of the components in type.
template<class Type>
inline constexpr size_t
    vertex_attrib_length = VertexAttribTypeTraits<Type>::length();

#define MAKE_SCALAR_BASE_TYPE_(GLtype, GL_TYPE)              \
  template<>                                                 \
  struct VertexAttribTypeTraits<GLtype> {                    \
    using storage_type = GLtype;                             \
    [[nodiscard]] static constexpr GLenum type() noexcept {  \
      return GL_TYPE;                                        \
    }                                                        \
    [[nodiscard]] static constexpr GLint length() noexcept { \
      return 1;                                              \
    }                                                        \
  };

#define MAKE_SCALAR_CAST_TYPE_(constrain, GLtype)            \
  template<constrain Type>                                   \
  struct VertexAttribTypeTraits<Type> {                      \
    using storage_type = GLtype;                             \
    [[nodiscard]] static constexpr GLenum type() noexcept {  \
      return vertex_attrib_type<GLtype>;                     \
    }                                                        \
    [[nodiscard]] static constexpr GLint length() noexcept { \
      return 1;                                              \
    }                                                        \
  };

// Scalar integral types.
MAKE_SCALAR_BASE_TYPE_(GLbyte, GL_BYTE)
MAKE_SCALAR_BASE_TYPE_(GLubyte, GL_UNSIGNED_BYTE)
MAKE_SCALAR_BASE_TYPE_(GLshort, GL_SHORT)
MAKE_SCALAR_BASE_TYPE_(GLushort, GL_UNSIGNED_SHORT)
MAKE_SCALAR_BASE_TYPE_(GLint, GL_INT)
MAKE_SCALAR_BASE_TYPE_(GLuint, GL_UNSIGNED_INT)
MAKE_SCALAR_CAST_TYPE_(std::signed_integral, GLint)
MAKE_SCALAR_CAST_TYPE_(std::unsigned_integral, GLuint)

// Scalar floating-point types.
MAKE_SCALAR_BASE_TYPE_(GLfloat, GL_FLOAT)
MAKE_SCALAR_BASE_TYPE_(GLdouble, GL_DOUBLE)
MAKE_SCALAR_CAST_TYPE_(std::floating_point, GLdouble)

#undef MAKE_SCALAR_BASE_TYPE_
#undef MAKE_SCALAR_CAST_TYPE_

// GLM vector type.
template<vertex_attrib_valid_type Type, //
         glm::length_t Length, glm::qualifier Qualifier>
struct VertexAttribTypeTraits<glm::vec<Length, Type, Qualifier>> {
  using storage_type = glm::vec<Length, vertex_attrib_storage_type_t<Type>,
                                glm::qualifier::packed_highp>;
  [[nodiscard]] static constexpr GLenum type() noexcept {
    return vertex_attrib_type<Type>;
  }
  [[nodiscard]] static constexpr GLint length() noexcept {
    return static_cast<GLint>(Length) * vertex_attrib_length<Type>;
  }
};

// GLM matrix type.
template<vertex_attrib_valid_type Type, //
         glm::length_t Rows, glm::length_t Cols, glm::qualifier Qualifier>
struct VertexAttribTypeTraits<glm::mat<Rows, Cols, Type, Qualifier>> {
  using storage_type = glm::mat<Rows, Cols, vertex_attrib_storage_type_t<Type>,
                                glm::qualifier::packed_highp>;
  [[nodiscard]] static constexpr GLenum type() noexcept {
    return vertex_attrib_type<Type>;
  }
  [[nodiscard]] static constexpr GLint length() noexcept {
    return static_cast<GLint>(Rows * Cols) * vertex_attrib_length<Type>;
  }
};

/// @brief OpenGL host buffer.
class HostBuffer {
private:

  GLenum value_type_{};
  GLint value_length_{};
  GLsizei value_stride_{};
  std::vector<std::byte> host_data_{};

public:

  /// @brief Construct an empty host buffer.
  HostBuffer() = default;

  /// @brief Construct a host buffer with values.
  template<vertex_attrib_type_range Range>
  explicit HostBuffer(Range&& values) {
    assign(std::forward<Range>(values));
  }

  /// @brief Buffer value type.
  [[nodiscard]] GLenum value_type() const noexcept {
    return value_type_;
  }

  /// @brief Buffer value length.
  [[nodiscard]] GLint value_length() const noexcept {
    return value_length_;
  }

  /// @brief Buffer value stride.
  [[nodiscard]] GLsizei value_stride() const noexcept {
    return value_stride_;
  }

  /// @brief Buffer size in bytes.
  [[nodiscard]] size_t size() const noexcept {
    return host_data_.size();
  }

  /// @brief Pointer to the host buffer data.
  [[nodiscard]] const std::byte* data() const noexcept {
    return host_data_.data();
  }

  /// @brief Buffer size in values.
  [[nodiscard]] GLsizei num_values() const noexcept {
    return host_data_.size() / value_stride_;
  }

  /// @brief Assign the buffer values.
  template<vertex_attrib_type_range Range>
  void assign(Range&& values) {
    using Type = std::ranges::range_value_t<Range>;
    using StorageType = vertex_attrib_storage_type_t<Type>;
    // Assign the type properties.
    value_type_ = vertex_attrib_type<Type>;
    value_length_ = vertex_attrib_length<Type>;
    value_stride_ = sizeof(StorageType);
    // Copy range bytes into the host buffer.
    host_data_.clear();
    if constexpr (std::ranges::sized_range<Range>) {
      host_data_.reserve(sizeof(StorageType) * values.size());
    }
    std::ranges::copy( //
        values | std::views::transform([](const Type& value) {
          union {
            StorageType value;
            std::array<std::byte, sizeof(StorageType)> value_bytes;
          } union_cast = {.value = static_cast<StorageType>(value)};
          return union_cast.value_bytes;
        }) | std::views::join,
        std::back_inserter(host_data_));
  }

}; // class HostBuffer

/// @brief OpenGL host-device buffer.
class HostDeviceBuffer final : public HostBuffer, public DeviceBuffer {
public:

  /// @brief Construct an empty host-device buffer.
  HostDeviceBuffer() = default;

  /// @brief Move-construct a host-device buffer.
  HostDeviceBuffer(HostDeviceBuffer&&) = default;

  /// @brief Move-assign the host-device buffer.
  HostDeviceBuffer& operator=(HostDeviceBuffer&&) = default;

  /// @brief Construct a host-device buffer.
  template<vertex_attrib_type_range Range>
  explicit HostDeviceBuffer(Range&& values) {
    assign(std::forward<Range>(values));
  }

  /// @brief Assign the buffer values.
  template<vertex_attrib_type_range Range>
  void assign(Range&& values) {
    HostBuffer::assign(values);
    glBindBuffer(GL_ARRAY_BUFFER, id());
    glBufferData(GL_ARRAY_BUFFER, size(), data(), GL_STATIC_DRAW);
  }

}; // class HostDeviceBuffer

/// @brief OpenGL mesh.
class Mesh final {
private:

  VertexArray vertex_array_{};
  std::vector<const HostDeviceBuffer*> p_vertex_arrays_{};
  const HostDeviceBuffer* p_index_array_{};

public:

  /// @brief Push back the vertex array.
  void push_vertex_array(const HostDeviceBuffer& vertex_array) {
    STORM_INFO_("Pushing the vertex array: "
                "value_type = {:#x}, value_length = {}, num_values = {}.",
                vertex_array.value_type(), vertex_array.value_length(),
                vertex_array.num_values());
    p_vertex_arrays_.push_back(&vertex_array);
  }

  /// @brief Set the index array.
  void set_index_array(const HostDeviceBuffer& index_array) {
    STORM_INFO_("Setting the index array: "
                "value_type = {:#x}, num_values = {}.",
                index_array.value_type(), index_array.num_values());
    p_index_array_ = &index_array;
  }

  /// @brief Build the mesh.
  void build() {
    STORM_INFO_("Building the mesh: "
                "it has {} vertex buffers and {} index buffers",
                p_vertex_arrays_.size(), p_index_array_ != nullptr ? 1 : 0);
    BindVertexArray bind_vertex_array{vertex_array_};
    if (p_index_array_ != nullptr) {
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, *p_index_array_);
    }
    for (size_t i = 0; i < p_vertex_arrays_.size(); ++i) {
      const auto index = static_cast<GLuint>(i);
      glBindBuffer(GL_ARRAY_BUFFER, *p_vertex_arrays_[i]);
      glVertexAttribPointer(index, //
                            p_vertex_arrays_[i]->value_length(),
                            p_vertex_arrays_[i]->value_type(),
                            /*normalized*/ GL_FALSE, //
                            p_vertex_arrays_[i]->value_stride(),
                            /*offset*/ nullptr);
      glEnableVertexAttribArray(index);
    }
  }

  /// @brief Draw the mesh.
  void draw(GLenum mode = GL_TRIANGLES) const {
    STORM_ASSERT_(mode == GL_LINES || mode == GL_TRIANGLES,
                  "Unsupported draw mode {:#x}!", mode);
    BindVertexArray bind_vertex_array{vertex_array_};
    if (p_index_array_ != nullptr) {
      glDrawElements(mode, //
                     p_index_array_->num_values(), p_index_array_->value_type(),
                     /*indices*/ nullptr);
    } else {
      glDrawArrays(mode, //
                   /*first*/ 0, p_vertex_arrays_.front()->num_values());
    }
  }

}; // class Mesh

} // namespace Storm::Vulture::gl
