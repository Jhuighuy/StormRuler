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

#include <Storm/Utils/Containers.hpp>

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef> // std::byte
#include <optional>
#include <ranges>
#include <string>
#include <string_view>
#include <vector>

#include <glm/glm.hpp>
// clang-format off
#include <GL/glew.h>
#include <GLFW/glfw3.h>
// clang-format on


namespace Storm::gl {

/// @brief RAII OpenGL debug output.
class DebugOutput {
public:

  /// @brief Enable OpenGL debug output.
  /// Sticks to GL_ARB_debug_output extension since debug output
  /// is inside OpenGL since 4.3, and we are using 3.3.
  DebugOutput() noexcept {
    if (GLEW_ARB_debug_output != GL_TRUE) { return; }
    glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS_ARB);
    glDebugMessageCallbackARB(&on_message_, nullptr);
    STORM_INFO_("OpenGL debug output enabled.");
  }

  /// @brief Disable OpenGL debug output.
  ~DebugOutput() noexcept {
    if (GLEW_ARB_debug_output != GL_TRUE) { return; }
    glDisable(GL_DEBUG_OUTPUT_SYNCHRONOUS_ARB);
    STORM_INFO_("OpenGL debug output disabled.");
  }

private:

  static void on_message_(GLenum source, GLenum type, GLuint id,
                          GLenum severity, [[maybe_unused]] GLsizei length,
                          const GLchar* message,
                          [[maybe_unused]] const void* user_param) {
    const char* debug_error_source;
    switch (source) {
      case GL_DEBUG_SOURCE_API_ARB: //
        debug_error_source = "API call";
        break;
      case GL_DEBUG_SOURCE_WINDOW_SYSTEM_ARB:
        debug_error_source = "window system API all";
        break;
      case GL_DEBUG_SOURCE_SHADER_COMPILER_ARB:
        debug_error_source = "shader compiler";
        break;
      case GL_DEBUG_SOURCE_THIRD_PARTY_ARB:
        debug_error_source = "third party API";
        break;
      case GL_DEBUG_SOURCE_APPLICATION_ARB:
        debug_error_source = "application";
        break;
      case GL_DEBUG_SOURCE_OTHER_ARB: //
        debug_error_source = "other";
        break;
      default: //
        debug_error_source = "unknown source";
        break;
    }

    const char* debug_type;
    switch (type) {
      default:
      case GL_DEBUG_TYPE_OTHER_ARB: //
        debug_type = "other issue";
        break;
      case GL_DEBUG_TYPE_ERROR_ARB: //
        debug_type = "error";
        break;
      case GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR_ARB:
        debug_type = "deprecated behavior";
        break;
      case GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR_ARB:
        debug_type = "undefined behavior";
        break;
      case GL_DEBUG_TYPE_PORTABILITY_ARB: //
        debug_type = "portability issue";
        break;
      case GL_DEBUG_TYPE_PERFORMANCE_ARB: //
        debug_type = "performance issue";
        break;
    }

    switch (severity) {
      case GL_DEBUG_SEVERITY_LOW_ARB:          //
        STORM_INFO_("OpenGL: {} {} {:#x}: {}", //
                    debug_error_source, debug_type, id, message);
        break;
      case GL_DEBUG_SEVERITY_MEDIUM_ARB:
        STORM_WARNING_("OpenGL: {} {} {:#x}: {}", //
                       debug_error_source, debug_type, id, message);
        break;
      default:
      case GL_DEBUG_SEVERITY_HIGH_ARB:
        STORM_ERROR_("OpenGL: {} {} {:#x}: {}", //
                     debug_error_source, debug_type, id, message);
        break;
    }
  }

}; // class DebugOutput

/// @brief OpenGL object.
class Object {
private:

  GLuint id_ = 0;

public:

  /// @brief Construct an object.
  Object() = default;

  /// @brief Move-construct an object.
  constexpr Object(Object&& other) noexcept
      : id_{std::exchange(other.id_, 0)} {}
  /// @brief Move-assign operator.
  constexpr Object& operator=(Object&& other) noexcept {
    id_ = std::exchange(other.id_, 0);
    return *this;
  }

  Object(const Object&) = delete;
  Object& operator=(const Object&) = delete;

  /// @brief Get the object ID.
  [[nodiscard]] GLuint id() const noexcept {
    return id_;
  }
  /// @brief Set a new ID.
  void set_id(GLuint id) noexcept {
    id_ = id;
  }

  /// @brief Cast to object ID operator.
  [[nodiscard]] operator GLuint() const noexcept {
    return id_;
  }

}; // class Object

/// @brief OpenGL buffer.
class Buffer final : public Object {
public:

  /// @brief Construct a buffer.
  Buffer() noexcept {
    GLuint buffer_id;
    glGenBuffers(1, &buffer_id);
    set_id(buffer_id);
  }

  /// @brief Move-construct a buffer.
  Buffer(Buffer&&) = default;

  /// @brief Move-assign the buffer.
  Buffer& operator=(Buffer&&) = default;

  /// @brief Destruct the buffer.
  ~Buffer() noexcept {
    const GLuint buffer_id = id();
    glDeleteBuffers(1, &buffer_id);
  }

}; // class Buffer

/// @brief OpenGL vertex array.
class VertexArray final : public Object {
public:

  /// @brief Construct a vertex array.
  VertexArray() noexcept {
    GLuint vertex_array_id;
    glGenVertexArrays(1, &vertex_array_id);
    set_id(vertex_array_id);
  }

  /// @brief Move-construct a vertex array.
  VertexArray(VertexArray&&) = default;

  /// @brief Move-assign the vertex array.
  VertexArray& operator=(VertexArray&&) = default;

  /// @brief Destruct the vertex array.
  ~VertexArray() noexcept {
    glBindVertexArray(0);
    const GLuint vertex_array_id = id();
    glDeleteVertexArrays(1, &vertex_array_id);
  }

}; // class VertexArray

/// @brief RAII binder of an OpenGL vertex array.
class BindVertexArray final {
public:

  /// @brief Bind the @p vertex_array.
  BindVertexArray(const VertexArray& vertex_array) noexcept {
    STORM_ASSERT_(vertex_array != 0, "Invalid vertex array!");
    glBindVertexArray(vertex_array);
  }

  /// @brief Unbind the vertex_array.
  ~BindVertexArray() noexcept {
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
MAKE_SCALAR_CAST_TYPE_(index, GLuint)

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

/// @brief OpenGL mesh.
class Mesh final {
private:

  struct HostBuffer {
    GLenum value_type{};
    GLint value_length{};
    GLsizei value_stride{};
    std::vector<std::byte> host_data{};
    [[nodiscard]] GLsizei num_values() const noexcept {
      return host_data.size() / value_stride;
    }
  };

  struct HostDeviceBuffer : HostBuffer {
    Buffer device_buffer{};
  };

  VertexArray vertex_array_{};
  std::optional<HostDeviceBuffer> element_hd_buffer_{};
  std::vector<HostDeviceBuffer> vertex_hd_buffers_{};

public:

  /// @brief Push back the vertex attributes range.
  // clang-format off
  template<std::ranges::range Range>
    requires vertex_attrib_valid_type<std::ranges::range_value_t<Range>>
  void push_attribs(Range&& values) {
    // clang-format on
    HostDeviceBuffer vertex_hd_buffer(
        make_host_buffer_(std::forward<Range>(values)));
    STORM_INFO_("Pushing the vertex index array: "
                "value_type = {:#x}, value_length = {}, num_values = {}.",
                vertex_hd_buffer.value_type, vertex_hd_buffer.value_length,
                vertex_hd_buffer.num_values());
    glBindBuffer(GL_ARRAY_BUFFER, vertex_hd_buffer.device_buffer);
    glBufferData(GL_ARRAY_BUFFER, vertex_hd_buffer.host_data.size(),
                 vertex_hd_buffer.host_data.data(), GL_STATIC_DRAW);
    vertex_hd_buffers_.push_back(std::move(vertex_hd_buffer));
  }

  /// @brief Set the element indices range.
  // clang-format off
  template<std::ranges::range Range>
    requires vertex_attrib_valid_type<std::ranges::range_value_t<Range>>
  void set_element_indices(Range&& indices) {
    // clang-format on
    HostDeviceBuffer element_hd_buffer(
        make_host_buffer_(std::forward<Range>(indices)));
    STORM_INFO_("Setting the element index array: "
                "value_type = {:#x}, num_values = {}.",
                element_hd_buffer.value_type, element_hd_buffer.num_values());
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, element_hd_buffer.device_buffer);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, element_hd_buffer.host_data.size(),
                 element_hd_buffer.host_data.data(), GL_STATIC_DRAW);
    element_hd_buffer_ = std::move(element_hd_buffer);
  }

  /// @brief Build the mesh.
  void build() noexcept {
    STORM_INFO_("Building the mesh: "
                "it has {} vertex buffers and {} index buffers",
                vertex_hd_buffers_.size(),
                element_hd_buffer_.has_value() ? 1 : 0);
    BindVertexArray bind_vertex_array{vertex_array_};
    if (element_hd_buffer_.has_value()) {
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, element_hd_buffer_->device_buffer);
    }
    for (size_t i = 0; i < vertex_hd_buffers_.size(); ++i) {
      const auto index = static_cast<GLuint>(i);
      glBindBuffer(GL_ARRAY_BUFFER, vertex_hd_buffers_[i].device_buffer);
      glVertexAttribPointer(index, //
                            vertex_hd_buffers_[i].value_length,
                            vertex_hd_buffers_[i].value_type,
                            /*normalized*/ GL_FALSE, //
                            vertex_hd_buffers_[i].value_stride,
                            /*offset*/ nullptr);
      glEnableVertexAttribArray(index);
    }
  }

  /// @brief Draw the mesh.
  void draw(GLenum mode = GL_TRIANGLES) noexcept {
    BindVertexArray bind_vertex_array{vertex_array_};
    if (element_hd_buffer_.has_value()) {
      glDrawElements(mode, //
                     element_hd_buffer_->num_values(),
                     element_hd_buffer_->value_type,
                     /*indices*/ nullptr);
    } else {
      glDrawArrays(mode, //
                   /*first*/ 0, vertex_hd_buffers_.front().num_values());
    }
  }

private:

  // clang-format off
  template<std::ranges::range Range>
    requires vertex_attrib_valid_type<std::ranges::range_value_t<Range>>
  [[nodiscard]] static HostBuffer make_host_buffer_(Range&& values) {
    // clang-format on
    using Type = std::ranges::range_value_t<Range>;
    using StorageType = vertex_attrib_storage_type_t<Type>;
    // Create the host buffer.
    HostBuffer buffer{
        .value_type = vertex_attrib_type<Type>,
        .value_length = vertex_attrib_length<Type>,
        .value_stride = sizeof(StorageType),
    };
    // Copy range bytes into the host buffer.
    if constexpr (std::ranges::sized_range<Range>) {
      buffer.host_data.reserve(sizeof(StorageType) * values.size());
    }
    std::ranges::copy( //
        values | std::views::transform([](const Type& value) {
          union {
            StorageType value;
            std::array<std::byte, sizeof(StorageType)> value_bytes;
          } union_cast = {.value = static_cast<StorageType>(value)};
          return union_cast.value_bytes;
        }) | std::views::join,
        std::back_inserter(buffer.host_data));
    return buffer;
  }

}; // class Mesh

/// @brief OpenGL shader.
class Shader final : public Object {
public:

  /// @brief Construct a shader.
  Shader() = default;

  /// @brief Move-construct a shader.
  Shader(Shader&&) = default;

  /// @brief Move-assign the shader.
  Shader& operator=(Shader&&) = default;

  /// @brief Destruct the shader.
  ~Shader() noexcept {
    unload();
  }

  /// @brief Load the shader of type @p shader_type from @p shader_source.
  [[nodiscard]] bool load(GLenum shader_type,
                          std::string_view shader_source) noexcept {
    // Create a shader.
    STORM_ASSERT_(shader_type == GL_VERTEX_SHADER ||
                      shader_type == GL_FRAGMENT_SHADER,
                  "Invalid shader type!");
    unload();
    set_id(glCreateShader(shader_type));

    // Upload the shader source code.
    const auto shader_source_data = shader_source.data();
    const auto shader_source_size = static_cast<GLint>(shader_source.size());
    glShaderSource(id(), 1, &shader_source_data, &shader_source_size);

    // Try to compile the shader and check result.
    glCompileShader(id());
    GLint status;
    glGetShaderiv(id(), GL_COMPILE_STATUS, &status);
    std::string info_log(1024, '\0');
    GLsizei info_log_length;
    glGetShaderInfoLog(id(), static_cast<GLsizei>(info_log.size()),
                       &info_log_length, info_log.data());
    if (status != GL_TRUE) {
      STORM_ERROR_("Failed to compile shader: {}", info_log.c_str());
      return false;
    }
    if (info_log_length != 0) {
      STORM_WARNING_("Shader compilation message: {}", info_log.c_str());
    }
    return true;
  }

  /// @brief Unload the shader.
  void unload() noexcept {
    glDeleteShader(id());
    set_id(0);
  }

}; // class Shader

/// @brief OpenGL shader program.
class Program final : public Object {
public:

  /// @brief Construct a program.
  Program() = default;

  /// @brief Move-construct a program.
  Program(Program&&) = default;

  /// @brief Move-assign the program.
  Program& operator=(Program&&) = default;

  /// @brief Destruct the program.
  ~Program() noexcept {
    unload();
  }

  /// @brief Load the program from @p vertex_shader @p fragment_shader.
  [[nodiscard]] bool load(const Shader& vertex_shader,
                          const Shader& fragment_shader) noexcept {
    // Create a program.
    unload();
    set_id(glCreateProgram());

    // Attach the vertex and fragment shaders.
    STORM_ASSERT_(vertex_shader != 0, "Invalid vertex shader!");
    glAttachShader(id(), vertex_shader);
    STORM_ASSERT_(fragment_shader != 0, "Invalid fragment shader!");
    glAttachShader(id(), fragment_shader);

    // Try to link the program and check result.
    glLinkProgram(id());
    GLint status;
    glGetProgramiv(id(), GL_LINK_STATUS, &status);
    std::string info_log(1024, '\0');
    GLsizei info_log_length;
    glGetProgramInfoLog(id(), static_cast<GLsizei>(info_log.size()),
                        &info_log_length, info_log.data());
    if (status != GL_TRUE) {
      STORM_ERROR_("Failed to link program: {}", info_log.c_str());
      return false;
    }
    if (info_log_length != 0) {
      STORM_WARNING_("Program linking message: {}", info_log.c_str());
    }
    return true;
  }

  /// @brief Unload the program.
  void unload() noexcept {
    glUseProgram(0);
    glDeleteProgram(id());
    set_id(0);
  }

}; // class Program

/// @brief RAII binder of an OpenGL program.
class BindProgram final {
public:

  /// @brief Bind the @p program.
  BindProgram(const Program& program) noexcept {
    STORM_ASSERT_(program != 0, "Invalid program!");
    glUseProgram(program);
  }

  /// @brief Unbind the program.
  ~BindProgram() noexcept {
    glUseProgram(0);
  }

}; // class BindProgram

} // namespace Storm::gl
