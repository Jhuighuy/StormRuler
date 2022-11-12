// Copyright © 2020 - 2023 Oleg Butakov
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

#include <Storm/Utils/Meta.hpp>

#include <algorithm>
#include <concepts>
#include <ranges>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include <GL/glew.h>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

// A macro, used to spawn shaders in '.glsl' files.
#define STORM_VULTURE_SHADER_(type, source)              \
  []() {                                                 \
    return ::Storm::Vulture::gl::Shader{                 \
        ::Storm::Vulture::gl::ShaderType::type, source}; \
  }(),

namespace Storm::Vulture::gl {

/// @brief OpenGL shader type.
enum class ShaderType : GLenum {
  vertex = GL_VERTEX_SHADER,
  geometry = GL_GEOMETRY_SHADER,
  fragment = GL_FRAGMENT_SHADER,
}; // enum class ShaderType

/// @brief OpenGL shader.
class Shader final : detail_::noncopyable_ {
private:

  GLuint shader_id_;

public:

  /// @brief Construct a shader.
  explicit Shader(ShaderType type) {
    shader_id_ = glCreateShader(static_cast<GLenum>(type));
  }
  /// @brief Construct a shader with @p source.
  Shader(ShaderType type, std::string_view source) : Shader{type} {
    assign(source);
  }

  /// @brief Move-construct a shader.
  Shader(Shader&&) = default;
  /// @brief Move-assign the shader.
  Shader& operator=(Shader&&) = default;

  /// @brief Destruct the shader.
  ~Shader() {
    glDeleteShader(shader_id_);
  }

  /// @brief Cast to shader ID.
  [[nodiscard]] constexpr operator GLuint() const noexcept {
    return shader_id_;
  }

  /// @brief Load the shader of type @p type from @p source.
  void assign(std::string_view source) {
    // Upload the shader source code.
    const auto source_data = source.data();
    const auto source_size = static_cast<GLint>(source.size());
    glShaderSource(shader_id_, 1, &source_data, &source_size);

    // Try to compile the shader and check result.
    glCompileShader(shader_id_);
    GLint status;
    glGetShaderiv(shader_id_, GL_COMPILE_STATUS, &status);
    std::string info_log(1024, '\0');
    GLsizei info_log_length;
    glGetShaderInfoLog(shader_id_, static_cast<GLsizei>(info_log.size()),
                       &info_log_length, info_log.data());
    if (status != GL_TRUE) {
      STORM_THROW_GL_("Failed to compile shader: {}", info_log.c_str());
    }
    if (info_log_length != 0) {
      STORM_WARNING_("Shader compilation message: {}", info_log.c_str());
    }
  }

}; // class Shader

// -----------------------------------------------------------------------------

/// @brief OpenGL shader program.
class Program final : detail_::noncopyable_ {
private:

  GLuint program_id_ = 0;

public:

  /// @brief Construct a program.
  Program() {
    program_id_ = glCreateProgram();
  }
  /// @brief Construct a program with @p shaders.
  template<std::same_as<Shader>... Shader>
  explicit Program(const Shader&... shaders) : Program() {
    assign(shaders...);
  }

  /// @brief Move-construct a program.
  Program(Program&&) = default;
  /// @brief Move-assign the program.
  Program& operator=(Program&&) = default;

  /// @brief Destruct the program.
  ~Program() noexcept {
    glDeleteProgram(program_id_);
  }

  /// @brief Cast to shader ID.
  [[nodiscard]] constexpr operator GLuint() const noexcept {
    return program_id_;
  }

  /// @brief Load the program with @p shaders.
  template<std::same_as<Shader>... Shader>
  void assign(const Shader&... shaders) {
    // Attach the shaders.
    (glAttachShader(program_id_, shaders), ...);

    // Try to link the program and check result.
    glLinkProgram(program_id_);
    GLint status;
    glGetProgramiv(program_id_, GL_LINK_STATUS, &status);
    std::string info_log(1024, '\0');
    GLsizei info_log_length;
    glGetProgramInfoLog(program_id_, static_cast<GLsizei>(info_log.size()),
                        &info_log_length, info_log.data());
    if (status != GL_TRUE) {
      STORM_THROW_GL_("Failed to link program: {}", info_log.c_str());
    }
    if (info_log_length != 0) {
      STORM_WARNING_("Program linking message: {}", info_log.c_str());
    }
  }

  /// @brief Get program uniform location.
  /// @{
  [[nodiscard]] GLint operator[](const std::string& name) const {
    return (*this)[name.c_str()];
  }
  [[nodiscard]] GLint operator[](const char* name) const {
    STORM_ASSERT_(name != nullptr, "Invalid uniform name!");
    return glGetUniformLocation(program_id_, name);
  }
  /// @}

}; // class Program

/// @brief Binder for an OpenGL program.
class BindProgram final {
public:

  /// @brief Bind the @p program.
  explicit BindProgram(const Program& program) {
    STORM_ASSERT_(program != 0, "Invalid program!");
    glUseProgram(program);
  }

  /// @brief Set the scalar uniform value.
  /// @{
  void set_uniform(GLint location, GLfloat value) {
    STORM_ASSERT_(location >= 0, "Invalid location!");
    glUniform1f(location, value);
  }
  void set_uniform(GLint location, GLint value) {
    STORM_ASSERT_(location >= 0, "Invalid location!");
    glUniform1i(location, value);
  }
  void set_uniform(GLint location, GLuint value) {
    STORM_ASSERT_(location >= 0, "Invalid location!");
    glUniform1ui(location, value);
  }
  /// @}

  /// @brief Set the vector uniform value.
  /// @{
  template<glm::length_t Length, glm::qualifier Qualifier>
  void set_uniform(GLint location,
                   const glm::vec<Length, GLfloat, Qualifier>& value) {
    STORM_ASSERT_(location >= 0, "Invalid location!");
    if constexpr (Length == 2) {
      glUniform2f(location, value.x, value.y);
    } else if constexpr (Length == 3) {
      glUniform3f(location, value.x, value.y, value.z);
    } else if constexpr (Length == 4) {
      glUniform4f(location, value.x, value.y, value.z, value.w);
    } else {
      static_assert(
          meta::always_false<std::integral_constant<glm::length_t, Length>>);
    }
  }
  template<glm::length_t Length, glm::qualifier Qualifier>
  void set_uniform(GLint location,
                   const glm::vec<Length, GLint, Qualifier>& value) {
    STORM_ASSERT_(location >= 0, "Invalid location!");
    if constexpr (Length == 2) {
      glUniform2i(location, value.x, value.y);
    } else if constexpr (Length == 3) {
      glUniform3i(location, value.x, value.y, value.z);
    } else if constexpr (Length == 4) {
      glUniform4i(location, value.x, value.y, value.z, value.w);
    } else {
      static_assert(
          meta::always_false<std::integral_constant<glm::length_t, Length>>);
    }
  }
  template<glm::length_t Length, glm::qualifier Qualifier>
  void set_uniform(GLint location,
                   const glm::vec<Length, GLuint, Qualifier>& value) {
    STORM_ASSERT_(location >= 0, "Invalid location!");
    if constexpr (Length == 2) {
      glUniform2ui(location, value.x, value.y);
    } else if constexpr (Length == 3) {
      glUniform3ui(location, value.x, value.y, value.z);
    } else if constexpr (Length == 4) {
      glUniform4ui(location, value.x, value.y, value.z, value.w);
    } else {
      static_assert(
          meta::always_false<std::integral_constant<glm::length_t, Length>>);
    }
  }
  /// @}

  /// @brief Set the matrix uniform value.
  template<glm::length_t Rows, glm::length_t Cols, glm::qualifier Qualifier>
  void set_uniform(GLint location,
                   const glm::mat<Rows, Cols, GLfloat, Qualifier>& value) {
    STORM_ASSERT_(location >= 0, "Invalid location!");
    const auto packed_value =
        static_cast<glm::mat<Rows, Cols, GLfloat, glm::packed>>(value);
    if constexpr (Rows == Cols) {
      if constexpr (Rows == 2) {
        glUniformMatrix2fv(location, 1, GL_FALSE, glm::value_ptr(packed_value));
      } else if constexpr (Rows == 3) {
        glUniformMatrix3fv(location, 1, GL_FALSE, glm::value_ptr(packed_value));
      } else if constexpr (Rows == 4) {
        glUniformMatrix4fv(location, 1, GL_FALSE, glm::value_ptr(packed_value));
      } else {
        static_assert(
            meta::always_false<std::integral_constant<glm::length_t, Rows>>);
      }
    } else if constexpr (Rows == 2 && Cols == 3) {
      glUniformMatrix2x3fv(location, 1, GL_FALSE, glm::value_ptr(packed_value));
    } else if constexpr (Rows == 3 && Cols == 2) {
      glUniformMatrix3x2fv(location, 1, GL_FALSE, glm::value_ptr(packed_value));
    } else if constexpr (Rows == 3 && Cols == 4) {
      glUniformMatrix3x4fv(location, 1, GL_FALSE, glm::value_ptr(packed_value));
    } else if constexpr (Rows == 4 && Cols == 3) {
      glUniformMatrix4x3fv(location, 1, GL_FALSE, glm::value_ptr(packed_value));
    } else {
      static_assert(
          meta::always_false<std::integer_sequence<glm::length_t, Rows, Cols>>);
    }
  }

#if 0
  /// @brief Set the scalar uniform array value.
  /// @{
  // clang-format off
  template<std::ranges::contiguous_range Range>
    requires std::same_as<std::ranges::range_value_t<Range>, GLfloat> 
  void set_uniform_array(GLint location, Range&& values) {
    // clang-format on
    STORM_ASSERT_(location >= 0, "Invalid location!");
    glUniform1fv(location, static_cast<GLsizei>(values.size()), values.data());
  }
  // clang-format off
  template<std::ranges::contiguous_range Range>
    requires std::same_as<std::ranges::range_value_t<Range>, GLint> 
  void set_uniform_array(GLint location, Range&& values) {
    // clang-format on
    STORM_ASSERT_(location >= 0, "Invalid location!");
    glUniform1iv(location, static_cast<GLsizei>(values.size()), values.data());
  }
  // clang-format off
  template<std::ranges::contiguous_range Range>
    requires std::same_as<std::ranges::range_value_t<Range>, GLuint> 
  void set_uniform_array(GLint location, Range&& values) {
    // clang-format on
    STORM_ASSERT_(location >= 0, "Invalid location!");
    glUniform1uiv(location, static_cast<GLsizei>(values.size()), values.data());
  }
  /// @}

  /// @brief Set the matrix uniform array value.
  // clang-format off
  template<std::ranges::sized_range Range>
    requires requires { 
      []<glm::length_t Rows, glm::length_t Cols, glm::qualifier Qualifier>(
          meta::type<glm::mat<Rows, Cols, GLfloat, Qualifier>>) {
      }(meta::type_v<std::ranges::range_value_t<Range>>);
    } 
  void set_uniform_array(GLint location, Range&& values) {
    // clang-format on
    STORM_ASSERT_(location >= 0, "Invalid location!");
    const auto size = static_cast<GLsizei>(values.size());
    [&]<glm::length_t Rows, glm::length_t Cols, glm::qualifier Qualifier>(
        meta::type<glm::mat<Rows, Cols, GLfloat, Qualifier>>) {
      decltype(auto) packed_values = [&]() -> decltype(auto) {
        if constexpr (std::ranges::contiguous_range<Range> &&
                      detail_::one_of_(Qualifier, glm::packed_lowp,
                                       glm::packed_mediump,
                                       glm::packed_highp)) {
          return std::forward<Range>(values);
        } else {
          return values | std::views::transform([](const auto& value) {
                   return static_cast< //
                       glm::mat<Rows, Cols, GLfloat, glm::packed_highp>>(value);
                 });
        }
      }();
      if constexpr (Rows == Cols) {
        if constexpr (Rows == 2) {
          glUniformMatrix2fv(location, size, GL_FALSE,
                             glm::value_ptr(*packed_values.data()));
        } else if constexpr (Rows == 3) {
          glUniformMatrix3fv(location, size, GL_FALSE,
                             glm::value_ptr(*packed_values.data()));
        } else if constexpr (Rows == 4) {
          glUniformMatrix4fv(location, size, GL_FALSE,
                             glm::value_ptr(*packed_values.data()));
        } else {
          static_assert(
              meta::always_false<std::integral_constant<glm::length_t, Rows>>);
        }
      } else if constexpr (Rows == 2 && Cols == 3) {
        glUniformMatrix2x3fv(location, size, GL_FALSE,
                             glm::value_ptr(*packed_values.data()));
      } else if constexpr (Rows == 3 && Cols == 2) {
        glUniformMatrix3x2fv(location, size, GL_FALSE,
                             glm::value_ptr(*packed_values.data()));
      } else if constexpr (Rows == 3 && Cols == 4) {
        glUniformMatrix3x4fv(location, size, GL_FALSE,
                             glm::value_ptr(*packed_values.data()));
      } else if constexpr (Rows == 4 && Cols == 3) {
        glUniformMatrix4x3fv(location, size, GL_FALSE,
                             glm::value_ptr(*packed_values.data()));
      } else {
        static_assert(meta::always_false<
                      std::integer_sequence<glm::length_t, Rows, Cols>>);
      }
    }(meta::type_v<std::ranges::range_value_t<Range>>);
  }
#endif

}; // class BindProgram

} // namespace Storm::Vulture::gl
