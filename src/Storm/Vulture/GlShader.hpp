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

#include <string>
#include <string_view>

#include <GL/glew.h>
#include <glm/glm.hpp>

namespace Storm::Vulture::gl {

/// @brief OpenGL shader.
class Shader final : public IdHolder {
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
  void load(GLenum type, std::string_view source) {
    // Create a shader.
    STORM_ASSERT_(type == GL_VERTEX_SHADER || //
                      type == GL_GEOMETRY_SHADER || type == GL_FRAGMENT_SHADER,
                  "Unsupported shader type {:#x}!", type);
    unload();
    set_id(glCreateShader(type));

    // Upload the shader source code.
    const auto source_data = source.data();
    const auto source_size = static_cast<GLint>(source.size());
    glShaderSource(id(), 1, &source_data, &source_size);

    // Try to compile the shader and check result.
    glCompileShader(id());
    GLint status;
    glGetShaderiv(id(), GL_COMPILE_STATUS, &status);
    std::string info_log(1024, '\0');
    GLsizei info_log_length;
    glGetShaderInfoLog(id(), static_cast<GLsizei>(info_log.size()),
                       &info_log_length, info_log.data());
    if (status != GL_TRUE) {
      STORM_THROW_GL_("Failed to compile shader: {}", info_log.c_str());
    }
    if (info_log_length != 0) {
      STORM_WARNING_("Shader compilation message: {}", info_log.c_str());
    }
  }

  /// @brief Unload the shader.
  void unload() noexcept {
    glDeleteShader(id());
    set_id(0);
  }

}; // class Shader

/// @brief OpenGL shader program.
class Program final : public IdHolder {
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

  /// @brief Load the program from shaders.
  template<std::same_as<Shader>... Shader>
  void load(const Shader&... shaders) {
    // Create a program.
    unload();
    set_id(glCreateProgram());

    // Attach the shaders.
    const auto attach_shader = [&](const auto& shader) {
      glAttachShader(id(), shader);
    };
    (attach_shader(shaders), ...);

    // Try to link the program and check result.
    glLinkProgram(id());
    GLint status;
    glGetProgramiv(id(), GL_LINK_STATUS, &status);
    std::string info_log(1024, '\0');
    GLsizei info_log_length;
    glGetProgramInfoLog(id(), static_cast<GLsizei>(info_log.size()),
                        &info_log_length, info_log.data());
    if (status != GL_TRUE) {
      STORM_THROW_GL_("Failed to link program: {}", info_log.c_str());
    }
    if (info_log_length != 0) {
      STORM_WARNING_("Program linking message: {}", info_log.c_str());
    }
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
  explicit BindProgram(const Program& program) {
    STORM_ASSERT_(program != 0, "Invalid program!");
    glUseProgram(program);
  }

  /// @brief Unbind the program.
  ~BindProgram() {
    glUseProgram(0);
  }

}; // class BindProgram

} // namespace Storm::Vulture::gl
