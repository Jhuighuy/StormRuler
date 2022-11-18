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

#include <GL/glew.h>

namespace Storm::Vulture::gl {

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
      default:
        STORM_DEBUG_("OpenGL: {} {} {:#x}: {}", //
                     debug_error_source, debug_type, id, message);
        break;
      case GL_DEBUG_SEVERITY_LOW_ARB:
        STORM_INFO_("OpenGL: {} {} {:#x}: {}", //
                    debug_error_source, debug_type, id, message);
        break;
      case GL_DEBUG_SEVERITY_MEDIUM_ARB:
        STORM_WARNING_("OpenGL: {} {} {:#x}: {}", //
                       debug_error_source, debug_type, id, message);
        break;
      case GL_DEBUG_SEVERITY_HIGH_ARB:
        STORM_ERROR_("OpenGL: {} {} {:#x}: {}", //
                     debug_error_source, debug_type, id, message);
        break;
    }
  }

}; // class DebugOutput

} // namespace Storm::Vulture::gl
