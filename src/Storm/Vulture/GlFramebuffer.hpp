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

#include <Storm/Vulture/GlTexture.hpp>

#include <array>

#include <GL/glew.h>

namespace Storm::Vulture::gl {

/// @brief OpenGL framebuffer.
class Framebuffer : detail_::noncopyable_ {
private:

  GLuint framebuffer_id_;

public:

  /// @brief Construct a framebuffer.
  Framebuffer() {
    glGenFramebuffers(1, &framebuffer_id_);
  }

  /// @brief Construct a framebuffer with @p textures.
  template<pixel... Pixels>
  explicit Framebuffer(const Texture2D<Pixels>&... textures) : Framebuffer{} {
    assign(textures...);
  }

  /// @brief Move-construct a framebuffer.
  Framebuffer(Framebuffer&&) = default;
  /// @brief Move-assign the framebuffer.
  Framebuffer& operator=(Framebuffer&&) = default;

  /// @brief Destruct the framebuffer.
  ~Framebuffer() {
    glDeleteFramebuffers(1, &framebuffer_id_);
  }

  /// @brief Cast to framebuffer ID.
  [[nodiscard]] constexpr operator GLuint() const noexcept {
    return framebuffer_id_;
  }

  /// @brief Build the framebuffer with @p textures.
  template<pixel... Pixels>
  void assign(const Texture2D<Pixels>&... textures) {
    glBindFramebuffer(GL_FRAMEBUFFER, framebuffer_id_);
    std::array<GLenum, sizeof...(textures)> attachments{};
    const auto attach_texture =
        [&]<class Pixel>(size_t index, const Texture2D<Pixel>& texture) {
          /// @todo Depth attachments!
          attachments[index] = GL_COLOR_ATTACHMENT0 + index;
          glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0 + index,
                               texture, /*level*/ 0);
        };
    size_t index = 0;
    (attach_texture(index++, textures), ...);
    glDrawBuffers(GLsizei{sizeof...(textures)}, attachments.data());
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
      STORM_THROW_GL_("Failed to initialize framebuffer!");
    }
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
  }

  /// @brief Draw to framebuffer.
  template<std::invocable DrawFn>
  void draw_into(DrawFn draw_fn) {
    glBindFramebuffer(GL_FRAMEBUFFER, framebuffer_id_);
    draw_fn();
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
  }

}; // class Framebuffer

} // namespace Storm::Vulture::gl
