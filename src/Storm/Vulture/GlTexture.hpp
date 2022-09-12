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

#include <Storm/Vulture/GlBuffer.hpp>

#include <GL/glew.h>
#include <glm/glm.hpp>

namespace Storm::Vulture::gl {

/// @brief OpenGL texture binding target.
enum class TextureTarget : GLenum {
  texture_buffer = GL_TEXTURE_BUFFER,
}; // enum class TextureTarget

/// @brief OpenGL texture.
class Texture : detail_::noncopyable_ {
private:

  GLuint texture_id_;

public:

  /// @brief Construct a texture.
  Texture() {
    glGenTextures(1, &texture_id_);
  }

  /// @brief Move-construct a texture.
  Texture(Texture&&) = default;
  /// @brief Move-assign the texture.
  Texture& operator=(Texture&&) = default;

  /// @brief Destruct the texture.
  ~Texture() {
    glDeleteTextures(1, &texture_id_);
  }

  /// @brief Cast to texture ID.
  [[nodiscard]] constexpr operator GLuint() const noexcept {
    return texture_id_;
  }

  /// @brief Bind the texture to @p target.
  /// @{
  void bind(TextureTarget target) const {
    glBindTexture(static_cast<GLenum>(target), texture_id_);
  }
  void bind(TextureTarget target, size_t slot) const {
    glActiveTexture(GL_TEXTURE0 + slot);
    bind(target);
  }
  /// @}

}; // class Texture

/// @brief OpenGL texture buffer.
class TextureBuffer : public Texture {
public:

  /// @brief Construct a texture.
  template<class Type>
  TextureBuffer(const Buffer<Type>& buffer, GLenum intf) {
    Texture::bind(TextureTarget::texture_buffer);
    glTexBuffer(GL_TEXTURE_BUFFER, intf, buffer);
  }

  /// @brief Move-construct a texture.
  TextureBuffer(TextureBuffer&&) = default;
  /// @brief Move-assign the texture.
  TextureBuffer& operator=(TextureBuffer&&) = default;

  /// @brief Destruct the texture.
  ~TextureBuffer() = default;

  /// @brief Bind the texture buffer to @p slot.
  void bind(size_t slot) const {
    Texture::bind(TextureTarget::texture_buffer, slot);
  }

}; // class TextureBuffer

} // namespace Storm::Vulture::gl
