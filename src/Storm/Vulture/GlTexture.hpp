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
  texture2D = GL_TEXTURE_2D,
  texture_buffer = GL_TEXTURE_BUFFER,
}; // enum class TextureTarget

/// @brief OpenGL pixel description type.
struct pixel_desc_t {
  /// @brief Pixel integral format.
  GLenum internal_format;
  /// @brief Pixel format.
  GLenum format;
  /// @brief Pixel component type.
  GLenum type;
}; // struct pixel_desc_t

/// @brief OpenGL default pixel description.
/// Copied straightly from the OpenGL ES 3.0's `glTexImage2D`
/// sized internal format table.
template<class Pixel>
inline constexpr auto default_pixel_desc_v = []() -> pixel_desc_t {
  if constexpr (std::same_as<Pixel, GLbyte>)
    return {GL_R8I, GL_RED_INTEGER, GL_BYTE};
  if constexpr (std::same_as<Pixel, GLubyte>)
    return {GL_R8UI, GL_RED_INTEGER, GL_UNSIGNED_BYTE};
  if constexpr (std::same_as<Pixel, GLshort>)
    return {GL_R16I, GL_RED_INTEGER, GL_SHORT};
  if constexpr (std::same_as<Pixel, GLushort>)
    return {GL_R16UI, GL_RED_INTEGER, GL_UNSIGNED_SHORT};
  if constexpr (std::same_as<Pixel, GLint>)
    return {GL_R32I, GL_RED_INTEGER, GL_INT};
  if constexpr (std::same_as<Pixel, GLuint>)
    return {GL_R32UI, GL_RED_INTEGER, GL_UNSIGNED_INT};
  if constexpr (std::same_as<Pixel, GLfloat>)
    return {GL_R32F, GL_RED, GL_FLOAT};

  if constexpr (std::same_as<Pixel, glm::tvec2<GLbyte, glm::packed>>)
    return {GL_RG8I, GL_RG_INTEGER, GL_BYTE};
  if constexpr (std::same_as<Pixel, glm::tvec2<GLubyte, glm::packed>>)
    return {GL_RG8UI, GL_RG_INTEGER, GL_UNSIGNED_BYTE};
  if constexpr (std::same_as<Pixel, glm::tvec2<GLshort, glm::packed>>)
    return {GL_RG16I, GL_RG_INTEGER, GL_SHORT};
  if constexpr (std::same_as<Pixel, glm::tvec2<GLushort, glm::packed>>)
    return {GL_RG16UI, GL_RG_INTEGER, GL_UNSIGNED_SHORT};
  if constexpr (std::same_as<Pixel, glm::tvec2<GLint, glm::packed>>)
    return {GL_RG32I, GL_RG_INTEGER, GL_INT};
  if constexpr (std::same_as<Pixel, glm::tvec2<GLuint, glm::packed>>)
    return {GL_RG32UI, GL_RG_INTEGER, GL_UNSIGNED_INT};
  if constexpr (std::same_as<Pixel, glm::tvec2<GLfloat, glm::packed>>)
    return {GL_RG32F, GL_RG, GL_FLOAT};

  if constexpr (std::same_as<Pixel, glm::tvec3<GLbyte, glm::packed>>)
    return {GL_RGB8I, GL_RGB_INTEGER, GL_BYTE};
  if constexpr (std::same_as<Pixel, glm::tvec3<GLubyte, glm::packed>>)
    return {GL_RGB8UI, GL_RGB_INTEGER, GL_UNSIGNED_BYTE};
  if constexpr (std::same_as<Pixel, glm::tvec3<GLshort, glm::packed>>)
    return {GL_RGB16I, GL_RGB_INTEGER, GL_SHORT};
  if constexpr (std::same_as<Pixel, glm::tvec3<GLushort, glm::packed>>)
    return {GL_RGB16UI, GL_RGB_INTEGER, GL_UNSIGNED_SHORT};
  if constexpr (std::same_as<Pixel, glm::tvec3<GLint, glm::packed>>)
    return {GL_RGB32I, GL_RGB_INTEGER, GL_INT};
  if constexpr (std::same_as<Pixel, glm::tvec3<GLuint, glm::packed>>)
    return {GL_RGB32UI, GL_RGB_INTEGER, GL_UNSIGNED_INT};
  if constexpr (std::same_as<Pixel, glm::tvec3<GLfloat, glm::packed>>)
    return {GL_RGB32F, GL_RGB, GL_FLOAT};

  if constexpr (std::same_as<Pixel, glm::tvec4<GLbyte, glm::packed>>)
    return {GL_RGBA8I, GL_RGBA_INTEGER, GL_BYTE};
  if constexpr (std::same_as<Pixel, glm::tvec4<GLubyte, glm::packed>>)
    return {GL_RGBA8UI, GL_RGBA_INTEGER, GL_UNSIGNED_BYTE};
  if constexpr (std::same_as<Pixel, glm::tvec4<GLshort, glm::packed>>)
    return {GL_RGBA16I, GL_RGBA_INTEGER, GL_SHORT};
  if constexpr (std::same_as<Pixel, glm::tvec4<GLushort, glm::packed>>)
    return {GL_RGBA16UI, GL_RGBA_INTEGER, GL_UNSIGNED_SHORT};
  if constexpr (std::same_as<Pixel, glm::tvec4<GLint, glm::packed>>)
    return {GL_RGBA32I, GL_RGBA_INTEGER, GL_INT};
  if constexpr (std::same_as<Pixel, glm::tvec4<GLuint, glm::packed>>)
    return {GL_RGBA32UI, GL_RGBA_INTEGER, GL_UNSIGNED_INT};
  if constexpr (std::same_as<Pixel, glm::tvec4<GLfloat, glm::packed>>)
    return {GL_RGBA32F, GL_RGBA, GL_FLOAT};

  return {};
}();

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

/// @brief OpenGL 2D texture.
template<class Pixel, pixel_desc_t PixelDesc = default_pixel_desc_v<Pixel>>
class Texture2D : public Texture {
public:

  /// @brief Construct a 2D texture.
  Texture2D() = default;

  /// @brief Bind the texture buffer to @p slot.
  void bind(size_t slot) const {
    Texture::bind(TextureTarget::texture2D, slot);
  }

  /// @brief Assign the texture width and height.
  void assign(GLsizei width, GLsizei height) {
    Texture::bind(TextureTarget::texture2D);
    glTexImage2D(GL_TEXTURE_2D,
                 /*level*/ 0, PixelDesc.internal_format, width, height,
                 /*border*/ 0, PixelDesc.format, PixelDesc.type,
                 /*data*/ nullptr);
  }

}; // class Texture2D

/// @brief OpenGL texture buffer.
class TextureBuffer : public Texture {
public:

  /// @brief Construct a texture.
  // clang-format off
  template<class Pixel>
  TextureBuffer(const Buffer<Pixel>& buffer) {
    // clang-format on
    assign(buffer);
  }

  /// @brief Bind the texture buffer to @p slot.
  void bind(size_t slot) const {
    Texture::bind(TextureTarget::texture_buffer, slot);
  }

  /// @brief Assign the texture buffer underlying @p buffer.
  template<class Pixel, pixel_desc_t PixelDesc = default_pixel_desc_v<Pixel>>
  void assign(const Buffer<Pixel>& buffer) {
    Texture::bind(TextureTarget::texture_buffer);
    glTexBuffer(GL_TEXTURE_BUFFER, PixelDesc.internal_format, buffer);
  }

}; // class TextureBuffer

} // namespace Storm::Vulture::gl
