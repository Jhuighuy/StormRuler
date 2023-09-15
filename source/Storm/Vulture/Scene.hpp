// Copyright (C) 2020-2023 Oleg Butakov
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
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <Storm/Base.hpp>

#include <Storm/Vulture/GlShader.hpp>
#include <Storm/Vulture/GlVertexArray.hpp>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/quaternion.hpp>

namespace Storm::Vulture::scene {

/// @brief Scene entity transform.
class Transform final {
private:

  glm::vec3 _position{};
  glm::quat _rotation{glm::vec3{}}; // initializes zero rotation.
  glm::vec3 _scale{1.0, 1.0, 1.0};

public:

  /// @brief Construct the transform.
  constexpr Transform() = default;

  /// @brief Transform position.
  constexpr const glm::vec3& position() const noexcept {
    return _position;
  }

  /// @brief Transform rotation.
  /// @{
  constexpr const glm::quat& rotation() const noexcept {
    return _rotation;
  }
  glm::vec3 rotation_degrees() const noexcept {
    return glm::degrees(glm::eulerAngles(_rotation));
  }
  /// @}

  /// @brief Transform scale.
  constexpr const glm::vec3& scale() const noexcept {
    return _scale;
  }

  /// @brief Translate the transform.
  constexpr void translate(const glm::vec3& delta) noexcept {
    _position += delta;
  }

  /// @brief Rotate the transform.
  /// @{
  constexpr void rotate(const glm::quat& delta) noexcept {
    _rotation *= delta;
  }
  constexpr void rotate_degrees(const glm::vec3& delta_degrees) noexcept {
    _rotation *= glm::quat{glm::radians(delta_degrees)};
  }
  /// @}

  /// @brief Rescale the transform.
  /// @{
  constexpr void rescale(float factor) noexcept {
    _scale *= factor;
  }
  constexpr void rescale(const glm::vec3& factor) noexcept {
    _scale *= factor;
  }
  /// @}

  /// @brief Transform model matrix.
  /// @todo GLM's matmul operator is not constexpr?
  glm::mat4 model_matrix() const noexcept {
    auto model_matrix = glm::translate(glm::mat4(1.0f), _position) *
                        glm::toMat4(_rotation) *
                        glm::scale(glm::mat4(1.0f), _scale);
    return model_matrix;
  }

}; // class Transform

// -----------------------------------------------------------------------------

/// @brief Scene camera.
class Camera final {
private:

  float _near = 0.01f, _far = 0.99f;
  float _orbit = 1.0f;
  Transform _transform{};
  glm::mat4 _projection_matrix{};

public:

  /// @brief Construct the camera.
  constexpr Camera() = default;

  /// @brief Construct the camera with @p parent.
  constexpr explicit Camera(const Transform& parent) : _transform{parent} {}

  /// @brief Camera transform.
  /// @{
  constexpr Transform& transform() noexcept {
    return _transform;
  }
  constexpr const Transform& transform() const noexcept {
    return _transform;
  }
  /// @}

  /// @brief Camera orbit.
  /// @{
  constexpr float orbit() const noexcept {
    return _orbit;
  }
  constexpr void set_orbit(float orbit) noexcept {
    _orbit = std::clamp(orbit, _near, _far);
  }
  /// @}

  /// @brief Set perspective projection matrix.
  void set_perspective(float aspect_ratio, float fov_degrees = 60.0f,
                       float the_near = 0.001f,
                       float the_far = 1000.0f) noexcept {
    // Beware: GLM's documentations says that FOV is in degress,
    // but actually it is in radians!
    _projection_matrix = glm::perspective(glm::radians(fov_degrees),
                                          aspect_ratio, the_near, the_far);
    _near = the_near, _far = the_far;
  }

  /// @brief Set orthographic projection matrix.
  void set_ortographic(float aspect_ratio, float height = 1.0f,
                       float the_near = 0.001f,
                       float the_far = 1000.0f) noexcept {
    const float width = aspect_ratio * height;
    _projection_matrix =
        glm::ortho(-0.5f * width, +0.5f * width, //
                   -0.5f * height, +0.5f * height, the_near, the_far);
    _near = the_near, _far = the_far;
  }

  /// @brief Camera view-projection matrix.
  /// @todo GLM's matmul operator is not constexpr?
  glm::mat4 view_projection_matrix() const noexcept {
    const glm::mat4 view_matrix = glm::inverse(
        _transform.model_matrix() *
        glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f, _orbit)));
    return _projection_matrix * view_matrix;
  }

}; // class Camera

// -----------------------------------------------------------------------------

/// @brief Scene mesh renderer.
class MeshRenderer final {
private:

  Transform _transform{};
  // gl::Mesh _mesh{};
  gl::Program _program{};

public:

  /// @brief Construct the mesh renderer.
  MeshRenderer() = default;

  /// @brief Construct the mesh renderer with @p parent.
  explicit MeshRenderer(const Transform& parent) : _transform{parent} {}

  /// @brief Mesh renderer transform.
  /// @{
  constexpr Transform& transform() noexcept {
    return _transform;
  }
  constexpr const Transform& transform() const noexcept {
    return _transform;
  }
  /// @}

  /// @brief Mesh renderer mesh.
  /// @{
  // constexpr gl::Mesh& mesh() noexcept {
  //  return _mesh;
  //}
  // constexpr const gl::Mesh& mesh() const noexcept {
  //  return _mesh;
  //}
  /// @}

  /// @brief Mesh renderer program.
  /// @{
  constexpr gl::Program& program() noexcept {
    return _program;
  }
  constexpr const gl::Program& program() const noexcept {
    return _program;
  }
  /// @}

  /// @brief Draw the mesh.
  /// @{
  // void draw(const Camera& camera, const gl::Program& program,
  //           GLenum mode = GL_TRIANGLES) const {
  //   // gl::BindProgram bind_program{program};
  //   _mesh.draw(mode);
  // }
  // void draw(const Camera& camera, GLenum mode = GL_TRIANGLES) const {
  //   draw(camera, _program, mode);
  // }
  /// @}

}; // class Camera

} // namespace Storm::Vulture::scene
