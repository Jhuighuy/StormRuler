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
#include <Storm/Vulture/GlShader.hpp>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/quaternion.hpp>

namespace Storm::Vulture::scene {

/// @brief Scene entity transform.
class Transform final {
private:

  const Transform* parent_{};
  glm::vec3 position_{};
  glm::quat rotation_{};
  glm::vec3 scale_{1.0, 1.0, 1.0};

public:

  /// @brief Construct the transform.
  constexpr Transform() = default;

  /// @brief Construct the transform with @p parent.
  constexpr explicit Transform(const Transform& parent) : parent_{&parent} {}

  /// @brief Transform position.
  [[nodiscard]] constexpr const glm::vec3& position() const noexcept {
    return position_;
  }

  /// @brief Transform rotation.
  [[nodiscard]] constexpr const glm::quat& rotation() const noexcept {
    return rotation_;
  }

  /// @brief Transform scale.
  [[nodiscard]] constexpr const glm::vec3& scale() const noexcept {
    return scale_;
  }

  /// @brief Translate the transform.
  constexpr void translate(const glm::vec3& delta) noexcept {
    position_ += delta;
  }

  /// @brief Rotate the transform.
  constexpr void rotate(const glm::vec3& delta_degrees) noexcept {
    rotation_ *= glm::quat{glm::radians(delta_degrees)};
  }

  /// @brief Rescale the transform.
  /// @{
  constexpr void rescale(float factor) noexcept {
    scale_ *= factor;
  }
  constexpr void rescale(const glm::vec3& factor) noexcept {
    scale_ *= factor;
  }
  /// @}

  /// @brief Transform model matrix.
  [[nodiscard]] /*constexpr*/ glm::mat4 //
  model_matrix() const noexcept {
    auto model_matrix = glm::translate(glm::mat4(1.0f), position_) *
                        glm::toMat4(rotation_) *
                        glm::scale(glm::mat4(1.0f), scale_);
    if (parent_ != nullptr) {
      model_matrix = parent_->model_matrix() * model_matrix;
    }
    return model_matrix;
  }

}; // class Transform

/// @brief Scene camera.
class Camera final {
private:

  Transform transform_{};
  glm::mat4 projection_matrix_{};

public:

  /// @brief Construct the camera.
  constexpr Camera() = default;

  /// @brief Construct the camera with @p parent.
  constexpr explicit Camera(const Transform& parent) : transform_{parent} {}

  /// @brief Camera transform.
  /// @{
  [[nodiscard]] Transform& transform() noexcept {
    return transform_;
  }
  [[nodiscard]] const Transform& transform() const noexcept {
    return transform_;
  }
  /// @}

  /// @brief Set perspective projection matrix.
  void set_perspective(float aspect_ratio, float fov_degrees = 60.0f, //
                       float near = 0.01f, float far = 100.0f) noexcept {
    projection_matrix_ = glm::perspective(fov_degrees, aspect_ratio, near, far);
  }

  /// @brief Camera view-projection matrix.
  [[nodiscard]] /*constexpr*/ glm::mat4
  view_projection_matrix() const noexcept {
    const glm::mat4 view_matrix = glm::inverse(transform_.model_matrix());
    return projection_matrix_ * view_matrix;
  }

}; // class Camera

/// @brief Scene mesh renderer.
class MeshRenderer final {
private:

  Transform transform_{};
  gl::Mesh mesh_{};
  gl::Program program_{};

public:

  /// @brief Construct the mesh renderer.
  MeshRenderer() = default;

  /// @brief Construct the mesh renderer with @p parent.
  explicit MeshRenderer(const Transform& parent) : transform_{parent} {}

  /// @brief Mesh renderer transform.
  /// @{
  [[nodiscard]] Transform& transform() noexcept {
    return transform_;
  }
  [[nodiscard]] const Transform& transform() const noexcept {
    return transform_;
  }
  /// @}

  /// @brief Mesh renderer mesh.
  /// @{
  [[nodiscard]] gl::Mesh& mesh() noexcept {
    return mesh_;
  }
  [[nodiscard]] const gl::Mesh& mesh() const noexcept {
    return mesh_;
  }
  /// @}

  /// @brief Mesh renderer program.
  /// @{
  [[nodiscard]] gl::Program& program() noexcept {
    return program_;
  }
  [[nodiscard]] const gl::Program& program() const noexcept {
    return program_;
  }
  /// @}

  /// @brief Draw the mesh.
  /// @{
  void draw(const Camera& camera, const gl::Program& program,
            GLenum mode = GL_TRIANGLES) const {
    // gl::BindProgram bind_program{program};
    mesh_.draw(mode);
  }
  void draw(const Camera& camera, GLenum mode = GL_TRIANGLES) const {
    draw(camera, program_, mode);
  }
  /// @}

}; // class Camera

} // namespace Storm::Vulture::scene
