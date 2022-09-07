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

#include <Storm/Mallard/Entity.hpp>
#include <Storm/Mallard/Mesh.hpp>

#include <Storm/Vulture/GlBuffer.hpp>
#include <Storm/Vulture/GlDebug.hpp>
#include <Storm/Vulture/GlShader.hpp>
#include <Storm/Vulture/GlWindow.hpp>
#include <Storm/Vulture/Scene.hpp>

#include <concepts>
#include <cstddef> // std::byte
#include <optional>
#include <ranges>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#include <GL/glew.h>
#include <glm/glm.hpp>

namespace Storm::Vulture {

template<mesh Mesh>
void visualize_mesh(const Mesh& mesh) {
  gl::Framework framework{};

  constexpr static const char* window_title = "Strom::Vulture Visualizer";
  constexpr static size_t window_width = 1600;
  constexpr static size_t window_height = 900;
  gl::Window window{};
  window.load(framework, window_title, window_width, window_height);
  gl::BindWindow bind_window{window};
  glViewport(0, 0, window_width, window_height);

  // Initialize.
  gl::DebugOutput debug_output{};
  gl::Shader vertex_shader{};
  vertex_shader.load(GL_VERTEX_SHADER, R"(
      #version 330 core
      layout(location = 0) in vec3 positionMS;
      uniform mat4 view_projection;
      void main() {
        gl_Position = vec4(positionMS - vec3(1.0, 0.5, 0.0), 1.0);
        gl_Position = view_projection * gl_Position;
      }
    )");

  gl::Shader fragment_shader{};
  fragment_shader.load(GL_FRAGMENT_SHADER, R"(
      #version 330 core
      out vec4 color;
      uniform vec3 col;
      void main() {
        color = vec4(col, 1.0);
      }
    )");

  gl::Program program;
  program.load(std::move(vertex_shader), std::move(fragment_shader));

  const gl::HostDeviceBuffer mesh_nodes(
      nodes(mesh) | std::views::transform([](NodeView<const Mesh> node) {
        return node.position();
      }));
  const gl::HostDeviceBuffer mesh_edges(
      edges(mesh) | std::views::transform([](EdgeView<const Mesh> edge) {
        return edge.nodes() |
               std::views::transform([](NodeView<const Mesh> node) {
                 return static_cast<size_t>(node.index());
               });
      }) |
      std::views::join);
  const gl::HostDeviceBuffer mesh_cells(
      cells(mesh) | std::views::transform([](CellView<const Mesh> cell) {
        return cell.nodes() |
               std::views::transform([](NodeView<const Mesh> node) {
                 return static_cast<size_t>(node.index());
               });
      }) |
      std::views::join);

  scene::MeshRenderer edge_renderer{};
  edge_renderer.mesh().push_vertex_array(mesh_nodes);
  edge_renderer.mesh().set_index_array(mesh_edges);
  edge_renderer.mesh().build();

  scene::MeshRenderer cell_renderer{};
  cell_renderer.mesh().push_vertex_array(mesh_nodes);
  cell_renderer.mesh().set_index_array(mesh_cells);
  cell_renderer.mesh().build();

  // Setup camera.
  scene::Camera camera{};
  camera.transform().translate(glm::vec3{0.0f, 0.0f, 1.0f});
  camera.set_perspective(16.0 / 9.0);
  window.on_key_down(
      [&]() {
        camera.transform().translate(glm::vec3{0.0, +0.05, 0.0});
      },
      gl::Key::w);
  window.on_key_down(
      [&]() {
        camera.transform().translate(glm::vec3{-0.05, 0.0, 0.0});
      },
      gl::Key::a);
  window.on_key_down(
      [&]() {
        camera.transform().translate(glm::vec3{0.0, -0.05, 0.0});
      },
      gl::Key::s);
  window.on_key_down(
      [&]() {
        camera.transform().translate(glm::vec3{+0.05, 0.0, 0.0});
      },
      gl::Key::d);
  window.on_key_down(
      [&]() { //
        camera.transform().translate(glm::vec3{0.0, 0.0, +0.05});
      },
      gl::Key::q);
  window.on_key_down(
      [&]() { //
        camera.transform().translate(glm::vec3{0.0, 0.0, -0.05});
      },
      gl::Key::e);

  while (!glfwWindowShouldClose(window)) {
    glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    gl::BindProgram bind_program{program};
    const auto vp = camera.view_projection_matrix();
    glUniformMatrix4fv(glGetUniformLocation(program, "view_projection"), 1,
                       GL_FALSE, &vp[0][0]);
    glm::vec3 col{0.9f, 0.9f, 0.9f};
    glUniform3f(glGetUniformLocation(program, "col"), col.r, col.g, col.b);
    cell_renderer.draw(camera, program);
    col = glm::vec3{0.1f, 0.1f, 0.9f};
    glUniform3f(glGetUniformLocation(program, "col"), col.r, col.g, col.b);
    edge_renderer.draw(camera, program, GL_LINES);

    window.update();
  }
}

} // namespace Storm::Vulture
