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

#include <Storm/Vulture/OpenGL.hpp>
#include <Storm/Vulture/Window.hpp>

#include <concepts>
#include <cstddef> // std::byte
#include <optional>
#include <ranges>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#include <glm/glm.hpp>
// clang-format off
#include <GL/glew.h>
#include <GLFW/glfw3.h>
// clang-format on

namespace Storm {

template<mesh Mesh>
void visualize_mesh(const Mesh& mesh) {
  glfwInit();

  constexpr static size_t window_width = 1600;
  constexpr static size_t window_height = 900;
  constexpr static const char* window_title = "Mallard Mesh Visualizer";
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GL_TRUE);
  glfwWindowHint(GLFW_DOUBLEBUFFER, GL_TRUE);
  glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
  const auto window = glfwCreateWindow(window_width, window_height,
                                       window_title, nullptr, nullptr);
  if (window == nullptr) {
    STORM_ERROR_("Failed to create a GLFW window!");
    return;
  }
  glfwMakeContextCurrent(window);

  glewExperimental = GL_TRUE;
  if (GLenum status = glewInit(); status != GLEW_OK) {
    STORM_ERROR_("Failed to initialize GLEW: {}!", glewGetErrorString(status));
    return;
  }
  glViewport(0, 0, window_width, window_height);

  // Initialize.
  gl::DebugOutput debug_output{};
  gl::Shader vertex_shader{};
  vertex_shader.load(GL_VERTEX_SHADER, R"(
      #version 330 core
      layout(location = 0) in vec3 positionMS;
      uniform vec2 pos;
      uniform float scale;
      void main() {
        gl_Position = vec4(positionMS - vec3(1.0, 0.5, 0.0), 1.0);
        gl_Position.xy += pos;
        gl_Position.xyz *= scale;
        gl_Position.x /= 16.0 / 9.0;
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

  gl::Mesh edge_mesh;
  edge_mesh.push_vertex_array(mesh_nodes);
  edge_mesh.set_index_array(mesh_edges);
  edge_mesh.build();

  gl::Mesh cell_mesh;
  cell_mesh.push_vertex_array(mesh_nodes);
  cell_mesh.set_index_array(mesh_cells);
  cell_mesh.build();

  float scale = 1.01;
  glm::vec2 pos{};
  glm::vec3 col{0.1f, 0.1f, 0.9f};

  while (!glfwWindowShouldClose(window)) {
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) { pos.y += 0.05; }
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) { pos.x -= 0.05; }
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) { pos.y -= 0.05; }
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) { pos.x += 0.05; }
    if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS) { scale *= 1.01; }
    if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS) { scale /= 1.01; }

    // render frame
    {
      glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
      glClear(GL_COLOR_BUFFER_BIT);
      gl::BindProgram bind_program{program};
      glUniform2f(glGetUniformLocation(program, "pos"), pos.x, pos.y);
      glUniform1f(glGetUniformLocation(program, "scale"), scale);
      col = glm::vec3{0.9f, 0.9f, 0.9f};
      glUniform3f(glGetUniformLocation(program, "col"), col.r, col.g, col.b);
      cell_mesh.draw();
      col = glm::vec3{0.1f, 0.1f, 0.9f};
      glUniform3f(glGetUniformLocation(program, "col"), col.r, col.g, col.b);
      edge_mesh.draw(GL_LINES);
    }
    glfwSwapBuffers(window);
    glfwPollEvents();
  }

  glfwTerminate();
}

} // namespace Storm
