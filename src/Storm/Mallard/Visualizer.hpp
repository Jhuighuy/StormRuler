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
  gl::Shader vertex_shader{gl::ShaderType::vertex};
  vertex_shader.load(R"(
      #version 330 core
      layout(location = 0) in vec3 positionMS;
      uniform mat4 view_projection;
      void main() {
        gl_Position = vec4(positionMS - vec3(1.0, 0.5, 0.0), 1.0);
        gl_Position = view_projection * gl_Position;
      }
    )");

  gl::Shader geometry_shader{gl::ShaderType::geometry};
  geometry_shader.load(R"(
      #version 330 core
      layout(points) in;
      layout(triangle_strip, max_vertices = 4) out;
      void main() {
        vec4 pos = gl_in[0].gl_Position;
        gl_Position = pos + vec4(-0.002, -0.002, 0.0, 0.0);
        EmitVertex();   
        gl_Position = pos + vec4( 0.002, -0.002, 0.0, 0.0);
        EmitVertex();
        gl_Position = pos + vec4(-0.002,  0.002, 0.0, 0.0);
        EmitVertex();
        gl_Position = pos + vec4( 0.002,  0.002, 0.0, 0.0);
        EmitVertex();
        EndPrimitive();
      }
    )");

  gl::Shader fragment_shader{gl::ShaderType::fragment};
  fragment_shader.load(R"(
      #version 330 core
      out vec4 color;
      uniform vec3 col;
      void main() {
        color = vec4(col, 1.0);
        color.x = gl_PrimitiveID/79672.0f;
      }
    )");

  gl::Program program;
  program.load(vertex_shader, fragment_shader);

  gl::Program node_program{};
  node_program.load(vertex_shader, geometry_shader, fragment_shader);

  const gl::Buffer nodes_buffer(
      nodes(mesh) | std::views::transform([](NodeView<const Mesh> node) {
        return static_cast<glm::vec2>(node.position());
      }));
  const gl::Buffer<GLuint> edge_nodes_buffer(
      edges(mesh) | std::views::transform([](EdgeView<const Mesh> edge) {
        return edge.nodes() |
               std::views::transform([](NodeView<const Mesh> node) {
                 return static_cast<GLuint>(node.index());
               });
      }) |
      std::views::join);
  const gl::Buffer<GLuint> cell_nodes_buffer(
      cells(mesh) | std::views::transform([](CellView<const Mesh> cell) {
        return cell.nodes() |
               std::views::transform([](NodeView<const Mesh> node) {
                 return static_cast<GLuint>(node.index());
               });
      }) |
      std::views::join);

  gl::VertexArray node_vertex_array{};
  node_vertex_array.build(nodes_buffer);

  gl::VertexArray edge_vertex_array{};
  edge_vertex_array.build_indexed(edge_nodes_buffer, nodes_buffer);

  gl::VertexArray cell_vertex_array{};
  cell_vertex_array.build_indexed(cell_nodes_buffer, nodes_buffer);

  // Setup camera.
  scene::Camera camera{};
  camera.transform().translate(glm::vec3{0.0f, 0.0f, 1.0f});
  camera.set_perspective(16.0 / 9.0);
  window.on_key_down(
      [&] {
        camera.transform().translate(glm::vec3{0.0, +0.05, 0.0});
      },
      gl::Key::w);
  window.on_key_down(
      [&] {
        camera.transform().translate(glm::vec3{-0.05, 0.0, 0.0});
      },
      gl::Key::a);
  window.on_key_down(
      [&] {
        camera.transform().translate(glm::vec3{0.0, -0.05, 0.0});
      },
      gl::Key::s);
  window.on_key_down(
      [&] {
        camera.transform().translate(glm::vec3{+0.05, 0.0, 0.0});
      },
      gl::Key::d);
  window.on_key_down(
      [&] {
        camera.transform().rotate_degrees(glm::vec3{0.0, +0.5, 0.0});
      },
      gl::Key::q);
  window.on_key_down(
      [&] {
        camera.transform().rotate_degrees(glm::vec3{0.0, -0.5, 0.0});
      },
      gl::Key::e);

  while (!glfwWindowShouldClose(window)) {
    glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    gl::BindProgram bind_program{program};

    bind_program.set_uniform(program.uniform_location("view_projection"),
                             camera.view_projection_matrix());
    bind_program.set_uniform(program.uniform_location("col"),
                             glm::vec3{0.9f, 0.9f, 0.9f});
    glBindVertexArray(cell_vertex_array.vertex_array_id_);
    glDrawElements(GL_TRIANGLES, num_cells(mesh) * 3, GL_UNSIGNED_INT,
                   /*indices*/ nullptr);
    // cell_renderer.draw(camera, program);

    bind_program.set_uniform(program.uniform_location("col"),
                             glm::vec3{0.1f, 0.1f, 0.9f});
    glBindVertexArray(edge_vertex_array.vertex_array_id_);
    glDrawElements(GL_LINES, num_edges(mesh) * 2, GL_UNSIGNED_INT,
                   /*indices*/ nullptr);
    //  edge_renderer.draw(camera, program, GL_LINES);

    {
      gl::BindProgram bind_program{node_program};
      bind_program.set_uniform(program.uniform_location("view_projection"),
                               camera.view_projection_matrix());
      bind_program.set_uniform(program.uniform_location("col"),
                               glm::vec3{0.9f, 0.1f, 0.1f});
      glBindVertexArray(node_vertex_array.vertex_array_id_);
      glDrawArrays(GL_POINTS, 0, num_nodes(mesh));
    }

    window.update();
  }
}

} // namespace Storm::Vulture
