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
#include <Storm/Vulture/GlVertexArray.hpp>
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
  IndexedVector<GLfloat, CellIndex<Mesh>> cell_data{};
  std::ranges::copy(
      cells(mesh) | std::views::transform([](CellView<const Mesh> cell) {
        const real_t v = std::sin(3.0 * cell.barycenter_position().x) *
                         std::cos(7.0 * cell.barycenter_position().y);
        return static_cast<GLfloat>((v + 1.0) / 2.0);
      }),
      std::back_inserter(cell_data));

  gl::Framework framework{};

  constexpr static const char* window_title = "Strom::Vulture Visualizer";
  constexpr static size_t window_width = 1600;
  constexpr static size_t window_height = 900;
  gl::Window window{};
  window.load(framework, window_title, window_width, window_height);
  gl::BindWindow bind_window{window};
  glViewport(0, 0, window_width, window_height);
  gl::DebugOutput debug_output{};

  // Initialize.
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

  gl::Shader node_geometry_shader{gl::ShaderType::geometry};
  node_geometry_shader.load(R"(
      #version 330 core
      layout(points) in;
      layout(triangle_strip, max_vertices = 4) out;
      out vec2 position;
      void main() {
        
        vec4 pos = gl_in[0].gl_Position;
        const vec2 sz = vec2(0.002*9.0/16.0, 0.002);

        position = vec2(-1, -1);
        gl_Position = pos + vec4(sz * position, 0.0, 0.0);
        EmitVertex();   

        position = vec2(+1, -1);
        gl_Position = pos + vec4(sz * position, 0.0, 0.0);
        EmitVertex();

        position = vec2(-1, +1);
        gl_Position = pos + vec4(sz * position, 0.0, 0.0);
        EmitVertex();

        position = vec2(+1, +1);
        gl_Position = pos + vec4(sz * position, 0.0, 0.0);
        EmitVertex();

        EndPrimitive();
      }
    )");

  gl::Shader node_fragment_shader{gl::ShaderType::fragment};
  node_fragment_shader.load(R"(
      #version 330 core
      in vec2 position;
      out vec4 color;
      void main() {
        if (length(position) > 1.0f) discard;
        color = vec4(0.9f, 0.9f, 0.9f, 1.0f);
      }
    )");

  gl::Shader edge_fragment_shader{gl::ShaderType::fragment};
  edge_fragment_shader.load(R"(
      #version 330 core
      out vec4 color;
      void main() {
        color = vec4(0.8f, 0.8f, 0.8f, 1.0f);
      }
    )");

  gl::Shader cell_fragment_shader{gl::ShaderType::fragment};
  cell_fragment_shader.load(R"(
      #version 330 core
      uniform samplerBuffer values;
      out vec4 color;
      void main() {
        float value = texelFetch(values, gl_PrimitiveID).r;
        value = clamp(value, 0.0f, 1.0f);
        color.rgb = mix(vec3(0.0f, 0.0f, 1.0f), vec3(1.0f, 0.0f, 0.0f), value);
        color.a = 1.0f;
      }
    )");

  gl::Program node_program{};
  node_program.load(vertex_shader, node_geometry_shader, node_fragment_shader);

  gl::Program edge_program{};
  edge_program.load(vertex_shader, edge_fragment_shader);

  gl::Program cell_program{};
  cell_program.load(vertex_shader, cell_fragment_shader);

  const gl::Buffer nodes_buffer(
      nodes(mesh) | std::views::transform([](NodeView<const Mesh> node) {
        return static_cast<glm::vec2>(node.position());
      }));
  const gl::Buffer edge_nodes_buffer(
      edges(mesh) | std::views::transform([](EdgeView<const Mesh> edge) {
        return edge.nodes() |
               std::views::transform([](NodeView<const Mesh> node) {
                 return static_cast<GLuint>(node.index());
               });
      }) |
      std::views::join);
  const gl::Buffer cell_nodes_buffer(
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

  // Setup data.
  const gl::Buffer cell_data_buffer(cell_data);
  GLuint texture_id_;
  glGenTextures(1, &texture_id_);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_BUFFER, texture_id_);
  glTexBuffer(GL_TEXTURE_BUFFER, GL_R32F, cell_data_buffer);

  // Setup camera.
  scene::Camera camera{};
  camera.set_perspective(16.0 / 9.0);
  // Reset the camera.
  window.on_key_up({gl::Key::r},
                   [&] { camera.transform() = scene::Transform{}; });
  // Translate the camera.
  window.on_key_down(gl::Key::w, [&] {
    camera.transform().translate(glm::vec3{0.0, +0.05, 0.0});
  });
  window.on_key_down(gl::Key::a, [&] {
    camera.transform().translate(glm::vec3{-0.05, 0.0, 0.0});
  });
  window.on_key_down(gl::Key::s, [&] {
    camera.transform().translate(glm::vec3{0.0, -0.05, 0.0});
  });
  window.on_key_down(gl::Key::d, [&] {
    camera.transform().translate(glm::vec3{+0.05, 0.0, 0.0});
  });
  // Rotate the camera.
  window.on_key_up(gl::Key::q, [&] {
    camera.transform().rotate_degrees(glm::vec3{0.0f, +90.0f, 0.0f});
  });
  window.on_key_up(gl::Key::e, [&] {
    camera.transform().rotate_degrees(glm::vec3{0.0f, -90.0f, 0.0f});
  });
  window.on_key_down(gl::Key::up, [&] {
    camera.transform().rotate_degrees(glm::vec3{+0.5, 0.0, 0.0});
  });
  window.on_key_down(gl::Key::down, [&] {
    camera.transform().rotate_degrees(glm::vec3{-0.5, 0.0, 0.0});
  });
  window.on_key_down(gl::Key::left, [&] {
    camera.transform().rotate_degrees(glm::vec3{0.0, +0.5, 0.0});
  });
  window.on_key_down(gl::Key::right, [&] {
    camera.transform().rotate_degrees(glm::vec3{0.0, -0.5, 0.0});
  });
  window.on_scroll([&](const glm::dvec2& offset) {
    camera.set_orbit(camera.orbit() - 0.25 * offset.y);
  });

  bool draw_nodes = false;
  window.on_key_up(gl::Key::n, [&] { draw_nodes = !draw_nodes; });
  bool draw_edges = false;
  window.on_key_up(gl::Key::m, [&] { draw_edges = !draw_edges; });

  while (!glfwWindowShouldClose(window)) {
    glClearColor(0.2f, 0.2f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    const auto view_projection = camera.view_projection_matrix();

    {
      gl::BindProgram bind_program{cell_program};
      bind_program.set_uniform(cell_program.uniform_location("values"), 0);
      bind_program.set_uniform(cell_program.uniform_location("view_projection"),
                               view_projection);
      cell_vertex_array.draw_indexed(gl::DrawMode::triangles,
                                     num_cells(mesh) * 3);
    }

    if (draw_edges) {
      gl::BindProgram bind_program{edge_program};
      bind_program.set_uniform(edge_program.uniform_location("view_projection"),
                               view_projection);
      edge_vertex_array.draw_indexed(gl::DrawMode::lines, num_edges(mesh) * 2);
    }

    if (draw_nodes) {
      gl::BindProgram bind_program{node_program};
      bind_program.set_uniform(node_program.uniform_location("view_projection"),
                               view_projection);
      node_vertex_array.draw(gl::DrawMode::points, num_nodes(mesh));
    }

#if 0
    bind_program.set_uniform(program.uniform_location("col"),
                             glm::vec3{0.1f, 0.1f, 0.9f});
    glBindVertexArray(edge_vertex_array.vertex_array_id_);
    glDrawElements(GL_LINES, num_edges(mesh) * 2, GL_UNSIGNED_INT,
                   /*indices*/ nullptr);

    {
      gl::BindProgram bind_program{node_program};
      bind_program.set_uniform(program.uniform_location("view_projection"),
                               camera.view_projection_matrix());
      bind_program.set_uniform(program.uniform_location("col"),
                               glm::vec3{0.9f, 0.1f, 0.1f});
      glBindVertexArray(node_vertex_array.vertex_array_id_);
      glDrawArrays(GL_POINTS, 0, num_nodes(mesh));
    }
#endif

    window.update();
  }
}

} // namespace Storm::Vulture
