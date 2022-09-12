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
#include <Storm/Vulture/GlFramebuffer.hpp>
#include <Storm/Vulture/GlShader.hpp>
#include <Storm/Vulture/GlTexture.hpp>
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

  // Setup window.
  constexpr static const char* window_title = "Strom::Vulture Visualizer";
  constexpr static size_t window_width = 1600;
  constexpr static size_t window_height = 900;
  gl::Window window{};
  window.load(framework, window_title, window_width, window_height);
  gl::BindWindow bind_window{window};
  glViewport(0, 0, window_width, window_height);
  gl::DebugOutput debug_output{};

  // Setup framebuffer.
  gl::Texture2D<glm::vec4> color_texture{window_width, window_height};
  gl::Framebuffer framebuffer{color_texture};
  gl::Buffer screen_quad_buffer{
      std::array{glm::vec2{+1.0f, +1.0f}, glm::vec2{+1.0f, -1.0f},
                 glm::vec2{-1.0f, +1.0f}, glm::vec2{+1.0f, -1.0f},
                 glm::vec2{-1.0f, -1.0f}, glm::vec2{-1.0f, +1.0f}}};
  gl::VertexArray screen_quad_vertex_array{screen_quad_buffer};
  gl::Program screen_quad_program{
#include <Storm/Vulture/ShaderScreenQuad.glsl>
  };

  // Setup nodes.
  const gl::Buffer node_positions_buffer(
      nodes(mesh) | std::views::transform([](NodeView<const Mesh> node) {
        return static_cast<glm::vec2>(node.position());
      }));
  gl::VertexArray node_vertex_array{node_positions_buffer};
  gl::Program node_program{
#include <Storm/Vulture/ShaderNodes.glsl>
  };
  // Setup node data.
  gl::Buffer<GLuint> node_states_buffer(num_nodes(mesh));
  const gl::TextureBuffer node_states_texture_buffer{node_states_buffer};

  // Setup edges.
  const gl::Buffer edge_nodes_buffer(
      edges(mesh) | std::views::transform([](EdgeView<const Mesh> edge) {
        return edge.nodes() |
               std::views::transform([](NodeView<const Mesh> node) {
                 return static_cast<GLuint>(node.index());
               });
      }) |
      std::views::join);
  gl::VertexArray edge_vertex_array{gl::indexed_v, //
                                    edge_nodes_buffer, node_positions_buffer};
  gl::Program edge_program{
#include <Storm/Vulture/ShaderEdges.glsl>
  };
  // Setup edge data.
  gl::Buffer<GLuint> edge_states_buffer(num_edges(mesh));
  const gl::TextureBuffer edge_states_texture_buffer{edge_states_buffer};

  // Setup cells.
  const gl::Buffer cell_nodes_buffer(
      cells(mesh) | std::views::transform([](CellView<const Mesh> cell) {
        return cell.nodes() |
               std::views::transform([](NodeView<const Mesh> node) {
                 return static_cast<GLuint>(node.index());
               });
      }) |
      std::views::join);
  gl::VertexArray cell_vertex_array{gl::indexed_v, //
                                    cell_nodes_buffer, node_positions_buffer};
  gl::Program cell_program{
#include <Storm/Vulture/ShaderCells.glsl>
  };
  // Setup cell data.
  const gl::Buffer cell_data_buffer(cell_data);
  const gl::TextureBuffer cell_data_texture_buffer{cell_data_buffer};
  gl::Buffer<GLuint> cell_states_buffer(num_cells(mesh));
  const gl::TextureBuffer cell_states_texture_buffer{cell_states_buffer};

  NodeView node(mesh, NodeIndex{13});
  node_states_buffer.set(static_cast<size_t>(node.index()), 1);
  node.for_each_edge([&](auto edge) { //
    edge_states_buffer.set(static_cast<size_t>(edge.index()), 1);
  });

  // Setup camera.
  scene::Camera camera{};
  camera.set_perspective(16.0 / 9.0);
  const auto reset_camera = [&]() {
    camera.transform() = scene::Transform{};
    const auto center = mesh.aabb().center();
    camera.transform().translate(
        glm::vec3(static_cast<glm::vec2>(center), 0.0f));
    const auto extents = mesh.aabb().extents();
    float orbit = static_cast<float>(0.5 * glm::max(extents.x, extents.y));
    camera.set_orbit(orbit);
  };
  reset_camera();
  // Reset the camera.
  window.on_key_up({gl::Key::r}, reset_camera);
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
  window.on_scroll([&](const glm::dvec2& offset) {
    camera.set_orbit(camera.orbit() - 0.25 * offset.y);
  });
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

  bool draw_nodes = false;
  window.on_key_up(gl::Key::n, [&] { draw_nodes = !draw_nodes; });
  bool draw_edges = false;
  window.on_key_up(gl::Key::m, [&] { draw_edges = !draw_edges; });
  bool draw_cells = true;
  window.on_key_up(gl::Key::b, [&] { draw_cells = !draw_cells; });

  window.main_loop([&] {
    framebuffer.draw_into([&]() {
      const auto view_projection_matrix = camera.view_projection_matrix();

      if (draw_cells) {
        glClearColor(0.2f, 0.2f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        gl::BindProgram bind_program{cell_program};
        cell_states_texture_buffer.bind(0);
        cell_data_texture_buffer.bind(1);
        bind_program.set_uniform(cell_program["cell_states"], 0);
        bind_program.set_uniform(cell_program["cell_values"], 1);
        bind_program.set_uniform(cell_program["view_projection_matrix"],
                                 view_projection_matrix);
        cell_vertex_array.draw_indexed(gl::DrawMode::triangles,
                                       num_cells(mesh) * 3);
      }

      if (draw_edges) {
        gl::BindProgram bind_program{edge_program};
        edge_states_texture_buffer.bind(0);
        bind_program.set_uniform(edge_program["edge_states"], 0);
        bind_program.set_uniform(edge_program["view_projection_matrix"],
                                 view_projection_matrix);
        edge_vertex_array.draw_indexed(gl::DrawMode::lines,
                                       num_edges(mesh) * 2);
      }

      if (draw_nodes) {
        gl::BindProgram bind_program{node_program};
        node_states_texture_buffer.bind(0);
        bind_program.set_uniform(node_program["node_states"], 0);
        bind_program.set_uniform(node_program["view_projection_matrix"],
                                 view_projection_matrix);
        bind_program.set_uniform(node_program["point_size"],
                                 glm::vec2{0.002f * 9.0f / 16.0f, 0.002f});
        node_vertex_array.draw(gl::DrawMode::points, num_nodes(mesh));
      }
    });

    {
      gl::BindProgram bind_program{screen_quad_program};
      color_texture.bind(0);
      bind_program.set_uniform(node_program["frame_texture"], 0);
      screen_quad_vertex_array.draw(gl::DrawMode::triangles, 6);
    }
  });
}

} // namespace Storm::Vulture
