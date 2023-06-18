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

#include <Storm/Utils/FilteringStreambuf.hpp>
#include <Storm/Utils/Meta.hpp>

#include <Storm/Mallard/Mesh.hpp>
#include <Storm/Mallard/Shape.hpp>

#include <filesystem>
#include <fstream>
#include <type_traits>
#include <vector>

namespace Storm {

// -----------------------------------------------------------------------------

/// @brief Read the @p mesh, generated by TetGen (Triangle in 2D).
/// @param path Path to the one of TetGen output file.
/// @see https://wias-berlin.de/software/tetgen/1.5/doc/manual/manual.pdf
/// @warning Unstable API!
template<mesh Mesh>
  requires (_detail::_in_range(mesh_dim_v<Mesh>, size_t{2}, size_t{3}))
void read_mesh_from_tetgen(Mesh& mesh, std::filesystem::path path) {
  constexpr bool mesh3D = mesh_dim_v<Mesh> == 3;

  STORM_INFO("Loading the TetGen/Triangle mesh from path '{}'...",
             path.replace_extension("*").string());

  { // Read the nodes.
    path.replace_extension(".node");
    std::ifstream node_stream(path);
    if (!node_stream.is_open()) {
      STORM_THROW_IO("Cannot open the node file '{}'!", path.string());
    }
    FilteringStreambuf<char, '#', '\n'> filter_comments(node_stream);

    // Read the nodes header.
    size_t num_nodes, dim, num_node_attribs;
    bool nodes_have_labels;
    node_stream >> num_nodes >> dim >> num_node_attribs >> nodes_have_labels;
    if (node_stream.bad()) {
      STORM_THROW_IO("Cannot read the node file '{}' header!", path.string());
    }
    if (dim != mesh_dim_v<Mesh>) {
      STORM_THROW_IO("Unexpected number of the dimensions in node file '{}' "
                     "header! Expected {}, got {}.",
                     path.string(), mesh_dim_v<Mesh>, dim);
    }
    if (num_node_attribs != 0) {
      STORM_WARNING("Header of the node file '{}' specifies {} attributes per "
                    "node, which are ignored!",
                    path.string(), num_node_attribs);
    }

    // Read the nodes.
    mesh.reserve(num_nodes, meta::type_v<NodeIndex>);
    for (size_t node_entry = 0; node_entry < num_nodes; ++node_entry) {
      size_t node_index, node_label = 0;
      mesh_vec_t<Mesh> node_pos{};
      std::vector<real_t> node_attribs(num_node_attribs);
      node_stream >> node_index >> node_pos(0) >> node_pos(1);
      if constexpr (mesh3D) node_stream >> node_pos(2);
      for (real_t& node_attrib : node_attribs) node_stream >> node_attrib;
      if (nodes_have_labels) node_stream >> node_label;
      if (node_stream.bad()) {
        STORM_THROW_IO("Cannot read the node # {} from file '{}'!", //
                       node_entry, path.string());
      }
      mesh.insert(node_pos, meta::type_v<NodeIndex>);
      /// @todo node_label, node_attribs!
    }

    STORM_INFO( //
        "Done reading {} nodes from file '{}'.", num_nodes, path.string());
  }

  std::vector<Label> edge_labels{};
  { // Read the edges.
    path.replace_extension(".edge");
    std::ifstream edge_stream(path);
    if (!edge_stream.is_open()) {
      STORM_THROW_IO("Cannot open the edge file '{}'!", path.string());
    }
    FilteringStreambuf<char, '#', '\n'> filter_comments(edge_stream);

    // Read the edges header.
    size_t num_edges;
    bool edges_have_labels;
    edge_stream >> num_edges >> edges_have_labels;
    if (edge_stream.bad()) {
      STORM_THROW_IO("Cannot read the edge file '{}' header!", path.string());
    }

    // Read the edges.
    mesh.reserve(num_edges, meta::type_v<EdgeIndex>);
    if (edges_have_labels) edge_labels.reserve(num_edges);
    for (size_t edge_entry = 0; edge_entry < num_edges; ++edge_entry) {
      size_t edge_index;
      Label edge_label;
      shapes::Seg edge{};
      edge_stream >> edge_index >> edge.n1 >> edge.n2;
      if (edges_have_labels) edge_stream >> edge_label;
      if (edge_stream.bad()) {
        STORM_THROW_IO("Cannot read the edge # {} from file '{}'!", //
                       edge_entry, path.string());
      }
      mesh.insert(edge, meta::type_v<EdgeIndex>);
      if (edges_have_labels) edge_labels.push_back(edge_label);
    }

    STORM_INFO( //
        "Done reading {} edges from file '{}'.", num_edges, path.string());
  }

  [[maybe_unused]] std::vector<Label> face_labels{};
  if constexpr (mesh3D) {
    // Read the faces.
    path.replace_extension(".face");
    std::ifstream face_stream(path);
    if (!face_stream.is_open()) {
      STORM_THROW_IO("Cannot open the face file '{}'!", path.string());
    }
    FilteringStreambuf<char, '#', '\n'> filter_comments(face_stream);

    // Read the faces header.
    size_t num_faces;
    bool faces_have_labels;
    face_stream >> num_faces >> faces_have_labels;
    if (face_stream.bad()) {
      STORM_THROW_IO("Cannot read the face file '{}' header!", path.string());
    }

    // Read the faces.
    mesh.reserve(num_faces, meta::type_v<FaceIndex<Mesh>>);
    if (faces_have_labels) face_labels.reserve(num_faces);
    for (size_t face_entry = 0; face_entry < num_faces; ++face_entry) {
      size_t face_index;
      Label face_label;
      shapes::Triangle face{};
      face_stream >> face_index >> face.n1 >> face.n2 >> face.n3;
      if (faces_have_labels) face_stream >> face_label;
      if (face_stream.bad()) {
        STORM_THROW_IO("Cannot read the face # {} from file '{}'!", //
                       face_entry, path.string());
      }
      mesh.insert(face, meta::type_v<FaceIndex<Mesh>>);
      if (faces_have_labels) face_labels.push_back(face_label);
    }

    STORM_INFO( //
        "Done reading {} faces from file '{}'.", num_faces, path.string());
  }

  { // Read the cells.
    path.replace_extension(".ele");
    std::ifstream cell_stream(path);
    if (!cell_stream.is_open()) {
      STORM_THROW_IO("Cannot open the cell file '{}'!", path.string());
    }
    FilteringStreambuf<char, '#', '\n'> filter_comments(cell_stream);

    // Read the cells header.
    size_t num_cells, num_cell_nodes;
    bool cells_have_attribs;
    cell_stream >> num_cells >> num_cell_nodes >> cells_have_attribs;
    if (cell_stream.bad()) {
      STORM_THROW_IO("Cannot read the cell file '{}' header!", path.string());
    }
    if (num_cell_nodes != mesh_dim_v<Mesh> + 1) {
      STORM_THROW_IO("Unexpected number of the nodes per cell in the cell "
                     "file '{}' header! Expected {}, got {}.",
                     path.string(), mesh_dim_v<Mesh> + 1, num_cell_nodes);
    }
    if (cells_have_attribs) {
      STORM_WARNING("Header of the cell file '{}' specifies a regional "
                    "attribute per cell, which is ignored!",
                    path.string());
    }

    // Read the cells.
    mesh.reserve(num_cells, meta::type_v<CellIndex<Mesh>>);
    for (size_t cell_entry = 0; cell_entry < num_cells; ++cell_entry) {
      size_t cell_index, cell_attrib = 0;
      std::conditional_t<mesh3D, shapes::Tetrahedron, shapes::Triangle> cell{};
      cell_stream >> cell_index >> cell.n1 >> cell.n2 >> cell.n3;
      if constexpr (mesh3D) cell_stream >> cell.n4;
      if (cells_have_attribs) cell_stream >> cell_attrib;
      if (cell_stream.bad()) {
        STORM_THROW_IO("Cannot read the cell # {} from file '{}'!", //
                       cell_entry, path.string());
      }
      mesh.insert(cell, meta::type_v<CellIndex<Mesh>>); /// @todo cell_attrib!
    }

    STORM_INFO( //
        "Done reading {} cells from file '{}'.", num_cells, path.string());
  }

  // Assign the labels.
  // (This should be done at the end, since
  //  TetGen may not generate all the edges/faces.)
  if (!edge_labels.empty()) {
    mesh.assign_labels(edge_labels, meta::type_v<EdgeIndex>);
    STORM_INFO("Done assigning edge labels.");
  }
  if constexpr (mesh3D) {
    if (!face_labels.empty()) {
      mesh.assign_labels(face_labels, meta::type_v<FaceIndex<Mesh>>);
      STORM_INFO("Done assigning face labels.");
    }
  }
}

// -----------------------------------------------------------------------------

} // namespace Storm
