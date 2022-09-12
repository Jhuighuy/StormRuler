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

STORM_VULTURE_SHADER_(vertex, R"(
#version 330 core

layout(location = 0) in vec4 position_world_space;

uniform mat4 view_projection_matrix;

void main() {
  gl_Position = view_projection_matrix * position_world_space;
}
)")

STORM_VULTURE_SHADER_(fragment, R"(
#version 330 core

in vec2 position_model_space;

out vec4 fragment_color;

uniform usamplerBuffer edge_states;

const vec4 regular_color = vec4(0.8, 0.8, 0.8, 1.0);
const vec4 selected_color = vec4(0.8, 0.2, 0.8, 1.0);

void main() {
  uint state = texelFetch(edge_states, gl_PrimitiveID).r;
  fragment_color = state == uint(0) ? regular_color : selected_color;
}
)")
