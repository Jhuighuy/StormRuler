/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// Copyright (C) 2022 Oleg Butakov
///
/// Permission is hereby granted, free of charge, to any person
/// obtaining a copy of this software and associated documentation
/// files (the "Software"), to deal in the Software without
/// restriction, including without limitation the rights  to use,
/// copy, modify, merge, publish, distribute, sublicense, and/or
/// sell copies of the Software, and to permit persons to whom the
/// Software is furnished to do so, subject to the following
/// conditions:
///
/// The above copyright notice and this permission notice shall be
/// included in all copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
/// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
/// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
/// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
/// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
/// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
/// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
/// OTHER DEALINGS IN THE SOFTWARE.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///

#include <GL/freeglut.h>
#include <cmath>
#include <utility>
#include <vector>

#include <Storm/Base.hpp>
#include <Storm/Blass/MatrixDense.hpp>

/// @todo Change `g_view_rotation_{x|y}` to `g_view_rotation.{x|y}`

namespace Storm::SpherePack::GL {

using Vec3D = Storm::DenseVector<real_t, 3>;

struct Sphere {
  Vec3D x_curr, x_prev;
};

std::vector<Sphere> g_spheres;

static real_t g_cylinder_radius{1.25};
static real_t g_cylinder_height{1.50};

static real_t g_view_distance = 3.0;
static real_t g_view_rotation_x{-30.0}, g_view_rotation_y{0.0};

static void on_idle() {
  static int kkk = 0;
  if (kkk % 5 == 0 && g_spheres.size() < 200) {
    g_spheres.emplace_back();
    const real_t r = (g_cylinder_radius * rand()) / RAND_MAX;
    const real_t phi = (2.0 * M_PI * rand()) / RAND_MAX;
    g_spheres.back().x_curr(0) = r * std::cos(phi);
    g_spheres.back().x_curr(2) = r * std::sin(phi);
    g_spheres.back().x_curr(1) = g_cylinder_height * 0.5;
    g_spheres.back().x_prev = g_spheres.back().x_curr;
  }
  kkk += 1;
  for (size_t ss{0}; ss < 8; ++ss) {
    // Apply the cylinder constrains.
    for (Sphere& sphere : g_spheres) {
      Vec3D& x = sphere.x_curr;
      // Clamp on floor and ceil.
      x(1) = std::clamp(x(1), //
                        -0.5 * g_cylinder_height + 0.1,
                        +0.5 * g_cylinder_height - 0.1);
      // Clamp on walls.
      auto w = select_rows(x, 0, 2);
      if (real_t d = norm_2(w); d > g_cylinder_radius - 0.1) {
        w *= (g_cylinder_radius - 0.1) / d;
      }
    }

    // Apply the collision, O(n²).
    for (size_t i{0}; i < g_spheres.size(); ++i) {
      Vec3D& x_i{g_spheres[i].x_curr};
      for (size_t j{i + 1}; j < g_spheres.size(); ++j) {
        Vec3D& x_j{g_spheres[j].x_curr};
        Vec3D axis;
        axis <<= x_i - x_j;
        if (const real_t dist = norm_2(axis); dist < 0.2) {
          Vec3D n;
          n <<= axis / dist;
          real_t delta{0.2 - dist};
          x_i += 0.5 * delta * n;
          x_j -= 0.5 * delta * n;
        }
      }
    }

    real_t en{};
    Vec3D down{};
    down(1) = -1.0;
    for (Sphere& sphere : g_spheres) {
      Vec3D new_r_previous{sphere.x_curr};
      sphere.x_curr <<= 2.0 * sphere.x_curr - sphere.x_prev + 0.00001 * down;
      sphere.x_prev <<= new_r_previous;
      en += norm_2(sphere.x_curr - sphere.x_prev);
    }
    en /= g_spheres.size();
    printf("en = %f\n", en);
  }

  // Require the redisplay.
  glutPostRedisplay();
}

static void on_display() {
  // Clear the screen with some pretty blue color.
  glClearColor(0.1f, 0.1f, 0.6f, 0.0f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Setup the projection matrix.
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(60.0, 1024.0 / 768.0, 0.1, 1000.0);

  // Setup the view matrix.
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0.0, 0.0, -g_view_distance, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
  glRotated(+g_view_rotation_x, 1.0, 0.0, 0.0);
  glRotated(-g_view_rotation_y, 0.0, 1.0, 0.0);

  // Draw the cylinder.
  glPushMatrix();
  glTranslated(0.0, 0.5 * g_cylinder_height, 0.0);
  glRotated(90.0, 1.0, 0.0, 0.0);
  glutWireCylinder(g_cylinder_radius, g_cylinder_height, 60, 1);
  glPopMatrix();

  // Draw the spheres.
  for (const Sphere& sphere : g_spheres) {
    glPushMatrix();
    glTranslated(sphere.x_curr(0), sphere.x_curr(1), sphere.x_curr(2));
    glutSolidSphere(0.1, 10, 5);
    glPopMatrix();
  }

  // Post the scene.
  glFlush();
  glutSwapBuffers();
}

static void on_key_pressed(unsigned char key, //
                           [[maybe_unused]] int x, [[maybe_unused]] int y) {
  constexpr real_t view_rotation_speed{0.5};
  switch (key) {
    case 'w': g_view_rotation_x += view_rotation_speed; break;
    case 'a': g_view_rotation_y -= view_rotation_speed; break;
    case 's': g_view_rotation_x -= view_rotation_speed; break;
    case 'd': g_view_rotation_y += view_rotation_speed; break;
  }
}

static int main(int argc, char** argv) {
  srand(123);
  glutInit(&argc, argv);

  glutInitContextVersion(3, 3);
  glutInitContextProfile(GLUT_COMPATIBILITY_PROFILE);

  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(1024, 768);
  glutInitWindowPosition(100, 100);
  glutCreateWindow("SpherePack");

  glutIdleFunc(on_idle);
  glutDisplayFunc(on_display);
  glutKeyboardFunc(on_key_pressed);

  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
  glShadeModel(GL_SMOOTH);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glViewport(0, 0, 1024, 768);

  glutMainLoop();

  return 0;
}

} // namespace Storm::SpherePack::GL

int main(int argc, char** argv) {
  return Storm::SpherePack::GL::main(argc, argv);
}
