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

#include <concepts>
#include <functional>

#include <glm/glm.hpp>
// clang-format off
#include <GL/glew.h>
#include <GLFW/glfw3.h>
// clang-format on

namespace Storm::gl {

/// @brief Input modifiers.
enum class Modifiers : int {
  shift = GLFW_MOD_SHIFT,
  control = GLFW_MOD_CONTROL,
  alt = GLFW_MOD_ALT,
  super = GLFW_MOD_SUPER,
  caps_lock = GLFW_MOD_CAPS_LOCK,
  num_lock = GLFW_MOD_NUM_LOCK,
}; // enum class Modifiers

[[nodiscard]] constexpr inline int operator&(Modifiers a,
                                             Modifiers b) noexcept {
  return static_cast<int>(a) & static_cast<int>(b);
}

/// @brief Keyboard key.
enum class Key : int {
  space = GLFW_KEY_SPACE,
  apostrophe = GLFW_KEY_APOSTROPHE,
  comma = GLFW_KEY_COMMA,
  minus = GLFW_KEY_MINUS,
  period = GLFW_KEY_PERIOD,
  slash = GLFW_KEY_SLASH,
  _0 = GLFW_KEY_0,
  _1 = GLFW_KEY_1,
  _2 = GLFW_KEY_2,
  _3 = GLFW_KEY_3,
  _4 = GLFW_KEY_4,
  _5 = GLFW_KEY_5,
  _6 = GLFW_KEY_6,
  _7 = GLFW_KEY_7,
  _8 = GLFW_KEY_8,
  _9 = GLFW_KEY_9,
  semicolon = GLFW_KEY_SEMICOLON,
  equal = GLFW_KEY_EQUAL,
  a = GLFW_KEY_A,
  b = GLFW_KEY_B,
  c = GLFW_KEY_C,
  d = GLFW_KEY_D,
  e = GLFW_KEY_E,
  f = GLFW_KEY_F,
  g = GLFW_KEY_G,
  h = GLFW_KEY_H,
  i = GLFW_KEY_I,
  j = GLFW_KEY_J,
  k = GLFW_KEY_K,
  l = GLFW_KEY_L,
  m = GLFW_KEY_M,
  n = GLFW_KEY_N,
  o = GLFW_KEY_O,
  p = GLFW_KEY_P,
  q = GLFW_KEY_Q,
  r = GLFW_KEY_R,
  s = GLFW_KEY_S,
  t = GLFW_KEY_T,
  u = GLFW_KEY_U,
  v = GLFW_KEY_V,
  w = GLFW_KEY_W,
  x = GLFW_KEY_X,
  y = GLFW_KEY_Y,
  z = GLFW_KEY_Z,
  left_bracket = GLFW_KEY_LEFT_BRACKET,
  backslash = GLFW_KEY_BACKSLASH,
  right_bracket = GLFW_KEY_RIGHT_BRACKET,
  grave_accent = GLFW_KEY_GRAVE_ACCENT,
  world_1 = GLFW_KEY_WORLD_1,
  world_2 = GLFW_KEY_WORLD_2,
  escape = GLFW_KEY_ESCAPE,
  enter = GLFW_KEY_ENTER,
  tab = GLFW_KEY_TAB,
  backspace = GLFW_KEY_BACKSPACE,
  insert = GLFW_KEY_INSERT,
  del = GLFW_KEY_DELETE,
  right = GLFW_KEY_RIGHT,
  left = GLFW_KEY_LEFT,
  down = GLFW_KEY_DOWN,
  up = GLFW_KEY_UP,
  page_up = GLFW_KEY_PAGE_UP,
  page_down = GLFW_KEY_PAGE_DOWN,
  home = GLFW_KEY_HOME,
  end = GLFW_KEY_END,
  caps_lock = GLFW_KEY_CAPS_LOCK,
  scroll_lock = GLFW_KEY_SCROLL_LOCK,
  num_lock = GLFW_KEY_NUM_LOCK,
  print_screen = GLFW_KEY_PRINT_SCREEN,
  pause = GLFW_KEY_PAUSE,
  f1 = GLFW_KEY_F1,
  f2 = GLFW_KEY_F2,
  f3 = GLFW_KEY_F3,
  f4 = GLFW_KEY_F4,
  f5 = GLFW_KEY_F5,
  f6 = GLFW_KEY_F6,
  f7 = GLFW_KEY_F7,
  f8 = GLFW_KEY_F8,
  f9 = GLFW_KEY_F9,
  f10 = GLFW_KEY_F10,
  f11 = GLFW_KEY_F11,
  f12 = GLFW_KEY_F12,
  f13 = GLFW_KEY_F13,
  f14 = GLFW_KEY_F14,
  f15 = GLFW_KEY_F15,
  f16 = GLFW_KEY_F16,
  f17 = GLFW_KEY_F17,
  f18 = GLFW_KEY_F18,
  f19 = GLFW_KEY_F19,
  f20 = GLFW_KEY_F20,
  f21 = GLFW_KEY_F21,
  f22 = GLFW_KEY_F22,
  f23 = GLFW_KEY_F23,
  f24 = GLFW_KEY_F24,
  f25 = GLFW_KEY_F25,
  kp_0 = GLFW_KEY_KP_0,
  kp_1 = GLFW_KEY_KP_1,
  kp_2 = GLFW_KEY_KP_2,
  kp_3 = GLFW_KEY_KP_3,
  kp_4 = GLFW_KEY_KP_4,
  kp_5 = GLFW_KEY_KP_5,
  kp_6 = GLFW_KEY_KP_6,
  kp_7 = GLFW_KEY_KP_7,
  kp_8 = GLFW_KEY_KP_8,
  kp_9 = GLFW_KEY_KP_9,
  kp_decimal = GLFW_KEY_KP_DECIMAL,
  kp_divide = GLFW_KEY_KP_DIVIDE,
  kp_multiply = GLFW_KEY_KP_MULTIPLY,
  kp_subtract = GLFW_KEY_KP_SUBTRACT,
  kp_add = GLFW_KEY_KP_ADD,
  kp_enter = GLFW_KEY_KP_ENTER,
  kp_equal = GLFW_KEY_KP_EQUAL,
  left_shift = GLFW_KEY_LEFT_SHIFT,
  left_control = GLFW_KEY_LEFT_CONTROL,
  left_alt = GLFW_KEY_LEFT_ALT,
  left_super = GLFW_KEY_LEFT_SUPER,
  right_shift = GLFW_KEY_RIGHT_SHIFT,
  right_control = GLFW_KEY_RIGHT_CONTROL,
  right_alt = GLFW_KEY_RIGHT_ALT,
  right_super = GLFW_KEY_RIGHT_SUPER,
  menu = GLFW_KEY_MENU,
}; // enum class Key

/// @brief Mouse button.
enum class MouseButton : int {
  _1 = GLFW_MOUSE_BUTTON_1,
  _2 = GLFW_MOUSE_BUTTON_2,
  _3 = GLFW_MOUSE_BUTTON_3,
  _4 = GLFW_MOUSE_BUTTON_4,
  _5 = GLFW_MOUSE_BUTTON_5,
  _6 = GLFW_MOUSE_BUTTON_6,
  _7 = GLFW_MOUSE_BUTTON_7,
  _8 = GLFW_MOUSE_BUTTON_8,
  left = GLFW_MOUSE_BUTTON_LEFT,
  right = GLFW_MOUSE_BUTTON_RIGHT,
  middle = GLFW_MOUSE_BUTTON_MIDDLE,
}; // enum class MouseButton

/// @brief OpenGL window.
class Window {
private:

  GLFWwindow* window_;
  std::function<void()> on_close_fn_{};
  std::function<void(size_t, size_t)> on_resize_fn_{};
  std::function<void(Key, Modifiers)> on_key_down_fn_{};
  std::function<void(Key, Modifiers)> on_key_up_fn_{};
  std::function<void(MouseButton, Modifiers)> on_mouse_button_down_fn_{};
  std::function<void(MouseButton, Modifiers)> on_mouse_button_up_fn_{};
  std::function<void(glm::dvec2)> on_scroll_fn_{};
  std::function<void(glm::dvec2)> on_set_cursor_pos_fn_{};

public:

  void init() {
    // Set the callbacks.
    glfwSetWindowUserPointer(window_, this);
    glfwSetWindowCloseCallback(window_, &on_close_);
    glfwSetWindowSizeCallback(window_, &on_resize_);
    glfwSetKeyCallback(window_, &on_key_);
    glfwSetMouseButtonCallback(window_, &on_mouse_button_);
    glfwSetScrollCallback(window_, &on_scroll_);
    glfwSetCursorPosCallback(window_, &on_set_cursor_pos_);
  }

public:

  /// @brief Set handler for the window close event.
  template<std::invocable CloseFn>
  void on_close(CloseFn on_close_fn) {
    on_close_fn_ = std::move(on_close_fn);
  }

  /// @brief Set handler for the window resize event.
  template<std::invocable<size_t, size_t> ResizeFn>
  void on_close(ResizeFn on_resize_fn) {
    on_resize_fn_ = std::move(on_resize_fn);
  }

  /// @brief Set handler for the key pressed event.
  template<std::invocable<int, int> KeyDownFn>
  void on_key_down(KeyDownFn on_key_down_fn) {
    on_key_down_fn_ = std::move(on_key_down_fn);
  }

  /// @brief Set handler for the key released event.
  template<std::invocable<int, int> KeyUpFn>
  void on_key_up(KeyUpFn on_key_up_fn) {
    on_key_up_fn_ = std::move(on_key_up_fn);
  }

  /// @brief Set handler for the mouse button pressed event.
  template<std::invocable<int, int> MouseButtonDownFn>
  void on_mouse_button_down(MouseButtonDownFn on_mouse_button_down_fn) {
    on_mouse_button_down_fn_ = std::move(on_mouse_button_down_fn);
  }

  /// @brief Set handler for the mouse button released event.
  template<std::invocable<int, int> KeyUpFn>
  void on_mouse_button_up(KeyUpFn on_mouse_button_up_fn) {
    on_mouse_button_up_fn_ = std::move(on_mouse_button_up_fn);
  }

  /// @brief Set handler for the scroll event.
  template<std::invocable<glm::dvec2> ScrollFn>
  void on_scroll(ScrollFn on_scroll_fn) {
    on_scroll_fn_ = std::move(on_scroll_fn);
  }

  /// @brief Set handler for the set cursor position event.
  template<std::invocable<glm::dvec2> SetCursorPosFn>
  void on_set_cursor_pos(SetCursorPosFn on_set_cursor_pos_fn) {
    on_set_cursor_pos_fn_ = std::move(on_set_cursor_pos_fn);
  }

private:

  static const Window* get_self_(GLFWwindow* window) noexcept {
    const auto self =
        static_cast<const Window*>(glfwGetWindowUserPointer(window));
    STORM_ASSERT_(self != nullptr, "Pointer to self was not set!");
    return self;
  }

  static void on_close_(GLFWwindow* window) noexcept {
    const auto self = get_self_(window);
    if (self->on_close_fn_ != nullptr) { self->on_close_fn_(); }
  }

  static void on_resize_(GLFWwindow* window, int width, int height) noexcept {
    const auto self = get_self_(window);
    if (self->on_resize_fn_ != nullptr) {
      self->on_resize_fn_(static_cast<size_t>(width),
                          static_cast<size_t>(height));
    }
  }

  static void on_key_(GLFWwindow* window, int key,
                      [[maybe_unused]] int scancode, int action,
                      int mods) noexcept {
    const auto self = get_self_(window);
    switch (action) {
      case GLFW_PRESS:
        if (self->on_key_down_fn_ != nullptr) {
          self->on_key_down_fn_(static_cast<Key>(key),
                                static_cast<Modifiers>(mods));
        }
        break;
      case GLFW_RELEASE:
        if (self->on_key_up_fn_ != nullptr) { //
          self->on_key_up_fn_(static_cast<Key>(key),
                              static_cast<Modifiers>(mods));
        }
        break;
      default: //
        STORM_WARNING_("GLFW: unsupported key action {:#x}!", action);
        break;
    }
  }

  static void on_mouse_button_(GLFWwindow* window, int mouse_button, int action,
                               int mods) noexcept {
    const auto self = get_self_(window);
    switch (action) {
      case GLFW_PRESS:
        if (self->on_mouse_button_down_fn_ != nullptr) {
          self->on_mouse_button_down_fn_(static_cast<MouseButton>(mouse_button),
                                         static_cast<Modifiers>(mods));
        }
        break;
      case GLFW_RELEASE:
        if (self->on_mouse_button_up_fn_ != nullptr) {
          self->on_mouse_button_up_fn_(static_cast<MouseButton>(mouse_button),
                                       static_cast<Modifiers>(mods));
        }
        break;
      default: //
        STORM_WARNING_("GLFW: unsupported mouse button action {:#x}!", action);
        break;
    }
  }

  static void on_scroll_(GLFWwindow* window, //
                         double x_offset, double y_offset) noexcept {
    const auto self = get_self_(window);
    if (self->on_scroll_fn_ != nullptr) {
      self->on_scroll_fn_({x_offset, y_offset});
    }
  }

  static void on_set_cursor_pos_(GLFWwindow* window, //
                                 double x_pos, double y_pos) noexcept {
    const auto self = get_self_(window);
    if (self->on_set_cursor_pos_fn_ != nullptr) {
      self->on_set_cursor_pos_fn_({x_pos, y_pos});
    }
  }

}; // class Window

} // namespace Storm::gl
