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

#include <algorithm>
#include <concepts>
#include <functional>
#include <map>
#include <utility>
#include <vector>

#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

namespace Storm::Vulture::gl {

/// @brief Input modifiers.
enum class Modifiers : int {
  none = 0,
  shift = GLFW_MOD_SHIFT,
  control = GLFW_MOD_CONTROL,
  alt = GLFW_MOD_ALT,
  super = GLFW_MOD_SUPER,
  caps_lock = GLFW_MOD_CAPS_LOCK,
  num_lock = GLFW_MOD_NUM_LOCK,
}; // enum class Modifiers

[[nodiscard]] constexpr inline auto operator&(Modifiers a,
                                              Modifiers b) noexcept {
  return static_cast<int>(a) & static_cast<int>(b);
}
[[nodiscard]] constexpr inline auto operator|(Modifiers a,
                                              Modifiers b) noexcept {
  return static_cast<Modifiers>(static_cast<int>(a) | static_cast<int>(b));
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

/// @brief OpenGL framework.
class Framework final {
private:

  friend class Window;
  bool loaded_ = false;
  bool glew_loaded_ = false;

public:

  /// @brief Construct a framework.
  Framework() = default;

  /// @brief Move-construct a framework.
  Framework(Framework&&) = default;
  /// @brief Move-assign the framework.
  Framework& operator=(Framework&& window) = default;

  Framework(const Framework&) = delete;
  Framework& operator=(const Framework& window) = delete;

  /// @brief Destroy the framefowork.
  ~Framework() {
    unload();
  }

  /// @brief Load the framework.
  void load() {
    if (loaded_) { return; }
    glfwSetErrorCallback(&on_error_);
    glfwInit();
    STORM_INFO_("GLFW initialized, version '{}'!", glfwGetVersionString());
    loaded_ = true;
  }

  /// @brief Unload the framework.
  void unload() {
    if (!loaded_) { return; }
  }

private:

  static void on_error_(int error, const char* description) {
    STORM_THROW_GL_("GLFW: Error {:#x}: {}", error, description);
  }

}; // class Framework

/// @brief OpenGL window.
class Window final {
private:

  GLFWwindow* underlying_{};
  std::vector<std::function<void()>> on_close_fns_{};
  std::vector<std::function<void(size_t, size_t)>> on_resize_fns_{};
  std::multimap<std::pair<Key, Modifiers>, std::function<void()>>
      on_key_down_fns_{}, on_key_up_fns_{};
  std::multimap<std::pair<MouseButton, Modifiers>, std::function<void()>>
      on_mouse_button_down_fns_{}, on_mouse_button_up_fns_{};
  std::vector<std::function<void(glm::dvec2)>> on_scroll_fns_{};
  std::vector<std::function<void(glm::dvec2)>> on_set_cursor_pos_fns_{};

public:

  /// @brief Construct a window.
  Window() = default;

  /// @brief Move-construct a window.
  Window(Window&&) = default;
  /// @brief Move-assign the window.
  Window& operator=(Window&&) = default;

  Window(const Window&) = delete;
  Window& operator=(const Window& window) = delete;

  /// @brief Destroy the window.
  ~Window() {
    glfwDestroyWindow(underlying_);
  }

  /// @brief Underlying window pointer.
  [[nodiscard]] constexpr GLFWwindow* underlying() const noexcept {
    return underlying_;
  }
  /// @brief Cast to underlying window pointer.
  [[nodiscard]] constexpr operator GLFWwindow*() const noexcept {
    return underlying_;
  }

  /// @brief Get window width.
  [[nodiscard]] size_t width() const {
    STORM_ASSERT_(underlying_ != nullptr, "Window is not loaded!");
    int width;
    glfwGetWindowSize(underlying_, &width, nullptr);
    return static_cast<size_t>(width);
  }

  /// @brief Get window height.
  [[nodiscard]] size_t height() const {
    STORM_ASSERT_(underlying_ != nullptr, "Window is not loaded!");
    int height;
    glfwGetWindowSize(underlying_, nullptr, &height);
    return static_cast<size_t>(height);
  }

  /// @brief Load the window.
  void load(Framework& framework, //
            const char* title, size_t width, size_t height) {
    framework.load();

    // Set the window hints.
    glfwDefaultWindowHints();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GL_TRUE);
    glfwWindowHint(GLFW_SAMPLES, 4);

    // Create the window.
    STORM_ASSERT_(title != nullptr, "Invalid window title!");
    STORM_ASSERT_(width > 0 && height > 0, "Invalid window size!");
    underlying_ =
        glfwCreateWindow(static_cast<int>(width), static_cast<int>(height),
                         title, /*monitor*/ nullptr, /*share*/ nullptr);

    // Load GLEW.
    if (!framework.glew_loaded_) {
      glfwMakeContextCurrent(underlying_);
      glewExperimental = GL_TRUE;
      if (GLenum status = glewInit(); status == GLEW_OK) {
        STORM_INFO_("GLEW initialized!");
      } else {
        STORM_THROW_GL_("GLEW: failed to initialize, {:#x}: {}!", //
                        status, glewGetErrorString(status));
      }
      framework.glew_loaded_ = true;
      glfwMakeContextCurrent(nullptr);
    }

    // Set the callbacks.
    glfwSetWindowUserPointer(underlying_, this);
    glfwSetWindowCloseCallback(underlying_, &on_close_);
    glfwSetWindowSizeCallback(underlying_, &on_resize_);
    glfwSetKeyCallback(underlying_, &on_key_);
    glfwSetMouseButtonCallback(underlying_, &on_mouse_button_);
    glfwSetScrollCallback(underlying_, &on_scroll_);
    glfwSetCursorPosCallback(underlying_, &on_set_cursor_pos_);
  }

  /// @brief Run the window main loop.
  template<std::invocable RenderFn>
  void main_loop(RenderFn render_fn) {
    while (!glfwWindowShouldClose(underlying_)) {
      render_fn();
      glfwSwapBuffers(underlying_);
      glfwPollEvents();
    }
  }

  /// @brief Set handler for the window close event.
  template<std::invocable CloseFn>
  void on_close(CloseFn on_close_fn) {
    on_close_fns_.push_back(std::move(on_close_fn));
  }

  /// @brief Set handler for the window resize event.
  template<std::invocable<size_t, size_t> ResizeFn>
  void on_resize(ResizeFn on_resize_fn) {
    on_resize_fns_.push_back(std::move(on_resize_fn));
  }

  /// @brief Set handler for the key pressed event.
  /// @{
  template<std::invocable KeyDownFn>
  void on_key_down(Key key, KeyDownFn on_key_down_fn) {
    on_key_down(key, Modifiers::none, std::move(on_key_down_fn));
  }
  template<std::invocable KeyDownFn>
  void on_key_down(Key key, Modifiers mods, KeyDownFn on_key_down_fn) {
    on_key_down_fns_.emplace(std::pair{key, mods}, std::move(on_key_down_fn));
  }
  /// @}

  /// @brief Set handler for the key released event.
  /// @{
  template<std::invocable KeyUpFn>
  void on_key_up(Key key, KeyUpFn on_key_up_fn) {
    on_key_up(key, Modifiers::none, std::move(on_key_up_fn));
  }
  template<std::invocable KeyUpFn>
  void on_key_up(Key key, Modifiers mods, KeyUpFn on_key_up_fn) {
    on_key_up_fns_.emplace(std::pair{key, mods}, std::move(on_key_up_fn));
  }
  /// @}

  /// @brief Set handler for the mouse button pressed event.
  /// @{
  template<std::invocable MouseButtonDownFn>
  void on_mouse_button_down(MouseButton mouse_button,
                            MouseButtonDownFn on_mouse_button_down_fn) {
    on_mouse_button_down(mouse_button, Modifiers::none,
                         std::move(on_mouse_button_down_fn));
  }
  template<std::invocable MouseButtonDownFn>
  void on_mouse_button_down(MouseButton mouse_button, Modifiers mods,
                            MouseButtonDownFn on_mouse_button_down_fn) {
    on_mouse_button_down_fns_.emplace(std::pair{mouse_button, mods},
                                      std::move(on_mouse_button_down_fn));
  }
  /// @}

  /// @brief Set handler for the mouse button released event.
  /// @{
  template<std::invocable MouseButtonDownFn>
  void on_mouse_button_up(MouseButton mouse_button,
                          MouseButtonDownFn on_mouse_button_up_fn) {
    on_mouse_button_up(mouse_button, Modifiers::none,
                       std::move(on_mouse_button_up_fn));
  }
  template<std::invocable MouseButtonUpFn>
  void on_mouse_button_up(MouseButton mouse_button, Modifiers mods,
                          MouseButtonUpFn on_mouse_button_up_fn) {
    on_mouse_button_up_fns_.emplace(std::pair{mouse_button, mods},
                                    std::move(on_mouse_button_up_fn));
  }

  /// @brief Set handler for the scroll event.
  template<std::invocable<glm::dvec2> ScrollFn>
  void on_scroll(ScrollFn on_scroll_fn) {
    on_scroll_fns_.push_back(std::move(on_scroll_fn));
  }

  /// @brief Set handler for the set cursor position event.
  template<std::invocable<glm::dvec2> SetCursorPosFn>
  void on_set_cursor_pos(SetCursorPosFn on_set_cursor_pos_fn) {
    on_set_cursor_pos_fns_.push_back(std::move(on_set_cursor_pos_fn));
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
    for (const auto& close_fn : self->on_close_fns_) {
      close_fn();
    }
  }

  static void on_resize_(GLFWwindow* window, int width, int height) noexcept {
    const auto self = get_self_(window);
    for (const auto& on_resize_fn : self->on_resize_fns_) {
      on_resize_fn(static_cast<size_t>(width), static_cast<size_t>(height));
    }
  }

  static void on_key_(GLFWwindow* window, int key,
                      [[maybe_unused]] int scancode, int action,
                      int mods) noexcept {
    const auto self = get_self_(window);
    switch (action) {
      case GLFW_PRESS:
      case GLFW_REPEAT: {
        const auto [first, last] = self->on_key_down_fns_.equal_range(
            {static_cast<Key>(key), static_cast<Modifiers>(mods)});
        std::for_each(first, last, [](const auto& pair) { pair.second(); });
        break;
      }
      case GLFW_RELEASE: {
        const auto [first, last] = self->on_key_up_fns_.equal_range(
            {static_cast<Key>(key), static_cast<Modifiers>(mods)});
        std::for_each(first, last, [](const auto& pair) { pair.second(); });
        break;
      }
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
      case GLFW_REPEAT: {
        const auto [first, last] = self->on_mouse_button_down_fns_.equal_range(
            {static_cast<MouseButton>(mouse_button),
             static_cast<Modifiers>(mods)});
        std::for_each(first, last, [](const auto& pair) { pair.second(); });
        break;
      }
      case GLFW_RELEASE: {
        const auto [first, last] = self->on_mouse_button_up_fns_.equal_range(
            {static_cast<MouseButton>(mouse_button),
             static_cast<Modifiers>(mods)});
        std::for_each(first, last, [](const auto& pair) { pair.second(); });
        break;
      }
      default: //
        STORM_WARNING_("GLFW: unsupported mouse button action {:#x}!", action);
        break;
    }
  }

  static void on_scroll_(GLFWwindow* window, //
                         double x_offset, double y_offset) noexcept {
    const auto self = get_self_(window);
    for (const auto& on_scroll_fn : self->on_scroll_fns_) {
      on_scroll_fn({x_offset, y_offset});
    }
  }

  static void on_set_cursor_pos_(GLFWwindow* window, //
                                 double x_pos, double y_pos) noexcept {
    const auto self = get_self_(window);
    for (const auto& on_set_cursor_pos_fn : self->on_set_cursor_pos_fns_) {
      on_set_cursor_pos_fn({x_pos, y_pos});
    }
  }

}; // class Window

/// @brief RAII binder of an OpenGL window.
class BindWindow final {
public:

  /// @brief Bind the @p window.
  explicit BindWindow(const Window& window) {
    glfwMakeContextCurrent(window);
  }

  /// @brief Unbind the window.
  ~BindWindow() {
    glfwMakeContextCurrent(nullptr);
  }

}; // class BindWindow

} // namespace Storm::Vulture::gl
