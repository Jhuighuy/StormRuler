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

#pragma once

#include <concepts>
#include <string_view>

#include <stormBase.hxx>

namespace Storm {

/// @brief Spawns a reflectible enumeration body.
/// @param Class The current class name.
#define StormEnum_(Class)                                                      \
private:                                                                       \
  using Class_ = Class;                                                        \
  friend class Enum<Class>;                                                    \
  static constexpr size_t BaseCounter_{__COUNTER__};                           \
                                                                               \
  template<size_t Index_, class Func>                                          \
  struct ForEachValueImpl_ {                                                   \
    static constexpr void call_(Func const&) noexcept {}                       \
  };                                                                           \
  template<class Func>                                                         \
  static constexpr void forEachValue_(Func const& func) noexcept {             \
    ForEachValueImpl_<0, Func>::call_(func);                                   \
  }                                                                            \
                                                                               \
public:                                                                        \
  /** @brief Construct the enumeration from the @p value. */                   \
  template<std::integral Integer>                                              \
  constexpr explicit Class(Integer value) : Enum<Class>(value) {}

/// @brief Spawns a reflectible enumeration value.
#define StormEnumValue_(Name, ...) StormEnumValueS_(Name, ##__VA_ARGS__, #Name)

/// @brief Spawns a names reflectible enumeration value.
#define StormEnumValueS_(Name, String, ...)                                    \
  static constexpr Enum<Class_> Name{__COUNTER__ - BaseCounter_ - 1};          \
                                                                               \
private:                                                                       \
  template<class Func>                                                         \
  struct ForEachValueImpl_<static_cast<size_t>(Name), Func> {                  \
    static constexpr void call_(Func const& func) noexcept {                   \
      func(Name, String);                                                      \
      ForEachValueImpl_<static_cast<size_t>(Name) + 1, Func>::call_(func);     \
    }                                                                          \
  };                                                                           \
                                                                               \
public:

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Reflectible enumeration.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Derived, std::integral Underlying = int>
class Enum {
private:
  Underlying Value_;

public:
  /// @brief Construct the enumeration from @p value.
  template<std::integral Integer>
  constexpr explicit Enum(Integer value = {}) noexcept :
      Value_(static_cast<Underlying>(value)) {}

  /// @brief Convert the enumeration into the derived type.
  /// @{
  constexpr operator Derived&() noexcept {
    return static_cast<Derived&>(*this);
  }
  constexpr operator Derived const&() const noexcept {
    return static_cast<Derived const&>(*this);
  }
  /// @}

  /// @brief Convert the enumeration into the integral type.
  template<std::integral Integer>
  constexpr explicit operator Integer() const noexcept {
    return static_cast<Integer>(Value_);
  }

  /// @brief Enumeration comparison operators.
  constexpr auto operator<=>(Enum const& other) const noexcept = default;

  /// @brief Convert the enumeration into the string.
  constexpr std::string_view toString() const noexcept;

  /// @brief Convert the specified @p string to the enumeration.
  static constexpr Enum fromString(std::string_view string);

}; // class Enum

template<class Derived, std::integral Underlying>
constexpr std::string_view
Enum<Derived, Underlying>::toString() const noexcept {
  std::string_view result{};

  Derived::forEachValue_([&](Enum<Derived> value, std::string_view string) {
    if (*this == value) result = string;
  });

  return result;

} // Enum::ToString

template<class Derived, std::integral Underlying>
constexpr Enum<Derived, Underlying>
Enum<Derived, Underlying>::fromString(std::string_view string) {
  bool found = false;
  Enum<Derived, Underlying> result;

  Derived::forEachValue_([&](Enum<Derived> value, std::string_view theString) {
    if (string == theString) result = value, found = true;
  });
  if (!found) { throw std::invalid_argument("Invalid enum string"); }

  return result;

} // Enum::FromString

} // namespace Storm
