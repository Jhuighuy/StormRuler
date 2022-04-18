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

#include <string_view>
#include <type_traits>

#include <stormBase.hxx>

namespace Storm {

/// ----------------------------------------------------------------- ///
/// @brief Macro for spawning a reflectible enumeration body.
/// ----------------------------------------------------------------- ///
#define STORM_ENUM_(Class) \
  private: \
    \
    using Class_ = Class; \
    friend class Enum<Class>; \
    static constexpr size_t BaseCounter_{__COUNTER__}; \
    \
    template<size_t Index_, class Func> \
    struct ReflectionChain_ { \
      static constexpr void Invoke_(Func const&) noexcept {} \
    }; \
    template<class Func> \
    static constexpr void Reflect_(Func const& func) noexcept { \
      ReflectionChain_<0, Func>::Invoke_(func); \
    } \
    \
  public: \
    \
    /** @brief Construct the enumeration from the integer. */ \
    template<class Integer> \
      requires (std::is_integral_v<Integer>) \
    constexpr explicit Class(Integer value) : Enum<Class>(value) { \
    }

/// ----------------------------------------------------------------- ///
/// @brief Macro for spawning a reflectible enumeration value body.
/// ----------------------------------------------------------------- ///
#define STORM_ENUM_VALUE_(Name, ...) \
    STORM_ENUM_VALUE__(Name, ##__VA_ARGS__, #Name)

#define STORM_ENUM_VALUE__(Name, String, ...) \
    static constexpr Enum<Class_> Name{__COUNTER__ - BaseCounter_ - 1}; \
    \
  private: \
    \
    template<class Func> \
    struct ReflectionChain_<static_cast<size_t>(Name), Func> { \
      static constexpr void Invoke_(Func const& func) noexcept { \
        func(Name, String); \
        ReflectionChain_<static_cast<size_t>(Name) + 1, Func>::Invoke_(func); \
      } \
    }; \
    \
  public:

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Reflectible enumeration.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Derived, class Underlying = int>
  requires (std::is_integral_v<Underlying>)
class Enum {
private:
  Underlying Value_;

public:

  /// @brief Construct the enumeration from the integer.
  template<class Integer = Underlying>
    requires (std::is_integral_v<Integer>)
  constexpr explicit Enum(Integer value = {}) : 
      Value_(static_cast<Underlying>(value)) {
  }

  /// @brief Convert the enumeration into the derived type.
  constexpr operator Derived() const noexcept {
    return static_cast<Derived const&>(*this);
  }

  /// @brief Convert the enumeration into the integral type.
  template<class Integer>
    requires (std::is_integral_v<Integer>)
  constexpr explicit operator Integer() const noexcept {
    return static_cast<Integer>(Value_);
  }

  /// @brief Enumeration comparison operators.
  /// @{
  constexpr bool operator==(Enum const& other) const noexcept {
    return Value_ == other.Value_;
  }
  constexpr bool operator!=(Enum const& other) const noexcept {
    return Value_ != other.Value_;
  }
  /// @}

  /// @brief Convert the enumeration into the string.
  constexpr std::string_view ToString() const noexcept;

  /// @brief Convert the specified @p string to the enumeration. 
  static constexpr Enum<Derived, Underlying> FromString(std::string_view string);

}; // class Enum<...>

template<class Derived, class Underlying>
constexpr std::string_view 
    Enum<Derived, Underlying>::ToString() const noexcept {

  std::string_view result{};

  Derived::Reflect_(
    [&](Enum<Derived> value, std::string_view string) {
      if (*this == value) result = string;
    });

  return result;

} // Enum<...>::ToString

template<class Derived, class Underlying>
constexpr Enum<Derived, Underlying> 
    Enum<Derived, Underlying>::FromString(std::string_view string) {

  bool found = false;
  Enum<Derived, Underlying> result;

  Derived::Reflect_(
    [&](Enum<Derived> value, std::string_view theString) {
      if (string == theString) result = value, found = true;
    });
  if (!found) {
    throw std::invalid_argument("Invalid enum string");
  }

  return result;

} // Enum<...>::FromString

} // namespace Storm
