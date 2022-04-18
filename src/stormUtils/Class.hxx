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

#include <memory>
#include <complex>
#include <string_view>
#include <string>
#include <vector>
#include <type_traits>

#include <stormBase.hxx>
#include <stormUtils/Enum.hxx>

namespace Storm {

/// ----------------------------------------------------------------- ///
/// @brief Macro for spawning the reflectible class body.
/// ----------------------------------------------------------------- ///
#define STORM_CLASS_(Class, ...) \
  private: \
    \
    using Class_ = Class; \
    static constexpr size_t BaseCounter_{__COUNTER__}; \
    \
    template<size_t Index, class Visitor> \
    struct ReflectionChain_ { \
      static constexpr void Invoke_(Class_&, Visitor const&) {} \
    }; \
    \
  public: \
    //\
    //void Reflect(Visitor const& visitor) noexcept override { \
    //  __VA_ARGS__::Reflect(visitor); \
    //  ReflectionChain_<0, Visitor>::Invoke_(*this, visitor); \
    //}

/// ----------------------------------------------------------------- ///
/// @brief Macro for spawning a reflectible class field body.
/// ----------------------------------------------------------------- ///
#define STORM_FIELD_(Type, Name, ...) \
    Type Name __VA_ARGS__; \
    \
  private: \
    \
    static constexpr size_t Name##Index_{__COUNTER__ - BaseCounter_ - 1}; \
    \
    template<class Visitor> \
    struct ReflectionChain_<Name##Index_, Visitor> { \
      static constexpr void Invoke_(Class_& self, Visitor const& visitor) { \
        visitor(#Name, self.Name); \
        ReflectionChain_<Name##Index_ + 1, Visitor>::Invoke_(self, visitor); \
      } \
    }; \
    \
  public:

//struct FieldInfo;
using FieldInfo = std::string_view;

/// ----------------------------------------------------------------- ///
/// @brief Base reflectible object.
/// ----------------------------------------------------------------- ///
class Object {
public:

  virtual ~Object() = default;

  /// @brief Reflect the object with @p visitor.
  //virtual void Reflect(class Visitor const&) noexcept {}

}; // class Object

/// ----------------------------------------------------------------- ///
/// @brief Visitor for object reflection.
/// ----------------------------------------------------------------- ///
class Visitor {
public:

  /// @brief Visit a boolean field.
  void operator()(FieldInfo const& fieldInfo, bool& value) const;

  /// @brief Visit an integer field.
  template<class Integer>
    requires (std::is_integral_v<Integer>)
  void operator()(FieldInfo const& fieldInfo, Integer& value) const;

  /// @brief Visit a floating point field.
  template<class Number>
    requires (std::is_floating_point_v<Number>)
  void operator()(FieldInfo const& fieldInfo, Number& value) const;

  /// @brief Visit a complex floating point field.
  template<class Number>
    requires (std::is_floating_point_v<Number>)
  void operator()(FieldInfo const& fieldInfo, 
                  std::complex<Number>& value) const;

  /// @brief Visit a enumeration field.
  template<class Derived, class Underlying>
  void operator()(FieldInfo const& fieldInfo, 
                  Enum<Derived, Underlying>& value) const;

  /// @brief Visit a string field.
  void operator()(FieldInfo const& fieldInfo, std::string& value) const;

  /// @brief Visit a vector field.
  template<class Value>
  void operator()(FieldInfo const& fieldInfo, 
                  std::vector<Value>& value) const;

  /// @brief Visit an object field.
  /// @{
  template<class Derived>
    requires (std::is_base_of_v<Object, Derived>)
  void operator()(FieldInfo const& fieldInfo, Derived& value) const;
  template<class Derived>
    requires (std::is_base_of_v<Object, Derived>)
  void operator()(FieldInfo const& fieldInfo, 
                  std::unique_ptr<Derived>& value) const;
  template<class Derived>
    requires (std::is_base_of_v<Object, Derived>)
  void operator()(FieldInfo const& fieldInfo, 
                  std::shared_ptr<Derived>& value) const;
  /// @}

}; // class Visitor

} // namespace Storm
