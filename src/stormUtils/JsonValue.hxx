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
#include <string>
#include <string_view>
#include <vector>
#include <ranges>
#include <map>

#include <stormBase.hxx>
#include <stormUtils/Enum.hxx>

namespace Storm {

/// ----------------------------------------------------------------- ///
/// @brief Abstract JSON value. 
/// ----------------------------------------------------------------- ///
class JsonValue : public std::enable_shared_from_this<JsonValue> {
public:

  /// @brief Default virtual destructor.
  virtual ~JsonValue() = default;

}; // class JsonValue

/// ----------------------------------------------------------------- ///
/// @brief Boolean JSON value.
/// ----------------------------------------------------------------- ///
class JsonBoolValue : public JsonValue {
public:

  /// @brief Boolean value.
  bool Value{};

  /// @brief Construct a boolean value.
  explicit JsonBoolValue(bool value) : Value(value) {}

}; // class JsonBoolValue

inline std::shared_ptr<JsonValue> MakeJsonValue(bool value) {

  return std::make_shared<JsonBoolValue>(value);

} // MakeJsonValue<bool>

/// ----------------------------------------------------------------- ///
/// @brief Numeric JSON value.
/// ----------------------------------------------------------------- ///
class JsonNumericValue : public JsonValue {
public:

  /// @brief Numeric value.
  long double Value{};

  /// @brief Construct a numeric value.
  /// @{
  explicit JsonNumericValue(long double value) : Value(value) {}
  template<class Integer>
    requires (std::is_integral_v<Integer>)
  explicit JsonNumericValue(Integer value); /// @todo
  /// @}

}; // class JsonNumericValue

inline std::shared_ptr<JsonValue> MakeJsonValue(long double value) {

  return std::make_shared<JsonNumericValue>(value);

} // MakeJsonValue<long double>

template<class Integer>
  requires (std::is_integral_v<Integer>)
inline std::shared_ptr<JsonValue> MakeJsonValue(Integer value) {

  return std::make_shared<JsonNumericValue>(value);

} // MakeJsonValue<Integer>

/// ----------------------------------------------------------------- ///
/// @brief String JSON value.
/// ----------------------------------------------------------------- ///
class JsonStringValue : public JsonValue {
public:

  /// @brief String value.
  std::string Value;

  /// @brief Construct a string value.
  explicit JsonStringValue(std::string_view value) : Value(value) {}

}; // class JsonStringValue

inline std::shared_ptr<JsonValue> MakeJsonValue(std::string_view value) {

  return std::make_shared<JsonStringValue>(value);

} // MakeJsonValue<std::string_view>

template<class Derived, class Underlying>
inline std::shared_ptr<JsonValue> MakeJsonValue(Enum<Derived, Underlying> value) {

  return std::make_shared<JsonStringValue>(value.ToString());

} // MakeJsonValue<Enum>

/// ----------------------------------------------------------------- ///
/// @brief Array JSON value.
/// ----------------------------------------------------------------- ///
class JsonArrayValue : public JsonValue {
public:

  /// @brief Array of inner values.
  std::vector<std::shared_ptr<JsonValue>> Values;

}; // class JsonArrayValue

template<class Range>
inline std::shared_ptr<JsonValue> MakeJsonValue(Range const& range) {

  auto array = std::make_shared<JsonArrayValue>();
  for (auto const& value : range) {
    array->Values.push_back(MakeJsonValue(value));
  }
  return array;

} // MakeJsonValue<Enum>

/// ----------------------------------------------------------------- ///
/// @brief Null JSON value.
/// ----------------------------------------------------------------- ///
class JsonNullValue : public JsonValue {};

/// ----------------------------------------------------------------- ///
/// @brief Object JSON value.
/// ----------------------------------------------------------------- ///
class JsonObjectValue : public JsonValue {
public:

  /// @brief Dictionary of fields.
  std::map<std::string, std::shared_ptr<JsonValue>> Fields;

}; // class JsonObjectValue

} // namespace Storm
