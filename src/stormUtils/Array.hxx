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

#include <vector>

#include <stormBase.hxx>

namespace Storm {

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Dynamic array class, a wrapper for @c std::vector \
///   class with a support of the custom index type.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value, class Index = size_t>
class Array {
private:
  std::vector<Value> Vector_;

public:

  Array() noexcept = default;

  Array(Array&&) noexcept = default;

  Array(Array const&) = default;

  explicit Array(Index size) : 
      Vector_(static_cast<size_t>(size)) {
  }
  Array(Index size, Value const& value) : 
      Vector_(static_cast<size_t>(size), value) {
  }

  ~Array() noexcept = default;

  template<class Iterator>
  Array(Iterator first, Iterator last) : 
      Vector_(first, last) {
  }

  Array(std::initializer_list<Value> list) : 
      Vector_(list) {
  }

  Array& operator=(Array&&) noexcept = default;
  
  Array& operator=(Array const&) = default;

  Array& operator=(std::initializer_list<Value> list) {
    return Vector_ = list, *this;
  }

  void Assign(Index size, Value const& value) {
    Vector_.assign(static_cast<size_t>(size), value);
  }

  template<class Iterator>
  void Assign(Iterator first, Iterator last) {
    Vector_.assign(first, last);
  }

  void Assign(std::initializer_list<Value> list) {
    Vector_.assign(list);
  }

  auto Begin() noexcept {
    return Vector_.begin();
  }
  auto Begin() const noexcept {
    return Vector_.begin();
  }

  auto End() noexcept {
    return Vector_.end();
  }
  auto End() const noexcept {
    return Vector_.end();
  }

  Value& Front() noexcept {
    return Vector_.front();
  }
  Value const& Front() const noexcept {
    return Vector_.front();
  }

  Value& Back() noexcept {
    return Vector_.back();
  }
  Value const& Back() const noexcept {
    return Vector_.back();
  }

  Value& operator[](Index index) noexcept {
    return Vector_[static_cast<size_t>(index)];
  }
  Value const& operator[](Index index) const noexcept {
    return Vector_[static_cast<size_t>(index)];
  }

  Value* Data() noexcept {
    return Vector_.data();
  }
  Value const* Data() const noexcept {
    return Vector_.data();
  }

  bool Empty() const noexcept {
    return Vector_.empty();
  }

  Index Size() const noexcept {
    return Index(Vector_.size());
  }

  Index Capacity() const noexcept {
    return Index(Vector_.capacity());
  }

  void Resize(Index size) {
    Vector_.resize(static_cast<size_t>(size));
  }

  void Reserve(Index capacity) {
    Vector_.reserve(static_cast<size_t>(capacity));
  }

  void Clear() {
    Vector_.clear();
  }

  void ShrinkToFit() {
    Vector_.shrink_to_fit();
  }

  void PushBack(Value&& value) {
    Vector_.push_back(std::forward<Value>(value));
  }
  void PushBack(Value const& value) {
    Vector_.push_back(value);
  }

  template<class... Args>
  Value& EmplaceBack(Args&&... args) {
    return Vector_.emplace_back(std::forward<Args>(args)...);
  }

  void PopBack() {
    Vector_.pop_back();
  }

  /// @todo Insert and Erase are missing, \
  ///   implement them on demand.

}; // class Array<...>

} // namespace Storm
