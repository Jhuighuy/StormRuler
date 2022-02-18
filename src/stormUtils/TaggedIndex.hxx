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

#include <stormBase.hxx>

namespace Storm {

template<class Tag, class Index = size_t>
using TaggedIndex = Index;

#if 0
template<class Tag, class Index = size_t>
class TaggedIndex {
private:
  Index Index_ = 0;

public:

  constexpr TaggedIndex() noexcept = default;
  
  constexpr TaggedIndex(Index index) noexcept :
    Index_{index} {
  }

  constexpr Index Get() const {
    return Index_;
  }

  constexpr auto operator<=>(TaggedIndex const& other) const noexcept = default;

  constexpr auto operator-(TaggedIndex const& other) const noexcept {
    return Index_ - other.Index_;
  }

  constexpr TaggedIndex& operator++() noexcept {
    return ++Index_, *this;
  }
  constexpr TaggedIndex& operator--() noexcept {
    return --Index_, *this;
  }

  constexpr auto operator++(int) noexcept {
    TaggedIndex copy(*this);
    return ++Index_, copy;
  }
  constexpr auto operator--(int) noexcept {
    TaggedIndex copy(*this);
    return --Index_, copy;
  }

}; // class TaggedIndex<...>
#endif

template<class>
class BaseIterator {

};

template<class Value, class Index = size_t>
class Array {
public:

  Index Size() const;

  Value const& operator[](Index index) const;

};

template<class Value, class ColumnIndex>
class TableRow {
public:

  ColumnIndex const* Begin() const;
  ColumnIndex const* End() const;

};

template<class Value, class RowIndex = size_t, class ColumnIndex = RowIndex>
class Table {
public:

  TableRow<Value, ColumnIndex> const& operator[](RowIndex index) const;

};

template<size_t Dim>
class Vec {};

} // namespace Storm
