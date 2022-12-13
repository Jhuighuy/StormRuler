// Copyright (C) 2020-2023 Oleg Butakov
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR Allocator PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT
// SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
// OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

#pragma once

#include <Storm/Base.hpp>

#include <Storm/Utils/Meta.hpp>

#include <istream>
#include <streambuf>

namespace Storm
{

// -----------------------------------------------------------------------------

/// @brief Filtering stream buffer.
/// @tparam BeginFilter First filtering character, typically '#'.
/// @tparam EndFilter Last filtering character, typically '\n'.
template<class Char, Char BeginFilter, Char EndFilter>
class FilteringStreambuf final : public std::basic_streambuf<Char>
{
private:

  std::istream* p_stream_;
  std::streambuf* p_streambuf_;
  Char buffer_;

  static constexpr auto begin_filter_ =
      std::basic_streambuf<Char>::traits_type::to_int_type(BeginFilter);
  static constexpr auto end_filter_ =
      std::basic_streambuf<Char>::traits_type::to_int_type(EndFilter);
  static constexpr auto eof_ = std::basic_streambuf<Char>::traits_type::eof();

public:

  /// @brief Construct a filtering stream buffer.
  explicit FilteringStreambuf(std::basic_istream<Char>& stream)
      : p_stream_{&stream}, p_streambuf_{p_stream_->rdbuf()}
  {
    p_stream_->rdbuf(this);
  }

  /// @brief Destroy a filtering stream buffer.
  ~FilteringStreambuf() final
  {
    p_stream_->rdbuf(p_streambuf_);
  }

protected:

  typename std::basic_streambuf<Char>::int_type underflow() final
  {
    auto bumped = p_streambuf_->sbumpc();
    if (bumped == begin_filter_) {
      while (bumped != eof_ && bumped != end_filter_) {
        bumped = p_streambuf_->sbumpc();
      }
    }
    if (bumped != eof_) {
      buffer_ = bumped;
      this->setg(&buffer_, &buffer_, &buffer_ + 1);
    }
    return bumped;
  }

}; // class FilteringStreambuf

// -----------------------------------------------------------------------------

} // namespace Storm
