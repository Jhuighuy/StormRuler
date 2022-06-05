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

#include <array>
#include <initializer_list>
#include <ranges>
#include <utility>

#include <stormBase.hxx>
#include <stormSolvers/Matrix.hxx>

namespace Storm {

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Statically-sized matrix.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Value, size_t NumRows, size_t NumCols>
class StaticMatrix : public BaseMatrix<StaticMatrix<Value, NumRows, NumCols>> {
private:

  std::array<std::array<Value, NumCols>, NumRows> Coeffs_;

public:

  /// @brief Default constructor.
  constexpr StaticMatrix() = default;

  /// @brief Construct a matrix with an @p initList.
  constexpr StaticMatrix(std::initializer_list<Value> initList) noexcept {
    StormAssert(initList.size() == size());
    std::ranges::copy(initList, data());
  }

  template<class T>
  constexpr explicit StaticMatrix(BaseMatrixView<T> const& mat) noexcept {
    *this = mat;
  }

  template<class T>
  constexpr auto& operator=(BaseMatrixView<T> const& mat) noexcept {
    this->assign(mat);
    return *this;
  }

  /// @brief Size of the matrix data.
  constexpr static auto size() noexcept {
    return std::integral_constant<size_t, NumRows * NumCols>{};
  }

  /// @brief Shape of the matrix.
  constexpr static auto shape() noexcept {
    return std::pair(std::integral_constant<size_t, NumRows>{},
                     std::integral_constant<size_t, NumCols>{});
  }

  /// @brief Pointer to the beginning of the matrix data.
  /// @{
  constexpr Value* data() noexcept {
    return Coeffs_[0].data();
  }
  constexpr Value const* data() const noexcept {
    return Coeffs_[0].data();
  }
  /// @}

  /// @brief Reference to the coefficient at @p rowIndex and @p colIndex.
  /// @{
  constexpr Value& operator()(size_t rowIndex, size_t colIndex = 0) noexcept {
    StormAssert(rowIndex < NumRows && colIndex < NumCols);
    return (Coeffs_[rowIndex])[colIndex];
  }
  constexpr Value const& operator()(size_t rowIndex,
                                    size_t colIndex = 0) const noexcept {
    return const_cast<StaticMatrix&>(*this)(rowIndex, colIndex);
  }
  /// @}

}; // class StaticMatrix

} // namespace Storm
