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
template<class Value, size_t num_rows, size_t num_cols>
class StaticMatrix :
    public BaseMatrix<StaticMatrix<Value, num_rows, num_cols>> {
private:

  std::array<std::array<Value, num_cols>, num_rows> Coeffs_;

public:

  /// @brief Default constructor.
  constexpr StaticMatrix() = default;

  /// @brief Construct a matrix with an @p initList.
  constexpr StaticMatrix(std::initializer_list<Value> initList) noexcept {
    STORM_ASSERT_(initList.size() == size());
    std::ranges::copy(initList, data());
  }

  template<class T>
  constexpr explicit StaticMatrix(const BaseMatrixView<T>& mat) noexcept {
    *this = mat;
  }

  template<class T>
  constexpr auto& operator=(const BaseMatrixView<T>& mat) noexcept {
    this->assign(mat);
    return *this;
  }

  /// @brief Size of the matrix data.
  constexpr static auto size() noexcept {
    return std::integral_constant<size_t, num_rows * num_cols>{};
  }

  /// @brief Shape of the matrix.
  constexpr static auto shape() noexcept {
    return std::pair(std::integral_constant<size_t, num_rows>{},
                     std::integral_constant<size_t, num_cols>{});
  }

  /// @brief Pointer to the beginning of the matrix data.
  /// @{
  constexpr Value* data() noexcept {
    return Coeffs_[0].data();
  }
  constexpr const Value* data() const noexcept {
    return Coeffs_[0].data();
  }
  /// @}

  /// @brief Reference to the coefficient at @p row_index and @p col_index.
  /// @{
  constexpr Value& operator()(size_t row_index, size_t col_index = 0) noexcept {
    STORM_ASSERT_(row_index < num_rows && col_index < num_cols);
    return (Coeffs_[row_index])[col_index];
  }
  constexpr const Value& operator()(size_t row_index,
                                    size_t col_index = 0) const noexcept {
    return const_cast<StaticMatrix&>(*this)(row_index, col_index);
  }
  /// @}

}; // class StaticMatrix

} // namespace Storm
