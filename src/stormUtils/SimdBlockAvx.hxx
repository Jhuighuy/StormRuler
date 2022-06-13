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

#ifndef STORM_SIMD_BLOCK_GUARD_
#error #include <stormUtils/SimdBlock.hxx> instead
#endif

#include <immintrin.h>

#include <stormBase.hxx>
#include <stormUtils/Math.hxx>
#include <stormUtils/SimdBlock.hxx>

namespace Storm {

/// ----------------------------------------------------------------- ///
/// @brief @c double[4] SIMD block.
/// ----------------------------------------------------------------- ///
template<>
class alignas(32) SimdBlock<double, 4> {
public:

  __m256d values;

  STORM_FORCE_INLINE_
  auto& operator=(double a) noexcept {
    return values = _mm256_set1_pd(a), *this;
  }

}; // SimdBlock<double, 4>

template<>
class SimdBlockMulExpr<double, 4> {
public:

  SimdBlock<double, 4> a_arg, b_arg;

  STORM_FORCE_INLINE_
  operator SimdBlock<double, 4>() const noexcept {
    return {_mm256_mul_pd(a_arg.values, b_arg.values)};
  }

}; // SimdBlockMulExpr<double, 4>

STORM_FORCE_INLINE_
auto operator+(const SimdBlock<double, 4>& x) noexcept {
  return x;
}

STORM_FORCE_INLINE_
auto operator-(const SimdBlock<double, 4>& x) noexcept {
  const auto sign_mask{_mm256_castsi256_pd(_mm256_set1_epi64x(1ull << 63))};
  return SimdBlock<double, 4>{_mm256_xor_pd(sign_mask, x.values)};
}

/// @{
STORM_FORCE_INLINE_
auto operator*(double a, const SimdBlock<double, 4>& x) noexcept {
  return SimdBlockMulExpr<double, 4>{{_mm256_set1_pd(a)}, x};
}
STORM_FORCE_INLINE_
auto operator*(const SimdBlock<double, 4>& x, double a) noexcept {
  return SimdBlockMulExpr<double, 4>{{_mm256_set1_pd(a)}, x};
}
STORM_FORCE_INLINE_
auto& operator*=(SimdBlock<double, 4>& x, double a) noexcept {
  return x = x * a;
}
/// @}

/// @{
STORM_FORCE_INLINE_
auto operator/(const SimdBlock<double, 4>& x, double a) noexcept {
  return SimdBlockMulExpr<double, 4>{{_mm256_set1_pd(1.0 / a)}, x};
}
STORM_FORCE_INLINE_
auto& operator/=(SimdBlock<double, 4>& x, double a) noexcept {
  return x = x / a;
}
/// @}

/// @{
STORM_FORCE_INLINE_
auto operator+(const SimdBlock<double, 4>& x,
               const SimdBlock<double, 4>& y) noexcept {
  return SimdBlock<double, 4>{_mm256_add_pd(x.values, y.values)};
}
STORM_FORCE_INLINE_
auto operator+(const SimdBlock<double, 4>& x,
               const SimdBlockMulExpr<double, 4>& y_expr) noexcept {
  return SimdBlock<double, 4>{
      _mm256_fmadd_pd(y_expr.a_arg.values, y_expr.b_arg.values, x.values)};
}
STORM_FORCE_INLINE_
auto operator+(const SimdBlockMulExpr<double, 4>& x_expr,
               const SimdBlock<double, 4>& y) noexcept {
  return SimdBlock<double, 4>{
      _mm256_fmadd_pd(x_expr.a_arg.values, x_expr.b_arg.values, y.values)};
}
STORM_FORCE_INLINE_
auto operator+(const SimdBlockMulExpr<double, 4>& x_expr,
               const SimdBlockMulExpr<double, 4>& y_expr) noexcept {
  return static_cast<SimdBlock<double, 4>>(x_expr) + y_expr;
}
STORM_FORCE_INLINE_
auto& operator+=(SimdBlock<double, 4>& x, const auto& y_any) noexcept {
  return x = x + y_any;
}
/// @}

/// @{
STORM_FORCE_INLINE_
auto operator-(const SimdBlock<double, 4>& x,
               const SimdBlock<double, 4>& y) noexcept {
  return SimdBlock<double, 4>{_mm256_sub_pd(x.values, y.values)};
}
STORM_FORCE_INLINE_
auto operator-(const SimdBlock<double, 4>& x,
               const SimdBlockMulExpr<double, 4>& y_expr) noexcept {
  return SimdBlock<double, 4>{
      _mm256_fnmadd_pd(y_expr.a_arg.values, y_expr.b_arg.values, x.values)};
}
STORM_FORCE_INLINE_
auto operator-(const SimdBlockMulExpr<double, 4>& x_expr,
               const SimdBlock<double, 4>& y) noexcept {
  return SimdBlock<double, 4>{
      _mm256_fmsub_pd(x_expr.a_arg.values, x_expr.b_arg.values, y.values)};
}
STORM_FORCE_INLINE_
auto operator-(const SimdBlockMulExpr<double, 4>& x_expr,
               const SimdBlockMulExpr<double, 4>& y_expr) noexcept {
  return static_cast<SimdBlock<double, 4>>(x_expr) - y_expr;
}
STORM_FORCE_INLINE_
auto& operator-=(SimdBlock<double, 4>& x,
                 const SimdBlock<double, 4>& y) noexcept {
  return x = x - y;
}
/// @}

namespace math {

  STORM_FORCE_INLINE_
  auto abs(const SimdBlock<double, 4>& x) noexcept {
    const auto sign_mask{_mm256_castsi256_pd(_mm256_set1_epi64x(1ull << 63))};
    return SimdBlock<double, 4>{_mm256_andnot_pd(sign_mask, x.values)};
  }

  STORM_FORCE_INLINE_
  auto sqrt(const SimdBlock<double, 4>& x) noexcept {
    return SimdBlock<double, 4>{_mm256_sqrt_pd(x.values)};
  }

} // namespace math

} // namespace Storm
