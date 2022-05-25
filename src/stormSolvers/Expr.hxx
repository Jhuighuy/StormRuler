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

#include <tuple>
#include <type_traits>

namespace Storm {

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Expression tree node.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Callable, class... Operands>
class Expr {
private:
  Callable Callable_;
  std::tuple<Operands const&...> Operands_;

public:

  /// @brief Construct an expression with \
  ///   @p callable and @p operadands.
  constexpr explicit Expr(Callable&& callable,
                          Operands const&... operands) noexcept :
    Callable_{std::forward<Callable>(callable)}, Operands_{operands...} {
  }

  /// @brief Evaluate the expression at the @p indices.
  template<class... Indices>
  constexpr auto operator()(Indices... indices) const noexcept {
    return std::apply(
      [this, indices...](Operands const&... operands) {
        return Callable_(operands(indices...)...);
      }, Operands_);
  }

}; // class Expr

template<class Any>
struct IsExpr_t : std::false_type {};
template<class Callable, class... Operands>
struct IsExpr_t<Expr<Callable, Operands...>> : std::true_type {};

template<class Any>
constexpr bool IsExpr = IsExpr_t<Any>::value;

#define StormUnaryExpr_(OP) \
  template<class E_> \
    requires (IsExpr<E_>) \
  constexpr auto operator OP(E_ const& e) noexcept { \
    return Expr{[](auto const& e) { return OP e; }, e}; \
  }

#define StormBinaryExpr_(OP) \
  /** @{ */ \
  template<class L_, class R_> \
    requires (IsExpr<L_> || IsExpr<R_>) \
  constexpr auto operator OP(L_ const& l, R_ const& r) noexcept { \
    if constexpr (!IsExpr<L_>) { \
      return Expr([l](auto const& r) { return l OP r; }, r); \
    } else if constexpr (!IsExpr<R_>) { \
      return Expr([r](auto const& l) { return l OP r; }, l); \
    } else { \
      return Expr([](auto const& l, auto const& r) { return l OP r; }, l, r); \
    } \
  } \
  /** @} */

/// @brief "Unary plus" an expression @p e. 
StormUnaryExpr_(+)

/// @brief Negate an expression @p e. 
StormUnaryExpr_(-)

/// @brief Component-wise add expressions @p l and @p r.
StormBinaryExpr_(+)

/// @brief Component-wise subtract expressions @p l and @p r.
StormBinaryExpr_(-)

/// @brief Component-wise multiply expressions @p l and @p r.
StormBinaryExpr_(*)

/// @brief Component-wise divide expressions @p l and @p r.
StormBinaryExpr_(/)

#undef StormUnaryExpr_
#undef StormBinaryExpr_

} // namespace Storm
